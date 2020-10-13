
# runtime: ['intake', 'zarr', 'numpy', 'pandas', 'scikit-allel', 'requests', 'aiohttp', 'gcsfs']
# tests: ['nose >= 1.3', 'ddt']


import argparse
import io
import unittest
from io import StringIO
from unittest.mock import Mock, patch, MagicMock, call

import allel
import pandas
import zarr
from ddt import ddt, data

from scripts.vgc import ConcordanceResult, classify, classify_sample, classify_chromosome, \
    compute_hom_alt_mismatch, HOMOZIGOUS_ALTERNATE_CATEGORY, MISMATCH_CATEGORY, compute_het_mismatch, \
    HETEROZIGOUS_CATEGORY, Callset, FilteringCallset, Commands, ArgumentParserBuilder, Summarizer, to_filtered_callset

A_URL = "A_URL"
A_FILTER_PATH = "A_FILTER_PATH"

SAMPLE_A = 'A'
SAMPLE_B = 'B'
SOME_SAMPLES = [SAMPLE_A, SAMPLE_B]

CHROMOSOME_1 = "YL"
CHROMOSOME_2 = "YR"
SOME_CHROMOSOMES = [CHROMOSOME_1, CHROMOSOME_2]

CHROMOSOME_1_GT = [[0, 1], [1, 1]]
CHROMOSOME_2_GT = [[0, 0], [-1, -1], [2, 0]]
SOME_CHROMOSOMES_GT = {CHROMOSOME_1: CHROMOSOME_1_GT, CHROMOSOME_2: CHROMOSOME_2_GT}

CHROMOSOME_1_SITE_FILTERS = [False, True]
CHROMOSOME_2_SITE_FILTERS = [True, False, True]
SOME_CHROMOSOMES_FILTERS = {CHROMOSOME_1: CHROMOSOME_1_SITE_FILTERS, CHROMOSOME_2: CHROMOSOME_2_SITE_FILTERS}


class TestSampleConcordanceResult(unittest.TestCase):

    def test_add_and_print(self):
        r = ConcordanceResult()
        recorder = r.record_chromosome_statistics("SAMPLE", "CHROMOSOME")
        recorder.record("LEFT", "RIGHT", 20)
        recorder.record("LEFT1", "RIGHT1", 201)
        recorder = r.record_chromosome_statistics("SAMPLE2", "OTHER_CHROMOSOME")
        recorder.record("OLEFT", "ORIGHT", 10)

        stream = io.StringIO()
        r.print(stream)
        self.assertEqual("""sample,chromosome,control,test,count
SAMPLE,CHROMOSOME,LEFT,RIGHT,20
SAMPLE,CHROMOSOME,LEFT1,RIGHT1,201
SAMPLE2,OTHER_CHROMOSOME,OLEFT,ORIGHT,10
""", stream.getvalue())


class TestClassify(unittest.TestCase):
    def setUp(self) -> None:
        self.control_creator = Mock()
        self.control_creator.return_value = MagicMock()
        self.test_creator = Mock()
        self.test_creator.return_value = MagicMock()
        self.results = MagicMock()

    @patch('scripts.vgc.classify_sample')
    def test_no_sample(self, classify_mock):
        classify([], SOME_CHROMOSOMES, self.control_creator, self.test_creator, self.results)
        self.control_creator.assert_not_called()
        self.test_creator.assert_not_called()
        classify_mock.assert_not_called()

    @patch('scripts.vgc.classify_sample')
    def test_with_sample_and_chromosomes(self, classify_mock):
        classify(SOME_SAMPLES, SOME_CHROMOSOMES, self.control_creator, self.test_creator, self.results)
        self.assertCountEqual(classify_mock.call_args_list,
                              [call(SAMPLE_A, SOME_CHROMOSOMES, self.control_creator.return_value,
                                    self.test_creator.return_value, self.results),
                               call(SAMPLE_B, SOME_CHROMOSOMES, self.control_creator.return_value,
                                    self.test_creator.return_value, self.results)])

    @patch('scripts.vgc.classify_sample')
    def test_no_chromosomes(self, classify_mock):
        classify(SOME_SAMPLES, [], self.control_creator, self.test_creator, self.results)
        self.assertCountEqual(classify_mock.call_args_list,
                              [call(SAMPLE_A, [], self.control_creator.return_value, self.test_creator.return_value,
                                    self.results),
                               call(SAMPLE_B, [], self.control_creator.return_value, self.test_creator.return_value,
                                    self.results)])


class TestClassifySample(unittest.TestCase):
    def setUp(self) -> None:
        self.control = MagicMock()
        self.test = MagicMock()
        self.results = MagicMock()
        self.record_chromosome_statistics_1 = MagicMock()
        self.record_chromosome_statistics_2 = MagicMock()

        def side_effect(sample, chromosome):
            if sample == SAMPLE_A and chromosome == CHROMOSOME_1:
                return self.record_chromosome_statistics_1
            if sample == SAMPLE_A and chromosome == CHROMOSOME_2:
                return self.record_chromosome_statistics_2
            return None

        self.results.record_chromosome_statistics.side_effect = side_effect
        self.control_gt_chromosome_1 = MagicMock()
        self.control_gt_chromosome_2 = MagicMock()
        self.test_gt_chromosome_1 = MagicMock()
        self.test_gt_chromosome_2 = MagicMock()
        control_map = {CHROMOSOME_1: self.control_gt_chromosome_1, CHROMOSOME_2: self.control_gt_chromosome_2}
        self.control.gt.side_effect = lambda chrom: control_map[chrom]
        test_map = {CHROMOSOME_1: self.test_gt_chromosome_1, CHROMOSOME_2: self.test_gt_chromosome_2}
        self.test.gt.side_effect = lambda chrom: test_map[chrom]

    @patch('scripts.vgc.classify_chromosome')
    def test_no_chromosomes(self, classify_mock):
        classify_sample(SAMPLE_A, [], self.control, self.test, self.results)
        classify_mock.assert_not_called()

    @patch('scripts.vgc.classify_chromosome')
    def test_with_chromosomes(self, classify_mock):
        classify_sample(SAMPLE_A, SOME_CHROMOSOMES, self.control, self.test, self.results)
        self.assertCountEqual(classify_mock.call_args_list,
                              [call(control=self.control_gt_chromosome_1, test=self.test_gt_chromosome_1,
                                    recorder=self.record_chromosome_statistics_1),
                               call(control=self.control_gt_chromosome_2, test=self.test_gt_chromosome_2,
                                    recorder=self.record_chromosome_statistics_2)])


class TestClassifyChromosome(unittest.TestCase):
    def setUp(self) -> None:
        self.control = MagicMock()
        self.test = MagicMock()
        self.recorder = MagicMock()
        self.het_mismatch = 12
        self.hom_alt_mismatch = 9

    @patch('scripts.vgc._classify_non_mismatch')
    @patch('scripts.vgc.compute_hom_alt_mismatch')
    @patch('scripts.vgc.compute_het_mismatch')
    def test_classify_chromosome(self, het_mock, hom_alt_mock, non_mismatch_mock):
        het_mock.return_value = self.het_mismatch
        hom_alt_mock.return_value = self.hom_alt_mismatch
        classify_chromosome(control=self.control, test=self.test, recorder=self.recorder)
        self.assertCountEqual(het_mock.call_args_list,
                              [call(self.control, self.test, self.recorder)])
        self.assertCountEqual(hom_alt_mock.call_args_list,
                              [call(self.control, self.test, self.recorder)])
        self.assertCountEqual(non_mismatch_mock.call_args_list,
                              [call(self.het_mismatch, self.hom_alt_mismatch, self.control, self.test, self.recorder)])


@ddt
class TestMismatchChromosome(unittest.TestCase):
    def setUp(self) -> None:
        self.recorder = MagicMock()

    @data(
        ([-1, -1], [0, 1], 0),
        ([0, 0], [0, 1], 0),
        ([1, 1], [0, 1], 0),
        ([0, 1], [-1, -1], 0),
        ([0, 1], [0, 0], 0),
        ([0, 1], [1, 1], 0),
        ([0, 1], [0, 1], 0),
        ([0, 1], [0, 2], 1),
        ([0, 1], [0, 3], 1),
        ([0, 2], [0, 1], 1),
        ([0, 2], [0, 2], 0),
        ([1, 2], [1, 2], 0),
        ([1, 2], [1, 3], 1),
    )
    def test_compute_het_mismatch(self, value):
        control, test, expected = value
        lgt = allel.GenotypeVector([control])
        rgt = allel.GenotypeVector([test])
        actual = compute_het_mismatch(control=lgt, test=rgt, recorder=self.recorder)
        self.assertEqual(expected, actual)
        self.assertCountEqual(self.recorder.record.call_args_list,
                              [call(HETEROZIGOUS_CATEGORY, MISMATCH_CATEGORY, expected)])

    @data(
        ([0, 0], [0, 0], 0),
        ([-1, -1], [1, 1], 0),
        ([0, 0], [1, 1], 0),
        ([0, 1], [1, 1], 0),
        ([1, 1], [-1, -1], 0),
        ([1, 1], [0, 0], 0),
        ([1, 1], [0, 1], 0),
        ([1, 1], [1, 1], 0),
        ([2, 2], [2, 2], 0),
        ([3, 3], [3, 3], 0),
        ([1, 1], [2, 2], 1),
        ([1, 1], [3, 3], 1),
        ([2, 2], [1, 1], 1),
        ([2, 2], [3, 3], 1),
        ([3, 3], [1, 1], 1),
        ([3, 3], [2, 2], 1),
    )
    def test_compute_hom_alt_mismatch(self, value):
        control, test, expected = value
        lgt = allel.GenotypeVector([control])
        rgt = allel.GenotypeVector([test])
        actual = compute_hom_alt_mismatch(control=lgt, test=rgt, recorder=self.recorder)
        self.assertEqual(expected, actual)
        self.assertCountEqual(self.recorder.record.call_args_list,
                              [call(HOMOZIGOUS_ALTERNATE_CATEGORY, MISMATCH_CATEGORY, expected)])


class TestCallset(unittest.TestCase):

    def test_call_set(self):
        zarr_group = generate_gt_data(SAMPLE_A, SOME_CHROMOSOMES_GT)
        under_test = Callset(sample=SAMPLE_A, callset=zarr_group)
        actual = under_test.gt(CHROMOSOME_1)
        self.assertTrue((allel.GenotypeVector(CHROMOSOME_1_GT) == actual).all())
        actual = under_test.gt(CHROMOSOME_2)
        self.assertTrue((allel.GenotypeVector(CHROMOSOME_2_GT) == actual).all())

    @patch('scripts.vgc.zarr.open')
    def test_new_instance(self, mock):
        expected = MagicMock()

        def side_effect(path, mode='z'):
            if mode == 'r' and path == SAMPLE_A:
                return expected
            return None

        mock.side_effect = side_effect

        actual = Callset.new_instance(sample=SAMPLE_A, file_path_format="{sample}")
        self.assertEqual(expected, actual.callset)
        self.assertEqual(SAMPLE_A, actual.sample)


class TestFilteredCallset(unittest.TestCase):

    def setUp(self) -> None:
        self.callset = MagicMock()
        self.callset.gt.side_effect = lambda chromosome: allel.GenotypeVector(SOME_CHROMOSOMES_GT[chromosome])

    def test_call_filtered_set(self):
        site_filters = generate_filter(SOME_CHROMOSOMES_FILTERS)
        under_test = FilteringCallset(site_filters=site_filters, callset=self.callset)
        actual = under_test.gt(CHROMOSOME_1)
        self.assertTrue((allel.GenotypeVector([[1, 1]]) == actual).all())
        actual = under_test.gt(CHROMOSOME_2)
        self.assertTrue((allel.GenotypeVector([[0, 0], [2, 0]]) == actual).all())

    @patch('scripts.vgc.zarr.open')
    def test_new_instance(self, mock):
        site_filter_path = "A_PATH"
        callset = MagicMock()
        expected = MagicMock()

        def side_effect(path, mode='z'):
            if mode == 'r' and path == site_filter_path:
                return expected
            return None

        mock.side_effect = side_effect
        actual = FilteringCallset.new_test_instance(site_filter_path=site_filter_path, callset=callset)
        self.assertEqual(callset, actual.callset)
        self.assertEqual(expected, actual.site_filters)


class TestToFilteredCallset(unittest.TestCase):

    def setUp(self) -> None:
        self.callset_mock = MagicMock()
        self.filtered_callset_mock = MagicMock()
        self.a_file_path_format = "A_FILE_PATH_FORMAT"

    @patch('scripts.vgc.Callset.new_instance')
    @patch('scripts.vgc.FilteringCallset.new_instance')
    def test_to_filtered_callset(self, filtered_callset_mock, callset_mock):
        callset_mock.side_effect = self._callset_mock_side_effect()
        filtered_callset_mock.side_effect = self.filtered_callset_mock_side_effect()
        actual = to_filtered_callset(A_URL, A_FILTER_PATH, self.a_file_path_format, SAMPLE_A)
        self.assertEqual(self.filtered_callset_mock, actual)

    def filtered_callset_mock_side_effect(self):
        return lambda url, path, callset: self.filtered_callset_mock \
            if callset == self.callset_mock and url == A_URL and path == A_FILTER_PATH else None

    def _callset_mock_side_effect(self):
        return lambda sample='', file_path_format='': self.callset_mock \
            if sample == SAMPLE_A and file_path_format == self.a_file_path_format else None


class TestCountParser(unittest.TestCase):

    def setUp(self):
        self.under_test = ArgumentParserBuilder.new_instance(lambda: ErrorRaisingArgumentParser()) \
            .with_count() \
            .build()

    def test_count(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url',
             'url', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_single_sample(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url',
             'url', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_single_chromosome(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', 'sample2', '--chromosomes', 'chromosome1', '--filter-catalog-url',
             'url', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_output_is_defaulted(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--samples',
             'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url',
             'url', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output=None,
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_chromosomes_is_defaulted(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', 'sample2', '--filter-catalog-url',
             'url', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['2L', '2R', '3L', '3R', 'X'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_url_is_defaulted(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-path', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='https://malariagen.github.io/intake/gcs.yml', path='path',
                                            command=Commands.COUNT))

    def test_path_is_defaulted(self):
        args = self.under_test.parse_args(
            ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--samples',
             'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url', 'url'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='url', path='ag3.site_filters_dt_20200416_gamb_colu_arab',
                                            command=Commands.COUNT))

    def test_count_short_options(self):
        args = self.under_test.parse_args(
            ['count', '-c', 'control', '-t', 'test', '-o', 'output', '-s', 'sample1', 'sample2', '--chromosomes',
             'chromosome1', 'chromosome2', '--filter-catalog-url', 'url', '-f', 'path'])
        self.assertEqual(args,
                         argparse.Namespace(control='control', test='test', output='output',
                                            samples=['sample1', 'sample2'], chromosomes=['chromosome1', 'chromosome2'],
                                            url='url', path='path', command=Commands.COUNT))

    def test_fail_if_no_control(self):
        with self.assertRaises(ValueError) as cm:
            self.under_test.parse_args(
                ['count', '--test', 'test', '--output', 'output', '--samples',
                 'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url',
                 'url', '--filter-path', 'path'])
        self.assertEqual(cm.exception.args[0], 'the following arguments are required: --control/-c')

    def test_fail_if_no_test(self):
        with self.assertRaises(ValueError) as cm:
            self.under_test.parse_args(
                ['count', '--control', 'control', '--output', 'output', '--samples',
                 'sample1', 'sample2', '--chromosomes', 'chromosome1', 'chromosome2', '--filter-catalog-url',
                 'url', '--filter-path', 'path'])
        self.assertEqual(cm.exception.args[0], 'the following arguments are required: --test/-t')

    def test_fail_if_no_samples(self):
        with self.assertRaises(ValueError) as cm:
            self.under_test.parse_args(
                ['count', '--control', 'control', '--test', 'test', '--output', 'output', '--chromosomes',
                 'chromosome1', 'chromosome2', '--filter-catalog-url', 'url', '--filter-path', 'path'])
        self.assertEqual(cm.exception.args[0], 'the following arguments are required: --samples/-s')


class TestSummarizeParser(unittest.TestCase):

    def setUp(self):
        self.under_test = ArgumentParserBuilder.new_instance(lambda: ErrorRaisingArgumentParser()) \
            .with_summarize() \
            .build()

    def test_summarize(self):
        args = self.under_test.parse_args(
            ['summarize', '--inputs', 'input1', 'input2', '--output', 'output'])
        self.assertEqual(args,
                         argparse.Namespace(output='output', inputs=['input1', 'input2'], command=Commands.SUMMARIZE))

    def test_single_input(self):
        args = self.under_test.parse_args(
            ['summarize', '--inputs', 'input1', '--output', 'output'])
        self.assertEqual(args,
                         argparse.Namespace(output='output', inputs=['input1'], command=Commands.SUMMARIZE))

    def test_output_is_defaulted(self):
        args = self.under_test.parse_args(
            ['summarize', '--inputs', 'input1', 'input2'])
        self.assertEqual(args,
                         argparse.Namespace(output=None, inputs=['input1', 'input2'], command=Commands.SUMMARIZE))

    def test_summarize_short_options(self):
        args = self.under_test.parse_args(
            ['summarize', '-i', 'input1', 'input2', '-o', 'output'])
        self.assertEqual(args,
                         argparse.Namespace(output='output', inputs=['input1', 'input2'], command=Commands.SUMMARIZE))

    def test_fail_if_no_input(self):
        with self.assertRaises(ValueError) as cm:
            self.under_test.parse_args(
                ['summarize', '--output', 'output'])
        self.assertEqual(cm.exception.args[0], 'the following arguments are required: --inputs/-i')


def mock_function():
    pass


class ErrorRaisingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def generate_gt_data(sample: str, chromosome_gt: dict, chunks=1):
    root = zarr.group()
    sample_group = root.create_group(sample)
    sample_group.create_groups(*chromosome_gt.keys())
    for i in sample_group:
        gt = chromosome_gt[i]
        sites = len(gt)
        gt_data = [[x] for x in gt]
        calldata = sample_group[i].create_group("calldata")
        calldata.create_dataset('GT', shape=(sites, 1, 2), chunks=(chunks, 1, 2), dtype='int8',
                                data=gt_data)

    return root


def generate_filter(chromosome_bool: dict):
    root = zarr.group()
    root.create_groups(*chromosome_bool.keys())
    for i in root:
        bools = chromosome_bool[i]
        variants = root[i].create_group('variants')
        variants.create_dataset('filter_pass', shape=(len(bools),), chunks=(len(bools),), dtype='bool', data=bools)
    return root


SUMMARY_AA0052 = """sample,chromosome,control,test,count
AA0052-C,2L,het,mismatch,0
AA0052-C,2L,hom_alt,mismatch,0
AA0052-C,2L,hom_ref,hom_ref,31905424
AA0052-C,2L,hom_ref,het,11
AA0052-C,2L,hom_ref,hom_alt,0
AA0052-C,2L,hom_ref,missing,35
AA0052-C,2L,het,hom_ref,15
AA0052-C,2L,het,het,317310
AA0052-C,2L,het,hom_alt,16
AA0052-C,2L,het,missing,0
AA0052-C,2L,hom_alt,hom_ref,0
AA0052-C,2L,hom_alt,het,9
AA0052-C,2L,hom_alt,hom_alt,248559
AA0052-C,2L,hom_alt,missing,0
AA0052-C,2L,missing,hom_ref,55
AA0052-C,2L,missing,het,0
AA0052-C,2L,missing,hom_alt,1
AA0052-C,2L,missing,missing,58548
AA0052-C,2R,het,mismatch,0
AA0052-C,2R,hom_alt,mismatch,0
AA0052-C,2R,hom_ref,hom_ref,39944849
AA0052-C,2R,hom_ref,het,12
AA0052-C,2R,hom_ref,hom_alt,0
AA0052-C,2R,hom_ref,missing,27
AA0052-C,2R,het,hom_ref,12
AA0052-C,2R,het,het,371735
AA0052-C,2R,het,hom_alt,3
AA0052-C,2R,het,missing,0
AA0052-C,2R,hom_alt,hom_ref,0
AA0052-C,2R,hom_alt,het,4
AA0052-C,2R,hom_alt,hom_alt,175609
AA0052-C,2R,hom_alt,missing,0
AA0052-C,2R,missing,hom_ref,38
AA0052-C,2R,missing,het,0
AA0052-C,2R,missing,hom_alt,1
AA0052-C,2R,missing,missing,69377
"""

SUMMARY_AA0053 = """sample,chromosome,control,test,count
AA0053-C,2L,het,mismatch,0
AA0053-C,2L,hom_alt,mismatch,0
AA0053-C,2L,hom_ref,hom_ref,25437963
AA0053-C,2L,hom_ref,het,9
AA0053-C,2L,hom_ref,hom_alt,0
AA0053-C,2L,hom_ref,missing,46
AA0053-C,2L,het,hom_ref,7
AA0053-C,2L,het,het,268531
AA0053-C,2L,het,hom_alt,3
AA0053-C,2L,het,missing,0
AA0053-C,2L,hom_alt,hom_ref,0
AA0053-C,2L,hom_alt,het,2
AA0053-C,2L,hom_alt,hom_alt,113941
AA0053-C,2L,hom_alt,missing,0
AA0053-C,2L,missing,hom_ref,18
AA0053-C,2L,missing,het,0
AA0053-C,2L,missing,hom_alt,0
AA0053-C,2L,missing,missing,48865
AA0053-C,2R,het,mismatch,0
AA0053-C,2R,hom_alt,mismatch,0
AA0053-C,2R,hom_ref,hom_ref,32809697
AA0053-C,2R,hom_ref,het,11
AA0053-C,2R,hom_ref,hom_alt,0
AA0053-C,2R,hom_ref,missing,55
AA0053-C,2R,het,hom_ref,19
AA0053-C,2R,het,het,365106
AA0053-C,2R,het,hom_alt,8
AA0053-C,2R,het,missing,4
AA0053-C,2R,hom_alt,hom_ref,0
AA0053-C,2R,hom_alt,het,1
AA0053-C,2R,hom_alt,hom_alt,151264
AA0053-C,2R,hom_alt,missing,0
AA0053-C,2R,missing,hom_ref,26
AA0053-C,2R,missing,het,0
AA0053-C,2R,missing,hom_alt,1
AA0053-C,2R,missing,missing,64670
"""


class TestSummarization(unittest.TestCase):

    def setUp(self) -> None:
        self.df_AA0052 = Summarizer.to_dataframe(StringIO(SUMMARY_AA0052))
        self.df_AA0053 = Summarizer.to_dataframe(StringIO(SUMMARY_AA0053))
        self.under_test = Summarizer([self.df_AA0052, self.df_AA0053])

    def test_total_concordance(self):
        actual = self.under_test.compute_total_concordance()
        self.assertEqual('Total', actual.iloc[0]['description'])
        self.assertEqual(132351448, actual.iloc[0]['match'])
        self.assertAlmostEqual(132351897, actual.iloc[0]['total'], delta=0.000001)
        self.assertAlmostEqual(99.999661, actual.iloc[0]['concordance'], delta=0.000001)
        self.assertEqual('Total no missing no homref_homref', actual.iloc[1]['description'])
        self.assertEqual(2012197, actual.iloc[1]['total'])
        self.assertEqual(2012055, actual.iloc[1]['match'])
        self.assertAlmostEqual(99.992943, actual.iloc[1]['concordance'], delta=0.000001)

    def test_total_concordance_per_category(self):
        actual = self.under_test.compute_concordance_per_category()
        self.assertEqual(18, len(actual))
        self._assert_concordance_per_category_row(actual.iloc[0], 'het', 'het', 1322769, 1322682, 99.993423)
        self._assert_concordance_per_category_row(actual.iloc[1], 'het', 'hom_alt', 1322769, 30, 0.002268)
        self._assert_concordance_per_category_row(actual.iloc[2], 'het', 'hom_ref', 1322769, 53, 0.004007)
        self._assert_concordance_per_category_row(actual.iloc[3], 'het', 'mismatch', 1322769, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[4], 'het', 'missing', 1322769, 4, 0.000302)
        self._assert_concordance_per_category_row(actual.iloc[5], 'hom_alt', 'het', 689389, 16, 0.002321)
        self._assert_concordance_per_category_row(actual.iloc[6], 'hom_alt', 'hom_alt', 689389, 689373, 99.997679)
        self._assert_concordance_per_category_row(actual.iloc[7], 'hom_alt', 'hom_ref', 689389, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[8], 'hom_alt', 'mismatch', 689389, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[9], 'hom_alt', 'missing', 689389, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[10], 'hom_ref', 'het', 130098139, 43, 0.000033)
        self._assert_concordance_per_category_row(actual.iloc[11], 'hom_ref', 'hom_alt', 130098139, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[12], 'hom_ref', 'hom_ref', 130098139, 130097933,
                                                  99.999842)
        self._assert_concordance_per_category_row(actual.iloc[13], 'hom_ref', 'missing', 130098139, 163, 0.000125)
        self._assert_concordance_per_category_row(actual.iloc[14], 'missing', 'het', 241600, 0, 0)
        self._assert_concordance_per_category_row(actual.iloc[15], 'missing', 'hom_alt', 241600, 3, 0.001242)
        self._assert_concordance_per_category_row(actual.iloc[16], 'missing', 'hom_ref', 241600, 137, 0.056705)
        self._assert_concordance_per_category_row(actual.iloc[17], 'missing', 'missing', 241600, 241460, 99.942053)

    def test_compute_concordance_per_sample_category(self):
        actual = self.under_test.compute_concordance_per_sample_category()
        self.assertEqual(36, len(actual))
        self._assert_sample_concordance(actual.iloc[0], 'AA0052-C', 'het', 'het', 689091, 689045, 99.993325)
        self._assert_sample_concordance(actual.iloc[1], 'AA0052-C', 'het', 'hom_alt', 689091, 19, 0.002757)
        self._assert_sample_concordance(actual.iloc[2], 'AA0052-C', 'het', 'hom_ref', 689091, 27, 0.003918)
        self._assert_sample_concordance(actual.iloc[3], 'AA0052-C', 'het', 'mismatch', 689091, 0, 0)
        self._assert_sample_concordance(actual.iloc[4], 'AA0052-C', 'het', 'missing', 689091, 0, 0)
        self._assert_sample_concordance(actual.iloc[5], 'AA0052-C', 'hom_alt', 'het', 424181, 13, 0.003065)
        self._assert_sample_concordance(actual.iloc[6], 'AA0052-C', 'hom_alt', 'hom_alt', 424181, 424168, 99.996935)
        self._assert_sample_concordance(actual.iloc[7], 'AA0052-C', 'hom_alt', 'hom_ref', 424181, 0, 0)
        self._assert_sample_concordance(actual.iloc[8], 'AA0052-C', 'hom_alt', 'mismatch', 424181, 0, 0)
        self._assert_sample_concordance(actual.iloc[9], 'AA0052-C', 'hom_alt', 'missing', 424181, 0, 0)
        self._assert_sample_concordance(actual.iloc[10], 'AA0052-C', 'hom_ref', 'het', 71850358, 23, 0.000032)
        self._assert_sample_concordance(actual.iloc[11], 'AA0052-C', 'hom_ref', 'hom_alt', 71850358, 0, 0)
        self._assert_sample_concordance(actual.iloc[12], 'AA0052-C', 'hom_ref', 'hom_ref', 71850358, 71850273,
                                        99.999882)
        self._assert_sample_concordance(actual.iloc[13], 'AA0052-C', 'hom_ref', 'missing', 71850358, 62, 0.000086)
        self._assert_sample_concordance(actual.iloc[14], 'AA0052-C', 'missing', 'het', 128020, 0, 0)
        self._assert_sample_concordance(actual.iloc[15], 'AA0052-C', 'missing', 'hom_alt', 128020, 2, 0.001562)
        self._assert_sample_concordance(actual.iloc[16], 'AA0052-C', 'missing', 'hom_ref', 128020, 93, 0.072645)
        self._assert_sample_concordance(actual.iloc[17], 'AA0052-C', 'missing', 'missing', 128020, 127925, 99.925793)
        self._assert_sample_concordance(actual.iloc[18], 'AA0053-C', 'het', 'het', 633678, 633637, 99.99353)
        self._assert_sample_concordance(actual.iloc[19], 'AA0053-C', 'het', 'hom_alt', 633678, 11, 0.001736)
        self._assert_sample_concordance(actual.iloc[20], 'AA0053-C', 'het', 'hom_ref', 633678, 26, 0.004103)
        self._assert_sample_concordance(actual.iloc[21], 'AA0053-C', 'het', 'mismatch', 633678, 0, 0)
        self._assert_sample_concordance(actual.iloc[22], 'AA0053-C', 'het', 'missing', 633678, 4, 0.000631)
        self._assert_sample_concordance(actual.iloc[23], 'AA0053-C', 'hom_alt', 'het', 265208, 3, 0.001131)
        self._assert_sample_concordance(actual.iloc[24], 'AA0053-C', 'hom_alt', 'hom_alt', 265208, 265205, 99.998869)
        self._assert_sample_concordance(actual.iloc[25], 'AA0053-C', 'hom_alt', 'hom_ref', 265208, 0, 0)
        self._assert_sample_concordance(actual.iloc[26], 'AA0053-C', 'hom_alt', 'mismatch', 265208, 0, 0)
        self._assert_sample_concordance(actual.iloc[27], 'AA0053-C', 'hom_alt', 'missing', 265208, 0, 0)
        self._assert_sample_concordance(actual.iloc[28], 'AA0053-C', 'hom_ref', 'het', 58247781, 20, 0.000034)
        self._assert_sample_concordance(actual.iloc[29], 'AA0053-C', 'hom_ref', 'hom_alt', 58247781, 0, 0)
        self._assert_sample_concordance(actual.iloc[30], 'AA0053-C', 'hom_ref', 'hom_ref', 58247781, 58247660,
                                        99.999792)
        self._assert_sample_concordance(actual.iloc[31], 'AA0053-C', 'hom_ref', 'missing', 58247781, 101, 0.000173)
        self._assert_sample_concordance(actual.iloc[32], 'AA0053-C', 'missing', 'het', 113580, 0, 0)
        self._assert_sample_concordance(actual.iloc[33], 'AA0053-C', 'missing', 'hom_alt', 113580, 1, 0.00088)
        self._assert_sample_concordance(actual.iloc[34], 'AA0053-C', 'missing', 'hom_ref', 113580, 44, 0.038739)
        self._assert_sample_concordance(actual.iloc[35], 'AA0053-C', 'missing', 'missing', 113580, 113535, 99.96038)

    def test_compute_concordance_per_chromosome_category(self):
        actual = self.under_test.compute_concordance_per_chromosome_category()
        self.assertEqual(36, len(actual))
        self._assert_chromosome_concordance(actual.iloc[0], '2L', 'het', 'het', 585882, 585841, 99.993002)
        self._assert_chromosome_concordance(actual.iloc[1], '2L', 'het', 'hom_alt', 585882, 19, 0.003243)
        self._assert_chromosome_concordance(actual.iloc[2], '2L', 'het', 'hom_ref', 585882, 22, 0.003755)
        self._assert_chromosome_concordance(actual.iloc[3], '2L', 'het', 'mismatch', 585882, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[4], '2L', 'het', 'missing', 585882, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[5], '2L', 'hom_alt', 'het', 362511, 11, 0.003034)
        self._assert_chromosome_concordance(actual.iloc[6], '2L', 'hom_alt', 'hom_alt', 362511, 362500, 99.996966)
        self._assert_chromosome_concordance(actual.iloc[7], '2L', 'hom_alt', 'hom_ref', 362511, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[8], '2L', 'hom_alt', 'mismatch', 362511, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[9], '2L', 'hom_alt', 'missing', 362511, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[10], '2L', 'hom_ref', 'het', 57343488, 20, 0.000035)
        self._assert_chromosome_concordance(actual.iloc[11], '2L', 'hom_ref', 'hom_alt', 57343488, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[12], '2L', 'hom_ref', 'hom_ref', 57343488, 57343387, 99.999824)
        self._assert_chromosome_concordance(actual.iloc[13], '2L', 'hom_ref', 'missing', 57343488, 81, 0.000141)
        self._assert_chromosome_concordance(actual.iloc[14], '2L', 'missing', 'het', 107487, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[15], '2L', 'missing', 'hom_alt', 107487, 1, 0.00093)
        self._assert_chromosome_concordance(actual.iloc[16], '2L', 'missing', 'hom_ref', 107487, 73, 0.067915)
        self._assert_chromosome_concordance(actual.iloc[17], '2L', 'missing', 'missing', 107487, 107413, 99.931154)

        self._assert_chromosome_concordance(actual.iloc[18], '2R', 'het', 'het', 736887, 736841, 99.993758)
        self._assert_chromosome_concordance(actual.iloc[19], '2R', 'het', 'hom_alt', 736887, 11, 0.001493)
        self._assert_chromosome_concordance(actual.iloc[20], '2R', 'het', 'hom_ref', 736887, 31, 0.004207)
        self._assert_chromosome_concordance(actual.iloc[21], '2R', 'het', 'mismatch', 736887, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[22], '2R', 'het', 'missing', 736887, 4, 0.000543)
        self._assert_chromosome_concordance(actual.iloc[23], '2R', 'hom_alt', 'het', 326878, 5, 0.00153)
        self._assert_chromosome_concordance(actual.iloc[24], '2R', 'hom_alt', 'hom_alt', 326878, 326873, 99.99847)
        self._assert_chromosome_concordance(actual.iloc[25], '2R', 'hom_alt', 'hom_ref', 326878, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[26], '2R', 'hom_alt', 'mismatch', 326878, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[27], '2R', 'hom_alt', 'missing', 326878, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[28], '2R', 'hom_ref', 'het', 72754651, 23, 0.000032)
        self._assert_chromosome_concordance(actual.iloc[29], '2R', 'hom_ref', 'hom_alt', 72754651, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[30], '2R', 'hom_ref', 'hom_ref', 72754651, 72754546, 99.999856)
        self._assert_chromosome_concordance(actual.iloc[31], '2R', 'hom_ref', 'missing', 72754651, 82, 0.000113)
        self._assert_chromosome_concordance(actual.iloc[32], '2R', 'missing', 'het', 134113, 0, 0)
        self._assert_chromosome_concordance(actual.iloc[33], '2R', 'missing', 'hom_alt', 134113, 2, 0.001491)
        self._assert_chromosome_concordance(actual.iloc[34], '2R', 'missing', 'hom_ref', 134113, 64, 0.047721)
        self._assert_chromosome_concordance(actual.iloc[35], '2R', 'missing', 'missing', 134113, 134047, 99.950788)

    def _assert_concordance_per_category_row(self, local, control, test, total, match, concordance):
        self.assertEqual(control, local['control'])
        self.assertEqual(test, local['test'])
        self.assertEqual(match, local['match'])
        self.assertAlmostEqual(total, local['total'], delta=0.000001)
        self.assertAlmostEqual(concordance, local['concordance'], delta=0.000001)

    def _assert_sample_concordance(self, local, sample, control, test, total, match, concordance):
        self.assertEqual(sample, local['sample'])
        self.assertEqual(control, local['control'])
        self.assertEqual(test, local['test'])
        self.assertEqual(match, local['match'])
        self.assertAlmostEqual(total, local['total'], delta=0.000001)
        self.assertAlmostEqual(concordance, local['concordance'], delta=0.000001)

    def _assert_chromosome_concordance(self, local, chromosome, control, test, total, match, concordance):
        self.assertEqual(chromosome, local['chromosome'])
        self.assertEqual(control, local['control'])
        self.assertEqual(test, local['test'])
        self.assertEqual(match, local['match'])
        self.assertAlmostEqual(total, local['total'], delta=0.000001)
        self.assertAlmostEqual(concordance, local['concordance'], delta=0.000001)
