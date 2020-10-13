#!/usr/bin/env python3

# runtime: ['intake', 'zarr', 'numpy', 'pandas', 'scikit-allel', 'requests', 'aiohttp', 'gcsfs']
# tests: ['nose >= 1.3', 'ddt']
# Example use:
# count: vgc.py count -c '/path/to/zarrs/{sample}/{sample}.genotypes.zarr.zip' -t '/path/to/zarrs/{sample}.zarr.zip'  -s AA0052-C AB0252-C AC0010-C AJ0037-C AN0131-C AN0280-Cx AN0326-C AR0078-C AZ0156-C BL0358-C -o /some/dir/results.csv 
# summary: vgc.py summarize -i results.*.csv -o /some/dir/10_samples

import argparse
import sys
from enum import Enum
from functools import partial
from sys import stdout

import allel
import intake
import numpy as np
import pandas
import zarr

# ZARR paths
FILTER_PASS_PATH_FORMAT = '{chromosome}/variants/filter_pass'
GT_PATH_FORMAT = '{sample}/{chromosome}/calldata/GT'

# Genotypes categories
MISMATCH_CATEGORY = 'mismatch'
MISSING_CATEGORY = 'missing'
HOMOZIGOUS_ALTERNATE_CATEGORY = 'hom_alt'
HETEROZIGOUS_CATEGORY = 'het'
HOMOZEGOUS_REFERENCE_CATEGORY = 'hom_ref'
NO_MISMATCH_CATEGORIES = [HOMOZEGOUS_REFERENCE_CATEGORY, HETEROZIGOUS_CATEGORY, HOMOZIGOUS_ALTERNATE_CATEGORY,
                          MISSING_CATEGORY]

# Genotype count output headers
SAMPLE_HEADER = 'sample'
CHROMOSOME_HEADER = 'chromosome'
CONTROL_HEADER = 'control'
TEST_HEADER = 'test'
COUNT = 'count'

"""
A class holding the statistics (ie count) for concordance calculation
"""


class ChromosomeConcordanceResultRecorder:
    """
    Internal class allowing to record statistics for a chromosome
    """

    def __init__(self, parent_result, sample: str, chromosome: str):
        """

        Parameters
        ----------
        parent_result: ConcordanceResult
            The parent ConcordanceResult
        sample: str
            The sample to records statistics on
        chromosome: str
            The chromosone to record statistics on
        """
        self.results = parent_result
        self.sample = sample
        self.chromosome = chromosome

    def record(self, control: str, test: str, count: int):
        """
        Records statistics
        Parameters
        ----------
        control: str
            the genotype category for the control side
        test: str
            the genotype category for the test side
        count
            the number of genotypes

        """
        series = pandas.Series([self.sample, self.chromosome, control, test, count],
                               index=[SAMPLE_HEADER, CHROMOSOME_HEADER, CONTROL_HEADER, TEST_HEADER, COUNT])
        self.results.df = self.results.df.append(series, ignore_index=True)


class ConcordanceResult:
    def __init__(self):
        self.df = pandas.DataFrame({
            SAMPLE_HEADER: pandas.Series([], dtype='str'),
            CHROMOSOME_HEADER: pandas.Series([], dtype='str'),
            CONTROL_HEADER: pandas.Series([], dtype='str'),
            TEST_HEADER: pandas.Series([], dtype='str'),
            COUNT: pandas.Series([], dtype='int'),
        })

    def record_chromosome_statistics(self, sample: str, chromosome: str) -> ChromosomeConcordanceResultRecorder:
        """
        Build a recorder of statistics for a sample/chromosoe pair
        Parameters
        ----------
        sample: str
            The sample to record statistics on
        chromosome: str
            The chromosome to record statistics on

        Returns: ChromosomeConcordanceResultRecorder
            A corresponding ChromosomeConcordanceResultRecorder
        -------
        """
        return ChromosomeConcordanceResultRecorder(self, sample, chromosome)

    def print(self, path_or_stream):
        """
        Prints the recorded statistics as csv
        Parameters
        ----------
        path_or_stream:
            Either a path to a file or a stream (stdout for example)


        """
        self.df.to_csv(path_or_buf=path_or_stream, index=False)


class Callset:
    """
    Wrapper class around a genotype callset stored in zarr
    """

    @staticmethod
    def new_instance(*, sample: str, file_path_format: str):
        """

        Parameters
        ----------
        sample: str
            The sample
        file_path_format: str
            The format of the path to the zarr file

        Returns: Callset
            The corresponding Callset
        -------

        """
        callset = zarr.open(file_path_format.format(sample=sample), mode='r')
        return Callset(sample=sample, callset=callset)

    def __init__(self, *, sample, callset):
        """
        Constructor, used for testing, use factory method
        Parameters
        ----------
        sample: str
            the sample
        callset:
            the root zarr group
        """
        self.sample = sample
        self.callset = callset

    def gt(self, chromosome: str) -> allel.GenotypeVector:
        """
        Builds a genotype vector for a chroosome
        Parameters
        ----------
        chromosome: str
            The chromosome
        Returns: allel.GenotypeVector
            A genotype vector for the chromosome
        -------
        """
        gt = self.callset[GT_PATH_FORMAT.format(sample=self.sample, chromosome=chromosome)]
        return allel.GenotypeVector(gt[:, 0])


class FilteringCallset:
    """
    Wrapper class around callset that applies a pass site filter
    """

    @staticmethod
    def new_test_instance(*, site_filter_path: str, callset: Callset):
        """
        A factory method used to build filter pass from local file
        Parameters
        ----------
        site_filter_path: str
            the path to the file containing the site filter
        callset: Callset
            the genotype callset to filter

        Returns: FilteringCallset
            The corresponding FilteringCallset
        -------
        """
        site_filters = zarr.open(site_filter_path, mode='r')
        return FilteringCallset(site_filters=site_filters, callset=callset)

    @staticmethod
    def new_instance(*, url: str, path: str, callset: Callset):
        """
        Factory method of FilteringCallset from an online catalog
        Parameters
        ----------
        url: str
            the url to the catalog
        path: str
            the path in the catalog to the site filter
        callset: Callset
            the genotype callset to filter

        Returns: FilteringCallset
            The corresponding FilteringCallset
        -------

        """
        cat = intake.open_catalog(url)
        site_filters = cat[path].to_zarr()
        return FilteringCallset(site_filters=site_filters, callset=callset)

    def __init__(self, *, site_filters, callset):
        """
        Constructor
        Parameters
        ----------
        site_filters: zarr group
            The site filter
        callset: Callset
            The genotype callset to filter
        """
        self.site_filters = site_filters
        self.callset = callset

    def _pass_filter_for_chromosome(self, chromosome: str):
        """
        Return pass filter for a chromosome
        Parameters
        ----------
        chromosome: str
            The chromosome to filter

        Returns
            The pass filter
        -------

        """
        path = FILTER_PASS_PATH_FORMAT.format(chromosome=chromosome)
        return self.site_filters[path][:] if path in self.site_filters else None

    def gt(self, chromosome: str) -> allel.GenotypeVector:
        """

        Parameters
        ----------
        chromosome: str
          The chromosome

        Returns: allel.GenotypeVector
            The filtered genotype vector
        -------

        """
        filter_pass = self._pass_filter_for_chromosome(chromosome)
        unfiltered_gt = self.callset.gt(chromosome)
        return unfiltered_gt.compress(filter_pass, axis=0) if filter_pass is not None else unfiltered_gt


def to_filtered_callset(url: str, filter_path: str, file_path_format: str, sample: str):
    """
    A factory method of filtered callset
    Parameters
    ----------
    url: str
        The url to the catalog
    filter_path: str
        The path to the filter in the catalog
    file_path_format: str
        The path to the zarr, as a format based on sample
    sample: str
        The sample

    Returns: FilteringCallset
        The corresponding FilteringCallset
    -------

    """
    callset = Callset.new_instance(sample=sample, file_path_format=file_path_format)
    return FilteringCallset.new_instance(url=url, path=filter_path, callset=callset)


def classify_chromosome(*, control: allel.GenotypeVector, test: allel.GenotypeVector,
                        recorder: ChromosomeConcordanceResultRecorder):
    """
    Classifies a chromosomes
    Parameters
    ----------
    control: allel.GenotypeVector
        the genotype vector for control side
    test: allel.GenotypeVector
        the genotype vector for test side
    recorder: ChromosomeConcordanceResultRecorder
        the recorder of statistics

    -------

    """
    het_mismatch = compute_het_mismatch(control, test, recorder)
    hom_alt_mismatch = compute_hom_alt_mismatch(control, test, recorder)
    _classify_non_mismatch(het_mismatch, hom_alt_mismatch, control, test, recorder)


def _classify_non_mismatch(het_mismatch, hom_alt_mismatch, control, test, recorder):
    """
    Classifies a chromosomes considering het_het and hom_alt_hom_alt mismatches

    Parameters
    ----------
    het_mismatch:
        number of het-het mismatches
    hom_alt_mismatch
        number of hom-alt_hom-alt mismatches
    control: allel.GenotypeVector
        the genotype vector for control side
    test: allel.GenotypeVector
        the genotype vector for test side
    recorder: ChromosomeConcordanceResultRecorder
        the recorder of statistics

    """
    mismatches = np.diag([0, het_mismatch, hom_alt_mismatch, 0])
    x1 = [control.is_hom_ref(), control.is_het(), control.is_hom_alt(), control.is_missing()]
    x2 = [test.is_hom_ref(), test.is_het(), test.is_hom_alt(), test.is_missing()]
    for i, t1 in enumerate(x1):
        for j, t2 in enumerate(x2):
            recorder.record(NO_MISMATCH_CATEGORIES[i], NO_MISMATCH_CATEGORIES[j],
                            np.count_nonzero(t1 & t2) - mismatches[i][j])


def compute_hom_alt_mismatch(control, test, recorder):
    """
    Compute hom_alt mismatches, for example [1,1] vs [2,2]
    Parameters
    ----------
    control: allel.GenotypeVector
        the genotype vector for control side
    test: allel.GenotypeVector
        the genotype vector for test side
    recorder: ChromosomeConcordanceResultRecorder
        the recorder of statistics

    Returns:
        the number of hom_alt mismatches
    -------

    """
    hom_alt_mismatch = _mismatch_chromosome(lambda gtf: gtf.is_hom_alt(), control, test)
    recorder.record(HOMOZIGOUS_ALTERNATE_CATEGORY, MISMATCH_CATEGORY, hom_alt_mismatch)
    return hom_alt_mismatch


def compute_het_mismatch(control, test, recorder):
    """
    Compute het mismatches, for example [0,1] vs [0,2]
    Parameters
    ----------
    control: allel.GenotypeVector
        the genotype vector for control side
    test: allel.GenotypeVector
        the genotype vector for test side
    recorder: ChromosomeConcordanceResultRecorder
        the recorder of statistics

    Returns:
        the number of hety mismatches
    -------

    """
    het_mismatch = _mismatch_chromosome(lambda gtf: gtf.is_het(), control, test)
    recorder.record(HETEROZIGOUS_CATEGORY, MISMATCH_CATEGORY, het_mismatch)
    return het_mismatch


def _mismatch_chromosome(filter_function, control, test):
    """
    Generic mismatch recorder
    Parameters
    ----------
    filter_function
      the function to apply to the genotype vector to select relevant genotype
    control: allel.GenotypeVector
        the genotype vector for control side
    test: allel.GenotypeVector
        the genotype vector for test side

    Returns
        The number of mismatches
    -------

    """
    combined_filter = filter_function(control) & filter_function(test)
    filtered_control = control.compress(combined_filter, axis=0)
    filtered_test = test.compress(combined_filter, axis=0)
    return len(filtered_control) - np.sum((filtered_control == filtered_test).all(1))


def classify_sample(sample, chromosomes, control, test, results):
    """
    Record statistics and categories for a sample
    Parameters
    ----------
    sample: str
        the sample
    chromosomes: list(str)
        the chromosomes to classify
    control:
        the genotype vector of the control side
    test
        the genotype vector of the test side
    results
        the results to record the stats into
    -------

    """
    for chromosome in chromosomes:
        print("Processing sample {}, chromosome {}".format(sample, chromosome), file=sys.stderr)
        classify_chromosome(control=control.gt(chromosome), test=test.gt(chromosome),
                            recorder=results.record_chromosome_statistics(sample, chromosome))


def classify(samples, chromosomes, control_callset_factory, test_callset_factory, results):
    """
    Classify genotypes in categories across samples/chromosomes
    Parameters
    ----------
    samples: list(str)
        the samples to classify
    chromosomes: list(str)
        the chromosomes to classify
    control_callset_factory:
        the factory of genotype callset for the control side
    test_callset_factory
        the factory of genotype callset for the test side
    results
        the results to record the stats into
    """
    for sample in samples:
        control_callset = control_callset_factory(sample)
        test_callset = test_callset_factory(sample)
        classify_sample(sample, chromosomes, control_callset, test_callset, results)


class Commands(Enum):
    """
    Enum of commands that the script can take
    """
    COUNT = 1
    SUMMARIZE = 2


class ArgumentParserBuilder:
    """ Helper class to build and test the argument parsing"""

    @staticmethod
    def new_instance(factory=lambda: argparse.ArgumentParser()):
        return ArgumentParserBuilder(factory())

    def __init__(self, parser):
        self.parser = parser
        self.subparsers = self.parser.add_subparsers(help='sub-command help')

    def with_count(self):
        count = self.subparsers.add_parser("count", help='Count the gt by categories')
        count.add_argument('--control', '-c', dest='control', required=True,
                           help='A format that points to the control zarr file, example /control/{sample}.zar')
        count.add_argument('--test', '-t', dest='test', required=True,
                           help='A format that points to the test zarr file, example /test/{sample}.zar')
        count.add_argument('--samples', '-s', nargs='+', dest='samples', required=True,
                           help='The list of samples')
        count.add_argument('--output', '-o', dest='output', required=False, default=None,
                           help='The output file, by default prints to stdout')
        count.add_argument('--chromosomes', nargs='+', dest='chromosomes', required=False,
                           default=['2L', '2R', '3L', '3R', 'X'], help='The list of chromosomes')
        count.add_argument('--filter-catalog-url', dest='url', required=False,
                           default='https://malariagen.github.io/intake/gcs.yml',
                           help='The url of the catalog containing the filters')
        count.add_argument('--filter-path', '-f', dest='path', required=False,
                           default='ag3.site_filters_dt_20200416_gamb_colu_arab',
                           help='The path in the catalog to the filter required, '
                                'ie ag3.site_filters_dt_20200416_gamb_colu_arab')
        count.set_defaults(command=Commands.COUNT)
        return self

    def with_summarize(self):
        count = self.subparsers.add_parser("summarize", help='Summarize results of several count output')
        count.add_argument('--inputs', '-i', nargs='+', dest='inputs', required=True,
                           help='The list of files to summarize')
        count.add_argument('--output', '-o', dest='output', required=False, default=None,
                           help='The output file, by default prints to stdout')
        count.set_defaults(command=Commands.SUMMARIZE)
        return self

    def build(self):
        return self.parser


class Summarizer:
    """Summarizer manipulates the statistics collected and prints concordance based on various groupings"""

    @staticmethod
    def to_dataframe(stream):
        """
        Helper function reading a csv into a pandas dataframes
        Parameters
        ----------
        stream:
            the stream or filename

        Returns
            the pandas dataframe
        -------

        """
        return pandas.read_csv(stream, dtype={SAMPLE_HEADER: 'str', CHROMOSOME_HEADER: 'str', CONTROL_HEADER: 'str',
                                              TEST_HEADER: 'str', COUNT: 'int64'})

    def __init__(self, dataframes):
        """
        Constructors
        Parameters
        ----------
        dataframes:
            list of pandas dataframe to be merged
        """
        self.dataframe = pandas.concat(dataframes, ignore_index=True, sort=False)
        self.dataframe['mismatch'] = self.dataframe[CONTROL_HEADER] != self.dataframe[TEST_HEADER]
        self.dataframe['category'] = self.dataframe[[CONTROL_HEADER, TEST_HEADER]].agg('_'.join, axis=1)

    def compute_total_concordance(self):
        """
        Computes total concordance
        Returns
        -------
            A dataframe containing total, match and concordance
        """
        total = self._compute_total_concordance()
        no_missing_homref_homref = self._compute_concordance_no_missing_no_homref_homref()
        result = pandas.concat([total, no_missing_homref_homref], ignore_index=True, sort=False)
        result['description'] = ['Total', 'Total no missing no homref_homref']
        return result.reindex(columns=['description', 'total', 'match', 'concordance'])

    def compute_concordance_per_category(self):
        """
        Computes concordance per category
        Returns
        -------
            A corresponding dataframe containing total, match and concordance per category
        """
        return self._concordance_per_category(self.dataframe)

    def compute_concordance_per_sample_category(self):
        """
        Computes concordance per sample and category
        Returns
        -------
            A corresponding dataframe containing total, match and concordance per sample and category
        """
        return self.dataframe.groupby([TEST_HEADER, CONTROL_HEADER, SAMPLE_HEADER]) \
            .agg(count=('count', sum)) \
            .groupby(level=2) \
            .apply(lambda x: self._concordance_per_category(x)) \
            .reset_index() \
            .drop(columns='level_1') \
            .sort_values([SAMPLE_HEADER, CONTROL_HEADER, TEST_HEADER])

    def compute_concordance_per_chromosome_category(self):
        """
        Computes concordance per chromosome and category
        Returns
        -------
            A corresponding dataframe containing total, match and concordance per chromosome and category
        """
        return self.dataframe.groupby([TEST_HEADER, CONTROL_HEADER, CHROMOSOME_HEADER]) \
            .agg(count=('count', sum)) \
            .groupby(level=2) \
            .apply(lambda x: self._concordance_per_category(x)) \
            .reset_index() \
            .drop(columns='level_1') \
            .sort_values([CHROMOSOME_HEADER, CONTROL_HEADER, TEST_HEADER])

    def print(self, *, category, chromosome, count, sample, total):
        """
        Writes summaries to corresponding streams
        Parameters
        ----------
        category:
            category stream for a per category concordance
        chromosome
            chromosome stream for per chromosome per category concordance
        count
            count stream, for a single merged file of count
        sample
            sample stream for per sample per category concordance
        total
            total stream for total concordance
        """
        self.print_dataframe(count, self.dataframe)
        self.print_dataframe(total, self.compute_total_concordance())
        self.print_dataframe(category, self.compute_concordance_per_category())
        self.print_dataframe(sample, self.compute_concordance_per_sample_category())
        self.print_dataframe(chromosome, self.compute_concordance_per_chromosome_category())

    def _compute_total_concordance(self):
        grouped = self.dataframe.groupby('mismatch') \
            .agg(total=('count', sum), match=('count', sum), concordance=('count', sum)) \
            .apply(lambda x:
                   100 * x / x.sum() if x.name == 'concordance'
                   else x * x.sum() / x if x.name == 'total'
                   else x) \
            .reset_index() \
            .astype({"total": int})
        return grouped[grouped.mismatch == False].drop(columns='mismatch')

    def _compute_concordance_no_missing_no_homref_homref(self):
        actual = self.compute_concordance_per_category()
        actual['category'] = actual[[CONTROL_HEADER, TEST_HEADER]].agg('_'.join, axis=1)
        actual = actual.reset_index()
        match = self._sum_matches(actual, actual.category.isin(['het_het', 'hom_alt_hom_alt']))
        total = self._sum_matches(actual, (actual.test != 'missing') & (actual.control != 'missing') & (
                actual.category != 'hom_ref_hom_ref'))
        return pandas.DataFrame({"total": total, "match": match, "concordance": 100 * match / total}, index=[0])

    @staticmethod
    def _sum_matches(actual, condition):
        return actual[condition].agg({'match': ['sum']})['match'][0]

    @staticmethod
    def _concordance_per_category(dataframe):
        return dataframe.groupby([TEST_HEADER, CONTROL_HEADER]) \
            .agg(count=('count', sum)) \
            .groupby(level=1) \
            .apply(lambda x: Summarizer._concordance_grouping(x)) \
            .reset_index() \
            .drop(columns='level_1') \
            .sort_values([CONTROL_HEADER, TEST_HEADER])

    @staticmethod
    def _concordance_grouping(dataframe):
        return dataframe.groupby(TEST_HEADER) \
            .agg(total=('count', sum), match=('count', sum), concordance=('count', sum)) \
            .apply(lambda x:
                   100 * x / x.sum() if x.name == 'concordance'
                   else (x + 1) * x.sum() / (x + 1) if x.name == 'total'
                   else x) \
            .reset_index() \
            .astype({"total": int})

    @staticmethod
    def print_dataframe(output, dataframe):
        print(dataframe.to_string(index=False), file=output, flush=True)


def main():
    parser = ArgumentParserBuilder.new_instance().with_count().with_summarize().build()
    arguments = dict(vars(parser.parse_args()))
    if arguments.get('command', None) == Commands.COUNT:
        control_callset_factory = partial(to_filtered_callset, arguments['url'], arguments['path'],
                                          arguments['control'])
        test_callset_factory = partial(to_filtered_callset, arguments['url'], arguments['path'], arguments['test'])
        results = ConcordanceResult()
        classify(arguments['samples'], arguments['chromosomes'], control_callset_factory, test_callset_factory, results)
        results.print(stdout if arguments['output'] is None else arguments['output'])
    elif arguments.get('command', None) == Commands.SUMMARIZE:
        dataframes = [Summarizer.to_dataframe(i) for i in arguments['inputs']]
        summarizer = Summarizer(dataframes)
        summarizer.print(**to_streams(arguments['output']))


def to_streams(out):
    """
    Work out the streams to use for summarize output based on the output argument
    Parameters
    ----------
    out: str
        the output argument
    Returns
    -------
        A dictionary of streams to use with keys: category, chromosome, count, sample, total
    """
    def to_stream(name_format, file):
        return open(name_format.format(file), 'w')

    if out is None:
        return {'category': stdout, 'chromosome': stdout, 'count': stdout, 'sample': stdout, 'total': stdout}
    return {'category': to_stream("{}.categories.txt", out),
            'chromosome': to_stream("{}.chromosomes.txt", out),
            'count': to_stream("{}.count.csv", out),
            'sample': to_stream("{}.samples.txt", out),
            'total': to_stream("{}.totals.txt", out)
            }


if __name__ == '__main__':
    main()
