import itertools
import pathlib
import subprocess
import tempfile
import unittest
from functools import partial

import numpy
import zarr
from numpy import int8, vectorize, fromfunction


class TestScript(unittest.TestCase):

    def setUp(self) -> None:
        self.project_dir = str(pathlib.Path(__file__).parent.parent.absolute())

    def test_count(self):
        with tempfile.TemporaryDirectory() as directory:
            self._generate_data(directory)
            generated_output_file = directory + '/S1.count.csv'
            self._run_process(self._generate_count_command_line(directory, generated_output_file))
            expected_output_file = self.project_dir + "/fixture/expected.S1.count.csv"
            self._compare_output_files(expected_output_file, generated_output_file)

    def test_summarize(self):
        with tempfile.TemporaryDirectory() as directory:
            self._run_process(self._generate_summarize_command_line(directory))
            for file in ['S1.output.categories.txt', 'S1.output.chromosomes.txt', 'S1.output.count.csv',
                         'S1.output.samples.txt', 'S1.output.totals.txt']:
                self._compare_output_files(directory + '/' + file, self.project_dir + "/fixture/" + file)

    def _generate_data(self, directory):
        value_generator = GTValueGenerator()
        generate_zarr(directory + '/S1.genotypes.zarr', value_generator.generate_left)
        generate_zarr(directory + '/S1.zarr', value_generator.generate_right)

    def _generate_count_command_line(self, directory, generated_output_file):
        return ["python", (self.project_dir + "/scripts/vector_genotype_concordance.py"), 'count',
                '-c', directory + '/{sample}.genotypes.zarr',
                '-t', directory + '/{sample}.zarr',
                '-s', 'S1',
                '--chromosomes', '2L', '2R',
                '-o', generated_output_file
                ]

    def _generate_summarize_command_line(self, directory):
        return ["python", (self.project_dir + "/scripts/vector_genotype_concordance.py"), 'summarize',
                '-i', self.project_dir + "/fixture/expected.S1.count.csv",
                '-o', directory + '/S1.output'
                ]

    def _compare_output_files(self, expected_output_file, generated_output_file):
        with open(generated_output_file, mode="r") as output, \
                open(expected_output_file, mode="r") as expected:
            output_lines = output.readlines()
            expected_lines = expected.readlines()
            for el, ol in zip(expected_lines, output_lines):
                self.assertEqual(el, ol)

    def _run_process(self, process_to_run):
        print(" ".join(process_to_run))
        with subprocess.Popen(process_to_run,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              universal_newlines=True) as proc:
            proc.wait()
            print(proc.stderr.read())
            print(proc.stdout.read())
            if proc.returncode != 0:
                self.fail("The subprocess failed")


class GTValueGenerator:

    def __init__(self):
        nucleotides = [int8(0), int8(1), int8(2), int8(3)]
        valid_genotypes = [numpy.array([a, b]) for (a, b) in itertools.product(nucleotides, nucleotides)]
        valid_genotypes.append(numpy.array([int8(-1), int8(-1)]))
        possible_combinations = [x for x in itertools.product(valid_genotypes, valid_genotypes)]
        self.left_side = [a for (a, _) in possible_combinations]
        self.right_side = [a for (_, a) in possible_combinations]
        self.len = len(possible_combinations)

    def __left_value(self, i, j, k):
        return self.left_side[(i % self.len)][k]

    def __right_value(self, i, j, k):
        return self.right_side[i % self.len][k]

    def generate_left(self, size):
        return self.__generate(self.__left_value, size)

    def __generate(self, side, size):
        return fromfunction(vectorize(side), (size, 1, 2), dtype='int8')

    def generate_right(self, size):
        return self.__generate(self.__right_value, size)


def _generate_gts(sample, chrom, size):
    group = sample.create_group(chrom)
    calldata = group.create_group('calldata')
    calldata['GT'] = zarr.zeros((size, 1, 2), dtype='int8')
    return calldata['GT']


def _generate_zarr(scratch_name, modifier):
    dest = zarr.open_group(scratch_name, mode='w')
    sample = dest.create_group('S1')
    _generate_gts(sample, '2R', 60132453)
    gts = _generate_gts(sample, '2L', 48525747)
    modifier(gts)


def _modify_gt_array(generator_function, gt_array):
    step = 100_000
    for j in range(0, 1_000_000, step):
        gt_array[j:j + step, :, :] = generator_function(step)


def generate_zarr(output, generator):
    _generate_zarr(output, partial(_modify_gt_array, generator))


if __name__ == '__main__':
    unittest.main()
