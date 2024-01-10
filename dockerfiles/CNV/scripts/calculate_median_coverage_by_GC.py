#!/usr/bin/python3

# This script calculates the median coverage by GC window for a list of samples from Ag1000G

import time  # required to output the time at which the script was run
from sys import stdout  # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv  # this import is needed in order for the script to handle command line arguments
import socket  # this import allows us to access the name of the computer that the script is being run on
import numpy as np
import pandas as pd
from re import *


# Write a function to normalised the coverage for each GC content bin
def normalise_coverage_by_GC(coverage, median_coverage_by_GC, ploidy=2):
	output = coverage.copy()
	# For each counts value in the output object, associate it with the median coverage for its GC bin
	output['expcov'] = [median_coverage_by_GC.loc[x] for x in output['GC']]
	# Now divide each counts value by the coverage associated with its GC bin, we multiply by 2 because
	# that was the ploidy of our median coverage.
	output['Normalised_coverage'] = ploidy * output['Counts total'] / output['expcov']
	# In very rare cases, the expected coverage ends up being 0 (because none of the few windows with a given
	# GC content have any coverage). In these cases, we manually set the coverage to be 0 (rather than NaN)
	which_zeros = (output['expcov'] == 0)
	output.loc[(which_zeros, 'Normalised_coverage')] = 0
	# Return the new dataframe 
	return output


# Write a function to load and mask the counts data for a given chromosome in a given sample
def load_and_mask(working_folder, chrom, sample_name, gc, filterpass):
	this_filename = working_folder + '/' + chrom + '/counts_for_HMM_' + sample_name + '_' + chrom + '_output.csv'
	print('\tLoading file ' + this_filename)
	stdout.flush()
	these_counts = pd.read_csv(this_filename, sep='\t')
	# Check that the number of bins for the coverage is correct (we have already checked that accessibility
	# and GC have the same number of bins).
	if these_counts.shape[0] != gc.shape[0]:
		raise Exception(
			'Fail. Coverage and GC_content on chromosome ' + chrom + ' should have the same number of bins.')
	# Add the GC content information to the counts table. 
	these_counts['GC'] = gc.iloc[:, 2].values
	these_counts['Chrom'] = chrom
	# Get the accessibility-masked coverage
	output_masked_counts = these_counts.loc[np.array(filterpass), :]
	return (output_masked_counts)


# Check how many arguments you have and do anything necessary
if (len(argv) == 9):
	accessibility_threshold = float(argv[1])
	accessibility_mask_file = argv[2]
	mapq0_threshold = float(argv[3])
	mapq0_file = argv[4]
	sample_manifest = argv[5]
	gc_content_filename = argv[6]
	working_folder = argv[7]
	output_file_key = argv[8]
else:
	raise Exception(
		"Fail. There should be eight command line arguments (accessibility_threshold, accessibility_mask_file, mapq0_threshold, mapq0_file, sample_manifest, gc_content_filename, working_folder, output_file_key)"
	)

print('Running ' + argv[0] + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime(
	"%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
print('Input arguments:')
print('\taccessibility_threshold: ' + str(accessibility_threshold))
print('\taccessibility_mask_file: ' + accessibility_mask_file)
print('\tmapq0_threshold: ' + str(mapq0_threshold))
print('\tmapq0_file: ' + mapq0_file)
print('\tsample_manifest: ' + sample_manifest)
print('\tgc_content_filename: ' + gc_content_filename)
print('\tworking_folder: ' + working_folder)
print('\toutput_file_key: ' + output_file_key + '\n')
stdout.flush()

# Set the chromosomes 
autosomes = ['2L', '2R', '3L', '3R']
sex_chrom = 'X'
chroms = autosomes + [sex_chrom]

with open(sample_manifest, 'r') as f:
	sample_names = [x.rstrip('\n') for x in f.readlines()]

mapq0_prop = pd.read_csv(mapq0_file, sep='\t')
mapq0_prop['prop_mapq0'] = mapq0_prop['Count mapq = 0'] / (mapq0_prop['Count mapq > 0'] + mapq0_prop['Count mapq = 0'])
# Create the mask for mapq0. There are windows where prop_mapq0 is NA (because there are no reads mapping in that window. 
# They get removed later anyway because they are inaccessible, but even in the event that they are not removed, the 
# output of prop_mapq0 <= mapq0_threshold will be False in these cases. 
mapq0_prop['filterpass'] = mapq0_prop.prop_mapq0 <= mapq0_threshold
mapq0_prop_grouped = mapq0_prop.groupby('Chrom')

# Load up the GC composition for the Agam genome. 
gc_all = pd.read_csv(gc_content_filename, sep='\t', header=None)
# Take the floor of the percentage GC content
gc_all.iloc[:, 2] = (gc_all.iloc[:, 2] * 100).astype(int)
gc_grouped = gc_all.groupby(0)

# Load up the accessibility data
acc_all = pd.read_csv(accessibility_mask_file, sep='\t')
acc_all['filterpass'] = acc_all.Mean_accessibility >= accessibility_threshold
acc_all['doublefilterpass'] = acc_all.filterpass & mapq0_prop.filterpass
acc_grouped = acc_all.groupby('Chrom')

# Check that the GC content and accessibility have the same number of bins for each chromosome
for chrom in chroms:
	if gc_grouped.get_group(chrom).shape[0] != acc_grouped.get_group(chrom).shape[0]:
		raise Exception('Fail. GC content and accessibility on chromosome ' + chrom + ' should have the same number of bins.')

# Get the number of each GC bin after filtering
gc_mainchroms = gc_all.loc[gc_all[0].apply(lambda x: x in chroms), :].copy()
gc_mainchroms['filterpass'] = np.array(acc_all.doublefilterpass)
# I don't know why, but "query" doesn't work here, so need to do this in two steps. 
output_table = pd.DataFrame(gc_mainchroms.groupby(2).apply(lambda x: sum(x.filterpass)), columns=['bin_freq'])# .query('bin_freq > 0')
output_table = output_table.loc[output_table.bin_freq > 0, :].copy()
output_table.index = output_table.index.rename('GC')

# Here is where we will store the variance data
output_variance = pd.DataFrame(0, columns=chroms + ['autosomes'], index=sample_names)

# Now for each sample get the counts data for all the chromosomes
for this_sample in sample_names:
	print('\nAnalysing sample ' + this_sample + '.')
	stdout.flush()
	# Load and mask the data for each chromosome
	autosomal_masked_raw_counts = pd.concat(map(lambda c: load_and_mask(working_folder,
	                                                                    c,
	                                                                    this_sample,
	                                                                    gc_grouped.get_group(c),
	                                                                    acc_grouped.get_group(c).doublefilterpass),
												autosomes))
	sexchrom_masked_raw_counts = load_and_mask(working_folder, sex_chrom, this_sample,
	                                           gc_grouped.get_group(sex_chrom),
	                                           acc_grouped.get_group(sex_chrom).doublefilterpass)
		
	# Now we have combined the results from all of the autosomes, we can group the coverage by GC bin
	# and then compute the median
	counts_by_GC = autosomal_masked_raw_counts.groupby('GC')
	median_counts_by_GC = counts_by_GC['Counts total'].median()
	# Add these to the global data
	output_table[this_sample] = median_counts_by_GC
	
	# Get the normalised coverage 
	autosomal_masked_counts = normalise_coverage_by_GC(autosomal_masked_raw_counts, median_counts_by_GC)
	# We still use the autosomal median_counts_by_GC for the X, because we want to normalised by the diploid expectation
	sexchrom_masked_counts = normalise_coverage_by_GC(sexchrom_masked_raw_counts, median_counts_by_GC)
	
	# Now we want to calculate the variance in coverage for this sample, both for each chromosome
	# and overall. We filter out the 1% of windows with the highest coverage before calculating
	# the variance. 
	print('\tCalculating normalised coverage and variance.')
	stdout.flush()
	masked_counts_by_autosome = autosomal_masked_counts.groupby('Chrom')
	
	def filter_by_quantile(x, thresh):
		return (x[x < np.quantile(x, thresh)])
	
	quantile_threshold = 0.99
	for chrom in autosomes:
		output_variance.loc[this_sample, chrom] = np.var(filter_by_quantile(masked_counts_by_autosome.get_group(chrom).Normalised_coverage, quantile_threshold))
	output_variance.loc[this_sample, sex_chrom] = np.var(filter_by_quantile(sexchrom_masked_counts.Normalised_coverage, quantile_threshold))
	# Now calculate the variance across the autosomes
	output_variance.loc[this_sample, 'autosomes'] = np.var(filter_by_quantile(autosomal_masked_counts.Normalised_coverage, quantile_threshold))

# Write the output to file
output_filename = working_folder + '/median_coverage_by_GC_masked_' + sub('\.', '', str(accessibility_threshold)) + '_' + sub('\.', '', str(mapq0_threshold)) + '_' + output_file_key + '.csv'
output_variance_filename = working_folder + '/coverage_variance_masked_' + sub('\.', '', str(accessibility_threshold)) + '_' + sub('\.', '', str(mapq0_threshold)) + '_' + output_file_key + '.csv'

print('\nSaving coverage output to file ' + output_filename)
stdout.flush()
output_table.to_csv(output_filename, sep='\t')

print('\nSaving variance output to file ' + output_variance_filename)
stdout.flush()
output_variance.to_csv(output_variance_filename, sep='\t')

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')

