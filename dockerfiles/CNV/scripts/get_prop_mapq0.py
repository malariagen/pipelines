#!/usr/bin/python3
#get_prop_mapq0.py

# This script takes the coverage data for all individuals and calculates the proportion of reads that have mapq=0
# across all individuals in the dataset. 

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from re import *
import pandas as pd
import os

if (len(argv) == 3):
	list_of_samples = argv[1]
	output_filename = argv[2]
	coverage_counts_directory = '.'
elif (len(argv) == 4):
	list_of_samples = argv[1]
	output_filename = argv[2]
	coverage_counts_directory = argv[3]
else:
	raise Exception("Fail. There should be one or two command line arguments (list_of_samples, output_filename [, coverage_counts_directory])")


print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
stdout.flush()

print('Input arguments:')
print('\tlist_of_samples: ' + list_of_samples)
print('\toutput_filename: ' + output_filename)
print('\tcoverage_counts_directory: ' + coverage_counts_directory + '\n')
stdout.flush()

# Get file objects for all of the files
with open(list_of_samples, 'r') as f:
	samples = [x.rstrip('\n') for x in f.readlines()]

# Identify the files relating to the sample manifest
all_files = os.listdir(coverage_counts_directory)
files_to_load = [os.path.join(coverage_counts_directory, x) for x in all_files for y in samples if search(y + '.*output.csv', x)]
# Check all of the samples are represented:
if len(files_to_load) != len(samples):
	raise Exception('Manifest does not match available files in folder')

# Get the data from the first file
print('Loading ' + files_to_load[0])
stdout.flush()
ftable = pd.read_csv(files_to_load[0], sep = '\t')
positions = ftable['Position']
total_mapqpositive = ftable['Counts mapq >= 10'] + ftable['Counts 0 < mapq < 10']
total_mapq0 = ftable['Counts mapq0']
counts_table = pd.concat([positions, total_mapqpositive.copy(), total_mapq0.copy()], 1)
del(ftable)

def add_to_table(input_table, filename):
	print('Loading ' + filename)
	stdout.flush()
	ftable = pd.read_csv(filename, sep = '\t')
	if not all(ftable['Position'] == input_table['Position']):
		raise Exception('Positions do not match between input table and new table')
	total_mapqpositive = ftable['Counts mapq >= 10'] + ftable['Counts 0 < mapq < 10']
	total_mapq0 = ftable['Counts mapq0']
	input_table.iloc[:, 1:] = input_table.iloc[:, 1:] + pd.concat([total_mapqpositive, total_mapq0], 1)
	return(input_table)
	

# Add the data from subsequent files
for f in files_to_load[1:]:
	counts_table = add_to_table(counts_table, f)

# Change the column names to match previous versions of the script and subsequent scripts in the pipeline
counts_table.columns = ['Position', 'Count_mapq_gr0', 'Count_mapq_eq0']

print('Writing output file')
stdout.flush()
counts_table.to_csv(output_filename, sep = '\t', index = False)

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()

