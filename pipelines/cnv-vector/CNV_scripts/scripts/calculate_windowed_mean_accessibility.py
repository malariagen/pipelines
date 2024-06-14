#!/usr/bin/python

import time # required to output the time at which the script was run
import os
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import numpy as np
import pandas as pd
from re import *
import zarr
import allel


if (len(argv) == 4):
	input_accessibility_filename = argv[1]
	positions_filename = argv[2]
	output_filename = argv[3]
	window_size = 300
if (len(argv) == 5):
	input_accessibility_filename = argv[1]
	positions_filename = argv[2]
	output_filename = argv[3]
	window_size = int(argv[4])
else:
	raise Exception("Fail. There should be three or four command line arguments (input_accessibility_filename, positions_filename, output_filename, [, window_size])")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')
print('Input arguments:\n\tinput_accessibility_filename: ' + input_accessibility_filename + '\n\tpositions_filename: ' + positions_filename + '\n\toutput_filename: ' + output_filename + '\n\twindow_size: ' + str(window_size) + '\n\n')
stdout.flush()

chroms = ['2L', '2R', '3L', '3R', 'X']
acc = zarr.open(os.path.expanduser(input_accessibility_filename), mode = 'r')
pos = zarr.open(os.path.expanduser(positions_filename), mode = 'r')

outfile = open(output_filename, 'w')
outfile.write('Chrom\tPosition\tMean_accessibility\n')
for chrom in chroms:
	print('Calculating average accessibility over chromosome ' + chrom + '.\n')
	this_acc = acc[chrom + '/variants/filter_pass']
	these_pos = allel.SortedIndex(pos[chrom + '/variants/POS'])
	# Do all but the last window in the for loop
	for p in range(1,these_pos[-1] - 299, window_size):
		endp = p + window_size - 1
		# If there are no windows in the range with a value, then there are no accessible positions in the 
		# window
		try:
			this_window_acc = this_acc[these_pos.locate_range(p, endp)]
		except KeyError:
			this_window_acc = 0
		# We always divide by the window size, even if there are fewer than window_size accessibility
		# values in this interval. That's because any missing values would be N positions, which should
		# not be counted as accessible. 
		outfile.write(chrom + '\t' + str(p - 1) + '\t' + str(float(np.sum(this_window_acc))/window_size) + '\n')
	# Do the last window on its own, because it may be smaller than window_size
	p = endp + 1
	this_window_acc = this_acc[these_pos.locate_range(p, p + window_size - 1)]
	window_len = len(this_window_acc)
	outfile.write(chrom + '\t' + str(p - 1) + '\t' + str(float(np.sum(this_window_acc)) / window_len) + '\n')

outfile.close()

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()

