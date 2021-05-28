# This script will calculate the coverage at every x bp along a genomic region for a given bamfile

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pysam

# Check how many arguments you have and do anything necessary
if (len(argv) == 5):
	input_filename = argv[1]
	chromosome = argv[2]
	interval = argv[3]
	window_size = argv[4]
	minqual = 10
	output_filename = 'counts_for_HMM_python3.txt'
elif (len(argv) == 6):
	input_filename = argv[1]
	chromosome = argv[2]
	interval = argv[3]
	window_size = argv[4]
	minqual = argv[5]
	output_filename = 'counts_for_HMM_python3.txt'
elif (len(argv) == 7):
	input_filename = argv[1]
	chromosome = argv[2]
	interval = argv[3]
	window_size = argv[4]
	minqual = argv[5]
	output_filename = argv[6]
else:
	raise Exception("Fail. There should be between four and six command line arguments (input_filename, chromosome, interval, window_size [, minqual [, output_filename]]]]])")


print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n')
print('Input arguments:\n\tinput_filename: ' + input_filename + '\n\tchromosome: ' + chromosome + '\n\tinterval: ' + str(interval) + '\n\twindow_size: ' + str(window_size) + '\n\tminqual: ' + str(minqual) + '\n\toutput_filename: ' + output_filename + '\n\n')
stdout.flush()

print('Loading alignment file ' + input_filename)
stdout.flush()
bambam = pysam.Samfile(input_filename, 'rb')

# Check that the requested chromosome is a valid reference and get its length
if chromosome not in bambam.references:
	raise Exception('Fail. Requested chromosome was not present in the reference')
else:
	chrom_id = bambam.get_tid(chromosome)
	chromosome_length = bambam.lengths[chrom_id]

print('Calculating coverage along chromosome ' + chromosome + ' from start to end (' + str(chromosome_length) + ') in intervals of ' + str(interval) + 'bp.')
stdout.flush()

# We start from position 1 and move up through the positions in steps of size "interval", calculating
# the coverage at each step. 
positions = range(0, chromosome_length-1, int(interval))
output = []
output_mapq0 = []
output_belowT = []
output_total = []
for i in positions:
	alignments = bambam.fetch(chromosome, i, i+int(window_size))
	count = 0
	count_mapq0 = 0
	count_belowT = 0
	for al in alignments:
		if not al.is_unmapped:
			# We count coverage as reads that start in the current interval
			if al.reference_start >= i:
				if al.mapq >= float(minqual):
					count += 1
				else:
					if al.mapq == 0:
						count_mapq0 += 1
					else:
						count_belowT += 1
	output += [str(count)]
	output_mapq0 += [str(count_mapq0)]
	output_belowT += [str(count_belowT)]
	output_total += [str(count + count_mapq0 + count_belowT)]
outputfile = open(output_filename, 'w')
outputfile.write('Position\tCounts mapq >= ' + str(minqual) + '\tCounts 0 < mapq < ' + str(minqual) + '\tCounts mapq0\tCounts total\n')
for i in range(len(positions)):
	outputfile.write(str(positions[i]) + '\t' + output[i] + '\t' + output_belowT[i] + '\t' + output_mapq0[i] + '\t' + output_total[i] + '\n')
outputfile.close()

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()

