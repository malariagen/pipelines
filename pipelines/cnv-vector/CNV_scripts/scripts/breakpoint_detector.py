# breakpoint_detector.py

# This script goes through an alignment file and records the positions at which soft_clipping is detected
# in the aligned reads

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pysam

# Run the script:

if (len(argv) == 4):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename_root = 'breakpoint_detector'
	minqual = 10
elif (len(argv) == 5):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename_root = argv[4]
	minqual = 10
elif (len(argv) == 6):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename_root = argv[4]
	minqual = argv[5]
else:
	raise Exception("Fail. There should be between 3 and  5 command line arguments (input_filename, chromosome, region [, output_filename_root [,minqual]])")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n')
print('Input arguments:\n\tinput_filename: ' + input_filename + '\n\tchromosome: ' + chromosome + '\n\tregion: ' + region + '\n\tminqual: ' + str(minqual) + '\n\toutput_filename_root: ' + output_filename_root + '\n\n')
stdout.flush()

print('Loading alignment file ' + input_filename)
stdout.flush()
bambam = pysam.Samfile(input_filename, 'rb')

output_filename = output_filename_root + '.csv'
pre_clipping_fastq = output_filename_root + '_preclipping.fastq'
pre_clipping_fastq_file = open(pre_clipping_fastq, 'w')
post_clipping_fastq = output_filename_root + '_postclipping.fastq'
post_clipping_fastq_file = open(post_clipping_fastq, 'w')

# Check that the requested chromosome is a valid reference and get its length
if chromosome not in bambam.references:
	raise Exception('Fail. Requested chromosome was not present in the reference')
else:
	chrom_id = bambam.get_tid(chromosome)
	chromosome_length = bambam.lengths[chrom_id]

print('Going through chromosome, looking for soft-clipped reads')
stdout.flush()

soft_clipping_start_points = []
soft_clipping_start_points_mapq0 = []
soft_clipping_start_points_belowT = []
soft_clipping_end_points = []
soft_clipping_end_points_mapq0 = []
soft_clipping_end_points_belowT = []
count = 0
count_mapq0 = 0
count_belowT = 0

# Get the info from the bamfile
if region == ':':
	alignments = bambam.fetch(chromosome)
else:
	region_start, region_end = [int(x) for x in region.split(':')]
	alignments = bambam.fetch(chromosome, region_start, region_end)

for al in alignments:
	# Set the possible states as False
	soft_clipping_start_point = False
	soft_clipping_end_point = False
	# We check that the read is mapped
	if al.is_unmapped:
		continue  
	# The deal with the weird situation that required this change to the previous version of this script
	# (see top of the page) we check that the positions vector is not empty
	if len(al.positions) == 0:
		continue
	cigar = al.cigar
	# Check the first and last entries of the cigar to see if it is soft clipped. If there is only
	# one element in the cigar, it can't be soft clipped
	if len(cigar) == 1:
		continue
	# We require soft-clipping to be at least 10 bases long
	if (cigar[0][0] == 4) & (cigar[0][1] >= 10):
		# We have soft-clipping at the start of the read. Record it as a soft-clipping end point
		soft_clipping_end_point = al.positions[0]
		# We will output a fastq file containing the soft-clipped reads. 
		pre_clipping_fastq_file.write('@' + al.query_name + '_bases_soft_clipped_left_of_pos_' + str(soft_clipping_end_point) + '\n')
		pre_clipping_fastq_file.write(al.seq[:cigar[0][1]] + '\n')
		pre_clipping_fastq_file.write('+\n')
		pre_clipping_fastq_file.write(al.qual[:cigar[0][1]] + '\n')
	if (cigar[-1][0] == 4) & (cigar[-1][1] >= 10):
		# We have soft-clipping at the end of the read. Record it as a soft-clipping start point
		soft_clipping_start_point = al.positions[-1]+2 #(we're adding 1 because of python indexing, and 1 because al.positions[-1] is the last base that matches, not the first base that is soft-clipped
		# We will output a fastq file containing the soft-clipped reads. 
		post_clipping_fastq_file.write('@' + al.query_name + '_bases_soft_clipped_right_of_pos_' + str(soft_clipping_start_point) + '\n')
		post_clipping_fastq_file.write(al.seq[(100-cigar[-1][1]):] + '\n')
		post_clipping_fastq_file.write('+\n')
		post_clipping_fastq_file.write(al.qual[(100-cigar[-1][1]):] + '\n')

	# We check the mapq value. 
	if al.mapq >= float(minqual):
		count += 1
		if soft_clipping_start_point:
			soft_clipping_start_points += [[soft_clipping_start_point, al.seq[(100-cigar[-1][1]):]]]
		if soft_clipping_end_point:
			soft_clipping_end_points += [[soft_clipping_end_point, al.seq[:cigar[0][1]]]]
	else:
		if al.mapq == 0:
			count_mapq0 += 1
			if soft_clipping_start_point:
				soft_clipping_start_points_mapq0 += [[soft_clipping_start_point, al.seq[(100-cigar[-1][1]):]]]
			if soft_clipping_end_point:
				soft_clipping_end_points_mapq0 += [[soft_clipping_end_point, al.seq[:cigar[0][1]]]]
		else:
			count_belowT += 1
			if soft_clipping_start_point:
				soft_clipping_start_points_belowT += [[soft_clipping_start_point, al.seq[(100-cigar[-1][1]):]]]
			if soft_clipping_end_point:
				soft_clipping_end_points_belowT += [[soft_clipping_end_point, al.seq[:cigar[0][1]]]]
post_clipping_fastq_file.close()
pre_clipping_fastq_file.close()

print('Writing output to file ' + output_filename)
stdout.flush()
outputfile = open(output_filename, 'w')
outputfile.write('#Total mapq >= ' + str(minqual) + ': ' + str(count) + '\n')
outputfile.write('#Total mapq < ' + str(minqual) + ': ' + str(count_belowT) + '\n')
outputfile.write('#Total mapq0: ' + str(count_mapq0) + '\n')
outputfile.write('Type\tPosition\tClipped_sequence\n')
for pos in soft_clipping_start_points:
	outputfile.write('soft clipping start point mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in soft_clipping_start_points_mapq0:
	outputfile.write('soft clipping start point0\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in soft_clipping_start_points_belowT:
	outputfile.write('soft clipping start point mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in soft_clipping_end_points:
	outputfile.write('soft clipping end point mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in soft_clipping_end_points_mapq0:
	outputfile.write('soft clipping end point0\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in soft_clipping_end_points_belowT:
	outputfile.write('soft clipping end point mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
outputfile.close()

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()



