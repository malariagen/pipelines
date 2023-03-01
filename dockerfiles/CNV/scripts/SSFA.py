# This script goes through an alignment file and records the positions of reads within a specified region
# whose mates map to a different chromosome or discordantly on the same chromosome

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pysam

if (len(argv) == 4):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename = 'counts_for_mHMM_python2_v3.txt'
	minqual = 10
elif (len(argv) == 5):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename = argv[4]
	minqual = 10
elif (len(argv) == 6):
	input_filename = argv[1]
	chromosome = argv[2]
	region = argv[3]
	output_filename = argv[4]
	minqual = argv[5]
else:
	raise Exception("Fail. There should be between 3 and  5 command line arguments (input_filename, chromosome, region [, output_filename [,minqual]])")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n')
print('Input arguments:\n\tinput_filename: ' + input_filename + '\n\tchromosome: ' + chromosome + '\n\tregion: ' + region + '\n\tminqual: ' + str(minqual) + '\n\toutput_filename: ' + output_filename + '\n\n')
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

print('Going through chromosome, looking for Cross-chrom, Face Away, Same Strand and Far Mapped reads')
stdout.flush()

output_crosschrom = []
output_crosschrom_mapq0 = []
output_crosschrom_belowT = []
output_FA = []
output_FA_mapq0 = []
output_FA_belowT = []
output_SS = []
output_SS_mapq0 = []
output_SS_belowT = []
output_FM = []
output_FM_mapq0 = []
output_FM_belowT = []
count = 0
count_mapq0 = 0
count_belowT = 0

# Get the info from the bamfile
if region == ':':
	alignments = bambam.fetch(chromosome)
else:
	region_start, region_end = [int(x) for x in region.split(':')]
	alignments = bambam.fetch(chromosome, region_start, region_end)

# To avoid counting the same pair twice, we record the these reads we have already counted.
done_reads = []
for al in alignments:
	# If the reads are a proper pair, just move to the next one
	if al.is_proper_pair:
		continue
	# We check that the read and its mate are both mapped
	if al.is_unmapped or al.mate_is_unmapped:
		continue  
	# Set all statuses to False
	pair_is_crosschrom = False
	pair_is_SS = False
	pair_is_FA = False
	pair_is_FM = False
	# Check we haven't already counted this pair (the two mates have the same name, so the following
	# check is correct).
	if al.query_name not in done_reads:
		done_reads += [al.query_name]
		# Look for genes that map discordantly.
		this_mate = bambam.mate(al)
		if bambam.references[this_mate.rname] != chromosome:
			pair_is_crosschrom = True
		else:
			# If the reads map to the same chromosome, look at whether they are FA, SS or FM.
			# We check whether the read pair is same strand
			if al.is_reverse == al.mate_is_reverse:
				pair_is_SS = True
			else:
				# If the positions of the read and its mate are the same, we have no information about FA so we 
				# record it as the default False. There is clearly no FM.
				if al.pos < al.mpos:
					# We check whether the read pair is face-away
					if al.is_reverse:
						pair_is_FA = True
					else:
						# We check how far the reads are apart from each other. If it's more than 1000bp, we consider it FM
						if al.mpos > (al.pos + 1000):
							pair_is_FM = True
						else:
							continue
				elif al.pos > al.mpos:
					# We check whether the read pair is face-away
					if al.is_reverse:
						# We check how far the reads are apart from each other. If it's more than 1000bp, we consider it FM
						if al.pos > (al.mpos + 1000):
							pair_is_FM = True
						else:
							continue
					else:
						pair_is_FA = True
	else:
		continue
	# We check the mapq value. In order to pass a certain threshold, the mapq of both the read and its mate need to be
	# high enough
	lowest_mapq = min([al.mapq, this_mate.mapq])
	if lowest_mapq >= float(minqual):
		count += 1
		if pair_is_crosschrom:
			output_crosschrom += [[al.pos, bambam.references[this_mate.rname], al.mpos]]
		if pair_is_SS:
			output_SS += [[al.pos, al.mpos]]
		if pair_is_FA:
			if al.pos == al.mpos:
				raise Exception('This shouldn\'t be happening')
			output_FA += [[al.pos, al.mpos]]
		if pair_is_FM:
			if al.pos == al.mpos:
				raise Exception('This shouldn\'t be happening')
			output_FM += [[al.pos, al.mpos]]
	else:
		if lowest_mapq == 0:
			count_mapq0 += 1
			if pair_is_crosschrom:
				output_crosschrom_mapq0 += [[al.pos, bambam.references[this_mate.rname], al.mpos]]
			if pair_is_SS:
				output_SS_mapq0 += [[al.pos, al.mpos]]
			if pair_is_FA:
				if al.pos == al.mpos:
					raise Exception('This shouldn\'t be happening')
				output_FA_mapq0 += [[al.pos, al.mpos]]
			if pair_is_FM:
				if al.pos == al.mpos:
					raise Exception('This shouldn\'t be happening')
				output_FM_mapq0 += [[al.pos, al.mpos]]
		else:
			count_belowT += 1
			if pair_is_crosschrom:
				output_crosschrom_belowT += [[al.pos, bambam.references[this_mate.rname], al.mpos]]
			if pair_is_SS:
				output_SS_belowT += [[al.pos, al.mpos]]
			if pair_is_FA:
				if al.pos == al.mpos:
					raise Exception('This shouldn\'t be happening')
				output_FA_belowT += [[al.pos, al.mpos]]
			if pair_is_FM:
				if al.pos == al.mpos:
					raise Exception('This shouldn\'t be happening')
				output_FM_belowT += [[al.pos, al.mpos]]

print('Writing output to file ' + output_filename)
stdout.flush()
outputfile = open(output_filename, 'w')
outputfile.write('#Total mapq >= ' + str(minqual) + ': ' + str(count) + '\n')
outputfile.write('#Total mapq < ' + str(minqual) + ': ' + str(count_belowT) + '\n')
outputfile.write('#Total mapq0: ' + str(count_mapq0) + '\n')
outputfile.write('Type\tPosition\tMate position\n')
for pos in output_crosschrom:
	outputfile.write('crosschrom mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + pos[1] + ':' + str(pos[2]) + '\n')
for pos in output_crosschrom_mapq0:
	outputfile.write('crosschrom mapq0\t' + str(pos[0]) + '\t' + pos[1] + ':' + str(pos[2]) + '\n')
for pos in output_crosschrom_belowT:
	outputfile.write('crosschrom mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + pos[1] + ':' + str(pos[2]) + '\n')
for pos in output_SS:
	outputfile.write('SS mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_SS_mapq0:
	outputfile.write('SS mapq0\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_SS_belowT:
	outputfile.write('SS mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FA:
	outputfile.write('FA mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FA_mapq0:
	outputfile.write('FA mapq0\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FA_belowT:
	outputfile.write('FA mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FM:
	outputfile.write('FM mapq >= ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FM_mapq0:
	outputfile.write('FM mapq0\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
for pos in output_FM_belowT:
	outputfile.write('FM mapq < ' + str(minqual) + '\t' + str(pos[0]) + '\t' + str(pos[1]) + '\n')
outputfile.close()

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
stdout.flush()


