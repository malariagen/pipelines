import argparse
import os


"""
Given a file with the length of each chromosome, 
creates a dictionary mapping chromosome name (as a string) 
to chromosome length (as an int)
"""
def get_chromomsome_lengths(chromosome_lengths_file):
    with open(chromosome_lengths_file) as f:
        chromosome_lengths = {l[0]: int(l[1]) for l in [line.strip().split('\t') for line in f]}
    return chromosome_lengths


"""
Given dictionary of chromosome name to chromosome length, 
creates a dictionary of chromosome name to a list of intervals. 
The size of the intervals and overlap are specified as arguments. 
"""
def get_intervals(chromosome_lengths, interval_size, overlap_size):
    intervals = {}
    for ch in chromosome_lengths:
        length = chromosome_lengths[ch]
        intervals[ch] = ["{}:{}-{}".format(ch, start, stop)
                         if stop <= length 
                         else "{}:{}-{}".format(ch, start, length) 
                         for start, stop in 
                         [(x, x + interval_size + overlap_size)  
                         for x in range(0, length + 1, interval_size)]]
    return intervals


"""
Write an interval list for each chromosome to a given output directory. 
"""
def write_interval_lists(intervals, output_dir):
    for ch in intervals:
        filename = os.path.join(output_dir, ch + "_interval.txt")
        with open(filename, 'w') as f:
            f.writelines("{}\n".format(interval) for interval in intervals[ch]) 


"""
Creates interval lists in the specified output directory for the chromosomes listed in the chromosome lenghts file.
An example of a chromosome lengths file is at:
gs://malariagen/references/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.chr_lengths.txt
Default interval size is 2Mb and default overlap size is 400Kb (note that arguments to the script are in bp).
"""
def main():
    parser = argparse.ArgumentParser(description="Creat an interval list for each chromosome")
    parser.add_argument('--chromosome-lengths-file',
                        dest='chromosome_lengths_file',
                        required=True,
                        help="Path to chromosome lengths file")
    parser.add_argument('--interval-size',
                        dest='interval_size',
                        default=2000000,
                        type=int,
                        help="Size of intervals in bp")
    parser.add_argument('--overlap-size',
                        dest='overlap_size',
                        default=400000,
                        type=int,
                        help="Size of overlap in bp")
    parser.add_argument('--output-dir',
                        dest='output_dir',
                        required=True,
                        help="Directory to write interval lists to")
    

    args = parser.parse_args()

    chromosome_lengths_file = args.chromosome_lengths_file
    interval_size = args.interval_size
    overlap_size = args.overlap_size
    output_dir = args.output_dir
    
    # Check that output_dir exists and is writeable
    try:
        assert(os.access(output_dir, os.W_OK))
    except AssertionError:
        print("The output directory specifed does not exist or is not writeable.")
        sys.exit(1)

    chromosome_lengths = get_chromomsome_lengths(chromosome_lengths_file)
    intervals = get_intervals(chromosome_lengths, interval_size, overlap_size)
    write_interval_lists(intervals, output_dir)


if __name__ == '__main__':
    main()

