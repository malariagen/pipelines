workingfolder=$1
manifest=$2
accessibility_file=$3
GC_content_file=$4
sample_group_id=$5
scriptsfolder=~/scripts/CNV_scripts/scripts
allchrom=(2L 2R 3L 3R X)

source activate cnv37 

startingfolder=$(pwd)

# We will combine all of the mapq0 prop files into a single file
mapqprop_output_filepath=$workingfolder/mapq_proportions_allchrom_${sample_group_id}.csv
echo -e 'Chrom\tPosition\tCount mapq > 0\tCount mapq = 0' > $mapqprop_output_filepath

# Calculate the windowed proportion of mapq0 reads across samples
echo "Calculating windowed proportion of mapq0 reads for sample group ${sample_group_id}"
for chrom in ${allchrom[@]}
do
	cd $workingfolder/$chrom
	python $scriptsfolder/get_prop_mapq0.py \
	       $manifest \
	       mapq_proportions_${chrom}_${sample_group_id}.csv \
		   . \
	       > get_prop_mapq0_${chrom}_${sample_group_id}.log 2>&1
	cd $startingfolder
	tail -n +2 $workingfolder/$chrom/mapq_proportions_${chrom}_${sample_group_id}.csv | sed "s/^/${chrom}\t/" >> $mapqprop_output_filepath
done

# Now calculate median coverage by GC 
echo "Calculating median and variance of coverage by GC bin for sample group ${sample_group_id}"
python $scriptsfolder/calculate_median_coverage_by_GC.py 0.9 \
                                                         $accessibility_file \
                                                         0.5 \
                                                         $mapqprop_output_filepath \
                                                         $manifest \
                                                         $GC_content_file \
                                                         $workingfolder \
                                                         $sample_group_id \
                                                         > $workingfolder/calculate_mean_coverage_by_GC_09_05_${sample_group_id}.log 2>&1


