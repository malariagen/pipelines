# This shell script was originally written by Eric Lucas and has been modified by Kevin Palis <palis@broadinstitute.org> to run in Terra
# as part of the CNV pipeline.


# bamfilefolder=$1
# bamfilepaths=$2
# samplenum=$3 
# outputfolder=$4

# scriptsfolder=~/scripts/CNV_scripts/scripts

# samplename=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f1))
# bampath=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f2))
# bamfile=${bamfilefolder}/${samplename}.bam

bamfile = $1
samplename = $2
outputfolder = $3


#echo Downloading bam file for $samplename from $bampath to $bamfile

# Removing with the assumption that in Terra, the input files are already localized for us
# Download the required bam
#curl ${bampath} > ${bamfile} 2> ${bamfilefolder}/${samplename}_download.log
#curl ${bampath}.bai > ${bamfile}.bai 2>> ${bamfilefolder}/${samplename}_download.log

# Start the conda env
source activate cnv37 

# Get the coverage data
echo Calculating coverage
allchrom=(2L 2R 3L 3R X)
coveragefolder=$outputfolder/coverage
for chrom in ${allchrom[@]}
do
	#create the directory structure needed (this was likely pre-provisioned in the baremetal version of this pipeline)
	mkdir -p ${coveragefolder}/${chrom}/coveragelogs/
	python ${scriptsfolder}/counts_for_HMM.py \
	       $bamfile \
	       $chrom \
	       300 300 10 \
	       ${coveragefolder}/${chrom}/counts_for_HMM_${samplename}_${chrom}_output.csv \
	       > ${coveragefolder}/${chrom}/coveragelogs/counts_for_HMM_${samplename}_${chrom}.log 2>&1
done

# Commenting out all the code for discordant reads out as we've identified that we don't need this for this pipeline.
# TODO: Remove the code block entirely once this sub-pipeline is in working order.
# echo Identifying discordant reads
# Get the discordant reads
# SSFA_script=$scriptsfolder/SSFA.py
# SSFAfolder=$outputfolder/diagnostic_reads/SSFA
# python $SSFA_script $bamfile 2R 3425000:3650000 ${SSFAfolder}/2R/Ace1_region/${samplename}_Ace1_SSFA_output.csv 10 > ${SSFAfolder}/2R/Ace1_region/SSFAlogs/${samplename}_Ace1_SSFA_output.log 2>&1
# python $SSFA_script $bamfile 2R 28460000:28570000 ${SSFAfolder}/2R/Cyp6_region/${samplename}_CYP6_SSFA_output.csv 10 > ${SSFAfolder}/2R/Cyp6_region/SSFAlogs/${samplename}_CYP6_SSFA_output.log 2>&1
# python $SSFA_script $bamfile 3R 6900000:7000000 ${SSFAfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_SSFA_output.csv 10 > ${SSFAfolder}/3R/Cyp6zm_region/SSFAlogs/${samplename}_CYP6ZM_SSFA_output.log 2>&1
# python $SSFA_script $bamfile 3R 28570000:28620000 ${SSFAfolder}/3R/Gste_region/${samplename}_GST_SSFA_output.csv 10 > ${SSFAfolder}/3R/Gste_region/SSFAlogs/${samplename}_GST_SSFA_output.log 2>&1
# python $SSFA_script $bamfile X 15220000:15255000 ${SSFAfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_SSFA_output.csv 10 > ${SSFAfolder}/X/Cyp9k1_region/SSFAlogs/${samplename}_CYP9K1_SSFA_output.log 2>&1

# # Get the soft clipped reads
# breakpoints_script=$scriptsfolder/breakpoint_detector.py
# breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
# python $breakpoints_script $bamfile 2R 3425000:3650000 ${breakpointsfolder}/2R/Ace1_region/${samplename}_Ace1_breakpoints_output 10 > ${breakpointsfolder}/2R/Ace1_region/breakpointlogs/${samplename}_Ace1_breakpoints_output.log 2>&1
# python $breakpoints_script $bamfile 2R 28460000:28570000 ${breakpointsfolder}/2R/Cyp6_region/${samplename}_CYP6_breakpoints_output 10 > ${breakpointsfolder}/2R/Cyp6_region/breakpointlogs/${samplename}_CYP6_breakpoints_output.log 2>&1
# python $breakpoints_script $bamfile 3R 6900000:7000000 ${breakpointsfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_breakpoints_output 10 > ${breakpointsfolder}/3R/Cyp6zm_region/breakpointlogs/${samplename}_CYP6ZM_breakpoints_output.log 2>&1
# python $breakpoints_script $bamfile 3R 28570000:28620000 ${breakpointsfolder}/3R/Gste_region/${samplename}_GST_breakpoints_output 10 > ${breakpointsfolder}/3R/Gste_region/breakpointlogs/${samplename}_GST_breakpoints_output.log 2>&1
# python $breakpoints_script $bamfile X 15220000:15255000 ${breakpointsfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_breakpoints_output 10 > ${breakpointsfolder}/X/Cyp9k1_region/breakpointlogs/${samplename}_CYP9K1_breakpoints_output.log 2>&1

# Delete the bam file 
# echo Deleting bam file
# rm ${bamfile}*

echo "Windowed Coverage finished successfully."

