bamfilefolder=$1
manifest=$2
samplenum=$3 
outputfolder=$4
scriptsfolder=~/scripts/CNV_scripts/scripts

samplename=($(head -n$samplenum $manifest | tail -n1))
bamfile=${bamfilefolder}/${samplename}.bam

echo $bamfile

source activate cnv37 

SSFA_script=$scriptsfolder/SSFA.py
SSFAfolder=$outputfolder/SSFA
python $SSFA_script $bamfile 2R 3425000:3650000 ${SSFAfolder}/2R/Ace1_region/${samplename}_Ace1_SSFA_output.csv 10 > ${SSFAfolder}/2R/Ace1_region/SSFAlogs/${samplename}_Ace1_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2R 28460000:28570000 ${SSFAfolder}/2R/Cyp6_region/${samplename}_CYP6_SSFA_output.csv 10 > ${SSFAfolder}/2R/Cyp6_region/SSFAlogs/${samplename}_CYP6_SSFA_output.log 2>&1
python $SSFA_script $bamfile 3R 6900000:7000000 ${SSFAfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_SSFA_output.csv 10 > ${SSFAfolder}/3R/Cyp6zm_region/SSFAlogs/${samplename}_CYP6ZM_SSFA_output.log 2>&1
python $SSFA_script $bamfile 3R 28570000:28620000 ${SSFAfolder}/3R/Gste_region/${samplename}_GST_SSFA_output.csv 10 > ${SSFAfolder}/3R/Gste_region/SSFAlogs/${samplename}_GST_SSFA_output.log 2>&1
python $SSFA_script $bamfile X 15220000:15255000 ${SSFAfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_SSFA_output.csv 10 > ${SSFAfolder}/X/Cyp9k1_region/SSFAlogs/${samplename}_CYP9K1_SSFA_output.log 2>&1

breakpoints_script=$scriptsfolder/breakpoint_detector.py
breakpointsfolder=$outputfolder/breakpoints
python $breakpoints_script $bamfile 2R 3425000:3650000 ${breakpointsfolder}/2R/Ace1_region/${samplename}_Ace1_breakpoints_output 10 > ${breakpointsfolder}/2R/Ace1_region/breakpointlogs/${samplename}_Ace1_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2R 28460000:28570000 ${breakpointsfolder}/2R/Cyp6_region/${samplename}_CYP6_breakpoints_output 10 > ${breakpointsfolder}/2R/Cyp6_region/breakpointlogs/${samplename}_CYP6_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 3R 6900000:7000000 ${breakpointsfolder}/3R/Cyp6zm_region/${samplename}_CYP6ZM_breakpoints_output 10 > ${breakpointsfolder}/3R/Cyp6zm_region/breakpointlogs/${samplename}_CYP6ZM_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 3R 28570000:28620000 ${breakpointsfolder}/3R/Gste_region/${samplename}_GST_breakpoints_output 10 > ${breakpointsfolder}/3R/Gste_region/breakpointlogs/${samplename}_GST_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile X 15220000:15255000 ${breakpointsfolder}/X/Cyp9k1_region/${samplename}_CYP9K1_breakpoints_output 10 > ${breakpointsfolder}/X/Cyp9k1_region/breakpointlogs/${samplename}_CYP9K1_breakpoints_output.log 2>&1

