samplelist=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
bamfilefolder=$rootfolder/bamlinks/phase3
outputfolder=$rootfolder/$study/breakpoints
logfolder=$outputfolder/logfolders/breakpoint_detection
errorfolder=$outputfolder/errorfolders/breakpoint_detection

mkdir -p $outputfolder
echo "Identifying breakpoint reads for samples from ${samplelist} on `date`." >> $outputfolder/added_samples.log

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder/2R/Ace1_region/breakpointlogs
mkdir -p $outputfolder/2R/Cyp6_region/breakpointlogs
mkdir -p $outputfolder/3R/Cyp6zm_region/breakpointlogs
mkdir -p $outputfolder/3R/Gste_region/breakpointlogs
mkdir -p $outputfolder/X/Cyp9k1_region/breakpointlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $samplelist))

bsub -J "coverageArray[1-$numbams]" \
     -q small \
     -o $logfolder/breakpoint_output_%J.%I.txt \
     -e $errorfolder/breakpoint_error_%J.%I.txt \
     ' '${scriptsfolder}'/breakpoints.sh '${bamfilefolder}' '${samplelist}' ${LSB_JOBINDEX} '$outputfolder

