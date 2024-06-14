samplelist=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
bamfilefolder=$rootfolder/bamlinks/phase3
outputfolder=$rootfolder/$study/diagnostic_reads
logfolder=$outputfolder/logfolders/diagnostic_reads_detection
errorfolder=$outputfolder/errorfolders/diagnostic_reads_detection
SSFAfolder=$outputfolder/SSFA
breakpointsfolder=$outputfolder/breakpoints

mkdir -p $outputfolder
echo "Identifying diagnostic reads for samples from ${samplelist} on `date`." >> $outputfolder/added_samples.log

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $SSFAfolder/2R/Ace1_region/SSFAlogs
mkdir -p $SSFAfolder/2R/Cyp6_region/SSFAlogs
mkdir -p $SSFAfolder/3R/Cyp6zm_region/SSFAlogs
mkdir -p $SSFAfolder/3R/Gste_region/SSFAlogs
mkdir -p $SSFAfolder/X/Cyp9k1_region/SSFAlogs
mkdir -p $breakpointsfolder/2R/Ace1_region/breakpointlogs
mkdir -p $breakpointsfolder/2R/Cyp6_region/breakpointlogs
mkdir -p $breakpointsfolder/3R/Cyp6zm_region/breakpointlogs
mkdir -p $breakpointsfolder/3R/Gste_region/breakpointlogs
mkdir -p $breakpointsfolder/X/Cyp9k1_region/breakpointlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $samplelist))

bsub -J "coverageArray[1-$numbams]" \
     -R"select[mem>300] rusage[mem=300]" \
     -M300 \
     -o $logfolder/diagnostic_reads_output_%J.%I.txt \
     -e $errorfolder/diagnostic_reads_error_%J.%I.txt \
     ' '${scriptsfolder}'/get_diagnostic_reads.sh '${bamfilefolder}' '${samplelist}' ${LSB_JOBINDEX} '$outputfolder

