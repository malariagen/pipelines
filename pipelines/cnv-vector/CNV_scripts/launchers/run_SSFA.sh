samplelist=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
bamfilefolder=$rootfolder/bamlinks/phase3
outputfolder=$rootfolder/$study/SSFA
logfolder=$outputfolder/logfolders/SSFA_detection
errorfolder=$outputfolder/errorfolders/SSFA_detection

mkdir -p $outputfolder
echo "Identifying SSFA reads for samples from ${samplelist} on `date`." >> $outputfolder/added_samples.log

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder/2R/Ace1_region/SSFAlogs
mkdir -p $outputfolder/2R/Cyp6_region/SSFAlogs
mkdir -p $outputfolder/3R/Cyp6zm_region/SSFAlogs
mkdir -p $outputfolder/3R/Gste_region/SSFAlogs
mkdir -p $outputfolder/X/Cyp9k1_region/SSFAlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $samplelist))

bsub -J "coverageArray[1-$numbams]" \
     -R"select[mem>300] rusage[mem=300]" \
     -M300 \
     -o $logfolder/SSFA_output_%J.%I.txt \
     -e $errorfolder/SSFA_error_%J.%I.txt \
     ' '${scriptsfolder}'/SSFA.sh '${bamfilefolder}' '${samplelist}' ${LSB_JOBINDEX} '$outputfolder

