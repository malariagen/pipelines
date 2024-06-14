samplelist=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
bamfilefolder=$rootfolder/bamlinks/phase3
outputfolder=$rootfolder/$study/coverage
logfolder=$outputfolder/logfolders/coverage_calculation
errorfolder=$outputfolder/errorfolders/coverage_calculation
allchrom=(2L 2R 3L 3R X)

mkdir -p $outputfolder
echo "Calculating coverage for samples from ${samplelist} on `date`." >> $outputfolder/added_samples.log

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder
for chrom in ${allchrom[@]}
do
	mkdir -p $outputfolder/$chrom/coveragelogs
done

# Get the number of bamfiles that need processing
numbams=($(wc -l $samplelist))

bsub -J "coverageArray[1-$numbams]" \
     -R"select[mem>300] rusage[mem=300]" \
     -M300 \
     -o $logfolder/coverage_output_%J.%I.txt \
     -e $errorfolder/coverage_error_%J.%I.txt \
     ' '${scriptsfolder}'/get_windowed_coverage.sh '${bamfilefolder}' '${samplelist}' ${LSB_JOBINDEX} '$outputfolder

