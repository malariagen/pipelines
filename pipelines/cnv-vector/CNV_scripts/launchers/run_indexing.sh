samplelist=$1

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
bamfilefolder=$rootfolder/bamlinks/phase3
logfolder=$rootfolder/logfolders/bam_indexing
errorfolder=$rootfolder/errorfolders/bam_indexing

mkdir -p $logfolder
mkdir -p $errorfolder

numbams=($(wc -l $samplelist))

bsub -J "samindexArray[1-$numbams]" \
     -R"select[mem>300] rusage[mem=300]" \
     -M300 \
     -o $logfolder/test_output_%J.%I.txt \
     -e $errorfolder/test_error_%J.%I.txt \
     ' '${scriptsfolder}'/indexwrapper.sh '${bamfilefolder}' '${samplelist}' ${LSB_JOBINDEX}'

