logfolder=/lustre/scratch118/malaria/team112/personal/el10/logfolders/accessibility_phase3
errorfolder=/lustre/scratch118/malaria/team112/personal/el10/errorfolders/accessibility_phase3
scriptsfolder=~/scripts/CNV_scripts/scripts

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/accessibility_phase3_log.txt \
     -e $errorfolder/accessibility_phase3_error.txt \
     ${scriptsfolder}/calculate_windowed_accessibility.sh 

