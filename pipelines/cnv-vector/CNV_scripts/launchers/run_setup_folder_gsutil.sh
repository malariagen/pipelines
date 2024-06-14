release=$1
speciescalls_folder=$2
sampleset=$3

scriptsfolder=~/scripts/CNV_scripts/scripts
study=${release}_${sampleset}
gsutil_path=gs://vo_agam_release/$release
workingfolder=/lustre/scratch118/malaria/team112/personal/el10/$study
logfolder=$workingfolder/folder_setup_log

mkdir -p $logfolder

bsub -o $logfolder/folder_setup_log_%J.txt \
     -e $logfolder/folder_setup_error_%J.txt \
     -q long \
     -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
     -M500 \
     $scriptsfolder/setup_folder_gsutil.sh $workingfolder \
                                           $gsutil_path \
                                           $speciescalls_folder \
                                           $sampleset
