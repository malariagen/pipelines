release=$1
study_id=$2

workingfolder=~/personal/${release}_${study_id}
gsutil_path=gs://vo_agam_release/${release}
logfolder=$workingfolder/folder_setup_log 
mkdir -p $logfolder 
speciescalls_folder=species_calls_20200422 
~/scripts/CNV_scripts/scripts/setup_folder_gsutil.sh $workingfolder $gsutil_path $speciescalls_folder $study_id > $logfolder/setup_folder_gsutil.log 2>&1 & 
