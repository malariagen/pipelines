samplelist=$1
species_id_file=$2
metadata_file=$3
study=$4

scriptsfolder=~/scripts/CNV_scripts/scripts
workingfolder=/lustre/scratch118/malaria/team112/personal/el10/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
coverage_variance_file=$workingfolder/coverage/coverage_variance_masked_09_05_all.csv
ncores=5
logfolder=$diagnostic_reads_folder/logfolders/target_regions_analysis
errorfolder=$diagnostic_reads_folder/errorfolders/target_regions_analysis

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_analysis_output_%J.txt \
     -e $errorfolder/target_regions_analysis_error_%J.txt \
     -q long \
     -n $ncores \
     -R"select[mem>2000] rusage[mem=2000] span[hosts=1]" \
     -M2000 \
     $scriptsfolder/target_regions_analysis.sh $workingfolder \
                                               $samplelist \
                                               $species_id_file \
                                               $coverage_variance_file \
                                               $metadata_file \
                                               $ncores
