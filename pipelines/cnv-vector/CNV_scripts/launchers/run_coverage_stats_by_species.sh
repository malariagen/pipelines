study=$1
specieslist=${@:2}

rootfolder=/lustre/scratch118/malaria/team112/personal/el10
logfolder=$rootfolder/phase3_cnv/coverage/logfolders/stats_by_species
errorfolder=$rootfolder/phase3_cnv/coverage/errorfolders/stats_by_species
workingfolder=$rootfolder/$study/coverage
datafolder=$rootfolder/phase3_data
scriptsfolder=~/scripts/CNV_scripts/scripts

mkdir -p $logfolder
mkdir -p $errorfolder

for species in $specieslist
do
	bsub -o $logfolder/${species}_coverage_stats_log.txt \
		 -e $errorfolder/${species}_coverage_stat_error.txt \
		 ${scriptsfolder}/get_coverage_stats_by_sample_set.sh $workingfolder \
															  $datafolder/phase3_manifest_${species}.txt \
															  $datafolder/windowed_accessibility/mean_accessibility_${species}.csv \
															  $datafolder/tables/Agam_genome_GC_content.csv \
															  $species 
done

