study=$1

rootfolder=/lustre/scratch118/malaria/team112/personal/el10
logfolder=$rootfolder/$study/coverage/logfolders/stats_by_species
errorfolder=$rootfolder/$study/coverage/errorfolders/stats_by_species
workingfolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
datafolder=$rootfolder/phase3_data
scriptsfolder=~/scripts/CNV_scripts/scripts

mkdir -p $logfolder
mkdir -p $errorfolder

# Find out what species are present in this sample set
s=($(ls $manifestfolder/sample_manifest_*.txt))
ss=(${s[@]#$manifestfolder\/sample_manifest_})
specieslist=(${ss[@]%.txt})

for species in ${specieslist[@]}
do
	bsub -o $logfolder/${species}_coverage_stats_log.txt \
		 -e $errorfolder/${species}_coverage_stat_error.txt \
		 ${scriptsfolder}/get_coverage_stats_by_sample_set_vobs.sh $workingfolder \
															       $manifestfolder/sample_manifest_${species}.txt \
															       $datafolder/windowed_accessibility/mean_accessibility_${species}.csv \
															       $datafolder/mapq_proportions/mapq_proportions_allchrom_${species}.csv \
															       $datafolder/tables/Agam_genome_GC_content.csv \
															       $species 
done

