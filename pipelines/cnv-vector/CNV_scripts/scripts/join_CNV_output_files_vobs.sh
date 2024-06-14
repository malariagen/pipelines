study=$1

rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
outputfolder=$rootfolder/$study/coverage_CNVs

mkdir -p $outputfolder

all_CNV_folders=($(ls -d $coveragefolder/2L/HMM_output/*_CNV))
ss=(${all_CNV_folders[@]#$coveragefolder/2L/HMM_output/})
specieslist=(${ss[@]%_CNV})

for species in ${specieslist[@]}
do
	combined_full_CNV_table=$outputfolder/full_coverage_CNV_table_${species}.csv
	combined_full_raw_CNV_table=$outputfolder/full_raw_CNV_table_${species}.csv

	cat $coveragefolder/2L/HMM_output/${species}_CNV/full_coverage_CNV_table_2L.csv > $combined_full_CNV_table
	tail -n +2 $coveragefolder/2R/HMM_output/${species}_CNV/full_coverage_CNV_table_2R.csv >> $combined_full_CNV_table
	tail -n +2 $coveragefolder/3L/HMM_output/${species}_CNV/full_coverage_CNV_table_3L.csv >> $combined_full_CNV_table
	tail -n +2 $coveragefolder/3R/HMM_output/${species}_CNV/full_coverage_CNV_table_3R.csv >> $combined_full_CNV_table
	tail -n +2 $coveragefolder/X/HMM_output/${species}_CNV/full_coverage_CNV_table_X.csv >> $combined_full_CNV_table

	cat $coveragefolder/2L/HMM_output/${species}_CNV/full_raw_CNV_table_2L.csv > $combined_full_raw_CNV_table
	tail -n +2 $coveragefolder/2R/HMM_output/${species}_CNV/full_raw_CNV_table_2R.csv >> $combined_full_raw_CNV_table
	tail -n +2 $coveragefolder/3L/HMM_output/${species}_CNV/full_raw_CNV_table_3L.csv >> $combined_full_raw_CNV_table
	tail -n +2 $coveragefolder/3R/HMM_output/${species}_CNV/full_raw_CNV_table_3R.csv >> $combined_full_raw_CNV_table
	tail -n +2 $coveragefolder/X/HMM_output/${species}_CNV/full_raw_CNV_table_X.csv >> $combined_full_raw_CNV_table
done
