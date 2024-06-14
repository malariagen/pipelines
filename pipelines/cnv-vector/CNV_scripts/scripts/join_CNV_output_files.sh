study=$1

coveragefolder=/lustre/scratch118/malaria/team112/personal/el10/$study/coverage
outputfolder=/lustre/scratch118/malaria/team112/personal/el10/$study/coverage_CNVs
combined_full_CNV_table_gambcolu=$outputfolder/full_coverage_CNV_table_gambcolu.csv
combined_full_raw_CNV_table_gambcolu=$outputfolder/full_raw_CNV_table_gambcolu.csv
combined_full_CNV_table_arabiensis=$outputfolder/full_coverage_CNV_table_arabiensis.csv
combined_full_raw_CNV_table_arabiensis=$outputfolder/full_raw_CNV_table_arabiensis.csv

mkdir -p $outputfolder

cat $coveragefolder/2L/HMM_output/gambcolu_CNV/full_coverage_CNV_table_2L.csv > $combined_full_CNV_table_gambcolu
tail -n +2 $coveragefolder/2R/HMM_output/gambcolu_CNV/full_coverage_CNV_table_2R.csv >> $combined_full_CNV_table_gambcolu
tail -n +2 $coveragefolder/3L/HMM_output/gambcolu_CNV/full_coverage_CNV_table_3L.csv >> $combined_full_CNV_table_gambcolu
tail -n +2 $coveragefolder/3R/HMM_output/gambcolu_CNV/full_coverage_CNV_table_3R.csv >> $combined_full_CNV_table_gambcolu
tail -n +2 $coveragefolder/X/HMM_output/gambcolu_CNV/full_coverage_CNV_table_X.csv >> $combined_full_CNV_table_gambcolu

cat $coveragefolder/2L/HMM_output/gambcolu_CNV/full_raw_CNV_table_2L.csv > $combined_full_raw_CNV_table_gambcolu
tail -n +2 $coveragefolder/2R/HMM_output/gambcolu_CNV/full_raw_CNV_table_2R.csv >> $combined_full_raw_CNV_table_gambcolu
tail -n +2 $coveragefolder/3L/HMM_output/gambcolu_CNV/full_raw_CNV_table_3L.csv >> $combined_full_raw_CNV_table_gambcolu
tail -n +2 $coveragefolder/3R/HMM_output/gambcolu_CNV/full_raw_CNV_table_3R.csv >> $combined_full_raw_CNV_table_gambcolu
tail -n +2 $coveragefolder/X/HMM_output/gambcolu_CNV/full_raw_CNV_table_X.csv >> $combined_full_raw_CNV_table_gambcolu

cat $coveragefolder/2L/HMM_output/arabiensis_CNV/full_coverage_CNV_table_2L.csv > $combined_full_CNV_table_arabiensis
tail -n +2 $coveragefolder/2R/HMM_output/arabiensis_CNV/full_coverage_CNV_table_2R.csv >> $combined_full_CNV_table_arabiensis
tail -n +2 $coveragefolder/3L/HMM_output/arabiensis_CNV/full_coverage_CNV_table_3L.csv >> $combined_full_CNV_table_arabiensis
tail -n +2 $coveragefolder/3R/HMM_output/arabiensis_CNV/full_coverage_CNV_table_3R.csv >> $combined_full_CNV_table_arabiensis
tail -n +2 $coveragefolder/X/HMM_output/arabiensis_CNV/full_coverage_CNV_table_X.csv >> $combined_full_CNV_table_arabiensis

cat $coveragefolder/2L/HMM_output/arabiensis_CNV/full_raw_CNV_table_2L.csv > $combined_full_raw_CNV_table_arabiensis
tail -n +2 $coveragefolder/2R/HMM_output/arabiensis_CNV/full_raw_CNV_table_2R.csv >> $combined_full_raw_CNV_table_arabiensis
tail -n +2 $coveragefolder/3L/HMM_output/arabiensis_CNV/full_raw_CNV_table_3L.csv >> $combined_full_raw_CNV_table_arabiensis
tail -n +2 $coveragefolder/3R/HMM_output/arabiensis_CNV/full_raw_CNV_table_3R.csv >> $combined_full_raw_CNV_table_arabiensis
tail -n +2 $coveragefolder/X/HMM_output/arabiensis_CNV/full_raw_CNV_table_X.csv >> $combined_full_raw_CNV_table_arabiensis

