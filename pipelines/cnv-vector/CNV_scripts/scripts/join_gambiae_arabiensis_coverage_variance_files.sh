study=$1

coveragefolder=/lustre/scratch118/malaria/team112/personal/el10/$study/coverage
arabiensis_coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_arabiensis.csv
gambcolu_coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_gambcolu.csv
combined_coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv

(cat $arabiensis_coverage_variance_file; tail -n +2 $gambcolu_coverage_variance_file) > $combined_coverage_variance_file

