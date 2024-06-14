study=$1

# turn on extended globbing
shopt -s extglob

coveragefolder=/lustre/scratch118/malaria/team112/personal/el10/$study/coverage
coverage_variance_files=($(ls $coveragefolder/coverage_variance_masked_09_05_!(all).csv))
combined_coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv

if [ ${#coverage_variance_files[@]} -eq 1 ]; then
	cp $coverage_variance_files $combined_coverage_variance_file
else
	(cat ${coverage_variance_files[0]}; sed -s 1d ${coverage_variance_files[@]:1}) > $combined_coverage_variance_file
fi
