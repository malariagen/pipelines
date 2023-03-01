workingfolder=$1
manifest=$2
species_id_file=$3
coverage_variance_file=$4
metadata=$5
ncores=$6

coveragefolder=$workingfolder/coverage
diagnostic_reads_folder=$workingfolder/diagnostic_reads
gene_coordinates_file=~/personal/phase3_data/tables/gene_regions.csv
scriptsfolder=~/scripts/CNV_scripts/scripts
plotting_functions_file=$scriptsfolder/plotting_functions.r

cd $workingfolder

mkdir -p target_regions_analysis

export R_LIBS_USER="~/R-modules:$R_LIBS_USER"

R-3.6.1 --version

R-3.6.1 --slave -f $scriptsfolder/target_regions_analysis.r --args $manifest \
                                                                   $gene_coordinates_file \
                                                                   $metadata \
                                                                   $species_id_file \
                                                                   $coverage_variance_file \
                                                                   $coveragefolder \
                                                                   $diagnostic_reads_folder \
                                                                   $plotting_functions_file \
                                                                   $ncores \
                                                                   > target_regions_analysis/target_regions_analysis.log 2>&1


