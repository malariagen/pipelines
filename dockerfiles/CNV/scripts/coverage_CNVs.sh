coveragefolder=$1
manifest=$2
chrom=$3
output_name=$4
coverage_variance_file=$5
ncores=$6
metadata=$7

gene_coordinates_file=~/personal/phase3_data/tables/gene_regions.csv
detox_genes_file=~/personal/phase3_data/tables/detox_genes.txt
workingfolder=$coveragefolder/$chrom/HMM_output
outputfolder=$workingfolder/$output_name
scriptsfolder=~/scripts/CNV_scripts/scripts

R --version

R --slave -f $scriptsfolder/CNV_analysis.r --args $chrom \
                                                  $manifest \
                                                  $coverage_variance_file \
                                                  $gene_coordinates_file \
                                                  $detox_genes_file \
                                                  $workingfolder \
                                                  $ncores \
                                                  $outputfolder \
                                                  $metadata \
												  > $coveragefolder/$chrom/CNV_analysis_logs/CNV_analysis_${output_name}.log 2>&1

