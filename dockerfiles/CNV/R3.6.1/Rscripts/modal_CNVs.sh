coveragefolder=$1
manifest=$2
output_name=$3
coverage_variance_file=$4
metadata=$5
outputfolder=$6

gene_coordinates_file=~/personal/phase3_data/tables/gene_regions.csv
scriptsfolder=~/scripts/CNV_scripts/scripts

R --version

R --slave -f $scriptsfolder/modal_cnv.r --args $manifest \
                                               $gene_coordinates_file \
                                               $metadata \
                                               $output_name \
                                               $coverage_variance_file \
                                               $coveragefolder \
                                               $outputfolder \
                                               > $outputfolder/modal_CNV_${output_name}.log 2>&1

