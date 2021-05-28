coveragefolder=$1
manifest=$2
chrom=$3
species_id_file=$4
GC_content_file=$5
coverage_by_GC_file=$coveragefolder/$6
coverage_variance_file=$coveragefolder/$7
mapq_prop_file=$coveragefolder/$8

scriptsfolder=~/scripts/CNV_scripts/scripts

source activate cnv37 

python ${scriptsfolder}/HMM_process.py \
       $manifest \
       $chrom \
       $coveragefolder \
       $GC_content_file \
       $coverage_by_GC_file \
       $coverage_variance_file \
       $mapq_prop_file \
       0.5 \
       > ${coveragefolder}/${chrom}/HMM_logs_${species_id_file}/HMM_${chrom}.log 2>&1

