metadata_file=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/phase3_data/tables/gene_regions.csv
detox_genes_file=$rootfolder/phase3_data/tables/detox_genes.txt
ncores=2
logfolder=$coveragefolder/logfolders/CNV_analysis
errorfolder=$coveragefolder/errorfolders/CNV_analysis

allchrom=(2L 2R 3L 3R X)

mkdir -p $logfolder
mkdir -p $errorfolder

all_manifests=($(ls $manifestfolder/sample_manifest_*.txt))

for sample_manifest in ${all_manifests[@]}
do
	# Pull out the species from the manifest names
	s=${sample_manifest#$manifestfolder\/sample_manifest_}
	species=${s%.txt}
	if [ $species == "known_species" ]; then
		continue
	fi

	output_id=${species}_CNV

	for chrom in ${allchrom[@]}
	do
		mkdir -p $coveragefolder/$chrom/CNV_analysis_logs

		bsub -o $logfolder/CNV_analysis_output_${chrom}_${species}_%J.txt \
			 -e $errorfolder/CNV_analysis_error_${chrom}_${species}_%J.txt \
			 -n $ncores \
			 -q long \
			 -R"select[mem>200] rusage[mem=200] span[hosts=1]" \
			 -M200 \
			 ${scriptsfolder}/coverage_CNVs_vobs.sh $coveragefolder \
											        $sample_manifest \
											        $chrom \
											        $output_id \
											        $coverage_variance_file \
											        $ncores \
											        $metadata_file \
											        $gene_coordinates_file \
											        $detox_genes_file
	done
done
