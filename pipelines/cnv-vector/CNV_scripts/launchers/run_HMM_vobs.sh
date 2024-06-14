study=$1

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
logfolder=$coveragefolder/logfolders/HMM
errorfolder=$coveragefolder/errorfolders/HMM

Agam_GC_file=$rootfolder/phase3_data/tables/Agam_genome_GC_content.csv

allchrom=(2L 2R 3L 3R X)

# Create the outputfolders if necessary
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
	echo "Running HMM for samples from $sample_manifest on `date`." >> $coveragefolder/HMM_added_samples.log

	# Get some of the necessary input filenames
	coverage_by_GC_file=median_coverage_by_GC_masked_09_05_${species}.csv
	coverage_variance_file=coverage_variance_masked_09_05_${species}.csv
	mapq_prop_file=$rootfolder/phase3_cnv/coverage/mapq_proportions_allchrom_${species}.csv

	# Run the HMM script on each species on each chromosome
	for chrom in ${allchrom[@]}
	do
		mkdir -p $coveragefolder/$chrom/HMM_output
		mkdir -p $coveragefolder/$chrom/HMM_logs_$species

		bsub -o $logfolder/HMM_output_${species}_%J.txt \
			 -e $errorfolder/HMM_error_${species}_%J.txt \
			 ${scriptsfolder}/coverage_HMM_vobs.sh $coveragefolder \
											       $sample_manifest \
											       $chrom \
											       $species \
											       $Agam_GC_file \
											       $coverage_by_GC_file \
											       $coverage_variance_file \
											       $mapq_prop_file
	done
done
