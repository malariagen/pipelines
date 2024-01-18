metadata_file=$1
study=$2

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/phase3_data/tables/gene_regions.csv
outputfolder=$rootfolder/$study/modal_CNVs
logfolder=$coveragefolder/logfolders/modal_CNV
errorfolder=$coveragefolder/errorfolders/modal_CNV

mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder

all_manifests=($(ls $manifestfolder/sample_manifest_*.txt))

for sample_manifest in ${all_manifests[@]}
do
	# Pull out the species from the manifest names
	s=${sample_manifest#$manifestfolder\/sample_manifest_}
	species=${s%.txt}
	if [ $species == "known_species" ]; then
		continue
	fi

bsub -o $logfolder/modal_CNVs_${species}_%J.txt \
     -e $errorfolder/modal_CNVs_${species}_%J.txt \
     -q long \
     -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
     -M500 \
     ${scriptsfolder}/modal_CNVs.sh $coveragefolder \
	                                $sample_manifest \
		                            $species \
	                                $coverage_variance_file \
	                                $metadata_file \
	                                $outputfolder

done
