arabiensis_samplelist=$1
gambcolu_samplelist=$2
study=$3

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
logfolder=$coveragefolder/logfolders/HMM
errorfolder=$coveragefolder/errorfolders/HMM

Agam_GC_file=$rootfolder/phase3_data/tables/Agam_genome_GC_content.csv
arabiensis_coverage_by_GC_file=median_coverage_by_GC_masked_09_05_arabiensis.csv
arabiensis_coverage_variance_file=coverage_variance_masked_09_05_arabiensis.csv
arabiensis_mapq_prop_file=mapq_proportions_allchrom_arabiensis.csv
gambcolu_coverage_by_GC_file=median_coverage_by_GC_masked_09_05_gambcolu.csv
gambcolu_coverage_variance_file=coverage_variance_masked_09_05_gambcolu.csv
gambcolu_mapq_prop_file=mapq_proportions_allchrom_gambcolu.csv

allchrom=(2L 2R 3L 3R X)

echo "Running HMM for samples from $arabiensis_samplelist and $gamcolu_samplelist on `date`." >> $coveragefolder/HMM_added_samples.log

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder
for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/HMM_logs_arabiensis
	mkdir -p $coveragefolder/$chrom/HMM_logs_gambcolu
	mkdir -p $coveragefolder/$chrom/HMM_output

	# arabiensis
	bsub -o $logfolder/HMM_output_arabiensis_%J.txt \
		 -e $errorfolder/HMM_error_arabiensis_%J.txt \
		 ${scriptsfolder}/coverage_HMM.sh $coveragefolder \
	                                      $arabiensis_samplelist \
		                                  $chrom \
		                                  arabiensis \
	                                      $Agam_GC_file \
	                                      $arabiensis_coverage_by_GC_file \
	                                      $arabiensis_coverage_variance_file \
	                                      $arabiensis_mapq_prop_file

	# gambiae coluzzii
	bsub -o $logfolder/HMM_output_gambcolu_%J.txt \
		 -e $errorfolder/HMM_error_gambcolu_%J.txt \
		 ${scriptsfolder}/coverage_HMM.sh $coveragefolder \
	                                      $gambcolu_samplelist \
		                                  $chrom \
		                                  gambcolu \
	                                      $Agam_GC_file \
	                                      $gambcolu_coverage_by_GC_file \
	                                      $gambcolu_coverage_variance_file \
	                                      $gambcolu_mapq_prop_file
done
