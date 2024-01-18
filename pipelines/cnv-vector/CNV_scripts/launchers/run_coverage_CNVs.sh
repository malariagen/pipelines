gambcolu_samplelist=$1
arabiensis_samplelist=$2
metadata_file=$3
study=$4

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
ncores=2
logfolder=$coveragefolder/logfolders/CNV_analysis
errorfolder=$coveragefolder/errorfolders/CNV_analysis

allchrom=(2L 2R 3L 3R X)

mkdir -p $logfolder
mkdir -p $errorfolder
for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/CNV_analysis_logs

	bsub -o $logfolder/CNV_analysis_output_${chrom}_gambcolu_%J.txt \
		 -e $errorfolder/CNV_analysis_error_${chrom}_gambcolu_%J.txt \
         -n $ncores \
         -q long \
         -R"select[mem>200] rusage[mem=200] span[hosts=1]" \
         -M200 \
		 ${scriptsfolder}/coverage_CNVs.sh $coveragefolder \
	                                       $gambcolu_samplelist \
		                                   $chrom \
		                                   gambcolu_CNV \
	                                       $coverage_variance_file \
	                                       $ncores \
	                                       $metadata_file

	bsub -o $logfolder/CNV_analysis_output_${chrom}_arabiensis_%J.txt \
		 -e $errorfolder/CNV_analysis_error_${chrom}_arabiensis_%J.txt \
         -n $ncores \
         -q long \
         -R"select[mem>200] rusage[mem=200] span[hosts=1]" \
         -M200 \
		 ${scriptsfolder}/coverage_CNVs.sh $coveragefolder \
	                                       $arabiensis_samplelist \
		                                   $chrom \
		                                   arabiensis_CNV \
	                                       $coverage_variance_file \
	                                       $ncores \
	                                       $metadata_file
done
