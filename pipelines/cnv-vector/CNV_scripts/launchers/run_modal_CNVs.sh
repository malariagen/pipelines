gambcolu_samplelist=$1
arabiensis_samplelist=$2
metadata_file=$3
study=$4

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
coveragefolder=$rootfolder/$study/coverage
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
outputfolder=$rootfolder/$study/modal_CNVs
logfolder=$coveragefolder/logfolders/modal_CNV
errorfolder=$coveragefolder/errorfolders/modal_CNV

mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder

bsub -o $logfolder/modal_CNVs_gambcolu_%J.txt \
     -e $errorfolder/modal_CNVs_gambcolu_%J.txt \
     -q long \
     -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
     -M500 \
     ${scriptsfolder}/modal_CNVs.sh $coveragefolder \
	                                $gambcolu_samplelist \
		                            gambcolu \
	                                $coverage_variance_file \
	                                $metadata_file \
	                                $outputfolder

bsub -o $logfolder/modal_CNVs_arabiensis_%J.txt \
     -e $errorfolder/modal_CNVs_arabiensis_%J.txt \
     -q long \
     -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
     -M500 \
     ${scriptsfolder}/modal_CNVs.sh $coveragefolder \
	                                $arabiensis_samplelist \
		                            arabiensis \
	                                $coverage_variance_file \
	                                $metadata_file \
	                                $outputfolder

