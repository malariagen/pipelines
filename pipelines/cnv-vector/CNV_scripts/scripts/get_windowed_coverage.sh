bamfilefolder=$1
manifest=$2
samplenum=$3 
outputfolder=$4
scriptsfolder=~/scripts/CNV_scripts/scripts
allchrom=(2L 2R 3L 3R X)

samplename=($(head -n$samplenum $manifest | tail -n1))
bamfile=${bamfilefolder}/${samplename}.bam

echo $bamfile

source activate cnv37 

for chrom in ${allchrom[@]}
do
	python ${scriptsfolder}/counts_for_HMM.py \
	       $bamfile \
	       $chrom \
	       300 300 10 \
	       ${outputfolder}/${chrom}/counts_for_HMM_${samplename}_${chrom}_output.csv \
	       > ${outputfolder}/${chrom}/coveragelogs/counts_for_HMM_${samplename}_${chrom}.log 2>&1
done


