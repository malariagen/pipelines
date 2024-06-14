bamfilefolder=$1
manifest=$2
samplenum=$3 

samplename=($(head -n$samplenum $manifest | tail -n1))
bamfile=${bamfilefolder}/${samplename}.bam

echo $bamfile

module load -s samtools/1.9

samtools index $bamfile
