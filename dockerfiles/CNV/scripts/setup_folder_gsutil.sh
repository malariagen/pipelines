# This script prepares a folder in which a given sample set will be run
folder_location=$1
gsutil_path=$2
speciescalls_folder=$3
sampleset=$4

# turn on extended globbing
shopt -s extglob

# Create the folder
mkdir -p $folder_location/data
# Download the metadata to the folder
sample_metadata=$gsutil_path/metadata/general/$sampleset/samples.meta.csv
gsutil cp $sample_metadata $folder_location/data
metadata_filename="$(basename -- $sample_metadata)"
# Create a sample manifest from the metadata
sed -e '1d' -e '2,$s/,.*//' $folder_location/data/$metadata_filename > $folder_location/data/sample_manifest.txt

# Download the table linking sample names with bam files
sample_bampaths=$gsutil_path/metadata/general/$sampleset/wgs_snp_data.csv
gsutil cp $sample_bampaths $folder_location/data
bampaths_filename="$(basename -- $sample_bampaths)"

# Download the species calls
sample_speciescalls=$gsutil_path/metadata/$speciescalls_folder/$sampleset/samples.species_aim.csv
gsutil cp $sample_speciescalls $folder_location/data
speciescalls_filename="$(basename -- $sample_speciescalls)"
# Some of the scripts will expect tab-delimited metadata and species calls
sed -e "s/,/\t/g" $folder_location/data/$metadata_filename > $folder_location/data/sample_metadata_tabs.tsv
sed -e "s/,/\t/g" $folder_location/data/$speciescalls_filename > $folder_location/data/sample_speciescalls_tabs.tsv

# To create species-specific manifests, we want to pull out the gambcolu/arabiensis
# species call data. We could do this by column number, but I don't know whether the
# columns numbers will be consistent. So we do it by column names. The first step for
# this is to cut down the species calls file to just two columns (sample name and 
# species call). Then we can be sure that the species call will be column 2.
species_manifest (){ 
	oldIFS=$IFS
	IFS=","
	read names
	while read $names; do
		echo -e ${sample_id}\\t${species_gambcolu_arabiensis/_/}
	done
	IFS=$oldIFS
}

# So now we apply that function, and then split the sample names into two files based
# on column 2 (which is guaranteed to be the species calls, as long as the column names
# are always the same). 
cat $folder_location/data/$speciescalls_filename | species_manifest | awk -v output_folder="$folder_location/data" '{print $1 > output_folder"/sample_manifest_"$2".txt"}'
cat $folder_location/data/sample_manifest_!(intermediate).txt > $folder_location/data/sample_manifest_known_species.txt

# We want to pull out the paths to the bamfile for each sample
bamfiles (){ 
	oldIFS=$IFS
	IFS=","
	read names
	while read $names; do
		echo -e ${sample_id}\\t${alignments_bam}
	done
	IFS=$oldIFS
}

cat $folder_location/data/$bampaths_filename | bamfiles > $folder_location/data/bampaths.csv

# Create a folder for the bamfiles
mkdir -p $folder_location/bamlinks
# Download all the required bams
cd $folder_location/bamlinks
while IFS=$'\t' read sample bampath
do 
	curl ${bampath} > ${sample}.bam
	curl ${bampath}.bai > ${sample}.bam.bai
done < ../data/bampaths.csv




