# This script prepares a folder in which a given sample set will be run
folder_location=$1
sample_metadata=$2
sample_speciescalls=$3
bamfile_folder=$4
bamindex_folder=$5

# turn on extended globbing
shopt -s extglob

# Create the folder
mkdir -p $folder_location/data
# Copy the metadata to the folder
cp $sample_metadata $folder_location/data
metadata_filename="$(basename -- $sample_metadata)"
# Create a sample manifest from the metadata
sed -e '1d' -e '2,$s/,.*//' $folder_location/data/$metadata_filename > $folder_location/data/sample_manifest.txt
# Copy the species calls
cp $sample_speciescalls $folder_location/data
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
# are always the same. 
cat $sample_speciescalls | species_manifest | awk -v output_folder="$folder_location/data" '{print $1 > output_folder"/sample_manifest_"$2".txt"}'
cat $folder_location/data/sample_manifest_!(known_species).txt > $folder_location/data/sample_manifest_known_species.txt

# Create a folder for the bamlinks
mkdir -p $folder_location/bamlinks
# Make symlinks to all the required bams
cd $folder_location/bamlinks
while read p
do 
	ln -s $bamfile_folder/${p}.bam
	ln -s $bamindex_folder/${p}.bam.bai
done < $folder_location/data/sample_manifest.txt




