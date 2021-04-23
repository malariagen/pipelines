#!/usr/bin/env bash

module load samtools/1.9
file_size_threshold_bytes=2000000000
batch=$1

declare -A expected_rows_per_chromosome
expected_rows_per_chromosome=(['2R']=60132453 ['3R']=52226568 ['2L']=48525747 ['UNKN']=27274988 \
                              ['3L']=40758473 ['X']=23385349 ['Y_unplaced']=135155 ['Mt']=15363)

for path in $(s3cmd ls s3://$batch | grep '.*\.vcf.gz$' | awk -v threshold=${file_size_threshold_bytes} '$3 < threshold {print $4}'); do 
	filename=${path##*/}
	echo "Downloading $filename ..."
	s3cmd get ${path} > /dev/null 2>&1
	echo "Checking number of rows per chromosome in VCF ..."
	tmpfile=$(mktemp -p $PWD check.XXXXXX)
	bgzip -dc $filename | cut -f 1 | uniq -cd > "$tmpfile"
	for chr in "${!expected_rows_per_chromosome[@]}"; do
		rows=$(grep "$chr" "$tmpfile" | awk '{print $1}')
		if [ -z "$rows" ]; then
			echo "No rows found for $chr in ${filename}!"
		elif [ "$rows" != "${expected_rows_per_chromosome[$chr]}" ]; then
			echo "An unexpected number of rows ($rows) found for $chr in ${filename}. Expected ${expected_rows_per_chromosome[$chr]}"
		fi
	done
	echo
	echo "Summary of rows per chromosome for $filename:" 
	cat "$tmpfile"
	echo
	rm "$tmpfile"
done
