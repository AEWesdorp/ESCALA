#!/bin/bash

# Specify the directory containing the .fastq.gz files
DIRECTORY="../../../output/cfspi/samples/results/clean_fastq/"  # Change this to your folder path

# Loop through each .fastq.gz file in the directory

# Loop through each .fastq.gz file in the directory
for file in "$DIRECTORY"/*.fastq.gz; do
	# Extract the base filename without the extension
	base_name=$(basename "$file" "_trimmed_truncated.fastq.gz")

	# Construct the new filename
	new_file_name="${base_name}_endMotif_count.txt"

	echo "Processing file: $(basename "$file")"
	echo "New file name: $new_file_name"

	# Use zcat to decompress and process the file, and redirect output to the new file
	zcat "$file" | awk 'NR % 4 == 2 {print substr($0, 1, 3)}' | \
		sort --parallel=32 -T ../../../output/nucleotideFreq_QCreads/temp/ | uniq -c | \
		awk 'BEGIN {{ OFS=","; print "EndMotif,Count" }} {{print $2,$1}}' >  "../../../output/nucleotideFreq_QCreads/${new_file_name}"
done

