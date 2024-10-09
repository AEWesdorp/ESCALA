#!/bin/bash

# Define variables
output_dir="../../output/sizeSelection/tmp"
output_file="../../output/sizeSelection/BAL_sizeSelection_output.txt"
samples_file="../../scripts/01_PreProcessing/config/samples.txt"
taxID_list=(1131492 5052 746128)
size_ranges=(50 75 100 125 150 1000)

echo -e "sample_id\tdatabase\tCT\ttaxID\tFilter\tQC_number_reads\ttaxID_number_reads" > "$output_file" 
# Function to check if a file exists
check_file_exists() {
	if [[ ! -f "$1" ]]; then
		echo "Error: File not found: $1"
		#exit 1
	fi
}

# Function to filter Kraken2 output by taxID
filter_by_taxID() {
	local taxID="$1"
	local sample_id="$2"
	local database="$3"
	local CT="$4"
	local in_k2_rep="$5"
	local in_k2_out="$6"
	local tax_names_incl_daughters="${output_dir}/${sample_id}_taxid_names_incl_daughters_${taxID}.txt"
	local out_k2_filtered="${output_dir}/${sample_id}_${database}_conf${CT}_filtered_${taxID}.output"

	echo "Processing taxID: $taxID"

	./../resources/extract_kraken2_taxids_edit.py -r "$in_k2_rep" -t "$taxID" --include-children --output_file "$tax_names_incl_daughters"
	check_file_exists "$tax_names_incl_daughters"

	awk 'NR == FNR { keywords[$1]=1; next; } { if ($3 in keywords) print $0; }' "$tax_names_incl_daughters" "$in_k2_out" > "$out_k2_filtered"
	check_file_exists "$out_k2_filtered"

	local taxID_number_reads=$(wc -l < "$out_k2_filtered")
	echo -e "${sample_id}\t${database}\tconf${CT}\t${taxID}\tnoFilter\t${QC_number_reads}\t${taxID_number_reads}" >> "$output_file"
}

# Function to filter reads by size range
filter_by_size_range() {
	local max_value="$1"
	local sample_id="$2"
	local merged_fq="$3"
	local out_filtered_fq="$4"

	echo "Filtering reads for size range <= $max_value"
	awk -v max_value="$max_value" '($5 == "P" || $5 == "I") && ($4 >= 35 && $4 <= max_value) {print $0}' "$merged_fq" | sort -k1,1 --temporary-directory="${output_dir}" --parallel=8 --buffer-size=25% --compress-program=gzip --batch-size=1000 > "$out_filtered_fq"
	check_file_exists "$out_filtered_fq"
}

# Main processing loop for each sample
for sample_id in $(cut -f 1 "$samples_file" | grep -v "sample_name" | grep "B" ); do
	echo "Processing sample: $sample_id"

	# Define input/output files for the sample
	database="EPRSc2"
	CT="0.4"
	in_fq1="../../output/cfspi/samples/results/clean_fastq/${sample_id}_R1_trimmed_truncated.fastq.gz"
	in_fq2="../../output/cfspi/samples/results/clean_fastq/${sample_id}_R2_trimmed_truncated.fastq.gz"
	in_k2_out="../../output/cfspi/samples/results/kraken2_output/after_host_mapping/${sample_id}_${database}_conf${CT}.output"
	in_k2_rep="../../output/cfspi/samples/results/kraken2_report/after_host_mapping/${sample_id}_${database}_conf${CT}.report"
	out_merged_fq="${output_dir}/${sample_id}_merged.fastq"
	out_merged_fq_size_range="${output_dir}/${sample_id}_merged_size_range_${taxID}.fastq"

	# Merging files (uncomment this if `bbmerge.sh` is available)
	#bbmerge.sh in="$in_fq1" in2="$in_fq2" outinsert="$out_merged_fq" threads=16

	# Check merged FASTQ file exists and count number of reads
	check_file_exists "$out_merged_fq"
	QC_number_reads=$(wc -l < "$out_merged_fq")
	echo "Number of reads in merged FASTQ: $QC_number_reads"

	# Process each taxID
	for taxID in "${taxID_list[@]}"; do
		filter_by_taxID "$taxID" "$sample_id" "$database" "$CT" "$in_k2_rep" "$in_k2_out"
	done

	# Filter by size range
	for max_value in "${size_ranges[@]}"; do
		filter_by_size_range "$max_value" "$sample_id" "$out_merged_fq" "$out_merged_fq_size_range"

		for taxID in "${taxID_list[@]}"; do
			echo "Processing size range and taxID: $max_value, $taxID"
			out_k2_filtered="${output_dir}/${sample_id}_${database}_conf${CT}_filtered_${taxID}.output"
			sorted_read_names=$(awk '{print $2}' "$out_k2_filtered" | sort -u)

			echo "Number of reads in range <= $max_value"
			QC_number_reads=$(wc -l < "$out_merged_fq_size_range")

			echo "Number of reads in range <= $max_value, of taxID $taxID"
			#taxID_number_reads=$(grep -F -f <(echo "$sorted_read_names") "$out_merged_fq_size_range" | wc -l)
			temp_sorted_read_names=$(mktemp)
			echo "$sorted_read_names" > "$temp_sorted_read_names"
			taxID_number_reads=$(grep -F -f "$temp_sorted_read_names" "$out_merged_fq_size_range" | wc -l)
			rm "$temp_sorted_read_names"
			echo -e "${sample_id}\t${database}\tconf${CT}\t${taxID}\t35to${max_value}\t${QC_number_reads}\t${taxID_number_reads}" >> "$output_file"
		done
	done
done
