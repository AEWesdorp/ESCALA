import sys
import pysam

def parse_bam(bam_file):
    paired_reads = {}

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate through each read in the BAM file
        for read in bam:
            # Check if the read is mapped and its mate is also mapped
            if not read.is_unmapped and not read.mate_is_unmapped:
                # If it's R1, store its start position
                if read.is_read1:
                    read_name = read.query_name
                    if read_name not in paired_reads:
                        paired_reads[read_name] = {'R1': read.reference_start}
                    else:
                        paired_reads[read_name]['R1'] = read.reference_start
                # If it's R2, calculate TLEN and print the information
                elif read.is_read2:
                    read_name = read.query_name
                    if read_name in paired_reads:
                        paired_reads[read_name]['R2'] = read.reference_start
                        paired_reads[read_name]['TLEN'] = abs(read.template_length)
                        print(f"Start of R1: {paired_reads[read_name]['R1']}, Start of R2: {paired_reads[read_name]['R2']}, TLEN: {paired_reads[read_name]['TLEN']}")

if __name__ == "__main__":
    # Check if the user provided a BAM file as input
    if len(sys.argv) != 2:
        print("Usage: python parse_bam.py input.bam")
        sys.exit(1)

    # Get the BAM file path from the command line argument
    bam_file = sys.argv[1]

    # Call the parse_bam function with the provided BAM file
    parse_bam(bam_file)
