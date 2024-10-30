#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
	 echo "Usage: $0 genome.fa"
	  exit 1
  fi

  # Input genome file
  filename=$1

  ###nonMT
  output_nonMT_name="../../../output/nucleotideFreqGenomes/$(basename "$filename" | sed 's/\.[^.]*$//')_nuclFreq_nonMT.txt"
  echo "Generating... "
  echo $output_nonMT_name
  bioawk -c fastx '
  /^>/ { 
  if ($name == "chrM") next;  # Skip sequences named "chrM"
  }
  {
	  # Convert sequence to uppercase
	  newseq = toupper($seq);

	  # Loop through the sequence to count nucleotide frequencies
	  for (i = 1; i <= length(newseq) - 1; i++) {
		  nuc = substr(newseq, i, 1);
		  nuc_freq[nuc]++;
	  }
  }

  END {
  # Print nucleotide frequencies after processing all sequences
  for (nuc in nuc_freq) {
	  print nuc, nuc_freq[nuc];
  }
  }' $filename > $output_nonMT_name

  ###MT 
  output_MT_name="../../../output/nucleotideFreqGenomes/$(basename "$filename" | sed 's/\.[^.]*$//')_nuclFreq_MT.txt"
  echo "Generating... "
  echo $output_MT_name
  bioawk -c fastx '
  /^>/ { 
  if ($name != "chrM") next;  # Skip sequences NOT named "chrM"
  }
  {
	  # Convert sequence to uppercase
	  newseq = toupper($seq);

	  # Loop through the sequence to count nucleotide frequencies
	  for (i = 1; i <= length(newseq) - 1; i++) {
		  nuc = substr(newseq, i, 1);
		  nuc_freq[nuc]++;
	  }
  }

  END {
  # Print nucleotide frequencies after processing all sequences
  for (nuc in nuc_freq) {
	  print nuc, nuc_freq[nuc];
  }
  }' $filename > $output_MT_name

