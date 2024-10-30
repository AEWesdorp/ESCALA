#!/bin/bash
#Check if input file is provided
if [ "$#" -ne 1 ]; then
	 echo "Usage: $0 genome.fa"
	  exit 1
  fi

  # Input genome file
  filename=$1

  ###nonMT
  output_nonMT_name="../../../output/nucleotideFreqGenomes/$(basename "$filename" | sed 's/\.[^.]*$//')_dinuclFreq_nonMT.txt"
  echo "Generating... "
  echo $output_nonMT_name
  bioawk -c fastx '
  /^>/ { 
  if ($name == "chrM") next;  # Skip sequences named "chrM"
  }
  {
	  # Convert sequence to uppercase
	  newseq = toupper($seq);

	  # Loop through the sequence to count dinucleotide frequencies
	  for (i = 1; i <= length(newseq) - 1; i++) {
		  dinuc = substr(newseq, i, 2);
		  dinuc_freq[dinuc]++;
	  }
  }

  END {
  # Print dinucleotide frequencies after processing all sequences
  for (dinuc in dinuc_freq) {
	  print dinuc, dinuc_freq[dinuc];
  }
  }' $filename > $output_nonMT_name

  ###MT 
  output_MT_name="../../../output/nucleotideFreqGenomes/$(basename "$filename" | sed 's/\.[^.]*$//')_dinuclFreq_MT.txt"
  echo "Generating... "
  echo $output_MT_name
  bioawk -c fastx '
  /^>/ { 
  if ($name != "chrM") next;  # Skip sequences NOT named "chrM"
  }
  {
	  # Convert sequence to uppercase
	  newseq = toupper($seq);

	  # Loop through the sequence to count dinucleotide frequencies
	  for (i = 1; i <= length(newseq) - 1; i++) {
		  dinuc = substr(newseq, i, 2);
		  dinuc_freq[dinuc]++;
	  }
  }

  END {
  # Print dinucleotide frequencies after processing all sequences
  for (dinuc in dinuc_freq) {
	  print dinuc, dinuc_freq[dinuc];
  }
  }' $filename > $output_MT_name

