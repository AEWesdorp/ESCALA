## TaxID of Aspergillaceae (1131492), Aspergillus (5052), Aspergillus fumigatus (746128)
for taxid in $("1131492", "5052", "746128"); do 
for sample in $(cut -f 1  ../../01_PreProcessing/cfspi/config/samples.txt | grep -v sample_name); do
	echo ${sample}
	./extract_kraken2_taxids_edit.py -r ../../../output/cfspi/samples/results/kraken2_report/after_host_mapping/${sample}_EPRSc2_conf0.4.report -t ${taxid}  \
		--include-children --output_file ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_taxIDs.txt

	awk 'NR == FNR { keywords[$1]=1; next; } { if ($3 in keywords) print $2; }' ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_taxIDs.txt ../../../output/cfspi/samples/results/kraken2_output/after_host_mapping/${sample}_EPRSc2_conf0.4.output > ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt

	nr_lines=$( wc -l  ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt | cut -d" " -f 1 )
	echo $nr_lines

	##mapping for reads classified as taxID
	if [ $nr_lines != "0" ]; then
		echo "mapping reads classified to this taxa"
		filterbyname.sh in1=../../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r1.fq in2=../../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r2.fq out1=../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq out2=../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq names=../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt include=TRUE overwrite=TRUE

		for AspGen in $(ls ../resources/EPRSc2_Aspergillus/*fna ); do
			echo $AspGen
			sh_AspGen=$(echo $AspGen | sed 's|.*FungiDB-46_||g' | sed 's|_Genome_cleaned_final.fna||g')
			mp_AspGen=$(echo $AspGen | sed 's|.fna$||g')


			output_stem="../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}"

			bowtie2 -p 8 -x ${mp_AspGen} -1 ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq -2 ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq > ${output_stem}.sam;
			samtools view -b ${output_stem}.sam -o ${output_stem}.bam;
			samtools sort ${output_stem}.bam -o ${output_stem}_srt.bam
			samtools index ${output_stem}_srt.bam

			rm ${output_stem}.sam
			rm ${output_stem}.bam

			##get TLEN en motifs
			samtools view -q 30 -F 256 -F 2048 -o ${output_stem}_srt_tmp.bam ${output_stem}_srt.bam;
			samtools view -c ${output_stem}_srt_tmp.bam -L ../../01_PreProcessing/resources/nonMT.txt > ${output_stem}_readCount.txt;

			samtools view -f 64 -o ${output_stem}_srt_tmp_R1.bam ${output_stem}_srt_tmp.bam;

			samtools view ${output_stem}_srt_tmp_R1.bam | \
				awk '{{print $9, substr($10, 1, 3)}}' | sort |  uniq -c | \
				awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R1_TLEN_EndMotif.txt;

			samtools view -f 128 -o ${output_stem}_srt_tmp_R2.bam ${output_stem}_srt_tmp.bam;

			samtools view ${output_stem}_srt_tmp_R2.bam | \
				awk '{{print $9, substr($10, 1, 3)}}' | sort |  uniq -c | \
				awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R2_TLEN_EndMotif.txt;

		done

		rm ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq
		rm ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq
		#rm ${output_stem}_srt_tmp.bam
		rm ${output_stem}_srt_tmp_R1.bam
		rm ${output_stem}_srt_tmp_R2.bam
	fi

	rm ../../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt
done
done
