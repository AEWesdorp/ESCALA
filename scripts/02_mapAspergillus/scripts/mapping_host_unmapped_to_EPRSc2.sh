#mapping of all non-host-mapped reads to Aspergillus genomes 
for sample in $(cut -f 1  ../../01_PreProcessing/cfspi/config/samples.txt | grep -v sample_name); do
	echo ${sample}
	for AspGen in $(ls ../resources/EPRSc2_Aspergillus/*fna ); do
		echo $AspGen

		if [ -e ${output_stem}_srt.bam ]; then
			echo "File ${output_stem}_srt.bam  exists, skipping..."
		else

		sh_AspGen=$(echo $AspGen | sed 's|.*FungiDB-46_||g' | sed 's|_Genome_cleaned_final.fna||g')
		mp_AspGen=$(echo $AspGen | sed 's|.fna$||g')

		output_stem="../../../output/mapAspergillus/${sample}_unmappedHost___${sh_AspGen}"

		bowtie2 -p 8 -x ${mp_AspGen} -1 ../../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r1.fq -2 ../../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r2.fq > ${output_stem}.sam;
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
			awk '$9 > 0 {{print $9, substr($10, 1, 3)}}' | sort |  uniq -c | \
			awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R1pos_TLEN_EndMotif.txt;

		samtools view ${output_stem}_srt_tmp_R1.bam | \
			awk '$9 < 0 {{print $9, substr($10, length($10) - 2 , length($10) )}}' | sort |  uniq -c | \
			awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R1neg_TLEN_EndMotif.txt;

		samtools view -f 128 -o ${output_stem}_srt_tmp_R2.bam ${output_stem}_srt_tmp.bam;

		samtools view ${output_stem}_srt_tmp_R2.bam | \
			awk '$9 > 0 {{print $9, substr($10, 1, 3)}}' | sort |  uniq -c | \
			awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R2pos_TLEN_EndMotif.txt;

		samtools view ${output_stem}_srt_tmp_R2.bam | \
			awk '$9 < 0 {{print $9, substr($10, length($10) - 2 , length($10) )}}' | sort |  uniq -c | \
			awk 'BEGIN {{ OFS=","; print "TLEN,EndMotif,Count" }} {{print $2,$3,$1}}' > ${output_stem}_R2neg_TLEN_EndMotif.txt;

		rm ${output_stem}_srt_tmp_R1.bam
		rm ${output_stem}_srt_tmp_R2.bam
		fi
	done
done

