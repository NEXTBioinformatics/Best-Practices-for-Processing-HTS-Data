Processing of DNA-seq
========
- Alignment

How to map paired end reads to a reference genome using BWA MEM::

	/data/apps/bwa-0.7.12/bwa mem \
	-M \
	-t8 \
	-R"@RG\tID:exome_1\tPL:ILLUMINA\tPU:exome.11\tLB:$id\tSM:$id\tCN:AAU" \
	$assembly \
	$rawDataDir/$R1 \
	$rawDataDir/$R2 \
	> aligned_reads.sam
	
- Marking Duplicates

Aligned reads are sorted, converted to BAM, and PCR duplicates are marked::

	## Sort and convert to bam
	$java -jar $picard SortSam \
	INPUT=aligned_reads.sam \
	OUTPUT=sorted_reads.bam \
	SORT_ORDER=coordinate

	## Mark PCR duplicates
	$java -jar $picard MarkDuplicates \
	INPUT=sorted_reads.bam \
	OUTPUT=dedup_reads.bam \
	METRICS_FILE=metrics.txt

	## Index bam file
	$java -jar $picard BuildBamIndex \
	INPUT=dedup_reads.bam
	\end{lstlisting}
	
- Base Quality Recalibration

Quality scores in BAM files are recalibrated to adjust for bias in the quality scores which might affect the final variant calls::

	## Analyze patterns of covariation in the sequence dataset
	$java -jar $gatk \
	-T BaseRecalibrator \
	-R $assembly \
	-nct 8 \
	-I dedup_reads.bam \
	-knownSites $known_snp \
	-knownSites $known_indels \
	-L $regions \
	-o recal_data.table

	## Do a second pass to analyze covariation remaining after recalibration
	$java -jar $gatk \
	-T BaseRecalibrator \
	-R $assembly \
	-nct 8 \
	-I dedup_reads.bam \
	-knownSites $known_snp \
	-knownSites $known_indels \
	-BQSR recal_data.table \
	-L $regions \
	-o post_recal_data.table

	## Apply recalibration
	$java -jar $gatk \
	-T PrintReads \
	-R $assembly \
	-nct 8 \
	-I dedup_reads.bam \
	-BQSR recal_data.table \
	-o recal_reads.bam
	
