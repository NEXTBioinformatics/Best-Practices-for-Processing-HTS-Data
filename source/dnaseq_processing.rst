Processing of DNA-seq
=====================
.. highlight:: bash

The workflow below requires specification of paths to a number of datasets and programs::
	
	## Datasets
	assembly="path to reference genome in .fasta format"
	known_snp="path to dbSNP variants in .vcf format"
	known_snp_hc="path to high confidence 1000G SNPs in .vcf format"
	known_indels="path to Mills and 1000g INDELs in .vcf format" 
	cosmic="path to .vcf file with known coding and non-coding somatic mutations in the COSMIC database"
	target_regions="path to bed file with target regions for capture kit"
	
	## Programs
	java="path to java 1.8 or newer"
	gatk="path to GATK .jar file"
	picard="path to picard .jar file"
	varscan2="path to varscan2 .jar file"
	samtools="path to samtools"
	
- Alignment

Paired end reads are mapped to a reference genome using BWA MEM::

	/data/apps/bwa-0.7.12/bwa mem \
	-M \
	-t28 \
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
	METRICS_FILE=dup_metrics.txt

	## Index bam file
	$java -jar $picard BuildBamIndex \
	INPUT=dedup_reads.bam
	\end{lstlisting}

- INDEL realignment

Although INDEL realignment is included in the haplotypecaller used in MUTECT2 the reads are realigned in the .bams, so that they look more reasonable when manually validating variants with IGV::

	## Create a target list of intervals to be realigned
	java -jar $gatk \
	-T RealignerTargetCreator \
	-R $assembly \
	-I dedup_reads.bam \
	-known $known_indels \
	-o realignment_targets.list

	## Perform realignment of the target intervals
	java -jar $gatk \
	-T IndelRealigner \
	-R $assembly \
	-I dedup_reads.bam \
	-targetIntervals realignment_targets.list \
	-known $known_indels \
	-o realigned_dedup_reads.bam

	
- Base Quality Recalibration

Quality scores in BAM files are recalibrated to adjust for bias in the quality scores which might affect the final variant calls::

	## Analyze patterns of covariation in the sequence dataset
	$java -jar $gatk \
	-T BaseRecalibrator \
	-R $assembly \
	-nct 28 \
	-I realigned_dedup_reads.bam \
	-knownSites $known_snp \
	-knownSites $known_snp_hc \
	-knownSites $known_indels \
	-o recal_data.table

	## Do a second pass to analyze covariation remaining after recalibration
	$java -jar $gatk \
	-T BaseRecalibrator \
	-R $assembly \
	-nct 8 \
	-I realigned_dedup_reads.bam \
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
	-I realigned_dedup_reads.bam \
	-BQSR recal_data.table \
	-o recal_realigned_dedup_reads.bam

- Alignment Metrics

quality metrics for alignment and duplication are calculated using picard tools::

	$java -jar $picard BedToIntervalList \
	I=$target_regions \
	O=target_Picard \
	SD=$dictionary

	$java -jar $picard CollectHsMetrics \
	I=recal_realigned_dedup_reads.bam \
	O=HSmetrics.txt \
	R=$assembly \
	TARGET_INTERVALS=target_Picard \
	BAIT_INTERVALS=target_Picard
	
		
- Variant calling for somatic mutations

Somatic variants are called using both Mutect2 and Varscan2, and variants are subsequently merged and filtered. A more detailed description is found in the Mutect2 pitfalls section.
Variant calling with Mutect2 can optionally be parallelized by chromosome using the -L parameter for faster runtimes::

	## Run Mutect2
	$java -jar $gatk \
	--analysis_type MuTect2 \
	--reference_sequence $assembly \
	--input_file:normal normal.bam \
	--input_file:tumor tumor.bam \
	--out $inTumor/somatic_variants.vcf \
	--max_alt_alleles_in_normal_count  1000000 \
	--max_alt_allele_in_normal_fraction 0.1 \
	--cosmic $cosmic \
	--dbsnp $known_snp \
	-nct 28
	
Variant calling with varscan2 requires an mpileup file which can be built with samtools using the aligned BAM files for tumor and normal samples::

	## Build mpileup with samtools
	$samtools mpileup \
	-f $assembly \
	-q 1 \
	-B normal.bam \
	tumor.bam > normal-tumor.mpileup

Variants may then be called with varscan2 and high confidence SNPs/INDELs can be extracted using the processSomatic command::

	## Run varscan2 somatic
	$java -jar $varscan2 \
	somatic \
	normal-tumor.mpileup \
	tumor_variants.varscan2 \
	--mpileup 1 \
	--min-var-freq 0.02 \
	--output-vcf

	## Process SNPs
	$java -jar $varscan2 \
	processSomatic \
	tumor_variants.varscan2.snp.vcf
	
	## Process INDELs
	$java -jar $varscan2 \
	processSomatic \
	tumor_variants.varscan2.indel.vcf

- Variant filtration

The final set of somatic SNPs / INDELS are found by combining and filtering outputs from Mutect2 and varscan2 as described in the Mutect2 pitfalls section.
Briefly, for a variant to pass filtering the following must be fulfilled::

	1) PASS in Mutect2 or called by MuTect2 + PASS in varscan2 HC
	2) Tumor AF > 4 * Normal AF
	3) QSS / AD > 25