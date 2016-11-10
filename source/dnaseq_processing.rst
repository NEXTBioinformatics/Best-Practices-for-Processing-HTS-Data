Processing of DNA-seq
========
- Alignment

Paired end reads are mapped to a reference genome using BWA MEM::

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
		
- Variant calling for germline variants

Germline variants are called using the Haplotypecaller in GATK::

	$java -Xmx10g -jar $gatk \
	--analysis_type HaplotypeCaller \
	-nct 8 \
	--reference_sequence $assembly \
	-I $inNormal/recal_reads.bam \
	-o $inNormal/variants.vcf \
	-L $regions \
	--dbsnp $known_snp
	
- Variant calling for somatic mutations

Somatic variants are called using Mutect2 which calls somatic SNPs and INDELs simultaneously::

	$java -Xmx10g -jar $gatk \
	--analysis_type MuTect2 \
	--reference_sequence $assembly \
	--input_file:normal $inNormal/recal_reads.bam \
	--input_file:tumor $inTumor/recal_reads.bam \
	--out $inTumor/somatic_variants.vcf \
	--cosmic $cosmic \
	--dbsnp $known_snp \
	-L $regions \
	-nct 8
	
Somatic variants may subsequently be annotated with e.g. cancer specific information using Oncotator::

	## Filter out variants with PASS
	/data/apps/vcftools_0.1.13/bin/vcftools \
	--vcf somatic_variants.vcf \
	--remove-filtered-all \
	--out somatic_variants.filtered \
	--recode

	## Start Virtual Machine
	source /data/users/rasmus/software/oncotator_vm_1.9/bin/activate

	## Run oncotator
	/data/users/rasmus/software/oncotator_vm_1.9/bin/oncotator \
	-i VCF \
	-o TCGAMAF \
	--db-dir /data/appdata/oncotator_v1_ds_Jan262014/ \
	somatic_variants.filtered.recode.vcf \
	somatic_variants_filtered_oncotator.maf \
	hg19