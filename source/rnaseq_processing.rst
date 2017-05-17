Processing of RNA-seq
=====================
- SDU
.. highlight:: bash

The SDU workflow for processing RNA-seq data is given as::

	#!/bin/sh

	#### PATHS TO REQUIRED SOFTWARE, INPUT & OUTPUT FOLDERS
	DAT="/path/to/data" ### PATH TO FOLDER CONTAINING GZ-COMPRESSED FASTQ FILES
	OUT="/path/to/output" ### PATH TO FOLDER DEPOSITING THE RESULTS
	REFT="/path/to/reference" ### PATH TO FOLDER TRANSCRIPTOME REFERENCE (gtf/gff file)
	REFG="/path/to/reference" ### PATH TO FOLDER GENOME REFERENCE (fa file and bowtie index files)  
	SAMTOOLS="/path/to/samtools" ### PATH TO SAMTOOLS 
	TOPHAT="/path/to/tophat2" ### PATH TO TOPHAT2
	BOWTIE="/path/to/bowtie" ### PATH TO BOWTIE 
	BBMAP="/path/to/BBMAP" ### PATH TO BBMAP
	PYTHON="/path/to/python"
	HTSEQ="/path/to/HTSeq"

	export PATH=$SAMTOOLS:$BOWTIE:$TOPHAT:$HTSEQ:$PATH
	echo $PATH

	### UNPACKING gz-compressed fastq files
	cd $DAT
	for i in *_R1.fastq.gz;
	do
	newfile=$(basename $i _R1.fastq.gz)
	gunzip /$DAT/${newfile}_R1.fastq.gz
	gunzip /$DAT/${newfile}_R2.fastq.gz
	done

	### Quality cleaning (adaptor removal, trimming of low quality bases and reads)
	for i in *_R1.fastq;
	do
	newfile=$(basename $i _R1.fastq)
	$BBMAP/bbduk.sh -Xmx20g \
			in1=/$DAT/${newfile}_R1.fastq \
			in2=/$DAT/${newfile}_R2.fastq \
			out1=/$DAT/${newfile}_clean_R1.fastq \
			out2=/$DAT/${newfile}_clean_R2.fastq \
			ref=$BBMAP/resources/adapters.fa \
			ktrim=r ktrim=l k=23 mink=11 hdist=1 tpe tbo \
			qtrim="rl" trimq=10 maq=10 minlen=25
	done

	#### Bowtie2-Tophat2 alignment
	for i in *_clean_R1.fastq;
	do
	newfile=$(basename $i _clean_R1.fastq)
	$TOPHAT -p 1 -G $REFT \
			--output-dir $OUT/${newfile} \
			$REFG $DAT/${newfile}_clean_R1.fastq \
			$DAT/${newfile}_clean_R2.fastq 
	$SAMTOOLS/samtools index $OUT/${newfile}/accepted_hits.bam
	mv $OUT/${newfile}/accepted_hits.bam $OUT/${newfile}/${newfile}_accepted_hits.bam
	mv $OUT/${newfile}/accepted_hits.bam.bai $OUT/${newfile}/${newfile}_accepted_hits.bam.bai

	### Count Matrix construction by HTSeq
	$python -m HTSeq.scripts.count \
			--format bam \
			--mode union \
			--stranded no \
			--minaqual 1 \
			--type gene \
			--idattr gene_id \
			$OUT/${newfile}/${newfile}_accepted_hits.bam $REFT \
			> $OUT/${newfile}_gene_read_counts_table.tsv

	done
	
- AAUH

 The spliced alignment of RNA-seq performed with tophat in the above script can altertively be done using a 2-pass alignment with `STAR <https://github.com/alexdobin/STAR>`_. 
 This is the suggested method in the GATK best practices
 
A script for doing this is shown below::

	## Create temp dir
	mkdir /scratch/$PBS_JOBID
	TMPDIR=/scratch/$PBS_JOBID
	cd $TMPDIR

	## Path to paired fastq files
	rawDataDir="Path to directory with QC fastq files"
	R1=$(ls $rawDataDir | grep $id | grep R1)
	R2=$(ls $rawDataDir | grep $id | grep R2)

	## Move fastq files to scratch
	cp $rawDataDir/$R1 $rawDataDir/$R2 $TMPDIR

	## Programs
	STAR="path to star"

	## Reference data
	assembly="Path to reference genome fasta"
	genomeDir="Path to STAR indexed reference genome"

	##########################
	#### Align using STAR ####
	##########################
	
	### Do 1st pass
	mkdir $TMPDIR/1pass
	cd $TMPDIR/1pass

	$STAR \
	--genomeDir $genomeDir \
	--readFilesIn ../$R1 ../$R2 \
	--readFilesCommand zcat \
	--runThreadN 8

	### Create new index using splice junction information from 1st pass
	mkdir $TMPDIR/b37_2pass

	$STAR \
	--runMode genomeGenerate \
	--genomeDir $TMPDIR/b37_2pass \
	--genomeFastaFiles $assembly \
	--sjdbFileChrStartEnd $TMPDIR/1pass/SJ.out.tab \
	--sjdbOverhang 75 \
	--genomeSAsparseD 2 \
	--runThreadN 8 \
	--limitGenomeGenerateRAM 20000000000
	
	### Do 2nd pass
	mkdir $TMPDIR/2pass
	cd $TMPDIR/2pass
	$STAR \
	--genomeDir $TMPDIR/b37_2pass \
	--readFilesIn ../$R1 ../$R2 \
	--readFilesCommand zcat \
	--runThreadN 8 \
	--outSAMstrandField intronMotif \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${id}_STAR_
	#--outSaMmapqUnique 60 \

	## Return output
	cp * $outDir

	## Clean up scratch
	cd /scratch
	rm -fr $PBS_JOBID
	
Transcipt level expression can then be inferred using `Cufflinks <https://github.com/cole-trapnell-lab/cufflinks>`_.
This is done using the script below::

	## Paths
	mkdir /scratch/$PBS_JOBID
	TMPDIR=/scratch/$PBS_JOBID
	bamFile="path to STAR aligned BAM file"
	cufflinks="path to cufflinks"
	gff="Path to gff file"
	outDir="Path for output files"
	
	cd $TMPDIR
	cp $bamFile $TMPDIR

	# Construct the mask file
	grep rRNA $gff > mask.gff3
	grep tRNA $gff >> mask.gff3

	# Run cufflinks
	echo Running cufflinks ...
	$cufflinks \
	--GTF-guide $gff \
	--mask-file mask.gff3 \
	--library-type fr-unstranded \
	--num-threads 8 \
	--output-dir cufflinks \
	--quiet \
	$bamFile

	echo moving files to home dir
	cp -fr $TMPDIR/cufflinks/* $outDir


