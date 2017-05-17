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

 