Structural Variants
=====================
.. highlight:: bash

- SDU
The following script developed at SDU might be used to infer breakpoints from DNA-seq data using Breakdancer::

	#!/bin/sh
	### Specify paths to requiered software, input and output folders
	BWA="/path/to/bwa"
	picard="/path/to/picard-tools"
	samtools="/path/to/samtools"
	DAT="/path/to/fastqfiles/"
	OUTPUT="/path/to/outputfolder"
	REF="/path/to/reference"
	BBMAP="/path/to/bbmap"
	BREAKDANCER="/path/to/breakdancer"

	###### Unpacking compressed fastq files
	cd $DAT
	for i in *_R1.fastq.gz;
	do
	newfile=$(basename $i _R1.fastq.gz)
	gunzip $DAT/${newfile}_R1.fastq.gz
	gunzip $DAT/${newfile}_R2.fastq.gz
	done

	###### Quality filtering
	for i in *_R1.fastq;
	do
	newfile=$(basename $i _R1.fastq)
	$BBMAP/bbduk.sh -Xmx20g in1=/$DAT/${newfile}_R1.fastq in2=/$DAT/${newfile}_R2.fastq out1=/$DAT/${newfile}_clean_R1.fastq out2=/$DAT/${newfile}_clean_R2.fastq ref=$BBMAP/resources/adapters.fa ktrim=r ktrim=l k=23 mink=11 hdist=1 tpe tbo qtrim="rl" trimq=10 maq=10 minlen=25
	done

	###### BWA alignment
	for i in *_clean_R1.fastq;
	do
	newfile=$(basename $i _clean_R1.fastq)
	$BWA mem -t 4 -M -R '@RG\tID:${newfile}.lib.run\tLB:${newfile}.lib\tPL:ILLUMINA\tSM:${newfile}' $REF/hg19.fa $DAT/${newfile}_clean_R1.fastq $DAT/${newfile}_clean_R2.fastq > $OUTPUT/${newfile}.sam

	### Samtools processing of aligned reads
	$samtools view -bS -@ 4 $OUTPUT/${newfile}.sam > $OUTPUT/${newfile}.bam 

	rm $OUTPUT/${newfile}.sam

	$samtools sort $OUTPUT/${newfile}.bam -o $OUTPUT/${newfile}_sorted.bam

	$samtools index $OUTPUT/${newfile}_sorted.bam $OUTPUT/${newfile}_sorted.bai

	### PCR duplicate removal
	java -jar $picard/MarkDuplicates.jar INPUT=$OUTPUT/${newfile}_sorted.bam OUTPUT=$OUTOUT/${newfile}_sorted_nodup.bam METRICS_FILE=${newfile}_dup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT
	$samtools index $OUTPUT/${newfile}_sorted_nodup.bam $OUTPUT/${newfile}_sorted_nodup.bai
	rm $OUTPUT/${newfile}_sorted.bam
	rm $OUTPUT/${newfile}_sorted.bai

	done

	### Run breakdancer
	cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
	cd $OUTPUT
	for i in *_sorted_nodup.bam;
	do
	newfile=$(basename $i _sorted_nodup.bam)
	perl $BREAKDANCER/bam2cfg.pl $OUTPUT/${newfile}_sorted_nodup.bam  > $OUTPUT/${newfile}.cfg
	perl $BREAKDANCER/BreakDancerMax.pl $OUTPUT/${newfile}.cfg -d -g -h -y > $OUTPUT/${newfile}.ctx

	done
	
- AAUH
Fusion genes may be inferred from RNA-seq data using the program `FusionMAP <http://www.arrayserver.com/wiki/index.php?title=FusionMap>`_. 
Running Fusionmap on a Linux machine requires installing Mono; instructions can be found on the FusionMap website.
Fusionmap Requires a control file, with the following specifications::

	MonoPath="path to mono"
	PairedEnd=True  //Automatically pair two fastq files as one sample to run fusion analysis
	RnaMode=True    //Detect fusion results
	ThreadNumber=8  //Possible values: 1-100. Default value=1
	FileFormat=FASTQ //Possible values: FASTQ, QSEQ, FASTA, BAM. Default value=FASTQ
	CompressionMethod=Gzip   //Gzip formatted input files
	Gzip=True        //Gzip
	OutputFusionReads=True //Out put Fusion reads as BAM files for genome browser. Default value=False

	<Output>
	TempPath=temp
	OutputPath=output
	OutputName=OUT.file
	
FusionMAP may be run using the script below::

	## Set paths
	mkdir /scratch/$PBS_JOBID
	TMPDIR=/scratch/$PBS_JOBID
	cd $TMPDIR
	rawDataDir="Path to QC fastQ files"
	outDir="Path for output"
	controlDir="path to control file"
	mono="path to Mono"
	fusionmap="path to fusionMap"

	## Paired fastq files
	R1=$(ls $rawDataDir | grep fastq | grep R1)
	R2=$(ls $rawDataDir | grep fastq | grep R2)

	## Move fastq files to scratch on node
	cp $rawDataDir/$R1 $rawDataDir/$R2 $TMPDIR

	cd $TMPDIR

	## Prepare control file
	sed s/R1.file/$R1/g $controlDir/fusionmap.control.txt |\
		sed s/R2.file/$R2/g | sed s/OUT.file/$id/g > fusion.control

	## Run FusionMap
	$mono $fusionmap \
		--semap /data/apps/software/FusionMap_2015-03-31 \
		Human.B37.3 RefGene \
		fusion.control

	mv -fr $TMPDIR/output $outDir
	rm -fr $TMPDIR

