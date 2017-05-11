#!/bin/sh

DAT="/path/to/data" ### PATH TO FOLDER CONTAINING GZ-COMPRESSED FASTQ FILES

OUT="/path/to/output" ### PATH TO FOLDER DEPOSITING THE RESULTS

REFT="/path/to/reference" ### PATH TO FOLDER TRANSCRIPTOME REFERENCE (fa file and bowtie index files)

REFG="/path/to/reference" ### PATH TO FOLDER TRANSCRIPTOME REFERENCE (GTF file)  

SAMTOOLS="/path/to/samtools" ### PATH TO SAMTOOLS 

TOPHAT="/path/to/tophat2" ### PATH TO SAMTOOLS

BOWTIE="/path/to/bowtie" ### PATH TO BOWTIE 

BBMAP="/path/to/BBMAP" ### PATH TO BBMAP

PYTHON="/path/to/python"

HTSEQ="/path/to/HTSeq"

export PATH=$SAMTOOLS:$BOWTIE:$TOPHAT:$HTSEQ:$PATH
echo $PATH

cd $DAT

for i in *_R1.fastq.gz;
do
newfile=$(basename $i _R1.fastq.gz)

gunzip /$DAT/${newfile}_R1.fastq.gz
gunzip /$DAT/${newfile}_R2.fastq.gz

done

for i in *_R1.fastq;
do
newfile=$(basename $i _R1.fastq)

/$BBMAP/bbduk.sh -Xmx20g in1=/$DAT/${newfile}_R1.fastq in2=/$DAT/${newfile}_R2.fastq out1=/$DAT/${newfile}_clean_R1.fastq out2=/$DAT/${newfile}_clean_R2.fastq ref=$BBMAP/resources/adapters.fa ktrim=r ktrim=l k=23 mink=11 hdist=1 tpe tbo qtrim="rl" trimq=10 maq=10 minlen=25

done

for i in *_clean_R1.fastq;
do
newfile=$(basename $i _clean_R1.fastq)

$TOPHAT -p 1 -G $REFT --output-dir $OUT/${newfile} $REFG $DAT/${newfile}_clean_R1.fastq $DAT/${newfile}_clean_R2.fastq 

$SAMTOOLS/samtools index $OUT/${newfile}/accepted_hits.bam

mv $OUT/${newfile}/accepted_hits.bam $OUT/${newfile}/${newfile}_accepted_hits.bam
mv $OUT/${newfile}/accepted_hits.bam.bai $OUT/${newfile}/${newfile}_accepted_hits.bam.bai


DAT="/data/mark/data/maria/poolM"
REF="/data/mark/tools/ref"
REFBWT="/data/mark/tools/ref/bwt/"
ALGTOPHAT="/data/mark/tools/alignments/tophat"
FASTX="/data/mark/tools/fastx/bin"
CUF="/data/mark/tools/cufflinks-2.2.1.Linux_x86_64"
SAM="/data/mark/tools/samtools-1.3.1/samtools"
BEDTOOLS="/data/mark/tools/bedtools225/bin/bedtools"


newfile="ms4_il_rep2"

$python -m HTSeq.scripts.count --format bam --mode union --stranded no --minaqual 1 --type gene --idattr gene_id $OUT/${newfile}/${newfile}_accepted_hits.bam $REFT > $OUT/${newfile}_gene_read_counts_table.tsv

done 
 





