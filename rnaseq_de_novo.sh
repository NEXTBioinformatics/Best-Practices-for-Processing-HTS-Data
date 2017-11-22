#!/bin/sh

DAT="/path/to/fastq" ### PATH TO FASTQ FILES
OUT="/path/to/outputfolder" PATH TO OUTPUT FOLDER
TRI="/path/to/trinity" PATH TO TRINITY
GMAP="/path/to/gmap" PATH TO GMAP
REFG="/path/to/reference_genome" ### PATH TO REFERENCE GENOME
SAM="/path/to/samtools" ## PATH TO SAMTOOLS

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

cd $OUT

$TRI/bin/Trinity --max_memory 510G --seqType fq --SS_lib_type FR --left $DAT/${newfile}_R1.fastq --right $DAT/${newfile}_R2.fastq --CPU 48 --trimmomatic

cd $OUT/trinity_out_dir

$TRIN/bin/TrinityStats.pl Trinity.fasta > ${newfile}_stats.txt

$GMAP/bin/gmap -n 0 -t 48 -D $REFG -d hg19 Trinity.fasta -f samse > ${newfile}_trinity_gmap.sam

$SAM/samtools view -Sb ${newfile}_trinity_gmap.sam > ${newfile}_trinity_gmap.bam

$SAM/samtools sort ${newfile}_trinity_gmap.bam ${newfile}_trinity_gmap

$SAM/samtools index ${newfile}_trinity_gmap.bam

done
 

