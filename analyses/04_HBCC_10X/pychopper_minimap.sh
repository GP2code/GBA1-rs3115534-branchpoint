#!/bin/sh



ml pychopper
ml minimap2
ml samtools


MAIN="./HBCC_82040_10X_GBA/"
mkdir $MAIN/pychopper/
mkdir $MAIN/stats/
mkdir $MAIN/mapped/

IN=$1


pychopper -k PCS111 -r $MAIN/pychopper/stats/${IN}.merged.pdf -S $MAIN/pychopper/stats/${IN}.merged.tsv $MAIN/cells_fastq/merged/${IN}_unique_reads.fastq $MAIN/pychopper/${IN}.merged.fastq


minimap2 -ax splice -k 14 -uf ./resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa $MAIN/pychopper/${IN}.merged.fastq - | samtools view -q 40 -F 2304 -Sb - | samtools sort -@ 10 - > $MAIN/mapped/${IN}.pychopper.sorted.bam

samtools index $MAIN/mapped/${IN}.pychopper.sorted.bam

