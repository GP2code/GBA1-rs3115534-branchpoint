#!/bin/sh
ml samtools
ml minimap2

MAIN="./data/HBCC/10X/HBCC_82040_10X_GBA"

IN=$1 
minimap2 -ax splice -k 14 -uf ./resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${IN}.fastq.gz - | samtools view -@ 10 -bh - | samtools sort -@ 10 - > ${IN}.hg38.sorted.bam

samtools index ${IN}.hg38.sorted.bam


