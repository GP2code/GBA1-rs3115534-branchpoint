#!/bin/sh



ml pychopper
ml minimap2
ml samtools


MAIN=./GBA1_CRISPR/RNA/
IN=$1
PREFIX=$2


pychopper -k PCS111 -r $MAIN/pychopper/stats/${IN}.merged.pdf -S $MAIN/pychopper/stats/${IN}.merged.tsv $MAIN/CRISPR_${PREFIX}_RNA/${IN}.fastq.gz $MAIN/pychopper/pychopper_fastqs/${PREFIX}.fastq.gz


minimap2 -ax splice -k 14 -uf ./resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa $MAIN/pychopper_fastqs/${PREFIX}.fastq.gz - | samtools view -q 40 -F 2304 -Sb - | samtools sort -@ 10 - > $MAIN/minimap2/${PREFIX}.sorted.bam

samtools index $MAIN/minimap2/${PREFIX}.sorted.bam
