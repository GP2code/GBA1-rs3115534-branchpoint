#!/bin/sh

MAIN="./10X/HBCC_82040_10X_GBA/"
mkdir $MAIN/cells_fastq/
IN=$1
grep -A 2 -B 1 -f $MAIN/cell_types/$1.txt $MAIN/gba_hg38.fastq > $MAIN/cells_fastq/grep/$1.fastq
sed -i '/^--$/d' $MAIN/cells_fastq/grep/$1.fastq
wc -l $MAIN/cells_fastq/grep/$1.fastq >> $MAIN/cells_fastq/grep/cell_numbers.txt


