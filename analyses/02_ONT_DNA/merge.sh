#!/bin/sh
MAIN=./GBA1_CRISPR/DNA/
ml samtools

PREFIX=$1

cd $MAIN/CRISPR_"$1"_DNA/mapped/

# Count the number of .bam files in the current directory
bam_count=$(ls -1 *.bam 2>/dev/null | wc -l)

# Check the value of bam_count
if [ $bam_count -eq 2 ]; then
  # If there are two .bam files, merge them into merged.bam
  samtools merge -o "$1"_merged.bam *.bam
  samtools index "$1"_merged.bam
  echo "Merged two .bam files into merged.bam"
elif [ $bam_count -eq 1 ]; then
  # If there is only one .bam file, rename it to merged.bam
  mv *.bam "$1"_merged.bam
  mv *.bam.bai "$1"_merged.bam.bai
  echo "Renamed the single .bam file to merged.bam"
else
  # If there are zero or more than two .bam files, display an error message
  echo "Error: There are $bam_count .bam files in the directory. Expected 1 or 2."
fi

