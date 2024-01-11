#!/bin/sh

# Set directory variables 
MAIN="./GBA1_CRISPR/RNA/stringtie/"

ml stringtie

OUT=$1
prefix=$2
input=$3

stringtie -v -L -G GBA1_GBAP1_ours_new.gtf ${input} > ${OUT}/${prefix}_ref.gtf



