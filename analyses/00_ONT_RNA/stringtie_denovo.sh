#!/bin/sh

# Set directory variables 
MAIN="./GBA1_CRISPR/RNA/stringtie/"


OUT=$1
prefix=$2
input=$3

module load stringtie
stringtie -L -R -m 50 -o ${OUT}/${prefix}_denovo.gtf ${input}