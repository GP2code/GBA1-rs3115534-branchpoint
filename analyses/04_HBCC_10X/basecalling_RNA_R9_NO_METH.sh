#!/bin/bash
module load guppy/6.1.2
FAST5_PATH=$1
OUT_PATH=$2

guppy_basecaller -i ${FAST5_PATH} -s ${OUT_PATH} -c dna_r9.4.1_450bps_sup_prom.cfg -x cuda:all -r --read_batch_size 50000 -q 50000 --chunks_per_runner 200

#how to run:
#sbatch --partition=gpu --cpus-per-task=14 --mem=200g --gres=gpu:v100x:2,lscratch:200 --time=4-0 --wrap="bash basecalling.sh /fast5_path /out_path"
