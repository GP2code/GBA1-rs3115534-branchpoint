#!/bin/bash
#module load guppy/6.3.8
#module load guppy/6.4.2
module load guppy/6.1.2 #withR9
FAST5_PATH=$1
OUT_PATH=$2

#dna_r10.4.1_e8.2_260bps_hac_prom.cfg
#dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg
guppy_basecaller -i ${FAST5_PATH} -s ${OUT_PATH} -c dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg \
-x cuda:all -r --read_batch_size 50000 -q 50000 --chunks_per_runner 195 --bam_out

# sbatch --partition=gpu --mail-type=END,TIME_LIMIT_80 --cpus-per-task=20 --mem=20g --gres=gpu:v100x:2,lscratch:200 --time=5-0 --wrap="bash $BASE/1_run_guppy_basecaller_less_readsize.sh $FAST5PATH $OUTPUT_DIR"