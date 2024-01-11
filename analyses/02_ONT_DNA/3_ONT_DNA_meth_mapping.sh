#!/bin/bash

#### For DNA with methylation tags !!!!


FASTQ_PATH=$1
BAM_FILE=$2
OUT_PATH=$3
SAMPLE_NAME=$4




echo "name ${SAMPLE_NAME}"
echo "Processing ${SAMPLE_NAME}"
echo "Input from ${FASTQ_PATH}"
echo "Output going to ${OUT_PATH}"

echo "map with minimap2"
module load minimap2
module load samtools

#R9
samtools fastq -TMm,Ml ${FASTQ_PATH}/${BAM_FILE} |  minimap2 -y -x map-ont -t 20 -a --eqx -k 17 -K 10g /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa - | samtools view -@ 10 -bh - | samtools sort -@ 10 - > ${OUT_PATH}/${SAMPLE_NAME}.sorted.bam

#R10
#samtools fastq -TMM,ML ${FASTQ_PATH}/${BAM_FILE} |  minimap2 -y -x map-ont -t 20 -a --eqx -k 17 -K 10g #/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa - | samtools view -@ 10 -bh - | samtools sort -@ 10 - > ${OUT_PATH}/${SAMPLE_NAME}.sorted.bam

samtools index ${OUT_PATH}/${SAMPLE_NAME}.sorted.bam

#how to run:
#### sbatch --cpus-per-task=24 --mem=300g --mail-type=BEGIN,END --time=4-0 3_mapping_sort.sh \
#### /data/Neuro_Longread/KOLF_meth_Pilar/MG_IPS_NEURON/unmapped/ \
#### merged_all.bam \
#### /data/Neuro_Longread/KOLF_meth_Pilar/MG_IPS_NEURON/mapped/ \
#### mg_ips_neuron_mapped 
