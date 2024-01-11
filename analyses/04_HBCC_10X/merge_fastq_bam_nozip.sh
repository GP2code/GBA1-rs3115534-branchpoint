#!/bin/bash
FASTQ_PATH=$1
SAMPLE_NAME=$2
module load samtools
module load seqkit
module load pycoQC
cd ${FASTQ_PATH}
echo "starting with fastq"
cat *.fastq > ${SAMPLE_NAME}.fastq; \
echo 'number of pass reads in sequencing_summary.txt:'; grep TRUE ../sequencing_summary.txt \
| wc -l; echo 'number of reads in merged fastq:'; awk 'END{print NR/4}' ${SAMPLE_NAME}.fastq; \
echo "done with fastq"
echo "starting with seqkit"
seqkit stats -j 10 -a ${SAMPLE_NAME}.fastq -T > stats.pass.tsv
echo "done with seqkit"
echo "starting with pycoQC"
pycoQC -f ../sequencing_summary.txt -o pycoQC_${SAMPLE_NAME}.html
echo "done with pycoQC"

#how to run:
# sbatch --cpus-per-task=10 --mem=200g --mail-type=END --time=12:00:00 merge_fastq_bam.sh PATH_TO_DATA SAMPLE_NAME
## pass
# sbatch --cpus-per-task=10 --mem=200g --mail-type=END --time=12:00:00 merge_fastq_bam.sh /data/CARD/projects/KOLF_genome/RNA/CELL_IPSC_RNA_TEST_PROM003_PROCESSING/pass CELL_IPSC_RNA_TEST_PROM003

## fail
# sbatch --cpus-per-task=10 --mem=200g --mail-type=END --time=12:00:00 merge_fastq_bam.sh /data/CARD/projects/KOLF_genome/RNA/CELL_IPSC_RNA_TEST_PROM003_PROCESSING/fail CELL_IPSC_RNA_TEST_PROM003
