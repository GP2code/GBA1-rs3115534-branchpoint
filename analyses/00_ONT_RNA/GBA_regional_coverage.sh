#!/bin/sh

# Load required modules
module load samtools

# Set up required variables
MAIN=./GBA1_CRISPR/RNA/depth/
SAMPLE_GENO_NAME=$1
BAM_PATH=$2
SAMPLE_NAME=$3

# Get coverage statistics
cat $MAIN/GBA1.bed | while read -r line
do
    chr=$(echo $line | cut -d" " -f1)
    start=$(echo $line | cut -d" " -f2)
    end=$(echo $line | cut -d" " -f3)
    region=$(echo $line | cut -d" " -f4)
    cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r $chr:$start-$end ${BAM_PATH}/"$id".sorted.bam | grep -v "#" >> $MAIN/$id.cov
    done 
done

# Add in ID and genotype
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
        awk -F '\t' -v OFS='\t' -v var="$id\t$geno" '{ $(NF+1) = var; print }' $MAIN/$id.cov > $MAIN/$id.id.cov
done



# Add regions
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    paste -d '\t' $MAIN/$id.id.cov $MAIN/GBA1_regions.txt > $MAIN/$id.id.regions.cov
done


# Get unique read counts
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    samtools view -c -F 260 ${BAM_PATH}/"$id"_GBA.bam >> $MAIN/uniq_reads.txt
done

# Transform to per million
while read -r line; do
  result=$(echo "scale=6; $line / 1000000" | bc)
  echo $result >> $MAIN/uniq_reads_million.txt
done < $MAIN/uniq_reads.txt

# Match them up with id
paste -d '\t' $MAIN/${SAMPLE_NAME} $MAIN/uniq_reads_million.txt > $MAIN/uniq_reads_ids.txt


# Add reads column to table
cat $MAIN/uniq_reads_ids.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/$id.id.regions.cov > $MAIN/$id.id.regions.reads.cov
done


# Cat files

cat $MAIN/*id.regions.reads.cov* > $MAIN/all_table_nh.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $7/$13}' $MAIN/all_table_nh.txt > $MAIN/all_table_nh_norm.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tGENOTYPE\tREGION\tUNIQREADS\tNORMALIZEDDEPTH" | cat - $MAIN/all_table_nh_norm.txt > $MAIN/regional_cov_all.txt

mkdir $MAIN/tmp/
mv *.cov $MAIN/tmp/
mv *nh* $MAIN/tmp/

echo "Done"
