#!/bin/sh

# Load required modules
module load samtools

# Set up required variables
MAIN=./GBA1_CRISPR/RNA/depth/
SAMPLE_GENO_NAME=$1
BAM_PATH=$2
SAMPLE_NAME=$3

# Get coverage statistics
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r chr1:155234452-155244627 ${BAM_PATH}/"$id".sorted.bam | grep -v "#" >> $MAIN/"$id"_whole.cov
    done 


# Add in ID and genotype
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
        awk -F '\t' -v OFS='\t' -v var="$id\t$geno" '{ $(NF+1) = var; print }' $MAIN/"$id"_whole.cov > $MAIN/"$id"_whole.id.cov
done



# Get unique reads
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    samtools view -c -F 260 ${BAM_PATH}/"$id"_GBA.bam >> $MAIN/uniq_reads_whole.txt
done

# Transform to per million
while read -r line; do
  result=$(echo "scale=6; $line / 1000000" | bc)
  echo $result >> $MAIN/uniq_reads_million.txt
done < $MAIN/uniq_reads.txt


# Match them up with id
paste -d '\t' $MAIN/${SAMPLE_NAME} $MAIN/uniq_reads_million.txt > $MAIN/uniq_reads_ids_whole.txt


# Add reads column to table
cat $MAIN/uniq_reads_ids_whole.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/"$id"_whole.id.cov > $MAIN/"$id"_whole.id.reads.cov
done


# Cat files

cat $MAIN/*_whole.id.reads.cov* > $MAIN/all_table_nh_whole.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $7/$12}' $MAIN/all_table_nh_whole.txt > $MAIN/all_table_nh_whole_norm.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tGENOTYPE\tUNIQREADS\tNORMALIZEDDEPTH" | cat - $MAIN/all_table_nh_whole_norm.txt > $MAIN/cov_all_whole.txt

mkdir $MAIN/tmp/
mv *.cov $MAIN/tmp/
mv *nh* $MAIN/tmp/

echo "Done"





# Same thing for GBAP1

# Get coverage statistics
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r chr1:155213825-155227534  ${BAM_PATH}/"$id".sorted.bam | grep -v "#" >> $MAIN/"$id"_whole.GBAP1.cov
    done 


# Add in ID and genotype
cat $MAIN/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
        awk -F '\t' -v OFS='\t' -v var="$id\t$geno" '{ $(NF+1) = var; print }' $MAIN/"$id"_whole.GBAP1.cov > $MAIN/"$id"_whole.id.GBAP1.cov
done



# Match them up with id
paste -d '\t' $MAIN/${SAMPLE_NAME} $MAIN/uniq_reads_million.txt > $MAIN/uniq_reads_ids_whole.txt
#head ./STAR/gba1/depth/coriell_ilm_uniq_reads_ids.txt


# Add reads column to table
cat $MAIN/uniq_reads_ids_whole.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/"$id"_whole.id.GBAP1.cov > $MAIN/"$id"_whole.id.reads.GBAP1.cov
done


# Cat files

cat $MAIN/*_whole.id.reads.GBAP1.cov* > $MAIN/all_table_nh_whole.GBAP1.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $7/$12}' $MAIN/all_table_nh_whole.GBAP1.txt > $MAIN/all_table_nh_whole_norm.GBAP1.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tGENOTYPE\tUNIQREADS\tNORMALIZEDDEPTH" | cat - $MAIN/all_table_nh_whole_norm.GBAP1.txt > $MAIN/cov_all_whole.GBAP1.txt

mkdir $MAIN/tmp/
mv *.cov $MAIN/tmp/
mv *nh* $MAIN/tmp/

echo "Done"

