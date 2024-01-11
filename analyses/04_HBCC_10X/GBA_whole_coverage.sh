#!/bin/sh

# Load required modules
module load samtools

# Set up required variables
MAIN="./10X/HBCC_82040_10X_GBA/"
SAMPLE_GENO_NAME=$1
BAM_PATH=$2
SAMPLE_NAME=$3

# Get coverage statistics
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r chr1:155234452-155244627 ${BAM_PATH}/$id.pychopper.sorted.bam | grep -v "#" >> $MAIN/depth/"$id"_whole.cov
    done 


# Add in ID and genotype
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
        awk -F '\t' -v OFS='\t' -v var="$id" '{ $(NF+1) = var; print }' $MAIN/depth/"$id"_whole.cov > $MAIN/depth/"$id"_whole.id.cov
done



# Get unique read counts
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    wc -l $MAIN/cell_types/$id.txt >> $MAIN/depth/uniq_cells_whole.txt
done


# Match them up with id
paste -d '\t' $MAIN/cell_types/${SAMPLE_NAME} $MAIN/depth/uniq_cells_whole.txt > $MAIN/depth/uniq_cells_ids_whole.txt
#head ./STAR/gba1/depth/coriell_ilm_uniq_reads_ids.txt


# Add reads column to table
cat $MAIN/depth/uniq_cells_ids_whole.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/depth/"$id"_whole.id.cov > $MAIN/depth/"$id"_whole.id.cells.cov
done


# Cat files

cat $MAIN/depth/*_whole.id.cells.cov* > $MAIN/depth/all_table_nh_whole.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $7/$11}' $MAIN/depth/all_table_nh_whole.txt > $MAIN/depth/all_table_nh_whole_norm.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tUNIQREADS\tNORMALIZEDDEPTH" | cat - $MAIN/depth/all_table_nh_whole_norm.txt > $MAIN/depth/cov_all_whole.txt

mkdir $MAIN/depth/tmp/
mv *.cov $MAIN/depth/tmp/
mv *nh* $MAIN/depth/tmp/
mv *uniq* $MAIN/depth/tmp/

echo "Done"


# Same thing for GBAP1
# Get coverage statistics
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r chr1:155213825-155227534 ${BAM_PATH}/$id.pychopper.sorted.bam | grep -v "#" >> $MAIN/depth/"$id"_whole_GBAP1.cov
    done 


# Add in ID and genotype
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
        awk -F '\t' -v OFS='\t' -v var="$id" '{ $(NF+1) = var; print }' $MAIN/depth/"$id"_whole_GBAP1.cov > $MAIN/depth/"$id"_whole_GBAP1.id.cov
done



# Get unique read counts
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    wc -l $MAIN/cell_types/$id.txt >> $MAIN/depth/uniq_cells_whole.txt
done


# Match them up with id
paste -d '\t' $MAIN/cell_types/${SAMPLE_NAME} $MAIN/depth/uniq_cells_whole.txt > $MAIN/depth/uniq_cells_ids_whole_GBAP1.txt
#head ./STAR/gba1/depth/coriell_ilm_uniq_reads_ids.txt


# Add reads column to table
cat $MAIN/depth/uniq_cells_ids_whole_GBAP1.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/depth/"$id"_whole_GBAP1.id.cov > $MAIN/depth/"$id"_whole_GBAP1.id.cells.cov
done


# Cat files

cat $MAIN/depth/*_whole_GBAP1.id.cells.cov* > $MAIN/depth/all_table_nh_whole_GBAP1.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $7/$11}' $MAIN/depth/all_table_nh_whole_GBAP1.txt > $MAIN/depth/all_table_nh_whole_GBAP1_norm.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tUNIQCELLS\tNORMALIZEDDEPTH" | cat - $MAIN/depth/all_table_nh_whole_GBAP1_norm.txt > $MAIN/depth/cov_all_whole_GBAP1.txt

mv *.cov $MAIN/depth/tmp/
mv *nh* $MAIN/depth/tmp/
mv *uniq* $MAIN/depth/tmp/

echo "Done"
