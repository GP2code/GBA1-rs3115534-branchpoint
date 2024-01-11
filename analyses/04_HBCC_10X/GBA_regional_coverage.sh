#!/bin/sh

# Load required modules
module load samtools

# Set up required variables
MAIN="./10X/HBCC_82040_10X_GBA/"
SAMPLE_GENO_NAME=$1
BAM_PATH=$2
SAMPLE_NAME=$3

# Get coverage statistics
cat $MAIN/depth/GBA1.bed | while read -r line
do
    chr=$(echo $line | cut -d" " -f1)
    start=$(echo $line | cut -d" " -f2)
    end=$(echo $line | cut -d" " -f3)
    region=$(echo $line | cut -d" " -f4)
    cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line
    do
        id=$(echo $line | cut -d" " -f1)
        geno=$(echo $line | cut -d" " -f2)
            samtools coverage -q5 -Q20 --ff UNMAP,SECONDARY,QCFAIL,DUP -r $chr:$start-$end ${BAM_PATH}/$id.pychopper.sorted.bam | grep -v "#" >> $MAIN/depth/$id.cov
    done 
done

# Add in ID and genotype
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
        awk -F '\t' -v OFS='\t' -v var="$id" '{ $(NF+1) = var; print }' $MAIN/depth/$id.cov > $MAIN/depth/$id.id.cov
done


# Add regions
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    geno=$(echo $line | cut -d" " -f2)
    paste -d '\t' $MAIN/depth/$id.id.cov $MAIN/depth/GBA1_regions.txt > $MAIN/depth/$id.id.regions.cov
done


# Get unique read counts
cat $MAIN/cell_types/${SAMPLE_GENO_NAME} | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    wc -l $MAIN/cell_types/$id.txt >> $MAIN/depth/uniq_cells.txt
done



# Match them up with id
paste -d '\t' $MAIN/cell_types/${SAMPLE_NAME} $MAIN/depth/uniq_cells.txt > $MAIN/depth/uniq_cells_ids.txt


# Add reads column to table
cat $MAIN/depth/uniq_cells_ids.txt | while read -r line 
do
    id=$(echo $line | cut -d" " -f1)
    reads=$(echo $line | cut -d" " -f2)
    awk -F '\t' -v OFS='\t' -v var="$reads" '{ $(NF+1) = var; print }' $MAIN/depth/$id.id.regions.cov > $MAIN/depth/$id.id.regions.cells.cov
done


# Cat files

cat $MAIN/depth/*id.regions.cells.cov* > $MAIN/depth/all_table_nh.txt

# Create normalized column
awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $7/$12}' $MAIN/depth/all_table_nh.txt > $MAIN/depth/all_table_nh_norm.txt

# Add in header
echo -e "CHR\tSTART\tEND\tNUMREADS\tCOVBASES\tCOVERAGE\tMEANDEPTH\tMEANBASEQ\tMEANMAPQ\tSAMPLE\tREGION\tUNIQCELLS\tNORMALIZEDDEPTH" | cat - $MAIN/depth/all_table_nh_norm.txt > $MAIN/depth/regional_cov_all.txt

mkdir $MAIN/depth/tmp/
mv *.cov $MAIN/depth/tmp/
mv *nh* $MAIN/depth/tmp/
mv *uniq* $MAIN/depth/tmp/

echo "Done"