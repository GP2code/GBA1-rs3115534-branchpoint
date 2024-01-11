#!/bin/sh

MAIN="./HBCC_82040_10X_GBA/"
PREFIX=$1


# Read and trim sequences from file A into an array
mapfile -t sequences_A < "$MAIN/cell_types/${PREFIX}.txt"
sequences_A=("${sequences_A[@]// /}")

# Read sequences from file B into a variable
sequences_B=$(< "$MAIN/cells_fastq/grep/${PREFIX}.fastq")

# Initialize an associative array to store counts
declare -A sequence_counts

# Loop through each query sequence from file A
for query_seq in "${sequences_A[@]}"; do
  # Use grep to count occurrences of the trimmed query sequence in file B
  count=$(grep -o -F -c "$query_seq" <<< "$sequences_B")
  
  # Store the count in the associative array
  sequence_counts["$query_seq"]=$count
done

# Print the counts
for query_seq in "${sequences_A[@]}"; do
  echo "$query_seq: ${sequence_counts[$query_seq]}" >> $MAIN/sequence_diversity/${PREFIX}_counts.txt
done

count=$(grep -c ': [^0]' $MAIN/9CYCLES/sequence_diversity/nopychopper/agrep/${PREFIX}_counts.txt)
total=$(wc -l < $MAIN/9CYCLES/sequence_diversity/nopychopper/agrep/${PREFIX}_counts.txt)
echo "${PREFIX}" >> $MAIN/sequence_diversity/total_sequence_diversity.txt
echo "Count of entries with nonzero values: $count" >> $MAIN/sequence_diversity/total_sequence_diversity.txt
echo "Total number of lines in the file: $total" >> $MAIN/sequence_diversity/total_sequence_diversity.txt
result=$(echo "scale=2; $count / $total" | bc) >> $MAIN/sequence_diversity/total_sequence_diversity.txt
echo "Result: $result" >> $MAIN/sequence_diversity/total_sequence_diversity.txt

