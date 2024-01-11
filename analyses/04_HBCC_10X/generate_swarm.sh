#!/bin/bash

# Input and output file paths
MAIN="./HBCC_82040_10X_GBA/cell_types/"
input_file="$1.txt"
output_swarm="$1.swarm"
output_fastq_dir="$MAIN/agrep_"$1"/"


# Loop through each line in the input file
while IFS= read -r line; do
    # Generate and print the grep command
    grep_command="agrep -1 "$line" ./HBCC_82040_10X_GBA/gba_hg38.fastq >> $output_fastq_dir/$line.txt"
    echo "$grep_command" >> "$output_swarm"
done < "$input_file"

# Output a message
echo "Swarm script generated: $output_swarm"
