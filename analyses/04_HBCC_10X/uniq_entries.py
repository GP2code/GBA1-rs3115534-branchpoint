from Bio import SeqIO
import sys
import os

if len(sys.argv) != 2:
    print("Usage: python uniq_entries.py <prefix>")
    sys.exit(1)

prefix = sys.argv[1]
fastq_file = f"{prefix}.gba.fastq"
output_file = f"{prefix}_unique_reads.fastq"

# Set to store unique read IDs
unique_read_ids = set()

with open(fastq_file, "r") as f_in, open(output_file, "w") as f_out:
    for record in SeqIO.parse(f_in, "fastq"):
        read_id = record.id
        if read_id not in unique_read_ids:
            # If read ID is unique, write the entry to the output file
            SeqIO.write(record, f_out, "fastq")
            unique_read_ids.add(read_id)

print(f"Unique reads written to: {output_file}")