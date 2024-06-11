# The-genome-architecture-

#all the codes are wriiten in bash scripts.

---------------------------------------------

#to add "speciesname" in front for every sequence accession number for easy identification.

#!/bin/bash

# Directory containing your TXT files
input_dir="/Users/apple/Desktop/dom/subject_db_copy"
output_dir="/Users/apple/Desktop/dom/new"
mkdir -p "$output_dir"  # Ensure the output directory exists

# Process each .txt file in the directory
for file in "$input_dir"/*.txt; do
    # Extract the base filename without the extension
    base_filename=$(basename "$file" .txt)
    
    # Construct the output filename
    output_file="$output_dir/${base_filename}.txt"

    # Process the file and write to the output directory
    sed -e "s/>/>${base_filename}_/" -e 's/|/_/g' "$file" > "$output_file"
    
    echo "Processed and output to $output_file"
done

echo "All files have been processed."

-----------------------------------------------

#to replace "speciesname" to ">speciesname_"in front for every sequence accession number for easy identification.

#!/bin/bash

# Directory containing the files to be processed
input_dir="/Users/apple/Desktop/dom/query_seqs"
output_dir="//Users/apple/Desktop/dom/new"
mkdir -p "$output_dir"  # Create output directory if it doesn't exist

# Loop over all .fa files in the input directory
for file in "$input_dir"/*.fa; do
    # Create an output file path
    output_file="$output_dir/$(basename "$file")"

    # Use sed to replace "human_" with ">human_" and "mouse_" with ">mouse_"
    sed -e 's/human_/>human_/g' -e 's/mouse_/>mouse_/g' "$file" > "$output_file"

    echo "Processed $file and saved modifications to $output_file"
done

echo "All files have been processed."
-----------------------------------------------

#to translate the DNA sequnces to amino acid sequnces

from Bio.Seq import Seq

def translate_dna(input_file, output_file):
    records = SeqIO.parse(input_file, "fasta")
    protein_records = []

    for record in records:
        seq_length = len(record.seq)
        trim_length = seq_length - (seq_length % 3)
        trimmed_seq = record.seq[:trim_length]
        protein_seq = trimmed_seq.translate()
        protein_record = SeqRecord(protein_seq, id=record.id, description="translated sequence")
        protein_records.append(protein_record)

    SeqIO.write(protein_records, output_file, "fasta")


import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_dna(input_file, output_file):
    records = SeqIO.parse(input_file, "fasta")
    protein_records = [SeqRecord(record.seq.translate(), id=record.id, description="translated sequence") for record in records]
    SeqIO.write(protein_records, output_file, "fasta")

def process_files(input_path, output_dir):
    if os.path.isdir(input_path):
        # Process each file in the directory
        for filename in os.listdir(input_path):
            file_path = os.path.join(input_path, filename)
            if filename.endswith(".txt"):  # Ensure processing only text files
                output_file = os.path.join(output_dir, filename.replace('.txt', '_protein.fasta'))
                translate_dna(file_path, output_file)
                print(f"Translated {file_path} to {output_file}")
    elif os.path.isfile(input_path):
        # Process a single file
        filename = os.path.basename(input_path)
        output_file = os.path.join(output_dir, filename.replace('.txt', '_protein.fasta'))
        translate_dna(input_path, output_file)
        print(f"Translated {input_path} to {output_file}")

if __name__ == "__main__":
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    process_files(input_path, output_dir)

import os

def ensure_directory_exists(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def process_files(input_path, output_dir):
    ensure_directory_exists(output_dir)
    # Your existing logic.
----------------------------------------------------

#to run above python script for translation

#!/bin/bash

input_dir="subject_db"
output_dir="new"

python ./translate.py "$input_dir" "$output_dir"

---------------------------------------------------------
#to merge all the translated sequnces in single fasta file.

#!/bin/bash

# Define the directory containing your FASTA files
directory="./aa_data"

# Define the output file where the merged FASTA will be saved
output_file="./aa_merged_x"

# Ensure the output file is empty
> "$output_file"

# Loop through all FASTA files in the specified directory, sorted alphabetically
for file in $(find "$directory" -name '*.fasta' | sort); do
    echo "Merging $file into $output_file"
    # Optionally, you can add a newline before each new file's content to ensure no run-ons from previous files
    echo "" >> "$output_file"
    cat "$file" >> "$output_file"
done

echo "All FASTA files have been merged into $output_file"

---------------------------------------------------------

#to replace the "X" in translated sequnces to "-"

from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

def replace_x_with_dash(fasta_file, output_file):
    records = list(SeqIO.parse(fasta_file, 'fasta'))  # Read all records from the FASTA file
    modified_records = []
    
    for record in records:
        # Replace 'X' with '-' directly in the sequence string and create a new Seq object
        modified_seq = str(record.seq).replace('X', '-')
        modified_record = SeqRecord(Seq(modified_seq), id=record.id, description=record.description)
        modified_records.append(modified_record)
    
    # Write the modified sequences back to a new file
    SeqIO.write(modified_records, output_file, 'fasta')

if __name__ == "__main__":
    fasta_file = './aa_merged_x'  # Path to the input FASTA file
    output_file = './finalaa.fasta'  # Path to the output FASTA file where modifications are saved
    replace_x_with_dash(fasta_file, output_file)

---------------------------------------------------------------

#to run above python script for translation

#!/bin/bash

input_dir="aa_merged_x"
output_dir="finalaa.fasta"

python ./xrep.py "$input_dir" "$output_dir"

-------------------------------------------------------------------

#to run blast on the trnaslated sequences against query sequences, using blast gompi.

#!/bin/bash

# Load the BLAST module
module load blast-uoneasy/2.14.1-gompi-2023a

# Define paths to your query and database
QUERY="./query_seqs/znf597_extracted.fa"
DATABASE="./mkdb/mkdb"  # Adjust 'finalaa' to your actual database prefix

# Define output file
OUTPUT="./blast_results/blast_znf597.txt"

# Run BLAST search
blastp -query $QUERY -db $DATABASE -out $OUTPUT -evalue 1e-7 -outfmt 6

echo "BLAST search has completed. Results are in $OUTPUT"

----------------------------------------------------------------------
