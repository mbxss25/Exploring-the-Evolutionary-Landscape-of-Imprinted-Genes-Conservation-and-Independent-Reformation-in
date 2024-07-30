## Evolution of Genomically Imprinted Gene Clusters in Mammals
Bioinformatics Individual Research Project (LIFE4137 UNUK) (SUM1 23-24)

NAME: Sahil Shelote

Student ID: 20595980

MSc Bioinformatics

School of Life Sciences

Supervisor: Mary O’connell

## **Overview** 
The repository created contains the data and scripts for my dissertation of MSc Bioinformatics, University of Nottingham. The project aimed at understanding the evolution of genomically imprinted gene clusters across the mammalian phylogenetic tree. The main goal is to investigate the conservation and divergence of imprinted gene clusters and explore the evolutionary rates of imprinted versus non-imprinted homologs.
The data of gene was sourced from the geneimprinting.com, this website confirmed the imprinting status of the gene in species. 

## Project Goals
**1**.Assemble Imprinted Gene Set:
collected a comprehensive set of imprinted genes from various mammalian species based on experimental evidence. Then Organized the data into a excel sheet  indicating imprinting status for each gene across species.

**2**.Extracting sequences and translation
Extracted query sequences for each imprinted gene and their homologs. This material was sourced from the Ensemble and biomart. For translation we used python package "biopython".

**3**.Finding homologous sequences: 
The BLAST was used  to findout the homologous sequences by programming a specific parameters. BLAST package used was "blast-uoneasy/2.14.1-gompi-2023a" module available on the ADA-HPC of university of Nottingham.

**4**.Sequence Alignment
The sequence alignment was accomplised with the help of MAFFT. The validation was done with the help of the matrix score generated while assembly.

**5**. Building phylogenetic trees
The trees were built with the help of the IQtree. Using IQtree allowed use to get standard outline of the phylogenetic trees of the speecies of the gene. We utilized only branch length present in the Newick format for the further analysis.

**6**.Vizualization of the phylogenetic trees
We find FigTree extremly helpful for vizualizing the data of the phylogenetic tree in high quality. The phylogenetic tree of each and every gene we studied are present on this github repository.

**7**.Statistical study
To determine which imprinted gene clusters are conserved and which have been disturbed we underwent few statistical studies. We Test the hypothesis if imprinted homologs evolve more rapidly than their non-imprinted counterparts. The part was done utilizing the R.

---------------------------------------------
## Acknowledgment 

I would like to express my sincere gratitude to all those who have supported and encouraged me throughout the course of my research and the writing of this thesis.First and foremost, I am deeply indebted to my supervisor, Mary O’connell, whose expertise, guidance, and patience were invaluable throughout this process. Her insightful feedback pushed me to sharpen my thinking and brought my work to a higher level.

I would also like to thank other staffs in laboratory of life sciences, James McInerney, Jonathan Fenn, Alan Beavan for their unwavering support and for always being available to discuss ideas and provide feedback.

I extend my gratitude to the staffs who helped us to resolve technical problems that occurred while undergoing this study. I would also like to acknowledge the High performace computing cluster of UoN i.e ada. and ChatGPT version 3, which helped me to rectify the code and setup this repository. 

---------------------------------------------
## Note

All the codes are wriiten in bash scripts. 

Output of each and every written code is present on the repository.

All the output data is submitted in .zip format,  Thus to view it,it should be downloaded.

The expected output for each submitted script will be tagged below it. 

---------------------------------------------
## Installation

Most  of the moducle utilized in this study were readily available on the Ada, thus there was not as such requirement of any specific installation. The module needsto be used were availed by executing code "module avail modulename".

Only MAFFT, IQtree and FigTree were the modules that were downloaded. The sources are mentioned below

MAFFT: conda install bioconda::mafft

IQtree: conda install -c bioconda iqtree

FigTree: http://tree.bio.ed.ac.uk/software/figtree/

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
#statistical analysis for comparing branch lengths within gene of different species. (chi-square test)

#Ampd3
  
  # Create a data frame with species and their branch lengths
  species <- c("asalmon", "catfish", "goldfish", "cow", "sheep", "rabbit", "Human", "macaque",
               "mouse", "rat", "pig", "opossum", "wallaby", "platypus", "platyfish", "sfightfish",
               "psolderfish", "imedaka", "mplatyfish", "platyfish_2", "sfightfish_2", "psolderfish_2")

branch_lengths <- c(0.3147155284, 0.1266752256, 0.2014479677, 0.0031934886, 0.0127031692, 0.0215495636,
                    0.0142581311, 0.00500846416, 0.0001841790, 0.0450726436, 0.0089775967, 0.0246567060,
                    0.0186123095, 0.0991088290, 0.1479684233, 0.3415666887, 0.0778627024, 0.1534971393,
                    0.0019947158, 0.0033699690, 0.2329132280, 0.0574678871)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#Asb4

# Create a data frame with species and their branch lengths
species <- c("catfish_1", "catfish_2", "cow", "sheep", "Human", "macaque", "rabbit", "mouse", "rat", "pig", 
             "dog", "hagfish", "reedfish", "asalmon", "imedaka", "jricefish", "mplatyfish", "platyfish", 
             "sfightfish", "psolderfish", "goldfish")

branch_lengths <- c(0.1931788056, 0.0581727292, 0.1891554181, 0.0044757454, 0.0000020219, 0.0106557828, 
                    0.0288790655, 0.0047532053, 0.0195811270, 0.0969453176, 0.0000020196, 1.2656155942, 
                    0.1748062969, 0.2777902469, 0.0311058711, 0.0270914320, 0.0060975198, 0.0116564547, 
                    0.1861661242, 0.1611768520, 0.2334015204)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#Ascl2

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "mouse", "rat", "Human", "wallaby", "macaque", "sheep")

branch_lengths <- c(0.0051102290, 0.0200226956, 0.0342312546, 0.0775635881, 0.0052534000, 
                    0.3947622722, 0.0123083594, 0.0217058550)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#Casd1

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "pig", 
             "sheep", "platypus", "wallaby", "reedfish", "asalmon", "imedaka", "jricefish", 
             "mplatyfish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.1851348947, 0.0027514162, 0.0217043375, 0.0045166718, 0.0052817374, 
                    0.0614995118, 0.0170436567, 0.0387399820, 0.0106050912, 0.0077293180, 
                    0.2055969372, 0.1566068830, 0.2609305303, 0.2822536152, 0.0082360888, 
                    0.0165311607, 0.0235498393, 0.0093565654, 0.1007548002, 0.1082811389, 
                    0.0689539943, 0.0547851827)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#cd81
# Create a data frame with species and their branch lengths
species <- c("catfish", "Human", "macaque", "wallaby", "rat", "pig", "sheep", "mouse")

branch_lengths <- c(1.0043826488, 0.0335098917, 0.1470924144, 0.2558775643, 0.0546522027, 
                    0.0738912833, 0.0716510339, 0.0000020899)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#cdh15
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "mouse", "rat", "opossum", 
             "wallaby", "platypus", "reedfish", "asalmon", "imedaka", "jricefish", "mplatyfish", 
             "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.2509165334, 0.0132993050, 0.0081423864, 0.0640112484, 0.0900800282, 
                    0.0149816539, 0.0490007644, 0.0116991378, 0.0325165102, 0.1769727597, 
                    0.0350557658, 0.1588057565, 0.5735342695, 0.3012832373, 0.0174192612, 
                    0.0275995115, 0.0450252652, 0.0142077159, 0.1602431030, 0.1215879011, 
                    0.0540397302, 0.0394150713)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#cdkn1c
# Create a data frame with species and their branch lengths
species <- c("cow", "Human_1", "Human_2", "macaque", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.1938272646, 0.0000027572, 0.0096730514, 0.0550106566, 
                    0.0489681178, 0.0000026236, 0.0620250762, 0.0072706210)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#commd1
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "rat", "pig", "Human", "macaque", "mouse", 
             "rabbit", "opossum", "wallaby", "platypus", "reedfish", "imedaka", "jricefish", 
             "mplatyfish", "platyfish", "psolderfish", "sfightfish", "goldfish_1", "goldfish_2", "goldfish_3")

branch_lengths <- c(0.1434924344, 0.0064036340, 0.2008323328, 0.1231680578, 0.2766745967, 
                    0.0168850978, 0.0057133563, 0.0054737796, 0.5669575617, 0.0346017916, 
                    0.0838978760, 0.0655284066, 0.1449379160, 0.2490364677, 0.0076581834, 
                    0.0090180629, 0.0054457508, 0.0000023046, 0.0812101115, 0.1504127491, 
                    0.0356161656, 0.1716286224, 0.1388049333)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#dcn
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit", "opossum", "wallaby", "platypus", 
             "asalmon", "imedaka", "jricefish", "psolderfish", "mplatyfish", "platyfish", 
             "sfightfish", "catfish", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4", 
             "reedfish", "pig", "sheep")

branch_lengths <- c(0.0152487971, 0.0521972698, 0.0233298584, 0.0334490431, 0.1570637106, 
                    0.2086109322, 0.1364125380, 0.1382093485, 0.5986987478, 0.1726030041, 
                    0.2116401333, 0.2220391274, 0.0000010830, 0.0000010830, 0.5879259668, 
                    0.3957401854, 0.0080024618, 0.0184104625, 0.0000010830, 0.0056724237, 
                    0.3987936391, 0.0425973991, 0.0076538467)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#dhcr7
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "Human", "macaque", "pig", "mouse", "rat", 
             "wallaby", "hagfish", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", 
             "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.1969023973, 0.0318870644, 0.0195957703, 0.0794013637, 0.0154399837, 
                    0.0190447704, 0.0553330437, 0.0169642277, 0.0376949159, 0.1634785768, 
                    0.5590368037, 0.1810968476, 0.0000028547, 0.2214952488, 0.1830180714, 
                    0.1458193336, 0.0634282189, 0.0591565570, 0.0975854895)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#dio3
# Create a data frame with species and their branch lengths
species <- c("cow", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4", "goldfish_5", "goldfish_6", "goldfish_7", "goldfish_8",
             "goldfish_9", "imedaka_1", "psolderfish", "imedaka_2", "wallaby", "Human", "rabbit", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0140382600, 0.0000000000, 0.0000000000, 0.0000023513, 0.0000020148, 
                    0.0000023513, 0.0000023513, 0.0000023513, 0.0152778019, 0.0166736659,
                    0.1207853984, 0.0600380674, 0.2655185674, 0.2557983962, 0.0036133745, 
                    0.7455425898, 0.0068906397, 0.0000023513, 0.0351084795, 0.0000023513)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#dlk1
# Create a data frame with species and their branch lengths
species <- c("asalmon", "catfish", "goldfish_1", "goldfish_2", "cow", "sheep", "pig", 
             "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "reedfish", 
             "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish")

branch_lengths <- c(0.4057310239, 0.3826195137, 0.0402946619, 0.0710395819, 0.0184243643, 
                    0.0234691899, 0.1176038598, 0.0444382782, 0.0275059345, 0.0621860475, 
                    0.1782185569, 0.0312280948, 0.0241757726, 0.3931132548, 0.5348869800, 
                    0.0182727458, 0.0109113360, 0.3020223949, 0.0000028229, 0.1198581614, 
                    0.1501241519)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)


#dscam
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "goldfish_1", "goldfish_2", "imedaka", "jricefish", "platyfish", "sfightfish", 
             "psolderfish", "opossum", "Human_1", "Human_2", "macaque", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0000010079, 0.0044170449, 0.0114301393, 0.0096477861, 0.0127407911, 
                    0.0098147807, 0.0249369016, 0.0158660176, 0.0149145804, 0.0269899974, 
                    0.0000010079, 0.0000010079, 0.0005101618, 0.0000010079, 0.0004977071, 
                    0.0065725829, 0.0014862782)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#gab1
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "rabbit", "mouse", "rat", 
             "Human", "macaque", "opossum", "wallaby", "platypus", "reedfish", 
             "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", 
             "goldfish_1", "goldfish_2")

branch_lengths <- c(0.1872110972, 0.0042481795, 0.0041316295, 0.0379862532, 0.0217405804, 
                    0.0435208401, 0.0379426528, 0.0537176696, 0.0058519801, 0.0050837760, 
                    0.0381628339, 0.0203924264, 0.0256492203, 0.2533123599, 0.0034330709, 
                    0.0337640006, 0.1296191896, 0.0508342445, 0.0432779283, 0.0732887035, 
                    0.0361420839)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#gatm
# Create a data frame with species and their branch lengths
species <- c("catfish_1", "catfish_2", "cow", "sheep", "dog", "opossum", "Human", "macaque", 
             "rabbit", "reedfish", "goldfish_1", "goldfish_2", "goldfish_3", "imedaka", 
             "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish", "catfish_3", "catfish_4")

branch_lengths <- c(0.0116263328, 0.0000010172, 0.0162133489, 0.0158547634, 0.0149086526, 
                    0.0368823341, 0.0088509362, 0.0171503419, 0.0488938241, 0.1130818122, 
                    0.0046721519, 0.0000010172, 0.0937337570, 0.0061363519, 0.0599774593, 
                    0.0000010172, 0.0000010172, 0.1142939837, 0.0508423918, 0.0000010172, 0.0044562039)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#glis3
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat",
             "opossum", "wallaby", "platypus", "imedaka", "jricefish", "mplatyfish", "platyfish",
             "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.3535785593, 0.0100974211, 0.0094605835, 0.0370199285, 0.0514208913, 
                    0.0329289929, 0.0090276525, 0.0404214959, 0.0210653072, 0.0274262935, 
                    0.0138939529, 0.0284853682, 0.0798407329, 0.0139139455, 0.0383598213, 
                    0.0704298497, 0.0078575163, 0.3130425228, 0.1804415493, 0.0025676831, 
                    0.0085933001)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#grb10
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "goldfish", "imedaka", "jricefish", "platyfish", "sfightfish", 
             "psolderfish", "hagfish", "opossum", "wallaby", "Human", "macaque", "rabbit", 
             "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0097308884, 0.0699152432, 0.0921884775, 0.0784976876, 0.0136011283, 
                    0.1085984878, 0.0595480677, 0.0228720876, 1.0950071459, 0.0163622940, 
                    0.0202290764, 0.0042553214, 0.0569503444, 0.1021281608, 0.0227580550, 
                    0.0293922435, 0.0505513568, 0.0188727870)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#htr2a
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "opossum", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4",
             "goldfish_5", "goldfish_6", "imedaka", "jricefish", "platyfish", "sfightfish",
             "psolderfish", "Human", "macaque", "mouse", "rat", "rabbit", "pig", "sheep")

branch_lengths <- c(0.0071940736, 0.0301051387, 0.2632995537, 0.0326559751, 0.1235121393, 
                    0.0380213620, 0.0042245463, 0.1312696812, 0.5056138548, 0.1362131867,
                    0.0034831404, 0.2167364877, 0.2292252639, 0.0882849903, 0.0000025144,
                    0.0076212171, 0.0096377885, 0.0204047070, 0.0508567287, 0.0285034451,
                    0.0247523659)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#igf2
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "Human", "macaque", "wallaby", "mouse", "rat", "sheep")

branch_lengths <- c(5.2499771172, 0.0000023631, 0.2501200858, 0.0653180309, 0.0000022539, 
                    1.0064636886, 0.0493037101, 0.0289216618, 0.0729959662)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#igf2r
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", 
             "mouse", "rat", "opossum", "platypus", "reedfish", "asalmon", "imedaka", 
             "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.2908769061, 0.0236149481, 0.0281474810, 0.0864078277, 0.0776954217, 
                    0.0217273885, 0.0300121360, 0.2009358115, 0.0328980534, 0.0408796312, 
                    0.1759768633, 0.2009510702, 0.4308277537, 0.3209193880, 0.0503666143, 
                    0.0402080160, 0.2429290927, 0.2083643943, 0.1700630678, 0.1217247684, 
                    0.0936552527)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#impact
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "pig", "Human", "macaque", "rabbit", 
             "mouse", "rat", "opossum", "reedfish", "imedaka", "jricefish", 
             "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.3645343694, 0.0673480048, 0.0455625425, 0.0544064007, 
                    0.0040610938, 0.0000019972, 0.2821983197, 0.1092647689, 
                    0.1001771118, 0.1447400460, 0.4472904332, 0.2440556192, 
                    0.1359002442, 0.1586463997, 0.1700433487, 0.0841683299, 
                    0.0166098663, 0.0089697999)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#lin28b
# Create a data frame with species and their branch lengths
species <- c("catfish_1", "catfish_2", "cow", "pig", "sheep", "dog", "Human", "mouse", 
             "rat", "macaque", "rabbit", "opossum", "wallaby", "reedfish", "imedaka", 
             "jricefish", "platyfish", "psolderfish", "sfightfish", "goldfish_1", 
             "goldfish_2", "goldfish_3", "goldfish_4")

branch_lengths <- c(0.2099123774, 0.5742520224, 0.0988597469, 0.0047952992, 0.8534256275, 
                    0.0455625425, 0.0022102646, 0.0284116869, 0.0414143863, 0.0305357474, 
                    0.0660834331, 0.0097311564, 0.0954293807, 0.4970692898, 0.0395855195, 
                    0.0360508139, 0.0514808765, 0.0529380138, 0.0492276856, 0.0078481891, 
                    0.0065363268, 0.0144290125, 0.0114238039)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#lrrtm1
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "Human", "macaque", "rat", "mouse", 
             "rabbit", "dog", "opossum", "platypus", "reedfish", "imedaka", "jricefish", 
             "platyfish", "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.0737786785, 0.0021673003, 0.0000025938, 0.0063561514, 0.0021678522, 
                    0.0021720593, 0.0115169572, 0.0042691969, 0.0822475890, 0.0125564962, 
                    0.0459871771, 0.0824201214, 0.0891323712, 0.0076169050, 0.0010838226, 
                    0.0422891854, 0.0242981219, 0.0205629284, 0.0141519619)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#naa60
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "pig", "rabbit", "mouse", 
             "macaque", "rat", "opossum", "sheep")

branch_lengths <- c(0.0000020286, 0.1536512310, 0.0203021797, 0.0044994776, 0.0078637490, 
                    0.0041069618, 0.7456022334, 0.1040223800, 0.0479091230, 0.0082257527)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#osbpl5
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", 
             "mouse", "rat", "opossum", "wallaby", "reedfish", "imedaka", 
             "jricefish", "platyfish", "sfightfish", "psolderfish", 
             "goldfish_1", "goldfish_2")

branch_lengths <- c(0.1458124884, 0.0057729919, 0.0168204792, 0.0388053212, 
                    0.0591310361, 0.0076157521, 0.0206628626, 0.0180450198, 
                    0.0161956445, 0.0585569429, 0.0599842536, 0.4413103892, 
                    0.0205741695, 0.0152159753, 0.1029129909, 0.0582269026, 
                    0.1816927301, 0.0348310027, 0.0786396244)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#plagl1
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "hagfish", "opossum", "Human", "macaque", 
             "rabbit", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0194626249, 0.0191125805, 1.1371154170, 0.2021487465, 
                    0.0000029958, 0.0110204225, 0.0618356602, 0.0950481189, 
                    0.0879394191, 0.0068137620, 0.0383106375)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#pon2
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "opossum", "wallaby", 
             "rabbit", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0055585162, 0.1605448571, 0.0092445340, 0.0132764349, 
                    0.2413293314, 0.0697122226, 0.0372865997, 0.0323480216, 
                    0.0149859150, 0.0190708493, 0.0143350995)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#pon3
# Create a data frame with species and their branch lengths
species <- c("cow", "Human", "macaque", "rabbit", "mouse", "rat", "pig", "sheep")

branch_lengths <- c(0.0797574974, 0.0177860978, 0.0176633739, 0.1974271084, 
                    0.0465122882, 0.1177570492, 0.1023570242, 0.0258231285)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#rb1
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "pig", "Human", "macaque", 
             "rabbit", "mouse", "rat", "platypus", "opossum", "wallaby", 
             "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", 
             "psolderfish", "goldfish_1", "goldfish_2")

branch_lengths <- c(0.2055812757, 0.0081434574, 0.0021595416, 0.0332168073, 
                    0.0196866621, 0.0066102568, 0.0000027648, 0.0297671795, 
                    0.0186874015, 0.0324412708, 0.1023786286, 0.0521032423, 
                    0.0305978487, 0.3856571933, 0.0213174560, 0.0170554119, 
                    0.1354484982, 0.1206411549, 0.1489856306, 0.0751267468, 
                    0.0240459072)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#sfmbt2
# Create a data frame with species and their branch lengths
species <- c("asalmon", "catfish", "goldfish_1", "goldfish_2", "cow", "sheep", "pig", "dog", 
             "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "platypus", "reedfish", 
             "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish")

branch_lengths <- c(0.3893371610, 0.2038967584, 0.0488864476, 0.0268946144, 0.0015119353, 
                    0.0147962157, 0.0634928125, 0.0430449771, 0.0091752679, 0.0142262131, 
                    0.0250512057, 0.1195188937, 0.1190891183, 0.0611817190, 0.0423816106, 
                    0.1418717385, 0.0149312441, 0.0295483094, 0.1616777144, 0.0000024080, 
                    0.1093360667, 0.1054227938)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#slc38a4
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "pig", "Human", "macaque", "mouse", "rat", 
             "rabbit", "opossum", "reedfish", "imedaka", "jricefish", "platyfish", 
             "psolderfish", "sfightfish", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4")

branch_lengths <- c(0.1192941451, 0.0069283680, 0.0164346867, 0.0753295420, 0.0529297067, 
                    0.0190784196, 0.0115024123, 0.0588143047, 0.0439372516, 0.2304188869, 
                    0.0937383792, 0.1607310087, 0.0222129042, 0.0188338567, 0.1495928077, 
                    0.1035944995, 0.0995584317, 0.0000010158, 0.0042735032, 0.0043250523, 
                    0.0044021198)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#Tssc4
# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "mouse", "rat", 
             "reedfish", "asalmon", "goldfish")

branch_lengths <- c(0.2412293671, 0.0650635270, 0.0531560164, 0.3414039653, 0.1929408882, 
                    0.0039909191, 0.0298890041, 0.0552319188, 0.0402110501, 0.7353518874, 
                    0.9254881734, 0.5227650953)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#zc3h12c
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "wallaby", 
             "platypus", "reedfish", "pig", "sheep")

branch_lengths <- c(0.0093644817, 0.0770295289, 0.0122680133, 0.0172141740, 0.0377101396, 
                    0.0188618148, 0.0162480470, 0.0228034250, 0.0292569244, 0.1484013199, 
                    0.5536163764, 0.0241887663, 0.0063339392)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#zfat
# Create a data frame with species and their branch lengths
species <- c("cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", 
             "wallaby", "reedfish", "goldfish1", "goldfish2", "goldfish3")

branch_lengths <- c(0.0085967938, 0.0240774901, 0.0540641387, 0.0900411568, 0.0082355482, 0.0216130280, 
                    0.1848942958, 0.0192151115, 0.0185209272, 0.0632621309, 0.0354126886, 0.3890100043, 
                    0.0016577287, 0.0803823496, 0.0670357840)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

#znf597
# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit")

branch_lengths <- c(0.3889007220, 0.1635983387, 0.0307547632, 0.0315059287, 0.2078778802)

data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Calculate the mean and standard deviation of branch lengths
mean_branch_length <- mean(data$Branch_Length)
sd_branch_length <- sd(data$Branch_Length)

# Categorize branch lengths into "Short", "Average", and "Long"
threshold_long <- mean_branch_length + sd_branch_length
threshold_short <- mean_branch_length - sd_branch_length

data$Category <- cut(data$Branch_Length, breaks = c(-Inf, threshold_short, threshold_long, Inf), labels = c("Short", "Average", "Long"))

# Observed frequencies
observed <- table(data$Category)

# Expected frequencies assuming uniform distribution
expected <- rep(length(branch_lengths) / length(observed), length(observed))

# Perform chi-square test
chi_square_test <- chisq.test(observed, p = expected / sum(expected))

# Print the results
print(data)
print(chi_square_test)

-------------------

#stastical analysis of imprinted/nonimprinted vs others.


```{r}
#Ampd3
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("asalmon", "catfish", "goldfish", "cow", "sheep", "rabbit", "Human", "macaque", "mouse", "rat", "pig", "opossum", "wallaby", "platypus", "platyfish", "sfightfish", "psolderfish", "imedaka", "mplatyfish", "platyfish_2", "sfightfish_2", "psolderfish_2")
branch_lengths <- c(0.3147155284, 0.1266752256, 0.2014479677, 0.0031934886, 0.0127031692, 0.0215495636, 0.0142581311, 0.00500846416, 0.0001841790, 0.0450726436, 0.0089775967, 0.0246567060, 0.0186123095, 0.0991088290, 0.1479684233, 0.3415666887, 0.0778627024, 0.1534971393, 0.0019947158, 0.0033699690, 0.2329132280, 0.0574678871)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#Asb4
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish_1", "catfish_2", "cow", "sheep", "Human", "macaque", "rabbit", "mouse", "rat", "pig", "dog", "hagfish", "reedfish", "asalmon", "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish", "goldfish")
branch_lengths <- c(0.1931788056, 0.0581727292, 0.1891554181, 0.0044757454, 0.0000020219, 0.0106557828, 0.0288790655, 0.0047532053, 0.0195811270, 0.0969453176, 0.0000020196, 1.2656155942, 0.1748062969, 0.2777902469, 0.0311058711, 0.0270914320, 0.0060975198, 0.0116564547, 0.1861661242, 0.1611768520, 0.2334015204)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "Human", "mouse", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#ascl4
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "mouse", "rat", "Human", "wallaby", "macaque", "sheep")
branch_lengths <- c(0.0051102290, 0.0200226956, 0.0342312546, 0.0775635881, 0.0052534000, 0.3947622722, 0.0123083594, 0.0217058550)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "mouse", "Human") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#casd1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "pig", "sheep", "platypus", "wallaby", "reedfish", "asalmon", "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.1851348947, 0.0027514162, 0.0217043375, 0.0045166718, 0.0052817374, 0.0614995118, 0.0170436567, 0.0387399820, 0.0106050912, 0.0077293180, 0.2055969372, 0.1566068830, 0.2609305303, 0.2822536152, 0.0082360888, 0.0165311607, 0.0235498393, 0.0093565654, 0.1007548002, 0.1082811389, 0.0689539943, 0.0547851827)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("pig", "sheep") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#cd81
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "Human", "macaque", "wallaby", "rat", "pig", "sheep", "mouse")
branch_lengths <- c(1.0043826488, 0.0335098917, 0.1470924144, 0.2558775643, 0.0546522027, 0.0738912833, 0.0716510339, 0.0000020899)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "pig", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#cdh15
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "mouse", "rat", "opossum", "wallaby", "platypus", "reedfish", "asalmon", "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.2509165334, 0.0132993050, 0.0081423864, 0.0640112484, 0.0900800282, 0.0149816539, 0.0490007644, 0.0116991378, 0.0325165102, 0.1769727597, 0.0350557658, 0.1588057565, 0.5735342695, 0.3012832373, 0.0174192612, 0.0275995115, 0.0450252652, 0.0142077159, 0.1602431030, 0.1215879011, 0.0540397302, 0.0394150713)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#cdkn1c
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "Human_1", "Human_2", "macaque", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.1938272646, 0.0000027572, 0.0096730514, 0.0550106566, 0.0489681178, 0.0000026236, 0.0620250762, 0.0072706210)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#commd1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "rat", "pig", "Human", "macaque", "mouse", "rabbit", "opossum", "wallaby", "platypus", "reedfish", "imedaka", "jricefish", "mplatyfish", "platyfish", "psolderfish", "sfightfish", "goldfish_1", "goldfish_2", "goldfish_3")
branch_lengths <- c(0.1434924344, 0.0064036340, 0.2008323328, 0.1231680578, 0.2766745967, 0.0168850978, 0.0057133563, 0.0054737796, 0.5669575617, 0.0346017916, 0.0838978760, 0.0655284066, 0.1449379160, 0.2490364677, 0.0076581834, 0.0090180629, 0.0054457508, 0.0000023046, 0.0812101115, 0.1504127491, 0.0356161656, 0.1716286224, 0.1388049333)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#dcn
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit", "opossum", "wallaby", "platypus", "asalmon", "imedaka", "jricefish", "psolderfish", "mplatyfish", "platyfish", "sfightfish", "catfish", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4", "reedfish", "pig", "sheep")
branch_lengths <- c(0.0152487971, 0.0521972698, 0.0233298584, 0.0334490431, 0.1570637106, 0.2086109322, 0.1364125380, 0.1382093485, 0.5986987478, 0.1726030041, 0.2116401333, 0.2220391274, 0.0000010830, 0.0000010830, 0.5879259668, 0.3957401854, 0.0080024618, 0.0184104625, 0.0000010830, 0.0056724237, 0.3987936391, 0.0425973991, 0.0076538467)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#dhcr7
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "Human", "macaque", "pig", "mouse", "rat", "wallaby", "hagfish", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.1969023973, 0.0318870644, 0.0195957703, 0.0794013637, 0.0154399837, 0.0190447704, 0.0553330437, 0.0169642277, 0.0376949159, 0.1634785768, 0.5590368037, 0.1810968476, 0.0000028547, 0.2214952488, 0.1830180714, 0.1458193336, 0.0634282189, 0.0591565570, 0.0975854895)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```



```{r}
#dio3
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4", "goldfish_5", "goldfish_6", "goldfish_7", "goldfish_8", "goldfish_9", "goldfish_10", "imedaka_1", "psolderfish", "imedaka_2", "wallaby", "Human", "rabbit", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.0140382600, 0.0000000000, 0.0000000000, 0.0000023513, 0.0000023513, 0.0000020148, 0.0000023513, 0.0000023513, 0.0152778019, 0.0166736659, 0.0889459754, 0.1207853984, 0.0600380674, 0.2655185674, 0.2557983962, 0.0036133745, 0.7455425898, 0.0068906397, 0.0000023513, 0.0351084795, 0.0000023513)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "wallaby", "Human", "mouse", "rat", "sheep", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#dlk1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("asalmon", "catfish", "goldfish_1", "goldfish_2", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "reedfish", "imedaka_1", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish")
branch_lengths <- c(0.4057310239, 0.3826195137, 0.0402946619, 0.0710395819, 0.0184243643, 0.0234691899, 0.1176038598, 0.0444382782, 0.0275059345, 0.0621860475, 0.1782185569, 0.0312280948, 0.0241757726, 0.3931132548, 0.5348869800, 0.0182727458, 0.0109113360, 0.3020223949, 0.0000028229, 0.1198581614, 0.1501241519)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "sheep", "pig", "Human", "macaque", "mouse", "rat", "opossum") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#dlx5
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "pig", "hagfish", "opossum", "wallaby", "Human", "macaque", "mouse", "rat", "rabbit", "sheep")
branch_lengths <- c(0.0000023732, 0.0069602623, 0.0104846260, 0.9348232823, 0.0355983414, 0.0427534505, 0.0033832821, 0.0069798239, 0.0000023732, 0.0104613447, 0.4645074289, 0.0000023732)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("pig", "Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#dscam
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "goldfish1", "goldfish2", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "opossum", "Human1", "Human2", "macaque", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.0000010079, 0.0044170449, 0.0114301393, 0.0096477861, 0.0127407911, 0.0098147807, 0.0249369016, 0.0158660176, 0.0149145804, 0.0269899974, 0.0000010079, 0.0000010079, 0.0005101618, 0.0000010079, 0.0004977071, 0.0065725829, 0.0014862782)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human1", "Human2", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#gab1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "rabbit", "mouse", "rat", "Human", "macaque", "opossum", "wallaby", "platypus", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish1", "goldfish2")
branch_lengths <- c(0.1872110972, 0.0042481795, 0.0041316295, 0.0379862532, 0.0217405804, 0.0435208401, 0.0379426528, 0.0537176696, 0.0058519801, 0.0050837760, 0.0381628339, 0.0203924264, 0.0256492203, 0.2533123599, 0.0034330709, 0.0337640006, 0.1296191896, 0.0508342445, 0.0432779283, 0.0732887035, 0.0361420839)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```



```{r}
#glis3
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "wallaby", "platypus", "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish", "goldfish1", "goldfish2")
branch_lengths <- c(0.3535785593, 0.0100974211, 0.0094605835, 0.0370199285, 0.0514208913, 0.0329289929, 0.0090276525, 0.0404214959, 0.0210653072, 0.0274262935, 0.0138939529, 0.0284853682, 0.0798407329, 0.0139139455, 0.0383598213, 0.0704298497, 0.0078575163, 0.3130425228, 0.1804415493, 0.0025676831, 0.0085933001)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)


```
```{r}
#grb10
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "goldfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "hagfish", "opossum", "wallaby", "Human", "macaque", "rabbit", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.0097308884, 0.0699152432, 0.0921884775, 0.0784976876, 0.0136011283, 0.1085984878, 0.0595480677, 0.0228720876, 1.0950071459, 0.0163622940, 0.0202290764, 0.0042553214, 0.0569503444, 0.1021281608, 0.0227580550, 0.0293922435, 0.0505513568, 0.0188727870)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("wallaby", "Human", "mouse", "sheep") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#htr2a
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "opossum", "goldfish1", "goldfish2", "goldfish3", "goldfish4", "goldfish5", "goldfish6", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "Human", "macaque", "mouse", "rat", "rabbit", "pig", "sheep")
branch_lengths <- c(0.0071940736, 0.0301051387, 0.2632995537, 0.0326559751, 0.1235121393, 0.0380213620, 0.0042245463, 0.1312696812, 0.0206405066, 0.1362131867, 0.0034831404, 0.2167364877, 0.2292252639, 0.0882849903, 0.0000025144, 0.0076212171, 0.0096377885, 0.0204047070, 0.0508567287, 0.0285034451, 0.0247523659)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("cow", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#igf2
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "Human", "macaque", "wallaby", "mouse", "rat", "sheep")
branch_lengths <- c(5.2499771172, 0.0000023631, 0.2501200858, 0.0653180309, 0.0000022539, 1.0064636886, 0.0493037101, 0.0289216618, 0.0729959662)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("catfish", "Human", "wallaby", "mouse", "rat") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#igf2r
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "platypus", "reedfish", "asalmon", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.2908769061, 0.0236149481, 0.0281474810, 0.0864078277, 0.0776954217, 0.0217273885, 0.0300121360, 0.2009358115, 0.0328980534, 0.0408796312, 0.1759768633, 0.2009510702, 0.4308277537, 0.3209193880, 0.0503666143, 0.0402080160, 0.2429290927, 0.2083643943, 0.1700630678, 0.1217247684, 0.0936552527)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("catfish", "Human", "rat", "cow", "sheep", "pig", "dog", "mouse", "rat", "opossum", "platypus", "reedfish", "asalmon", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```


```{r}
#impact
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "dog", "pig", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.3645343694, 0.0673480048, 0.0455625425, 0.0544064007, 0.0040610938, 0.0000019972, 0.2821983197, 0.1092647689, 0.1001771118, 0.1447400460, 0.4472904332, 0.2440556192, 0.1359002442, 0.1586463997, 0.1700433487, 0.0841683299, 0.0166098663, 0.0089697999)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("pig", "Human", "rabbit", "mouse", "rat") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#lin28b
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish_1", "catfish_2", "cow", "pig", "sheep", "dog", "Human", "mouse", "rat", "macaque", "rabbit", "opossum", "wallaby", "reedfish", "imedaka", "jricefish", "platyfish", "psolderfish", "sfightfish", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4")
branch_lengths <- c(0.2099123774, 0.5742520224, 0.0988597469, 0.0047952992, 0.8534256275, 0.1058963644, 0.0022102646, 0.0284116869, 0.0414143863, 0.0305357474, 0.0660834331, 0.0097311564, 0.0954293807, 0.4970692898, 0.0395855195, 0.0360508139, 0.0514808765, 0.0529380138, 0.0492276856, 0.0078481891, 0.0065363268, 0.0144290125, 0.0114238039)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#lrrtm1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "Human", "macaque", "rat", "mouse", "rabbit", "dog", "opossum", "platypus", "reedfish", "imedaka", "jricefish", "platyfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.0737786785, 0.0021673003, 0.0000025938, 0.0063561514, 0.0021678522, 0.0021720593, 0.0115169572, 0.0042691969, 0.0822475890, 0.0125564962, 0.0459871771, 0.0824201214, 0.0891323712, 0.0076169050, 0.0010838226, 0.0422891854, 0.0242981219, 0.0205629284, 0.0141519619)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "macaque", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#naa60


## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "pig", "rabbit", "mouse", "macaque", "rat", "opossum", "sheep")
branch_lengths <- c(0.0000020286, 0.1536512310, 0.0203021797, 0.0044994776, 0.0078637490, 0.0041069618, 0.7456022334, 0.1040223800, 0.0479091230, 0.0082257527)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#osbpl5
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "mouse", "rat", "opossum", "wallaby", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.1458124884, 0.0057729919, 0.0168204792, 0.0388053212, 0.0591310361, 0.0076157521, 0.0206628626, 0.0180450198, 0.0161956445, 0.0585569429, 0.0599842536, 0.4413103892, 0.0205741695, 0.0152159753, 0.1029129909, 0.0582269026, 0.1816927301, 0.0348310027, 0.0786396244)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```


```{r}
#pon2
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "opossum", "wallaby", "rabbit", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.0055585162, 0.1605448571, 0.0092445340, 0.0132764349, 0.2413293314, 0.0697122226, 0.0372865997, 0.0323480216, 0.0149859150, 0.0190708493, 0.0143350995)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "pig") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#pon3
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "Human", "macaque", "rabbit", "mouse", "rat", "pig", "sheep")
branch_lengths <- c(0.0797574974, 0.0177860978, 0.0176633739, 0.1974271084, 0.0465122882, 0.1177570492, 0.1023570242, 0.0258231285)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "sheep") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#rb1
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "sheep", "dog", "pig", "Human", "macaque", "rabbit", "mouse", "rat", "platypus", "opossum", "wallaby", "reedfish", "imedaka", "jricefish", "platyfish", "sfightfish", "psolderfish", "goldfish_1", "goldfish_2")
branch_lengths <- c(0.0081434574, 0.0021595416, 0.0332168073, 0.0196866621, 0.0066102568, 0.0000027648, 0.0297671795, 0.0186874015, 0.0324412708, 0.1023786286, 0.0521032423, 0.0305978487, 0.3856571933, 0.0213174560, 0.0170554119, 0.1354484982, 0.1206411549, 0.1489856306, 0.0751267468, 0.0240459072)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#sfmbt2
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("asalmon", "catfish", "goldfish_1", "goldfish_2", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "platypus", "reedfish", "imedaka", "jricefish", "mplatyfish", "platyfish", "sfightfish", "psolderfish")
branch_lengths <- c(0.3893371610, 0.2038967584, 0.0488864476, 0.0268946144, 0.0015119353, 0.0147962157, 0.0634928125, 0.0430449771, 0.0091752679, 0.0142262131, 0.0250512057, 0.1195188937, 0.1190891183, 0.0611817190, 0.0423816106, 0.1418717385, 0.0149312441, 0.0295483094, 0.1616777144, 0.0000024080, 0.1093360667, 0.1054227938)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse", "rat", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#slc38a4
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "dog", "pig", "Human", "macaque", "mouse", "rat", "rabbit", "opossum", "reedfish", "imedaka", "jricefish", "platyfish", "psolderfish", "sfightfish", "goldfish_1", "goldfish_2", "goldfish_3", "goldfish_4")
branch_lengths <- c(0.1192941451, 0.0069283680, 0.0164346867, 0.0753295420, 0.0529297067, 0.0190784196, 0.0115024123, 0.0588143047, 0.0439372516, 0.2304188869, 0.0937383792, 0.1607310087, 0.0222129042, 0.0188338567, 0.1495928077, 0.1035944995, 0.0995584317, 0.0000010158, 0.0042735032, 0.0043250523, 0.0044021198)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("pig", "mouse", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```


```{r}
#tssc4
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish", "cow", "sheep", "pig", "dog", "Human", "macaque", "mouse", "rat", "reedfish", "asalmon", "goldfish")
branch_lengths <- c(0.2412293671, 0.0650635270, 0.0531560164, 0.3414039653, 0.1929408882, 0.0039909191, 0.0298890041, 0.0552319188, 0.0402110501, 0.7353518874, 0.9254881734, 0.5227650953)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#zc3h12c
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "wallaby", "platypus", "reedfish", "pig", "sheep")
branch_lengths <- c(0.0093644817, 0.0770295289, 0.0122680133, 0.0172141740, 0.0377101396, 0.0188618148, 0.0162480470, 0.0228034250, 0.0292569244, 0.1484013199, 0.5536163764, 0.0241887663, 0.0063339392)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#zfat
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("catfish1", "catfish2", "cow", "sheep", "pig", "dog", "Human", "macaque", "rabbit", "mouse", "rat", "opossum", "wallaby", "reedfish", "goldfish1", "goldfish2", "goldfish3")
branch_lengths <- c(0.3021483855, 9.9999988669, 0.0085967938, 0.0240774901, 0.0540641387, 0.0900411568, 0.0082355482, 0.0216130280, 0.1848942958, 0.0192151115, 0.0185209272, 0.0632621309, 0.0354126886, 0.3890100043, 0.0016577287, 0.0803823496, 0.0670357840)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "mouse") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```
```{r}
#znf597
## Load necessary library
library(dplyr)

# Create a data frame with species and their branch lengths
species <- c("cow", "dog", "Human", "macaque", "rabbit")
branch_lengths <- c(0.3889007220, 0.1635983387, 0.0307547632, 0.0315059287, 0.2078778802)

# Create a data frame
data <- data.frame(Species = species, Branch_Length = branch_lengths)

# Assign groups
data <- data %>%
  mutate(Category = case_when(
    Species %in% c("Human", "cow") ~ "Imprinted/non-imprinted",
    TRUE ~ "Others"
  ))

# Perform ANOVA
anova_result <- aov(Branch_Length ~ Category, data = data)

# Print the results
summary(anova_result)

# If ANOVA is significant, perform post-hoc pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Function for pairwise t-tests
compare_groups <- function(group1, group2) {
  data1 <- data[data$Category == group1, ]
  data2 <- data[data$Category == group2, ]
  t.test(data1$Branch_Length, data2$Branch_Length)
}

# Perform t-tests for pairwise comparisons
imprinted_vs_others <- compare_groups("Imprinted/non-imprinted", "Others")

# Print t-test results
print(imprinted_vs_others)

```

