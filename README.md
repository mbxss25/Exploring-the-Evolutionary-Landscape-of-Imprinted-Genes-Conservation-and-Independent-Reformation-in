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

All the codes are wriiten in bash scripts, Thus chmod +x filename.sh is used for each an every script in sbatch to make it executable. 

Python scripts are saved in .py format, however are executed in sbatch format.

Output of each and every written code is present on the repository.

All the output data is submitted in .zip format,  Thus to view it,it should be downloaded.

The expected output for each submitted script will be tagged below it will file number. 

---------------------------------------------
## Installation

Most  of the moducle utilized in this study were readily available on the Ada, thus there was not as such requirement of any specific installation. The module needsto be used were availed by executing code "module avail modulename".

Only MAFFT, IQtree and FigTree were the modules that were downloaded. The sources are mentioned below

MAFFT: conda install bioconda::mafft

IQtree: conda install -c bioconda iqtree

FigTree: http://tree.bio.ed.ac.uk/software/figtree/

---------------------------------------------
## Uploaded files 

1. selected_genes.zip
2. query_seqs.zip
3_subject_dborg.zip
4_subject_db.zip
5_aa_data.zip
6_merged_all_aminoacids.fasta.zip
7. merged_clusters.zip
8. new_blastmkdata.zip
9. final_blast.zip
11. ready_clusters (homologous seqs).zip
12. aligned_sequences.zip
13. output_trees.zip
15. final_scripts.zip
branch_lenght_data.docx
chi-sq_test_data.csv
figtree_vizuals.zip
others_vs_impnonimp_ANOVA.csv
gene_set

