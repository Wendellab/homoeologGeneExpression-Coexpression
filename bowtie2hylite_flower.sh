# HyLiTE documentation: http://hylite.sourceforge.net

# log in to server biocrunch
mkdir mapping_hylite

# locate reference
ls ~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2hylite/

# Run bowtie2 mapping
module load python/2
mkdir mapping_hylite/sam
bash runBowtie2.sh >mapping_hylite/runBowtie2_073117.txt

# Run HyLite, note that bowtie2 2.7 is preferred
HyLiTE -v -S -f spf_flower.txt -r ~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2hylite/D5.transcripts.fa -n mapping_hylite/results0817 >mapping_hylite/cmd00817.log 2>&1
## Options as:
## -v turns on verbose runtime comments
## -f specifies the protocol file
## -r specifies the .fasta file containing the reference gene sequences
## -n allows the user to provide a name for the HyLiTE analysis, and creates a directory for output files
## -S use mapping .sam results instead of .fastq from protocol file
## -b use pre built ref 

# Inspect resutls
# HyLiTE output files in folder "results030917"
# Explanation of HyLiTE output formats: https://hylite.sourceforge.io/outformat.html#outformat
# --- total read counts: "results030917.expression.txt"
# --- allele specific counts: "results030917.AD.AD-10-R1.read.txt", etc.
