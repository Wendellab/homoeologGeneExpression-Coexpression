# HyLiTE documentation: http://hylite.sourceforge.net

# log in to server biocrunch
mkdir bowtie2hylite
cd bowtie2hylite

# First, prepare D5 transcript fasta as ref gene sequences
## get genome
cp ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.fasta . 
## get primary transcript annotation
cp ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf .
## check to make sure there are the right number genes 
sed '/Gorai[.]N/d' D5.primaryOnly.gtf | cut -f9 | sort | uniq |wc -l
## extract primary transcript sequences from genome
module load rsem
rsem-prepare-reference --gtf D5.primaryOnly.gtf --bowtie2 Dgenome2_13.fasta D5
# D5.transcripts.fa contains all extracted transcript sequences from "exon" features, and *.bt2 are Bowtie2 indices.

# Next, prepare fq files
mkdir fastq
cd fastq
cp ~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/eflen_seed/*cut.fq.gz .
## uncompress all gz files to fq
gunzip -c A2-10-R1.cut.fq.gz > A2-10-R1.cut.fq

# Run bowtie2 mapping
module load python/2
bowtie2-build D5.transcripts.fa D5.transcripts
mkdir sam
bowtie2 -x D5.transcripts -U fastq/A2-10-R1.cut.fq -S sam/A2-10-R1.sam -N 1 -p 8

# Run HyLite, note that bowtie2 2.7 is preferred
HyLiTE -v -S -f sam2_protocol_file.txt -r D5.transcripts.fa -n results030917 >cmd0309.log 2>&1
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
