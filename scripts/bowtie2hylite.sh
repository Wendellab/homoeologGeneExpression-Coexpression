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

# ---------- Run bowtie2 mapping to obtain SAM files--------------
# ----------   code skipped --------------------------------------

## 2018-1-12 HyLiTE has now been updated and rewritten in Python 3
# module load python/3.6.3-u4oaxsb samtools/1.9-k6deoga py-scipy/1.1.0-py3-pzig4lr
# python3 HyLiTE-2.0.1-py3-none-any.whl --help
HyLiTE -v -S -f sam2_protocol_file.txt -r D5.transcripts.fa -n resultsXXXXXX >cmdXXXXXX.log 2>&1
## Options as:
## -v turns on verbose runtime comments
## -f specifies the protocol file
## -r specifies the .fasta file containing the reference gene sequences
## -n allows the user to provide a name for the HyLiTE analysis, and creates a directory for output files
## -S use mapping .sam results instead of .fastq from protocol file
## -b use pre built ref 

# Inspect resutls
# HyLiTE output files in folder "resultsXXXXXX"
# Explanation of HyLiTE output formats: https://hylite.sourceforge.io/outformat.html#outformat
# --- total read counts: "resultsXXXXXX.expression.txt"
# --- allele specific counts: "resultsXXXXXX.AD.AD-10-R1.read.txt", etc.
