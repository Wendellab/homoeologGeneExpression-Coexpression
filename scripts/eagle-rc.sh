## install EAGLE https://github.com/tony-kuo/eagle/
cd /work/LAS/jfw-lab/hugj2006/tools/
module purge
module load gcc git
git clone https://github.com/tony-kuo/eagle.git
cd eagle
git clone https://github.com/samtools/htslib.git
make

# code modified from https://github.com/tony-kuo/eagle/blob/master/scripts/tetraploid_analysis.sh
# Example workflow using EAGLE-RC for tetraploid Arabidopsis kamchatica (Arabidopsis halleri + Arabidopsis lyrata)

cd /work/LAS/jfw-lab/hugj2006/eflen2020/eagle

# gffread at: https://github.com/gpertea/gffread.
# LAST at: http://last.cbrc.jp/
# STAR at: https://github.com/alexdobin/STAR
# featureCounts at: http://bioinf.wehi.edu.au/featureCounts/
module load last/869-56gezob star/2.5.3a-rdzqclb samtools/1.9-k6deoga

## Homeolog identification

# Genome sequences
ln -s ../bowtie2hylite/A2Du_13.fasta
ln -s ../bowtie2hylite/Dgenome2_13.fasta

# Transcript sequences
ln -s ../bowtie2hylite/A2.transcripts.fa
ln -s ../bowtie2hylite/D5.transcripts.fa

# GTF
ln -s ../bowtie2hylite/A2.gtf
ln -s ../bowtie2hylite/D5.primaryOnly.gtf

# Reciprocal best hit
lastdb -uNEAR -R01 A2_db A2.transcripts.fa
lastdb -uNEAR -R01 D5_db D5.transcripts.fa
# -P8 use 8 cpu
lastal A2_db -P8 -D10000000000 D5.transcripts.fa | last-map-probs -m 0.49 > A.D.maf
lastal D5_db -P8 -D10000000000 A2.transcripts.fa | last-map-probs -m 0.49 > D.A.maf

# Create VCFs based on genotype differences between homeologs
python scripts/homeolog_genotypes.py -o A.vs.D -f exon -g A2.gtf A.D.maf D.A.maf # coordinates based on A2
python scripts/homeolog_genotypes.py -o D.vs.A -f exon -g D5.primaryOnly.gtf D.A.maf A.D.maf # coordinates based on D5

# Subgenome unique transcripts
cat A2.gtf | perl -ne 'chomp; m/transcript_id "(.*?)";/; print "$1\n";' | sort | uniq > A2.all.list
cat D5.primaryOnly.gtf | perl -ne 'chomp; m/transcript_id "(.*?)";/; print "$1\n";' | sort | uniq > D5.all.list
python scripts/tablize.py -v0 A.vs.D.reciprocal_best A2.all.list > A2.only.list
python scripts/tablize.py -v0 D.vs.A.reciprocal_best D5.all.list > D5.only.list

## Origin specific alignment with STAR

mkdir -p A2
STAR --runMode genomeGenerate --genomeDir A2 --genomeFastaFiles A2Du_13.fasta --sjdbGTFfile A2.gtf --runThreadN 36

mkdir -p D5
STAR --runMode genomeGenerate --genomeDir D5 --genomeFastaFiles Dgenome2_13.fasta --sjdbGTFfile D5.primaryOnly.gtf --runThreadN 36

## run STAR alignment, EARGLE-rc classification, and feature count
sbatch runStarEagle-rc.slurm

