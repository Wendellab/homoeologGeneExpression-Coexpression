#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-0:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "StarD5"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=runEagle.starD5.out.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever # use sinfo to get queue names

#!usr/bin/bash
echo ''
echo "Starting Job on "
date

#module load bowtie2/2.3.4.1-py2-jl5zqym
#module load samtools/1.9-k6deoga
module load last/869-56gezob star/2.5.3a-rdzqclb samtools/1.9-k6deoga

## seed
for j in $( ls ../seedFastq );  do
    echo $j
    time STAR --genomeDir D5 --readFilesCommand zcat --readFilesIn ../seedFastq/$j \
        --outFileNamePrefix starD5/${j/.cut.fq.gz/}- --runThreadN 36 --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > starD5/${j/.cut.fq.gz/}.bam
    samtools sort -o starD5/${j/.cut.fq.gz/}.refsort.bam starD5/${j/.cut.fq.gz/}.bam
    samtools index -c starD5/${j/.cut.fq.gz/}.refsort.bam      
done

## flower
for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' ); do
    echo $j
    echo ${j%%_1.fq.gz}_2.fq.gz    
    time STAR --genomeDir D5 --readFilesCommand zcat --readFilesIn ../flowerFastq/$j ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz\
        --outFileNamePrefix starD5/${j%%_1.fq.gz}- --runThreadN 36 --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > starD5/${j%%_1.fq.gz}.bam
    samtools sort -o starD5/${j%%_1.fq.gz}.refsort.bam starD5/${j%%_1.fq.gz}.bam
    samtools index -c starD5/${j%%_1.fq.gz}.refsort.bam      
done

echo ''
echo "Ending  Job on "
date
