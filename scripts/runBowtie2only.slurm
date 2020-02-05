#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-0:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "bowtieOnly"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=runBowtie2only.out.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever # use sinfo to get queue names

#!usr/bin/bash
echo ''
echo "Starting Job on "
date

#### Bowtie2 mapping against transcript reference
module load bowtie2/2.3.4.1-py2-jl5zqym samtools/1.9-k6deoga salmon/0.9.1-py2-mpich-434m5zb

cd /work/LAS/jfw-lab/hugj2006/eflen2020/bowtie2only 
# ls A2D5.transcripts.fa
# cat A2D5.transcripts.fa |sed '/>.*D/,+1d' >A2.transcripts.fa
# cat A2D5.transcripts.fa |sed '/>.*A/,+1d' >D5.transcripts.fa
# bowtie2-build A2D5.transcripts.fa A2D5
# bowtie2-build A2.transcripts.fa A2
# bowtie2-build D5.transcripts.fa D5


# A2 reads mapped to A2.transcript.fa (derived from D5 transcripts and snp index v4.0) reference for expected count
# D5 reads mapped to D5.transcript.fa reference for expected count
# ADs reads mapped to catenated A2D5.transcript.fa reference to observed count


## seed

for j in $( ls ../seedFastq | grep 'A2' );  do
    echo $j
    time bowtie2 -p 36 -q --local -x A2 -U ../seedFastq/$j | samtools view -Shb -o bam/${j/.cut.fq.gz/}.bam -
    samtools view -b -q 20 bam/${j/.cut.fq.gz/}.bam | samtools sort -o bam/${j/.cut.fq.gz/}.q20.sort.bam -
    samtools index bam/${j/.cut.fq.gz/}.q20.sort.bam
    salmon quant -l A -a bam/${j/.cut.fq.gz/}.q20.sort.bam -t A2.transcripts.fa -o bam/${j/.cut.fq.gz/}
done

for j in $( ls ../seedFastq | grep 'D5' );  do
    echo $j
    time bowtie2 -p 36 -q --local -x D5 -U ../seedFastq/$j | samtools view -Shb -o bam/${j/.cut.fq.gz/}.bam -
    samtools view -b -q 20 bam/${j/.cut.fq.gz/}.bam | samtools sort -o bam/${j/.cut.fq.gz/}.q20.sort.bam -
    samtools index bam/${j/.cut.fq.gz/}.q20.sort.bam
    salmon quant -l A -a bam/${j/.cut.fq.gz/}.q20.sort.bam -t D5.transcripts.fa -o bam/${j/.cut.fq.gz/}
done

for j in $( ls ../seedFastq | grep 'ADs' );  do
    echo $j
    time bowtie2 -p 36 -q --local -x A2D5 -U ../seedFastq/$j | samtools view -Shb -o bam/${j/.cut.fq.gz/}.bam -
    samtools view -b -q 20 bam/${j/.cut.fq.gz/}.bam | samtools sort -o bam/${j/.cut.fq.gz/}.q20.sort.bam -
    samtools index bam/${j/.cut.fq.gz/}.q20.sort.bam
    salmon quant -l A -a bam/${j/.cut.fq.gz/}.q20.sort.bam -t A2D5.transcripts.fa -o bam/${j/.cut.fq.gz/}
done


## flower

for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' |grep '\-A2\-'); do
    echo $j
    echo ${j%%_1.fq.gz}_2.fq.gz    
    time bowtie2 -p 36 -q --local --no-mixed --no-discordant --no-unal --dovetail -x A2 -1 ../flowerFastq/$j -2 ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz | samtools view -Shb -o bam/${j%%_1.fq.gz}.bam -
    samtools view -b -q 20 bam/${j%%_1.fq.gz}.bam | samtools sort -o bam/${j%%_1.fq.gz}.q20.sort.bam -
    samtools index bam/${j%%_1.fq.gz}.q20.sort.bam
    salmon quant -l A -a bam/${j%%_1.fq.gz}.q20.sort.bam -t A2.transcripts.fa -o bam/${j%%_1.fq.gz}
done

for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' |grep '\-D5\-'); do
    echo $j
    echo ${j%%_1.fq.gz}_2.fq.gz    
    time bowtie2 -p 36 -q --local --no-mixed --no-discordant --no-unal --dovetail -x D5 -1 ../flowerFastq/$j -2 ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz | samtools view -Shb -o bam/${j%%_1.fq.gz}.bam -
    samtools view -b -q 20 bam/${j%%_1.fq.gz}.bam | samtools sort -o bam/${j%%_1.fq.gz}.q20.sort.bam -
    samtools index bam/${j%%_1.fq.gz}.q20.sort.bam
    salmon quant -l A -a bam/${j%%_1.fq.gz}.q20.sort.bam -t D5.transcripts.fa -o bam/${j%%_1.fq.gz}
done

for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' |grep '\-ADs\-'); do
    echo $j
    echo ${j%%_1.fq.gz}_2.fq.gz    
    time bowtie2 -p 36 -q --local --no-mixed --no-discordant --no-unal --dovetail -x A2D5 -1 ../flowerFastq/$j -2 ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz | samtools view -Shb -o bam/${j%%_1.fq.gz}.bam -
    samtools view -b -q 20 bam/${j%%_1.fq.gz}.bam | samtools sort -o bam/${j%%_1.fq.gz}.q20.sort.bam -
    samtools index bam/${j%%_1.fq.gz}.q20.sort.bam
    salmon quant -l A -a bam/${j%%_1.fq.gz}.q20.sort.bam -t A2D5.transcripts.fa -o bam/${j%%_1.fq.gz}
done


echo ''
echo "Ending  Job on "
date
