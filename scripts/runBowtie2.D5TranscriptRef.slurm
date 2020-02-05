#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-0:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Bowtie2 hylite D5"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=runBowtie2.D5TranscriptRef.out.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever # use sinfo to get queue names

#!usr/bin/bash
echo ''
echo "Starting Job on "
date

module load bowtie2/2.3.4.1-py2-jl5zqym
module load samtools/1.9-k6deoga

# seed
for j in $( ls ../seedFastq );	do
echo $j
time bowtie2 -x D5.transcripts -U ../seedFastq/$j -S samD5/${j/.cut.fq.gz/}.sam -N 1 -p 36 2>samD5/${j/.cut.fq.gz/}.log
done

for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' ); do
echo $j
echo ${j%%_1.fq.gz}_2.fq.gz
echo ${j%%_1.fq.gz}.sam
echo ${j%%_1.fq.gz}.log
time bowtie2 -q -p 36 -t -N 1 --no-mixed --no-discordant --no-unal --dovetail -x D5.transcripts -1 ../flowerFastq/$j -2 ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz -S samD5/${j%%_1.fq.gz}.sam 2>samD5/${j%%_1.fq.gz}.log
done

# hylite 
## 2018-1-12 HyLiTE has now been updated and rewritten in Python 3
module load python/3.6.3-u4oaxsb samtools/1.9-k6deoga py-scipy/1.1.0-py3-pzig4lr
python3 HyLiTE-2.0.1-py3-none-any.whl -v -S -f samD5_protocol_file.txt -r D5.transcripts.fa -n resultsD5 >cmd.hylite.D5.log 2>&1

echo ''
echo "Ending  Job on "
date
