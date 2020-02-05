hea#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --constraint=AVX2  # request run on a machine that has AVX2 instructions in the CPU, working like a charm for local compiled program
#SBATCH -t 7-0:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "eagle"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=runEagle.eagle.out.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever # use sinfo to get queue names

#!usr/bin/bash
echo ''
echo "Starting Job on "
date

#module load bowtie2/2.3.4.1-py2-jl5zqym
#module load samtools/1.9-k6deoga
#module load last/869-56gezob star/2.5.3a-rdzqclb samtools/1.9-k6deoga
module load py-numpy/1.15.2-py3-wwyx7ek py-scipy subread/1.6.0-ak6vxhs

## EAGLE-RC: read classification and the quantification with featureCounts

## likelhood using A2 reference
cd starA2
# flower PE
for i in `ls ../../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../A2Du_13.fasta -v ../A.vs.D.gtf.vcf --splice --rc 1> $F.A.vs.D.txt 2> $F.A.vs.D.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired -a $F.refsort.bam --listonly -o $F.A.vs.D -v $F.A.vs.D.txt $F.A.vs.D.readinfo.txt > $F.A.vs.D.list
done
cd ..

## likelhood using D5 reference
cd starD5
# flower PE
for i in `ls ../../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../Dgenome2_13.fasta -v ../D.vs.A.gtf.vcf --splice --rc 1> $F.D.vs.A.txt 2> $F.D.vs.A.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired -a $F.refsort.bam --listonly -o $F.D.vs.A -v $F.D.vs.A.txt $F.D.vs.A.readinfo.txt > $F.D.vs.A.list
done
cd ..

## use likelihood information to partition BAM and count reads
# mkdir -p eagleBam
for i in `ls ../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    echo $F
    python scripts/ref2_consensus.py --pe -u -o eagleBam/$F.ref \
        -A starA2/$F.A.vs.D.list \
        -B starD5/$F.D.vs.A.list
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starA2/$F.refsort.bam -o eagleBam/$F.A eagleBam/$F.ref.chrA.list >flush.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starD5/$F.refsort.bam -o eagleBam/$F.D eagleBam/$F.ref.chrB.list >flush.txt
    featureCounts -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.counts.txt eagleBam/$F.A.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.counts.txt eagleBam/$F.D.ref.bam
done

## generate read count table
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.A.counts.txt > eagle.A.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.D.counts.txt > eagle.D.tsv

# Homeolog counts in terms of halleri gene id, transcript level
python scripts/tablize.py -a A.vs.D.reciprocal_best eagle.A.tsv | cut -f 1,3- | sort -k1 > eagle.A.homeolog.tsv
python scripts/tablize.py -a D.vs.A.reciprocal_best eagle.D.tsv | cut -f 2,3- | sort -k1 > eagle.D.homeolog.tsv


echo ''
echo "Ending  Job on "
date
