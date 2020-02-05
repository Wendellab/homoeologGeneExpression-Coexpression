#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --constraint=AVX2 
#SBATCH -t 7-0:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 36   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "eagle SE"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=runEagle.eagleSE.out.txt
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
# seed
for i in `ls ../../seedFastq/`; do
    F=`basename $i .cut.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../A2Du_13.fasta -v ../A.vs.D.gtf.vcf --splice --rc 1> $F.A.vs.D.txt 2> $F.A.vs.D.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc -a $F.refsort.bam --listonly -o $F.A.vs.D -v $F.A.vs.D.txt $F.A.vs.D.readinfo.txt > $F.A.vs.D.list
done
cd ..

## likelhood using D5 reference
cd starD5
# seed
for i in `ls ../../seedFastq/`; do
    F=`basename $i .cut.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../Dgenome2_13.fasta -v ../D.vs.A.gtf.vcf --splice --rc 1> $F.D.vs.A.txt 2> $F.D.vs.A.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc -a $F.refsort.bam --listonly -o $F.D.vs.A -v $F.D.vs.A.txt $F.D.vs.A.readinfo.txt > $F.D.vs.A.list
done
cd ..

## use likelihood information to partition BAM and count reads
# mkdir -p eagleBam
for i in `ls ../seedFastq/`; do
    F=`basename $i .cut.fq.gz`
    python scripts/ref2_consensus.py --pe -u -o eagleBam/$F.ref \
        -A starA2/$F.A.vs.D.list \
        -B starD5/$F.D.vs.A.list
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starA2/$F.refsort.bam -o eagleBam/$F.A eagleBam/$F.ref.chrA.list >flush.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starD5/$F.refsort.bam -o eagleBam/$F.D eagleBam/$F.ref.chrB.list >flush.txt
    featureCounts -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.counts.txt eagleBam/$F.A.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.counts.txt eagleBam/$F.D.ref.bam
done


echo ''
echo "Ending  Job on "
date
