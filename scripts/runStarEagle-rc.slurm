#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

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

#### STAR Alignment
module load last/869-56gezob star/2.5.3a-rdzqclb samtools/1.9-k6deoga

###A2

## seed
for j in $( ls ../seedFastq );  do
    echo $j
    time STAR --genomeDir A2 --readFilesCommand zcat --readFilesIn ../seedFastq/$j \
        --outFileNamePrefix starA2/${j/.cut.fq.gz/}- --runThreadN 36 --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > starA2/${j/.cut.fq.gz/}.bam
    samtools sort -o starA2/${j/.cut.fq.gz/}.refsort.bam starA2/${j/.cut.fq.gz/}.bam
    samtools index -c starA2/${j/.cut.fq.gz/}.refsort.bam      
done

## flower
for j in $( ls ../flowerFastq/ | grep '_1.fq.gz' ); do
    echo $j
    echo ${j%%_1.fq.gz}_2.fq.gz    
    time STAR --genomeDir A2 --readFilesCommand zcat --readFilesIn ../flowerFastq/$j ../flowerFastq/${j%%_1.fq.gz}_2.fq.gz\
        --outFileNamePrefix starA2/${j%%_1.fq.gz}- --runThreadN 36 --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random \
        --outFilterType BySJout --outStd SAM | samtools view -Shb - > starA2/${j%%_1.fq.gz}.bam
    samtools sort -o starA2/${j%%_1.fq.gz}.refsort.bam starA2/${j%%_1.fq.gz}.bam
    samtools index -c starA2/${j%%_1.fq.gz}.refsort.bam   
done

###D5

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


## EAGLE-RC: read classification and the quantification with featureCounts
module load py-numpy/1.15.2-py3-wwyx7ek py-scipy subread/1.6.0-ak6vxhs

## likelhood using A2 reference
cd starA2
# seed
for i in `ls ../../seedFastq/`; do
    F=`basename $i .cut.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../A2Du_13.fasta -v ../A.vs.D.gtf.vcf --splice --rc 1> $F.A.vs.D.txt 2> $F.A.vs.D.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc -a $F.refsort.bam --listonly -o $F.A.vs.D -v $F.A.vs.D.txt $F.A.vs.D.readinfo.txt > $F.A.vs.D.list
done
# flower PE
for i in `ls ../../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../A2Du_13.fasta -v ../A.vs.D.gtf.vcf --splice --rc 1> $F.A.vs.D.txt 2> $F.A.vs.D.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired -a $F.refsort.bam --listonly -o $F.A.vs.D -v $F.A.vs.D.txt $F.A.vs.D.readinfo.txt > $F.A.vs.D.list
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
# flower PE
for i in `ls ../../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle -t 36 -a $F.refsort.bam -r ../Dgenome2_13.fasta -v ../D.vs.A.gtf.vcf --splice --rc 1> $F.D.vs.A.txt 2> $F.D.vs.A.readinfo.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired -a $F.refsort.bam --listonly -o $F.D.vs.A -v $F.D.vs.A.txt $F.D.vs.A.readinfo.txt > $F.D.vs.A.list
done
cd ..

## use likelihood information to partition BAM and count reads
mkdir -p eagleBam
# SE
for i in `ls ../seedFastq/`; do
    F=`basename $i .cut.fq.gz`
    echo $F
    python scripts/ref2_consensus.py --pe -u -o eagleBam/$F.ref \
        -A starA2/$F.A.vs.D.list \
        -B starD5/$F.D.vs.A.list
    echo 'run eagle-rc'
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starA2/$F.refsort.bam -o eagleBam/$F.A eagleBam/$F.ref.chrA.list >flush.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --refonly --readlist -a starD5/$F.refsort.bam -o eagleBam/$F.D eagleBam/$F.ref.chrB.list >flush.txt
    featureCounts -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.ref.counts.txt eagleBam/$F.A.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.alt.counts.txt eagleBam/$F.A.alt.bam
    featureCounts -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.unk.counts.txt eagleBam/$F.A.unk.bam
    featureCounts -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.ref.counts.txt eagleBam/$F.D.ref.bam
    featureCounts -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.alt.counts.txt eagleBam/$F.D.alt.bam
    featureCounts -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.unk.counts.txt eagleBam/$F.D.unk.bam
done
# PE
for i in `ls ../flowerFastq/ |  grep '_1.fq.gz' `; do
    F=`basename $i _1.fq.gz`
    echo $F
    python scripts/ref2_consensus.py --pe -u -o eagleBam/$F.ref \
        -A starA2/$F.A.vs.D.list \
        -B starD5/$F.D.vs.A.list
    echo 'run eagle-rc'
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired --refonly --readlist -a starA2/$F.refsort.bam -o eagleBam/$F.A eagleBam/$F.ref.chrA.list >flush.txt
    /work/LAS/jfw-lab/hugj2006/tools/eagle/eagle-rc --paired --refonly --readlist -a starD5/$F.refsort.bam -o eagleBam/$F.D eagleBam/$F.ref.chrB.list >flush.txt
    featureCounts -p -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.ref.counts.txt eagleBam/$F.A.ref.bam
    featureCounts -p -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.alt.counts.txt eagleBam/$F.A.alt.bam
    featureCounts -p -T 8 -t exon -g transcript_id -a A2.gtf -o eagleBam/$F.A.unk.counts.txt eagleBam/$F.A.unk.bam
    featureCounts -p -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.ref.counts.txt eagleBam/$F.D.ref.bam
    featureCounts -p -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.alt.counts.txt eagleBam/$F.D.alt.bam
    featureCounts -p -T 8 -t exon -g transcript_id -a D5.primaryOnly.gtf -o eagleBam/$F.D.unk.counts.txt eagleBam/$F.D.unk.bam
done

## generate read count table
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.A.ref.counts.txt > eagle.A.ref.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.D.ref.counts.txt > eagle.D.ref.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.A.alt.counts.txt > eagle.A.alt.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.D.alt.counts.txt > eagle.D.alt.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.A.unk.counts.txt > eagle.A.unk.tsv
python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagleBam/*.D.unk.counts.txt > eagle.D.unk.tsv

# Homeolog counts in terms of A2 gene id, transcript level
#python scripts/tablize.py -a A.vs.D.reciprocal_best eagle.A.tsv | cut -f 1,3- | sort -k1 > eagle.A.homeolog.Aref.tsv
#python scripts/tablize.py -a D.vs.A.reciprocal_best eagle.D.tsv | cut -f 2,3- | sort -k1 > eagle.D.homeolog.Aref.tsv

# Homeolog counts in terms of D5 gene id, transcript level
python scripts/tablize.py -a A.vs.D.reciprocal_best eagle.A.ref.tsv | cut -f 2,3- | sort -k1 > eagle.A.homeolog.ref.tsv
python scripts/tablize.py -a D.vs.A.reciprocal_best eagle.D.ref.tsv | cut -f 1,3- | sort -k1 > eagle.D.homeolog.ref.tsv
python scripts/tablize.py -a A.vs.D.reciprocal_best eagle.A.alt.tsv | cut -f 2,3- | sort -k1 > eagle.A.homeolog.alt.tsv
python scripts/tablize.py -a D.vs.A.reciprocal_best eagle.D.alt.tsv | cut -f 1,3- | sort -k1 > eagle.D.homeolog.alt.tsv
python scripts/tablize.py -a A.vs.D.reciprocal_best eagle.A.unk.tsv | cut -f 2,3- | sort -k1 > eagle.A.homeolog.unk.tsv
python scripts/tablize.py -a D.vs.A.reciprocal_best eagle.D.unk.tsv | cut -f 1,3- | sort -k1 > eagle.D.homeolog.unk.tsv


echo ''
echo "Ending  Job on "
date

###### Below code relevant but not needed, since primary transcripts were used here and we don't examine subgenome unique reads

# Homeolog counts in terms of halleri gene id, gene level by summation of transcript counts
# perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.A.homeolog.tsv > temp.txt
# python scripts/tablize.py -add temp.txt > eagle.A.homeolog.genelevel.tsv
# perl -ne 'chomp; @t=split(/\s+/); @i=split(/\./, $t[0]); $a=join("\t", @t[1..$#t]); print "$i[0]\t$a\n";' eagle.D.homeolog.tsv > temp.txt
# python scripts/tablize.py -add temp.txt > eagle.D.homeolog.genelevel.tsv

# Subgenome unique mapped reads
# for i in `ls *_R1.fastq.gz`; do
#     F=`basename $i _R1.fastq.gz`
#     echo "" > dummy.txt
#     eagle-rc --refonly --readlist -a hal/$F.refsort.bam -u lyr/$F.refsort.bam -o eagle/$F.H.only dummy.txt
#     eagle-rc --refonly --readlist -a lyr/$F.refsort.bam -u hal/$F.refsort.bam -o eagle/$F.L.only dummy.txt
#     featureCounts -T 8 -t exon -g transcript_id -a $HALGTF.gtf -o eagle/$F.H.only.counts.txt eagle/$F.H.only.ref.bam
#     featureCounts -T 8 -t exon -g transcript_id -a $LYRGTF.gtf -o eagle/$F.L.only.counts.txt eagle/$F.L.only.ref.bam
# done

# python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.H.only.counts.txt > eagle.H.only.tsv
# python scripts/tablize.py -skip 1 -a -i 0 -c 6 eagle/*.L.only.counts.txt > eagle.L.only.tsv

# Subgenome unique mapped reads in subgenome unique genes
# python scripts/tablize.py -a H.only.list eagle.H.only.tsv > eagle.H.only.unique.tsv
# python scripts/tablize.py -a L.only.list eagle.L.only.tsv > eagle.L.only.unique.tsv
