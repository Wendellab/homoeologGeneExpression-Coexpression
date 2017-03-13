#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

module load gsnap/20160816
module load samtools/1.2
module load python
module load bambam/1.3

# For Single end reads

for j in $( ls | grep 'fq.gz$' );	do
   
echo ''
echo "===== Running GSNAP SE for"
echo $j
gsnap --gunzip -n 1 -N 1 -Q -t 6 --merge-distant-samechr --use-sarray=0 -D ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/snpindex4.1 -d D5 -v D5snp4.1 -A sam $j > ${j/.fq.gz/}.sam
# option -n tells to report one best alignment only, -N looks for novel splicing, -t 2 tells to use 2 threads, -Q output protein seq.
# bigram: ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/snpindex4.1/
# biocrunch: /home/jfw-lab-local/gmapdb/D5/snpindex4.1

echo ''
echo "===== Running SAMtools for" 
echo ${j/.fq.gz/}.sam
samtools view -Sb ${j/.fq.gz/}.sam > ${j/.fq.gz/}.bam
samtools sort -n ${j/.fq.gz/}.bam ${j/.fq.gz/}.nsort
samtools index ${j/.fq.gz/}.nsort.bam
rm ${j/.fq.gz/}.sam
rm ${j/.fq.gz/}.bam

echo ''
echo "===== Running Polycat for"
echo ${j/.fq.gz/}.nsort.bam
polyCat -x 1 -p 0 -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/snpindex4.1/D13.snp4.1 ${j/.fq.gz/}.nsort.bam
echo ""

echo ''
echo "===== Running HTSEq-count for"
echo ${j/.fq.gz/}.nsort.A.bam
echo ${j/.fq.gz/}.nsort.A.bam
echo ${j/.fq.gz/}.nsort.N.bam
# without name sort, Maximum alignment buffer size exceeded error
htseq-count -f bam --stranded=no -r name ${j/.fq.gz/}.nsort.A.bam ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > ${j/.fq.gz/}.nsort.A.txt
htseq-count -f bam --stranded=no -r name ${j/.fq.gz/}.nsort.D.bam ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > ${j/.fq.gz/}.nsort.D.txt
htseq-count -f bam --stranded=no -r name ${j/.fq.gz/}.nsort.N.bam ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > ${j/.fq.gz/}.nsort.N.txt
echo ""


done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd






# Among 24,589,150 SAM alignment pairs processed.
# __no_feature    1,948,227(7.92%) ===== reads not assigned to any feature, aka those mapped to non-genic regions here
# __ambiguous       416,480(1.69%) ===== assigned to more than one feature, probably reads assigned to region with multiple (anti-sense?) annotation 
# __too_low_aQual 1,930,408(7.85%) ===== skip all reads with alignment quality lower than the given minimum value, default: 10
# __not_aligned   2,276,510(9.26%) ===== in the SAM file without alignment
# __alignment_not_unique  35,698 (0.15%) ==== reads (or read pairs) with more than one reported alignment, recognized from the NH optional SAM field tag. maybe the pair mapped to different chromosome??
# The rest of counts assigned to genes sum up to 17,981,827 accounting for 73.13% of trimmed pairs.
