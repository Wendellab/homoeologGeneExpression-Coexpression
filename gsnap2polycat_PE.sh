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

for j in $( ls | grep '_1.fq.gz$' );	do
   
echo ''
echo "===== Running GSNAP PE for"
echo $j
gsnap --gunzip -n 1 -N 1 -Q -t 10 --merge-distant-samechr --use-sarray=0 -D /home/jfw-lab-local/gmapdb/D5/snpindex4.0 -d D5 -v D5snp4.0 -A sam $j ${j/_1.fq.gz/}_2.fq.gz > ${j/_1.fq.gz/}.sam 2>>${j/_1.fq.gz/}.log

# option -n tells to report one best alignment only, -N looks for novel splicing, -t 2 tells to use 2 threads, -Q output protein seq.
# bigram: ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/snpindex4.1/
# biocrunch: /home/jfw-lab-local/gmapdb/D5/snpindex4.1

echo ''
echo "===== Pruning SAM for"
echo ${j/_1.fq.gz/}.sam
echo ".....Begin samtools pruning, from ... to ..." >> ${j/_1.fq.gz/}.log
echo $(wc -l ${j/_1.fq.gz/}.sam) >> ${j/_1.fq.gz/}.log
cat ${j/_1.fq.gz/}.sam | grep 'XO:Z:CU\|XO:Z:HU\|XO:Z:CT\|@SQ' | awk '{if ($7=="=" || $3=="=" || $7=="*" || $3=="*" || $1=="@SQ") print $0}' > ${j/_1.fq.gz/}.pruned.sam
echo $(wc -l ${j/_1.fq.gz/}.pruned.sam) >> ${j/_1.fq.gz/}.log

echo ''
echo "===== Running SAMtools for" 
echo ${j/_1.fq.gz/}.pruned.sam
echo ".....Begin samtools transformation" >> ${j/_1.fq.gz/}.log
samtools view -Sb -@ 10 ${j/_1.fq.gz/}.pruned.sam -o ${j/_1.fq.gz/}.pruned.bam 2>> ${j/_1.fq.gz/}.log
samtools sort -@ 10 -n ${j/_1.fq.gz/}.pruned.bam ${j/_1.fq.gz/}.pruned.namesort 2>> ${j/_1.fq.gz/}.log
samtools index ${j/_1.fq.gz/}.pruned.namesort.bam
rm ${j/_1.fq.gz/}.sam
rm ${j/_1.fq.gz/}.pruned.sam
rm ${j/_1.fq.gz/}.pruned.bam


echo ''
echo "===== Running Polycat for"
echo ${j/_1.fq.gz/}.pruned.namesort.bam
echo "Begin polycat partition" >> ${j/_1.fq.gz/}.log
polyCat -x 1 -p 1 -s /home/jfw-lab-local/gmapdb/D5/snpindex4.0/D13.snp4.0 ${j/_1.fq.gz/}.pruned.namesort.bam 2>> ${j/_1.fq.gz/}.log
# -p 1 for pairend
echo ""

echo ''
echo "===== Running HTSEq-count for"
echo ${j/_1.fq.gz/}.pruned.namesort.bam
echo ${j/_1.fq.gz/}.pruned.namesort.A.bam
echo ${j/_1.fq.gz/}.pruned.namesort.D.bam
echo ${j/_1.fq.gz/}.pruned.namesort.N.bam
echo ".....Begin counting reads for T, A, D and N" >> ${j/_1.fq.gz/}.log
# without name sort, Maximum alignment buffer size exceeded error
htseq-count -f bam --stranded=no -r name ${j/_1.fq.gz/}.pruned.namesort.bam D5.primaryOnly.gtf > ${j/_1.fq.gz/}.T.txt
htseq-count -f bam --stranded=no -r name ${j/_1.fq.gz/}.pruned.namesort.A.bam D5.primaryOnly.gtf > ${j/_1.fq.gz/}.A.txt
htseq-count -f bam --stranded=no -r name ${j/_1.fq.gz/}.pruned.namesort.D.bam D5.primaryOnly.gtf > ${j/_1.fq.gz/}.D.txt
htseq-count -f bam --stranded=no -r name ${j/_1.fq.gz/}.pruned.namesort.N.bam D5.primaryOnly.gtf > ${j/_1.fq.gz/}.N.txt
echo ""


done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd

