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

for j in $( ls | grep 'namesort.bam' );	do
echo ''
echo "===== Running Polycat for"
echo $j
echo "Begin polycat partition"
#polyCat -x 1 -p 1 -s /home/jfw-lab-local/gmapdb/D5/snpindex4.0/D13.snp4.0 $j 2> ${j/[.].*/}.log
# -p 1 for pairend
echo ""

echo ''
echo "===== Running HTSEq-count for"
echo $j
echo ${j/.bam/}.A.bam
echo ${j/.bam/}.D.bam
echo ${j/.bam/}.N.bam
echo ${j/.bam/}.AD.bam
echo ".....Begin counting reads for T, A, D and N"
# without name sort, Maximum alignment buffer size exceeded error
htseq-count -f bam --stranded=no -r name $j D5.primaryOnly.gtf > ${j/[.].*/}.T.txt
htseq-count -f bam --stranded=no -r name ${j/.bam/}.A.bam D5.primaryOnly.gtf > ${j/[.].*/}.A.txt
htseq-count -f bam --stranded=no -r name ${j/.bam/}.D.bam D5.primaryOnly.gtf > ${j/[.].*/}.D.txt
htseq-count -f bam --stranded=no -r name ${j/.bam/}.AD.bam D5.primaryOnly.gtf > ${j/[.].*/}.AD.txt
htseq-count -f bam --stranded=no -r name ${j/.bam/}.N.bam D5.primaryOnly.gtf > ${j/[.].*/}.N.txt
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
