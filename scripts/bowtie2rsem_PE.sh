#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

module load bowtie2
module load rsem

for j in $( ls fastq/ | grep '_1.fq.gz' );	do
   
echo ''
echo "==========Running bowtie2-RSEM for"
echo $j
echo ""

## rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
rsem-calculate-expression -p 8 --bowtie2 --time --paired-end fastq/$j fastq/${j%%_1.fq.gz}_2.fq.gz ~/jfw-lab/GenomicResources/pseudogenomes/A2D5 mapping_rsem/${j%%_1.fq.gz} >mapping_rsem/${j%%_1.fq.gz}.log 2>&1

done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd
