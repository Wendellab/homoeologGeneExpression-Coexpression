# cd to hylite folder, now in /lss/research/jfw-lab/Projects/Eflen/flowerTimeDataset/mapping_hylite/results082417/bams
# both flower and seed data here

## quickly check bam result in comparison to hylite output, same order?
module load samtools
samtools flagstat AD-10-R1.sorted.bam
# 39610274 + 0 mapped (79.46% : N/A)
cd..
wc -l results082417.ADs.seed10-ADs-R1.read.txt 
# 39610275 results082417.ADs.seed10-ADs-R1.read.txt 

## get all read ID from mapped AD reads, check in A2 and D5 read list to see where they come from; then use hylite output to see which genome they are assigned to
# get mapped only -> read id 
cd bams/
samtools view -F 4 AD-10-R1.sorted.bam | awk '{print $1}' >AD.id.txt
# reads truly from A2
samtools view A2-10-R1.sorted.bam | awk '{print $1}'|sort|uniq >A2.txt
# reads truly from D5
samtools view D5-10-R1.sorted.bam | awk '{print $1}'|sort|uniq >D5.txt

# label subgenome origin for AD ids
perl labelAD.pl AD.id.txt A2.txt D5.txt >origin.txt

paste <(head results082417.ADs.seed10-ADs-R1.read.txt)| awk '{print $1, $6, $7}') <(head origin.txt) 

################# List all input files ###############

cat <(ls bams|grep '\-ADs\-.*bam$') <(ls bams|grep 'AD\-.*sorted.bam$') > fileLists.AD.txt
cat <(ls bams|grep '\-A2\-.*bam$') <(ls bams|grep '^A2\-.*sorted.bam$') > fileLists.A2.txt
cat <(ls bams|grep '\-D5\-.*bam$') <(ls bams|grep '^D5\-.*sorted.bam$') > fileLists.D5.txt
paste fileLists.AD.txt fileLists.A2.txt fileLists.D5.txt <(ls|grep '\-ADs\-.*.read.txt$') > fileLists.txt



################## {getOrigin.sh} ############################
module load samtools

for j in $( ls| grep 'AD-.*sorted.bam$' );	do
echo $j
echo A2-${j#AD-}
echo D5-${j#AD-}
echo ${j%%.*}.origin.txt

# AD mapped read id in order 
samtools view -F 4 $j | awk '{print $1}' >AD.id.txt
# reads truly from A2
samtools view A2-${j#AD-} | awk '{print $1}' >A2.txt
# reads truly from D5
samtools view D5-${j#AD-} | awk '{print $1}' >D5.txt
# perl
perl labelAD.pl A2.txt D5.txt AD.id.txt >origin.txt

# generate my file
echo results082417.ADs.seed10-ADs-R1.read.txt
echo results031417.AD.${j%%.*}.read.txt
paste <(awk '{print $1, $6, $7}' results031417.AD.${j%%.*}.read.txt) origin.txt >origin.${j%%.*}.txt

done

rm origin.txt
rm AD.id.txt
rm A2.txt
rm D5.txt

################## {labelAD.pl} ############################
# print out lines unique to file1 when compared to file 2
# Usage: script.pl fileA.txt fileD.txt fileAD.txt > outfile.txt

use strict;
use warnings;
use autodie;

my $f1 = shift || "A.s.txt";
my $f2 = shift || "D.s.txt";
my $f3 = shift || "AD.s.txt";

print "origin\n"
my %results;
open my $file1, '<', $f1;
while (my $line = <$file1>) { $results{$line} = "A2" }
open my $file2, '<', $f2;
while (my $line = <$file2>) { $results{$line} = "D5" }
open my $file3, '<', $f3;
while (my $line = <$file3>) { print "$results{$line}\n" }

################## R code #####################
x<- read.table("fileLists.txt")
for(i in 1:nrow(x)){
    fileAD<-x[i,1]
    fileA<-x[i,2]
    fileD<-x[i,3]
    filehylite<-x[i,4]
    
    # get mapped only -> read id 
    cmd<-paste0("samtools view -F 4 bams/", fileAD," | awk '{print $1}' >AD.id.txt")
    message(cmd); system(cmd);
    # reads truly from A2
    cmd<-paste0("samtools view bams/", fileA, " | awk '{print $1}'|sort|uniq >A2.txt")
    message(cmd); system(cmd);
    # reads truly from D5
    cmd<-paste0("samtools view bams/", fileD, " | awk '{print $1}'|sort|uniq >D5.txt")
    message(cmd); system(cmd);
    # run perl to get true origin
    cmd<-"perl labelAD.pl A2.txt D5.txt AD.id.txt >origin.txt"
    message(cmd); system(cmd);
    # include hylite assignment results
    cmd<-paste0("awk '{print $1, $6, $7}' ",filehylite, ">hl.txt; paste hl.txt origin.txt >origin/", gsub(".*ADs[.]|.read.txt","",filehylite),".txt")
    message(cmd); system(cmd);
    message("---------")
}

######################################################################
###  Output files as "origin.AD-10-R1.txt"
###   <gene>	<assigned>	<whetherNewSNP>	<trueOrigin>
###   GENE	CAT	NEW	Origin
###   Next run R code below to get datasets   
###      for counts table with origin info
######################################################################


load("~/working/output/R-01-hyliteDatasets.RData")
# "A2.Total"  "D5.Total"  "ADs.Total" "ADs.At"    "ADs.Dt"    "ADs.N"

id<-rownames(ADs.Total)
setwd("./origin/")
files<-list.files()

# Dt is ADs.Dt, and At is ADs.At
At <- data.frame(gene=id)
Dt <- data.frame(gene=id)
At.T <- data.frame(gene=id)
At.F <- data.frame(gene=id)
Dt.T <- data.frame(gene=id)
Dt.F <- data.frame(gene=id)

for(f in files) {
# f="LD7-ADs-L1.txt"
print(f)
x<-read.table(f, header=TRUE,sep=" ")
x$Origin <-gsub(".*\t","",x$NEW.Origin)
x$GENE <- gsub("[.]1$","",x$GENE)
dim(x<-x[x$CAT=="D5"|x$CAT=="A2",])
print("Correct assignment as:")
print(table(x$CAT == x$Origin)/nrow(x))
# Dt is ADs.Dt, and At is ADs.At
flag=gsub(".txt","",f)
At[,flag] <- as.numeric(table(x$GENE[x$CAT=="A2"])[id] )
Dt[,flag] <- as.numeric(table(x$GENE[x$CAT=="D5"])[id] )
print(head(At))
print(head(Dt))
At.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$Origin=="A2"])[id] )
At.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$Origin=="D5"])[id] )
Dt.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$Origin=="D5"])[id] )
Dt.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$Origin=="A2"])[id] )
}

for(i in c("At", "Dt","At.T","At.F","Dt.T","Dt.F"))
{
   count<-get(i)
   count[is.na(count)]=0
   count<-count[,-1]
   rownames(count)<-id
   assign(i,count)
}

# check if my tables replicate hylite outputs
unique(ADs.At-At)
unique(ADs.Dt-Dt)


save(A2.Total,D5.Total,ADs.Total,ADs.At,ADs.Dt, At.T, Dt.T, At.F, Dt.F, file="~/working/output/R-01-hyliteDatasets.true.RData")
