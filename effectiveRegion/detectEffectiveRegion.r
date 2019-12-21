## --------------------------------------------------------
## This code is for detecting effective transcript regions that are diagnostic of homoeolog origins. Modified based on `/jfw-lab/Projects/Eflen/eflen_recheck/methods.txt` from Corrinne Grover.
## --------------------------------------------------------
## Guanjing Hu. 08/15/2017

## Description: For each transcript, only regions contraining homoeolog diagnostic SNPs could

## Required INPUT include:
# 1. Transcript (Exon) annotation GFF file: remove comment lines, and keep transcript names in column 9
exonFile <- "D5.exon_unnamed.gff"
# 2. Homoeolog SNP index file
snpFile <- "D13.snp4.0.txt"
# 3. RNA-seq sequencing read length, e.g. SE100 or PE50
readLen <- 100

## OUTPUT files will include
snpRangesFile <- paste0("output/snp.rl",readLen,".gff")
EffectiveRegionsFile <- paste0("output/D5.theoretical.rl",readLen,".gff")
eflenListFile <-paste0("output/eflenList.rl",readLen,".txt")

## Load libraries
library(GenomicRanges)
library(rtracklayer)
library(genomation)

## read exon to get actual transcript length
gr = gffToGRanges(exonFile )
len<-data.frame(id=gr$group, len=width(gr))
len<-aggregate(len$len,by=list(len$id),sum)
names(len)<-c("id","trueLen")

## read in SNP file BUT make sure that the comment character is removed from line 1 first
snp <- read.table(snpFile, sep = "\t")
names(snp) = c("Chr","Pos","A","D")
head(snp)
#     Chr Pos A D
#   Chr01  59 C T
#   Chr01  62 T A
#   Chr01  72 G A
#   Chr01 124 G A

## make the range around (and including) the SNP that a RNA-seq read could represent
snp$low <- snp$Pos - (readLen-1)
snp$high <- snp$Pos + (readLen-1)
# need to remove minus positions
snp[snp$low<0,]
snp$low[snp$low<0]=1

## convert to GenomicRanges object, and merge overlapping features
snpGrange <- makeGRangesFromDataFrame(snp, seqnames.field=c("Chr"), start.field=c("low"), end.field=c("high"))
snpGrange <- reduce(snpGrange)
export.gff(snpGrange, snpRangesFile)

## run Bedtools intersect to get effective regions on exons
cmd<-paste0("intersectBed -a ",exonFile, " -b ",snpRangesFile, " > ",EffectiveRegionsFile)
cmd
system(cmd)

## read in effective region results
gr = gffToGRanges(EffectiveRegionsFile )
efflen<-data.frame(id=gr$group, len=width(gr))
efflen<-aggregate(efflen$len,by=list(efflen$id),sum)
names(efflen)<-c("id","effectiveLen")

## combine and get %
LEN<- merge(len,efflen,all.x=TRUE, all.y=TRUE)
LEN$effectiveLen[is.na(LEN$effectiveLen)]=0
LEN$percentEffective <- LEN$effectiveLen/LEN$trueLen
write.table(LEN,file=eflenListFile,sep="\t", row.names=FALSE)
