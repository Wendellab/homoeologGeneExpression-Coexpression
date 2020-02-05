## --------------------------------------------------------
## This code is for detecting effective transcript regions that are diagnostic of homoeolog origins. Modified based on `/jfw-lab/Projects/Eflen/eflen_recheck/methods.txt` from Corrinne Grover.
## --------------------------------------------------------
## Guanjing Hu. 01/21/2020

#############
### INPUT ###
#############

# 1. Transcript annotation GFF file (coding regions only): remove comment lines, and keep transcript names in column 9
exonFile <- "D5.exon_unnamed.gff"

# 2. Variant information, such as Homoeolog SNP index file from polyCat or vcf from Eagle-rc.
snpFile <- "D13.snp4.0.txt"
vcfFile <- "../eagle/D.vs.A.gtf.vcf"
library(data.table)
## Check SNP index or vcf
head(fread(snpFile))
#     Chr Pos A D
#   Chr01  59 C T
#   Chr01  62 T A
#   Chr01  72 G A
#   Chr01 124 G A
head(fread(vcfFile))
#      V1    V2 V3   V4 V5
#  Chr01 14860  .    T  C
#  Chr01 15022  . CAAC  -
#  Chr01 15083  .    T  C
#  Chr01 15297  .    T  C
#  Chr01 15314  .    T  C
#  Chr01 23458  .    C  A
#                                                                V6
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
#  Gorai.001G000400.1;evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388

# 3. RNA-seq sequencing read length, e.g. SE100 or PE50
readLen <- 100


###########
### FUN ###
###########
# require Bedtools installed and can call intersectBed
getEflen = function(gff, variants, readLen, out="eflen.txt"){
    ## Load libraries
    require(GenomicRanges)
    require(rtracklayer)
    require(genomation)
    require(data.table)
    
    ## read exon to get actual transcript length
    gr = gffToGRanges(gff )
    len<-data.frame(id=gr$group, len=width(gr))
    len<-aggregate(len$len,by=list(len$id),sum)
    names(len)<-c("id","trueLen")
    
    ## read variants loci
    loci = fread(variants,select=1:2)
    names(loci)=c("Chr","Pos")
    
    ## make the range around (and including) the SNP that a RNA-seq read could represent
    loci$low <- loci$Pos - (readLen-1)
    loci$high <- loci$Pos + (readLen-1)
    # need to remove minus positions
    loci$low[loci$low<=0]=1
    
    ## convert to GenomicRanges object, and merge overlapping features
    lociGrange <- makeGRangesFromDataFrame(loci, seqnames.field=c("Chr"), start.field=c("low"), end.field=c("high"))
    lociGrange <- reduce(lociGrange)
    export.gff(lociGrange, "loci.gff")

    ## run Bedtools intersect to get effective regions on exons
    system(paste0("intersectBed -a ",exonFile, " -b loci.gff > efregion.gff"))
    
    ## read in effective region results
    gr = gffToGRanges("efregion.gff" )
    efflen<-data.frame(id=gr$group, len=width(gr))
    efflen<-aggregate(efflen$len,by=list(efflen$id),sum)
    names(efflen)<-c("id","effectiveLen")

    ## combine and get %
    LEN<- merge(len,efflen,all.x=TRUE, all.y=TRUE)
    LEN$effectiveLen[is.na(LEN$effectiveLen)]=0
    LEN$percentEffective <- LEN$effectiveLen/LEN$trueLen
    write.table(LEN,file=out,sep="\t", row.names=FALSE)
    
    system("rm loci.gff efregion.gff")
    message(paste0("Effective lengths writen to file: ", out))
}

###########
### RUN ###
###########

getEflen(gff=exonFile, variants=snpFile, readLen=50, out="eflen.snp.rl50.txt")
getEflen(gff=exonFile, variants=snpFile, readLen=100, out="eflen.snp.rl100.txt")
getEflen(gff=exonFile, variants=snpFile, readLen=200, out="eflen.snp.rl200.txt")
getEflen(gff=exonFile, variants=snpFile, readLen=300, out="eflen.snp.rl300.txt")

getEflen(gff=exonFile, variants=vcfFile, readLen=50,  out="eflen.vcf.rl50.txt")
getEflen(gff=exonFile, variants=vcfFile, readLen=100, out="eflen.vcf.rl100.txt")
getEflen(gff=exonFile, variants=vcfFile, readLen=200, out="eflen.vcf.rl200.txt")
getEflen(gff=exonFile, variants=vcfFile, readLen=300, out="eflen.vcf.rl300.txt")
