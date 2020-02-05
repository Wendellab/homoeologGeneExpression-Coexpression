#!/usr/bin/env Rscript
library(edgeR)
fl <- system("ls /work/LAS/jfw-lab/hugj2006/eflen2020/bowtie2only/bamQ10/*/quant.sf", intern=T)
fla <- grep("/A2|-A2-",fl, value=T)
fld <- grep("/D5|-D5-",fl, value=T)
flad <- grep("ADs",fl, value=T)

#count
A2 <- as.data.frame(readDGE(fla, columns=c(1,5))$count)
D5 <- as.data.frame(readDGE(fld, columns=c(1,5))$count)
ADs <- as.data.frame(readDGE(flad, columns=c(1,5))$count)
A2$gene_id = rownames(A2)
D5$gene_id = rownames(D5)
ADs$gene_id = rownames(ADs)
count<-merge(merge(A2,D5,all.x=T,all.y=T,by="gene_id"), ADs, all.x=T,all.y=T,by="gene_id")
names(count) = gsub(".*bamQ10/|/quant|.cut.fq.gz","",names(count))
names(count)[grep("0",names(count))]<-unlist(lapply(strsplit(names(count)[grep("0",names(count))],"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )
count[is.na(count)]=0
write.table(count, "table.count.bowtie.txt", sep="\t", col.names=TRUE, row.names=FALSE)

#tpm
A2 <- as.data.frame(readDGE(fla, columns=c(1,4))$count)
D5 <- as.data.frame(readDGE(fld, columns=c(1,4))$count)
ADs <- as.data.frame(readDGE(flad, columns=c(1,4))$count)
A2$gene_id = rownames(A2)
D5$gene_id = rownames(D5)
ADs$gene_id = rownames(ADs)
count<-merge(merge(A2,D5,all.x=T,all.y=T,by="gene_id"), ADs, all.x=T,all.y=T,by="gene_id")
names(count) = gsub(".*bamQ10/|/quant|.cut.fq.gz","",names(count))
names(count)[grep("0",names(count))]<-unlist(lapply(strsplit(names(count)[grep("0",names(count))],"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )
count[is.na(count)]=0
write.table(count, "table.tpm.bowtie.txt", sep="\t", col.names=TRUE, row.names=FALSE)

###########modified table count column names to ensure flower first then seed samples
# raw count
# all <- d.count$count
all <- read.table("table.count.bowtie.txt", sep="\t",header=TRUE)
rownames(all)<-gsub(".1.D","d", gsub(".1.A","a", all$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)[c(12:33,1:11)]]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE),      grep("A2",colnames(all),value=TRUE)[c(12:33,1:11)]]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)[c(12:33,1:11)]]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)[c(12:33,1:11)]]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)[c(12:33,1:11)]]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)[c(12:33,1:11)]]
ADs.Total <- ADs.Dt+ ADs.At
##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total) # 1 A2.At == A2.Total
colSums(D5.At)/colSums(D5.Total) # 0 D5.Dt == D5.Total
colSums(ADs.At)/colSums(ADs.Total) # ~0.5
colSums(ADs.Dt)/colSums(ADs.Total) # ~0.5
colSums(ADs.At)/colSums(A2.Total)
colSums(ADs.Dt)/colSums(D5.Total)
colSums(A2.Total)/colSums(D5.Total)

##### inspect grouping by tissue types
total <- cbind(cbind(A2.Total, D5.Total),ADs.Total)
info  <- data.frame(sample=colnames(total), lib_size=colSums(total) )
# info  <- merge(info, targets, by.x='sample', by.y='Library')
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"[.]"))))
names(info)[3:5] <-c("tissue","genome","rep")
# check samples complete
xtabs(~genome+tissue,data=info)
#      tissue
# genome LD7 LD9 LDM SD5 SD7 SDM SDP9 SDPF seed10 seed20 seed30 seed40
#    A2    3   3   3   2   3   2    3    3      3      2      3      3
#    ADs   3   3   3   2   3   2    3    3      3      2      3      3
#    D5    3   3   3   2   3   2    3    3      3      2      3      3

##### Plot grouping
norm<-sweep(total,2,info$lib_size,"/")*10^6
pca=prcomp(t(log2(norm+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.bowtie.log2cpm.pdf")
print( ggplot(aes(PC1, PC2, color=info$tissue, shape=info$genome),data=dat) + geom_point() )
library(scales)
library(ape)
hc<-hclust( dist(t(log2(norm+1))) )
tre <- as.phylo(hc)
tre$tip.label == info$sample
# try to match ggplot color: library(scale); show_col(col4<-hue_pal()(4))
tipCol <- info$tissue
levels(tipCol) <- hue_pal()(nlevels(tipCol))
plot(tre,tip.col =as.character(tipCol), type="unrooted",cex=0.6, no.margin=TRUE)
plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
dev.off()

# save
save(A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt, file = "R-01-bowtieDatasets.RData")

### prepare network datasets of A2D5, A2.D5.tech and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- paste(info[names(ADs.Total),3], gsub("^.","R",info[names(ADs.Total),5]),sep=".")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    ad        <- rbind (a,d)
    return(ad)
}
## make raw network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks  <- list(A2D5 = A2D5, ADs = ADs   )

## TPM
all <- read.table("table.tpm.bowtie.txt", sep="\t",header=TRUE)
rownames(all)<-gsub(".1.D","d", gsub(".1.A","a", all$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)[c(12:33,1:11)]]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)[c(12:33,1:11)]]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)[c(12:33,1:11)]]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)[c(12:33,1:11)]]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)[c(12:33,1:11)]]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)[c(12:33,1:11)]]
ADs.Total <- ADs.Dt+ ADs.At
##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total) # variable, rather weird
colSums(ADs.Dt)/colSums(ADs.Total)
##### inspect grouping by tissue types
total <- cbind(cbind(A2.Total, D5.Total),ADs.Total)
info  <- data.frame(sample=colnames(total), lib_size=colSums(total) )
# info  <- merge(info, targets, by.x='sample', by.y='Library')
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"[.]"))))
names(info)[3:5] <-c("tissue","genome","rep")
# check samples complete
xtabs(~genome+tissue,data=info)
#      tissue
# genome LD7 LD9 LDM SD5 SD7 SDM SDP9 SDPF seed10 seed20 seed30 seed40
#    A2    3   3   3   2   3   2    3    3      3      2      3      3
#    ADs   3   3   3   2   3   2    3    3      3      2      3      3
#    D5    3   3   3   2   3   2    3    3      3      2      3      3
pca=prcomp(t(log2(total+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.bowtie.log2tpm.pdf")
print( ggplot(aes(PC1, PC2, color=info$tissue, shape=info$genome),data=dat) + geom_point() )
library(scales)
library(ape)
hc<-hclust( dist(t(log2(norm+1))) )
tre <- as.phylo(hc)
tre$tip.label == info$sample
# try to match ggplot color: library(scale); show_col(col4<-hue_pal()(4))
tipCol <- info$tissue
levels(tipCol) <- hue_pal()(nlevels(tipCol))
plot(tre,tip.col =as.character(tipCol), type="unrooted",cex=0.6, no.margin=TRUE)
plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
dev.off()

A2D5 <- prepDatasets(A2.Total, D5.Total)
ADs  <- prepDatasets(ADs.At, ADs.Dt )
networks.tpm  <- list(A2D5 = A2D5, ADs = ADs   )

# rlog transformation
# need DESeq2
library(DESeq2)
rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    coldata<-data.frame( sample = colnames(count), tissue = gsub("[.]R.","", colnames(count)), rep = gsub(".*[.]","", colnames(count)) )
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    colnames(count.rld)<-colnames(count)
    return(count.rld)
}
networks.rld <- list(A2D5 = rlogTransformation(networks$A2D5), ADs = rlogTransformation(networks$ADs) )

# RPKM transformation
# RPKM stands for Reads Per Kilobase of transcript per Million mapped reads.
# RPKM= 10^9 * number of mapped reads in gene
#             ----------------------------------------------------
#             total reads *  exon length
# "/lss/LAS/jfw-lab/Projects/Eflen/eflen_recheck/D5.trueANDtheo.lengths"
geneLen<-read.table("D5.trueANDtheo.lengths",header=TRUE, sep="\t")
table(geneLen$true>=geneLen$theoretical)
quantile(geneLen$true-geneLen$theoretical)
# prepare RPKM tranformation
rpkm<-function(count, length) {
    # match gene order
    gene <-  gsub(".$","", rownames(count)[1:nrow(A2.Total)])
    # table(match(gene, lenth$gene)==1:37223)
    len <- length[ match(gene, length$gene),2 ]
    len <- c(len, len)
    librarySize <- colSums(count)
    
    # OR x <- t( t(10^9 * count[,-1]/ geneLen) / librarySize )
    x<-sweep(sweep(count*10^9, 2 ,librarySize, FUN="/"), 1, len, FUN ="/" )
    return(x)
}
networks.rpkm<-lapply(networks, function(x){rpkm(x, geneLen[,c("gene","true" )])})
names(networks.rpkm)
# "A2D5"      "ADs"
for(i in 1:2) {print(table(is.na(as.matrix(networks.rpkm[[i]])) ) )}
for(i in 1:2) {print(table(is.infinite(as.matrix(networks.rpkm[[i]])) ) )}


coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, networks.rpkm, networks.tpm, file = "R-01-bowtieNetworkDatasets.RData")

####### Output files are:
# "pca.bowtie.log2cpm.pdf"
# "table.count.bowtie.txt"
# "table.tpm.bowtie.txt"
# "R-01-bowtieDatasets.RData"
# "R-01-bowtieNetworkDatasets.RData"
