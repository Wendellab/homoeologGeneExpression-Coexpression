#!/usr/bin/env Rscript
DIR <- "/work/maa146/EffectiveLength/kallisto/"
targets <- local({
     con <- textConnection(
       "\"LD7-A2-L1\"	\"flower\"	\"LD7\"	\"A2\"	\"L1\"
\"LD7-A2-L2\"	\"flower\"	\"LD7\"	\"A2\"	\"L2\"
\"LD7-A2-L3\"	\"flower\"	\"LD7\"	\"A2\"	\"L3\"
\"LD7-ADs-L1\"	\"flower\"	\"LD7\"	\"ADs\"	\"L1\"
\"LD7-ADs-L2\"	\"flower\"	\"LD7\"	\"ADs\"	\"L2\"
\"LD7-ADs-L3\"	\"flower\"	\"LD7\"	\"ADs\"	\"L3\"
\"LD7-D5-L1\"	\"flower\"	\"LD7\"	\"D5\"	\"L1\"
\"LD7-D5-L2\"	\"flower\"	\"LD7\"	\"D5\"	\"L2\"
\"LD7-D5-L3\"	\"flower\"	\"LD7\"	\"D5\"	\"L3\"
\"LD9-A2-L1\"	\"flower\"	\"LD9\"	\"A2\"	\"L1\"
\"LD9-A2-L2\"	\"flower\"	\"LD9\"	\"A2\"	\"L2\"
\"LD9-A2-L3\"	\"flower\"	\"LD9\"	\"A2\"	\"L3\"
\"LD9-ADs-L1\"	\"flower\"	\"LD9\"	\"ADs\"	\"L1\"
\"LD9-ADs-L2\"	\"flower\"	\"LD9\"	\"ADs\"	\"L2\"
\"LD9-ADs-L3\"	\"flower\"	\"LD9\"	\"ADs\"	\"L3\"
\"LD9-D5-L1\"	\"flower\"	\"LD9\"	\"D5\"	\"L1\"
\"LD9-D5-L2\"	\"flower\"	\"LD9\"	\"D5\"	\"L2\"
\"LD9-D5-L3\"	\"flower\"	\"LD9\"	\"D5\"	\"L3\"
\"LDM-A2-L1\"	\"flower\"	\"LDM\"	\"A2\"	\"L1\"
\"LDM-A2-L2\"	\"flower\"	\"LDM\"	\"A2\"	\"L2\"
\"LDM-A2-L3\"	\"flower\"	\"LDM\"	\"A2\"	\"L3\"
\"LDM-ADs-L1\"	\"flower\"	\"LDM\"	\"ADs\"	\"L1\"
\"LDM-ADs-L2\"	\"flower\"	\"LDM\"	\"ADs\"	\"L2\"
\"LDM-ADs-L3\"	\"flower\"	\"LDM\"	\"ADs\"	\"L3\"
\"LDM-D5-L1\"	\"flower\"	\"LDM\"	\"D5\"	\"L1\"
\"LDM-D5-L2\"	\"flower\"	\"LDM\"	\"D5\"	\"L2\"
\"LDM-D5-L3\"	\"flower\"	\"LDM\"	\"D5\"	\"L3\"
\"SD5-A2-S1\"	\"flower\"	\"SD5\"	\"A2\"	\"S1\"
\"SD5-A2-S5\"	\"flower\"	\"SD5\"	\"A2\"	\"S5\"
\"SD5-ADs-S1\"	\"flower\"	\"SD5\"	\"ADs\"	\"S1\"
\"SD5-ADs-S3\"	\"flower\"	\"SD5\"	\"ADs\"	\"S3\"
\"SD5-D5-S1\"	\"flower\"	\"SD5\"	\"D5\"	\"S1\"
\"SD5-D5-S3\"	\"flower\"	\"SD5\"	\"D5\"	\"S3\"
\"SD7-A2-S1\"	\"flower\"	\"SD7\"	\"A2\"	\"S1\"
\"SD7-A2-S4\"	\"flower\"	\"SD7\"	\"A2\"	\"S4\"
\"SD7-A2-S5\"	\"flower\"	\"SD7\"	\"A2\"	\"S5\"
\"SD7-ADs-S1\"	\"flower\"	\"SD7\"	\"ADs\"	\"S1\"
\"SD7-ADs-S2\"	\"flower\"	\"SD7\"	\"ADs\"	\"S2\"
\"SD7-ADs-S3\"	\"flower\"	\"SD7\"	\"ADs\"	\"S3\"
\"SD7-D5-S1\"	\"flower\"	\"SD7\"	\"D5\"	\"S1\"
\"SD7-D5-S2\"	\"flower\"	\"SD7\"	\"D5\"	\"S2\"
\"SD7-D5-S3\"	\"flower\"	\"SD7\"	\"D5\"	\"S3\"
\"SDM-A2-S5\"	\"flower\"	\"SDM\"	\"A2\"	\"S5\"
\"SDM-A2-S6\"	\"flower\"	\"SDM\"	\"A2\"	\"S6\"
\"SDM-ADs-S1\"	\"flower\"	\"SDM\"	\"ADs\"	\"S1\"
\"SDM-ADs-S2\"	\"flower\"	\"SDM\"	\"ADs\"	\"S2\"
\"SDM-D5-S1\"	\"flower\"	\"SDM\"	\"D5\"	\"S1\"
\"SDM-D5-S2\"	\"flower\"	\"SDM\"	\"D5\"	\"S2\"
\"SDP9-A2-S1\"	\"flower\"	\"SDP9\"	\"A2\"	\"S1\"
\"SDP9-A2-S3\"	\"flower\"	\"SDP9\"	\"A2\"	\"S3\"
\"SDP9-A2-S4\"	\"flower\"	\"SDP9\"	\"A2\"	\"S4\"
\"SDP9-ADs-S1\"	\"flower\"	\"SDP9\"	\"ADs\"	\"S1\"
\"SDP9-ADs-S2\"	\"flower\"	\"SDP9\"	\"ADs\"	\"S2\"
\"SDP9-ADs-S3\"	\"flower\"	\"SDP9\"	\"ADs\"	\"S3\"
\"SDP9-D5-S1\"	\"flower\"	\"SDP9\"	\"D5\"	\"S1\"
\"SDP9-D5-S2\"	\"flower\"	\"SDP9\"	\"D5\"	\"S2\"
\"SDP9-D5-S3\"	\"flower\"	\"SDP9\"	\"D5\"	\"S3\"
\"SDPF-A2-S1\"	\"flower\"	\"SDPF\"	\"A2\"	\"S1\"
\"SDPF-A2-S3\"	\"flower\"	\"SDPF\"	\"A2\"	\"S3\"
\"SDPF-A2-S6\"	\"flower\"	\"SDPF\"	\"A2\"	\"S6\"
\"SDPF-ADs-S1\"	\"flower\"	\"SDPF\"	\"ADs\"	\"S1\"
\"SDPF-ADs-S2\"	\"flower\"	\"SDPF\"	\"ADs\"	\"S2\"
\"SDPF-ADs-S3\"	\"flower\"	\"SDPF\"	\"ADs\"	\"S3\"
\"SDPF-D5-S1\"	\"flower\"	\"SDPF\"	\"D5\"	\"S1\"
\"SDPF-D5-S2\"	\"flower\"	\"SDPF\"	\"D5\"	\"S2\"
\"SDPF-D5-S3\"	\"flower\"	\"SDPF\"	\"D5\"	\"S3\"
\"A2-10-R1\"	\"seed\"	\"10\"	\"A2\"	\"R1\"
\"A2-10-R2\"	\"seed\"	\"10\"	\"A2\"	\"R2\"
\"A2-10-R3\"	\"seed\"	\"10\"	\"A2\"	\"R3\"
\"A2-20-R1\"	\"seed\"	\"20\"	\"A2\"	\"R1\"
\"A2-20-R2\"	\"seed\"	\"20\"	\"A2\"	\"R2\"
\"A2-20-R3\"	\"seed\"	\"20\"	\"A2\"	\"R3\"
\"A2-30-R1\"	\"seed\"	\"30\"	\"A2\"	\"R1\"
\"A2-30-R2\"	\"seed\"	\"30\"	\"A2\"	\"R2\"
\"A2-30-R3\"	\"seed\"	\"30\"	\"A2\"	\"R3\"
\"A2-40-R1\"	\"seed\"	\"40\"	\"A2\"	\"R1\"
\"A2-40-R2\"	\"seed\"	\"40\"	\"A2\"	\"R2\"
\"A2-40-R3\"	\"seed\"	\"40\"	\"A2\"	\"R3\"
\"ADs-10-R1\"	\"seed\"	\"10\"	\"ADs\"	\"R1\"
\"ADs-10-R2\"	\"seed\"	\"10\"	\"ADs\"	\"R2\"
\"ADs-10-R3\"	\"seed\"	\"10\"	\"ADs\"	\"R3\"
\"ADs-20-R1\"	\"seed\"	\"20\"	\"ADs\"	\"R1\"
\"ADs-20-R3\"	\"seed\"	\"20\"	\"ADs\"	\"R3\"
\"ADs-30-R1\"	\"seed\"	\"30\"	\"ADs\"	\"R1\"
\"ADs-30-R2\"	\"seed\"	\"30\"	\"ADs\"	\"R2\"
\"ADs-30-R3\"	\"seed\"	\"30\"	\"ADs\"	\"R3\"
\"ADs-40-R1\"	\"seed\"	\"40\"	\"ADs\"	\"R1\"
\"ADs-40-R2\"	\"seed\"	\"40\"	\"ADs\"	\"R2\"
\"ADs-40-R3\"	\"seed\"	\"40\"	\"ADs\"	\"R3\"
\"D5-10-R1\"	\"seed\"	\"10\"	\"D5\"	\"R1\"
\"D5-10-R2\"	\"seed\"	\"10\"	\"D5\"	\"R2\"
\"D5-10-R3\"	\"seed\"	\"10\"	\"D5\"	\"R3\"
\"D5-20-R1\"	\"seed\"	\"20\"	\"D5\"	\"R1\"
\"D5-20-R2\"	\"seed\"	\"20\"	\"D5\"	\"R2\"
\"D5-20-R3\"	\"seed\"	\"20\"	\"D5\"	\"R3\"
\"D5-30-R1\"	\"seed\"	\"30\"	\"D5\"	\"R1\"
\"D5-30-R2\"	\"seed\"	\"30\"	\"D5\"	\"R2\"
\"D5-30-R3\"	\"seed\"	\"30\"	\"D5\"	\"R3\"
\"D5-40-R1\"	\"seed\"	\"40\"	\"D5\"	\"R1\"
\"D5-40-R2\"	\"seed\"	\"40\"	\"D5\"	\"R2\"
\"D5-40-R3\"	\"seed\"	\"40\"	\"D5\"	\"R3\""
     )
     res <- utils::read.table(
       con,
       header    = FALSE,
       row.names = NULL,
       sep       = "\t",
       as.is     = TRUE
     )
     close(con)
     res
   })
setwd(DIR)
ROOT=system('git rev-parse --show-toplevel', intern=T)
.libPaths(c(file.path(ROOT, "Rlib"), .libPaths()))

library(edgeR)

colnames(targets) <- c("Library", "Set", "Tissue", "Genome", "Rep")
row.names(targets) <- targets$Library
targets$files = file.path(targets$Set, targets$Library, 'abundance.tsv')

d.count <- readDGE(targets, columns=c(1,4))
d.tpm   <- readDGE(targets, columns=c(1,5))

write.table(d.count$count, "table.count.kallisto.txt", sep="\t", col.names=NA)
write.table(d.tpm$count,   "table.tpm.kallisto.txt",   sep="\t", col.names=NA)

###########Jing modified table count column names, and edited below 9-5-18
# raw count
# all <- d.count$count
# raw count
all <- read.table("table.count.kallisto.txt", sep="\t",header=TRUE)
all <- all[,-grep("seed20.*R2",names(all))]
rownames(all)<-gsub(".1.D","d", gsub(".1.A","a", rownames(all)))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)]
ADs.Total <- ADs.Dt+ ADs.At

##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
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


##### Plot grouping
norm<-sweep(total,2,info$lib_size,"/")*10^6
pca=prcomp(t(log2(norm+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.kallisto.log2cpm.pdf")
print( ggplot(aes(PC1, PC2, color=info$tissue, shape=info$genome),data=dat) + geom_point() )
library(scales)
library(ape)
hc<-hclust( dist(t(log2(norm+1))) )
tre <- as.phylo(hc)
tre$tip.label == info$sample
# try to match ggplot color: library(scale); show_col(col4<-hue_pal()(4))
tipCol <- factor(info$tissue)
levels(tipCol) <- hue_pal()(nlevels(tipCol))
plot(tre,tip.col =as.character(tipCol), type="unrooted",cex=0.6, no.margin=TRUE)
plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
dev.off()

# save
save(A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt, file = "R-01-kallistoDatasets.RData")



### prepare network datasets of A2D5, A2.D5.tech and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- paste(info[info$genome=="ADs","tissue"], gsub("^.","R",info[info$genome=="ADs","rep"]),sep=".")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    ad        <- rbind (a,d)
    return(ad)
}
## make raw network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs   )


all <- read.table("table.tpm.kallisto.txt", sep="\t",header=TRUE)
all <- all[,-grep("seed20.*R2",names(all))]
rownames(all)<-gsub(".1.D","d", gsub(".1.A","a", rownames(all)))

A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",colnames(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("\\bD5\\b",colnames(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",colnames(all),value=TRUE)]
ADs.Total <- ADs.Dt+ ADs.At
##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
colSums(ADs.Dt)/colSums(ADs.Total)
##### inspect grouping by tissue types
total <- cbind(cbind(A2.Total, D5.Total),ADs.Total)
info  <- data.frame(sample=colnames(total), lib_size=colSums(total) )
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"[.]"))))
names(info)[3:5] <-c("tissue","genome","rep")
# check samples complete
xtabs(~genome+tissue,data=info)


A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs <- prepDatasets(ADs.At, ADs.Dt )

networks.tpm  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs   )


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
networks.rld <- list(A2D5 = rlogTransformation(networks$A2D5), A2D5.tech = rlogTransformation(networks$A2D5.tech), ADs = rlogTransformation(networks$ADs) )

coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, networks.tpm, file = "R-01-kallistoNetworkDatasets.RData")


####### Output files are:
# "pca.kallisto.log2cpm.pdf"
# "table.count.kallisto.txt"
# "table.tpm.kallisto.txt"
# "R-01-kallistoDatasets.RData"
# "R-01-kallistoNetworkDatasets.RData"
