############### Step 1a. Prepare EAGLE-RC network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# eaglerc were ran for both seed and flower datasets together
dir<- "/work/LAS/jfw-lab/hugj2006/eflen2020/eagle/"

# read count mapped to A2 reference, and classified as At
Aref<-read.table(paste0(dir,"eagle.A.homeolog.ref.tsv"))
Dref<-read.table(paste0(dir,"eagle.D.homeolog.ref.tsv"))
Aalt<-read.table(paste0(dir,"eagle.A.homeolog.alt.tsv"))
Dalt<-read.table(paste0(dir,"eagle.D.homeolog.alt.tsv"))
Aunk<-read.table(paste0(dir,"eagle.A.homeolog.unk.tsv"))
Dunk<-read.table(paste0(dir,"eagle.D.homeolog.unk.tsv"))
#function to rename columns and reorder
fixCol<-function(a){
    rownames(a)=a$V1
    a=a[,-1]
    names(a) = gsub("[.].*","",list.files(paste0(dir,"eagleBam"),pattern="A.ref.counts.txt$"))
    names(a)[grep("0",names(a))]<-unlist(lapply(strsplit(names(a)[grep("0",names(a))],"-" ),function(x)paste0("seed",x[2],"-",x[1],"-",x[3])))
    a=a[,c(34:99,1:33)]
    return(a)
}
Aref<-fixCol(Aref)
Dref<-fixCol(Dref)
Aalt<-fixCol(Aalt)
Dalt<-fixCol(Dalt)
Aunk<-fixCol(Aunk)
Dunk<-fixCol(Dunk)
# ADs
select=grep("ADs",colnames(Aref),value=TRUE)
colSums(Aref[, select]);
colSums(Dref[, select]);
colSums(Aalt[, select]); #0 likely taken cared by eagle consensus
colSums(Dalt[, select]); #0
colSums(Aunk[, select]);
colSums(Dunk[, select]);
ADs.At <- Aref[, select] + Dalt[,select]
ADs.Dt <- Dref[, select] + Aalt[,select]
ADs.N <- Aunk[,select] + Dunk[,select]
ADs.Total = ADs.At + ADs.Dt + ADs.N
# A2
select=grep("A2",colnames(Aref),value=TRUE)
colSums(Aref[, select]);
colSums(Dref[, select]); # small amount, reads of A2 origin assigned to D5 as error
colSums(Aalt[, select]); #0 likely taken cared by eagle consensus
colSums(Dalt[, select]); #0
colSums(Aunk[, select]);
colSums(Dunk[, select]);
A2.At <- Aref[,select] + Dalt[,select]
A2.Dt <- Aalt[,select] + Dref[,select]#0
A2.N  <- Aunk[,select] + Dunk[,select]
A2.Total <- A2.At+ A2.Dt + A2.N
# D5
select=grep("\\bD5\\b",colnames(Aref),value=TRUE)
colSums(Aref[, select]); # small amount, reads of A2 origin assigned to D5 as error
colSums(Dref[, select]);
colSums(Aalt[, select]); #0 likely taken cared by eagle consensus
colSums(Dalt[, select]); #0
colSums(Aunk[, select]);
colSums(Dunk[, select]);
D5.Dt <-Dref[,select] + Aalt[,select]
D5.At <-Dalt[,select] + Aref[,select]
D5.N <-Dunk[,select] + Aunk[,select]
D5.Total <- D5.At+ D5.Dt + D5.N

##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.Dt)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
colSums(ADs.Dt)/colSums(ADs.Total)
colSums(ADs.N)/colSums(ADs.Total) #0.4~ 14%
colSums(ADs.Total)/(colSums(A2.Total)+colSums(D5.Total))  # near identical

##### inspect grouping by tissue types
total <- cbind(cbind(A2.Total, D5.Total),ADs.Total)
info  <- data.frame(sample=colnames(total), lib_size=colSums(total) )
# info  <- merge(info, targets, by.x='sample', by.y='Library')
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"-"))))
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
pdf("pca.eaglerc.log2cpm.pdf")
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
save(A2.At, A2.Dt, A2.N, D5.At, D5.Dt, D5.N, A2.Total, D5.Total, ADs.Total, ADs.At, ADs.Dt, ADs.N, file = "R-01-eaglercDatasets.RData")

### prepare network datasets of A2D5 and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- paste(info[names(ADs.Total),3], gsub("^.","R",info[names(ADs.Total),5]),sep=".")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    rownames(a) <- gsub("[.]1$","a",rownames(a) )
    rownames(d) <- gsub("[.]1$","d",rownames(d) )
    ad        <- rbind (a,d)
    return(ad)
}
## make raw network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make raw network list
networks  <- list(A2D5      = A2D5,
ADs       = ADs   )

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

# rlog transformation
# need DESeq2
library(DESeq2)
rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    coldata<-data.frame( sample = names(count), dpa = gsub("dev|[.]R.","", names(count)), rep = gsub(".*[.]","", names(count)) )
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    names(count.rld)<-names(count)
    return(count.rld)
}
networks.rld <- list(A2D5 = rlogTransformation(A2D5), ADs = rlogTransformation(ADs) )

#save
coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, geneLen, networks.rpkm, file = "R-01-eaglercNetworkDatasets.RData")


###### Output files are:
# "pca.eaglerc.log2cpm.pdf"
# "table.count.eaglerc.txt"
# "R-01-eaglercDatasets.RData"
# "R-01-eaglercNetworkDatasets.RData"
