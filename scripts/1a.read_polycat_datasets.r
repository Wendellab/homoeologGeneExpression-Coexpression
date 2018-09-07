############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# PolyCat result files for seeds and flowers
seed_dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/seed_for_eflen_paper/polycatTests/htcountfiles/"
flower_dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/flowerTimeDataset/mapping_Gsnap/"

### seeds
dir<-seed_dir
fileL<-grep("txt",list.files(dir),value=TRUE)
 # 3 species x 11 samples x 5(total, A, D, N, AD)  = 165
remove(allcount)
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t")
    names(x) <- c("gene", file)
    if(!exists("allcount")) {allcount = x}
    else {allcount <- merge(allcount, x, by="gene")}
}
row.names(allcount) <- allcount$gene
names(allcount)<-gsub("sort.bam", "sort.T.bam", names(allcount) )
names(allcount)<-gsub("bam.", "", names(allcount) )
names(allcount)[-1]<-unlist(lapply(strsplit(names(allcount[-1]),"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )

### flower
dir<-flower_dir
fileL<-grep("-.*txt",list.files(dir),value=TRUE)
# 3 species x 22 samples x 5(total, A, D, N, AD)  = 330
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t")
    names(x) <- c("gene", file)
    if(!exists("allcountF")) {allcountF = x}
    else {allcountF <- merge(allcountF, x, by="gene")}
}
# merge put flower first
allcount<-merge(allcountF,allcount, by="gene")

# Prepare read count tables including only Chr genes
allcount13 <- allcount[grep("Gorai.0",allcount$gene),]
dim(allcount13)
names(allcount13)<-gsub(".pruned.namesort.bam|sort.|.txt", "", names(allcount13) )
names(allcount13)<-gsub("-AD-", "-ADs-", names(allcount13) )
# write polycat count tables
write.table(allcount13, "table.count.polycat.txt",sep="\t", row.names=FALSE)

allcount13<-read.table("table.count.polycat.txt",header=TRUE,sep="\t")
rownames(allcount13)<-allcount13$gene
##################
A2.Total  <- allcount13[,grep('A2.*[.]T$',names(allcount13))]
ADs.Total <- allcount13[,grep('AD.*[.]T$',names(allcount13))]
D5.Total  <- allcount13[,grep('[.]D5.*[.]T$',names(allcount13))]
##################
A2.At  <- allcount13[,grep('A2.*[.]A$',names(allcount13))]
ADs.At <- allcount13[,grep('AD.*[.]A$',names(allcount13))]
D5.At  <- allcount13[,grep('[.]D5.*[.]A$',names(allcount13))]
##################
A2.Dt  <- allcount13[,grep('A2.*[.]D$',names(allcount13))]
ADs.Dt <- allcount13[,grep('AD.*[.]D$',names(allcount13))]
D5.Dt  <- allcount13[,grep('[.]D5.*[.]D$',names(allcount13))]
##################
A2.N  <- allcount13[,grep('A2.*[.]N$',names(allcount13))]
ADs.N <- allcount13[,grep('AD.*[.]N$',names(allcount13))]
D5.N  <- allcount13[,grep('[.]D5.*[.]N$',names(allcount13))]
##################
A2.AD  <- allcount13[,grep('A2.*[.]AD$',names(allcount13))]
ADs.AD <- allcount13[,grep('AD.*[.]AD$',names(allcount13))]
D5.AD  <- allcount13[,grep('[.]D5.*[.]AD$',names(allcount13))]

##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
colSums(ADs.Dt)/colSums(ADs.Total)

##### inspect grouping by tissue types
total<-cbind(cbind(A2.Total, D5.Total),ADs.Total)
info<- data.frame(sample=names(total), lib_size=colSums(total) )
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"[.]"))))
names(info)[3:5] <-c("tissue","genome","rep")
# check samples complete
xtabs(~genome+tissue,data=info)
#      tissue
# genome LD7 LD9 LDM SD5 SD7 SDM SDP9 SDPF seed10 seed20 seed30 seed40
#    A2    3   3   3   2   3   2    3    3      3      2      3      3
#    ADs   3   3   3   2   3   2    3    3      3      2      3      3
#    D5    3   3   3   2   3   2    3    3      3      2      3      3

## make N partition table
ADs.Aportion <-ADs.At/( ADs.At+ ADs.Dt )
colMeans(ADs.Aportion)   # NA means At=0 and (At+Dt)=0
ADs.Aportion[is.na(ADs.Aportion)]<-0
ADs.AtN<- ADs.Total * ADs.Aportion
##################
ADs.Dportion <-ADs.Dt/( ADs.At+ ADs.Dt)
colMeans(ADs.Dportion)   # NA means At=0 and (At+Dt)=0
ADs.Dportion[is.na(ADs.Dportion)]<-0
ADs.DtN<- ADs.Total * ADs.Dportion
##################

# double check files
(colSums(ADs.At) + colSums(ADs.Dt) +colSums(ADs.N) +colSums(ADs.AD))/colSums(ADs.Total)
(colSums(ADs.AtN) + colSums(ADs.DtN ))/colSums(ADs.Total)  # close to 1

##### Plot grouping
norm<-sweep(total,2,info$lib_size,"/")*10^6
pca=prcomp(t(log2(norm+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.polycat.log2cpm.pdf")
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
save(A2.Total, D5.Total, ADs.Total, A2.At, D5.At, ADs.At, A2.Dt, D5.Dt, ADs.Dt, ADs.AtN, ADs.DtN, A2.N, D5.N, ADs.N, A2.AD, D5.AD, ADs.AD, file = "R-01-polycatDatasets.RData")


### prepare network datasets of A2D5, A2.D5.tech and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- paste(info[names(ADs.Total),3], gsub("^.","R",gsub(".T","",info[names(ADs.Total),5])),sep=".")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    rownames(a) <- paste0(rownames(a) , "a")
    rownames(d) <- paste0(rownames(d) , "d")
    ad        <- rbind (a,d)
    return(ad)
}
## make raw network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
ADs.ncorrect <- prepDatasets(ADs.AtN,  ADs.DtN )

# make list
networks  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs,
ADs.ncorrect  = ADs.ncorrect)

# double check content: 74446    33
dim(A2D5)
dim(A2D5.tech)
dim(ADs)
dim(ADs.ncorrect)



### rlog transformation
# need DESeq2
library(DESeq2)
rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    coldata<-data.frame( sample = names(count), tissue = gsub("[.]R.","", names(count)), rep = gsub(".*[.]","", names(count)) )
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    names(count.rld)<-names(count)
    return(count.rld)
}
networks.rld <- list(A2D5      = rlogTransformation(A2D5),
                     A2D5.tech = rlogTransformation(A2D5.tech),
                     ADs       = rlogTransformation(ADs),
                  ADs.ncorrect = rlogTransformation(ADs.ncorrect)   )

# RPKM transformation
# RPKM stands for Reads Per Kilobase of transcript per Million mapped reads.
# RPKM= 10^9 * number of mapped reads in gene
#             ----------------------------------------------------
#             total reads *  exon length
geneLen<-read.table("~/jfw-lab/Projects/Eflen/eflen_recheck/D5.trueANDtheo.lengths",header=TRUE, sep="\t")
table(geneLen$true>=geneLen$theoretical)
quantile(geneLen$true-geneLen$theoretical)
# prepare RPKM tranformation
rpkm<-function(count, length) {
    # match gene order
    gene <-  gsub(".$","", rownames(count)[1:37223])
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
# "A2D5"            "A2D5.tech"       "ADs"             "ADs.ncorrect"
for(i in 1:4) {print(table(is.na(as.matrix(networks.rpkm[[i]])) ) )}
for(i in 1:4) {print(table(is.infinite(as.matrix(networks.rpkm[[i]])) ) )}

#save
coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, geneLen, networks.rpkm, file = "R-01-polycatNetworkDatasets.RData")

###### Output files are:
# "pca.polycat.log2cpm.pdf"
# "table.count.polycat.txt"
# "R-01-polycatDatasets.RData"
# "R-01-polycatNetworkDatasets.RData"
