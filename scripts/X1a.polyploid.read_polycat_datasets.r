############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# PolyCat result files for seeds and flowers, only using Yuc or TX2094, because flower dataset has no TM1
seed_dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/seed_AD1_2016Dec/"
flower_dir <- "/work/LAS/jfw-lab/hugj2006/eflen/flowerCount/"
# once all done move to "/lss/research/jfw-lab/Projects/Eflen/flowerTimeDataset/mapping_Gsnap/"

### seeds
dir<-seed_dir
fileL<-grep("AD1_Yuc.*txt",list.files(dir),value=TRUE)
# 1 accession x 12 samples x 3(A, D, N, )  = 36
remove(allcount)
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t")
    names(x) <- c("gene", file)
    if(!exists("allcount")) {allcount = x}
    else {allcount <- merge(allcount, x, by="gene")}
}
row.names(allcount) <- allcount$gene
names(allcount)<-gsub(".nsort|.txt|AD1_", "", names(allcount) )
names(allcount)[-1]<-unlist(lapply(strsplit(names(allcount[-1]),"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )

### flower
dir<-flower_dir
fileL<-grep("-.*txt",list.files(dir),value=TRUE)
# 1 species x 24 samples x 2(A, D)  = 48
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t")
    names(x) <- c("gene", file)
    if(!exists("allcountF")) {allcount = x}
    else {allcountF <- merge(allcount, x, by="gene")}
}
# merge put flower first
allcount<-merge(allcountF,allcount, by="gene")


# Prepare read count tables including only Chr genes
allcount13 <- allcount[grep("Gorai.0",allcount$gene),]
dim(allcount13)
names(allcount13)<-gsub("-TX2094-", "-Yuc-", gsub(".pruned|.namesort|-NF|.txt", "", names(allcount13) ) )
names(allcount13)<-gsub("-215", "b", names(allcount13) )
# write polycat count tables
write.table(allcount13, "table.count.polycat.Yuc.txt",sep="\t", row.names=FALSE)

#allcount13<-read.table("table.count.polycat.Yuc.txt",header=TRUE,sep="\t")
rownames(allcount13)<-allcount13$gene
##################
Yuc.At  <- allcount13[,grep('[.]A$',names(allcount13))]
Yuc.Dt  <- allcount13[,grep('[.]D$',names(allcount13))]
gsub("..$","",names(Yuc.At)) == gsub("..$","",names(Yuc.Dt))
# if polyploid samples are truely polyploid
colSums(Yuc.Dt)/colSums(Yuc.At)

##### inspect grouping by tissue types
Yuc.AnD <- Yuc.At + Yuc.Dt
info<- data.frame(sample=gsub("..$","",names(Yuc.At)), lib_size=colSums(Yuc.AnD) )
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
total<-Yuc.AnD
norm<-sweep(total,2,info$lib_size,"/")*10^6
pca=prcomp(t(log2(norm+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.polycat.log2cpm.Yuc.pdf")
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
save(Yuc.At, Yuc.Dt, file = "R-01-polycatDatasets.Yuc.RData")

### prepare network datasets of A2D5, A2.D5.tech and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- gsub("-",".", gsub("..$","",names(Yuc.At)))
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
Yuc          <- prepDatasets(Yuc.At,   Yuc.Dt  )
dim(Yuc)

# make list
networks  <- list(Yuc      = Yuc)

### rlog transformation
# need DESeq2
library(DESeq2)
rlogTransformation<-function(x)
{
    count <- round( x ) #have to round ADs.ncorrect
    coldata<-data.frame( sample = names(count), rep = gsub(".*[.]","", names(count)) )
    dds <- DESeqDataSetFromMatrix( countData = count, colData = coldata, design = ~ sample)
    # rlog transformation, note that with defalt blind=TRUE, design is actually ~1
    rld <- rlog(dds)
    count.rld <- as.data.frame(assay(rld))
    names(count.rld)<-names(count)
    return(count.rld)
}
networks.rld <- list(Yuc      = rlogTransformation(Yuc) )

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
# "Yuc"
for(i in 1:4) {print(table(is.na(as.matrix(networks.rpkm[[i]])) ) )}
for(i in 1:4) {print(table(is.infinite(as.matrix(networks.rpkm[[i]])) ) )}

#save
coldata<-info
save(coldata, networks, networks.rld, geneLen, networks.rpkm, file = "R-01-polycatNetworkDatasets.Yuc.RData")

###### Output files are:
# "pca.polycat.log2cpm.pdf"
# "table.count.polycat.txt"
# "R-01-polycatDatasets.RData"
# "R-01-polycatNetworkDatasets.RData"
