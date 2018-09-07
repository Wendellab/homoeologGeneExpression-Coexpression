############### Step 1b. Prepare RSEM network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# RSEM result files for seeds and flowers
seed_dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2RSEM/"
flower_dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/flowerTimeDataset/mapping_rsem/"

# read seed count tables
dir<-seed_dir
fileL<-grep("genes.results", list.files(dir),value=TRUE)
fileL
# 3 species x 11 samples  = 33
rm(count)
rm(rpkm)
for ( file in fileL ) {
    x <- read.table(paste0(dir,file), header=TRUE, sep="\t")
    temp.count <- x[,c("gene_id", "expected_count")]
    temp.rpkm  <- x[,c("gene_id", "FPKM")]
    nn<- gsub(".genes.results","",file)
    names(temp.count)[2]<-nn
    names(temp.rpkm )[2]<-nn
    if (!exists("count")) { count<-temp.count }
    else count <- merge(count, temp.count, by="gene_id")
    if (!exists("rpkm"))  { rpkm <-temp.rpkm }
    else rpkm <- merge(rpkm, temp.rpkm, by="gene_id")
}
names(count)[-1]<-unlist(lapply(strsplit(names(count[-1]),"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )
names(rpkm)[-1]<-unlist(lapply(strsplit(names(rpkm[-1]),"-"),function(x)paste(paste0("seed",x[2]),x[1],x[3],sep="-")) )
countS<-count
rpkmS<-rpkm

# read flower count tables
dir<-flower_dir
fileL<-grep("genes.results", list.files(dir),value=TRUE)
fileL<-fileL[-grep("SD5-D5-S2|SD5-ADs-S2|SD5-A2-S4",fileL)]  # remove un-used
# 3 species x 22 samples  = 66
rm(count)
rm(rpkm)
for ( file in fileL ) {
    x <- read.table(paste0(dir,file), header=TRUE, sep="\t")
    temp.count <- x[,c("gene_id", "expected_count")]
    temp.rpkm  <- x[,c("gene_id", "FPKM")]
    nn<- gsub(".genes.results","",file)
    names(temp.count)[2]<-nn
    names(temp.rpkm )[2]<-nn
    if (!exists("count")) { count<-temp.count }
    else count <- merge(count, temp.count, by="gene_id")
    if (!exists("rpkm"))  { rpkm <-temp.rpkm }
    else rpkm <- merge(rpkm, temp.rpkm, by="gene_id")
}
# combine make flower first
count <- merge(count, countS, by="gene_id")
rpkm <- merge(rpkm, rpkmS, by="gene_id")

# write rsem count tables
write.table(count, "table.count.rsem.txt",sep="\t", row.names=FALSE)
write.table(rpkm, "table.rpkm.rsem.txt",sep="\t", row.names=FALSE)


# raw count
all<-count
rownames(all)<-gsub(".D","d", gsub(".A","a", count$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("-D5-",names(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("-D5-",names(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",names(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",names(all),value=TRUE)]
ADs.Total <- ADs.Dt+ ADs.At

##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
colSums(ADs.Dt)/colSums(ADs.Total)

##### inspect grouping by tissue types
total<-cbind(cbind(A2.Total, D5.Total),ADs.Total)
info<- data.frame(sample=names(total), lib_size=colSums(total) )
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
pdf("pca.rsem.log2cpm.pdf")
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
save(A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt, file = "R-01-rsemDatasets.RData")



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
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs   )


## make rpkm network datasets
all<-rpkm
rownames(all)<-gsub(".D","d", gsub(".A","a", count$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("-D5-",names(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("-D5-",names(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("ADs",names(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("ADs",names(all),value=TRUE)]
ADs.Total <- ADs.Dt+ ADs.At
##### inspect if diploid samples are truely diploid
colSums(A2.At)/colSums(A2.Total)
colSums(D5.At)/colSums(D5.Total)
colSums(ADs.At)/colSums(ADs.Total)
colSums(ADs.Dt)/colSums(ADs.Total)
##### inspect grouping by tissue types
total<-cbind(cbind(A2.Total, D5.Total),ADs.Total)
info<- data.frame(sample=names(total), lib_size=colSums(total) )
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"-"))))
names(info)[3:5] <-c("tissue","genome","rep")
# check samples complete
xtabs(~genome+tissue,data=info)
#      tissue
# genome LD7 LD9 LDM SD5 SD7 SDM SDP9 SDPF seed10 seed20 seed30 seed40
#    A2    3   3   3   2   3   2    3    3      3      2      3      3
#    ADs   3   3   3   2   3   2    3    3      3      2      3      3
#    D5    3   3   3   2   3   2    3    3      3      2
# make network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks.rpkm  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs   )


# rlog transformation
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
networks.rld <- list(A2D5 = rlogTransformation(networks$A2D5), A2D5.tech = rlogTransformation(networks$A2D5.tech), ADs = rlogTransformation(networks$ADs) )
# save list

# coldata
coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, networks.rpkm, file = "R-01-rsemNetworkDatasets.RData")


####### Output files are:
# "pca.rsem.log2cpm.pdf"
# "table.count.rsem.txt"
# "table.rpkm.rsem.txt"
# "R-01-rsemDatasets.RData"
# "R-01-rsemNetworkDatasets.RData"