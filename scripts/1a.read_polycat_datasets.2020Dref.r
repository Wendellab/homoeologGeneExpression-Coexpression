############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################
# 2020 revision
# Read classification - old polycat [v1.3], new polycat [v2], old gsnap-polycat (original files)
# Counting - featureCounts, htseq, htseq primary alignments only

#FUN
readHTseqCounts<-function(dir){
    fileL<-grep("txt",list.files(dir),value=TRUE)
     # 3 species x 11 samples x 5(total, A, D, N, AD)  = 165
    for(file in fileL)
    {
        x <- read.table(paste0(dir,file), sep="\t")
        names(x) <- c("Geneid", file)
        if(!exists("allcount")) {allcount = x}
        else {allcount <- merge(allcount, x, by="Geneid")}
    }
    row.names(allcount) <- allcount$Geneid
    names(allcount)<-gsub("-",".",gsub(".txt", "", names(allcount) ))
    allcount<-allcount[grep("Gorai.0",rownames(allcount)),]
    return(allcount)
}


#### New GSNAP + New PolyCat + different counting
# gsnapD5.seed.newPolyCat.featureCounts = new gsnap against D5 ref using new polycat and feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.newPolyCat.featureCounts"
seed.NNF <- read.table(seed_file, sep="\t",header=T)
# htseq-count-new = new gsnap against D5 ref using new polycat and htseq (but not specifying --secondary-alignments=ignore)
seed_dir <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/htseq-count-new/"
seed.NNH=readHTseqCounts(seed_dir)
# htseq-count-new-primary-only = new gsnap against D5 ref using new polycat and htseq (specifying --secondary-alignments=ignore) * this should be more similar to featureCounts (edited)
seed_dir <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/htseq-count-new-primary-only/"
seed.NNHp=readHTseqCounts(seed_dir)


# gsnapD5.seed.oldPolyCat.featureCounts = new gsnap against D5 ref using old polycat version and feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.oldPolyCat.featureCounts"
seed.NOF <- read.table(seed_file, sep="\t",header=T)

# gsnapD5.seed.oldGSNAP.featureCounts = our original run; old gsnap against D5 with old polycat, but counted with feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.oldGSNAP.featureCounts"
# read files
seed.OOF <- read.table(seed_file, sep="\t",header=T)


##################
source("2.1.FUN.r")
quickEval=function(count){
    names(count) =gsub(".bam|.sort","",names(count))
    names(count)[grep("[1-9]$",names(count))] = paste0(grep("[1-9]$",names(count),value=T),".T")
    A2.Total  <- count[,grep('A2.*[.]T$',names(count))]
    ADs.Total <- count[,grep('AD.*[.]T$',names(count))]
    D5.Total  <- count[,grep('D5.*[.]T$',names(count))]
    ##################
    A2.At  <- count[,grep('A2.*[.]A$',names(count))]
    ADs.At <- count[,grep('AD.*[.]A$',names(count))]
    D5.At  <- count[,grep('D5.*[.]A$',names(count))]
    ##################
    A2.Dt  <- count[,grep('A2.*[.]D$',names(count))]
    ADs.Dt <- count[,grep('AD.*[.]D$',names(count))]
    D5.Dt  <- count[,grep('D5.*[.]D$',names(count))]
    ##
    total <- D5.Total + A2.Total
    totalA <- A2.Total
    totalD <- D5.Total
    assigned <- ADs.At + ADs.Dt
    assignedA <- ADs.At
    assignedD <- ADs.Dt
    efficiency <- assigned/total
    efficiencyA <- assignedA/totalA
    efficiencyD <- assignedD/totalD
    e <- c(sum(assigned)/sum(total), sum(assignedA)/sum(totalA), sum(assignedD)/sum(totalD))
    ##
    different <- abs(ADs.At-A2.Total)+abs(ADs.Dt-D5.Total)
    differentA <- abs(ADs.At-A2.Total)
    differentD <- abs(ADs.Dt-D5.Total)
    discrepancy <- different/total
    discrepancyA <- differentA/totalA
    discrepancyD <- differentD/totalD
    d <- c(sum(different)/sum(total), sum(differentA)/sum(totalA), sum(differentD)/sum(totalD))
    ##
    resBinA=binC(TP = A2.At, TN = D5.Dt, FP = D5.At ,FN = A2.Dt)
    resBinD=binC(TP = D5.Dt, TN = A2.At, FP = A2.Dt ,FN = D5.At)
    oa=unlist(resBinA$summary[1,])
    od=unlist(resBinD$summary[1,])
    ##
    me= c(e,d,oa[1],od[1],oa[2],od[2], oa[4],od[4],oa[3],oa[5])
    names(me) = c("Efficiency","Efficiency.At","Efficiency.Dt", "Discrepancy","Discrepancy.At","Discrepancy.Dt","Precision.At","Precision.Dt","Recall.At","Recall.Dt","F1.At","F1.Dt","Accuracy","MCC")
    return(me)
}
TotalLibSize =data.frame(OOF=colSums(seed.OOF[,grep("[1-9].sort.bam",names(seed.OOF))]), NOF=colSums(seed.NOF[,grep("[1-9]$",names(seed.NOF))]), NNH=colSums(seed.NNH[,grep("T$",names(seed.NNH))]), NNHp=colSums(seed.NNHp[,grep("T$",names(seed.NNHp))]))
k=c(rep("-",6),c(oa[1],od[1],oa[2],od[2], oa[4],od[4],oa[3],oa[5]))
eval=data.frame(OOF=quickEval(seed.OOF), NOF=quickEval(seed.NOF), NNH=quickEval(seed.NNH),NNHp=quickEval(seed.NNHp))
##################

# 2019 gsnap + old polycat + htseq
x<-read.table("table.count.polycat.2019.txt",header=TRUE,sep="\t")

dfA=data.frame(Nf=colSums(seed1[,grep("sort.A.bam",names(seed1))]), #Of=colSums(seed2[,grep("sort.A.bam",names(seed2))]),
Nh=colSums(seed3[,grep(".A$",names(seed3))]), Oh2019=colSums(x[,grep("seed.*A$",names(x))]))



flower_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/gsnapD5.flower.all.counts"
flower <- read.table(flower_file, sep="\t",header=T)

# merge put flower first
a<-merge(flower, seed3, by="Geneid")
names(a) <-gsub(".sort|.bam|bam[.]","",names(a))
# Prepare read count tables including only Chr genes
allcount13 <- a[grep("Gorai.0",a$Geneid),]
names(allcount13)[grep("[1-9]$",names(allcount13))]<-paste0(grep("[1-9]$",names(allcount13),value=T),".T")
names(allcount13)[grep("0",names(allcount13))]<-unlist(lapply(strsplit(names(allcount13)[grep("0",names(allcount13))],"[.]" ),function(x)paste0("seed",x[2],"-",x[1],"-",x[3],".",x[4])))
dim(allcount13) #37223 440
# write polycat count tables
write.table(allcount13, "table.count.polycat.txt",sep="\t", row.names=FALSE)

allcount13<-read.table("table.count.polycat.txt",header=TRUE,sep="\t")
# bring old count table from original analysis
x<-read.table("table.count.polycat.2019.txt",header=TRUE,sep="\t")
s2020=colSums(allcount13[,grep("T$",names(allcount13))])
s2019=colSums(x[,grep("T$",names(x))])
names(s2019)==names(s2020)
df=data.frame(s2019,s2020,diff =(s2020-s2019)/s2019)
colSums(allcount13[,grep("T$",names(allcount13))])/colSums(x[,grep("T$",names(x))])
# some mapping seems truncated LD7.A2.L1, LD7.D5.L3, LD9.A2.L2

rownames(allcount13)<-allcount13$Geneid
##################
A2.Total  <- allcount13[,grep('A2.*[.]T$',names(allcount13))]
ADs.Total <- allcount13[,grep('ADs.*[.]T$',names(allcount13))]
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
A2.N  <- allcount13[,grep('A2.*[.]unclass$',names(allcount13))]
ADs.N <- allcount13[,grep('AD.*[.]unclass$',names(allcount13))]
D5.N  <- allcount13[,grep('[.]D5.*[.]unclass$',names(allcount13))]
##################
A2.AD  <- allcount13[,grep('A2.*[.]chimera$',names(allcount13))]
ADs.AD <- allcount13[,grep('AD.*[.]chimera$',names(allcount13))]
D5.AD  <- allcount13[,grep('[.]D5.*[.]chimera$',names(allcount13))]


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
colMeans(ADs.Aportion)
ADs.AtN<- ADs.Total * ADs.Aportion
##################
ADs.Dportion <-ADs.Dt/( ADs.At+ ADs.Dt)
colMeans(ADs.Dportion)   # NA means At=0 and (At+Dt)=0
ADs.Dportion[is.na(ADs.Dportion)]<-0
colMeans(ADs.Dportion)
ADs.DtN<- ADs.Total * ADs.Dportion
##################

# double check files
(colSums(ADs.At) + colSums(ADs.Dt))/colSums(ADs.Total)
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
save(A2.Total, D5.Total, ADs.Total, A2.At, D5.At, ADs.At, A2.Dt, D5.Dt, ADs.Dt, ADs.AtN, ADs.DtN,  file = "R-01-polycatDatasets.RData")


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
geneLen<-read.table("D5.trueANDtheo.lengths",header=TRUE, sep="\t")
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
