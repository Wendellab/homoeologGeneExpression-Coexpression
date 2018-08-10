## Analysis by Guanjing Hu from Nov 5th, 2016
## updated from previous version: Eflen_Networks021616.r, Eflen_Networks051916.r, Eflen_Networks090816.r

ssh hugj2006@bigram.ent.iastate.edu
cd /home/hugj2006/jfw-lab/Projects/Eflen_networks
ln -s /home/hugj2006/jfw-lab/Papers/GroverEflen2015/count_files count_files021616
screen -S eflen
module load R
# Module name: R                          Version: 3.3.1
R
# start R analysis



## The analysis was coducted in following steps
## 1. Basic data processing and cleaning - input DESeq2 rlog table and trait table for sample clustering and detecting outliers.
## 2. Choosing the soft-thresholding power - default or powers making good fit of scale free topology.
## 3. Network Construction - single block, corType = "pearson" (not "bicor", see discussion below), networkType = "signed"
## 4. General network topology analysis - produce a few plots for exploration


############## Install WGCNA ####################
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
biocLite("impute")
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("WGCNA")
install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") )


sessionInfo()


############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################

#cd ~/jfw-lab/Projects/Eflen_networks/
getwd() # "/net/my.files.iastate.edu/ifs/isu/las/research/jfw-lab/Projects/Eflen_networks"
# get read count tables generated from HTseq-count, NOT counter
# htseq-count -f bam --stranded=no -r pos $j ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > htcount/$j.txt
BEGIN R
dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/seed_for_eflen_paper/polycatTests/htcountfiles/"
fileL<-list.files(dir)
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

# check countable reads
check<-cbind( 
colSums(allcount[,grep('.*sort.bam',names(allcount))] ),
colSums( allcount[grep("Gorai",allcount$gene),grep('.*sort.bam',names(allcount))] ),
t(allcount[grep("__",allcount$gene),grep('.*sort.bam',names(allcount))] ) )
colnames(check)[1:2]<-c("Mapped","Counted")
check<-as.data.frame(check)
check$percentageCount <- check$Counted/check$Mapped
check$Mapped[1:11]+check$Mapped[23:33] - check$Mapped[12:22]
write.table(check,file="/home/hugj2006/jfw-lab/Projects/Eflen/seed_for_eflen_paper/polycatTests/checkCounts.txt", sep="\t")

# Prepare read count tables including only Chr genes
allcount13 <- allcount[grep("Gorai.0",row.names(allcount)),]
dim(allcount13)

##################
A2.Total  <- allcount13[,grep('A2.*sort.bam',names(allcount13))]
ADs.Total <- allcount13[,grep('AD.*sort.bam',names(allcount13))]
D5.Total  <- allcount13[,grep('D5.*sort.bam',names(allcount13))]
##################
A2.At  <- allcount13[,grep('A2.*sort.A.bam',names(allcount13))]
ADs.At <- allcount13[,grep('AD.*sort.A.bam',names(allcount13))]
D5.At  <- allcount13[,grep('D5.*sort.A.bam',names(allcount13))]
##################
A2.Dt  <- allcount13[,grep('A2.*sort.D.bam',names(allcount13))]
ADs.Dt <- allcount13[,grep('AD.*sort.D.bam',names(allcount13))]
D5.Dt  <- allcount13[,grep('D5.*sort.D.bam',names(allcount13))]
##################
A2.N  <- allcount13[,grep('A2.*sort.N.bam',names(allcount13))]
ADs.N <- allcount13[,grep('AD.*sort.N.bam',names(allcount13))]
D5.N  <- allcount13[,grep('D5.*sort.N.bam',names(allcount13))]
##################
A2.AD  <- allcount13[,grep('A2.*sort.AD.bam',names(allcount13))]
ADs.AD <- allcount13[,grep('AD.*sort.AD.bam',names(allcount13))]
D5.AD  <- allcount13[,grep('D5.*sort.AD.bam',names(allcount13))]


##################
ADs.Aportion <-ADs.At/( ADs.At+ ADs.Dt )
summary(ADs.Aportion)   # NA means At=0 and (At+Dt)=0
ADs.Aportion[is.na(ADs.Aportion)]<-0
ADs.AtN<- ADs.Total * ADs.Aportion
##################
ADs.Dportion <-ADs.Dt/( ADs.At+ ADs.Dt)
summary(ADs.Dportion)   # NA means At=0 and (At+Dt)=0
ADs.Dportion[is.na(ADs.Dportion)]<-0
ADs.DtN<- ADs.Total * ADs.Dportion
##################

# double check files
colSums(A2.At)/colSums(A2.Total)
colSums(D5.Dt)/colSums(D5.Total)
(colSums(ADs.At) + colSums(ADs.Dt) +colSums(ADs.N) +colSums(ADs.AD))/colSums(ADs.Total) 
(colSums(ADs.AtN) + colSums(ADs.DtN ))/colSums(ADs.Total)  # close to 1

# save
save(A2.Total, D5.Total, ADs.Total, A2.At, D5.At, ADs.At, A2.Dt, D5.Dt, ADs.Dt, ADs.AtN, ADs.DtN, A2.N, D5.N, ADs.N, A2.AD, D5.AD, ADs.AD, file = "R-01-polycatDatasets.RData")


# prepare network datasets
name  <- c("dev10.R1", "dev10.R2", "dev10.R3", "dev20.R1", "dev20.R3", "dev30.R1", "dev30.R2", "dev30.R3", "dev40.R1", "dev40.R2", "dev40.R3")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    rownames(a) <- paste0(rownames(a) , "a")
    rownames(d) <- paste0(rownames(d) , "d")
    ad        <- rbind (a,d)
    return(ad)
}


A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
ADs.ncorrect <- prepDatasets(ADs.AtN,  ADs.DtN )


# double check content: 74446    12
dim(A2D5)
dim(A2D5.tech)
dim(ADs)
dim(ADs.ncorrect)

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
networks.rld <- list(A2D5      = rlogTransformation(A2D5),
                     A2D5.tech = rlogTransformation(A2D5.tech),
                     ADs       = rlogTransformation(ADs),
                  ADs.ncorrect = rlogTransformation(ADs.ncorrect)   )
# save list
networks  <- list(A2D5      = A2D5,
                  A2D5.tech = A2D5.tech,
                  ADs       = ADs,
               ADs.ncorrect = ADs.ncorrect   )

# coldata
coldata<-data.frame( sample = name, dpa = gsub("dev|[.]R.","", name), rep = gsub(".*[.]","", name) )
coldata
#    sample dpa rep
#  dev10.R1  10  R1
#  dev10.R2  10  R2
#  dev10.R3  10  R3
#  dev20.R1  20  R1
#  dev20.R3  20  R3
#  dev30.R1  30  R1
#  dev30.R2  30  R2
#  dev30.R3  30  R3
#  dev40.R1  40  R1
#  dev40.R2  40  R2
#  dev40.R3  40  R3

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
networks.rpkm$A2D5.tech.theo <- rpkm(networks$A2D5.tech, geneLen[,c("gene","theoretical" )])
networks.rpkm$ADs.theo       <- rpkm(networks$ADs,       geneLen[,c("gene","theoretical" )])
names(networks.rpkm)
# "A2D5"            "A2D5.tech"       "ADs"             "ADs.ncorrect"
# "ADs.theo"        "A2D5.tech.theo"
for(i in 1:6) {print(table(is.na(as.matrix(networks.rpkm[[i]])) ) )}
for(i in 1:6) {print(table(is.infinite(as.matrix(networks.rpkm[[i]])) ) )}
# emprical eflen corrected dataset has ~2% NA (=0/0) and ~0.03% Inf (e.g. =4/0), due to zero empirical length
# assign 0 to NA and Inf
for(i in 5:6) {
    networks.rpkm[[i]][is.na(as.matrix(networks.rpkm[[i]]))] <- 0
    networks.rpkm[[i]][is.infinite(as.matrix(networks.rpkm[[i]]))] <- 0  }

#save
save(coldata, networks, networks.rld, geneLen, networks.rpkm, file = "R-01-polycatNetworkDatasets.RData")



############### Step 1b. Prepare RSEM network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# get RSEM result files
dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2RSEM/"
fileL<-grep("genes.results", list.files(dir),value=TRUE)
fileL<-fileL[-grep("20-R2",fileL)]  # remove un-used
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
# save
save(count,rpkm, file = "R-01-rsemDatasets.RData")


# prepare network datasets
name  <- c("dev10.R1", "dev10.R2", "dev10.R3", "dev20.R1", "dev20.R3", "dev30.R1", "dev30.R2", "dev30.R3", "dev40.R1", "dev40.R2", "dev40.R3")
prepDatasets<-function(a,d)
{
    names(a) <- name
    names(d) <- name
    ad        <- rbind (a,d)
    return(ad)
}

# rpkm
all<-rpkm
rownames(all)<-gsub(".D","d", gsub(".A","a", count$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("D5",names(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("D5",names(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("AD",names(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("AD",names(all),value=TRUE)]
ADs.Total <- ADs.Dt+ ADs.At
# check
colSums(A2.At)/colSums(A2.Total)
colSums(D5.Dt)/colSums(D5.Total)
colSums(ADs.Dt)/colSums(ADs.Total)
# make network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks.rpkm  <- list(A2D5      = A2D5,
A2D5.tech = A2D5.tech,
ADs       = ADs   )


# raw
all<-count
rownames(all)<-gsub(".D","d", gsub(".A","a", count$gene_id))
# make datasets
A2.At <- all[grep("a$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("A2",names(all),value=TRUE)]
A2.Total <- A2.At+ A2.Dt
D5.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("D5",names(all),value=TRUE)]
D5.At <- all[grep("a$",rownames(all),value=TRUE), grep("D5",names(all),value=TRUE)]
D5.Total <- D5.Dt+ D5.At
ADs.At <- all[grep("a$",rownames(all),value=TRUE), grep("AD",names(all),value=TRUE)]
ADs.Dt <- all[grep("d$",rownames(all),value=TRUE), grep("AD",names(all),value=TRUE)]
# check
colSums(A2.At)/colSums(A2.Total)
colSums(D5.Dt)/colSums(D5.Total)
# make network datasets
A2D5         <- prepDatasets(A2.Total, D5.Total)
A2D5.tech    <- prepDatasets(A2.At,    D5.Dt   )
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make list
networks  <- list(A2D5      = A2D5,
                  A2D5.tech = A2D5.tech,
                  ADs       = ADs   )

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
networks.rld <- list(A2D5      = rlogTransformation(A2D5),
A2D5.tech = rlogTransformation(A2D5.tech),
ADs       = rlogTransformation(ADs)   )
# save list

# coldata
coldata<-data.frame( sample = name, dpa = gsub("dev|[.]R.","", name), rep = gsub(".*[.]","", name) )
coldata
#    sample dpa rep
#  dev10.R1  10  R1
#  dev10.R2  10  R2
#  dev10.R3  10  R3
#  dev20.R1  20  R1
#  dev20.R3  20  R3
#  dev30.R1  30  R1
#  dev30.R2  30  R2
#  dev30.R3  30  R3
#  dev40.R1  40  R1
#  dev40.R2  40  R2
#  dev40.R3  40  R3
save(coldata, networks, networks.rld, networks.rpkm, file = "R-01-rsemNetworkDatasets.RData")



############### Step 2. Check network datasets  ###############
## unix command: R CMD BATCH s2.R
################

polycat <- load("R-01-polycatNetworkDatasets.RData")
rld.polycat <- networks.rld
rpkm.polycat<- networks.rpkm
rsem    <- load("R-01-rsemNetworkDatasets.RData")
rld.rsem <- networks.rld
rpkm.rsem<- networks.rpkm




list2dataframe<-function(list)
{
    nn<- names(list);
    df<-cbind(as.numeric(as.matrix(list[[1]])), as.numeric(as.matrix(list[[2]])))
    for(i in 3:length(nn)){df<-cbind(df, as.numeric(as.matrix(list[[i]])))}
    df<-as.data.frame(df)
    names(df)<-nn
    return(df)
}
library(genefilter)
library(ggplot2)

# panel.smooth function is built in.
# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


# do some plots
pdf("s2.exploratory.plots.pdf")
for(i in c("rld.polycat","rpkm.polycat","rld.rsem","rpkm.rsem")){
    net.dat <- list2dataframe(get(i))
    pca = prcomp(t(net.dat))
    dat = as.data.frame(pca$x)
    dat$group= rownames(dat)
    print( ggplot(aes(PC1, PC2, color=group),data=dat) + geom_point() + ggtitle(i) )
    
    pairs(net.dat, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main=paste0("Expression Profile Scatterplot Matrix",i))
}
# so empir corrected profiles are more similiar to uncorrected rpkm than theoretical eflen corrected profiles??
dev.off()



############### Step 3. Prep for WGCNA  ###############
## unix command: R CMD BATCH s2.R
################

library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# enableWGCNAThreads();


rsem    <- load("R-01-rsemNetworkDatasets.RData")

################## RSEM rlog ##########################
# "coldata"       "networks"      "networks.rld"  "networks.rpkm"
nSets = 3
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
shortLabels = names(networks.rld) # "A2D5"         "A2D5.tech"    "ADs"
setLabels = paste(shortLabels, c("expected", "expected" ,"observed"), sep="_")
# "A2D5_expected"           "A2D5.tech_expected"   "ADs_observed"            "ADs.ncorrect_Npartition"
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = t(networks.rld[[1]]));
multiExpr[[2]] = list(data = t(networks.rld[[2]]));
multiExpr[[3]] = list(data = t(networks.rld[[3]]));
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 3
# $nGenes 74446
# $nSamples 11 11 11
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 8616 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Excluding genes from the calculation due to too many missing samples or zero variance.
# Update exprSize
checkSets(multiExpr)
# $nSets 3
# $nGenes 65159
# $nSamples 11 11 11
# $structureOK  TRUE
save(multiExpr,nSets, setLabels, shortLabels, file = "R-02-dataInput.rsem_rld.RData")


pdf(file = "s3.SampleClusteringS.rsem_rld.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;



################## RSEM rpkm, log2(x+1) ##########################
# "coldata"       "networks"      "networks.rld"  "networks.rpkm"
nSets = 3
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
shortLabels = names(networks.rpkm) # "A2D5"         "A2D5.tech"    "ADs"
setLabels = paste(shortLabels, c("expected", "expected" ,"observed"), sep="_")
# "A2D5_expected"           "A2D5.tech_expected"   "ADs_observed"            "ADs.ncorrect_Npartition"
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = t(log2(networks.rpkm[[1]]+1)));
multiExpr[[2]] = list(data = t(log2(networks.rpkm[[2]]+1)));
multiExpr[[3]] = list(data = t(log2(networks.rpkm[[3]]+1)));
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 3
# $nGenes 74446
# $nSamples 11 11 11
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 8616 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Excluding genes from the calculation due to too many missing samples or zero variance.
# Update exprSize
checkSets(multiExpr)
# $nSets 3
# $nGenes 66008
# $nSamples 11 11 11
# $structureOK  TRUE
save(multiExpr,nSets, setLabels, shortLabels, file = "R-02-dataInput.rsem_rpkm.RData")

pdf(file = "s3.SampleClusteringS.rsem_rpkm.log2.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;


polycat <- load("R-01-polycatNetworkDatasets.RData")
# "coldata"       "networks"      "networks.rld"  "geneLen"  "networks.rpkm"
################## polycat rlog ##########################
# "coldata"       "networks"      "networks.rld"  "networks.rpkm"
nSets = 4
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
shortLabels = names(networks.rld) # "A2D5"         "A2D5.tech"    "ADs"  "ADs.ncorrect"
setLabels = paste(shortLabels, c("expected", "expected" ,"observed", "Npartition"), sep="_")
# "A2D5_expected"           "A2D5.tech_expected"   "ADs_observed"            "ADs.ncorrect_Npartition"
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = t(networks.rld[[1]]));
multiExpr[[2]] = list(data = t(networks.rld[[2]]));
multiExpr[[3]] = list(data = t(networks.rld[[3]]));
multiExpr[[4]] = list(data = t(networks.rld[[4]]));
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 4
# $nGenes 74446
# $nSamples 11 11 11 11
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 8616 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Excluding genes from the calculation due to too many missing samples or zero variance.
# Update exprSize
checkSets(multiExpr)
# $nSets 4
# $nGenes 64276
# $nSamples 11 11 11 11
# $structureOK  TRUE
save(multiExpr,nSets, setLabels, shortLabels, file = "R-02-dataInput.polycat_rld.RData")

pdf(file = "s3.SampleClusteringS.polycat_rld.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;


################## polycat rpkm, log2(x+1) ##########################
# "coldata"       "networks"      "networks.rld"  "networks.rpkm"
nSets = 6
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
shortLabels = names(networks.rpkm) # "A2D5"         "A2D5.tech"    "ADs"  "ADs.ncorrect"
setLabels = paste(shortLabels, c("expected", "expected" ,"observed", "Npartition", "corrected","corrected"), sep="_")
# "A2D5_expected"           "A2D5.tech_expected"   "ADs_observed"            "ADs.ncorrect_Npartition" "A2D5.tech.theo_corrected" "ADs.theo_corrected"
# Form multi-set expression data
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = t(log2(networks.rpkm[[1]]+1)));
multiExpr[[2]] = list(data = t(log2(networks.rpkm[[2]]+1)));
multiExpr[[3]] = list(data = t(log2(networks.rpkm[[3]]+1)));
multiExpr[[4]] = list(data = t(log2(networks.rpkm[[4]]+1)));
multiExpr[[5]] = list(data = t(log2(networks.rpkm[[5]]+1)));
multiExpr[[6]] = list(data = t(log2(networks.rpkm[[6]]+1)));
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 6
# $nGenes 74446
# $nSamples 11 11 11 11
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 8616 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
# Excluding genes from the calculation due to too many missing samples or zero variance.
# Update exprSize
checkSets(multiExpr)
# $nSets 6
# $nGenes 64226
# $nSamples 11 11 11 11
# $structureOK  TRUE
save(multiExpr,nSets, setLabels, shortLabels, file = "R-02-dataInput.polycat_rpkm.RData")

pdf(file = "s3.SampleClusteringS.polycat_rpkm.log2.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;




#######################################
#  Choosing the soft-thresholding power: analysis of network topology

rm(list=ls())
rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

for(file in rdatafiles)
{
     print(file)
     load(file)
     tag<-gsub(".*Input[.]|[.].*","",file)

# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
# types<-c("unsigned", "signed", "signed hybrid")
    type <- "signed"
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type, blockSize=10000)[[2]])      }
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste0("s2.ChooseSoftThresholdPower_", tag,".pdf") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    dev.off()
    assign(paste0("powerTables.",tag),powerTables)
}

save( powerTables.polycat_rld, powerTables.polycat_rpkm, powerTables.rsem_rld, powerTables.rsem_rpkm, file = "R-02-chooseSoftThreshold.Rdata")
# blockSize(65830, rectangularBlocks = TRUE) can calculates suitable block size for Biocrunch
# I just use 10000 here

# Inspect "s2.ChooseSoftThresholdPower_???.pdf".
unlist(lapply(powerTables.rsem_rpkm, function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.8 & x$slope<0][1]} ) )
### polycat_rld: 24 28 24 24  =>28
### polycat_rpkm: 24 24 20 24 24 20  =>24
### rsem_rld: 30 28 26  => 30
### rsem_rpkm: 24,22,28 => 28

# ???? Why do RPKM datasets approximate power law well with b=1 ?????
# To reach over 80% scale free topology model fit, and also allows the connectivity stabilized.
Powers <- c(28,24,30,28)
names(Powers)<-c("polycat_rld","polycat_rpkm", "rsem_rld", "rsem_rpkm")


############### Step X.  Differential coexpression analysis  ###############
## nohup R CMD BATCH s4.networkConstruction.R &
################
library(DiffCorr)

rm(list=ls())
rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

rp <- rbind(c("Set", "X", "Y", "DCpairs" ) )
for(file in rdatafiles)
{
    print(file)
    load(file) # # multiExpr,nSets, setLabels, shortLabels
    tag<-gsub(".*Input[.]|[.].*","",file)
    
    # Analysis of differntial coexpression gene pairs
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2))
    {
        print(paste0("=== Comparing ", shortLabels[ pwset[1,i] ], " vs ", shortLabels[ pwset[2,i] ]))
        X<-t(multiExpr[[ pwset[1,i] ]]$data)
        Y<-t(multiExpr[[ pwset[2,i] ]]$data)
        outfile <- paste0("DC/",tag, ".", shortLabels[ pwset[1,i] ],"vs",shortLabels[ pwset[2,i] ],".res.txt" )
        
        # DiffCorr code cannot handle too big a correlation matrix
        # if( nrow(X) < 40000)
        # comp.2.cc.fdr(output.file=outfile, X, Y, threshold=0.05)
        res <- comp.2.big.cc(output.file=outfile,X,Y, 2000)
        rp <- rbind(rp, c( tag, shortLabels[ pwset[1,i] ], shortLabels[ pwset[2,i] ], nrow(res)   )  )
    }
    
}
DC<-rp
DC<-as.data.frame(rp[-1,])
names(DC)<-rp[1,]

n<-c()
for(file in rdatafiles)
{
    print(file)
    load(file) # # multiExpr,nSets, setLabels, shortLabels
    tag<-gsub(".*Input[.]|[.].*","",file)
    nGenes <- ncol(multiExpr[[1]]$data)
    
    # Analysis of differntial coexpression gene pairs
    pwset<-combn(nSets,2)
    for(i in 1:choose(nSets,2))
    {
        n<- c(n, nGenes)
    }
    
}
DC$nGenes <- n
DC$DCpercent <- as.numeric(as.character(DC$DCpairs))/(DC$nGenes*(DC$nGenes-1)/2)
write.table(DC, file<-"sX.differentialCoexpression.txt", row.names=FALSE, sep="\t")

DC[DC$X=="A2D5"&DC$Y=="ADs",]
#           Set    X   Y  DCpairs   DCpercent nGenes
#     rsem_rpkm A2D5 ADs 11504289 0.005280842  66008
#      rsem_rld A2D5 ADs  8665348 0.004082008  65159
#  polycat_rpkm A2D5 ADs  3277727 0.001589234  64226
#   polycat_rld A2D5 ADs  3059846 0.001481285  64276

#### 03072017

rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")
flag<-"A2D5vsADs.res."
for(file in rdatafiles)
{
    print(file)
    load(file) # # multiExpr,nSets, setLabels, shortLabels
    tag<-gsub(".*Input[.]|[.].*","",file)
    nGenes <- ncol(multiExpr[[1]]$data)
    print( table(gsub("Gorai.*00","",colnames(multiExpr[[1]]$data))) )
}
# "R-02-dataInput.rsem_rpkm.RData"
# a     d
# 32702 33306
# [1] "R-02-dataInput.rsem_rld.RData"
# a     d
# 32176 32983
# [1] "R-02-dataInput.polycat_rpkm.RData"
# a     d
# 32070 32156
# [1] "R-02-dataInput.polycat_rld.RData"
# a     d
# 32076 32200

# how many pairs with and between subgenome
[hugj2006@biocrunch DC]$ cat polycat_rld.A2D5vsADs.res.txt | awk '{print $1,$2}'|sed 's/Gorai...........//g'|sort | uniq -c
583208 a a
1381237 a d
1095401 d d
1 molecule X
[hugj2006@biocrunch DC]$ cat polycat_rpkm.A2D5vsADs.res.txt | awk '{print $1,$2}'|sed 's/Gorai...........//g'|sort | uniq -c
615534 a a
1487138 a d
1175055 d d
1 molecule X
[hugj2006@biocrunch DC]$ cat rsem_rpkm.A2D5vsADs.res.txt | awk '{print $1,$2}'|sed 's/Gorai...........//g'|sort | uniq -c
2969290 a a
6883220 a d
1651779 d d
1 molecule X
[hugj2006@biocrunch DC]$ cat rsem_rld.A2D5vsADs.res.txt | awk '{print $1,$2}'|sed 's/Gorai...........//g'|sort | uniq -c
2135851 a a
5255747 a d
1273750 d d
1 molecule X


############### Step 4.  Network Construction  ###############
## nohup R CMD BATCH s4.networkConstruction.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

ptm <- proc.time()

Powers <- c(28,24,30,28)
names(Powers)<-c("polycat_rld","polycat_rpkm", "rsem_rld", "rsem_rpkm")

rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

# rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.polycat_rpkm.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file)
    compr<-gsub("R-02-dataInput.|.RData","",file)
    powerEach = Powers[compr]
    print(paste0(file,", construct network with b = ",powerEach))
    print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    
    
    ###### calculate individual TOMs
    print("###### Calculate individual TOMs:")
    iTOMs = blockwiseIndividualTOMs(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    power = powerEach, networkType = "signed", TOMType = "signed",
    # Save individual TOMs?
    saveTOMs = TRUE,
    individualTOMFileNames = paste0(compr,".iTOM-%N-block.%b.RData")  )
    
    ###### calculate consensus modules
    print("###### Construct consensus networks:")
    cnet = blockwiseConsensusModules(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    power = powerEach, networkType = "signed", TOMType = "signed",
    # load previous TOMs
    individualTOMInfo = iTOMs,
    # Saving the consensus TOM
    saveConsensusTOMs = TRUE,
    consensusTOMFileNames = paste0(compr,".cTOM-Block%b.RData"),
    # Basic tree cut options
    deepSplit = 2,  #default, known to reasonable
    minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
    pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
    # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
    mergeCutHeight = 0.25,
    # others
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3)
    
    ###### calculate individual modules
    print("###### Construct individual networks:")
    # stupid blockwiseModules only load TOM rdata file with "TOM", not "tomDS"
    tomFiles<-grep(paste0(compr,".iTOM"),list.files(), value=TRUE)
    for(fl in tomFiles)
    {
        load(fl)
        TOM<-tomDS
        save(TOM, file=fl)
    }
    rm(TOM)
    collectGarbage()
    for(i in 1:nSets)
    {
        inet = blockwiseModules(
        # Input data
        multiExpr[[i]]$data,
        # Data checking options
        checkMissingData = TRUE,
        # Options for splitting data into blocks
        blocks =  iTOMs$blocks,
        #randomSeed = 12345,
        #maxBlockSize =  500,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
        # Network construction arguments: correlation options, use bicor instead of default pearson
        corType = "pearson",
        # Adjacency and topology overlap function options
        power = powerEach, networkType = "signed", TOMType = "signed",
        # load previous TOMs
        loadTOM = TRUE,
        saveTOMFileBase = paste0(compr,".iTOM-",i),
        # Basic tree cut options
        deepSplit = 2,  #default, known to reasonable
        minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
        pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
        # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
        mergeCutHeight = 0.25,
        # others
        reassignThreshold = 0,
        numericLabels = TRUE,
        verbose = 3)
        
        assign(paste0("inet",i), inet)
        
    }
    
    save(list=c( "iTOMs", "cnet", grep("inet.",ls(), value=TRUE)), file = paste0("R-04-buildNetwork.", compr,".RData"))
    collectGarbage()

}

proc.time() - ptm

###############bookmark
# With maxBlockSize = 20000, it took roughly a week



############### Step 5.  General network topology analysis  ###############
## nohup R CMD BATCH s5.networkTopology.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

plotRefModulePres<-function(mp, file= "")
{
    pdf(file)
    for(set in 2:nSets ){
        # specify the reference and the test networks
        ref=1; test = set
        Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
        Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
        # Look at the observed preservation statistics
        # Obs.PreservationStats
        # Z statistics from the permutation test analysis
        # Z.PreservationStats
        
        # Let us now visualize the data.
        modIDs = rownames(Obs.PreservationStats)
        modColors=labels2colors(order(as.numeric(modIDs) )-1 )
        moduleSize = Obs.PreservationStats$moduleSize
        # we will omit the grey module (background genes)
        # and the gold module (random sample of genes)
        selectModules = !(modColors %in% c("grey", "gold"))
        # Text labels for points
        point.label = modIDs[selectModules]
        # Composite preservation statistics
        medianRank=Obs.PreservationStats$medianRank.pres
        Zsummary=Z.PreservationStats$Zsummary.pres
        
        par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
        # plot medianRank versus module size
        plot(moduleSize[selectModules],medianRank[selectModules],col=1, bg=modColors[selectModules],
        pch = 21,main=paste("medianRank -",shortLabels[ref], "vs",shortLabels[test]),
        cex = 2, ylab ="medianRank",xlab="Module size", log="x")
        labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
        
        # plot Zsummary versus module size
        plot(moduleSize[selectModules],Zsummary[selectModules], col = 1, bg=modColors[selectModules],pch = 21,
        main=paste("Zsummary -",shortLabels[ref], "vs",shortLabels[test]),
        cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
        labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
        # Add threshold lines for Zsummary
        abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
    }
    dev.off()
}


ptm <- proc.time()

Powers <- c(28,24,30,28)
names(Powers)<-c("polycat_rld","polycat_rpkm", "rsem_rld", "rsem_rpkm")


rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file) # "multiExpr"   "nSets"       "setLabels"   "shortLabels"
    compr<-gsub("R-02-dataInput.|.RData","",file)
    powerEach = Powers[compr]
    print(paste0(file,", construct network with b = ",powerEach))
    print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    load(paste0("R-04-buildNetwork.", compr,".RData")) ->Names
    print(Names)  # "iTOMs" "cnet"  "inet1" "inet2" "inet3"
    
    ###  1.Comparing basic topological networks parameters
    ######################################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        gg<-net$goodGenes    #get the good genes descision
        
        ###  Get basic topological networks parameters
        print(paste("Start to extract network parameter for:",genome))
        # network concept
        adj <- adjacency(subDat, power = powerEach,type = "signed")   # calculate adjacency table
        # fundmentalNetworkConcept taking forever
        # net1C <- fundamentalNetworkConcepts(adj, GS = NULL)           # computes fundamental network concepts
        Size = dim(adj)[1]
        Connectivity = apply(adj, 2, sum)-1
        Density = sum(Connectivity)/(Size * (Size - 1))
        Centralization = Size * (max(Connectivity) - mean(Connectivity))/((Size -
        1) * (Size - 2))
        Heterogeneity = sqrt(Size * sum(Connectivity^2)/sum(Connectivity)^2 -
        1)
        # clustering coef takes forever to computate, skip
        # ClusterCoef = .ClusterCoef.fun(adj)
        fMAR = function(v) sum(v^2)/sum(v)
        MAR = apply(adj, 1, fMAR)
        ScaledConnectivity = Connectivity/max(Connectivity, na.rm = T)
        output = list(Connectivity = Connectivity, ScaledConnectivity = ScaledConnectivity,  MAR = MAR, Density = Density, Centralization = Centralization, Heterogeneity = Heterogeneity)
        param <- unlist(lapply(output,mean))
        
        ### Get module and dupllication related topology
        # how many modules detected??
        param["nModules"]<-length(unique(net$colors))
        # interal and crossing connections
        dt<-grep("d$",rownames(adj))
        size.dt<-length(dt)
        at<-grep("a$",rownames(adj))
        size.at<-length(at)
        param["at.density"] <- (sum(adj[at,at])-size.at)/(size.at*(size.at-1))
        param["dt.density"] <- (sum(adj[dt,dt])-size.dt)/(size.dt*(size.dt-1))
        param["ad.density"] <- mean(adj[at,dt])
        ### write results table
        if(!exists("Rtable"))
        {
            Rtable <- data.frame(param)
            names(Rtable)<-genome
        }else
        {Rtable[,genome] <- param }
        
        
        assign(paste0(genome,"Adj"),adj)
        assign(paste0(genome,"Concepts"),output)
        
    }
    #  Save basic topological networks parameters
    rownames(Rtable)[1:3]<- paste0("mean",rownames(Rtable)[1:3])
    write.table(Rtable, file = paste0("s5.parameters.",compr,".txt"),sep="\t")
    save( list=c( grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
    
    
    ###  2.Correlating node-specific network properties
    ###################################################
    pwset<-combn(nSets,2)
    
    # pdf(paste0("s5.correlatingParameters.",compr,".pdf") )
    # par(mfrow=c(2,2))
    # for( j in 1:3) # "Connectivity"       "ScaledConnectivity" "MAR"
    # {
    #    for(i in 1:ncol(pwset) )
    #    {   xnet <- get(paste0(shortLabels[ pwset[1,i]  ],"Concepts"))
    #        ynet <- get(paste0(shortLabels[ pwset[2,i]  ],"Concepts"))
    #        pp <- names(xnet)[j]
    #        corr <- cor.test(xnet[[j]], ynet[[j]])
    #        plot (xnet[[j]], ynet[[j]], main = paste(pp, ": cor=", round(corr$estimate,2), ", p=", corr$p.value, sep=""), xlab= shortLabels[pwset[1,i]], ylab = shortLabels[pwset[2,i]], type ="p", pch=16, col = rgb(0, 0, 0, 0.2) )
    #    }
    # }
    # dev.off()
    
    df_Connectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$Connectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$Connectivity)
    df_ScaledConnectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$ScaledConnectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$ScaledConnectivity)
    df_MAR <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$MAR, get(paste0(shortLabels[ 2  ],"Concepts"))$MAR)
    df_expr <- cbind(  as.numeric(as.matrix(multiExpr[[1]]$data)),  as.numeric(as.matrix(multiExpr[[2]]$data))  )
    for(i in 3:nSets)
    {
        df_Connectivity <- cbind( df_Connectivity, get(paste0(shortLabels[ i  ],"Concepts"))$Connectivity)
        df_ScaledConnectivity <- cbind( df_ScaledConnectivity , get(paste0(shortLabels[ i  ],"Concepts"))$ScaledConnectivity)
        df_MAR <- cbind( df_MAR, get(paste0(shortLabels[ i  ],"Concepts"))$MAR)
        df_expr <- cbind(  df_expr,  as.numeric(as.matrix(multiExpr[[i]]$data))  )
    }
    colnames(df_Connectivity) <-shortLabels
    colnames(df_ScaledConnectivity) <-shortLabels
    colnames(df_MAR) <-shortLabels
    colnames(df_expr) <-shortLabels
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    
    pdf(paste0("s5.correlatingParameters.",compr,".k.pdf") )
    pairs(df_Connectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".sk.pdf") )
    pairs(df_ScaledConnectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Scaled Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".mar.pdf") )
    pairs(df_MAR, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene MAR Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".expr.pdf") )
    pairs(df_expr, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Expression Scatterplot Matrix")
    dev.off()



    ###  3.WGCNA provided methods for comparison
    ############################################
    # get all anova p and put together
    dpa=as.factor(c(10,10,10,20,20,30,30,30,40,40,40))
    anovaP<-list()
    for(set in 1:nSets)
    {
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        MEs<-net$MEs
        pval<-apply(MEs,2,function(x){round(anova(aov(x~dpa) )$"Pr(>F)"[1],4)})
        pval<-as.data.frame(pval)
        pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
        pval$numeric<-as.numeric(substring(rownames(pval),3) )
        pval<-pval[order(pval$numeric),]
        pval$symbol[1]<-" "  # ME0 always meaningless
        anovaP[[set]]<-pval
    }
    names(anovaP)<-shortLabels
    
    # Module correspondence test
    pdf(paste0("s5.moduleCorrepondence.",compr,".pdf"),width=10,height=7)
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(8, 10.4, 2.7, 1)+0.3);
    # loop pairwise comparison
    for(i in 1:ncol(pwset))
    {
        print(paste0("Plot for pairwise set ",i))
        colnet<-get(paste0("inet",pwset[1,i]) )
        rownet<-get(paste0("inet",pwset[2,i]) )
        coln<-ncol(colnet$MEs )  # number of MEs
        rown<-ncol(rownet$MEs )
        # color list of MEs in the color of decreasing numbers of memebers
        colModule  <- labels2colors(as.numeric(names(table(colnet$colors)) ))
        rowModule  <- labels2colors(as.numeric(names(table(rownet$colors)) ))
        # colors for each gene
        colColors  <- labels2colors(colnet$colors )
        rowColors  <- labels2colors(rownet$colors )   # colors for each gene
        # anova significance sybol
        colP  <- anovaP[[pwset[1,i]]]$symbol
        rowP  <- anovaP[[pwset[2,i]]]$symbol
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # Execute all pairwaise comparisons
        for (rmod in 1:rown)
        for (cmod in 1:coln)
        {
            rMembers = (rowColors == rowModule[rmod] );
            cMembers = (colColors == colModule[cmod] );
            pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
            CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
        }
        
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rModTotals = apply(CountTbl, 1, sum)
        cModTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".",colModule, "_", cModTotals, colP),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".",rowModule, "_", rModTotals, rowP),
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
        
        # plotting only singficant modules
        labeledHeatmap( Matrix = pTable[rowP=="*",colP=="*"], colorLabels = TRUE,
        xLabels = paste(" ", colModule[colP=="*"]), yLabels = paste(" ", rowModule[rowP=="*"]),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".", colModule[colP=="*"], ": ", cModTotals[colP=="*"], colP[colP=="*"]),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".", rowModule[rowP=="*"], ": ", rModTotals[rowP=="*"], rowP[rowP=="*"]),
        textMatrix = CountTbl[rowP=="*",colP=="*"], colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules, significant only",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
    }
    dev.off()
    
    # Reference network based preservation test
    # set A2D5 as reference
    multiColor<-list()
    for(set in 1:nSets)
    {
        multiColor[[set]] <- get(paste0("inet",set) )$colors
    }
    names(multiColor) <-shortLabels
    # The number of permutations drives the computation time of the module preservation function. For a publication use 200 permutations.
    # But for brevity and testing, start with a small number
    nPermutations1=200
    # Set it to a low number (e.g. 3) if only the medianRank statistic and other observed statistics are needed.
    # Permutations are only needed for calculating Zsummary and other permutation test statistics.
    # set the random seed of the permutation test analysis
    set.seed(1)
    system.time({
        mp = modulePreservation(multiExpr, multiColor, networkType="signed", referenceNetworks = 1, testNetworks=NULL, nPermutations = nPermutations1,
         permutedStatisticsFile = paste0("s5.permutedStats-intrModule.", compr,".RData"),
        randomSeed = 1, quickCor = 0, verbose = 3)
    })
    # for 3 permutations, 4539.133s as roughly 1.26 hours
    # so it will probably takes 3.5 days to run 200 permutations
    save( list=c("anovaP", "mp", grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
    # plot preservation test results
    plotRefModulePres(mp, file= paste0("s5.refPreservation.", compr,".pdf"))
    
    
    ###  4.Homoeolog focused analysis
    ############################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        
        probes <- colnames(subDat )
        genes <- gsub(".$","",probes)
        subgenome <- gsub(".*00","",probes)
        homoeo <- c(length(intersect( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        colors <-net$colors
        inModuleNo<-aggregate(net$colors, by=list(genes), function(x)length(unique(x)))
        inSameModule <- length(which(inModuleNo$x==1))
        dominance<- as.data.frame(unclass(xtabs(~net$colors+subgenome) ) )
        dominance$chisq.p<-apply(dominance,1,function(x)chisq.test(x)$p.value)
        dominance$dominant <- apply(dominance,1,function(x)ifelse(x[3]>0.05,"-",ifelse(x[1]>x[2],"a","d")) )
        dominant <- factor(dominance$dominant, levels=c("a","b","-"))
        homoeo<-c(homoeo, tabulate(dominant) )
        names(homoeo)[4:6]<-levels(dominant)
        if(!exists("Rtableh"))
        { Rtableh<-data.frame(homoeo)
            names(Rtableh)<-genome
        }else {Rtableh[,genome]<-homoeo}

        assign(paste0(genome,"dominance"),dominance)
    }
    write.table(Rtableh, file = paste0("s5.homoeolog.",compr,".txt"),sep="\t")
    
    # save everything, AGAIN
    save( list=c("anovaP", "mp", grep(".dominance",ls(), value=TRUE), grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
}

proc.time() - ptm



############### Step END.  Extra analysis  ###############
## nohup R CMD BATCH s.networkTopology.R &
################

# Plot grapgh parameters with PCA
library(ggplot2)
pdf("s6.graphParametersPCA.pdf")
for(i in c("polycat_rld","polycat_rpkm","rsem_rld","rsem_rpkm"))
{
    f<-read.table(file = paste0("s5.parameters.",i,".txt"), header=TRUE, sep="\t" )
    pca = prcomp(t(f), scale.=T)
    
    # scale = 0 ensures that arrows are scaled to represent the loadings.
    # focus on the extreme ends (top, bottom, left, right), e.g. first principal component corresponds to a measure of left and right arrows.
    # For exact measure of a variable in a component, pca$rotation
    biplot(pca, scale = 0 , main=i)
    
    dat = as.data.frame(pca$x)
    dat$group= rownames(dat)
    portion.var <- as.numeric( summary(pca)$importance[2,] )
    print(
    ggplot(aes(PC1, PC2, color=group),data=dat) +
    geom_point() + ggtitle(i) +
    xlab(paste0("PC1: ",portion.var[1]*100,"% variance")) +
    ylab(paste0("PC2: ",portion.var[2]*100,"% variance"))
    )
}
# so empir corrected profiles are more similiar to uncorrected rpkm than theoretical eflen corrected profiles??
dev.off()



# Comparing parameters directly, get Z scores
rdatafiles<-c("R-05-networkTopology.rsem_rpkm.RData", "R-05-networkTopology.rsem_rld.RData", "R-05-networkTopology.polycat_rpkm.RData", "R-05-networkTopology.polycat_rld.RData")
Zresult <- rbind(c("Method","Obs","parameter","Zscore"))
for(f in rdatafiles)
{
    print(f)
    print("-----------------")
    
    rm(list=grep("Concepts",ls(),value=TRUE))
    ll<-load(f)
    
    ref= "A2D5Concepts"
    testG<-grep("Concepts",ls(), value=TRUE)
    testG<-testG[testG!=ref]
    nameG<-names(get(ref))
    for(test in testG)
    {
        result <- c(gsub(".RData|R-05-networkTopology.","",f), test)
        print(paste0(test," is compared to reference ", ref))
        for(i in 1:3)
        {
            x<- get(test)[[i]] - get(ref)[[i]]
            print(paste0("     ", nameG[i], ": mean - ", mean(x), "; sd - ", sd(x)))
            Zresult <- rbind(Zresult, c(result, nameG[i], x/sd(x)))

        }
    }
}
Z<-data.frame(Zresult[-1,])
names(Z) <- Zresult[1,]
Z$Obs<-gsub("Concepts","", Z$Obs)
Z$Zscore <- as.numeric(Z$Zscore)
Z$p.value <- pnorm(Z$Zscore)

write.table(Z, file = "s6.parameterZscore.txt",sep="\t", row.names=FALSE)


# "R-05-networkTopology.rsem_rpkm.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - 0.781465672186864; sd - 98.6522941737412"
# "     ScaledConnectivity: mean - 0.000947613568534769; sd - 0.0217144740877903"
# "     MAR: mean - 0.00450310922680791; sd - 0.0363252496455325"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -36.6065926997929; sd - 269.453978935194"
# "     ScaledConnectivity: mean - -0.00689189844026479; sd - 0.0593641940042516"
# "     MAR: mean - -0.00665182524599488; sd - 0.0741680503434127"
[1] "R-05-networkTopology.rsem_rpkm.RData"
[1] "-----------------"
[1] "A2D5.techConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - 0.785287565474176; sd - 57.5508927554109"
[1] "     ScaledConnectivity: mean - 0.00123342597671302; sd - 0.0199828561365023"
[1] "     MAR: mean - 0.00495548095961635; sd - 0.0432059016086636"
[1] "ADsConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -13.8622253567885; sd - 164.335781159567"
[1] "     ScaledConnectivity: mean - -0.00747326262845879; sd - 0.0568774292673812"
[1] "     MAR: mean - -0.00551956296440461; sd - 0.0883282894623805"


# "R-05-networkTopology.rsem_rld.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - 6.03432521331369; sd - 47.9331601209624"
# "     ScaledConnectivity: mean - 0.00385747478587089; sd - 0.0311011472237645"
# "     MAR: mean - 0.0108177203452433; sd - 0.0665238512411918"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -15.8272119289457; sd - 85.875353209918"
# "     ScaledConnectivity: mean - -0.00152867570149998; sd - 0.0564456719757156"
# "     MAR: mean - 0.000917384941257472; sd - 0.0845454737261381"


# "R-05-networkTopology.polycat_rpkm.RData"
# "-----------------"
# "A2D5.tech.theoConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -21.9512000285965; sd - 240.137834569537"
# "     ScaledConnectivity: mean - -0.00280878884515633; sd - 0.0466000052518468"
# "     MAR: mean - 0.00740427188747321; sd - 0.0679811075400122"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -21.9512000285965; sd - 240.137834569537"
# "     ScaledConnectivity: mean - -0.00280878884515633; sd - 0.0466000052518468"
# "     MAR: mean - 0.00740427188747321; sd - 0.0679811075400122"
# "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -24.3006652983133; sd - 202.22935432233"
# "     ScaledConnectivity: mean - -0.00126881936965774; sd - 0.0391716183150736"
# "     MAR: mean - 0.00560018426316752; sd - 0.060594947256564"
# "ADs.theoConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -23.3773020911422; sd - 247.50336056002"
# "     ScaledConnectivity: mean - -0.00385318188826418; sd - 0.0479840016318954"
# "     MAR: mean - 0.00369977480503477; sd - 0.0628610552950785"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -23.3773020911422; sd - 247.50336056002"
# "     ScaledConnectivity: mean - -0.00385318188826418; sd - 0.0479840016318954"
# "     MAR: mean - 0.00369977480503477; sd - 0.0628610552950785"
[1] "R-05-networkTopology.polycat_rpkm.RData"
[1] "-----------------"
[1] "A2D5.techConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -22.0765434162768; sd - 163.648523219436"
[1] "     ScaledConnectivity: mean - -0.00357149736531711; sd - 0.044308874996079"
[1] "     MAR: mean - 0.00925051099948535; sd - 0.0796255066248835"
[1] "A2D5.tech.theoConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -23.4762894834495; sd - 163.815572860369"
[1] "     ScaledConnectivity: mean - -0.00353351048896584; sd - 0.0443603130168054"
[1] "     MAR: mean - 0.00908386695385239; sd - 0.0799041069808706"
[1] "ADsConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -22.7441966541389; sd - 166.634032186635"
[1] "     ScaledConnectivity: mean - -0.00430959600795669; sd - 0.045165704424851"
[1] "     MAR: mean - 0.00491785060911022; sd - 0.073303344500602"
[1] "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -21.0670290871299; sd - 137.304169417926"
[1] "     ScaledConnectivity: mean - -0.00211100698492692; sd - 0.0370154162797463"
[1] "     MAR: mean - 0.00673809793888078; sd - 0.0712921811660626"
[1] "ADs.theoConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -24.1557787577917; sd - 166.67064086685"
[1] "     ScaledConnectivity: mean - -0.00426954496386292; sd - 0.0451968808000639"
[1] "     MAR: mean - 0.00474241482081538; sd - 0.0735389334852459"

# "R-05-networkTopology.polycat_rld.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -7.87451841852331; sd - 82.7567479726778"
# "     ScaledConnectivity: mean - -0.00128175364156708; sd - 0.0497757523553588"
# "     MAR: mean - 0.0126867937645758; sd - 0.0881149811226852"
# "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -4.52782204245409; sd - 73.9216718803422"
# "     ScaledConnectivity: mean - 0.00415615488438032; sd - 0.0450874147901161"
# "     MAR: mean - 0.00904897219457412; sd - 0.0806527881971959"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -11.1175246601803; sd - 81.7896357500773"
# "     ScaledConnectivity: mean - -0.00145226483274745; sd - 0.0492663618513508"
# "     MAR: mean - 0.00810894907197572; sd - 0.0825466056186121"



rdatafiles<-c("R-05-networkTopology.rsem_rpkm.RData", "R-05-networkTopology.rsem_rld.RData", "R-05-networkTopology.polycat_rpkm.RData", "R-05-networkTopology.polycat_rld.RData")
# view A2D5 vs A2D5.tech median rank
# mp$preservation$observed[[1]][[2]] ->mr
# mr[order(mr$medianRank.pres,decreasing=TRUE),1:2]

pdf("s6.sumZsummary.pdf")
Zlist<- c("Zsummary.pres", "Zdensity.pres", "Zconnectivity.pres")
upper<-c(100,160,60)
for(i in 1:3)
{
    z<-Zlist[i]
    u<-upper[i]
    load("R-05-networkTopology.polycat_rpkm.RData")
    plot(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="blue", ylim=c(0,u), xlab="Module Size", ylab="Zsummary", main=z)
    points(mp$preservation$Z[[1]][[4]]$moduleSize, mp$preservation$Z[[1]][[4]][,z], col="blue", pch = 2) #ncorrect
    # points(mp$preservation$Z[[1]][[6]]$moduleSize, mp$preservation$Z[[1]][[6]]$Zsummary.pres, col="blue", pch=0)
    
    load("R-05-networkTopology.polycat_rld.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="brown")
    points(mp$preservation$Z[[1]][[4]]$moduleSize, mp$preservation$Z[[1]][[4]][,z], col="brown", pch = 2) #ncorrect
    
    load("R-05-networkTopology.rsem_rpkm.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="green")
    load("R-05-networkTopology.rsem_rld.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="purple")
    
    abline(h=10,col="red")
    legend(100,u,legend=c("polycat_rpkm","polycat_rpkm: Npartition","polycat_rld","polycat_rld: Npartition","rsem_rpkm","rsem_rld"), col=c("blue","blue","brown","brown","green","purple"), pch=c(1,2,1,2,1,1))
    
}

#### Separability statistics performs poorly

load("R-05-networkTopology.polycat_rpkm.RData")
plot(mp$testSeparability$Z[[1]][[3]]$moduleSize,mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="blue", ylim=c(0,1500), xlab="Module Size", ylab="Zsummary", main="Z.separability.pres")
points(mp$testSeparability$Z[[1]][[4]]$moduleSize, mp$testSeparability$Z[[1]][[4]][,"Z.separability.pres"], col="blue", pch = 2) #ncorrect
# points(mp$preservation$Z[[1]][[6]]$moduleSize, mp$preservation$Z[[1]][[6]]$Zsummary.pres, col="blue", pch=0)

load("R-05-networkTopology.polycat_rld.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="brown")
points(mp$testSeparability$Z[[1]][[4]]$moduleSize, mp$testSeparability$Z[[1]][[4]][,"Z.separability.pres"], col="brown", pch = 2) #ncorrect

load("R-05-networkTopology.rsem_rpkm.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="green")
load("R-05-networkTopology.rsem_rld.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="purple")

legend(100,1500,legend=c("polycat_rpkm","polycat_rpkm: Npartition","polycat_rld","polycat_rld: Npartition","rsem_rpkm","rsem_rld"), col=c("blue","blue","brown","brown","green","purple"), pch=c(1,2,1,2,1,1))
abline(h=10,col="red")
dev.off()


# how homoeolog pairs in the same module

rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file) # "multiExpr"   "nSets"       "setLabels"   "shortLabels"
    compr<-gsub("R-02-dataInput.|.RData","",file)
    print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    load(paste0("R-04-buildNetwork.", compr,".RData")) ->Names
    print(Names)  # "iTOMs" "cnet"  "inet1" "inet2" "inet3"
    ###  4.Homoeolog focused analysis
    ############################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        
        probes <- colnames(subDat )
        genes <- gsub(".$","",probes)
        subgenome <- gsub(".*00","",probes)
        inPairs<-intersect( genes[subgenome =="a"], genes[subgenome=="d"] )
        homoeo <- c(length( inPairs ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        colors <-net$colors
        inModuleNo<-aggregate(net$colors, by=list(genes), function(x)length(unique(x)))
        rownames(inModuleNo)<-inModuleNo$Group.1
        inModuleNo<-inModuleNo[inPairs,]
        inSameModule <- length(which(inModuleNo$x==1))
        homoeo["ADinSameModule"]<-inSameModule
        if(set==1)
        {
            true<-inModuleNo;
            homoeo["ADinSM:sensitivity"]<-"-";
            homoeo["ADinSM:specificity"]<-"-"
        }else
        {
            truePair <- true$Group.1[true$x==1]
            obsPair  <- inModuleNo$Group.1[inModuleNo$x==1]
            # sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
            homoeo["ADinSM:sensitivity"]<- length(intersect(truePair,obsPair) )/length(truePair)
            
            truePairN <- true$Group.1[true$x>1]
            obsPairN  <- inModuleNo$Group.1[inModuleNo$x>1]
            # specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)
            homoeo["ADinSM:specificity"]<- length(intersect(truePairN,obsPairN) )/length(truePairN)
        }
        dominance<- as.data.frame(unclass(xtabs(~net$colors+subgenome) ) )
        dominance$chisq.p<-apply(dominance,1,function(x)chisq.test(x)$p.value)
        dominance$dominant <- apply(dominance,1,function(x)ifelse(x[3]>0.05,"-",ifelse(x[1]>x[2],"a","d")) )
        homoeo<-c(homoeo, unclass(table(dominance$dominant) ) )
        if(!exists("Rtableh"))
        { Rtableh<-data.frame(homoeo)
            names(Rtableh)<-genome
        }else {Rtableh[,genome]<-homoeo}
        
        assign(paste0(genome,"dominance"),dominance)
    }
    write.table(Rtableh, file = paste0("s5.homoeologN.",compr,".txt"),sep="\t")
    rm(Rtableh)
 }
# sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
# specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)


