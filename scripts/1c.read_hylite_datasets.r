############### Step 1a. Prepare hylite network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# HyLITE were ran for both seed and flower datasets together
dir<-"~/working/mapping_hylite/results082417/"
dir <- "/home/hugj2006/jfw-lab/Projects/Eflen/flowerTimeDataset/mapping_hylite/results082417/"

# get Total read count without partition
des <- "results082417"
total<-read.table(paste0(dir, des, ".expression.txt"),sep="\t",header=TRUE)
names(total)<-gsub("^A2-|^D5-|^ADs-","",gsub("[.]","-",names(total)))
rownames(total)<- total$GENE
total<-total[,-1]
# check samples
info<- data.frame(sample=names(total), lib_size=colSums(total) )
info<-cbind(info,t(as.data.frame(strsplit(as.character(info$sample),"-"))))
names(info)[3:5] <-c("tissue","genome","rep")


# check partitioned reads
fileL<-grep(paste0(des,".AD.*read.summary.txt"),list.files(dir),value=TRUE)
fileL
# 3 species x 11 samples  = 33
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t", header=TRUE)
    tag<-gsub("^ADs.","",gsub(paste0(des,"."),"",gsub(".read.summary.txt","",file)))
    print(tag)
    tt <- apply(x[,2:9],1,sum)
    # check if the sum of partition equals to total, somehow the last gene of total was not counted...
    print(unique( tt - total[,tag] ))
    print(unique( tt[1:37222] - total[1:37222,tag] ) )
}
# get total counts of last gene
BAMs<-grep("sorted.bam$",list.files(dir),value=TRUE)
t<-c()
for(bam in BAMs){
    t[bam]<-system(paste0("samtools idxstats ",dir, bam,"|grep 'Gorai.013G272700.1'|cut -f 3"),intern = TRUE)
}
names(t)<-gsub("^AD[.]","ADs.",gsub("-",".",gsub(".sorted.bam","",names(t))))
names(t)[grep("0",names(t))]<-unlist(lapply(strsplit(names(t)[grep("0",names(t))],"[.]" ),function(x)paste0("seed",x[2],".",x[1],".",x[3])))
totalN<-gsub("-",".",names(total))
total[37223,] <- as.numeric(t[totalN])

# prepare ADs.At, ADs.Dt, ADs.N and ADs.Total
rm(ADs.At)
rm(ADs.Dt)
rm(ADs.N)
rm(ADs.Total)
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t", header=TRUE)
    id <-gsub(".1$","",x$GENE)
    tag<-gsub("^ADs.","",gsub(paste0(des,"."),"",gsub(".read.summary.txt","",file)))
    print(tag)
    if(!exists("ADs.At")){
        ADs.At<-data.frame(x$A2 + x$A2.N)
        names(ADs.At)<-tag
        rownames(ADs.At)<-id}else{
        ADs.At<-cbind(ADs.At, x$A2 + x$A2.N)
        names(ADs.At)[ncol(ADs.At)]<-tag
        }
    if(!exists("ADs.Dt")){
        ADs.Dt<-data.frame(x$D5 + x$D5.N)
        names(ADs.Dt)<-tag
        rownames(ADs.Dt)<-id}else{
            ADs.Dt<-cbind(ADs.Dt, x$D5 + x$D5.N)
            names(ADs.Dt)[ncol(ADs.At)]<-tag
        }
    if(!exists("ADs.N")){
        ADs.N<-data.frame(rowSums(x[,-c(1:5)]))
        names(ADs.N)<-tag
        rownames(ADs.N)<-id}else{
            ADs.N<-cbind(ADs.N, rowSums(x[,-c(1:5)]))
            names(ADs.N)[ncol(ADs.N)]<-tag
        }
    if(!exists("ADs.Total")){
        ADs.Total<-data.frame(rowSums(x[,-1]))
        names(ADs.Total)<-tag
        rownames(ADs.Total)<-id}else{
        ADs.Total<-cbind(ADs.Total, rowSums(x[,-1]))
        names(ADs.Total)[ncol(ADs.Total)]<-tag
        }

}

# prepare total dataset
all<-total
row.names(all)<-gsub(".1$","",row.names(all))
A2.Total <- all[grep("A2",names(all),value=TRUE)]
D5.Total <- all[grep("-D5-",names(all),value=TRUE)]
ADs.Total0 <- all[grep("AD",names(all),value=TRUE)]

# check
colSums(A2.Total)/(colSums(D5.Total) + colSums(A2.Total) )
colSums(ADs.Total)/ colSums(ADs.Total0)
colSums(ADs.Total)/ ( colSums(ADs.At) + colSums(ADs.Dt) +colSums(ADs.N) )
colSums(ADs.Total0)/ (colSums(A2.Total) + colSums(D5.Total))
unique(ADs.Total ==(A2.Total+D5.Total) ) # yes!!!!!!!!!
(colSums(ADs.Dt)+ colSums(ADs.At))/colSums(ADs.Total)   # efficiency
colSums(ADs.At)/colSums(A2.Total)   # At efficiency
colSums(ADs.Dt)/colSums(D5.Total)   # Dt efficiency



##### Plot grouping
norm<-sweep(total,2,info$lib_size,"/")*10^6
pca=prcomp(t(log2(norm+1)))
dat = as.data.frame(pca$x)
library(ggplot2)
pdf("pca.hylite.log2cpm.pdf")
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
save(A2.Total, D5.Total, ADs.Total, ADs.At, ADs.Dt, ADs.N, file = "R-01-hyliteDatasets.RData")

### prepare network datasets of A2D5 and ADs, each with 74446 genes x 33 samples
# make consistent sample names for 33 samples
name <- paste(info[names(ADs.Total),3], gsub("^.","R",info[names(ADs.Total),5]),sep=".")
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
ADs          <- prepDatasets(ADs.At,   ADs.Dt  )
# make raw network list
networks  <- list(A2D5      = A2D5,
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
networks.rld <- list(A2D5 = rlogTransformation(A2D5), ADs = rlogTransformation(ADs) )

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
# "A2D5"      "ADs"
for(i in 1:2) {print(table(is.na(as.matrix(networks.rpkm[[i]])) ) )}
for(i in 1:2) {print(table(is.infinite(as.matrix(networks.rpkm[[i]])) ) )}

#save
coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, networks.rld, geneLen, networks.rpkm, file = "R-01-hyliteNetworkDatasets.RData")

# write hylite count tables - total 37223 x 99
write.table(total, "table.count.hylite.total.txt",sep="\t", row.names=TRUE)
# write hylite count tables - AD partitioned 74446 x 33
write.table(networks$ADs, "table.count.hylite.ADs.txt",sep="\t", row.names=TRUE)


###### Output files are:
# "pca.hylite.log2cpm.pdf"
# "table.count.hylite.total.txt"
# "table.count.hylite.ADs.txt"
# "R-01-hyliteDatasets.RData"
# "R-01-hyliteNetworkDatasets.RData"