############### Step 1a. Prepare hylite network datasets  ###############
## unix command: R CMD BATCH s1.R
################

# HyLITE were ran for both seed and flower datasets together
dir<-"/work/LAS/jfw-lab/hugj2006/eflen2020/bowtie2hylite/resultsA2_011820/"
bam.dir<-"/work/LAS/jfw-lab/hugj2006/eflen2020/bowtie2hylite/resultsA2/"
# get Total read count without partition
des <- "resultsA2_011820"
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
    print(unique( tt[length(tt)] - total[length(tt),tag] ) )
}
# get total counts of last gene
lastgene=rownames(total)[nrow(total)]
BAMs<-grep("sorted.bam$",list.files(bam.dir),value=TRUE)
t<-c()
for(bam in BAMs){
    t[bam]<-system(paste0("samtools idxstats ",bam.dir, bam,"|grep '",lastgene,"'|cut -f 3"),intern = TRUE)
}
names(t)<-gsub("-",".",gsub(".sorted.bam","",names(t)))
names(t)[grep("0",names(t))]<-unlist(lapply(strsplit(names(t)[grep("0",names(t))],"[.]" ),function(x)paste0("seed",x[2],".",x[1],".",x[3])))
totalN<-gsub("-",".",names(total))
total[nrow(total),] <- as.numeric(t[totalN])

# prepare ADs.At, ADs.Dt, ADs.N and ADs.Total
rm(ADs.At)
rm(ADs.Dt)
rm(ADs.N)
rm(ADs.Total)
for(file in fileL)
{
    x <- read.table(paste0(dir,file), sep="\t", header=TRUE)
    id <-gsub("[.]1$","",x$GENE)
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
row.names(all)<-gsub("[.]1$","",row.names(all))
A2.Total <- all[grep("A2",names(all),value=TRUE)]
D5.Total <- all[grep("-D5-",names(all),value=TRUE)]
ADs.Total0 <- all[grep("ADs",names(all),value=TRUE)]

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
pdf("pca.hyliteAref.log2cpm.pdf")
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
save(A2.Total, D5.Total, ADs.Total, ADs.At, ADs.Dt, ADs.N, file = "R-01-hyliteArefDatasets.RData")

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



#save
coldata<-data.frame( sample = name, tissue = gsub("[.]R.","", name), rep = gsub(".*[.]","", name) )
save(coldata, networks, file = "R-01-hyliteArefNetworkDatasets.RData")

# write hylite count tables - total 37223 x 99
write.table(total, "table.count.hyliteAref.total.txt",sep="\t", row.names=TRUE)
# write hylite count tables - AD partitioned 74446 x 33
write.table(networks$ADs, "table.count.hyliteAref.ADs.txt",sep="\t", row.names=TRUE)


###### Output files are:
# "pca.hylite.log2cpm.pdf"
# "table.count.hylite.total.txt"
# "table.count.hylite.ADs.txt"
# "R-01-hyliteDatasets.RData"
# "R-01-hyliteNetworkDatasets.RData"
