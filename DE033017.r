### Differential expression analysis
# [EBSeq](https://www.biostat.wisc.edu/~kendzior/EBSEQ/) is an empirical Bayesian DE analysis tool developed in UW-Madison, can take variance due to read mapping ambiguity into consideration by grouping isoforms with parent gene’s number of isoforms. In addition, it is more robust to outliers. RSEM also includes EBSeq in its folder named ‘EBSeq’.

library(EBSeq)
library(DESeq2)

pairwiseDE.ebseq<-function(expr, contrast, libSize=NULL, path="") {
    # for 2 conditions, e.g.
    # geneMat <- as.matrix( expr[,coldata$sample %in% c("A2.10", "D5.10" )] )
    select <- which(col_con %in% contrast)
    geneMat <- as.matrix( expr[,select])
    # library size factor
    if (is.null(libSize)) {
        Sizes = MedianNorm(geneMat)
    }else{
    Sizes <-libSize[select]
    }
    print(Sizes)
    EBOut = EBTest(Data = geneMat, Conditions =as.factor(col_con[select]), sizeFactors=Sizes, maxround=5 )
    # checking convergency, the differences between the 4th and 5th iterations need to be less than 0.01 or 0.001.
    print(EBOut$Alpha)
    print(EBOut$Beta)
    print(EBOut$P)
    #Checking the model fit and other diagnostics: QQP data points should lie on the y = x line for both conditions; DenNHist check the density plot of empirical q's vs the simulated q's, good model fit is expected.
    #par(mfrow=c(1,2))
    #QQP(EBOut)
    #DenNHist(EBOut)
    PP = GetPPMat(EBOut)
    
    EBDERes=GetDEResults(EBOut, FDR=0.05)
    # DEfound - DE gene IDs
    # PPEE and PPDE - the posterior probabilities of being EE or DE for each gene.
    # Status - each gene's status called by EBSeq.
    print( table(EBDERes$Status))
    # calculate fold change
    GeneFC=PostFC(EBOut)
    # str(GeneFC) # PostFC , RealFC , Direction
    # PlotPostVsRawFC(EBOut,GeneFC)
    res= data.frame( Status = EBDERes$Status,PostFC=NA, RealFC=NA, Direction=NA,PPDE=NA )
    res[names(GeneFC$PostFC),"PostFC"] <- GeneFC$PostFC
    res[names(GeneFC$PostFC),"RealFC"] <- GeneFC$RealFC
    res[names(GeneFC$PostFC),"Direction"] <- GeneFC$Direction
    res[rownames(PP),"PPDE"] <- PP[,2]
    
    write.table(res, file=paste(path, "ebseq.",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}



pairwiseDE.deseq2<-function(expr, contrast, libSize=NULL, path="")
{
    ## RSEM output expected_count is not integer, the counts need to be rounded for DESeq2 analysis.
    expr <- round(expr)
    select <- which(col_con %in% contrast)
    geneMat <- as.matrix( expr[,select])
    dds <- DESeqDataSetFromMatrix( countData = geneMat, colData = coldata[select,], design = ~ condition)
    # library size factor
    if (!is.null(libSize)) {
        sizeFactors(dds) <-libSize[select]
    }
    print(sizeFactors(dds))
    
    res <- results(DESeq(dds))
    print( summary(res,alpha=.05) ) # print results
    write.table(res, file=paste(path, "deseq.",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}




##################
##### Run DE #####
##################

setwd("~/jfw-lab/Projects/Eflen_networks/")
rdatafiles<-grep("Net",grep("R-01-",list.files(),value=TRUE),invert=TRUE,value=TRUE)
rdatafiles
# "R-01-hyliteDatasets.RData"
# "R-01-polycatDatasets.RData"
# "R-01-rsemDatasets.RData"

coln <- c("A-10-R1","A-10-R2","A-10-R3","A-20-R1","A-20-R3","A-30-R1","A-30-R2","A-30-R3","A-40-R1","A-40-R2","A-40-R3","D-10-R1","D-10-R2","D-10-R3","D-20-R1","D-20-R3","D-30-R1","D-30-R2","D-30-R3","D-40-R1","D-40-R2","D-40-R3")
col_con <- gsub("-R.","",coln)
coldata<-data.frame( sample= coln, condition =col_con )
batch <- rbind( c("A-10", "D-10" ), c("A-20", "D-20" ), c("A-30", "D-30" ), c("A-40", "D-40" ) )


for(file in rdatafiles)
{
    # file<-"R-01-polycatDatasets.RData"
    print(file)
    lnames<-load(file)
    flag<-gsub("R-01-|Datasets.RData","",file)

    ## true exprected
    expr.A2D5 <- cbind(A2.Total, D5.Total)
    expr.ADs <- cbind(ADs.At, ADs.Dt)
    
    ## ebseq
    ADs.libsize<-rep(MedianNorm(ADs.At + ADs.Dt),2)
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.A2D5,x, path=paste0("DE/",flag,".exp.") ))
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.ADs, x, path=paste0("DE/",flag,".obs."),libSize = ADs.libsize ))
    
    ## deseq
    libTotal<-colSums(ADs.At + ADs.Dt)
    ADs.libsize<-rep(libTotal/mean(libTotal),2)
    apply(batch,1,function(x) pairwiseDE.deseq2(expr.A2D5,x, path=paste0("DE/",flag,".exp.") ))
    apply(batch,1,function(x) pairwiseDE.deseq2(expr.ADs, x, path=paste0("DE/",flag,".obs."),libSize = ADs.libsize ))
}

##################
#### Summary #####
##################

# Compare all results
DEs<-c("mapping","type","DE_method","compr","A_up","D_up")
resL<-grep("vs.*txt$",list.files("DE/"),value=TRUE)
for(res in resL)
{
    labels<-unlist(strsplit(res,"[.]") )[1:4]
    x<-read.table(file=paste0("DE/", res), sep="\t",header=TRUE)
    # "A>D" and "D>A"
    if(labels[3]=="deseq"){de<-c(table(x$log2FoldChange[x$padj<0.05&!is.na(x$padj)]<0))}
    if(labels[3]=="ebseq"){de<-c(table(x$PostFC[x$Status=="DE"]<1))}
    DEs<-rbind(DEs,c(labels,de))
}
print(DEs)
options(stringsAsFactors = FALSE)
DEsummary<-data.frame(DEs[-1,])
names(DEsummary)<-DEs[1,]
DEsummary$sig<-as.numeric(DEsummary$A_up) + as.numeric(DEsummary$D_up)
write.table(DEsummary,"DE/DE.summary.txt",row.names=FALSE,sep="\t")

# Evaluate sensitivity and specificity
source("~/scripts/addTrans.r")
eval<-c("compr","mapping","DE_method","exp.sig","obs.sig","both.sig","sensitivity","specificity")
pdf("DE/DEevals.pdf")
par(mfrow=c(2,2))
for(compr in c("A-10vsD-10","A-20vsD-20","A-30vsD-30","A-40vsD-40"))
{
    for(mapping in c("polycat","rsem","hylite"))
    {
        for(method in c("deseq","ebseq"))
        {
            file.obs<- grep(paste(mapping,"obs",method,compr,sep="."),resL,value=TRUE)
            file.exp<- grep(paste(mapping,"exp",method,compr,sep="."),resL,value=TRUE)
            obs<-read.table(file=paste0("DE/", file.obs), sep="\t",header=TRUE)
            exp<-read.table(file=paste0("DE/", file.exp), sep="\t",header=TRUE)
            if(method=="deseq")
            {
                obs$sig<- obs$padj<0.05
                exp$sig<- exp$padj<0.05
            }
            if(method=="ebseq")
            {
                obs$sig<- obs$Status=="DE"
                obs$log2FoldChange <- log2(obs$PostFC)
                exp$sig<- exp$Status=="DE"
                exp$log2FoldChange <- log2(exp$PostFC)
            }
            obs.sig<-length(which(obs$sig==TRUE) )
            exp.sig<-length(which(exp$sig==TRUE) )
            both.sig<-length(which(obs$sig==TRUE & exp$sig==TRUE) )
            both.not.sig<-length(which(obs$sig!=TRUE & exp$sig!=TRUE) )
            # sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
            # specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)
            sensitivity <- both.sig/exp.sig
            specificity <- both.not.sig/(37223-exp.sig)
            eval<-rbind(eval,c(compr,mapping, method, exp.sig,obs.sig,both.sig,sensitivity,specificity))
            
            # plot log2 fold changes
            colors <- rep(addTrans("grey",50), length(row.names(exp)))
            # color DE genes: both - red, deseq2 - green, ebs - blue, esle=grey
            colors[ obs$sig==TRUE ] <- addTrans("blue", 10)
            colors[ exp$sig==TRUE ] <- addTrans("green",10)
            colors[obs$sig==TRUE & exp$sig==TRUE ] <- addTrans("red",10)
            # check gene order: unique(rownames(exp)==rownames(obs))
            plot(exp$log2FoldChange,obs$log2FoldChange, pch=".", col=colors, main =paste(mapping,method,compr,sep=".") )
            legend("bottomright", inset=.05, title="Signicant DEs", c("obs-only","exp-only","both"), fill=c("blue","green","red") )

        }
    }
}
dev.off()
print(eval)
DEeval<-data.frame(eval[-1,])
names(DEeval)<-eval[1,]
write.table(DEeval,"DE/DE.evaluation.txt",row.names=FALSE,sep="\t")



##################
### ROC curves ###
##################

library(ROCR)

resultAUC <- c("mapping","DE_method","compr","AUC")
pdf("DE/DE.ROC.pdf")
for(mapping in c("polycat","rsem","hylite"))
{
    for(method in c("deseq","ebseq"))
    {
        plot(0,0,xlim=c(0,1),ylim=c(0,1),pch='.',xlab="False Positive Rate",ylab="True Positive rate", main=paste0(mapping,": ", method))
        abline(a=0, b= 1)
        i=1
        aucs1<-c()
        for(compr in c("A-10vsD-10","A-20vsD-20","A-30vsD-30","A-40vsD-40"))
        {
            i=i+1
            file.obs<- grep(paste(mapping,"obs",method,compr,sep="."),resL,value=TRUE)
            file.exp<- grep(paste(mapping,"exp",method,compr,sep="."),resL,value=TRUE)
            obs<-read.table(file=paste0("DE/", file.obs), sep="\t",header=TRUE)
            exp<-read.table(file=paste0("DE/", file.exp), sep="\t",header=TRUE)
            if(method=="deseq")
            {
                trues <- obs$padj<0.05 & !is.na(obs$padj)
                predictions<- rank(1-exp$padj)
            }
            if(method=="ebseq")
            {
                trues <- obs$Status=="DE"
                predictions<- rank(exp$PPDE)
            }
            
            # ROC analysis
            pred <- prediction(predictions = predictions,labels= trues)
            roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
            auc.perf = performance(pred, measure = "auc")
            auc <- round(as.numeric(auc.perf@y.values ),3)
            resultAUC <- rbind(resultAUC, c(mapping, method, compr,auc))
            aucs1 <- c(aucs1, auc)
            #  plot(roc.perf,main=paste(mapping,compr) )
            lines(roc.perf@x.values[[1]], roc.perf@y.values[[1]], col = i)
        }
        legend(0.6,0.2,legend=paste0( c("A-10vsD-10","A-20vsD-20","A-30vsD-30","A-40vsD-40")," AUC: ",aucs1),col = 2:5, pch="_",cex=1)
        
    }
}
dev.off()
print(resultAUC)
DEresultAUC<-data.frame(resultAUC[-1,])
names(DEresultAUC)<-resultAUC[1,]
write.table(DEresultAUC,"DE/DE.ROC.txt",row.names=FALSE,sep="\t")


libray(ggplot2)
DEeval<-read.table("DE/DE.evaluation.txt",header=TRUE,sep="\t")
DEeval[,-c(1:3)] <- apply(DEeval[,-c(1:3)],2,as.numeric)
pdf("DE/DE.performance.pdf")
ggplot(data=DEeval, aes(x=mapping, y=sensitivity, fill=DE_method)) +
geom_boxplot()
ggplot(data=DEeval, aes(x=mapping, y=specificity, fill=DE_method)) +
geom_boxplot()
# low sensitivity from Hylite, specificity are similiar
# EBseq shows higher specificity than DEseq, probably due to higher stringency indicated by fewer DEs.
DEsummary<-read.table("DE/DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEsummary, aes(x=type, y=sig, fill=DE_method)) + geom_boxplot() + facet_grid(.~compr)
# significant DEs: mapping - hylite, polycat and rsem are similiar
ggplot(data=DEsummary, aes(x=type, y=sig, fill=mapping)) + geom_boxplot() + facet_grid(.~compr)
# exp vs obs: fewer significant DEs in observed
DEauc<-read.table("DE/DE.ROC.txt",header=TRUE,sep="\t")
DEauc[,4] <- as.numeric(DEauc[,4])
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEauc, aes(x=mapping, y=AUC, fill=DE_method)) + geom_boxplot()
dev.off()
# ggplot(data=x, aes(x=mapping, y=sensitivity, fill=DE_method)) + geom_bar(colour="black", stat="identity",position=position_dodge())


            
