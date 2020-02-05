### Differential expression analysis
# [EBSeq](https://www.biostat.wisc.edu/~kendzior/EBSEQ/) is an empirical Bayesian DE analysis tool developed in UW-Madison, can take variance due to read mapping ambiguity into consideration by grouping isoforms with parent gene’s number of isoforms. In addition, it is more robust to outliers. RSEM also includes EBSeq in its folder named ‘EBSeq’.

library(EBSeq)
library(DESeq2)

pairwiseDE.ebseq<-function(expr, contrast, coldata, libSize=NULL, path="") {
    # for 2 conditions, e.g.
    # geneMat <- as.matrix( expr[,coldata$sample %in% c("A2.10", "D5.10" )] )
    select <- which(coldata$condition %in% contrast)
    geneMat <- as.matrix( expr[,select])
    # library size factor
    if (is.null(libSize)) {
        Sizes = MedianNorm(geneMat)
    }else{
    Sizes <-libSize[select]
    }
    print(Sizes)
    EBOut = EBTest(Data = geneMat, Conditions =as.factor(coldata$condition[select]), sizeFactors=Sizes, maxround=5 )
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



pairwiseDE.deseq2<-function(expr, contrast, coldata, libSize=NULL, path="")
{
    ## RSEM output expected_count is not integer, the counts need to be rounded for DESeq2 analysis.
    expr <- round(expr)
    select <- which(coldata$condition %in% contrast)
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

multiDE.deseq2<-function(expr, coldata, libSize=NULL, path="")
{
    ## RSEM output expected_count is not integer, the counts need to be rounded for DESeq2 analysis.
    geneMat <- as.matrix( round(expr))
    dds <- DESeqDataSetFromMatrix( countData = geneMat, colData = coldata, design = ~ tissue+genome)
    # library size factor
    if (!is.null(libSize)) {
        sizeFactors(dds) <-libSize
    }
    print(sizeFactors(dds))
    ddsMF=DESeq(dds)
    res <- results(ddsMF, contrast=c("genome","A","D"))
    print( summary(res,alpha=.05) ) # print results
    write.table(res, file=paste0(path, "deseqMF.AvsD.txt"), sep="\t")
}


pairwiseDE.fisher<-function(expr, contrast, coldata, libSize=NULL, path="")
{
    expr <- round(expr)
    select <- which(coldata$condition %in% contrast)
    geneMat <- as.matrix( expr[,select])
    
    # remove noise: >=1 count in all sample
    filter<- apply(geneMat,1,min)>0
    
    # note that fishers exact test doesn't require replicates, two options i have: 1.combine replicates then test; 2. test individually then
    # to be completed
    
    # write.table(res, file=paste(path, "fisher.",paste(contrast, collapse="vs"),".txt", sep=""), sep="\t")
}



##################
##### Run DE #####
##################

# check rdata files in place
list.files(pattern="R-01-.*Datasets.RData")
# [1] "R-01-hyliteDatasets.RData"          "R-01-hyliteNetworkDatasets.RData"
# [3] "R-01-kallistoDatasets.RData"        "R-01-kallistoNetworkDatasets.RData"
# [5] "R-01-polycatDatasets.RData"         "R-01-polycatNetworkDatasets.RData"
# [7] "R-01-rsemDatasets.RData"            "R-01-rsemNetworkDatasets.RData"
# [9] "R-01-salmonDatasets.RData"          "R-01-salmonNetworkDatasets.RData"

methods = c("polycat","hylite","hyliteAref","eaglerc","rsem","kallisto","salmon","bowtie")

for(flag in methods)
{
    load(paste0("R-01-",flag, "Datasets.RData")) # need A2.Total, D5.Total, ADs.At, ADs.Dt
    load(paste0("R-01-",flag, "NetworkDatasets.RData")) # need coldata
    print(flag)
    
    # print(names(A2.Total))
    # print(names(D5.Total))
    ## true exprected
    expr.A2D5 <- cbind(A2.Total, D5.Total)
    expr.ADs <- cbind(ADs.At, ADs.Dt)
    info<-rbind(coldata,coldata)
    info$genome<-rep(c("A","D"), each=33)
    info$condition<-paste0(info$tissue,".", info$genome)
    batch<-cbind(paste0(unique(info$tissue),".A"),paste0(unique(info$tissue),".D"))

    ## ebseq
    print("Run ebseq")
    ADs.libsize<-rep(MedianNorm(ADs.At + ADs.Dt),2)
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.A2D5,x, coldata=info, path=paste0("DE/",flag,".exp.") ))
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.ADs, x, coldata=info, path=paste0("DE/",flag,".obs."),libSize = ADs.libsize ))
    
    ## deseq
    print("Run deseq2")
    libTotal<-colSums(ADs.At + ADs.Dt)
    ADs.libsize<-rep(libTotal/mean(libTotal),2)
    apply(batch,1,function(x) pairwiseDE.deseq2(expr.A2D5,x, coldata=info, path=paste0("DE/",flag,".exp.") ))
    apply(batch,1,function(x) pairwiseDE.deseq2(expr.ADs, x, coldata=info, path=paste0("DE/",flag,".obs."),libSize = ADs.libsize ))
    
    ## deseq multifactor design
    print("Run deseq2 MF")
    multiDE.deseq2(expr.A2D5, coldata=info, path=paste0("DE/",flag,".exp.") )
    multiDE.deseq2(expr.ADs, coldata=info, path=paste0("DE/",flag,".obs."),libSize = ADs.libsize )
    ## fisher's exact test
    # use above libsize factor for deseq
    # skip
}

######################
##### DE Summary #####
######################
# Compare all results
DEs<-c("mapping","type","DE_method","compr","A_up","D_up")
resL<-grep("vs.*txt$",list.files("DE/"),value=TRUE)
for(res in resL)
{
    labels<-unlist(strsplit(res,"[.]") )[1:4]
    x<-read.table(file=paste0("DE/", res), sep="\t",header=TRUE)
    # "A>D" and "D>A"
    if(labels[3]=="deseq"){de<-c(table(x$log2FoldChange[x$padj<0.05&!is.na(x$padj)]<0))}
    if(labels[3]=="deseqMF"){de<-c(table(x$log2FoldChange[x$padj<0.05&!is.na(x$padj)]<0))}
    if(labels[3]=="ebseq"){de<-c(table(x$PostFC[x$Status=="DE"]<1))}
    DEs<-rbind(DEs,c(labels,de))
}
print(DEs)
options(stringsAsFactors = FALSE)
DEsummary<-data.frame(DEs[-1,])
names(DEsummary)<-DEs[1,]
DEsummary[,5:6] = apply(DEsummary[,5:6], 2, as.numeric)
DEsummary$sig<-DEsummary$A_up + DEsummary$D_up
write.table(DEsummary,"s3.DE.summary.txt",row.names=FALSE,sep="\t")



# check bias
DEsummary$chisqP = apply(DEsummary[,5:6],1, function(x)chisq.test(x)$p.value)
DEsummary$balance = ifelse(DEsummary$A_up-DEsummary$D_up>0, "A","D")
DEsummary$balance[DEsummary$chisqP>0.05] = "N"
unique(DEsummary[DEsummary$type=="exp",c(1,3,4)] == DEsummary[DEsummary$type=="obs",c(1,3,4)]) #TRUE
x<-DEsummary[DEsummary$type=="exp",c(1,3,4)]
x$exp = DEsummary$balance[DEsummary$type=="exp"]
x$obs = DEsummary$balance[DEsummary$type=="obs"]
x$pattern=factor(paste0(x$exp,"/",x$obs),levels=c("A/A","D/D","N/N","A/D","D/A","A/N","D/N","N/A","N/D"))
x$dataset= paste0(x$mapping,".",x$DE_method)
xtabs(~dataset+pattern, data=x)
#                 pattern
# dataset              A/A D/D N/N A/D D/A A/N D/N N/A N/D
#   bowtie.deseq         7   0   0   4   0   0   0   1   0
#   bowtie.deseqMF       0   0   0   0   0   0   0   0   1
#   bowtie.ebseq         0   4   1   2   2   0   0   2   1
#   eaglerc.deseq        0   1   2   0   2   0   2   3   2
#   eaglerc.deseqMF      0   0   0   1   0   0   0   0   0
#   eaglerc.ebseq        3   0   1   2   0   1   0   1   4
#   hylite.deseq         8   0   0   0   0   2   0   2   0
#   hylite.deseqMF       0   0   0   0   0   0   0   0   1
#   hylite.ebseq         0   4   0   2   0   0   1   1   4
#   hyliteAref.deseq     0   9   0   0   3   0   0   0   0
#   hyliteAref.deseqMF   1   0   0   0   0   0   0   0   0
#   hyliteAref.ebseq     4   0   0   2   1   0   0   4   1
#   kallisto.deseq       7   0   0   3   0   1   0   1   0
#   kallisto.deseqMF     0   0   0   0   0   0   0   0   1
#   kallisto.ebseq       0   3   2   2   1   0   0   2   2
#   polycat.deseq        6   0   0   1   0   2   0   2   1
#   polycat.deseqMF      0   0   0   1   0   0   0   0   0
#   polycat.ebseq        1   3   2   1   0   0   0   3   2
#   rsem.deseq          12   0   0   0   0   0   0   0   0
#   rsem.deseqMF         0   0   0   0   0   0   0   0   1
#   rsem.ebseq           0   2   0   2   0   0   0   0   8
#   salmon.deseq         8   0   0   3   0   1   0   0   0
#   salmon.deseqMF       0   0   0   0   0   0   0   0   1
#   salmon.ebseq         0   3   2   1   1   0   0   2   3
##### Take-home message, the pattern of balance is quite random between obs and exp

# obs vs exp
R<-c("mapping","DE_method","compr","obs","exp","overlap","union")
obsL<-grep("obs.*vs.*txt$",list.files("DE/"),value=TRUE)
for(i in obsL){
    labels<-unlist(strsplit(i,"[.]") )[1:4]
    obs = read.table(file=paste0("DE/", i), sep="\t",header=TRUE)
    exp = read.table(file=paste0("DE/", gsub("obs","exp",i)), sep="\t",header=TRUE)
    if(labels[3]=="deseq"){
        obsDEGs = rownames(obs)[obs$padj<0.05&!is.na(obs$padj)]
        expDEGs = rownames(exp)[exp$padj<0.05&!is.na(exp$padj)]
    }
    if(labels[3]=="deseqMF"){
        obsDEGs = rownames(obs)[obs$padj<0.05&!is.na(obs$padj)]
        expDEGs = rownames(exp)[exp$padj<0.05&!is.na(exp$padj)]
    }
    if(labels[3]=="ebseq"){
        obsDEGs = rownames(obs)[obs$Status=="DE"]
        expDEGs = rownames(exp)[exp$Status=="DE"]
    }
    R<- rbind(R, c(labels[-2], length(obsDEGs),length(expDEGs),length(intersect(obsDEGs,expDEGs)),length(union(obsDEGs,expDEGs))))
}
print(R)
write.table(R,"s3.DE.summary2.txt",row.names=FALSE,col.names=FALSE,sep="\t")

x<-read.table("s3.DE.summary2.txt",head=TRUE,sep="\t")
t.test(x$obs,x$exp, paired=T)


#######################################
##### sensitivity and specificity #####
#######################################
# compare oobs and exp
DEsummary<-read.table("s3.DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
batch<-cbind(paste0(unique(DEsummary$compr),".A"),paste0(unique(DEsummary$compr),".D"))
batch=gsub("AvsD.","",batch)
source("addTrans.r")
pdf("s3.DE.evals.pdf")
eval<-c("compr","mapping","DE_method","exp.sig","obs.sig","both.sig","sensitivity","specificity","precision","accuracy","Fstat","MCC")
par(mfrow=c(3,2))
for(compr in paste0(batch[1:12,1],"vs",batch[1:12,2]) )
{
    for(mapping in methods)
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
            P <- length(which(exp$sig==TRUE) ) # condition positive
            TP <-as.numeric(length(which(obs$sig==TRUE & exp$sig==TRUE) ))
            TN <-as.numeric(length(which(obs$sig!=TRUE & exp$sig!=TRUE) ))
            FP <-as.numeric(length(which(obs$sig==TRUE & exp$sig!=TRUE) ))
            FN <-as.numeric(length(which(obs$sig!=TRUE & exp$sig==TRUE) ))
            # sensitivity, or recall, true-positive rate (TPR), TP/P
            sensitivity <- TP/P
            # specificity),true-negative rate (TNR), TN/TN+FP  both.not.sig / (37223-A2D5.not)
            specificity <- TN/(TN+FP)
            # precision, postivie predictive value, TP/(TP+FP)
            precision <- TP/(TP+FP)
            # accuracy
            accuracy <-(TP+TN)/(TP+TN+FP+FN)
            # F measure
            Fstat =2*TP/(2*TP+FP+FN)
            # MCC
            MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            eval<-rbind(eval,c(compr,mapping, method, P,TP+FP,TP,sensitivity,specificity,precision,accuracy,Fstat,MCC))
            
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
print(eval)
compr="AvsD"
method="deseqMF"
for(mapping in methods)
{
    file.obs<- grep(paste(mapping,"obs",method,compr,sep="."),resL,value=TRUE)
    file.exp<- grep(paste(mapping,"exp",method,compr,sep="."),resL,value=TRUE)
    obs<-read.table(file=paste0("DE/", file.obs), sep="\t",header=TRUE)
    exp<-read.table(file=paste0("DE/", file.exp), sep="\t",header=TRUE)
    obs$sig<- obs$padj<0.05
    exp$sig<- exp$padj<0.05
    P <- length(which(exp$sig==TRUE) ) # condition positive
    TP <-as.numeric(length(which(obs$sig==TRUE & exp$sig==TRUE) ))
    TN <-as.numeric(length(which(obs$sig!=TRUE & exp$sig!=TRUE) ))
    FP <-as.numeric(length(which(obs$sig==TRUE & exp$sig!=TRUE) ))
    FN <-as.numeric(length(which(obs$sig!=TRUE & exp$sig==TRUE) ))
    # sensitivity, or recall, true-positive rate (TPR), TP/P
    sensitivity <- TP/P
    # specificity),true-negative rate (TNR), TN/TN+FP  both.not.sig / (37223-A2D5.not)
    specificity <- TN/(TN+FP)
    # precision, postivie predictive value, TP/(TP+FP)
    precision <- TP/(TP+FP)
    # accuracy
    accuracy <-(TP+TN)/(TP+TN+FP+FN)
    # F measure
    Fstat =2*TP/(2*TP+FP+FN)
    # MCC
    MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    eval<-rbind(eval,c(compr,mapping, method, P,TP+FP,TP,sensitivity,specificity,precision,accuracy,Fstat,MCC))
          
    # plot log2 fold changes
    colors <- rep(addTrans("grey",50), length(row.names(exp)))
    # color DE genes: both - red, deseq2 - green, ebs - blue, esle=grey
    colors[ obs$sig==TRUE ] <- addTrans("blue", 10)
    colors[ exp$sig==TRUE ] <- addTrans("green",10)
    colors[obs$sig==TRUE & exp$sig==TRUE ] <- addTrans("red",10)
    # check gene order: unique(rownames(exp)==rownames(obs))
    plot(exp$log2FoldChange,obs$log2FoldChange, pch=".", col=colors, main=paste(mapping,method,compr,sep=".") )
    legend("bottomright", inset=.05, title="Signicant DEs", c("obs-only","exp-only","both"), fill=c("blue","green","red") )
}
print(eval)
dev.off()

DEeval<-data.frame(eval[-1,])
names(DEeval)<-eval[1,]
write.table(DEeval,"s3.DE.evals.txt",row.names=FALSE,sep="\t")

###########################
### ROC curves: p-value ###
###########################
library(ROCR)

resultAUC <- c("mapping","DE_method","compr","AUC")
pdf("s3.DE.ROC.pdf")
methods # "polycat"  "hylite"   "rsem"     "salmon"   "kallisto"
for(mapping in methods)
{
    for(method in c("deseq","ebseq"))
    {
        plot(0,0,xlim=c(0,1),ylim=c(0,1),pch='.',xlab="False Positive Rate",ylab="True Positive rate", main=paste0(mapping,": ", method))
        abline(a=0, b= 1)
        i=1
        aucs1<-c()
        for(compr in paste0(batch[1:12,1],"vs",batch[1:12,2]))
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
        if(method=="deseq"){
            i=13
            compr="AvsD"
            file.obs<- grep(paste(mapping,"obs","deseqMF",compr,sep="."),resL,value=TRUE)
            file.exp<- grep(paste(mapping,"exp","deseqMF",compr,sep="."),resL,value=TRUE)
            obs<-read.table(file=paste0("DE/", file.obs), sep="\t",header=TRUE)
            exp<-read.table(file=paste0("DE/", file.exp), sep="\t",header=TRUE)
            trues <- obs$padj<0.05 & !is.na(obs$padj)
            predictions<- rank(1-exp$padj)
            # ROC analysis
            pred <- prediction(predictions = predictions,labels= trues)
            roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
            auc.perf = performance(pred, measure = "auc")
            auc <- round(as.numeric(auc.perf@y.values ),3)
            resultAUC <- rbind(resultAUC, c(mapping, method, compr,auc))
            aucs1 <- c(aucs1, auc)
            #  plot(roc.perf,main=paste(mapping,compr) )
            lines(roc.perf@x.values[[1]], roc.perf@y.values[[1]], col = i)
            legend(0.5,0.55,legend=paste0( paste0(batch[,1],"vs",batch[,2])," AUC: ",aucs1),col = 2:i, pch="_",cex=1)
        }else{
            legend(0.5,0.5,legend=paste0( paste0(batch[1:12,1],"vs",batch[1:12,2])," AUC: ",aucs1),col = 2:i, pch="_",cex=1)
        }
    }
}
dev.off()
print(resultAUC)
DEresultAUC<-data.frame(resultAUC[-1,])
names(DEresultAUC)<-resultAUC[1,]
write.table(DEresultAUC,"s3.DE.ROC.txt",row.names=FALSE,sep="\t")


###########################
###### Summary Plots ######
###########################

    
# import analysis results
DEsummary<-read.table("s3.DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
DEsummary$mapping = factor(DEsummary$mapping, levels=methods)
DEeval<-read.table("s3.DE.evals.txt",header=TRUE,sep="\t")
DEeval[,-c(1:3)] <- apply(DEeval[,-c(1:3)],2,as.numeric)
DEeval$mapping = factor(DEeval$mapping, levels=methods)
DEauc<-read.table("s3.DE.ROC.txt",header=TRUE,sep="\t")
DEauc[,4] <- as.numeric(DEauc[,4])
DEauc$mapping = factor(DEauc$mapping, levels=methods)

# prep table
aggregate(AUC~mapping+DE_method, range,data=DEauc)
aggregate(precision~mapping+DE_method, range,data=DEeval)
aggregate(sensitivity~mapping+DE_method, range,data=DEeval)
aggregate(specificity~mapping+DE_method, range,data=DEeval)
aggregate(accuracy~mapping+DE_method, range,data=DEeval)
aggregate(Fstatus~mapping+DE_method, range,data=DEeval)
aggregate(MCC~mapping+DE_method, range,data=DEeval)

#

# plots
library(ggplot2)
library(reshape2)
pdf("s3.DE.performance.pdf")
# significant DEs: smaller numbers (p=0) in EBseq than DEseq, fewer observed than expected p=0.03, no significance difference between partitioning methods.
ggplot(data=DEsummary[DEsummary$DE_method!="deseqMF",], aes(x=mapping, y=sig, fill=type)) +scale_fill_brewer(palette=7) + geom_boxplot() + facet_grid(.~DE_method) + theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Number of DE genes")
ggplot(data=DEsummary[DEsummary$DE_method!="deseqMF",], aes(x=mapping, y=sig, fill=type)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Number of DE genes, check across sample condition")
ss=rbind(melt(DEeval[DEeval$DE_method!="deseqMF",],id.vars=c("mapping","DE_method","compr"),measure.vars = c("precision","sensitivity","accuracy","Fstat","MCC")), melt(DEauc,id.vars=c("mapping","DE_method","compr"),measure.vars = c("AUC")) )
ggplot(data=ss, aes(x=mapping, y=value, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Performance measures") + facet_wrap(~variable)
ggplot(data=ss, aes(x=mapping, y=value, fill=DE_method)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Performance measures") + facet_wrap(~variable)
dev.off()
# low sensitivity from Hylite, specificity are similiar
# EBseq shows higher specificity than DEseq, probably due to higher stringency indicated by fewer DEs.
## plot AUC
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
pdf("s3.DE.performance.fig.pdf")
ggplot(data=ss[ss$mapping!="hyliteAref",], aes(x=mapping, y=value, fill=DE_method)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Performance measures") + facet_wrap(~variable)
dev.off()


# more DE genes by DESeq2 than EBseq
x<-read.table("s3.DE.summary.txt",head=TRUE,sep="\t")
t.test(x$sig[x$DE_method=="deseq"],x$sig[x$DE_method=="ebseq"], paired=T)
# t = 12.529, df = 191, p-value < 2.2e-16
summary(x$sig[x$DE_method=="deseq"]/x$sig[x$DE_method=="ebseq"])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.615   1.027   1.204   1.159   1.301   1.498

## plot percentage DEs by comp
DEsummary$sigPerc= DEsummary$sig/372.23 # use 10-50 range
DEsummary$sigPerc[DEsummary$mapping=="eaglerc"]= DEsummary$sig[DEsummary$mapping=="eaglerc"]/284.01
summary(DEsummary$sigPerc)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 8.183  20.710  25.112  26.125  30.469  54.493
summary(DEsummary$sigPerc[DEsummary$DE_method!="deseqMF"])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 8.183  20.593  24.799  25.049  30.086  44.236


pdf("s3.DE.summary.pdf")
for(compr in unique(DEsummary$compr) )
{
    p<-ggplot(data=DEsummary[ DEsummary$compr==compr,], aes(x=type, y=sigPerc, fill=mapping)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() + scale_y_continuous(limits=c(10, 50))+ theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle(as.character(compr))
    print(p)
}
dev.off()

for(compr in unique(DEsummary$compr)){
    x=DEsummary[ DEsummary$compr==compr,]
    print(compr)
    print(t.test(x$sig[x$type=="obs"],x$sig[x$type=="exp"], paired=T))
}
x
#        mapping type DE_method compr  A_up  D_up   sig  sigPerc
# 13      bowtie  exp   deseqMF  AvsD 10125 10142 20267 54.44752
# 38      bowtie  obs   deseqMF  AvsD  9497  9878 19375 52.05115
# 63     eaglerc  exp   deseqMF  AvsD  8202  7104 15306 53.89247
# 88     eaglerc  obs   deseqMF  AvsD  6883  8444 15327 53.96641
# 113     hylite  exp   deseqMF  AvsD  9888  9732 19620 52.70935
# 138     hylite  obs   deseqMF  AvsD  6661 12292 18953 50.91744
# 163 hyliteAref  exp   deseqMF  AvsD  9860  7658 17518 47.06230
# 188 hyliteAref  obs   deseqMF  AvsD 10418  5965 16383 44.01311
# 213   kallisto  exp   deseqMF  AvsD 10039 10065 20104 54.00962
# 238   kallisto  obs   deseqMF  AvsD  8734 10400 19134 51.40370
# 263    polycat  exp   deseqMF  AvsD 10312  9700 20012 53.76246
# 288    polycat  obs   deseqMF  AvsD  9006 10384 19390 52.09145
# 313       rsem  exp   deseqMF  AvsD  9618  9843 19461 52.28219
# 338       rsem  obs   deseqMF  AvsD  6504 12946 19450 52.25264
# 363     salmon  exp   deseqMF  AvsD 10125 10159 20284 54.49319
# 388     salmon  obs   deseqMF  AvsD  8798 10497 19295 51.83623
t.test(x$sig[x$type=="obs"],x$sig[x$type=="exp"], paired=T)
#    Paired t-test
# data:  x$sig[x$type == "obs"] and x$sig[x$type == "exp"]
# t = -4.2073, df = 7, p-value = 0.003999
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -1028.0066  -288.2434
# sample estimates:
# mean of the differences
#                -658.125


# ANOVA tests

# DE
fit<- lm(sig~mapping+DE_method+type,data=DEsummary[DEsummary$DE_method!="deseqMF",])
anova(fit)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

# sensitivity
fit<- lm(sensitivity~mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

# specificity
fit<- lm(specificity~mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

# precision
fit<- lm(precision~mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

# Fstat
fit<- lm(Fstat~mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

# MCC
fit<- lm(MCC~mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)

## AUC
fit<- lm(AUC~mapping+DE_method,data=DEauc)
summary(a1<-aov(fit))
TukeyHSD(x=a1, conf.level=0.95)
