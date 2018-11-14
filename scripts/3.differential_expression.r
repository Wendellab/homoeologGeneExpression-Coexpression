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
grep("Datasets.RData",grep("R-01-",list.files(),value=TRUE),value=TRUE)
# [1] "R-01-hyliteDatasets.RData"          "R-01-hyliteNetworkDatasets.RData"
# [3] "R-01-kallistoDatasets.RData"        "R-01-kallistoNetworkDatasets.RData"
# [5] "R-01-polycatDatasets.RData"         "R-01-polycatNetworkDatasets.RData"
# [7] "R-01-rsemDatasets.RData"            "R-01-rsemNetworkDatasets.RData"
# [9] "R-01-salmonDatasets.RData"          "R-01-salmonNetworkDatasets.RData"

methods = c("polycat","hylite","rsem","salmon","kallisto")

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
    if(labels[3]=="ebseq"){de<-c(table(x$PostFC[x$Status=="DE"]<1))}
    DEs<-rbind(DEs,c(labels,de))
}
print(DEs)
options(stringsAsFactors = FALSE)
DEsummary<-data.frame(DEs[-1,])
names(DEsummary)<-DEs[1,]
DEsummary$sig<-as.numeric(DEsummary$A_up) + as.numeric(DEsummary$D_up)
write.table(DEsummary,"s3.DE.summary.txt",row.names=FALSE,sep="\t")


#######################################
##### sensitivity and specificity #####
#######################################
# compare oobs and exp
DEsummary<-read.table("s3.DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
# batch<-cbind(paste0(unique(DEsummary$compr),".A"),paste0(unique(DEsummary$compr),".D"))
source("addTrans.r")
eval<-c("compr","mapping","DE_method","exp.sig","obs.sig","both.sig","sensitivity","specificity","precision","Fstat","MCC")
pdf("s3.DE.evals.pdf")
par(mfrow=c(3,2))
for(compr in paste0(batch[,1],"vs",batch[,2]) )
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
            # F measure
            Fstat =2*TP/(2*TP+FP+FN)
            # MCC
            MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            eval<-rbind(eval,c(compr,mapping, method, P,TP+FP,TP,sensitivity,specificity,precision,Fstat,MCC))
            
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
        for(compr in paste0(batch[,1],"vs",batch[,2]))
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
        legend(0.5,0.5,legend=paste0( paste0(batch[,1],"vs",batch[,2])," AUC: ",aucs1),col = 2:i, pch="_",cex=1)
        
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
DEeval<-read.table("s3.DE.evals.txt",header=TRUE,sep="\t")
DEeval[,-c(1:3)] <- apply(DEeval[,-c(1:3)],2,as.numeric)
DEauc<-read.table("s3.DE.ROC.txt",header=TRUE,sep="\t")
DEauc[,4] <- as.numeric(DEauc[,4])

# plots
library(ggplot2)
pdf("s3.DE.performance.pdf")
# significant DEs: smaller numbers (p=0) in EBseq than DEseq, fewer observed than expected p=0.03, no significance difference between partitioning methods.
ggplot(data=DEsummary, aes(x=mapping, y=sig, fill=type)) +scale_fill_brewer(palette=7) + geom_boxplot() + facet_grid(.~DE_method) + theme(legend.position="bottom") + ggtitle("Number of DE genes")
ggplot(data=DEsummary, aes(x=mapping, y=sig, fill=type)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom") + ggtitle("Number of DE genes, check across sample condition")
# binary classification
ggplot(data=DEeval, aes(x=mapping, y=sensitivity, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom") + ggtitle("Sensitivity")
ggplot(data=DEeval, aes(x=mapping, y=specificity, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom") + ggtitle("Specificity")
ggplot(data=DEeval, aes(x=mapping, y=precision, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom") + ggtitle("Precision")
ggplot(data=DEeval, aes(x=mapping, y=Fstat, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom") + ggtitle("F1 score")
ggplot(data=DEeval, aes(x=mapping, y=MCC, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom") + ggtitle("MCC")
# low sensitivity from Hylite, specificity are similiar
# EBseq shows higher specificity than DEseq, probably due to higher stringency indicated by fewer DEs.
## plot AUC
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEauc, aes(x=mapping, y=AUC, fill=DE_method)) + geom_boxplot()+ theme(legend.position="bottom") + ggtitle("AUC")
dev.off()

## plot DE numbers by comp
pdf("s3.DE.summary.pdf")
summary(DEsummary$sig/372.23) # use 10-50 range
for(compr in unique(DEsummary$compr) )
{
    p<-ggplot(data=DEsummary[ DEsummary$compr==compr,], aes(x=type, y=sig/372.23, fill=mapping)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() + scale_y_continuous(limits=c(10, 50))+ theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle(as.character(compr))
    print(p)
}
dev.off()


# ANOVA tests

# DE
fit<- lm(sig~mapping+DE_method+type,data=DEsummary)
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





            
