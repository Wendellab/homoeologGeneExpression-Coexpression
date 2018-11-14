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
grep("Datasets.Yuc.RData",grep("R-01-",list.files(),value=TRUE),value=TRUE)
# "R-01-polycatNetworkDatasets.Yuc.RData"



for(flag in c("polycat"))
{
    load(paste0("R-01-",flag, "Datasets.Yuc.RData")) # need A2.Total, D5.Total, ADs.At, ADs.Dt
    load(paste0("R-01-",flag, "NetworkDatasets.Yuc.RData")) # need coldata
    print(flag)
    
    # print(names(A2.Total))
    # print(names(D5.Total))
    ## true exprected
    expr.ADs <- cbind(Yuc.At, Yuc.Dt)
    info<-rbind(coldata,coldata)
    info$genome<-rep(c("A","D"), each=ncol(Yuc.At))
    info$condition<-paste0(info$tissue,".", info$genome)
    rownames(info)<-names(expr.ADs)
    batch<-cbind(paste0(unique(info$tissue),".A"),paste0(unique(info$tissue),".D"))

    ## ebseq
    print("Run ebseq")
    ADs.libsize<-rep(MedianNorm(Yuc.At + Yuc.Dt),2)
    apply(batch,1,function(x) pairwiseDE.ebseq(expr.ADs, x, coldata=info, path=paste0("DE/Yuc.",flag,".obs."),libSize = ADs.libsize ))
    
    ## deseq
    print("Run deseq2")
    libTotal<-colSums(Yuc.At + Yuc.Dt)
    ADs.libsize<-rep(libTotal/mean(libTotal),2)
    
    apply(batch,1,function(x) pairwiseDE.deseq2(expr.ADs, x, coldata=info, path=paste0("DE/Yuc.",flag,".obs."),libSize = ADs.libsize ))
    
    ## fisher's exact test
    # use above libsize factor for deseq
    # skip
}

######################
##### DE Summary #####
######################

# Compare all results
DEs<-c("mapping","type","DE_method","compr","A_up","D_up")
resL<-grep("Yuc.*vs.*txt$",list.files("DE/"),value=TRUE)
for(res in resL)
{
    labels<-unlist(strsplit(res,"[.]") )[2:5]
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
write.table(DEsummary,"s3.DE.summary.Yuc.txt",row.names=FALSE,sep="\t")
---book
x<-read.table("s3.DE.summary.txt",sep="\t",header=TRUE)
fit<- lm(sig~compr+mapping+DE_method+type,data=x)
summary(aov(fit))
# Df    Sum Sq   Mean Sq F value Pr(>F)
# compr        11 803465226  73042293  70.732 <2e-16 ***
# mapping       2   9340192   4670096   4.522 0.0127 *
# DE_method     1 100793234 100793234  97.606 <2e-16 ***
# type          1    945594    945594   0.916 0.3404
# Residuals   128 132180023   1032656
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(a1<-aov(fit))
TukeyHSD(x=a1, 'DE_method', conf.level=0.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $DE_method
# diff       lwr       upr p adj
# ebseq-deseq -1673.264 -2008.384 -1338.144     0
#
TukeyHSD(x=a1, 'mapping', conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $mapping
#                     diff        lwr       upr     p adj
# polycat-hylite -478.3333  -970.2108  13.54413 0.0585413
# rsem-hylite    -585.9792 -1077.8566 -94.10171 0.0150617
# rsem-polycat   -107.6458  -599.5233 384.23163 0.8622411
#
TukeyHSD(x=a1, 'compr', conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $compr
#                              diff        lwr        upr     p adj
# LD9-LD7       -1713.75000 -3095.01022  -332.48978 0.0036122
# LDM-LD7        3634.41667  2253.15645  5015.67688 0.0000000
# SD5-LD7       -4335.25000 -5716.51022 -2953.98978 0.0000000
# SD7-LD7         563.66667  -817.59355  1944.92688 0.9690675
# SDM-LD7        -786.16667 -2167.42688   595.09355 0.7600944
# SDP9-LD7         90.08333 -1291.17688  1471.34355 1.0000000
# SDPF-LD7       2455.91667  1074.65645  3837.17688 0.0000018
# seed10-LD7     1320.66667   -60.59355  2701.92688 0.0751066
# seed20-LD7    -1143.08333 -2524.34355   238.17688 0.2123069
# seed30-LD7     1235.91667  -145.34355  2617.17688 0.1270631
# seed40-LD7    -4622.66667 -6003.92688 -3241.40645 0.0000000
# LDM-LD9        5348.16667  3966.90645  6729.42688 0.0000000
# SD5-LD9       -2621.50000 -4002.76022 -1240.23978 0.0000003
# SD7-LD9        2277.41667   896.15645  3658.67688 0.0000132
# SDM-LD9         927.58333  -453.67688  2308.84355 0.5272466
# SDP9-LD9       1803.83333   422.57312  3185.09355 0.0016056
# SDPF-LD9       4169.66667  2788.40645  5550.92688 0.0000000
# seed10-LD9     3034.41667  1653.15645  4415.67688 0.0000000
# seed20-LD9      570.66667  -810.59355  1951.92688 0.9661517
# seed30-LD9     2949.66667  1568.40645  4330.92688 0.0000000
# seed40-LD9    -2908.91667 -4290.17688 -1527.65645 0.0000000
# SD5-LDM       -7969.66667 -9350.92688 -6588.40645 0.0000000
# SD7-LDM       -3070.75000 -4452.01022 -1689.48978 0.0000000
# SDM-LDM       -4420.58333 -5801.84355 -3039.32312 0.0000000
# SDP9-LDM      -3544.33333 -4925.59355 -2163.07312 0.0000000
# SDPF-LDM      -1178.50000 -2559.76022   202.76022 0.1759587
# seed10-LDM    -2313.75000 -3695.01022  -932.48978 0.0000089
# seed20-LDM    -4777.50000 -6158.76022 -3396.23978 0.0000000
# seed30-LDM    -2398.50000 -3779.76022 -1017.23978 0.0000034
# seed40-LDM    -8257.08333 -9638.34355 -6875.82312 0.0000000
# SD7-SD5        4898.91667  3517.65645  6280.17688 0.0000000
# SDM-SD5        3549.08333  2167.82312  4930.34355 0.0000000
# SDP9-SD5       4425.33333  3044.07312  5806.59355 0.0000000
# SDPF-SD5       6791.16667  5409.90645  8172.42688 0.0000000
# seed10-SD5     5655.91667  4274.65645  7037.17688 0.0000000
# seed20-SD5     3192.16667  1810.90645  4573.42688 0.0000000
# seed30-SD5     5571.16667  4189.90645  6952.42688 0.0000000
# seed40-SD5     -287.41667 -1668.67688  1093.84355 0.9999212
# SDM-SD7       -1349.83333 -2731.09355    31.42688 0.0619396
# SDP9-SD7       -473.58333 -1854.84355   907.67688 0.9920955
# SDPF-SD7       1892.25000   510.98978  3273.51022 0.0006995
# seed10-SD7      757.00000  -624.26022  2138.26022 0.8014051
# seed20-SD7    -1706.75000 -3088.01022  -325.48978 0.0038410
# seed30-SD7      672.25000  -709.01022  2053.51022 0.8984318
# seed40-SD7    -5186.33333 -6567.59355 -3805.07312 0.0000000
# SDP9-SDM        876.25000  -505.01022  2257.51022 0.6151423
# SDPF-SDM       3242.08333  1860.82312  4623.34355 0.0000000
# seed10-SDM     2106.83333   725.57312  3488.09355 0.0000818
# seed20-SDM     -356.91667 -1738.17688  1024.34355 0.9993590
# seed30-SDM     2022.08333   640.82312  3403.34355 0.0001949
# seed40-SDM    -3836.50000 -5217.76022 -2455.23978 0.0000000
# SDPF-SDP9      2365.83333   984.57312  3747.09355 0.0000050
# seed10-SDP9    1230.58333  -150.67688  2611.84355 0.1311048
# seed20-SDP9   -1233.16667 -2614.42688   148.09355 0.1291348
# seed30-SDP9    1145.83333  -235.42688  2527.09355 0.2093103
# seed40-SDP9   -4712.75000 -6094.01022 -3331.48978 0.0000000
# seed10-SDPF   -1135.25000 -2516.51022   246.01022 0.2210043
# seed20-SDPF   -3599.00000 -4980.26022 -2217.73978 0.0000000
# seed30-SDPF   -1220.00000 -2601.26022   161.26022 0.1394203
# seed40-SDPF   -7078.58333 -8459.84355 -5697.32312 0.0000000
# seed20-seed10 -2463.75000 -3845.01022 -1082.48978 0.0000016
# seed30-seed10   -84.75000 -1466.01022  1296.51022 1.0000000
# seed40-seed10 -5943.33333 -7324.59355 -4562.07312 0.0000000
# seed30-seed20  2379.00000   997.73978  3760.26022 0.0000043
# seed40-seed20 -3479.58333 -4860.84355 -2098.32312 0.0000000
# seed40-seed30 -5858.58333 -7239.84355 -4477.32312 0.0000000

TukeyHSD(x=a1, 'type', conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $type
#              diff       lwr       upr     p adj
# obs-exp -162.0694 -497.1894 173.0505 0.3404128

# plots confirmed the anova and tukeyHSD results
library(ggplot2)
DEsummary<-read.table("s3.DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEsummary, aes(x=type, y=sig, fill=DE_method)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom")
# significant DEs: mapping - hylite, polycat and rsem are similiar
ggplot(data=DEsummary, aes(x=type, y=sig, fill=mapping)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom")
# exp vs obs: fewer significant DEs in observed
DEauc<-read.table("DE/DE.ROC.txt",header=TRUE,sep="\t")
DEauc[,4] <- as.numeric(DEauc[,4])
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEauc, aes(x=mapping, y=AUC, fill=DE_method)) + geom_boxplot() + theme(legend.position="bottom")



#######################################
##### sensitivity and specificity #####
#######################################
# compare oobs and exp
source("addTrans.r")
eval<-c("compr","mapping","DE_method","exp.sig","obs.sig","both.sig","sensitivity","specificity")
pdf("s3.DE.evals.pdf")
par(mfrow=c(3,2))
for(compr in paste0(batch[,1],"vs",batch[,2]) )
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
write.table(DEeval,"s3.DE.evals.txt",row.names=FALSE,sep="\t")

fit<- lm(sensitivity~compr+mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
# Df  Sum Sq Mean Sq F value   Pr(>F)
# compr       11 0.20774 0.018885   26.02  < 2e-16 ***
# mapping      2 0.02955 0.014775   20.36 2.13e-07 ***
# DE_method    1 0.00990 0.009903   13.65 0.000496 ***
# Residuals   57 0.04136 0.000726
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
TukeyHSD(x=a1, 'DE_method', conf.level=0.95)  # ebseq > deseq, sig
TukeyHSD(x=a1, 'mapping', conf.level=0.95)    # polycat > hylite = rsem

fit<- lm(specificity~compr+mapping+DE_method,data=DEeval)
summary(a1<-aov(fit))
# Df  Sum Sq Mean Sq F value   Pr(>F)
# compr       11 0.0538  0.0049    3.729  0.00049 ***
# mapping      2 0.0313  0.0157   11.958 4.61e-05 ***
# DE_method    1 1.8412  1.8412 1404.948  < 2e-16 ***
# Residuals   57 0.0747  0.0013
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
TukeyHSD(x=a1, 'DE_method', conf.level=0.95)  # ebseq > deseq, sig
TukeyHSD(x=a1, 'mapping', conf.level=0.95)    # polycat > hylite = rsem


###########################
### ROC curves: p-value ###
###########################

library(ROCR)

resultAUC <- c("mapping","DE_method","compr","AUC")
pdf("s3.DE.ROC.pdf")
for(mapping in c("polycat","rsem","hylite"))
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

# What is affecting DE ROC
# DEresultAUC<-read.table("s3.DE.ROC.txt",sep="\t",header=TRUE)
fit<- lm(AUC~compr+mapping+DE_method,data=DEresultAUC)
summary(a1<-aov(fit))
# Df  Sum Sq Mean Sq F value   Pr(>F)
# compr       11 0.07076 0.00643    5.76 3.65e-06 ***
# mapping      2 0.04674 0.02337   20.93 1.53e-07 ***
# DE_method    1 0.15079 0.15079  135.02  < 2e-16 ***
# Residuals   57 0.06366 0.00112
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
TukeyHSD(x=a1, 'DE_method', conf.level=0.95)  # ebseq > deseq, sig
TukeyHSD(x=a1, 'mapping', conf.level=0.95)    # polycat > hylite > rsem


###########################
###### Summary Plots ######
###########################

libray(ggplot2)
pdf("s3.DE.performance.pdf")

## plot DE numbers relative to DE_method or mapping, with obs and exp in x.axis and comparisons in facet.
DEsummary<-read.table("s3.DE.summary.txt",header=TRUE,sep="\t")
DEsummary[,-c(1:4)] <- apply(DEsummary[,-c(1:4)],2,as.numeric)
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEsummary, aes(x=type, y=sig, fill=DE_method)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom")
# significant DEs: mapping - hylite, polycat and rsem are similiar
ggplot(data=DEsummary, aes(x=type, y=sig, fill=mapping)) + geom_boxplot() + facet_grid(.~compr) + theme(legend.position="bottom")


## plot sensitivity and specificity
DEeval<-read.table("s3.DE.evals.txt",header=TRUE,sep="\t")
DEeval[,-c(1:3)] <- apply(DEeval[,-c(1:3)],2,as.numeric)
ggplot(data=DEeval, aes(x=mapping, y=sensitivity, fill=DE_method)) + geom_boxplot()
ggplot(data=DEeval, aes(x=DE_method, y=sensitivity, fill=mapping)) + geom_boxplot()
ggplot(data=DEeval, aes(x=mapping, y=specificity, fill=DE_method)) + geom_boxplot()
ggplot(data=DEeval, aes(x=DE_method, y=specificity, fill=mapping)) + geom_boxplot()
# low sensitivity from Hylite, specificity are similiar
# EBseq shows higher specificity than DEseq, probably due to higher stringency indicated by fewer DEs.


## plot AUC
DEauc<-read.table("s3.DE.ROC.txt",header=TRUE,sep="\t")
DEauc[,4] <- as.numeric(DEauc[,4])
# significant DEs: DE_method - deseq2 detects more DEs than ebseq
ggplot(data=DEauc, aes(x=mapping, y=AUC, fill=DE_method)) + geom_boxplot()
ggplot(data=DEauc, aes(x=DE_method, y=AUC, fill=mapping)) + geom_boxplot()

dev.off()

## plot DE numbers by comp
pdf("s3.DE.summary.pdf")
summary(DEsummary$sig/372.23) # use 10-50 range
for(compr in unique(DEsummary$compr) )
{
    p<-ggplot(data=DEsummary[ DEsummary$compr==compr,], aes(x=type, y=sig/372.23, fill=mapping)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() + scale_y_continuous(limits=c(10, 50))+ theme(panel.border=element_blank(),axis.line.x=element_line()) + scale_fill_manual(values=c( "#b4b2bc", "black","#22aad6")) + ggtitle(as.character(compr))
    print(p)
}
dev.off()


            
