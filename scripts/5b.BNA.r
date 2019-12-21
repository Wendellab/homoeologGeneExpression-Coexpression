# convert correlation matrix to list
cor2list<-function(cor_mat)
{
    library(reshape2)
    cor_mat[lower.tri(cor_mat,diag=TRUE)]<-NA
    cor_list <- subset(melt(cor_mat, na.rm=TRUE))
    return(cor_list)
}
mat <- structure(c(0, 0, 0, 0, 0.85, 0, 0, 0, 0.45, 0.85, 0, 0, 0.96, 0.56, 0.45, 0), .Dim = c(4L, 4L), .Dimnames = list(c("A1", "B1", "C1", "D1"), c("A1", "B1", "C1", "D1")))
mat
cor2list(mat)

#############################
## Binary Network Analysis ##
#############################

options(stringsAsFactors=FALSE)
library(psych)
library(ROCR)
library(WGCNA) # cor is faster

# Binary network generated multiple ways:
# 1. thresholded to the top correlation coefficients of gene connections 0.5%
# 2. thresholded by a set of fisher transformed Zscore cutoffs (1.5, 2.0, 2.5, 3.0)

# For each network construction method, use ROC and AUC to compare the performances of homoelog expression eistimates:
# Due to the large gene number of polyploid networks (> 60,000 genes), run permutation test by sampling about 10% of genes (nPermute = 10, sample size = 6,000)
nPermutation = 10
# sampleSize =6000

# for a fair comparison, i need to keep the same set of genes in all network
rdatafiles<-grep("R-05-dataInput",list.files(),value=TRUE)
rdatafiles
load(rdatafiles[1])
use=colnames(multiExpr[[1]]$data)
print(length(use))
for(i in rdatafiles[-1])
{
    load(i);print(i);
    ids =colnames(multiExpr[[1]]$data)
    print(length(ids))
    use=intersect(use, ids)
}
print(length(use)) #62656


ptm<-proc.time()
resultAUC <- c("permutation","dataset","method","category","AUC")
par(mfrow=c(2,2))
for(file in rdatafiles)
{
    # file = "R-05-dataInput.salmon_rld.RData"
    message(file)
    start<-proc.time()
    flag <- gsub("R-05-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    
    # get expression datasets of exp and obs
    exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data[,use]))
    obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data[,use]))
    
    # test small sets, 10 permutation
    nGenes<-nrow(exp)
    sampleSize <- round(nGenes/10)
    for(permute in 1:nPermutation)
    {
        message(paste0(".permutation ",permute))
        s1<-proc.time()
        select<- sample(1:nGenes,sampleSize) # 10% of ~70,000
        e<-exp[select,]
        o<-obs[select,]
        
        # get correlation matrix
        system.time( exp_cor <- cor(t(e), use = 'pairwise.complete.obs', nThread=8) )
        system.time( obs_cor <- cor(t(o), use = 'pairwise.complete.obs', nThread=8) )
        
        # get correlation and ranks
        #tri<-upper.tri(exp_cor, diag=FALSE)
        #exp_p <- exp_cor[tri]
        #obs_p <- obs_cor[tri]
        #np<-length(exp_p)
        # percentage rank correlation
        #exp_rank<-rank(exp_p)/np
        #obs_rank<-rank(obs_p)/np
        
        # convert matrix to list
        exp_corl<-cor2list(exp_cor)
        obs_corl<-cor2list(obs_cor)
        
        # percentage rank correlation
        exp_corl$rank<-rank(exp_corl$value)/nrow(exp_corl)
        obs_corl$rank<-rank(obs_corl$value)/nrow(obs_corl)
        
        # categorize edges: aa, dd, ad
        group <- as.factor(paste0(gsub(".*00","",exp_corl$Var1),gsub(".*00","",exp_corl$Var2)))
        group[group=="da"]<-"ad"
        group <- factor(group)
        
        
        # ------------------ Method 1 -------------------------------
        message("...top percentage")
        # top percentage
        Ps <- c(0.01, 0.05, 0.005, 0.001)
        # corresponding ranking cutoffs, consider above cutoff
        Rs<- (1-Ps)
        
        # define true values and predictions
        trues<-list()
        predictions<-list()
        prediction <- obs_corl$rank
        for(i in 1:length(Rs)) {
            cut <- Rs[i]
            predictions[[i]] <- prediction
            trues[[i]] <- (exp_corl$rank > cut)
        }
        # ROC analysis
        pred <- prediction(predictions = predictions,labels= trues)
        #roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
        auc.perf = performance(pred, measure = "auc")
        aucs1 <- round(as.numeric(auc.perf@y.values ),3)
        resultAUC <- rbind(resultAUC, cbind(rep(permute,length(Rs)), rep(flag,length(Rs)), paste0("Top",Ps), rep("-",length(Rs)),c(aucs1) ) )
        
        # ROC loop for edge categories
        for(c in levels(group))
        {
            sub_predictions<-lapply(predictions,function(x)x[group==c])
            sub_trues<-lapply(trues,function(x)x[group==c])
            pred <- prediction(predictions = sub_predictions,labels= sub_trues)
            #roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
            auc.perf = performance(pred, measure = "auc")
            aucs1 <- round(as.numeric(auc.perf@y.values ),3)
            resultAUC <- rbind(resultAUC, cbind(rep(permute,length(Rs)), rep(flag,length(Rs)), paste0("Top",Ps), rep(c,length(Rs)),c(aucs1) ) )
        }
        
        # -----------------------------------------------------------
        
        # ------------------ Method 2 -------------------------------
        message("...fisher's Z")
        # chosen Zs cutoffs
        Zs <-c(1.5, 2.0, 2.5, 3.0)
        # corresponding r cutoffs are
        Rs<-fisherz2r(Zs)   # 0.9051483 0.9640276 0.9866143 0.9950548
        
        # define true values and predictions
        trues<-list()
        predictions<-list()
        prediction <- obs_corl$value
        for(i in 1:length(Rs)) {
            cut <- Rs[i]
            predictions[[i]] <- prediction
            trues[[i]] <- (exp_corl$value > cut)
        }
        # in case trues has no TRUE
        filter<- which(lapply(trues,function(x)nlevels(factor(x)))==1)
        if(length(filter)>0)
        {predictions = predictions[-filter]; trues= trues[-filter]; Rs= Rs[-filter]; Zs=Zs[-filter]}
        # ROC analysis
        pred <- prediction(predictions = predictions,labels= trues)
        #roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
        auc.perf = performance(pred, measure = "auc")
        aucs2 <- round(as.numeric(auc.perf@y.values ),3)
        resultAUC <- rbind(resultAUC, cbind(rep(permute,length(Rs)), rep(flag,length(Rs)), paste0("Zcutoff",Zs), rep("-",length(Rs)),c(aucs2) ) )
        # ROC loop for edge categories
        for(c in levels(group))
        {
            sub_predictions<-lapply(predictions,function(x)x[group==c])
            sub_trues<-lapply(trues,function(x)x[group==c])
            check<-sapply(sub_trues,any) # see if any has no true values
            pred <- prediction(predictions = sub_predictions[check],labels= sub_trues[check])
            # roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
            auc.perf = performance(pred, measure = "auc")
            aucs2 <- rep(NA,length(Rs))
            aucs2[check] <- round(as.numeric(auc.perf@y.values ),3)
            resultAUC <- rbind(resultAUC, cbind(rep(permute,length(Rs)), rep(flag,length(Rs)), paste0("Zcutoff",Zs), rep(c,length(Rs)),c(aucs2) ) )
        }
        
        message(paste0("...running permutation ",permute, " for ", (proc.time()-s1)[3]/60," minutes."))
    }
    # -----------------------------------------------------------
    ptm<-proc.time()-start
    message(paste0("...running time for ",flag," used ", round(ptm[3]/60,1), " minutes."))
    save(resultAUC, file = paste0("s5.bna.AUC.Rdata"))
    
}

#############################
## Plots and post-analysis ##
#############################
load("s5.bna.AUC.Rdata")
sumAUC <- data.frame(resultAUC[-1,])
names(sumAUC) <- resultAUC[1,]
sumAUC$AUC<-as.numeric(as.character(sumAUC$AUC))
head(sumAUC)
# permutation         dataset   method category   AUC
#           1 hylite_log2rpkm  Top0.01        - 0.997
#           1 hylite_log2rpkm  Top0.05        - 0.994
#           1 hylite_log2rpkm Top0.005        - 0.997
#           1 hylite_log2rpkm Top0.001        - 0.998
#           1 hylite_log2rpkm  Top0.01       aa 0.994
#           1 hylite_log2rpkm  Top0.05       aa 0.992
xtabs(~dataset+method+ category, sumAUC)
sumAUC$mapping=gsub("_.*","",sumAUC$dataset)
sumAUC$normalization=gsub(".*_","",sumAUC$dataset)
sumAUC$cutoffs =gsub("Top|Zcutoff","",sumAUC$method)
sumAUC$threshold =gsub('[0-9]|[.]',"",sumAUC$method)
head(sumAUC)
# permutation         dataset   method category   AUC mapping normalization cutoffs threshold
#           1 hylite_log2rpkm  Top0.01        - 0.997  hylite      log2rpkm    0.01       Top
#           1 hylite_log2rpkm  Top0.05        - 0.994  hylite      log2rpkm    0.05       Top
#           1 hylite_log2rpkm Top0.005        - 0.997  hylite      log2rpkm   0.005       Top
#           1 hylite_log2rpkm Top0.001        - 0.998  hylite      log2rpkm   0.001       Top
#           1 hylite_log2rpkm  Top0.01       aa 0.994  hylite      log2rpkm    0.01       Top
#           1 hylite_log2rpkm  Top0.05       aa 0.992  hylite      log2rpkm    0.05       Top

sumAUC$method = factor(sumAUC$method, levels= c("Top0.05", "Top0.01","Top0.005","Top0.001","Zcutoff1.5", "Zcutoff2","Zcutoff2.5","Zcutoff3"))
summary(sumAUC)

# anova and TukeyHSD
fit<-lm(AUC~mapping+normalization+method,data=sumAUC[sumAUC$category=="-",] )
summary(fit)
# Call:
# lm(formula = AUC ~ mapping + normalization + method, data = sumAUC[sumAUC$category =="-", ])
#
# Residuals:
# Min        1Q    Median        3Q       Max
# -0.170409 -0.002325 -0.000107  0.005243  0.026696
#
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)       0.9971193  0.0017724 562.596  < 2e-16 ***
# mappingkallisto  -0.0032508  0.0015579  -2.087   0.0372 *
# mappingpolycat   -0.0016747  0.0015655  -1.070   0.2851
# mappingrsem      -0.0094672  0.0015556  -6.086 1.82e-09 ***
# mappingsalmon    -0.0022173  0.0015629  -1.419   0.1564
# normalizationrld  0.0001054  0.0009876   0.107   0.9150
# methodTop0.01     0.0028500  0.0019613   1.453   0.1466
# methodTop0.005    0.0034600  0.0019613   1.764   0.0781 .
# methodTop0.001    0.0040900  0.0019613   2.085   0.0374 *
# methodZcutoff1.5  0.0041000  0.0019613   2.090   0.0369 *
# methodZcutoff2    0.0031900  0.0019613   1.626   0.1043
# methodZcutoff2.5 -0.0009800  0.0019613  -0.500   0.6175
# methodZcutoff3   -0.0143484  0.0020162  -7.116 2.52e-12 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.01387 on 777 degrees of freedom
# Multiple R-squared:  0.1841,    Adjusted R-squared:  0.1715
# F-statistic: 14.62 on 12 and 777 DF,  p-value: < 2.2e-16
summary(a<-aov(fit))
# Df  Sum Sq  Mean Sq F value   Pr(>F)
# mapping         4 0.00891 0.002227   11.58 3.91e-09 ***
# normalization   1 0.00005 0.000046    0.24    0.624
# method          7 0.02478 0.003540   18.41  < 2e-16 ***
# Residuals     777 0.14945 0.000192
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 TukeyHSD(x=a, conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $mapping
#                           diff          lwr           upr     p adj
# kallisto-hylite  -0.0033437226 -0.007603345  0.0009158996 0.2016827
# polycat-hylite   -0.0014852321 -0.005765220  0.0027947558 0.8775894
# rsem-hylite      -0.0096518987 -0.013904881 -0.0053989163 0.0000000
# salmon-hylite    -0.0021232363 -0.006396360  0.0021498874 0.6544337
# polycat-kallisto  0.0018584906 -0.002414805  0.0061317865 0.7576818
# rsem-kallisto    -0.0063081761 -0.010554424 -0.0020619282 0.0005131
# salmon-kallisto   0.0012204863 -0.003045935  0.0054869074 0.9357717
# rsem-polycat     -0.0081666667 -0.012433344 -0.0038999892 0.0000021
# salmon-polycat   -0.0006380042 -0.004924759  0.0036487502 0.9942340
# salmon-rsem       0.0075286624  0.003268871  0.0117884542 0.0000160
#
# $normalization
#                     diff          lwr         upr     p adj
# rld-log2rpkm 0.000483279 -0.001454105 0.002420663 0.6245025
#
# $method
#                              diff         lwr           upr     p adj
# Top0.005-Top0.001     -0.00063000 -0.00659091  0.0053309097 0.9999830
# Top0.01-Top0.001      -0.00124000 -0.00720091  0.0047209097 0.9984278
# Top0.05-Top0.001      -0.00409000 -0.01005091  0.0018709097 0.4253652
# Zcutoff1.5-Top0.001    0.00001000 -0.00595091  0.0059709097 1.0000000
# Zcutoff2-Top0.001     -0.00090000 -0.00686091  0.0050609097 0.9998088
# Zcutoff2.5-Top0.001   -0.00507000 -0.01103091  0.0008909097 0.1629136
# Zcutoff3-Top0.001     -0.01840682 -0.02453108 -0.0122825708 0.0000000
# Top0.01-Top0.005      -0.00061000 -0.00657091  0.0053509097 0.9999864
# Top0.05-Top0.005      -0.00346000 -0.00942091  0.0025009097 0.6446957
# Zcutoff1.5-Top0.005    0.00064000 -0.00532091  0.0066009097 0.9999810
# Zcutoff2-Top0.005     -0.00027000 -0.00623091  0.0056909097 1.0000000
# Zcutoff2.5-Top0.005   -0.00444000 -0.01040091  0.0015209097 0.3151885
# Zcutoff3-Top0.005     -0.01777682 -0.02390108 -0.0116525708 0.0000000
# Top0.05-Top0.01       -0.00285000 -0.00881091  0.0031109097 0.8318791
# Zcutoff1.5-Top0.01     0.00125000 -0.00471091  0.0072109097 0.9983446
# Zcutoff2-Top0.01       0.00034000 -0.00562091  0.0063009097 0.9999998
# Zcutoff2.5-Top0.01    -0.00383000 -0.00979091  0.0021309097 0.5148536
# Zcutoff3-Top0.01      -0.01716682 -0.02329108 -0.0110425708 0.0000000
# Zcutoff1.5-Top0.05     0.00410000 -0.00186091  0.0100609097 0.4220237
# Zcutoff2-Top0.05       0.00319000 -0.00277091  0.0091509097 0.7342129
# Zcutoff2.5-Top0.05    -0.00098000 -0.00694091  0.0049809097 0.9996627
# Zcutoff3-Top0.05      -0.01431682 -0.02044108 -0.0081925708 0.0000000
# Zcutoff2-Zcutoff1.5   -0.00091000 -0.00687091  0.0050509097 0.9997942
# Zcutoff2.5-Zcutoff1.5 -0.00508000 -0.01104091  0.0008809097 0.1610210
# Zcutoff3-Zcutoff1.5   -0.01841682 -0.02454108 -0.0122925708 0.0000000
# Zcutoff2.5-Zcutoff2   -0.00417000 -0.01013091  0.0017909097 0.3989089
# Zcutoff3-Zcutoff2     -0.01750682 -0.02363108 -0.0113825708 0.0000000
# Zcutoff3-Zcutoff2.5   -0.01333682 -0.01946108 -0.0072125708 0.0000000
library(agricolae)
out <- duncan.test(a,trt="mapping"); out
# $statistics
#       MSerror  Df      Mean       CV
#   0.000192344 777 0.9943089 1.394819
#
# $parameters
#     test  name.t ntr alpha
#   Duncan mapping   5  0.05
#
# $duncan
# NULL
#
# $means
#                AUC         std   r   Min Max   Q25   Q50    Q75
# hylite   0.9976519 0.003307831 158 0.965   1 0.997 0.998 0.9990
# kallisto 0.9943082 0.012253106 159 0.916   1 0.995 0.997 0.9985
# polycat  0.9961667 0.003625774 156 0.975   1 0.995 0.997 0.9990
# rsem     0.9880000 0.028789848 160 0.803   1 0.995 0.997 0.9980
# salmon   0.9955287 0.009810095 157 0.911   1 0.995 0.998 0.9990
#
# $comparison
# NULL
#
# $groups
#                AUC groups
# hylite   0.9976519      a
# polycat  0.9961667     ab
# salmon   0.9955287     ab
# kallisto 0.9943082      b
# rsem     0.9880000      c
#
# attr(,"class")
# [1] "group"
out <- duncan.test(a,trt="method"); out
# $statistics
#       MSerror  Df      Mean       CV
#   0.000192344 777 0.9943089 1.394819
#
# $parameters
#     test name.t ntr alpha
#   Duncan method   8  0.05
#
# $duncan
# NULL
#
# $means
#                  AUC          std   r   Min   Max     Q25    Q50   Q75
# Top0.001   0.9979400 0.0011531249 100 0.995 1.000 0.99700 0.9980 0.999
# Top0.005   0.9973100 0.0009287224 100 0.995 0.999 0.99700 0.9970 0.998
# Top0.01    0.9967000 0.0010777830 100 0.994 0.999 0.99600 0.9970 0.997
# Top0.05    0.9938500 0.0019036117 100 0.989 0.997 0.99300 0.9950 0.995
# Zcutoff1.5 0.9979500 0.0011666667 100 0.995 1.000 0.99700 0.9980 0.999
# Zcutoff2   0.9970400 0.0033179037 100 0.985 1.000 0.99500 0.9980 0.999
# Zcutoff2.5 0.9928700 0.0127688502 100 0.915 1.000 0.99000 0.9985 1.000
# Zcutoff3   0.9793111 0.0396321979  90 0.803 1.000 0.97725 1.0000 1.000
#
# $comparison
# NULL
#
# $groups
#                  AUC groups
# Zcutoff1.5 0.9979500      a
# Top0.001   0.9979400      a
# Top0.005   0.9973100      a
# Zcutoff2   0.9970400     ab
# Top0.01    0.9967000     ab
# Top0.05    0.9938500     ab
# Zcutoff2.5 0.9928700      b
# Zcutoff3   0.9793111      c
#
# attr(,"class")
# [1] "group"



# why doesn hylite shows higher AUC than other pipelines? Bigger libSize?
# also hylite networks contain fewest genes
for(i in rdatafiles)
{load(i);print(i);print(rowSums(multiExpr[[1]]$data)[1])}

library(ggplot2)
pdf(paste0("s5.bna.plotAUC.sample", sampleSize,".permutation",nPermutation,".pdf"))
# one factor, AUC ~ mapping + normalization + method, data = sumAUC[sumAUC$category =="-", ])
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=mapping, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=normalization, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
# show only the significant: method and mapping
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC: all gene pairs")
#+ facet_grid(normalization~.)
# compare connections within network
ggplot(data=sumAUC[sumAUC$method=="Top0.05", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Top0.01", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Top0.005", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Zcutoff1.5", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Zcutoff2", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Zcutoff2.5", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
ggplot(data=sumAUC[sumAUC$method=="Zcutoff3", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)
dev.off()

fit<-lm(AUC~mapping+category+method,data=sumAUC)
summary(fit)
out <- duncan.test(a,trt="category"); out
TukeyHSD(x=a, "category", conf.level=0.95)
# dd> ad/- > aa
