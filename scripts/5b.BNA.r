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

# anova and TukeyHSD
summary(sumAUC)
fit<-lm(AUC~mapping+normalization+method+category,data=sumAUC )
summary(fit)
#
# Call:
# lm(formula = AUC ~ mapping + normalization + method + category,
#     data = sumAUC)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.18127 -0.01089  0.00399  0.01645  0.07763
#
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)
# (Intercept)       1.0047417  0.0020690 485.610  < 2e-16 ***
# mappingkallisto  -0.0279788  0.0016468 -16.990  < 2e-16 ***
# mappingpolycat   -0.0053881  0.0016468  -3.272  0.00108 **
# mappingrsem      -0.0266522  0.0016468 -16.185  < 2e-16 ***
# mappingsalmon    -0.0258631  0.0016468 -15.705  < 2e-16 ***
# normalizationrld  0.0145381  0.0010365  14.026  < 2e-16 ***
# methodTop0.005    0.0009750  0.0020664   0.472  0.63708
# methodTop0.01     0.0004600  0.0020664   0.223  0.82386
# methodTop0.05    -0.0026700  0.0020664  -1.292  0.19643
# methodZcutoff1.5  0.0003050  0.0020664   0.148  0.88267
# methodZcutoff2   -0.0162700  0.0020664  -7.873  4.7e-15 ***
# methodZcutoff2.5 -0.0513475  0.0020664 -24.848  < 2e-16 ***
# methodZcutoff3   -0.0637215  0.0020929 -30.447  < 2e-16 ***
# categoryaa       -0.0132642  0.0014645  -9.057  < 2e-16 ***
# categoryad       -0.0008953  0.0014654  -0.611  0.54127
# categorydd        0.0132722  0.0014635   9.069  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.02922 on 3165 degrees of freedom
#   (15 observations deleted due to missingness)
# Multiple R-squared:  0.5096,	Adjusted R-squared:  0.5072
# F-statistic: 219.2 on 15 and 3165 DF,  p-value: < 2.2e-16

summary(a<-aov(fit))
#                 Df Sum Sq Mean Sq F value Pr(>F)
# mapping          4 0.4811 0.12028   140.8 <2e-16 ***
# normalization    1 0.1804 0.18035   211.2 <2e-16 ***
# method           7 1.8666 0.26666   312.2 <2e-16 ***
# category         3 0.2804 0.09347   109.4 <2e-16 ***
# Residuals     3165 2.7030 0.00085
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 15 observations deleted due to missingness

TukeyHSD(x=a, conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $mapping
#                           diff          lwr          upr     p adj
# kallisto-hylite  -0.0292709415 -0.033763762 -0.024778121 0.0000000
# polycat-hylite   -0.0066803165 -0.011173137 -0.002187496 0.0004858
# rsem-hylite      -0.0279443790 -0.032437199 -0.023451559 0.0000000
# salmon-hylite    -0.0271553165 -0.031648137 -0.022662496 0.0000000
# polycat-kallisto  0.0225906250  0.018131781  0.027049469 0.0000000
# rsem-kallisto     0.0013265625 -0.003132282  0.005785407 0.9270263
# salmon-kallisto   0.0021156250 -0.002343219  0.006574469 0.6943407
# rsem-polycat     -0.0212640625 -0.025722907 -0.016805218 0.0000000
# salmon-polycat   -0.0204750000 -0.024933844 -0.016016156 0.0000000
# salmon-rsem       0.0007890625 -0.003669782  0.005247907 0.9889256
#
# $normalization
#                    diff        lwr       upr p adj
# rld-log2rpkm 0.01505878 0.01302686 0.0170907     0
#
# $method
#                              diff          lwr          upr     p adj
# Top0.005-Top0.001      0.00097500 -0.005292349  0.007242349 0.9997714
# Top0.01-Top0.001       0.00046000 -0.005807349  0.006727349 0.9999987
# Top0.05-Top0.001      -0.00267000 -0.008937349  0.003597349 0.9020538
# Zcutoff1.5-Top0.001    0.00030500 -0.005962349  0.006572349 0.9999999
# Zcutoff2-Top0.001     -0.01627000 -0.022537349 -0.010002651 0.0000000
# Zcutoff2.5-Top0.001   -0.05134750 -0.057614849 -0.045080151 0.0000000
# Zcutoff3-Top0.001     -0.06357276 -0.069917763 -0.057227753 0.0000000
# Top0.01-Top0.005      -0.00051500 -0.006782349  0.005752349 0.9999971
# Top0.05-Top0.005      -0.00364500 -0.009912349  0.002622349 0.6447438
# Zcutoff1.5-Top0.005   -0.00067000 -0.006937349  0.005597349 0.9999820
# Zcutoff2-Top0.005     -0.01724500 -0.023512349 -0.010977651 0.0000000
# Zcutoff2.5-Top0.005   -0.05232250 -0.058589849 -0.046055151 0.0000000
# Zcutoff3-Top0.005     -0.06454776 -0.070892763 -0.058202753 0.0000000
# Top0.05-Top0.01       -0.00313000 -0.009397349  0.003137349 0.7997689
# Zcutoff1.5-Top0.01    -0.00015500 -0.006422349  0.006112349 1.0000000
# Zcutoff2-Top0.01      -0.01673000 -0.022997349 -0.010462651 0.0000000
# Zcutoff2.5-Top0.01    -0.05180750 -0.058074849 -0.045540151 0.0000000
# Zcutoff3-Top0.01      -0.06403276 -0.070377763 -0.057687753 0.0000000
# Zcutoff1.5-Top0.05     0.00297500 -0.003292349  0.009242349 0.8387588
# Zcutoff2-Top0.05      -0.01360000 -0.019867349 -0.007332651 0.0000000
# Zcutoff2.5-Top0.05    -0.04867750 -0.054944849 -0.042410151 0.0000000
# Zcutoff3-Top0.05      -0.06090276 -0.067247763 -0.054557753 0.0000000
# Zcutoff2-Zcutoff1.5   -0.01657500 -0.022842349 -0.010307651 0.0000000
# Zcutoff2.5-Zcutoff1.5 -0.05165250 -0.057919849 -0.045385151 0.0000000
# Zcutoff3-Zcutoff1.5   -0.06387776 -0.070222763 -0.057532753 0.0000000
# Zcutoff2.5-Zcutoff2   -0.03507750 -0.041344849 -0.028810151 0.0000000
# Zcutoff3-Zcutoff2     -0.04730276 -0.053647763 -0.040957753 0.0000000
# Zcutoff3-Zcutoff2.5   -0.01222526 -0.018570263 -0.005880253 0.0000002
#
# $category
#                diff          lwr          upr     p adj
# aa--  -0.0132634700 -0.017027584 -0.009499356 0.0000000
# ad--  -0.0008943375 -0.004660834  0.002872159 0.9288308
# dd--   0.0132726228  0.009510882  0.017034364 0.0000000
# ad-aa  0.0123691326  0.008596737  0.016141528 0.0000000
# dd-aa  0.0265360928  0.022768445  0.030303740 0.0000000
# dd-ad  0.0141669603  0.010396932  0.017936989 0.0000000


# why doesn hylite shows higher AUC than other pipelines? Bigger libSize?
# also hylite networks contain fewest genes
for(i in rdatafiles)
{load(i);print(i);print(rowSums(multiExpr[[1]]$data)[1])}


library(ggplot2)
pdf(paste0("s5.bna.plotAUC.sample", sampleSize,".permutation",nPermutation,".pdf"))
# single factor
ggplot(data=sumAUC, aes(x=permutation, y=AUC)) + geom_boxplot() + ggtitle("Boxplot of AUC: all gene pairs")
ggplot(data=sumAUC, aes(x=dataset, y=AUC)) + geom_boxplot() + ggtitle("Boxplot of AUC: all gene pairs")
ggplot(data=sumAUC, aes(x=method, y=AUC)) + geom_boxplot() + ggtitle("Boxplot of AUC: all gene pairs")
ggplot(data=sumAUC, aes(x=category, y=AUC)) + geom_boxplot() + ggtitle("Boxplot of AUC: all gene pairs")
# two factors
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC, fill=dataset )) + geom_boxplot() + ggtitle("Boxplot of AUC: all gene pairs")
ggplot(data=sumAUC[sumAUC$category=="aa",], aes(x=method, y=AUC, fill=dataset )) + geom_boxplot() + ggtitle("Boxplot of AUC:  gene pairs in A subgenome")
ggplot(data=sumAUC[sumAUC$category=="dd",], aes(x=method, y=AUC, fill=dataset )) + geom_boxplot() + ggtitle("Boxplot of AUC:  gene pairs in D subgenome")
ggplot(data=sumAUC[sumAUC$category=="ad",], aes(x=method, y=AUC, fill=dataset )) + geom_boxplot() + ggtitle("Boxplot of AUC:  gene pairs between A&D subgenomes")
#
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC, fill=mapping)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() + scale_y_continuous(limits=c(10, 50))+ theme(panel.border=element_blank(),axis.line.x=element_line())
dev.off()