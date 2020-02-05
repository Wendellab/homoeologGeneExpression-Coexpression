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
rdatafiles<-grep("eaglerc",list.files(pattern="R-05-dataInput"),value=TRUE, invert=T)
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
print(length(use)) # 62192
sampleSize=round(length(use)/10)

ptm<-proc.time()
resultAUC <- c("permutation","dataset","method","category","AUC")
par(mfrow=c(2,2))
rdatafiles<-list.files(pattern="R-05-dataInput")
for(file in rdatafiles)
{
    # file = "R-05-dataInput.salmon_rld.RData"
    message(file)
    start<-proc.time()
    flag <- gsub("R-05-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    
    # get expression datasets of exp and obs
    if(grep("eaglerc",flag)){
        exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
        obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs" )]]$data))
    }else{
        exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data[,use]))
        obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data[,use]))
    }

    
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
# 14 dataset x 8 netType x 10 permu x 4 cat
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
#
# Call:
# lm(formula = AUC ~ mapping + normalization + method, data = sumAUC[sumAUC$category ==
#     "-", ])
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.283080 -0.002513 -0.000316  0.004838  0.023037
#
# Coefficients:
#                    Estimate Std. Error t value             Pr(>|t|)
# (Intercept)       0.9948366  0.0017426 570.900 < 0.0000000000000002 ***
# mappingeaglerc   -0.0058151  0.0016856  -3.450             0.000582 ***
# mappinghylite     0.0021349  0.0016856   1.267             0.205592
# mappingkallisto  -0.0011737  0.0016881  -0.695             0.487042
# mappingpolycat    0.0007424  0.0016961   0.438             0.661702
# mappingrsem      -0.0057271  0.0016881  -3.393             0.000717 ***
# mappingsalmon     0.0014233  0.0016989   0.838             0.402335
# normalizationrld  0.0001169  0.0009030   0.129             0.896988
# methodTop0.01     0.0029357  0.0017932   1.637             0.101895
# methodTop0.005    0.0036786  0.0017932   2.051             0.040469 *
# methodTop0.001    0.0045357  0.0017932   2.529             0.011567 *
# methodZcutoff1.5  0.0045500  0.0017932   2.537             0.011309 *
# methodZcutoff2    0.0041071  0.0017932   2.290             0.022190 *
# methodZcutoff2.5 -0.0007500  0.0017932  -0.418             0.675853
# methodZcutoff3   -0.0121466  0.0018437  -6.588      0.0000000000691 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.015 on 1091 degrees of freedom
# Multiple R-squared:  0.1389,    Adjusted R-squared:  0.1279
# F-statistic: 12.57 on 14 and 1091 DF,  p-value: < 0.00000000000000022
#
summary(a<-aov(fit))
#                 Df  Sum Sq  Mean Sq F value               Pr(>F)
# mapping          6 0.01088 0.001813   8.055         0.0000000158 ***
# normalization    1 0.00006 0.000057   0.253                0.615
# method           7 0.02868 0.004097  18.203 < 0.0000000000000002 ***
# Residuals     1091 0.24558 0.000225
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(x=a, conf.level=0.95)
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
#
# Fit: aov(formula = fit)
#
# $mapping
#                           diff           lwr           upr     p adj
# eaglerc-bowtie   -0.0060626194 -0.0110410405 -0.0010841984 0.0061863
# hylite-bowtie     0.0018873806 -0.0030910405  0.0068658016 0.9221124
# kallisto-bowtie  -0.0013397028 -0.0063258714  0.0036464659 0.9856035
# polycat-bowtie    0.0008269639 -0.0041829683  0.0058368962 0.9990102
# rsem-bowtie      -0.0058931619 -0.0108793305 -0.0009069932 0.0090231
# salmon-bowtie     0.0015936306 -0.0034244015  0.0066116626 0.9663220
# hylite-eaglerc    0.0079500000  0.0029951922  0.0129048078 0.0000497
# kallisto-eaglerc  0.0047229167 -0.0002396756  0.0096855090 0.0742271
# polycat-eaglerc   0.0068895833  0.0019031151  0.0118760516 0.0009440
# rsem-eaglerc      0.0001694575 -0.0047931348  0.0051320498 0.9999999
# salmon-eaglerc    0.0076562500  0.0026616439  0.0126508561 0.0001342
# kallisto-hylite  -0.0032270833 -0.0081896756  0.0017355090 0.4668120
# polycat-hylite   -0.0010604167 -0.0060468849  0.0039260516 0.9959013
# rsem-hylite      -0.0077805425 -0.0127431348 -0.0028179502 0.0000831
# salmon-hylite    -0.0002937500 -0.0052883561  0.0047008561 0.9999977
# polycat-kallisto  0.0021666667 -0.0028275367  0.0071608701 0.8605008
# rsem-kallisto    -0.0045534591 -0.0095238237  0.0004169055 0.0976483
# salmon-kallisto   0.0029333333 -0.0020689953  0.0079356620 0.5944906
# rsem-polycat     -0.0067201258 -0.0117143292 -0.0017259224 0.0014587
# salmon-polycat    0.0007666667 -0.0042593492  0.0057926825 0.9993698
# salmon-rsem       0.0074867925  0.0024844638  0.0124891211 0.0002181
#
# $normalization
#                      diff          lwr         upr     p adj
# rld-log2rpkm 0.0004541999 -0.001316323 0.002224722 0.6148152
#
# $method
#                                 diff           lwr           upr     p adj
# Top0.01-Top0.05        0.00293571429 -0.0025099502  0.0083813788 0.7276102
# Top0.005-Top0.05       0.00367857143 -0.0017670931  0.0091242359 0.4475853
# Top0.001-Top0.05       0.00453571429 -0.0009099502  0.0099813788 0.1841365
# Zcutoff1.5-Top0.05     0.00455000000 -0.0008956645  0.0099956645 0.1809076
# Zcutoff2-Top0.05       0.00410714286 -0.0013385216  0.0095528073 0.2997558
# Zcutoff2.5-Top0.05    -0.00075000000 -0.0061956645  0.0046956645 0.9998978
# Zcutoff3-Top0.05      -0.01211164188 -0.0177065303 -0.0065167535 0.0000000
# Top0.005-Top0.01       0.00074285714 -0.0047028073  0.0061885216 0.9999042
# Top0.001-Top0.01       0.00160000000 -0.0038456645  0.0070456645 0.9867728
# Zcutoff1.5-Top0.01     0.00161428571 -0.0038313788  0.0070599502 0.9860625
# Zcutoff2-Top0.01       0.00117142857 -0.0042742359  0.0066170931 0.9980657
# Zcutoff2.5-Top0.01    -0.00368571429 -0.0091313788  0.0017599502 0.4449260
# Zcutoff3-Top0.01      -0.01504735617 -0.0206422446 -0.0094524678 0.0000000
# Top0.001-Top0.005      0.00085714286 -0.0045885216  0.0063028073 0.9997496
# Zcutoff1.5-Top0.005    0.00087142857 -0.0045742359  0.0063170931 0.9997204
# Zcutoff2-Top0.005      0.00042857143 -0.0050170931  0.0058742359 0.9999978
# Zcutoff2.5-Top0.005   -0.00442857143 -0.0098742359  0.0010170931 0.2096447
# Zcutoff3-Top0.005     -0.01579021331 -0.0213851017 -0.0101953249 0.0000000
# Zcutoff1.5-Top0.001    0.00001428571 -0.0054313788  0.0054599502 1.0000000
# Zcutoff2-Top0.001     -0.00042857143 -0.0058742359  0.0050170931 0.9999978
# Zcutoff2.5-Top0.001   -0.00528571429 -0.0107313788  0.0001599502 0.0644834
# Zcutoff3-Top0.001     -0.01664735617 -0.0222422446 -0.0110524678 0.0000000
# Zcutoff2-Zcutoff1.5   -0.00044285714 -0.0058885216  0.0050028073 0.9999972
# Zcutoff2.5-Zcutoff1.5 -0.00530000000 -0.0107456645  0.0001456645 0.0630615
# Zcutoff3-Zcutoff1.5   -0.01666164188 -0.0222565303 -0.0110667535 0.0000000
# Zcutoff2.5-Zcutoff2   -0.00485714286 -0.0103028073  0.0005885216 0.1209883
# Zcutoff3-Zcutoff2     -0.01621878474 -0.0218136731 -0.0106238963 0.0000000
# Zcutoff3-Zcutoff2.5   -0.01136164188 -0.0169565303 -0.0057667535 0.0000000
#
library(agricolae)
out <- duncan.test(a,trt="mapping"); out
# $statistics
#        MSerror   Df      Mean       CV
#   0.0002250962 1091 0.9947025 1.508311
#
# $parameters
#     test  name.t ntr alpha
#   Duncan mapping   7  0.05
#
# $duncan
# NULL
#
# $means
#                AUC         std   r   Min Max     Q25   Q50   Q75
# bowtie   0.9960064 0.004050161 157 0.985   1 0.99400 0.997 0.999
# eaglerc  0.9899438 0.016406107 160 0.915   1 0.99200 0.997 0.998
# hylite   0.9978938 0.003339485 160 0.965   1 0.99775 0.999 1.000
# kallisto 0.9946667 0.016532158 159 0.825   1 0.99600 0.998 0.999
# polycat  0.9968333 0.004460339 156 0.961   1 0.99600 0.998 0.999
# rsem     0.9901132 0.033745882 159 0.694   1 0.99600 0.997 0.998
# salmon   0.9976000 0.002542113 155 0.987   1 0.99700 0.998 0.999
#
# $comparison
# NULL
#
# $groups
#                AUC groups
# hylite   0.9978938      a
# salmon   0.9976000      a
# polycat  0.9968333      a
# bowtie   0.9960064      a
# kallisto 0.9946667      a
# rsem     0.9901132      b
# eaglerc  0.9899438      b
#
# attr(,"class")
# [1] "group"
out <- duncan.test(a,trt="normalization"); out
# $statistics
#        MSerror   Df      Mean       CV
#   0.0002250962 1091 0.9947025 1.508311
#
# $parameters
#     test        name.t ntr alpha
#   Duncan normalization   2  0.05
#
# $duncan
# NULL
#
# $means
#                AUC        std   r   Min Max   Q25   Q50   Q75
# log2rpkm 0.9944982 0.01311365 560 0.825   1 0.996 0.997 0.999
# rld      0.9949121 0.01862325 546 0.694   1 0.996 0.998 0.999
#
# $comparison
# NULL
#
# $groups
#                AUC groups
# rld      0.9949121      a
# log2rpkm 0.9944982      a
#
# attr(,"class")
# [1] "group"
out <- duncan.test(a,trt="method"); out
# $statistics
#        MSerror   Df      Mean       CV
#   0.0002250962 1091 0.9947025 1.508311
#
# $parameters
#     test name.t ntr alpha
#   Duncan method   8  0.05
#
# $duncan
# NULL
#
# $means
#                  AUC          std   r   Min   Max     Q25   Q50   Q75
# Top0.001   0.9982286 0.0009394845 140 0.995 1.000 0.99800 0.998 0.999
# Top0.005   0.9973714 0.0012485552 140 0.993 0.999 0.99700 0.998 0.998
# Top0.01    0.9966286 0.0015141674 140 0.992 0.999 0.99600 0.997 0.998
# Top0.05    0.9936929 0.0028281637 140 0.986 0.997 0.99300 0.995 0.996
# Zcutoff1.5 0.9982429 0.0009281484 140 0.996 1.000 0.99800 0.998 0.999
# Zcutoff2   0.9978000 0.0026881033 140 0.988 1.000 0.99675 0.999 1.000
# Zcutoff2.5 0.9929429 0.0122948029 140 0.936 1.000 0.99075 0.999 1.000
# Zcutoff3   0.9813810 0.0430812919 126 0.694 1.000 0.98475 1.000 1.000
#
# $comparison
# NULL
#
# $groups
#                  AUC groups
# Zcutoff1.5 0.9982429      a
# Top0.001   0.9982286      a
# Zcutoff2   0.9978000      a
# Top0.005   0.9973714     ab
# Top0.01    0.9966286    abc
# Top0.05    0.9936929     bc
# Zcutoff2.5 0.9929429      c
# Zcutoff3   0.9813810      d
#
# attr(,"class")
# [1] "group"

sumAUC$mapping =factor(sumAUC$mapping,levels=c("polycat","hylite","eaglerc","rsem","kallisto","salmon","bowtie"))
pdf(paste0("s5.bna.plotAUC.sample", sampleSize,".permutation",nPermutation,".pdf"))
# one factor, AUC ~ mapping + normalization + method, data = sumAUC[sumAUC$category =="-", ])
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=mapping, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=normalization, y=AUC )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line())
# show only the significant: method and mapping
ggplot(data=sumAUC[sumAUC$category=="-",], aes(x=method, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC: all gene pairs")+ scale_fill_brewer(palette="Set1")
#+ facet_grid(normalization~.)
# compare connections within network
ggplot(data=sumAUC[sumAUC$method=="Top0.05", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Top0.01", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Top0.005", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Zcutoff1.5", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Zcutoff2", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Zcutoff2.5", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
ggplot(data=sumAUC[sumAUC$method=="Zcutoff3", ], aes(x=category, y=AUC, fill=mapping )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line()) + ggtitle("Boxplot of AUC") + facet_grid(method~.)+ scale_fill_brewer(palette="Set1")
dev.off()

fit<-lm(AUC~mapping+category+method,data=sumAUC)
summary(fit)
out <- duncan.test(a,trt="category"); out
TukeyHSD(x=a, "category", conf.level=0.95)
# dd> ad/- > aa

q()
n
