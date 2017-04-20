library(DiffCorr)
library(corrplot)
library(DGCA)
library(WGCNA) # cor is faster

##############
# FUN: Get DGCA class when given two different matrix of correlation and sample size
##############
getClass<-function(exp_corM, obs_corM, nSample)
{
    ## deal with upper tri only
    # id <- row.names(exp_corM)
    upper<- upper.tri(exp_corM,diag=F)
    # ind <- which( upper , arr.ind = TRUE )
    exp.cor<- exp_corM[upper]
    obs.cor<- obs_corM[upper]
    
    ## calculate individual p and diff.p, then classification
    exp.p <- cor2.test(nSample, exp.cor )
    obs.p <- cor2.test(nSample, obs.cor )
    diff <- compcorr(nSample, exp.cor, nSample, obs.cor)
    # qval <- p.adjust(diff$pval,"BH")
    class<-dCorClass(exp.cor, exp.p, obs.cor, obs.p, 1, sigThresh = 1.01, corSigThresh = 0.05, convertClasses = TRUE)
    aclass <- c(table(class))
    sig<- which(diff$pval<0.05)
    if( length(sig)>0)
    {
        pclass <- c(table(class[sig]))
    }else{pclass=0}
    return(list(p=pclass,a=aclass))
}


# Differential coexpression analysis for all pairs
# 1. How many pairs were differentially coexpressed between A2D5 and ADs
# 2. What are the patterns of DC

####################
## DC calculation ##
####################

setwd("~/jfw-lab/Projects/Eflen_networks/")
rdatafiles<-grep("R-02-dataInput",list.files(),value=TRUE)
rdatafiles
# "R-02-dataInput.hylite_rld.RData"   "R-02-dataInput.hylite_rpkm.RData"
# "R-02-dataInput.polycat_rld.RData"  "R-02-dataInput.polycat_rpkm.RData"
# "R-02-dataInput.rsem_rld.RData"     "R-02-dataInput.rsem_rpkm.RData"

options(scipen=999) # disable scientifc format

chunk_size =1000  #chunk huge matrix to relax memory usage

cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
resS <-c("dataset","genes","genePairs","sigP",paste0("P",cat))
resA <-c("dataset","genes","genePairs","allP",paste0("P",cat))

for(file in rdatafiles)
{
    # file = "R-02-dataInput.polycat_rld.RData"
    print(file)
    flag <- gsub("R-02-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    
    # get expression datasets of exp and obs
    exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    
    # matrix size
    size=nrow(obs)
    nSample=ncol(obs)  # 11
    
    # split the big matrix into chunks, and loop
    schunk <- chunk_size
    n <- floor(size/schunk)
    classP<-rep(0,10)  # store classes of gene pairs with significant DC (p<0.05)
    classA<-rep(0,10)  # store classes of all gene pairs regardless of DC, as background for class enrichment
    for(i in 0:n)
    {
        select1 <- c( (i*schunk+1): min(size,(i+1)*schunk) )
        select2 <- c( (i*schunk+1): size )
        print(paste0("....Calculating gene ",min(select1), " to gene ", max(select1)))
        ss<-proc.time()
        
        ## get correlation matrix, use WGCNA cor, 10min each
        exp_corM <- cor(t(exp[select1,]), t(exp[select2,]), use = 'pairwise.complete.obs', nThread=8)
        obs_corM <- cor(t(obs[select1,]), t(obs[select2,]), use = 'pairwise.complete.obs', nThread=8)
        
        ## summarize class
        pclass <- getClass(exp_corM, obs_corM, nSample)
        classP <- pclass$p+classP
        classA <- pclass$a+classA
        ## count time
        ee<-proc.time()
        run.time  <- ee-ss
        cat("\n ........running minutes: ", run.time[3]/60,"\n\n")
    }
    resS<-rbind(resS,c(flag, size, size*(size-1)/2,sum(classP),classP[-1]))
    print(resS)
    resA<-rbind(resA,c(flag, size, size*(size-1)/2,sum(classA),classA[-1]))
    print(resA)
    save(resS,resA,file="~/DC.all.Rdata")  # save during loop in case error
    gc()
}


res<-data.frame(as.matrix(resA[-1,]))
names(res)<-resA[1,]
resA<-res

res<-data.frame(as.matrix(resS[-1,]))
names(res)<-resS[1,]
resS<-res
save(resS,resA,file="~/DC.all.Rdata")

#####################################
## Enrichment analysis of DC class ##
#####################################

#load("~/DC.all.Rdata")

resA[,-1]<-apply(resA[,-1],2,as.numeric)
resS[,-1]<-apply(resS[,-1],2,as.numeric)
resA<-resA[c(3,4,1,2,5,6),]
resS<-resS[c(3,4,1,2,5,6),]

# Given the huge numbers of AllOut, fisher.test report integer overflow - use sum(as.numeric(.)), use chi square instead
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = 6, ncol = 9);
for(set in 1:6)
{
    for(cat in 1:9)
    {
        sigIn <- resS[set,4+cat]
        allIn <- resA[set,4+cat]
        sigOut <- resS$sigP[set]-sigIn
        allOut <- resA$allP[set]-allIn
        # print(categoryzation)
        categoryzation <- matrix(c(sigIn,allIn, sigOut, allOut), nrow = 2, dimnames = list(c("Sig", "Other"), c("inCat", "outCat")))
        chisq<-chisq.test(categoryzation)
        pTable[set, cat] = ifelse( chisq$residuals[1]>0, -log10(chisq$p.value),0);
        # pTable[set, cat] = -log10(fisher.test(categoryzation, alternative = "greater")$p.value);
    }
}
# pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
pdf("~/differentialCorClasses.all.pdf")
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
yLabels = resS$dataset, xLabels = names(resS)[5:13],
textMatrix = resS[,5:13], colors = blueWhiteRed(100)[50:100],
main = "Differential correlation classes and enrichment",
cex.text = 0.6, cex.lab = 1.5, setStdMargins = FALSE      )


## compile chisq residue table
resiTable = matrix(0, nrow = 6, ncol = 9);
sig<-c()
for(set in 1:6)
{
    tt<-data.frame(sig = as.numeric(resS[set, 5:13]), not.sig=as.numeric(resA[set,5:13])-as.numeric(resS[set, 5:13]))
    rownames(tt)<-names(resA)[5:13]
    chisq <- chisq.test(tt)
    sig[set]<-ifelse(chisq$p.value<0.05,"*"," ")
    resiTable[set,]<-chisq$residuals[,1]
}
colnames(resiTable)<-names(resS)[5:13]
rownames(resiTable)<-paste(sig,resS$dataset)
corrplot(resiTable, is.cor = FALSE, main="Chi-square test residuals for each dataset",mar=c(0,0,4,0))
dev.off()

quit(save="no")

############
## Extras ##
############
# chisq vs fisher's exact
# To test independence, use "Pearson's Chi-squared test with Yates' continuity correction" instead, but chisq p-value doesn't indicate direction of difference, need to check the sign of residue
library(corrplot)
set=1
tt<-data.frame(sig = as.numeric(resS[set, 5:13]), not.sig=as.numeric(resA[set,5:13])-as.numeric(resS[set, 5:13]))
rownames(tt)<-names(resA)[5:13]
chisq <- chisq.test(tt)
chisq
chisq$observed
round(chisq$expected,2)
# visualize Pearson residuals using the package corrplot:
corrplot(chisq$residuals, is.cor = FALSE)
# The sign of the standardized residuals is also very important to interpret the association between rows and columns as explained in the block below.
# Positive residuals are in blue. Positive values in cells specify an attraction (positive association) between the corresponding row and column variables.
# Negative residuals are in red. This implies a repulsion (negative association) between the corresponding row and column variables.
# Contibution in percentage (%)
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
# Visualize the contribution
corrplot(contrib, is.cor = FALSE)
