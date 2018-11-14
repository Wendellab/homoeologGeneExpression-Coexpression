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

# setwd("/work/LAS/jfw-lab/hugj2006/eflen/output")
rdatafiles<-grep("R-01.*.NetworkDatasets.RData",list.files(),value=TRUE)
rdatafiles
# "R-01-hyliteNetworkDatasets.RData"   "R-01-kallistoNetworkDatasets.RData"
# "R-01-polycatNetworkDatasets.RData"  "R-01-rsemNetworkDatasets.RData"
# "R-01-salmonNetworkDatasets.RData"

options(scipen=999) # disable scientifc format

chunk_size =1000  #chunk huge matrix to relax memory usage

cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
resS <-c("dataset","norm","genes","genePairs","sigP",paste0("P",cat))
resA <-c("dataset","norm","genes","genePairs","allP",paste0("P",cat))

for(file in rdatafiles)
{
    message(file)
    flag <- gsub("R-01-|NetworkDatasets.RData","",file)
    fn<-load(file) #  "coldata"      "networks"     "networks.rld" "networks.tpm"
    for(norm in c("rld","log2rpkm")){
        message(norm)
        if(norm=="rld"){
            multiExpr = networks.rld
            # get expression datasets of exp and obs
            exp<-as.data.frame(multiExpr[["A2D5"]])
            obs<-as.data.frame(multiExpr[["ADs"]])
        }else{
            if(flag%in%c("salmon","kallisto")){
                multiExpr = networks.tpm; message("tpm")
            }else{
                multiExpr = networks.rpkm
            }
            # get expression datasets of exp and obs
            exp<-as.data.frame(log2(multiExpr[["A2D5"]]+1))
            obs<-as.data.frame(log2(multiExpr[["ADs"]]+1))
        }
        # matrix size
        size=nrow(obs)
        nSample=ncol(obs)  # 33
        
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
        resS<-rbind(resS,c(flag, norm, size, size*(size-1)/2,sum(classP),classP[-1]))
        print(resS)
        resA<-rbind(resA,c(flag, norm, size, size*(size-1)/2,sum(classA),classA[-1]))
        print(resA)
        save(resS,resA,file="s4.DC.all.Rdata")  # save during loop in case error
        gc()
    }
}

res<-data.frame(as.matrix(resA[-1,]))
names(res)<-resA[1,]
resA<-res

res<-data.frame(as.matrix(resS[-1,]))
names(res)<-resS[1,]
resS<-res

resA[,-c(1,2)]<-apply(resA[,-c(1,2)],2,as.numeric)
resS[,-c(1,2)]<-apply(resS[,-c(1,2)],2,as.numeric)

rownames(resA)=NULL
rownames(resS)=NULL
save(resS,resA,file="s4.DC.all.Rdata")

#####################################
## Enrichment analysis of DC class ##
#####################################

load("s4.DC.all.Rdata")

# Given the huge numbers of AllOut, fisher.test report integer overflow - use sum(as.numeric(.)), use chi square instead
# Initialize tables of p-values and of the corresponding counts
cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
pTable = matrix(0, nrow = nrow(resA), ncol = length(cat));
for(set in 1:nrow(resA))
{
    for(c in 1:length(cat))
    {
        sigIn <- resS[set,5+c]
        allIn <- resA[set,5+c]
        sigOut <- resS$sigP[set]-sigIn
        allOut <- resA$allP[set]-allIn
        # print(categoryzation)
        categoryzation <- matrix(c(sigIn,allIn, sigOut, allOut), nrow = 2, dimnames = list(c("Sig", "Other"), c("inCat", "outCat")))
        chisq<-chisq.test(categoryzation)
        pTable[set, c] = ifelse( chisq$residuals[1]>0, -log10(chisq$p.value),0);
        # pTable[set, cat] = -log10(fisher.test(categoryzation, alternative = "greater")$p.value);
    }
}
pTable
# pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
pdf("s4.DC.all.pdf")
library(gplots)
textplot(resA)
textplot(resS)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
yLabels = resS$dataset, xLabels = names(resS)[5:13],
textMatrix = resS[,5:13], colors = blueWhiteRed(100)[50:100],
main = "Differential correlation classes and enrichment",
cex.text = 0.6, cex.lab = 1.5, setStdMargins = FALSE      )


## compile chisq residue table
resiTable = matrix(0, nrow = nrow(resA), ncol = length(cat));
sig<-c()
for(set in 1:nrow(resA))
{
    tt<-data.frame(sig = as.numeric(resS[set, 6:14]), not.sig=as.numeric(resA[set,6:14])-as.numeric(resS[set, 6:14]))
    rownames(tt)<-names(resA)[6:14]
    chisq <- chisq.test(tt)
    sig[set]<-ifelse(chisq$p.value<0.05,"*"," ")
    resiTable[set,]<-chisq$residuals[,1]
}
colnames(resiTable)<-names(resS)[6:14]
rownames(resiTable)<-paste(sig,resS$dataset,resS$norm)
corrplot(resiTable, is.cor = FALSE, main="Chi-square test residuals for each dataset",mar=c(0,0,4,0))
dev.off()

quit(save="no")


##############################################
## Identify differential coexpression genes ##
##############################################
options(scipen=999) # disable scientifc format
# p - the percentage of DC among all possible gene pairs, measures the extent of differential co-expression
# n - the number of all possible pairs for a given gene
# k - the number of DC pairs for a given gene
# P - the probability for a given gene to be differentially co-expressed, considering that the distribution of k significant in n possible pairs follows the binomial distribution model
# P = pbinom(k, size=n, prob=p, lower.tail=FALSE)

# get n and p from early Dc analysis
load("s4.DC.all.Rdata") # resS, resA
n.list <- resA$genes
names(n.list) <- paste(resA$dataset,resA$norm,sep="_")
p.list <- resS$sigP/resS$genePairs
names(p.list) <- paste(resS$dataset,resS$norm,sep="_")

# setwd("/work/LAS/jfw-lab/hugj2006/eflen/output")
rdatafiles<-grep("R-01.*.NetworkDatasets.RData",list.files(),value=TRUE)
rdatafiles
# "R-01-hyliteNetworkDatasets.RData"   "R-01-kallistoNetworkDatasets.RData"
# "R-01-polycatNetworkDatasets.RData"  "R-01-rsemNetworkDatasets.RData"
# "R-01-salmonNetworkDatasets.RData"

# get k and P for each dataset
for(file in rdatafiles)
{
    message(file)
    flag <- gsub("R-01-|NetworkDatasets.RData","",file)
    fn<-load(file) #  "coldata"      "networks"     "networks.rld" "networks.tpm"
    for(norm in c("rld","log2rpkm"))
    {
        message(norm)
        if(norm=="rld"){
            multiExpr = networks.rld
            # get expression datasets of exp and obs
            exp<-as.data.frame(multiExpr[["A2D5"]])
            obs<-as.data.frame(multiExpr[["ADs"]])
        }else
        {
            if(flag%in%c("salmon","kallisto")){
                multiExpr = networks.tpm; message("tpm")
            }else{
                multiExpr = networks.rpkm
            }
            # get expression datasets of exp and obs
            exp<-as.data.frame(log2(multiExpr[["A2D5"]]+1))
            obs<-as.data.frame(log2(multiExpr[["ADs"]]+1))
        }
        
        # n, p
        item=paste(flag,norm,sep="_")
        n<-n.list[item]
        p<-p.list[item]
        nSample <- ncol(exp)
        if(n!=nrow(obs) ) break
        
        # now get k, loop by gene
        require(svMisc)
        dcg <- data.frame(genes = rownames(exp), P=NA)
        for (i in 1:n)
        {
            ## get correlation matrix, use WGCNA cor
            exp_corL <- cor(t(exp[i,]), t(exp[-i,]), use = 'pairwise.complete.obs', nThread=8)
            obs_corL <- cor(t(obs[i,]), t(obs[-i,]), use = 'pairwise.complete.obs', nThread=8)
            
            ## calculate individual p and diff.p, then classification
            # exp_corL.p <- cor2.test(nSample, exp_corL )
            # obs_corL.p <- cor2.test(nSample, obs_corL )
            diff <- compcorr(nSample, exp_corL, nSample, obs_corL)
            k <- length(which(as.numeric(diff$pval)<0.05))
            dcg$P[i] = pbinom(k, size=n-1, prob=p, lower.tail=FALSE)
            
            ## report loop progress
            progress(i, max.value=n)
            if (i == n) cat("Done!\n")
        }
        
        # display results
        print(item)
        print(table(dcg$P<0.05))
        dcg$q.bh<- p.adjust(dcg$P, "BH")
        print(table(dcg$q.bh<0.05))
        assign(paste0(item,".dcg"),dcg)

    }
}

DCG=list(hylite_log2rpkm.dcg, hylite_rld.dcg, kallisto_log2rpkm.dcg, kallisto_rld.dcg, polycat_log2rpkm.dcg, polycat_rld.dcg, rsem_log2rpkm.dcg, rsem_rld.dcg, salmon_log2rpkm.dcg, salmon_rld.dcg)
names(DCG) = c("hylite_log2rpkm", "hylite_rld", "kallisto_log2rpkm", "kallisto_rld", "polycat_log2rpkm", "polycat_rld", "rsem_log2rpkm", "rsem_rld", "salmon_log2rpkm", "salmon_rld")
lapply(DCG, function(x)table(x$q.bh<0.05) )

resA$DCgene=lapply(DCG,function(x)length(which(x$q.bh<0.05)))[paste(resA$dataset,resA$norm,sep="_")]
resA$DCgenePerc =lapply(DCG,function(x)length(which(x$q.bh<0.05))/nrow(x))[paste(resA$dataset,resA$norm,sep="_")]

save(resS, resA, DCG, file="s4.DC.all.Rdata")

---------book 10/31/18
###########################
## Possible causes of DC ##
###########################
m<-load("s5.DC.all.Rdata")
m # "resS" "resA" "DCG"

# E, A, D
load("s2.assign_eval.hylite.Rdata")
load("s2.assign_eval.polycat.Rdata")
load("s2.assign_eval.rsem.Rdata")
source("addTrans.r")

# library(VennDiagram)
library(gplots)
library(scatterplot3d)
library(scales)
library(RColorBrewer)
show_col(hue_pal()(9))
show_col(brewer.pal(n = 9, name = "Paired"))
colors<-c(addTrans("grey", 150), addTrans("blue", 150), addTrans("grey", 150), addTrans("red", 150))
show_col(colors)

for(i in names(DCG))
{
    x<-DCG[[i]]
    x$class <- "NonSig"
    x$class[x$q.bh<0.05]<-"Sig"
    x$genome <- "At"
    x$genome[grep("d$",x$genes)] <- "Dt"
    x$gorai <- gsub(".$","",x$genes)
    
    # combine DC results with asignment metrics
    metrics <-get(gsub("_.*","",i))
    metrics$gene <-gsub("a$|d$","",metrics$gene)
    tt<-aggregate(metrics[,-c(1:4)],by=list(metrics$gene), function(x)mean(x,na.rm=TRUE))
    rownames(tt)<-tt$Group.1
    y<- cbind(x,tt[x$gorai,])
    
    #save image
    pdf(paste0("s5.DC.all.metrics.",i,".pdf"))
    
    # plot venn for overlap between At and Dt DC genes
    venn(list(At=x$gorai[x$genome=="At"&x$class=="Sig"], Dt=x$gorai[x$genome=="Dt"&x$class=="Sig"]))
    #venn.diagram(x=list(At=x$gorai[x$genome=="At"&x$class=="Sig"], Dt=x$gene[x$genome=="Dt"&x$class=="Sig"]), filename="test.tiff",fill = c("light blue", "pink"))
    
    # boxplot by DC class
    boxplot(geneLenM~genome+class,data=y, main ="Gene Length")
    boxplot(percentageEffectM~genome+class,data=y, main="Average percentage effiective region")
    boxplot(expression~genome+class,data=y, main="Average expression")
    boxplot(expression.log2~genome+class,data=y, main="Average log2 expression")
    boxplot(efficiency~genome+class,data=y, main="Efficiency")
    boxplot(accuracy~genome+class,data=y, main="Accuracy")
    boxplot(discrepancy0~genome+class,data=y, main="Discrepancy")
    
    # 3d scatter plot of significant changes
    attach(y)
    cc<-as.factor(paste(genome,class))
    levels(cc)
    scatterplot3d(accuracy, efficiency, discrepancy0, color=colors[cc], type="h", pch=20)
    legend("right", legend = levels(cc), col = colors, pch = 16, bg="white")
    scatterplot3d(geneLenM, percentageEffectM, expression.log2, color=colors[cc], type="h",pch=20, xlim=c(0,3000))
    legend("right", legend = levels(cc), col = colors, pch = 16, bg="white")
    detach(y)
    dev.off()
}

