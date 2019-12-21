library(DiffCorr)
library(DGCA)
library(WGCNA)
library(corrplot)
options(scipen=999) # disable scientifc format


# Differential coexpression analysis for homoeolog pairs
# 1. How many pairs were differentially coexpressed between A2D5 and ADs
# 2. What are the patterns of DC


# setwd("/work/LAS/jfw-lab/hugj2006/eflen/output")
rdatafiles<-grep("R-01.*.NetworkDatasets.RData",list.files(),value=TRUE)
rdatafiles
# "R-01-hyliteNetworkDatasets.RData"   "R-01-kallistoNetworkDatasets.RData"
# "R-01-polycatNetworkDatasets.RData"  "R-01-rsemNetworkDatasets.RData"
# "R-01-salmonNetworkDatasets.RData"

ss<-proc.time()
nSample=33
m=c()
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
        #
        id<-row.names(exp)
        At<-grep("a$",id)
        Dt<-grep("d$",id)
        id<- gsub(".$","",id)
        idP<- intersect(id[At],id[Dt])
        # homoeolog pairs
        res <- data.frame(At=paste0(idP,"a"),Dt=paste0(idP,"d"))
        exp.cor <- apply(res,1,function(x) cor(as.numeric(exp[x[1],]),as.numeric(exp[x[2],])))
        exp.p <- cor2.test(nSample, exp.cor )
        obs.cor <- apply(res,1,function(x) cor(as.numeric(obs[x[1],]),as.numeric(obs[x[2],])))
        obs.p <- cor2.test(nSample, obs.cor )
        diff <- compcorr(nSample, exp.cor, nSample, obs.cor)
        pValDiff.adjust <- p.adjust(diff$pval,"BH")
        class<-dCorClass(exp.cor, exp.p, obs.cor, obs.p, diff$pval, sigThresh = 1.01, corSigThresh = 0.05, convertClasses = TRUE)
        resF<-data.frame( At=res$At, Dt=res$Dt, exp.cor, exp.p, obs.cor, obs.p, zScoreDiff=diff$diff,pValDiff=diff$pval,pValDiff.adjust,class)
        
        print("Overall classification: ")
        print( round(table(class)/length(class)*100,1))
        
        assign(paste0(flag,"_",norm),resF)
        m=c(m,paste0(flag,"_",norm))
    }
}
ee<-proc.time()
run.time  <- ee-ss
cat("\n Number of minutes running:", run.time[3]/60,"\n\n")
save(list=m, file = paste0("s4.DC.homoeoPair.Rdata"))

###################
## Post Analysis ##
###################m

## Analysis of homoeolog pairs
mm<-load("s4.DC.homoeoPair.Rdata")
mm # "hylite_rld"   "hylite_rpkm"  "polycat_rld"  "polycat_rpkm" "rsem_rld"   "rsem_rpkm"

# summarize sig results
cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
res<-c("dataset","genes","sigP",paste0("P",cat))
for(i in mm)
{
    x<-get(i)
    temp<-c(i, nrow(x), length(which(x$pValDiff<0.05)), as.numeric(table(x$class[x$pValDiff<0.05]))[-1])
    res<-rbind(res,temp)
}
res
resS<-data.frame(as.matrix(res[-1,]))
colnames(resS) <-res[1,]
rownames(resS)=NULL
resS[,-1]<-apply(resS[,-1],2,as.numeric)

# summarize Not sig results
cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
res<-c("dataset","genes","allP",paste0("A",cat))
for(i in m)
{
    x<-get(i)
    temp<-c(i, nrow(x), nrow(x), as.numeric(table(x$class))[-1])
    res<-rbind(res,temp)
}
res
resA<-data.frame(as.matrix(res[-1,]))
colnames(resA) <-res[1,]
rownames(resA)=NULL
resA[,-1]<-apply(resA[,-1],2,as.numeric)


# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nrow(resS), ncol = length(cat));
for(set in 1:nrow(resS))
{
    for(c in 1:length(cat))
    {
        sigIn <- resS[set,3+c]
        allIn <- resA[set,3+c]
        sigOut <- resS$sigP[set]-sigIn
        allOut <- resA$allP[set]-allIn
        categoryzation <- matrix(c(sigIn,allIn, sigOut, allOut), nrow = 2, dimnames = list(c("Sig", "Other"), c("inCat", "outCat")))
        pTable[set, c] = -log10(fisher.test(categoryzation, alternative = "greater")$p.value);
    }
}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
# pTable[pTable>50 ] = 50 ;
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
pdf("s4.DC.homoeoPair.pdf")
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
labeledHeatmap( Matrix = pTable, colorLabels = TRUE, yLabels = resS$dataset, xLabels = names(resS)[4:12],
textMatrix = resS[,4:12], colors = blueWhiteRed(100)[50:100], main = "Differential correlation classes and enrichment", cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
#rld
labeledHeatmap( Matrix = pTable[grep("rld",resS$dataset),], colorLabels = TRUE, yLabels = resS$dataset[grep("rld",resS$dataset)], xLabels = names(resS)[4:12],
textMatrix = resS[grep("rld",resS$dataset),4:12], colors = blueWhiteRed(100)[50:100], main = "Differential correlation classes and enrichment, rld", cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
#rpkm
labeledHeatmap( Matrix = pTable[grep("log2rpkm",resS$dataset),], colorLabels = TRUE, yLabels = resS$dataset[grep("log2rpkm",resS$dataset)], xLabels = names(resS)[4:12],
textMatrix = resS[grep("log2rpkm",resS$dataset),4:12], colors = blueWhiteRed(100)[50:100], main = "Differential correlation classes and enrichment, log2rpkm", cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )

resiTable = matrix(0, nrow = nrow(resS), ncol = length(cat));
sig<-c()
for(set in 1:nrow(resS))
{
    tt<-data.frame(sig = as.numeric(resS[set, 4:12]), not.sig=as.numeric(resA[set,4:12])-as.numeric(resS[set, 4:12]))
    rownames(tt)<-names(resA)[4:12]
    chisq <- chisq.test(tt)
    sig[set]<-ifelse(chisq$p.value<0.05,"*"," ")
    resiTable[set,]<-chisq$residuals[,1]
}
colnames(resiTable)<-names(resS)[4:12]
rownames(resiTable)<-paste(sig,resS$dataset)
resiTable[is.na(resiTable)]<-0
corrplot(resiTable, is.cor = FALSE, main="Chi-square test residuals for each dataset",mar=c(0,0,4,0))
corrplot(resiTable[grep("log2rpkm",resS$dataset),], is.cor = FALSE, main="Chi-square test residuals for each dataset, log2rpkm",mar=c(0,0,4,0))
corrplot(resiTable[grep("rld",resS$dataset),], is.cor = FALSE, main="Chi-square test residuals for each dataset, rld",mar=c(0,0,4,0))
dev.off()




quit(save="no")
