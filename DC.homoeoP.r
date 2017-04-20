library(DiffCorr)
library(DGCA)
library(WGCNA)
library(corrplot)


# Differential coexpression analysis for homoeolog pairs
# 1. How many pairs were differentially coexpressed between A2D5 and ADs
# 2. What are the patterns of DC


setwd("~/jfw-lab/Projects/Eflen_networks/")
rdatafiles<-grep("R-02-dataInput",list.files(),value=TRUE)
rdatafiles
# "R-02-dataInput.hylite_rld.RData"   "R-02-dataInput.hylite_rpkm.RData"
# "R-02-dataInput.polycat_rld.RData"  "R-02-dataInput.polycat_rpkm.RData"
# "R-02-dataInput.rsem_rld.RData"     "R-02-dataInput.rsem_rpkm.RData"


ss<-proc.time()

nSample=11

for(file in rdatafiles)
{
    # file = "R-02-dataInput.polycat_rld.RData"
    print(file)
    flag <- gsub("R-02-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    
    # get expression datasets of exp and obs
    exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    
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
    print( table(class))
    print( table(class)/length(class) )
    print("Differentially co-expressed: ")
    print( table(pValDiff.adjust<0.05)/length(class)  )
    print( table(class[pValDiff.adjust<0.05])/length(class[pValDiff.adjust<0.05])  )
    print( table(diff$pval<0.05)/length(class)  )
    print( table(class[diff$pval<0.05])/length(class[diff$pval<0.05])  )
    
    assign(flag,resF)
}
ee<-proc.time()
run.time  <- ee-ss
cat("\n Number of minutes running:", run.time[3]/60,"\n\n")


save(list=gsub("R-02-dataInput.|.RData","",rdatafiles), file = paste0("~/DC.homoeoPair0331.Rdata"))

#### Post Analysis
## Analysis of homoeolog pairs
m<-load("~/DC.homoeoPair0331.Rdata")
m # "hylite_rld"   "hylite_rpkm"  "polycat_rld"  "polycat_rpkm" "rsem_rld"   "rsem_rpkm"

# summarize sig results
cat<-c("+/+","+/0","+/-","0/+", "0/0", "0/-", "-/+",  "-/0", "-/-")
res<-c("dataset","genes","sigP",paste0("P",cat))
for(i in m)
{
    x<-get(i)
    temp<-c(i, nrow(x), length(which(x$pValDiff<0.05)), as.numeric(table(x$class[x$pValDiff<0.05]))[-1])
    res<-rbind(res,temp)
}
res
resS<-data.frame(as.matrix(res[-1,]))
colnames(resS) <-res[1,]

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

#
resA[,-1]<-apply(resA[,-1],2,as.numeric)
resS[,-1]<-apply(resS[,-1],2,as.numeric)
resA<-resA[c(3,4,1,2,5,6),]
resS<-resS[c(3,4,1,2,5,6),]

# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = 6, ncol = 9);
for(set in 1:6)
{
    for(cat in 1:9)
    {
        sigIn <- resS[set,3+cat]
        allIn <- resA[set,3+cat]
        sigOut <- resS$sigP[set]-sigIn
        allOut <- resA$allP[set]-allIn
        categoryzation <- matrix(c(sigIn,allIn, sigOut, allOut), nrow = 2, dimnames = list(c("Sig", "Other"), c("inCat", "outCat")))
        pTable[set, cat] = -log10(fisher.test(categoryzation, alternative = "greater")$p.value);
    }
}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
pdf("~/differentialCorClasses.homoeolog.pdf")
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
yLabels = resS$dataset, xLabels = names(resS)[4:12],
textMatrix = resS[,4:12], colors = blueWhiteRed(100)[50:100],
main = "Differential correlation classes and enrichment",
cex.text = 1, cex.lab = 1.5, setStdMargins = FALSE      )

resiTable = matrix(0, nrow = 6, ncol = 9);
sig<-c()
for(set in 1:6)
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
dev.off()




quit(save="no")
