################
#### bowtie ####
################
# Additional piplines for 2020 revision

load("R-01-bowtieDatasets.RData")
# A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt,
load("R-01-bowtieNetworkDatasets.RData")
# coldata, networks, networks.rld, networks.tpm

s2o=function(x){unlist(lapply(strsplit(x,"[.]"),function(x)paste(paste0(x[2]),gsub("seed","",x[1]),x[3],sep="-")) )}

library(data.table)

At.T = ADs.At
At.F = ADs.At
Dt.T = ADs.Dt
Dt.F = ADs.Dt
for(i in 23:33){
    flag= names(ADs.At)[i]
    print(flag)
    ADflag =s2o(flag)
    A2flag=s2o(names(A2.Total)[i])
    D5flag=s2o(names(D5.Total)[i])
    # output ADs read ID
    # samtools view -F 4 ADs-10-R1.q10.sort.bam | cut -f1,3 >AD.id.txt; cut -f1 AD.id.txt >AD.txt
    system(paste0("samtools view -F 4 ../bowtie2only/bamQ10/", ADflag, ".q10.sort.bam | cut -f1,3 >AD.id.txt; cut -f1 AD.id.txt >AD.txt"))
    # out put A2 read ID
    # zcat ../../seedFastq/A2-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>A2.txt
    # zcat ../../seedFastq/D5-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>D5.txt
    system(paste0("zcat ../seedFastq/",A2flag,".cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>A2.txt"))
    system(paste0("zcat ../seedFastq/",D5flag,".cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>D5.txt"))
    # label AD reads
    system("perl labelAD.pl A2.txt D5.txt AD.txt | paste AD.id.txt - >AD.origin.txt")
    #
    x=fread("AD.origin.txt",select=c(2,3), header=FALSE)
    names(x)=c("id","origin")
    x$ref = gsub(".*[.]","",x$id)
    x$id=gsub(".1.D","d", gsub(".1.A","a",x$id))
    At.T[,flag]  = as.numeric( table(x$id[x$ref=="A" & x$origin=="A2"])[rownames(A2.Total)] )
    At.F[,flag]  = as.numeric( table(x$id[x$ref=="A" & x$origin=="D5"])[rownames(A2.Total)] )
    Dt.T[,flag]  = as.numeric( table(x$id[x$ref=="D" & x$origin=="D5"])[rownames(D5.Total)] )
    Dt.F[,flag]  = as.numeric( table(x$id[x$ref=="D" & x$origin=="A2"])[rownames(D5.Total)] )
    # check if make sense, should be zero if correct
    print(sum(At.T[,flag]+At.F[,flag] - ADs.At[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag] - ADs.Dt[,flag], na.rm=T))
}
for(i in 1:22){
    flag= names(ADs.At)[i]
    print(flag)
    ADflag =gsub("[.]","-",flag)
    A2flag=gsub("[.]","-",names(A2.Total)[i])
    D5flag=gsub("[.]","-",names(D5.Total)[i])
    # output ADs read ID
    # samtools view -F 4 ADs-10-R1.q10.sort.bam | cut -f1,3 >AD.id.txt; cut -f1 AD.id.txt >AD.txt
    system(paste0("samtools view -F 4 ../bowtie2only/bamQ10/", ADflag, ".q10.sort.bam | cut -f1,3 >AD.id.txt; cut -f1 AD.id.txt >AD.txt"))
    # out put A2 read ID
    # zcat ../../seedFastq/A2-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>A2.txt
    # zcat ../../seedFastq/D5-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>D5.txt
    system(paste0("zcat ../flowerFastq/",A2flag,"_1.fq.gz |grep '^@K'|sed 's/ .*//g' |sed 's/@//g'>A2.txt"))
    system(paste0("zcat ../flowerFastq/",D5flag,"_1.fq.gz |grep '^@K'|sed 's/ .*//g' |sed 's/@//g'>D5.txt"))
    # label AD reads
    system("perl labelAD.pl A2.txt D5.txt AD.txt | paste AD.id.txt - >AD.origin.txt")
    #
    x=fread("AD.origin.txt",select=c(2,3), header=FALSE)
    names(x)=c("id","origin")
    x$ref = gsub(".*[.]","",x$id)
    x$id=gsub(".1.D","d", gsub(".1.A","a",x$id))
    At.T[,flag]  = as.numeric( table(x$id[x$ref=="A" & x$origin=="A2"])[rownames(A2.Total)] )
    At.F[,flag]  = as.numeric( table(x$id[x$ref=="A" & x$origin=="D5"])[rownames(A2.Total)] )
    Dt.T[,flag]  = as.numeric( table(x$id[x$ref=="D" & x$origin=="D5"])[rownames(D5.Total)] )
    Dt.F[,flag]  = as.numeric( table(x$id[x$ref=="D" & x$origin=="A2"])[rownames(D5.Total)] )
    # check if make sense: PE reads counted likely twice ofshould be zero if correct
    print(sum(At.T[,flag]+At.F[,flag], na.rm=T)/sum(ADs.At[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag], na.rm=T)/sum(ADs.Dt[,flag], na.rm=T))
}


At.T[is.na(At.T)] = 0
At.F[is.na(At.F)] = 0
Dt.T[is.na(Dt.T)] = 0
Dt.F[is.na(Dt.F)] = 0
save( At.T, Dt.T, At.F, Dt.F, file="R-01-bowtieDatasets.true.RData")
q()
n

################################
## Read assignment evaluation ##
################################
source("2.1.FUN.r")

load("R-01-bowtieDatasets.RData")->l;l
# A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt,

load("R-01-bowtieNetworkDatasets.RData")->l;l
# coldata, networks, networks.rld, networks.tpm

load("R-01-bowtieDatasets.true.RData")->l;l
# At.T, Dt.T, At.F, Dt.F

# see if meets expectation
colSums(At.F+At.T-ADs.At) # SE 0, PE is not 0
colSums(Dt.F+Dt.T-ADs.Dt) # SE 0, PE is not 0, as expected

## calculate assignment efficiency - percentage of total reads assigned to At and Dt homoeologs
total <- D5.Total + A2.Total
totalA <- A2.Total
totalD <- D5.Total

assigned <- ADs.At + ADs.Dt
assignedA <- ADs.At
assignedD <- ADs.Dt

efficiency <- assigned/total
efficiencyA <- assignedA/totalA
efficiencyD <- assignedD/totalD

quantile(as.numeric(as.matrix(efficiency)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(efficiencyA)),na.rm=TRUE) #NaN from 0/0, no way to learn; Inf indicating incorrect assignment to zero total
quantile(as.numeric(as.matrix(efficiencyD)),na.rm=TRUE) #NaN from 0/0, no way to learn; Inf indicating incorrect assignment

# summary table
overall <- c(sum(assigned)/sum(total), sum(assignedA)/sum(totalA), sum(assignedD)/sum(totalD))
bySamples<-cbind(apply(assigned,2,sum)/apply(total,2,sum),apply(assignedA,2,sum)/apply(totalA,2,sum),apply(assignedD,2,sum)/apply(totalD,2,sum))
resEffi <- as.data.frame(rbind(overall, bySamples))
names(resEffi) <-c("total","At","Dt")
resEffi

## calculate assignment discrepancy - percentage of absolute diffence between observed and expected values  among assigned reads
# discrepancy is an estimate very similar to inaccuracy, but reflect the observed differences in read counts, while inaccuracy reflects the inherent assignment actions. Inaccurate read assignments both for At and Dt can sometimes cancel out each other and do not lead to discrepancy between expected and observed count numbers.
different <- abs(ADs.At-A2.Total)+abs(ADs.Dt-D5.Total)
differentA <- abs(ADs.At-A2.Total)
differentD <- abs(ADs.Dt-D5.Total)

discrepancy <- different/total
discrepancyA <- differentA/totalA
discrepancyD <- differentD/totalD

quantile(as.numeric(as.matrix(discrepancy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(discrepancyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(discrepancyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

summary(as.numeric(as.matrix(discrepancy)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.00    0.01    0.05    0.16    0.20    2.00  193731
# NaN meaning no discrepancy
summary(as.numeric(as.matrix(discrepancyA)))
summary(as.numeric(as.matrix(discrepancyD))) # inf means wrong assignment to zero expression

# summary table
overall <- c(sum(different)/sum(total), sum(differentA)/sum(totalA), sum(differentD)/sum(totalD))
bySamples<-cbind(apply(different,2,sum)/apply(total,2,sum),apply(differentA,2,sum)/apply(totalA,2,sum),apply(differentD,2,sum)/apply(totalD,2,sum))
resDisc <- as.data.frame(rbind(overall, bySamples))
names(resDisc) <-c("total","At","Dt")
resDisc


## binary classification evaluation
# A as TRUE
resBinA=binC(TP = At.T, TN = Dt.T, FP = At.F ,FN = Dt.F)
# D as TRUE
resBinD=binC(TP = Dt.T, TN = At.T, FP = Dt.F ,FN = At.F)
# summary
resBinA$summary
resBinD$summary
###############

## prepare explainatory variables for analysis
info<-coldata
info$type <- c(rep("PE",22), rep("SE",11))
info$ADsLibSize <- colSums(ADs.At+ADs.Dt)

# geneLenM, output from "detectEffectiveRegion.r"
system("cut -f3,4 eflen.snp.rl300.txt|paste eflen.snp.rl100.txt - >eflenList.snp.txt")
geneL <- read.table("eflenList.snp.txt",header=TRUE,sep="\t")
len<- geneL$trueLen
names(len)<-geneL$id
# length for 37223 genes
len<-len[gsub("a$","",rownames(A2.Total))]

# PercentageEffectM
ratio100<- geneL$percentEffective
names(ratio100)<-geneL$id
ratio300<- geneL$percentEffective.1
names(ratio300)<-geneL$id
# effective_length/true_length ratio for 37223 genes
ratio100<-ratio100[gsub("a$","",rownames(A2.Total))]
ratio300<-ratio300[gsub("a$","",rownames(A2.Total))]
# length corresponding to 37223 genes X 11 samples matrix
info$effectRegionRatio<-NA
info$effectRegionRatio[info$type=="PE"]<-"ratio300"
info$effectRegionRatio[info$type=="SE"]<-"ratio100"

## prepare explainatory variables -  expression
# total <- as.numeric(as.matrix(D5.Total + A2.Total))
expression<-total
expression.log2<-log2(expression+1)

# Estimation metrics of Efficiency, Accuracy and Discrepancy are data frames of 37223 genes x 33 samples
# Categorical explanatory variables for 33 samples include tissue (12), sequency type (SE&PE), and rep
# Continuous variables for 37223 genes include gene length and effective region ratio

## prepare explanatory variables as numeric data
gene <- rep(rownames(A2.Total),33)
sample <- rep(coldata$sample, each=37223)
tissue <- rep(coldata$tissue, each=37223)
rep <- rep(as.numeric(gsub("R","",coldata$rep)), each=37223)
# length corresponding to 37223 genes X 11 samples matrix
geneLenM <- rep(len, 33)
percentageEffectM <- c( rep(ratio300,22),rep(ratio100,11) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(resBinA$accuracy)), mcc=as.numeric(as.matrix(resBinA$MCC)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)),  as.numeric(as.matrix(resBinA$precision)), as.numeric(as.matrix(resBinD$precision)), as.numeric(as.matrix(resBinA$recall)), as.numeric(as.matrix(resBinD$recall)), as.numeric(as.matrix(resBinA$F1)), as.numeric(as.matrix(resBinD$F1)) )
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","precisionA", "precisionD", "recallA", "recallD", "F1A", "F1D")

## plot and save results
bowtie<-res
plotRes(bowtie, file="s2.eval.bowtie.pdf")

bowtieH<-resH
plotRes.homoeo(bowtieH, file="s2.eval.bowtie.homoeolog.pdf")

library(gplots)
pdf("s2.eval.bowtie.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()

bowtieSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)
save(bowtie, bowtieH, bowtieSummary, file="s2.assign_eval.bowtie.Rdata")
