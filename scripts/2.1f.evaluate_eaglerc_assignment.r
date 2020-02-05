# Additional piplines for 2020 revision

############### FUNC ###############
source("2.1.FUN.r")

################
#### eaglerc ####
################

lnames <- load('R-01-eaglercDatasets.RData')
lnames
#  "A2.At"     "A2.Dt"     "A2.N"      "D5.At"     "D5.Dt"     "D5.N"  "A2.Total"  "D5.Total"  "ADs.Total" "ADs.At"    "ADs.Dt"    "ADs.N"
lnames <- load('R-01-eaglercNetworkDatasets.RData')
lnames
# "coldata"      "networks"     "networks.rld" "networks.tpm"

# see if meets expectation
quantile(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
colSums(ADs.At + ADs.Dt)-colSums(A2.At + A2.Dt + D5.At + D5.Dt)
quantile(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
colSums(ADs.Total)-colSums(A2.Total + D5.Total)
# not perfect 0, but close enough

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
resBinA=binC(TP = A2.At, TN = D5.Dt, FP = D5.At ,FN = A2.Dt)
# D as TRUE
resBinD=binC(TP = D5.Dt, TN = A2.At, FP = A2.Dt ,FN = D5.At)
# summary
resBinA$summary
resBinD$summary
###############

## prepare explainatory variables for analysis
info<-coldata
info$type <- c(rep("PE",22), rep("SE",11))
info$ADsLibSize <- colSums(ADs.Total)

# geneLenM, output from "detectEffectiveRegion.r"
system("cut -f3,4 eflen.vcf.rl300.txt|paste eflen.vcf.rl100.txt - >eflenList.vcf.txt")
geneL <- read.table("eflenList.vcf.txt",header=TRUE,sep="\t")
len<- geneL$trueLen
names(len)<-geneL$id
# length for 37223 genes
len<-len[gsub("[.]1$","",rownames(A2.Total))]

#\ percentageEffectM
ratio100<- geneL$percentEffective
names(ratio100)<-geneL$id
ratio300<- geneL$percentEffective.1
names(ratio300)<-geneL$id
# effective_length/true_length ratio for 37223 genes
ratio100<-ratio100[gsub("[.]1$","",rownames(A2.Total))]
ratio300<-ratio300[gsub("[.]1$","",rownames(A2.Total))]
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
sample <- rep(coldata$sample, each=nrow(A2.Total))
tissue <- rep(coldata$tissue, each=nrow(A2.Total))
rep <- rep(as.numeric(gsub("R","",coldata$rep)), each=nrow(A2.Total))
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
eaglerc<-res
plotRes(eaglerc, file="s2.eval.eaglerc.pdf")

eaglercH<-resH
plotRes.homoeo(eaglercH, file="s2.eval.eaglerc.homoeolog.pdf")

library(gplots)
pdf("s2.eval.eaglerc.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()

eaglercSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)
save(eaglerc, eaglercH, eaglercSummary, file="s2.assign_eval.eaglerc.Rdata")
