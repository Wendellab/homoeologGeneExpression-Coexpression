# Additional piplines for 2020 revision

############### FUNC ###############
## binary classification evaluation
binC = function(TP,TN,FP,FN){
    precision=TP/(TP+FP)
    recall=TP/(TP+FN)
    accuracy = (TP+TN)/(TP+TN+FP+FN)
    F1=2*precision*recall/(precision+recall)
    MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    overall <- c(mean(as.numeric(as.matrix(precision)), na.rm=TRUE), mean(as.numeric(as.matrix(recall)), na.rm=TRUE),  mean(as.numeric(as.matrix(accuracy)), na.rm=TRUE), mean(as.numeric(as.matrix(F1)), na.rm=TRUE),  mean(as.numeric(as.matrix(MCC)), na.rm=TRUE) )
    bySamples<-cbind( apply(precision,2,function(x)mean(x,na.rm=TRUE)), apply(recall,2,function(x)mean(x,na.rm=TRUE)), apply(accuracy,2,function(x)mean(x,na.rm=TRUE)), apply(F1,2,function(x)mean(x,na.rm=TRUE)), apply(MCC,2,function(x)mean(x,na.rm=TRUE)))
    res <- as.data.frame(rbind(overall, bySamples))
    names(res) <-c( "Precision","Recall","Accuracy","F1","MCC")
    return(list(summary=res, precision=precision, recall=recall, accuracy=accuracy, F1=F1, MCC=MCC))
}
##
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
## generare histogram for estimation and quantitative explanatory variables, pairwise correlation matrix
plotRes <- function(res, file="", sampleN = 20000)
{
    pdf(file)
    par(mfrow=c(2,3))
    hist(res$efficiency)
    hist(res$discrepancy)
    hist(res$accuracy)
    hist(res$geneLenM, main="true_length")
    hist(res$percentageEffectM, main="ratio effective_length/true_length")
    hist(res$expression.log2)
    finite <- which(is.finite(res$discrepancy)&is.finite(res$efficiency))
    ss<-sample(finite,sampleN,replace=FALSE)
    par(mfrow=c(1,1))
    pairs( res[ss,c("efficiency", "discrepancy", "accuracy", "mcc", "geneLenM", "percentageEffectM", "expression", "expression.log2")], lower.panel=panel.smooth, upper.panel=panel.cor, pch='.', col = rgb(0, 0, 0, 0.05) )
    # exp.log2
    ss<-sample(1:nrow(networks$ADs)*33,sampleN*2,replace=FALSE)
    obs = log2(as.numeric(as.matrix(networks$ADs))[ss]+1)
    exp = log2(as.numeric(as.matrix(networks$A2D5))[ss]+1)
    r<-cor(exp, obs)
    plot(exp,obs,pch=".", col = rgb(0, 0, 0, 0.05), main=paste0("Exp vs Obs (log2 expression), r=",round(r,3)))
    lines(stats::lowess(exp, obs),  col = "red")
    abline(lm(obs~exp),  col = "blue")
    dev.off()
}

# generate scatter plots for At vs Dt
plotRes.homoeo<-function(resH, file="", sampleN = 20000)
{
    pdf(file)
    ss<-sample(1:nrow(A2.Total)*ncol(A2.Total),sampleN,replace=FALSE)
    A = log2(as.numeric(as.matrix(A2.Total))[ss]+1)
    D = log2(as.numeric(as.matrix(D5.Total))[ss]+1)
    r<-cor(A,D)
    plot(A,D,pch=".", col = rgb(0, 0, 0, 0.05), main=paste0("Log2 Expression: diploid A2 vs D5, r=",round(r,3)))
    lines(stats::lowess(A,D),  col = "red")
    abline(lm(D~A),  col = "blue")
    
    A = log2(as.numeric(as.matrix(ADs.At))[ss]+1)
    D = log2(as.numeric(as.matrix(ADs.Dt))[ss]+1)
    r<-cor(A,D)
    plot(A,D,pch=".", col = rgb(0, 0, 0, 0.05), main=paste0("Log2 Expression: polyploid At vs Dt, r=",round(r,3)))
    lines(stats::lowess(A,D),  col = "red")
    abline(lm(D~A),  col = "blue")
    
    finite <- which(is.finite(resH$efficiencyA)&is.finite(resH$efficiencyD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$efficiencyA[only],resH$efficiencyD[only])
    plot(resH$efficiencyA[only],resH$efficiencyD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Efficiency: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$efficiencyD[only], resH$efficiencyA[only]),  col = "red")
    abline(lm(resH$efficiencyD[only]~resH$efficiencyA[only]),  col = "blue")
    
    finite <- which(is.finite(resH$discrepancyA)&is.finite(resH$discrepancyD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$discrepancyA[only],resH$discrepancyD[only])
    plot(resH$discrepancyA[only],resH$discrepancyD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("discrepancy: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$discrepancyD[only], resH$discrepancyA[only]),  col = "red")
    abline(lm(resH$discrepancyD[only]~resH$discrepancyA[only]),  col = "blue")
    
    finite <- which(is.finite(resH$precisionA)&is.finite(resH$precisionD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$precisionA[only],resH$precisionD[only])
    plot(resH$precisionA[only],resH$precisionD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Precision: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$precisionD[only], resH$precisionA[only]),  col = "red")
    abline(lm(resH$precisionD[only]~resH$precisionA[only]),  col = "blue")
    
    finite <- which(is.finite(resH$recallA)&is.finite(resH$recallD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$recallA[only],resH$recallD[only])
    plot(resH$recallA[only],resH$recallD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Recall: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$recallD[only], resH$recallA[only]),  col = "red")
    abline(lm(resH$recallD[only]~resH$recallA[only]),  col = "blue")

    finite <- which(is.finite(resH$F1A)&is.finite(resH$F1D))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$F1A[only],resH$F1D[only])
    plot(resH$F1A[only],resH$F1D[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("F1 score: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$F1D[only], resH$F1A[only]),  col = "red")
    abline(lm(resH$F1D[only]~resH$F1A[only]),  col = "blue")
    

dev.off() }


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


---------book
#######################
#### GSNAP-polycat ####
#######################

lnames <- load('R-01-polycatDatasets.RData')
lnames
#  "A2.Total"  "D5.Total"  "ADs.Total" "A2.At"     "D5.At"     "ADs.At"
# "A2.Dt"     "D5.Dt"     "ADs.Dt"    "ADs.AtN"   "ADs.DtN"
lnames <- load('R-01-polycatNetworkDatasets.RData')
lnames
# "coldata"       "networks"      "networks.rld"  "geneLen"       "networks.rpkm"

# see if meets expectation
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
# not perfect 0, but close enough

## calculate assignment efficiency - percentage of total reads assigned to At and Dt homoeologs
total <- D5.Total + A2.Total
unique(total - ADs.Total)
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
system("cut -f3,4 eflen.snp.rl300.txt|paste eflen.snp.rl100.txt - >eflenList.snp.txt")
geneL <- read.table("eflenList.snp.txt",header=TRUE,sep="\t")
len<- geneL$trueLen
names(len)<-geneL$id
# length for 37223 genes
len<-len[rownames(A2.Total)]

# PercentageEffectM
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
polycat<-res
plotRes(polycat, file="s2.eval.polycat.pdf")

polycatH<-resH
plotRes.homoeo(polycatH, file="s2.eval.polycat.homoeolog.pdf")

library(gplots)
pdf("s2.eval.polycat.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()

polycatSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)

save(polycat, polycatH, polycatSummary, file="s2.assign_eval.polycat.Rdata")
# save(efficency, efficencyA, efficencyD, accuracy, accuracyA, accuracyD, discrepancy, discrepancyA, discrepancyD, info, geneL, ratio100, ratio300, file="s2.polycat.estimation.rdata")


################
#### bowtie ####
################

lnames <- load('R-01-bowtieDatasets.RData')
lnames
#  "A2.Total" "A2.At"    "A2.Dt"    "D5.Total" "D5.At"    "D5.Dt"    "ADs.At"  "ADs.Dt"
lnames <- load('R-01-bowtieNetworkDatasets.RData')
lnames
# "coldata"      "networks"     "networks.rld" "networks.tpm"

# see if meets expectation
quantile(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
quantile(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
# 0, perfect

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

## prepare explainatory variables for analysis
info<-coldata
info$type <- c(rep("PE",22), rep("SE",11))
info$ADsLibSize <- colSums(ADs.Total)

# geneLenM, output from "detectEffectiveRegion.r"
geneL <- read.table("eflenList.txt",header=TRUE,sep="\t")
len<- geneL$true
names(len)<-geneL$gene
# length for 37223 genes
len<-len[gsub("a$","",rownames(A2.Total))]

#\ percentageEffectM
ratio100<- geneL$theoretical100 / geneL$true
names(ratio100)<-geneL$gene
ratio300<- geneL$theoretical300 / geneL$true
names(ratio300)<-geneL$gene
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
len<- geneL$trueLen
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





#######################
#### Bowtie2-rsem ####
#######################

lnames <- load('R-01-rsemDatasets.RData')
lnames
#  "A2.Total" "A2.At"    "A2.Dt"    "D5.Total" "D5.At"    "D5.Dt"    "ADs.At"  "ADs.Dt"
lnames <- load('R-01-rsemNetworkDatasets.RData')
lnames
# "coldata"      "networks"     "networks.rld" "networks.rpkm"

# see if meets expectation
ADs.Total = ADs.At + ADs.Dt
quantile(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
quantile(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
# 0, perfect

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
geneL <- read.table("eflenList.txt",header=TRUE,sep="\t")
len<- geneL$true
names(len)<-geneL$gene
# length for 37223 genes
len<-len[gsub("a$","",rownames(A2.Total))]

#\ percentageEffectM
ratio100<- geneL$theoretical100 / geneL$true
names(ratio100)<-geneL$gene
ratio300<- geneL$theoretical300 / geneL$true
names(ratio300)<-geneL$gene
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
rsem<-res
plotRes(rsem, file="s2.eval.rsem.pdf")

rsemH<-resH
plotRes.homoeo(rsemH, file="s2.eval.rsem.homoeolog.pdf")

library(gplots)
pdf("s2.eval.rsem.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()


rsemSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)

save(rsem, rsemH, rsemSummary, file="s2.assign_eval.rsem.Rdata")


########################
#### Bowtie2-HyLite ####
########################

lnames <- load('R-01-hyliteDatasets.true.RData')
lnames
#  "A2.Total"  "D5.Total"  "ADs.Total" "ADs.At"    "ADs.Dt"    "At.T" "Dt.T"      "At.F"      "Dt.F"
lnames <- load('R-01-hyliteNetworkDatasets.RData')
lnames
# "coldata"       "networks"      "networks.rld"  "geneLen"       "networks.rpkm"

# see if meets expectation
quantile(as.numeric(as.matrix((ADs.Total) - (ADs.At + ADs.Dt + ADs.N) ) ) )
quantile(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
# 0, perfect
unique(as.numeric(as.matrix((ADs.At + ADs.Dt) - (At.T +At.F + Dt.T + Dt.F) ) ) )

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
resBinA=binC(TP = At.T, TN = Dt.T, FP = At.F, FN = Dt.F)
# D as TRUE
resBinD=binC(TP = Dt.T, TN = At.T, FP = Dt.F ,FN = At.F)
# summary
resBinA$summary
resBinD$summary


###############

## prepare explainatory variables for analysis
info<-coldata
info$type <- c(rep("PE",22), rep("SE",11))
info$ADsLibSize <- colSums(ADs.Total)

# geneLenM, output from "detectEffectiveRegion.r"
geneL <- read.table("eflenList.txt",header=TRUE,sep="\t")
len<- geneL$true
names(len)<-geneL$gene
# length for 37223 genes
len<-len[gsub("a$","",rownames(A2.Total))]

#\ percentageEffectM
ratio100<- geneL$theoretical100 / geneL$true
names(ratio100)<-geneL$gene
ratio300<- geneL$theoretical300 / geneL$true
names(ratio300)<-geneL$gene
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
hylite<-res
plotRes(hylite, file="s2.eval.hylite.pdf")

hyliteH<-resH
plotRes.homoeo(hyliteH, file="s2.eval.hylite.homoeolog.pdf")

library(gplots)
pdf("s2.eval.hylite.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()

hyliteSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)

save(hylite, hyliteH, hyliteSummary, file="s2.assign_eval.hylite.Rdata")

