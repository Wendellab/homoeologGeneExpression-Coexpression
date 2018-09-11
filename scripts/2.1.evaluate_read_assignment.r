# generare histogram for estimation and quantitative explanatory variables, pairwise correlation matrix
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
    pairs( res[ss,c("efficiency", "discrepancy", "accuracy", "precision", "geneLenM", "percentageEffectM", "expression", "expression.log2")], lower.panel=panel.smooth, upper.panel=panel.cor, pch='.', col = rgb(0, 0, 0, 0.05) )
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
    
    finite <- which(is.finite(resH$accuracyA)&is.finite(resH$accuracyD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$accuracyA[only],resH$accuracyD[only])
    plot(resH$accuracyA[only],resH$accuracyD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Accuracy (or Recall): At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$accuracyD[only], resH$accuracyA[only]),  col = "red")
    abline(lm(resH$accuracyD[only]~resH$accuracyA[only]),  col = "blue")
    
    finite <- which(is.finite(resH$precisionA)&is.finite(resH$precisionD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$precisionA[only],resH$precisionD[only])
    plot(resH$precisionA[only],resH$precisionD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Precision: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$precisionD[only], resH$precisionA[only]),  col = "red")
    abline(lm(resH$precisionD[only]~resH$precisionA[only]),  col = "blue")

    finite <- which(is.finite(resH$FstatA)&is.finite(resH$FstatD))
    only<-sample(finite, sampleN, replace=FALSE)
    r<-cor(resH$FstatA[only],resH$FstatD[only])
    plot(resH$FstatA[only],resH$FstatD[only],pch=".",xlim=c(0,1), ylim=c(0,1), col = rgb(0, 0, 0, 0.05), main=paste0("Fstat: At vs Dt, r=",round(r,3)))
    lines(stats::lowess(resH$FstatD[only], resH$FstatA[only]),  col = "red")
    abline(lm(resH$FstatD[only]~resH$FstatA[only]),  col = "blue")
    

dev.off()
}



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



#######################
#### GSNAP-polycat ####
#######################

lnames <- load('R-01-polycatDatasets.RData')
lnames
#  "A2.Total"  "D5.Total"  "ADs.Total" "A2.At"     "D5.At"     "ADs.At"    "A2.Dt"     "D5.Dt"     "ADs.Dt"
# "ADs.AtN"   "ADs.DtN"   "A2.N"      "D5.N"      "ADs.N"     "A2.AD"     "D5.AD"     "ADs.AD"
lnames <- load('R-01-polycatNetworkDatasets.RData')
lnames
# "coldata"       "networks"      "networks.rld"  "geneLen"       "networks.rpkm"

# see if meets expectation
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.Total) - (A2.Total + D5.Total) ) ) )
# 0, perfect

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


## calculate assignment accuracy - percentage of correctly partitioned reads among assigned reads
# mis-assigned reads cannot be directly counted from ADs reads, so using mis-assigned from independent A2 and D5 mapping to approximete.
# This approximation can be validated by comparing assigned reads between observation (ADs.At + ADs.Dt) and expectation (A2.At + A2.Dt + D5.At + D5.Dt), as:
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
unique(as.numeric(as.matrix((ADs.At) - (A2.At + D5.At) ) ) )
# given the equal amounts of assigned reads between ADs and A2+D5, we accecpt that polycat partitioning worked in a consistent manner, that mis-assigned reads during A2+D5 will stay mis-assigned in ADs
true <- D5.Dt + A2.At
trueA <- A2.At
trueD <- D5.Dt

diploid <- A2.At + A2.Dt + D5.Dt + D5.At
diploidA <- A2.At + A2.Dt
diploidD <- D5.Dt + D5.At

accuracy <- true/diploid
accuracyA <- trueA/diploidA
accuracyD <- trueD/diploidD

quantile(as.numeric(as.matrix(accuracy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

# summary table
overall <- c(sum(true)/sum(diploid), sum(trueA)/sum(diploidA), sum(trueD)/sum(diploidD))
bySamples<-cbind(apply(true,2,sum)/apply(diploid,2,sum),apply(trueA,2,sum)/apply(diploidA,2,sum),apply(trueD,2,sum)/apply(diploidD,2,sum))
resAccu <- as.data.frame(rbind(overall, bySamples))
names(resAccu) <-c("total","At","Dt")
resAccu

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

## Additionally from the perspective of Precision and recall, we obtained TP, TN, FP, FN from the diploid datasets for calculation. The recall is equivalent to previous Accuracy calculated for homoeologs.
recallA <- A2.At/(A2.At + A2.Dt)
recallD <- D5.Dt/(D5.Dt + D5.At)
recall<- accuracy
# precision
precisionA <- A2.At/(A2.At + D5.At)
precisionD <- D5.Dt/(D5.Dt + A2.Dt)
# or
positive <- A2.At + D5.At + D5.Dt + A2.Dt  #same as diploid
positiveA <- A2.At + D5.At
positiveD <- D5.Dt + A2.Dt
precision <- true/positive
precisionA <- trueA/positiveA
precisionD <- trueD/positiveD
# summary table
overall <- c(sum(true)/sum(positive), sum(trueA)/sum(positiveA), sum(trueD)/sum(positiveD))
bySamples<-cbind(apply(true,2,sum)/apply(positive,2,sum),apply(trueA,2,sum)/apply(positiveA,2,sum),apply(trueD,2,sum)/apply(positiveD,2,sum))
resPrec <- as.data.frame(rbind(overall, bySamples))
names(resPrec) <-c("total","At","Dt")
resPrec
# same values for paired reads between recall and precision, all as Accuracy
resPrec$total - resAccu$total
# F= 2x(precision x recall)/(precision + recall)
FstatA <- 2*(recallA*precisionA)/(recallA+precisionA)
FstatD <- 2*(recallD*precisionD)/(recallD+precisionD)
Fstat <- 2*(recall*precision)/(recall+precision)  # same as recall or precision or accuracy, not really meaningful

# summary table
overall <- c(mean(as.numeric(as.matrix(FstatA)), na.rm=TRUE), mean(as.numeric(as.matrix(FstatD)), na.rm=TRUE))
bySamples<-cbind( apply(FstatA,2,function(x)mean(x, na.rm=TRUE)),apply(FstatD,2,function(x)mean(x, na.rm=TRUE)))
resFstat <- as.data.frame(rbind(overall, bySamples))
names(resFstat) <-c( "At","Dt")
resFstat

## MCC
TP = A2.At
TN = D5.Dt
FP = D5.At
FN = A2.Dt
MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
TN = D5.Dt + D5.N
FN =  A2.Dt + A2.N
MCCAn = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
TP = D5.Dt
TN = A2.At + A2.N
FP = A2.Dt
FN = D5.At + D5.N
MCCDn = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

overall <- c(mean(as.numeric(as.matrix(MCC)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCAn)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCDn)), na.rm=TRUE) )
bySamples<-cbind( apply(MCC,2,function(x)mean(x, na.rm=TRUE)),apply(MCCAn,2,function(x)mean(x, na.rm=TRUE)),apply(MCCDn,2,function(x)mean(x, na.rm=TRUE)))
resMCC <- as.data.frame(rbind(overall, bySamples))
names(resMCC) <-c( "total", "At","Dt")
resMCC

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
len<-len[rownames(A2.Total)]

# percentageEffectM
ratio100<- geneL$theoretical100 / geneL$true
names(ratio100)<-geneL$gene
ratio300<- geneL$theoretical300 / geneL$true
names(ratio300)<-geneL$gene
# effective_length/true_length ratio for 37223 genes
ratio100<-ratio100[rownames(A2.Total)]
ratio300<-ratio300[rownames(A2.Total)]
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
percentageEffectM <- c( rep(ratio100,11),rep(ratio300,22) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(accuracy)), precision=as.numeric(as.matrix(precision)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)), as.numeric(as.matrix(accuracyA)), as.numeric(as.matrix(accuracyD)), as.numeric(as.matrix(precisionA)), as.numeric(as.matrix(precisionD)), as.numeric(as.matrix(FstatA)), as.numeric(as.matrix(FstatD)),as.numeric(as.matrix(MCC)), as.numeric(as.matrix(MCCAn)), as.numeric(as.matrix(MCCDn)) )
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","accuracyA", "accuracyD", "precisionA", "precisionD", "FstatA", "FstatD","MCC","MCCAn","MCCDn")

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
textplot(round(resAccu*100,2))
mtext("Accuracy (%)")
textplot(round(resPrec*100,2))
mtext("Precision (%)")
textplot(round(resFstat*100,2))
mtext("F measure (%)")
textplot(round(resMCC*100,2))
mtext("MCC (%)")
dev.off()

polycatSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, Accuracy=resAccu, Precision=resPrec, Fmeasure=resFstat, MCC= resMCC)

save(polycat, polycatH, polycatSummary, file="s2.assign_eval.polycat.Rdata")
# save(efficency, efficencyA, efficencyD, accuracy, accuracyA, accuracyD, discrepancy, discrepancyA, discrepancyD, info, geneL, ratio100, ratio300, file="s2.polycat.estimation.rdata")


################
#### Salmon ####
################

lnames <- load('R-01-salmonDatasets.RData')
lnames
#  "A2.Total" "A2.At"    "A2.Dt"    "D5.Total" "D5.At"    "D5.Dt"    "ADs.At"  "ADs.Dt"
lnames <- load('R-01-salmonNetworkDatasets.RData')
lnames
# "coldata"      "networks"     "networks.rld" "networks.tpm"

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


## calculate assignment accuracy - percentage of correctly partitioned reads among assigned reads
# mis-assigned reads cannot be directly counted from ADs reads, so using mis-assigned from independent A2 and D5 mapping to approximete.
# This approximation can be validated by comparing assigned reads between observation (ADs.At + ADs.Dt) and expectation (A2.At + A2.Dt + D5.At + D5.Dt), as:
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At) - (A2.At + D5.At) ) ) )
# given the equal amounts of assigned reads between ADs and A2+D5, we accecpt that salmon partitioning worked in a consistent manner, that mis-assigned reads during A2+D5 will stay mis-assigned in ADs
true <- D5.Dt + A2.At
trueA <- A2.At
trueD <- D5.Dt

diploid <- A2.At + A2.Dt + D5.Dt + D5.At
diploidA <- A2.At + A2.Dt
diploidD <- D5.Dt + D5.At

accuracy <- true/diploid
accuracyA <- trueA/diploidA
accuracyD <- trueD/diploidD

quantile(as.numeric(as.matrix(accuracy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

# summary table
overall <- c(sum(true)/sum(diploid), sum(trueA)/sum(diploidA), sum(trueD)/sum(diploidD))
bySamples<-cbind(apply(true,2,sum)/apply(diploid,2,sum),apply(trueA,2,sum)/apply(diploidA,2,sum),apply(trueD,2,sum)/apply(diploidD,2,sum))
resAccu <- as.data.frame(rbind(overall, bySamples))
names(resAccu) <-c("total","At","Dt")
resAccu

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

## Additionally from the perspective of Precision and recall, we obtained TP, TN, FP, FN from the diploid datasets for calculation. The recall is equivalent to previous Accuracy calculated for homoeologs.
recallA <- A2.At/(A2.At + A2.Dt)
recallD <- D5.Dt/(D5.Dt + D5.At)
recall<- accuracy
# precision
precisionA <- A2.At/(A2.At + D5.At)
precisionD <- D5.Dt/(D5.Dt + A2.Dt)
# or
positive <- A2.At + D5.At + D5.Dt + A2.Dt  #same as diploid
positiveA <- A2.At + D5.At
positiveD <- D5.Dt + A2.Dt
precision <- true/positive
precisionA <- trueA/positiveA
precisionD <- trueD/positiveD
# summary table
overall <- c(sum(true)/sum(positive), sum(trueA)/sum(positiveA), sum(trueD)/sum(positiveD))
bySamples<-cbind(apply(true,2,sum)/apply(positive,2,sum),apply(trueA,2,sum)/apply(positiveA,2,sum),apply(trueD,2,sum)/apply(positiveD,2,sum))
resPrec <- as.data.frame(rbind(overall, bySamples))
names(resPrec) <-c("total","At","Dt")
resPrec
# same values for paired reads between recall and precision, all as Accuracy
resPrec$total - resAccu$total
# F= 2x(precision x recall)/(precision + recall)
FstatA <- 2*(recallA*precisionA)/(recallA+precisionA)
FstatD <- 2*(recallD*precisionD)/(recallD+precisionD)
Fstat <- 2*(recall*precision)/(recall+precision)  # same as recall or precision or accuracy, not really meaningful

# summary table
overall <- c(mean(as.numeric(as.matrix(FstatA)), na.rm=TRUE), mean(as.numeric(as.matrix(FstatD)), na.rm=TRUE))
bySamples<-cbind( apply(FstatA,2,function(x)mean(x, na.rm=TRUE)),apply(FstatD,2,function(x)mean(x, na.rm=TRUE)))
resFstat <- as.data.frame(rbind(overall, bySamples))
names(resFstat) <-c( "At","Dt")
resFstat

## MCC
TP = A2.At
TN = D5.Dt
FP = D5.At
FN = A2.Dt
MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCCAn = MCC
MCCDn = MCC

overall <- c(mean(as.numeric(as.matrix(MCC)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCAn)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCDn)), na.rm=TRUE) )
bySamples<-cbind( apply(MCC,2,function(x)mean(x, na.rm=TRUE)),apply(MCCAn,2,function(x)mean(x, na.rm=TRUE)),apply(MCCDn,2,function(x)mean(x, na.rm=TRUE)))
resMCC <- as.data.frame(rbind(overall, bySamples))
names(resMCC) <-c( "total", "At","Dt")
resMCC


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
percentageEffectM <- c( rep(ratio100,11),rep(ratio300,22) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(accuracy)), precision=as.numeric(as.matrix(precision)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)), as.numeric(as.matrix(accuracyA)), as.numeric(as.matrix(accuracyD)), as.numeric(as.matrix(precisionA)), as.numeric(as.matrix(precisionD)), as.numeric(as.matrix(FstatA)), as.numeric(as.matrix(FstatD)),as.numeric(as.matrix(MCC)))
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","accuracyA", "accuracyD", "precisionA", "precisionD", "FstatA", "FstatD","MCC")

## plot and save results
salmon<-res
plotRes(salmon, file="s2.eval.salmon.pdf")

salmonH<-resH
plotRes.homoeo(salmonH, file="s2.eval.salmon.homoeolog.pdf")

library(gplots)
pdf("s2.eval.salmon.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resAccu*100,2))
mtext("Accuracy (%)")
textplot(round(resPrec*100,2))
mtext("Precision (%)")
textplot(round(resFstat*100,2))
mtext("F measure (%)")
textplot(round(resMCC*100,2))
mtext("MCC (%)")

dev.off()

salmonSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, Accuracy=resAccu, Precision=resPrec, Fmeasure=resFstat, MCC=resMCC)

save(salmon, salmonH, salmonSummary, file="s2.assign_eval.salmon.Rdata")

################
#### kallisto ####
################

lnames <- load('R-01-kallistoDatasets.RData')
lnames
#  "A2.Total" "A2.At"    "A2.Dt"    "D5.Total" "D5.At"    "D5.Dt"    "ADs.At"  "ADs.Dt"
lnames <- load('R-01-kallistoNetworkDatasets.RData')
lnames
# "coldata"      "networks"     "networks.rld" "networks.tpm"

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


## calculate assignment accuracy - percentage of correctly partitioned reads among assigned reads
# mis-assigned reads cannot be directly counted from ADs reads, so using mis-assigned from independent A2 and D5 mapping to approximete.
# This approximation can be validated by comparing assigned reads between observation (ADs.At + ADs.Dt) and expectation (A2.At + A2.Dt + D5.At + D5.Dt), as:
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At) - (A2.At + D5.At) ) ) )
# given the equal amounts of assigned reads between ADs and A2+D5, we accecpt that kallisto partitioning worked in a consistent manner, that mis-assigned reads during A2+D5 will stay mis-assigned in ADs
true <- D5.Dt + A2.At
trueA <- A2.At
trueD <- D5.Dt

diploid <- A2.At + A2.Dt + D5.Dt + D5.At
diploidA <- A2.At + A2.Dt
diploidD <- D5.Dt + D5.At

accuracy <- true/diploid
accuracyA <- trueA/diploidA
accuracyD <- trueD/diploidD

quantile(as.numeric(as.matrix(accuracy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

# summary table
overall <- c(sum(true)/sum(diploid), sum(trueA)/sum(diploidA), sum(trueD)/sum(diploidD))
bySamples<-cbind(apply(true,2,sum)/apply(diploid,2,sum),apply(trueA,2,sum)/apply(diploidA,2,sum),apply(trueD,2,sum)/apply(diploidD,2,sum))
resAccu <- as.data.frame(rbind(overall, bySamples))
names(resAccu) <-c("total","At","Dt")
resAccu

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

## Additionally from the perspective of Precision and recall, we obtained TP, TN, FP, FN from the diploid datasets for calculation. The recall is equivalent to previous Accuracy calculated for homoeologs.
recallA <- A2.At/(A2.At + A2.Dt)
recallD <- D5.Dt/(D5.Dt + D5.At)
recall<- accuracy
# precision
precisionA <- A2.At/(A2.At + D5.At)
precisionD <- D5.Dt/(D5.Dt + A2.Dt)
# or
positive <- A2.At + D5.At + D5.Dt + A2.Dt  #same as diploid
positiveA <- A2.At + D5.At
positiveD <- D5.Dt + A2.Dt
precision <- true/positive
precisionA <- trueA/positiveA
precisionD <- trueD/positiveD
# summary table
overall <- c(sum(true)/sum(positive), sum(trueA)/sum(positiveA), sum(trueD)/sum(positiveD))
bySamples<-cbind(apply(true,2,sum)/apply(positive,2,sum),apply(trueA,2,sum)/apply(positiveA,2,sum),apply(trueD,2,sum)/apply(positiveD,2,sum))
resPrec <- as.data.frame(rbind(overall, bySamples))
names(resPrec) <-c("total","At","Dt")
resPrec
# same values for paired reads between recall and precision, all as Accuracy
resPrec$total - resAccu$total
# F= 2x(precision x recall)/(precision + recall)
FstatA <- 2*(recallA*precisionA)/(recallA+precisionA)
FstatD <- 2*(recallD*precisionD)/(recallD+precisionD)
Fstat <- 2*(recall*precision)/(recall+precision)  # same as recall or precision or accuracy, not really meaningful

# summary table
overall <- c(mean(as.numeric(as.matrix(FstatA)), na.rm=TRUE), mean(as.numeric(as.matrix(FstatD)), na.rm=TRUE))
bySamples<-cbind( apply(FstatA,2,function(x)mean(x, na.rm=TRUE)),apply(FstatD,2,function(x)mean(x, na.rm=TRUE)))
resFstat <- as.data.frame(rbind(overall, bySamples))
names(resFstat) <-c( "At","Dt")
resFstat
## MCC
TP = A2.At
TN = D5.Dt
FP = D5.At
FN = A2.Dt
MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCCAn = MCC
MCCDn = MCC

overall <- c(mean(as.numeric(as.matrix(MCC)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCAn)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCDn)), na.rm=TRUE) )
bySamples<-cbind( apply(MCC,2,function(x)mean(x, na.rm=TRUE)),apply(MCCAn,2,function(x)mean(x, na.rm=TRUE)),apply(MCCDn,2,function(x)mean(x, na.rm=TRUE)))
resMCC <- as.data.frame(rbind(overall, bySamples))
names(resMCC) <-c( "total", "At","Dt")
resMCC


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
percentageEffectM <- c( rep(ratio100,11),rep(ratio300,22) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(accuracy)), precision=as.numeric(as.matrix(precision)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)), as.numeric(as.matrix(accuracyA)), as.numeric(as.matrix(accuracyD)), as.numeric(as.matrix(precisionA)), as.numeric(as.matrix(precisionD)), as.numeric(as.matrix(FstatA)), as.numeric(as.matrix(FstatD)), as.numeric(as.matrix(MCC)))
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","accuracyA", "accuracyD", "precisionA", "precisionD", "FstatA", "FstatD", "MCC")

## plot and save results
kallisto<-res
plotRes(kallisto, file="s2.eval.kallisto.pdf")

kallistoH<-resH
plotRes.homoeo(kallistoH, file="s2.eval.kallisto.homoeolog.pdf")

library(gplots)
pdf("s2.eval.kallisto.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resAccu*100,2))
mtext("Accuracy (%)")
textplot(round(resPrec*100,2))
mtext("Precision (%)")
textplot(round(resFstat*100,2))
mtext("F measure (%)")
textplot(round(resMCC*100,2))
mtext("MCC (%)")

dev.off()

kallistoSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, Accuracy=resAccu, Precision=resPrec, Fmeasure=resFstat, MCC=resMCC)

save(kallisto, kallistoH, kallistoSummary, file="s2.assign_eval.kallisto.Rdata")

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


## calculate assignment accuracy - percentage of correctly partitioned reads among assigned reads
# mis-assigned reads cannot be directly counted from ADs reads, so using mis-assigned from independent A2 and D5 mapping to approximete.
# This approximation can be validated by comparing assigned reads between observation (ADs.At + ADs.Dt) and expectation (A2.At + A2.Dt + D5.At + D5.Dt), as:
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt))/(ADs.At + ADs.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At + ADs.Dt) - (A2.At + A2.Dt + D5.At + D5.Dt) ) ) )
summary(as.numeric(as.matrix((ADs.At) - (A2.At + D5.At) ) ) )
# given the equal amounts of assigned reads between ADs and A2+D5, we accecpt that rsem partitioning worked in a consistent manner, that mis-assigned reads during A2+D5 will stay mis-assigned in ADs
true <- D5.Dt + A2.At
trueA <- A2.At
trueD <- D5.Dt

diploid <- A2.At + A2.Dt + D5.Dt + D5.At
diploidA <- A2.At + A2.Dt
diploidD <- D5.Dt + D5.At

accuracy <- true/diploid
accuracyA <- trueA/diploidA
accuracyD <- trueD/diploidD

quantile(as.numeric(as.matrix(accuracy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

# summary table
overall <- c(sum(true)/sum(diploid), sum(trueA)/sum(diploidA), sum(trueD)/sum(diploidD))
bySamples<-cbind(apply(true,2,sum)/apply(diploid,2,sum),apply(trueA,2,sum)/apply(diploidA,2,sum),apply(trueD,2,sum)/apply(diploidD,2,sum))
resAccu <- as.data.frame(rbind(overall, bySamples))
names(resAccu) <-c("total","At","Dt")
resAccu

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

## Additionally from the perspective of Precision and recall, we obtained TP, TN, FP, FN from the diploid datasets for calculation. The recall is equivalent to previous Accuracy calculated for homoeologs.
recallA <- A2.At/(A2.At + A2.Dt)
recallD <- D5.Dt/(D5.Dt + D5.At)
recall<- accuracy
# precision
precisionA <- A2.At/(A2.At + D5.At)
precisionD <- D5.Dt/(D5.Dt + A2.Dt)
# or
positive <- A2.At + D5.At + D5.Dt + A2.Dt  #same as diploid
positiveA <- A2.At + D5.At
positiveD <- D5.Dt + A2.Dt
precision <- true/positive
precisionA <- trueA/positiveA
precisionD <- trueD/positiveD
# summary table
overall <- c(sum(true)/sum(positive), sum(trueA)/sum(positiveA), sum(trueD)/sum(positiveD))
bySamples<-cbind(apply(true,2,sum)/apply(positive,2,sum),apply(trueA,2,sum)/apply(positiveA,2,sum),apply(trueD,2,sum)/apply(positiveD,2,sum))
resPrec <- as.data.frame(rbind(overall, bySamples))
names(resPrec) <-c("total","At","Dt")
resPrec
# same values for paired reads between recall and precision, all as Accuracy
resPrec$total - resAccu$total
# F= 2x(precision x recall)/(precision + recall)
FstatA <- 2*(recallA*precisionA)/(recallA+precisionA)
FstatD <- 2*(recallD*precisionD)/(recallD+precisionD)
Fstat <- 2*(recall*precision)/(recall+precision)  # same as recall or precision or accuracy, not really meaningful

# summary table
overall <- c(mean(as.numeric(as.matrix(FstatA)), na.rm=TRUE), mean(as.numeric(as.matrix(FstatD)), na.rm=TRUE))
bySamples<-cbind( apply(FstatA,2,function(x)mean(x, na.rm=TRUE)),apply(FstatD,2,function(x)mean(x, na.rm=TRUE)))
resFstat <- as.data.frame(rbind(overall, bySamples))
names(resFstat) <-c( "At","Dt")
resFstat

## MCC
TP = A2.At
TN = D5.Dt
FP = D5.At
FN = A2.Dt
MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
MCCAn = MCC
MCCDn = MCC

overall <- c(mean(as.numeric(as.matrix(MCC)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCAn)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCDn)), na.rm=TRUE) )
bySamples<-cbind( apply(MCC,2,function(x)mean(x, na.rm=TRUE)),apply(MCCAn,2,function(x)mean(x, na.rm=TRUE)),apply(MCCDn,2,function(x)mean(x, na.rm=TRUE)))
resMCC <- as.data.frame(rbind(overall, bySamples))
names(resMCC) <-c( "total", "At","Dt")
resMCC


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
percentageEffectM <- c( rep(ratio100,11),rep(ratio300,22) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(accuracy)), precision=as.numeric(as.matrix(precision)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)), as.numeric(as.matrix(accuracyA)), as.numeric(as.matrix(accuracyD)), as.numeric(as.matrix(precisionA)), as.numeric(as.matrix(precisionD)), as.numeric(as.matrix(FstatA)), as.numeric(as.matrix(FstatD)), as.numeric(as.matrix(MCC)))
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","accuracyA", "accuracyD", "precisionA", "precisionD", "FstatA", "FstatD","MCC")

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
textplot(round(resAccu*100,2))
mtext("Accuracy (%)")
textplot(round(resPrec*100,2))
mtext("Precision (%)")
textplot(round(resFstat*100,2))
mtext("F measure (%)")
textplot(round(resMCC*100,2))
mtext("MCC (%)")
dev.off()

rsemSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, Accuracy=resAccu, Precision=resPrec, Fmeasure=resFstat, MCC=resMCC)

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

## calculate assignment accuracy - percentage of correctly partitioned reads among assigned reads
# different from polycat, mis-assigned reads can be directly counted from ADs reads
summary(as.numeric(as.matrix(( (ADs.At + ADs.Dt) - (At.T + At.F + Dt.T + Dt.F))/(ADs.At + ADs.Dt) ) ) )
# given the equal amounts of assigned reads between ADs and A2+D5, we accecpt that hylite partitioning worked in a consistent manner, that mis-assigned reads during A2+D5 will stay mis-assigned in ADs
true <- At.T + Dt.T
trueA <- At.T
trueD <- Dt.T

diploid <- At.T + Dt.T + At.F + Dt.F
diploidA <- At.T + Dt.F
diploidD <- Dt.T + At.F

accuracy <- true/diploid
accuracyA <- trueA/diploidA
accuracyD <- trueD/diploidD

quantile(as.numeric(as.matrix(accuracy)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyA)),na.rm=TRUE) #NaN from 0/0, no way to learn
quantile(as.numeric(as.matrix(accuracyD)),na.rm=TRUE) #NaN from 0/0, no way to learn

# summary table
overall <- c(sum(true)/sum(diploid), sum(trueA)/sum(diploidA), sum(trueD)/sum(diploidD))
bySamples<-cbind(apply(true,2,sum)/apply(diploid,2,sum),apply(trueA,2,sum)/apply(diploidA,2,sum),apply(trueD,2,sum)/apply(diploidD,2,sum))
resAccu <- as.data.frame(rbind(overall, bySamples))
names(resAccu) <-c("total","At","Dt")
resAccu


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

## Additionally from the perspective of Precision and recall, we obtained TP, TN, FP, FN from the diploid datasets for calculation. The recall is equivalent to previous Accuracy calculated for homoeologs.
recallA <- At.T/(At.T + Dt.F)
recallD <- Dt.T/(Dt.T + At.F)
recall<- accuracy
# precision
precisionA <- At.T/(At.T + At.F)
precisionD <- Dt.T/(Dt.T + Dt.F)
# or
positive <- At.T + Dt.T + At.F + Dt.F  #same as diploid
positiveA <- At.T + At.F
positiveD <- Dt.T + Dt.F
precision <- true/positive
precisionA <- trueA/positiveA
precisionD <- trueD/positiveD
# summary table
overall <- c(sum(true)/sum(positive), sum(trueA)/sum(positiveA), sum(trueD)/sum(positiveD))
bySamples<-cbind(apply(true,2,sum)/apply(positive,2,sum),apply(trueA,2,sum)/apply(positiveA,2,sum),apply(trueD,2,sum)/apply(positiveD,2,sum))
resPrec <- as.data.frame(rbind(overall, bySamples))
names(resPrec) <-c("total","At","Dt")
resPrec
# same values for paired reads between recall and precision, all as Accuracy
resPrec$total - resAccu$total
# F= 2x(precision x recall)/(precision + recall)
FstatA <- 2*(recallA*precisionA)/(recallA+precisionA)
FstatD <- 2*(recallD*precisionD)/(recallD+precisionD)
Fstat <- 2*(recall*precision)/(recall+precision)  # same as recall or precision or accuracy, not really meaningful

# summary table
overall <- c(mean(as.numeric(as.matrix(FstatA)), na.rm=TRUE), mean(as.numeric(as.matrix(FstatD)), na.rm=TRUE))
bySamples<-cbind( apply(FstatA,2,function(x)mean(x, na.rm=TRUE)),apply(FstatD,2,function(x)mean(x, na.rm=TRUE)))
resFstat <- as.data.frame(rbind(overall, bySamples))
names(resFstat) <-c( "At","Dt")
resFstat

## MCC
TP = At.T
TN = Dt.T
FP = At.F
FN = Dt.F
MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
D5.N = A2.Total - At.T - Dt.F
A2.N = D5.Total - Dt.T - At.F
TN = Dt.T + D5.N
FN = Dt.F + A2.N
MCCAn = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
TP = Dt.T
TN = At.T + A2.N
FP = At.F
FN = Dt.F + D5.N
MCCDn = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

overall <- c(mean(as.numeric(as.matrix(MCC)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCAn)), na.rm=TRUE), mean(as.numeric(as.matrix(MCCDn)), na.rm=TRUE) )
bySamples<-cbind( apply(MCC,2,function(x)mean(x, na.rm=TRUE)),apply(MCCAn,2,function(x)mean(x, na.rm=TRUE)),apply(MCCDn,2,function(x)mean(x, na.rm=TRUE)))
resMCC <- as.data.frame(rbind(overall, bySamples))
names(resMCC) <-c( "total", "At","Dt")
resMCC


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
percentageEffectM <- c( rep(ratio100,11),rep(ratio300,22) )

## make working dataset
res <- data.frame(gene, sample, tissue, rep, geneLenM, percentageEffectM, expression=as.numeric(as.matrix(expression)), expression.log2=as.numeric(as.matrix(expression.log2)), efficiency=as.numeric(as.matrix(efficiency)), discrepancy=as.numeric(as.matrix(discrepancy)), accuracy=as.numeric(as.matrix(accuracy)), precision=as.numeric(as.matrix(precision)))
# discrepancy very tricky:
# res$discrepancy0 <- res$discrepancy
# res$discrepancy0[discrepancy>30] <- NA
resH <- data.frame(as.numeric(as.matrix(efficiencyA)), as.numeric(as.matrix(efficiencyD)), as.numeric(as.matrix(discrepancyA)), as.numeric(as.matrix(discrepancyD)), as.numeric(as.matrix(accuracyA)), as.numeric(as.matrix(accuracyD)), as.numeric(as.matrix(precisionA)), as.numeric(as.matrix(precisionD)), as.numeric(as.matrix(FstatA)), as.numeric(as.matrix(FstatD)), as.numeric(as.matrix(MCC)), as.numeric(as.matrix(MCCAn)), as.numeric(as.matrix(MCCDn)))
names(resH) <- c("efficiencyA", "efficiencyD", "discrepancyA", "discrepancyD","accuracyA", "accuracyD", "precisionA", "precisionD", "FstatA", "FstatD","MCC","MCCAn","MCCDn")

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
textplot(round(resAccu*100,2))
mtext("Accuracy (%)")
textplot(round(resPrec*100,2))
mtext("Precision (%)")
textplot(round(resFstat*100,2))
mtext("F measure (%)")
textplot(round(resMCC*100,2))
mtext("MCC (%)")
dev.off()

hyliteSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, Accuracy=resAccu, Precision=resPrec, Fmeasure=resFstat, MCC=resMCC)

save(hylite, hyliteH, hyliteSummary, file="s2.assign_eval.hylite.Rdata")

