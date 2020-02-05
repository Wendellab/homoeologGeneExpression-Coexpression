########################
#### Bowtie2-HyLite ####
########################

# get A2 or D5 origin for each mapped ADs reads
load("R-01-hyliteArefDatasets.RData")
# "A2.Total"  "D5.Total"  "ADs.Total" "ADs.At"    "ADs.Dt"    "ADs.N"

s2o=function(x){unlist(lapply(strsplit(x,"-"),function(x)paste(paste0(x[2]),gsub("seed","",x[1]),x[3],sep="-")) )}

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
    system(paste0("samtools view -F 4 ../bowtie2hylite/resultsA2/", ADflag, ".sorted.bam | cut -f1 >AD.txt"))
    # out put A2 read ID
    system(paste0("zcat ../seedFastq/",A2flag,".cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>A2.txt"))
    system(paste0("zcat ../seedFastq/",D5flag,".cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>D5.txt"))
    # label AD reads
    system("perl labelAD.pl A2.txt D5.txt AD.txt >origin.txt")
    ori =fread("origin.txt",header=FALSE)
    #
    x=fread(paste0("../bowtie2hylite/resultsA2_011820/resultsA2_011820.ADs.", flag,".read.txt"),select=c(1,6), header=TRUE)
    x$origin = ori
    # x$id=gsub(".1.D","d", gsub(".1.A","a",x$id))
    At.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$origin=="A2"])[rownames(A2.Total)] )
    At.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$origin=="D5"])[rownames(A2.Total)] )
    Dt.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$origin=="D5"])[rownames(D5.Total)] )
    Dt.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$origin=="A2"])[rownames(D5.Total)] )
    # check if make sense, should be zero if correct
    print(sum(At.T[,flag]+At.F[,flag] - ADs.At[,flag], na.rm=T))
    print(sum(At.T[,flag]+At.F[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag] - ADs.Dt[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag], na.rm=T))
}
for(i in 1:22){
    flag= names(ADs.At)[i]
    print(flag)
    ADflag =gsub("[.]","-",flag)
    A2flag=names(A2.Total)[i]
    D5flag=names(D5.Total)[i]
    # output ADs read ID
    # samtools view -F 4 ADs-10-R1.q10.sort.bam | cut -f1,3 >AD.id.txt; cut -f1 AD.id.txt >AD.txt
    system(paste0("samtools view -F 4 ../bowtie2hylite/resultsA2/", ADflag, ".sorted.bam | cut -f1 >AD.txt"))
    # out put A2 read ID
    # zcat ../../seedFastq/A2-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>A2.txt
    # zcat ../../seedFastq/D5-10-R1.cut.fq.gz |grep '^@HWI'|sed 's/ .*//g' |sed 's/@//g'>D5.txt
    system(paste0("zcat ../flowerFastq/",A2flag,"_1.fq.gz |grep '^@K'|sed 's/ .*//g' |sed 's/@//g'>A2.txt"))
    system(paste0("zcat ../flowerFastq/",D5flag,"_1.fq.gz |grep '^@K'|sed 's/ .*//g' |sed 's/@//g'>D5.txt"))
    # label AD reads
    system("perl labelAD.pl A2.txt D5.txt AD.txt >origin.txt")
    ori =fread("origin.txt",header=FALSE)
    #
    x=fread(paste0("../bowtie2hylite/resultsA2_011820/resultsA2_011820.ADs.", flag,".read.txt"),select=c(1,6), header=TRUE)
    x$origin = ori
    # x$id=gsub(".1.D","d", gsub(".1.A","a",x$id))
    At.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$origin=="A2"])[rownames(A2.Total)] )
    At.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="A2" & x$origin=="D5"])[rownames(A2.Total)] )
    Dt.T[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$origin=="D5"])[rownames(D5.Total)] )
    Dt.F[,flag]  = as.numeric( table(x$GENE[x$CAT=="D5" & x$origin=="A2"])[rownames(D5.Total)] )
    # check if make sense, should be zero if correct
    print(sum(At.T[,flag]+At.F[,flag] - ADs.At[,flag], na.rm=T))
    print(sum(At.T[,flag]+At.F[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag] - ADs.Dt[,flag], na.rm=T))
    print(sum(Dt.T[,flag]+Dt.F[,flag], na.rm=T))
}

At.T[is.na(At.T)] = 0
At.F[is.na(At.F)] = 0
Dt.T[is.na(Dt.T)] = 0
Dt.F[is.na(Dt.F)] = 0
save( At.T, Dt.T, At.F, Dt.F, file="R-01-hyliteArefDatasets.true.RData")

################################
## Read assignment evaluation ##
################################
source("2.1.FUN.r")

lnames <- load('R-01-hyliteArefDatasets.RData')
lnames
#  "A2.Total"  "D5.Total"  "ADs.Total" "ADs.At"    "ADs.Dt"
lnames <- load('R-01-hyliteArefDatasets.true.RData')
lnames
# "At.T" "Dt.T"      "At.F"      "Dt.F"
lnames <- load('R-01-hyliteArefNetworkDatasets.RData')
lnames
# "coldata"       "networks"

# see if meets expectation
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
library(gplots)
pdf("s2.eval.hyliteAref.summary.pdf")
textplot(round(resEffi*100,2))
mtext("Efficiency (%)")
textplot(round(resDisc*100,2))
mtext("Discrepancy (%)")
textplot(round(resBinA$summary*100,2))
mtext("At Metrics (%)")
textplot(round(resBinD$summary*100,2))
mtext("Dt Metrics (%)")
dev.off()

info<-coldata
info$type <- c(rep("PE",22), rep("SE",11))
info$ADsLibSize <- colSums(ADs.Total)

hyliteArefSummary<-list(info=info, Efficiency=resEffi, Discrepancy=resDisc, binA=resBinA$summary, binB=resBinD$summary)
save(hyliteArefSummary, file="s2.assign_eval.hyliteAref.Rdata")

