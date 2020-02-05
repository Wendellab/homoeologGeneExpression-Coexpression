options(scipen = 999)
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
#computation of the standard error of the mean
sem=function(x){sd(x)/sqrt(length(x))}
#95% confidence intervals of the mean
# ci=c(mean(x)-2*sem,mean(x)+2*sem)

#########################
######### summary #######
#########################
library(ggplot2)
list.files(pattern="s2.assign_eval.*Rdata")
load("s2.assign_eval.polycat.Rdata")
load("s2.assign_eval.rsem.Rdata")
load("s2.assign_eval.hylite.Rdata")
load("s2.assign_eval.salmon.Rdata")
load("s2.assign_eval.kallisto.Rdata")
# 2020 revision
load("s2.assign_eval.bowtie.Rdata")
load("s2.assign_eval.eaglerc.Rdata")
load("s2.assign_eval.hyliteAref.Rdata")


## check sample order
polycatSummary$info$sample == hyliteSummary$info$sample
rsemSummary$info$sample == hyliteSummary$info$sample
salmonSummary$info$sample == hyliteSummary$info$sample
kallistoSummary$info$sample == hyliteSummary$info$sample
polycatSummary$info$sample == bowtieSummary$info$sample
polycatSummary$info$sample == eaglercSummary$info$sample
polycatSummary$info$sample == hyliteArefSummary$info$sample


## double check libSize
size<-data.frame(polycat=polycatSummary$info$ADsLibSize, hylite=hyliteSummary$info$ADsLibSize, rsem=rsemSummary$info$ADsLibSize, salmon=salmonSummary$info$ADsLibSize, kallisto=kallistoSummary$info$ADsLibSize )
rownames(size)=polycatSummary$info$sample
pairs( size, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main ="Library size")
pairs( size[1:22,], lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main ="Library size: PE")
pairs( size[23:33,], lower.panel=panel.smooth, upper.panel=panel.cor, pch=20,  main ="Library size: SE")
# bar plots not very useful
# library(reshape)
# df<-size
# df$sample<-rownames(df)
# df<-melt(df,id="sample")
# ggplot(df, aes(x=sample,y=value, fill=variable))+ geom_bar(stat="identity", position=position_dodge())

## summary
methods = c("polycat","hylite","hyliteAref","eaglerc", "rsem", "salmon", "kallisto","bowtie")
if(exists("sumTbl")){rm(sumTbl)}
for(m in methods)
{
    x = get(paste0(m,"Summary"))
    # Efficiency
    e=cbind(t(x$Efficiency[1,]), apply(x$Efficiency,2,sem)); rownames(e)[1]="Efficiency"
    # Discrepancy
    d=cbind(t(x$Discrepancy[1,]), apply(x$Discrepancy,2,sem)); rownames(d)[1]="Discrepancy"
    # other metrics
    oa=cbind(t(x$binA[1,]), apply(x$binA,2,sem))
    od=cbind(t(x$binB[1,]), apply(x$binB,2,sem))
    # combine
    me= rbind(e,d,oa[1,],od[1,],oa[2,],od[2,], oa[4,],od[4,],oa[3,],oa[5,])
    rownames(me)[7:14]=c("Precision.At","Precision.Dt","Recall.At","Recall.Dt","F1.At","F1.Dt","Accuracy","MCC")
    colnames(me)=paste0(m,c(".overall",".se"))
    if(!exists("sumTbl")) {sumTbl=me}else {sumTbl=cbind(sumTbl,me)}
}
round(sumTbl*100,1)
write.table(round(sumTbl*100,1),file="s2.evaluation_summary.txt")
