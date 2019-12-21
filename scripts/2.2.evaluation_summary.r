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
c(mean(x)-2*sem,mean(x)+2*sem)

#########################
######### summary #######
#########################
library(ggplot2)
load("s2.assign_eval.polycat.Rdata")
load("s2.assign_eval.rsem.Rdata")
load("s2.assign_eval.hylite.Rdata")
load("s2.assign_eval.salmon.Rdata")
load("s2.assign_eval.kallisto.Rdata")

## check sample order
polycatSummary$info$sample == hyliteSummary$info$sample
rsemSummary$info$sample == hyliteSummary$info$sample
salmonSummary$info$sample == hyliteSummary$info$sample
kallistoSummary$info$sample == hyliteSummary$info$sample

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
methods = c("polycat","hylite","rsem", "salmon", "kallisto")
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

#########################
## Partition by %eflen ##
#########################

# geneLenM, output from "detectEffectiveRegion.r"
geneL <- read.table("eflenList.txt",header=TRUE,sep="\t")
geneL$ratio100 = geneL$theoretical100/geneL$true
geneL$ratio300 = geneL$theoretical300/geneL$true

methods = c("polycat","hylite","rsem", "salmon", "kallisto")
if(exists("sumTbl")){rm(sumTbl)}
for(m in methods)
{
    res = get(m)
    # bin genes by A to 6 classes
    breaks <- c(0, 0.80, 0.90, 0.95,1)
    res$c<-.bincode(x=res$percentageEffectM, b=breaks, FALSE)
    res$c[is.na(res$c)]=5
    cName<-c("[0-0.8)","[0.8-0.9)","[0.9-0.95)","[0.95-1)","[1]")
head(res)}

library(ggplot2)
library(dplyr)
# collect stats
load("s2.assign_eval.polycat.Rdata")
load("s2.assign_eval.rsem.Rdata")
load("s2.assign_eval.hylite.Rdata")
load("s2.assign_eval.salmon.Rdata")
load("s2.assign_eval.kallisto.Rdata")
polycat$pipeline="polycat"
hylite$pipeline="hylite"
rsem$pipeline="rsem"
salmon$pipeline="salmon"
kallisto$pipeline="kallisto"
df=rbind(polycat,hylite,rsem,salmon,kallisto)
# bin by %eflen
breaks <- c(0, 0.80, 0.90, 0.95,1)
df$c<-.bincode(x=df$percentageEffectM, b=breaks, FALSE)
df$c[is.na(df$c)]=5
cName<-c("[0-0.8)","[0.8-0.9)","[0.9-0.95)","[0.95-1)","[1]")
df$b = cName[df$c]
# clean
df$discrepancy[df$discrepancy>1]=1
df$efficiency[df$efficiency>1]=1
pdf("s2.eval_by_%eflen.pdf")
# Accuracy
df.summary <- df %>%
group_by(pipeline,b ) %>%
    summarise(
    sd = sd(accuracy,na.rm=T),
    mean = mean(accuracy,na.rm=T)
    )
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Accuracy")+theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# MCC
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(mcc,na.rm=T),
mean = mean(mcc,na.rm=T)
)
df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("MCC") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# Discrepancy
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(discrepancy,na.rm=T),
mean = mean(discrepancy,na.rm=T)
)
df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Discrepancy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# Efficency
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(efficiency,na.rm=T),
mean = mean(efficiency,na.rm=T)
)
df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Efficiency") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
dev.off()



---book--

## gene-wise metrics comparison, noting that sample order is different in hylite
message("Compare metrics between methods: ")
print("Compare metrics between methods: ")
pw<-combn(length(methods),2)
for(i in 1:ncol(pw)){
    print(paste0(methods[pw[1,i]]," vs ",methods[pw[2,i]]))
    res1<- get(methods[pw[1,i]])
    res2<- get(methods[pw[2,i]])
    for(m in c( "efficiency" , "discrepancy", "accuracy"))
    {
        select <- is.finite(res1[[m]]) & is.finite(res2[[m]])
        
        print(paste0("Testing ",m))
        print( t.test(res1[[m]][select],res2[[m]][select], paired=TRUE) )
        print( wilcox.test(res1[[m]][select],res2[[m]][select], paired=TRUE) )
    }
    
}

# compare At vs Dt
message("Compare metrics between homoeologs within each method: ")
print("Compare metrics between homoeologs within each method: ")
for(i in paste0(methods,"H")){
    resH<-get(i)
    for(m in c("efficiency","discrepancy","accuracy","precision","Fstat"))
    {
        # m<-"discrepancy"
        # res<-polycatH
        col<-grep(m,names(resH))
        select <- is.finite(resH[,col[1]]) & is.finite(resH[,col[2]])
        print(paste0(i,": ",m))
        print( t.test(resH[select,col[1]],resH[select, col[2]], paired=TRUE) )
        print( wilcox.test(resH[select,col[1]],resH[select, col[2]], paired=TRUE) )
    }
}

# multiple measurements as potential explainational variables to performance metrics

message("Linear regression analysis of method")
for(i in methods){
    print(paste0("Linear regression analysis of method: ",i))
    res<- get(i)
    
    for(m in c("accuracy","efficiency","discrepancy"))
    {
        message(paste0("...",m))
        res<-na.omit(res)
        
        model.null = glm(get(m) ~ 1,
        data=res)
        
        model.full = glm(get(m) ~ geneLenM + percentageEffectM + expression.log2 +sample,
        data=res)
        
        step(model.null,
        scope = list(upper=model.full),
        direction="both",
        test="Chisq",
        data=Data)
        
    }
}

# bar plots not very useful
library(reshape)
plotMetric=function(x, main=""){
    x$method=rownames(x)
    df<-melt(x,ids="method")
    df$method<- factor(df$method,levels = methods) # my order
    p<-ggplot(df, aes(x=method,y=value, fill=method))+ geom_bar(stat="identity", position=position_dodge()) +  geom_text( aes(label = round(value,4), y = value - 0.05), position = position_dodge(0.9), vjust = 0)+ facet_grid(variable~.) + ggtitle(main)
    print(p)
}
pdf("s2.assign_eval.summary.pdf")
for(i in names(sumTbl)){
    plotMetric(sumTbl[[i]],i)
}
# Library
library(fmsb)
ma = rep(1,4)
mi =rep(0.8,4)
a<- with(sumTbl,cbind(Precision$At,Accuracy$At, Fmeasure$At, MCC$total))
d<- with(sumTbl,cbind(Precision$Dt,Accuracy$Dt, Fmeasure$Dt, MCC$total))
rownames(a) = rownames(sumTbl$MCC)
rownames(d) = rownames(sumTbl$MCC)
colnames(a) = c("Precision","recall","F1 score", "MCC")
colnames(d)= c("Precision","recall","F1 score", "MCC")
a=as.data.frame(rbind(ma,mi,a))*100
d=as.data.frame(rbind(ma,mi,d))*100
#==================
# Plot 1: Default radar chart proposed by the library:
radarchart(a)
radarchart(d)
#==================
# Plot 2: Same plot with custom features
library(scales)
#show_col(hue_pal()(5))
colors = hue_pal()(5)
source("addTrans.r")
colors_border=addTrans(colors,240)
colors_in=addTrans(colors,40)
radarchart(a, axistype=1,
    #custom polygon
    pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(80,100,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 , title="At" )
legend(x=0.7, y=1, legend = rownames(a[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=1.2, pt.cex=3)
radarchart(d, axistype=1,
    #custom polygon
    pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(80,100,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 , title="Dt" )
legend(x=0.7, y=1, legend = rownames(d[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=1.2, pt.cex=3)

dev.off()

