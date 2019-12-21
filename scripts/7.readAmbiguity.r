options(scipen = 999)
###############
## Functions ##
###############
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
## Partition by %eflen ##
#########################
# bin criterion
breaks <- c(0, 0.80, 0.90, 0.95,1)
cName<-c("[0-0.8)","[0.8-0.9)","[0.9-0.95)","[0.95-1)","[1]")
# geneLenM, output from "detectEffectiveRegion.r"
geneL <- read.table("eflenList.txt",header=TRUE,sep="\t")
binEflen=function(x,breaks,cName){
    y<-.bincode(x, b=breaks, FALSE)
    y[is.na(y)]=5
    y = cName[y]
    return(y)
}
geneL$bin50<-binEflen(geneL$ratio50,breaks,cName)
geneL$bin100<-binEflen(geneL$ratio100,breaks,cName)
geneL$bin200<-binEflen(geneL$ratio200,breaks,cName)
geneL$bin300<-binEflen(geneL$ratio300,breaks,cName)
rownames(geneL)=geneL$id
# check
table(geneL$bin300)
round(table(geneL$bin300)/nrow(geneL)*100,1)
#    [0-0.8)  [0.8-0.9) [0.9-0.95)   [0.95-1)        [1] 
#      2561        345        220        344      34035 
#       6.8        0.9        0.6        0.9       90.7 
table(geneL$bin100)
round(table(geneL$bin100)/nrow(geneL)*100,1)
#   [0-0.8)  [0.8-0.9) [0.9-0.95)   [0.95-1)        [1] 
#      4342       2801       4849      10880      14633 
#      11.6        7.5       12.9       29.0       39.0 
save(geneL,file="s7.eflen.rdata")


############################################
######### Read assignment evaluation #######
############################################
library(ggplot2)
library(dplyr)
# load read assignment metrics, eflen info already incorporated
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
pdf("s7.eval_by_%eflen.pdf")
# Accuracy
df.summary <- df %>%
group_by(pipeline,b ) %>%
    summarise(
    sd = sd(accuracy,na.rm=T),
    mean = mean(accuracy,na.rm=T)
    )
accu=df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Accuracy")+theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# MCC
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(mcc,na.rm=T),
mean = mean(mcc,na.rm=T)
)
mcc=df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("MCC") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# Discrepancy
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(discrepancy,na.rm=T),
mean = mean(discrepancy,na.rm=T)
)
disc=df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Discrepancy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
# Efficency
df.summary <- df %>%
group_by(pipeline,b ) %>%
summarise(
sd = sd(efficiency,na.rm=T),
mean = mean(efficiency,na.rm=T)
)
effi=df.summary
ggplot(df.summary, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Efficiency") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
dev.off()

save(geneL, effi, disc, mcc, accu,file="s7.eflen.rdata")


###############################
######### DE evaluation #######
###############################
library(ggplot2)
library(dplyr)
library(ROCR)
load("s7.eflen.rdata")

# get batch for tissue samples
load("s2.assign_eval.polycat.Rdata")
info= unique(polycatSummary$info[,c("tissue","effectRegionRatio")])
info$batch =paste0(info$tissue,".Avs",info$tissue,".D")

# collect AUC for bins
rm(resultAUC))
methods = c("polycat","hylite","rsem","salmon","kallisto")
resL<-list.files("DE")
for(mapping in methods)
{
    for(method in c("deseq","ebseq"))
    {
        for(i in 1:nrow(info))
        {
            compr =info$batch[i]
            file.obs<- grep(paste(mapping,"obs",method,compr,sep="."),resL,value=TRUE)
            file.exp<- grep(paste(mapping,"exp",method,compr,sep="."),resL,value=TRUE)
            obs<-read.table(file=paste0("DE/", file.obs), sep="\t",header=TRUE)
            exp<-read.table(file=paste0("DE/", file.exp), sep="\t",header=TRUE)
            if(method=="deseq")
            {
                trues <- obs$padj<0.05 & !is.na(obs$padj)
                predictions<- rank(1-exp$padj)
            }
            if(method=="ebseq")
            {
                trues <- obs$Status=="DE"
                predictions<- rank(exp$PPDE)
            }
            # assign %eflen bin
            binType= gsub("ratio","bin", info$effectRegionRatio[i])
            bins = geneL[gsub("a$","",rownames(obs)),binType]
            # loop by Bin
            for(b in unique(bins))
            {
                select= which(bins==b)
                # ROC analysis
                pred <- prediction(predictions = predictions[select],labels= trues[select])
                roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
                auc.perf = performance(pred, measure = "auc")
                auc <- round(as.numeric(auc.perf@y.values ),3)
                temp = data.frame(pipeline=mapping, DEmethod=method, comparison=compr,Bin=b,AUC=auc)
                print(temp)
                if(exists("resultAUC")){resultAUC <- rbind(resultAUC, temp)}else{resultAUC=temp}
            }
        }          
    }
}
print(resultAUC)
resultAUC$Bin =factor(resultAUC$Bin, levels=c("[0-0.8)","[0.8-0.9)","[0.9-0.95)","[0.95-1)","[1]"))

# plot by bins
df.summary <- resultAUC %>%
group_by(pipeline,DEmethod, Bin ) %>%
    summarise(
    sd = sd(AUC,na.rm=T),
    mean = mean(AUC,na.rm=T)
    )
DE<-df.summary
pdf("s7.DE_by_eflen.pdf")
ggplot(df.summary[df.summary$DEmethod=="deseq",], aes(Bin, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("DE performance - AUC, DESeq2") + theme_bw()  + ylim(0.2,0.9) + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
ggplot(df.summary[df.summary$DEmethod=="ebseq",], aes(Bin, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("DE performance - AUC, EBSeq") + theme_bw()  + ylim(0.2,0.9) + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
ggplot(df.summary, aes(Bin, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("DE performance - AUC, EBSeq") + theme_bw()  + ylim(0.2,0.9) + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen") + facet_grid(. ~ DEmethod)
dev.off()

# save
save(geneL, effi, disc, mcc, accu,DE, file="s7.eflen.rdata")

###############################
######### DC evaluation #######
###############################
library(ggplot2)
library(dplyr)
load("s7.eflen.rdata")

## Analysis of homoeolog pairs
mm<-load("s4.DC.homoeoPair.Rdata")
mm # "hylite_rld"   "hylite_rpkm"  "polycat_rld"  "polycat_rpkm" "rsem_rld"   "rsem_rpkm"

# summarize sig results
rm(res)
for(i in mm)
{
    x<-get(i)
    bins = geneL[gsub("a$","",x$At),"bin300"]
    x$sig = ifelse(x$pValDiff<0.05 & !is.na(x$pValDiff), "DC","nonDC")
    temp<-data.frame(aggregate(x$sig, list(bins), function(x){length(which(x=="DC"))/length(x)}))
    temp$dataset=i
    if(exists("res")){res<-rbind(res,temp)}else{res=temp}
}
names(res)[1:2]=c("Bin","DC")
res$pipeline=gsub("_.*","",res$dataset)
res$transformation=gsub(".*_","",res$dataset)
res$Bin =factor(res$Bin, levels=c("[0-0.8)","[0.8-0.9)","[0.9-0.95)","[0.95-1)","[1]"))


# plot 
pdf("s7.DC_by_eflen.pdf")
ggplot(res, aes(Bin, DC)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + ggtitle("DC performance - %DC") + theme_bw()  + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen") + facet_grid(. ~ transformation)

ggplot(res[res$transformation=="rld",], aes(Bin, DC)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + ggtitle("DC performance - %DC, rld") + theme_bw()  + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen") 

ggplot(res[res$transformation=="log2rpkm",], aes(Bin, DC)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + ggtitle("DC performance - %DC, log2RPKM") + theme_bw()  + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
dev.off()

DC=res
# save
save(geneL, effi, disc, mcc, accu,DE,DC, file="s7.eflen.rdata")


####################################
######### Network evaluation #######
####################################
library(ggplot2)
library(dplyr)

load("s6.NC.rdata")
pdf("s7.NCcorr_by_eflen.pdf")
## per group
df.summary <- sumBinCorr %>%
group_by(pipeline,transformation, netType, bins ) %>%
    summarise(
    sd = sd(corr,na.rm=T),
    mean = mean(corr,na.rm=T)
    )
df.summary
ggplot(df.summary, aes(bins, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Node connectivity - corr") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen") + facet_wrap(. ~ transformation + netType)
## ignore normalization and netType
df.summary <- sumBinCorr %>%
group_by(bins,pipeline ) %>%
    summarise(
    sd = sd(corr,na.rm=T),
    mean = mean(corr,na.rm=T)
    )
df.summary
ggplot(df.summary, aes(bins, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Node connectivity - corr") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line()) + xlab("%Eflen")
dev.off()

NC=df.summary
# save
save(geneL, effi, disc, mcc, accu,DE,DC, NC,file="s7.eflen.rdata")



#################
## make Figure ##
#################
library(ggplot2)
library(ggpubr)
library(dplyr)
load("s7.eflen.rdata")

head(geneL)
size=nrow(geneL)
df50<-as.data.frame(table(geneL$bin50)/size); df50$seq ="50"
df100<-as.data.frame(table(geneL$bin100)/size); df100$seq ="100"
df200<-as.data.frame(table(geneL$bin200)/size); df200$seq ="200"
df300<-as.data.frame(table(geneL$bin300)/size); df300$seq ="300"
df=rbind(df50,df100,df200,df300)
df$seq=factor(df$seq,levels=c("50","100","200","300"))

abi =c( "(0.2-1]", "(0.1-0.2]", "(0.05-0.1]", "(0-0.05]", "[0]")
names(abi) = c( "[0-0.8)", "[0.8-0.9)", "[0.9-0.95)", "[0.95-1)", "[1]")


p1= ggplot(df, aes(Var1, Freq)) + geom_line(aes(group=seq,color=seq, linetype=seq)) +  geom_point(aes( color=seq)) + ggtitle("Bin size")+theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank(),legend.position="bottom")+ scale_color_grey(start = 0.8, end = 0.2) 

p2= ggplot(effi, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Efficiency") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

p3= ggplot(disc, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Discrepancy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

p4= ggplot(accu, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Accuracy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

p5= ggplot(mcc, aes(b, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("MCC") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

DE$pipeline=factor(DE$pipeline,levels=c("hylite","kallisto","polycat","rsem","salmon"))
p6= ggplot(DE, aes(Bin, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("DE - AUC") + theme_bw()  + ylim(0.2,0.9) + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ DEmethod)

p7 = ggplot(DC, aes(Bin, DC)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + ggtitle("DC - %change") + theme_bw()  + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ transformation)

p8= ggplot(NC, aes(bins, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("k - correlation") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))


p0=c()

pdf("s7.figure.pdf")
ggarrange(p1,p2,p3,p0,p4,p5,p6,p7,p8,
          labels = c("A", "B","C","", "D","E","F","G","H"),
          ncol = 3, nrow = 3,common.legend=TRUE) 
p1;p2;p3;p4;p5;p6;p7;p8;   
dev.off()


o=c("(0)", "(0-0.05)", "(0.05-0.1)", "(0.1-0.2)","(0.2-1)")
abi =c( "(0.2-1)", "(0.1-0.2)", "(0.05-0.1)", "(0-0.05)", "(0)")
names(abi) = c( "[0-0.8)", "[0.8-0.9)", "[0.9-0.95)", "[0.95-1)", "[1]")

df$Var2 =factor(abi[as.character(df$Var1)], levels = o)
p1= ggplot(df, aes(Var2, Freq)) + geom_line(aes(group=seq,color=seq, linetype=seq)) +  geom_point(aes( color=seq)) + ggtitle("Bin size")+theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank(),legend.position="bottom")+ scale_color_grey(start = 0.8, end = 0.2)

effi$b2 =factor(abi[as.character(effi$b)], levels = o)
p2= ggplot(effi, aes(b2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Efficiency") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

disc$b2 =factor(abi[as.character(disc$b)], levels = o)
p3= ggplot(disc, aes(b2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Discrepancy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

accu$b2 =factor(abi[as.character(accu$b)], levels = o)
p4= ggplot(accu, aes(b2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("Accuracy") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

mcc$b2 =factor(abi[as.character(mcc$b)], levels = o)
p5= ggplot(mcc, aes(b2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,  color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("MCC") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

DE$Bin2 =factor(abi[as.character(DE$Bin)], levels = o)
DE$pipeline=factor(DE$pipeline,levels=c("hylite","kallisto","polycat","rsem","salmon"))
p6= ggplot(DE, aes(Bin2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("DE - AUC") + theme_bw()  + ylim(0.2,0.9) + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ DEmethod)

DC$Bin2 =factor(abi[as.character(DC$Bin)], levels = o)
p7 = ggplot(DC, aes(Bin2, DC)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + ggtitle("DC - %change") + theme_bw()  + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ transformation)

NC$bins2 =factor(abi[as.character(NC$bins)], levels = o)
p8= ggplot(NC, aes(bins2, mean)) + geom_line(aes(group=pipeline,linetype=pipeline,   color=pipeline),position = position_dodge(0.3)) +  geom_point(aes( color=pipeline),position = position_dodge(0.3)) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = pipeline), position = position_dodge(0.3),width = 0.2) + ggtitle("k - correlation") + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))


p0=c()

pdf("s7.figure.v120819.pdf")
ggarrange(p1,p2,p3,p0,p4,p5,p6,p7,p8,
          labels = c("A", "B","C","", "D","E","F","G","H"),
          ncol = 3, nrow = 3,common.legend=TRUE)
p1;p2;p3;p4;p5;p6;p7;p8;
dev.off()
