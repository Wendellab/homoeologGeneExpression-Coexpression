library(ggplot2)
library(agricolae)
library(RColorBrewer)

options(scipen=999) # disable scientifc format
netType<-data.frame( net=1:11, edge=c(rep("binary",8),rep("weighted",3)),metric=c(rep("rank",4),rep("Zscore",4),rep("WGCNA",3)),cutoff=c(0.95,0.99,0.995,0.999,1.5,2,2.5,3,1,12,24))
netType$desc <- paste(netType$edge,netType$metric, netType$cutoff,sep=".")
netType
#  net     edge metric cutoff              desc
#    1   binary   rank  0.950  binary.rank.0.95
#    2   binary   rank  0.990  binary.rank.0.99
#    3   binary   rank  0.995 binary.rank.0.995
#    4   binary   rank  0.999 binary.rank.0.999
#    5   binary Zscore  1.500 binary.Zscore.1.5
#    6   binary Zscore  2.000   binary.Zscore.2
#    7   binary Zscore  2.500 binary.Zscore.2.5
#    8   binary Zscore  3.000   binary.Zscore.3
#    9 weighted  WGCNA  1.000  weighted.WGCNA.1
#   10 weighted  WGCNA 12.000 weighted.WGCNA.12
#   11 weighted  WGCNA 24.000 weighted.WGCNA.24

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
}

#######################
## Node connectivity ##
#######################

load("s6.NC.rdata") 
# sumCorr, 7 pipeline x 2 norm x 11 networks x 10 permu = 1540 rows
# sumBinCorr, 5 bin x 1540 = 5500
sumCorr$netType=factor(sumCorr$netType,levels=netType$desc)
sumCorr$pipeline=factor(sumCorr$pipeline,levels=c("polycat","hylite","eaglerc","rsem","kallisto","salmon","bowtie") )

#--------- NC correlation ---------#
pdf("s6.NC_corr.pdf") # FigS1B
p.cor=ggplot(data=sumCorr, aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Node Connectivity") + scale_fill_brewer(palette="Set1")
print(p.cor)
dev.off()
# check
fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value              Pr(>F)    
# pipeline          6  0.077  0.0128   1.243                0.281
# transformation    1  0.629  0.6286  60.833   0.0000000000000116 ***
# netType          10 14.333  1.4333 138.716 < 0.0000000000000002 ***
# Residuals      1483 15.323  0.0103
## pipeline is not robust
out <- duncan.test(a,trt="pipeline"); out$groups
#          correlation groups
# eaglerc    0.9456054      a
# salmon     0.9453349      a
# kallisto   0.9410813      a
# bowtie     0.9382415      a
# rsem       0.9370625      a
# polycat    0.9271907      a
# hylite     0.9269196      a
out <- duncan.test(a,trt="transformation"); out$groups
#          correlation groups
# log2rpkm   0.9572866      a
# rld        0.9163024      b
out <- duncan.test(a,trt="netType"); out$groups
#                   correlation groups
# weighted.WGCNA.1    0.9890500      a
# weighted.WGCNA.12   0.9876962      a
# weighted.WGCNA.24   0.9870191      a
# binary.rank.0.95    0.9849158      a
# binary.rank.0.99    0.9825111      a
# binary.rank.0.995   0.9803198      a
# binary.Zscore.1.5   0.9767840      a
# binary.rank.0.999   0.9715946      a
# binary.Zscore.2     0.9388603      b
# binary.Zscore.2.5   0.7751795      c
# binary.Zscore.3     0.6591825      d


#--------- Density ---------#
# Density A vs D vs A-D
library(ggpubr)
library(reshape2)
library(dplyr)
res=melt(data = sumCorr[,1:12], id.vars = names(sumCorr)[1:4], measure.vars = names(sumCorr)[5:12],value.name ="density")
res$subnetwork = "all"
res$subnetwork[grep(".A$",res$variable)] = "A"
res$subnetwork[grep("[.]D$",res$variable)] = "D"
res$subnetwork[grep(".interAD$",res$variable)] = "interAD"
res$dataset = factor(gsub("[.].*","", gsub("density[.]|", "", res$variable)))
res$group= paste(res$pipeline,res$transformation, res$netType, sep="_")
pdf("s6.NC_density.pdf")
for(each in unique(res$group)){
    print(each)  # 11x 5 x 2 =110
    df=res[res$group==each&res$subnetwork!="all",]
    title=paste0("pipeline-",unique(df$pipeline),", transformation-",unique(df$transformation),", netType-",unique(df$netType))
    # boxplot
    # p=ggplot(df, aes(x=subnetwork, y=density)) + geom_boxplot(aes(color = dataset)) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
    # line plot with errorbar
    df.summary2 <- df[,6:8]  %>%
       group_by(subnetwork, dataset) %>%
       summarise( sd = sd(density), density = mean(density) )
    df.summary2
    p=ggplot(df, aes(x=subnetwork, y=density)) + geom_jitter( aes(color = dataset), position = position_jitter(0.2)) + geom_line( aes(group = dataset, color = dataset), data = df.summary2) + geom_errorbar( aes(ymin = density-sd, ymax = density+sd, color = dataset), data = df.summary2, width = 0.1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
    print(p)
}
dev.off()
# make a summary table of density mean and sd as supplemental table
res.summary <- res  %>%
group_by(subnetwork, dataset, group) %>%
summarise( sd = sd(density), density = mean(density) )%>%
arrange(group)
res.summary
write.table(res.summary,file = "s6.NC_density.txt", sep="\t",row.names=FALSE )
# make wide form to report
w= data.frame(res.summary)
w$group2=paste(w$group,w$dataset)
all = w[w$subnetwork=="all",c("group2","density","sd")];names(all)=c("pipeline_transformation_network dataset","mean","sd")
A = w[w$subnetwork=="A",c("group2","density","sd")];names(A)[2:3]=c("A_mean","A_sd")
D = w[w$subnetwork=="D",c("group2","density","sd")];names(D)[2:3]=c("D_mean","D_sd")
AD = w[w$subnetwork=="interAD",c("group2","density","sd")];names(AD)[2:3]=c("interAD_mean","interAD_sd")
unique(all[,1]==AD$group2)
unique(A$group2==AD$group2)
unique(D$group2==AD$group2)
w= cbind(all, A[,-1],D[,-1], AD[,-1])
w$AvsD=ifelse(w$A_mean>w$D_mean,"A>D","A<D")
w$AvsD.pvalue=apply(w,1,function(x){y=as.numeric(x[4:7]);p=t.test2(m1=y[1],m2=y[3],s1=y[2],s2=y[4],n1=10,n2=10 )["p-value"]; return(p)})
write.table(w,file = "s6.NC_density2.txt", sep="\t",row.names=FALSE )
# ANOVA
fit<-lm(density~pipeline+transformation+netType+dataset+subnetwork, data= res)
print(summary(a<-aov(fit)))
#                   Df Sum Sq Mean Sq    F value               Pr(>F)
# pipeline           6   0.00   0.000      8.992       0.000000000794 ***
# transformation     1   0.06   0.060   1407.429 < 0.0000000000000002 ***
# netType           10 270.76  27.076 634515.909 < 0.0000000000000002 ***
# dataset            1   0.00   0.000      1.411                0.235
# subnetwork         3   0.04   0.014    323.697 < 0.0000000000000002 ***
# Residuals      12298   0.52   0.000
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.05857909      a
# rld      0.05416326      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                       auroc groups
# weighted.WGCNA.1  0.5227969729880      a
# binary.rank.0.95  0.0525930439248      b
# weighted.WGCNA.12 0.0212623248138      c
# binary.rank.0.99  0.0110134956537      d
# binary.rank.0.995 0.0055252181053      e
# weighted.WGCNA.24 0.0042932233319      f
# binary.Zscore.1.5 0.0013605687646      g
# binary.rank.0.999 0.0011678539711      g
# binary.Zscore.2   0.0000641870542      h
# binary.Zscore.2.5 0.0000051708123      h
# binary.Zscore.3   0.0000008712142      h
out <- duncan.test(a,trt="subnetwork"); print(out$groups)
#            density groups
# A       0.05920235      a
# D       0.05628238      b
# all     0.05591388      c
# interAD 0.05408609      d
# ----
## make figure
p=list()
res1=res[res$pipeline=="polycat"&res$transformation=="log2rpkm"&res$subnetwork!="all",]
for(i in 1:11){
    df=res1[res1$netType==netType$desc[i],]
    title=gsub("binary.|weighted.","",netType$desc[i])
    # line plot with errorbar
    df.summary2 <- df[,5:8]  %>%
       group_by(subnetwork, dataset) %>%
       summarise( sd = sd(density), density = mean(density) )
    df.summary2
    p[[i]]=ggplot(df, aes(x=subnetwork, y=density)) + geom_jitter( aes(color = dataset), position = position_jitter(0.2)) + geom_line( aes(group = dataset, color = dataset), data = df.summary2) + geom_errorbar( aes(ymin = density-sd, ymax = density+sd, color = dataset), data = df.summary2, width = 0.1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom") + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
}
q=list()
res1=res[res$pipeline=="eaglerc"&res$transformation=="log2rpkm"&res$subnetwork!="all",]
for(i in 1:11){
    df=res1[res1$netType==netType$desc[i],]
    title=gsub("binary.|weighted.","",netType$desc[i])
    # line plot with errorbar
    df.summary2 <- df[,5:8]  %>%
       group_by(subnetwork, dataset) %>%
       summarise( sd = sd(density), density = mean(density) )
    df.summary2
    q[[i]]=ggplot(df, aes(x=subnetwork, y=density)) + geom_jitter( aes(color = dataset), position = position_jitter(0.2)) + geom_line( aes(group = dataset, color = dataset), data = df.summary2) + geom_errorbar( aes(ymin = density-sd, ymax = density+sd, color = dataset), data = df.summary2, width = 0.1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() + theme(panel.border=element_blank(),axis.line.x=element_line(),axis.title.x=element_blank(),legend.position="bottom") + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
}
pdf("s6.figure6.pdf", height=5,width=8)
ggarrange(p[[2]],p[[6]],p[[10]],q[[2]],q[[6]],q[[10]],labels=c("A","B","C","E","F","G"),ncol = 3, nrow = 2,common.legend=TRUE)
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],ncol = 4, nrow = 3,common.legend=TRUE)
ggarrange(q[[1]],q[[2]],q[[3]],q[[4]],q[[5]],q[[6]],q[[7]],q[[8]],q[[9]],q[[10]],q[[11]],ncol = 4, nrow = 3,common.legend=TRUE)
dev.off()

#--------- Density ---------#
res=melt(data = sumCorr[,c(1:4,13:20)], id.vars = names(sumCorr)[1:4], measure.vars = names(sumCorr)[13:20],value.name ="path")
res=res[grep("WGCNA",res$netType,invert=T),]
res$subnetwork = "all"
res$subnetwork[grep(".A$",res$variable)] = "A"
res$subnetwork[grep("[.]D$",res$variable)] = "D"
res$subnetwork[grep(".interAD$",res$variable)] = "interAD"
res$dataset = factor(gsub("[.].*","", gsub("path[.]|", "", res$variable)))
res$group= paste(res$pipeline,res$transformation, res$netType, sep="_")
pdf("s6.NC_path.pdf")
for(each in unique(res$group)){
    print(each)  # 11x 5 x 2 =110
    df=res[res$group==each&res$subnetwork!="all",]
    title=paste0("pipeline-",unique(df$pipeline),", transformation-",unique(df$transformation),", netType-",unique(df$netType))
    # boxplot
    # p=ggplot(df, aes(x=subnetwork, y=density)) + geom_boxplot(aes(color = dataset)) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
    # line plot with errorbar
    df.summary2 <- df[,6:8]  %>%
       group_by(subnetwork, dataset) %>%
       summarise( sd = sd(path), path = mean(path) )
    df.summary2
    p=ggplot(df, aes(x=subnetwork, y=path)) + geom_jitter( aes(color = dataset), position = position_jitter(0.2)) + geom_line( aes(group = dataset, color = dataset), data = df.summary2) + geom_errorbar( aes(ymin =path-sd, ymax = path+sd, color = dataset), data = df.summary2, width = 0.1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
    print(p)
}
dev.off()
# make a summary table of density mean and sd as supplemental table
res.summary <- res  %>%
group_by(subnetwork, dataset, group) %>%
summarise( sd = sd(density), density = mean(density) )%>%
arrange(group)
res.summary
write.table(res.summary,file = "s6.NC_path.txt", sep="\t",row.names=FALSE )

#############################
## Functional connectivity ##
#############################

load("s6.FC.rdata") # sumCorr, 2 func x 7 pipeline x 2 norm x 11 net x 10 permu = 3080 rows
sumCorr$netType=factor(sumCorr$netType,levels=netType$desc)
sumCorr$pipeline=factor(sumCorr$pipeline,levels=c("polycat","hylite","eaglerc","rsem","kallisto","salmon","bowtie") )


fit<-lm(correlation~pipeline+transformation+netType+Func,data=sumCorr )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value               Pr(>F)
# pipeline          6   0.06  0.0105   0.421              0.86552
# transformation    1   0.26  0.2582  10.374              0.00129 **
# netType          10  26.50  2.6499 106.487 < 0.0000000000000002 ***
# Func              1   0.07  0.0690   2.773              0.09599 .
# Residuals      2906  72.32  0.0249

# GO 
fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr[sumCorr$Func=="go",] )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value               Pr(>F)
# pipeline          6   0.09  0.0149   0.584                0.743
# transformation    1   0.47  0.4653  18.186            0.0000213 ***
# netType          10  20.01  2.0006  78.187 < 0.0000000000000002 ***
# Residuals      1464  37.46  0.0256

# KEGG
fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr[sumCorr$Func=="kegg",] )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value              Pr(>F)
# pipeline          6   0.13  0.0213   0.941               0.464
# transformation    1   0.00  0.0010   0.044               0.834
# netType          10   8.88  0.8878  39.326 <0.0000000000000002 ***
# Residuals      1425  32.17  0.0226

pdf("s6.FC_corr.pdf")
p.go=ggplot(data=sumCorr[sumCorr$Func=="go",], aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Function Connectivity: go")+ scale_fill_brewer(palette="Set1")
p.kegg=ggplot(data=sumCorr[sumCorr$Func=="kegg",], aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Function Connectivity: kegg")+ scale_fill_brewer(palette="Set1")
print(p.go)
print(p.kegg)
dev.off()

# FC AUROC
pp1=ggplot(data=sumCorr[sumCorr$Func=="go",], aes(x=netType, y=exp, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Expected Functional Connectivity: ",i))+ scale_fill_brewer(palette="Set1")
pp2=ggplot(data=sumCorr[sumCorr$Func=="go",], aes(x=netType, y=obs, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Observed Functional Connectivity: ",i))+ scale_fill_brewer(palette="Set1")
pp3=ggplot(data=sumCorr[sumCorr$Func=="kegg",], aes(x=netType, y=exp, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Expected Functional Connectivity: ",i))+ scale_fill_brewer(palette="Set1")
pp4=ggplot(data=sumCorr[sumCorr$Func=="kegg",], aes(x=netType, y=obs, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Observed Functional Connectivity: ",i))+ scale_fill_brewer(palette="Set1")
pdf("s6.FC_auroc.pdf") # FigS2
print(pp1,pp2,pp3,pp4)
dev.off()
pdf("s6.figureS2.pdf", height=5,width=8)
ggarrange(pp1,pp2,pp3,pp4,ncol = 2, nrow = 2,common.legend=TRUE)
dev.off()

fit<-lm(obs~pipeline+transformation+netType,data= sumCorr[sumCorr$Func=="go",])
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq   F value               Pr(>F)
# pipeline          6  0.002  0.0003     5.979           0.00000344 ***
# transformation    1  0.016  0.0161   361.036 < 0.0000000000000002 ***
# netType          10  7.882  0.7882 17672.970 < 0.0000000000000002 ***
# Residuals      1522  0.068  0.0000
out <- duncan.test(a,trt="pipeline"); print(out$groups)
#                obs groups
# bowtie   0.5897346      a
# rsem     0.5886820     ab
# kallisto 0.5884983     ab
# hylite   0.5881554      b
# polycat  0.5881243      b
# salmon   0.5867275      c
# eaglerc  0.5866103      c
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.5759466      a
# rld      0.5711339      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                         obs groups
# weighted.WGCNA.12 0.6992980      a
# weighted.WGCNA.24 0.6965607      b
# binary.rank.0.95  0.6466008      c
# weighted.WGCNA.1  0.6444411      d
# binary.rank.0.99  0.6031047      e
# binary.rank.0.995 0.5807028      f
# binary.Zscore.1.5 0.5443248      g
# binary.rank.0.999 0.5409898      h
# binary.Zscore.2   0.5113248      i
# binary.Zscore.2.5 0.5013521      j
# binary.Zscore.3   0.5001369      j

fit<-lm(obs~pipeline+transformation+netType,data= sumCorr[sumCorr$Func=="kegg",])
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq  F value              Pr(>F)
# pipeline          6  0.001  0.0002    1.401               0.211
# transformation    1  0.032  0.0318  211.300 <0.0000000000000002 ***
# netType          10  9.182  0.9182 6101.442 <0.0000000000000002 ***
# Residuals      1522  0.229  0.0002
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
out <- duncan.test(a,trt="pipeline"); print(out$groups)
#                obs groups
# bowtie   0.5967204      a
# rsem     0.5966504      a
# kallisto 0.5965309      a
# hylite   0.5961833     ab
# eaglerc  0.5958560     ab
# polycat  0.5957867     ab
# salmon   0.5938914      b
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.6004896      a
# rld      0.5914015      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                         obs groups
# weighted.WGCNA.24 0.7162354      a
# weighted.WGCNA.12 0.7108080      b
# binary.rank.0.95  0.6664660      c
# weighted.WGCNA.1  0.6489908      d
# binary.rank.0.99  0.6143068      e
# binary.rank.0.995 0.5892329      f
# binary.Zscore.1.5 0.5496298      g
# binary.rank.0.999 0.5456207      h
# binary.Zscore.2   0.5128992      i
# binary.Zscore.2.5 0.5011782      j
# binary.Zscore.3   0.5000335      j


library(ggpubr)
pdf("s6.FC_auroc.pdf")
for(func in c("go","kegg")){
    print(func)
    rdatafiles = grep(func, grep("rdata",list.files("s6.FC"),value=TRUE),value=TRUE)
    length(rdatafiles)
    rm(sumRes)
    for(i in rdatafiles){
        title=paste0(func,": pipeline = ",gsub("_.*","",i),", transformation = ",gsub(".*_|[.].*","",i),", netType = ",gsub(".*rld.|.*rpkm.","",gsub(paste0(".",func,".*"),"",i)))
        l=load(paste0("s6.FC/",i))  #"obs.fc" "exp.fc"
        # collect
        obsA= data.frame(auc=obs.fc$A[,"auc"], genome="A",dataset="obs", term=rownames(obs.fc$A))
        obsD= data.frame(auc=obs.fc$D[,"auc"], genome="D",dataset="obs", term=rownames(obs.fc$D))
        expA= data.frame(auc=exp.fc$A[,"auc"], genome="A",dataset="exp", term=rownames(exp.fc$A))
        expD= data.frame(auc=exp.fc$D[,"auc"], genome="D",dataset="exp", term=rownames(exp.fc$D))
        # make df, plot box with signif A vs D for each dataset
        df=rbind(obsA,obsD,expA,expD)
        p=ggplot(df, aes(x=dataset, y=auc)) + geom_boxplot(aes(color = genome)) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = genome),label = "p.signif") + ggtitle(title)
        print(p)
    }
}
dev.off()

