library(ggplot2)
library(agricolae)

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

#######################
## Node connectivity ##
#######################

load("s6.NC.rdata") 
# sumCorr, 5 pipeline x 2 norm x 11 x 10 permu = 1100 rows
# sumBinCorr, 5 bin x 1100 = 5500

sumCorr$netType=factor(sumCorr$netType,levels=netType$desc)
pdf("s6.NC_corr.pdf")
ggplot(data=sumCorr, aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Node Connectivity")
#ggplot(data=sumCorr, aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Node Connectivity")+ facet_grid(transformation~.)
dev.off()

fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value              Pr(>F)    
# pipeline          4  0.037  0.0091   0.992               0.411    
# transformation    1  1.661  1.6610 180.189 <0.0000000000000002 ***
# netType          10 16.611  1.6611 180.201 <0.0000000000000002 ***
# Residuals      1084  9.992  0.0092                                
## pipeline is not robust

out <- duncan.test(a,trt="pipeline"); out$groups
#          correlation groups
# salmon     0.9270970      a
# kallisto   0.9246529      a
# rsem       0.9208376      a
# hylite     0.9145941      a
# polycat    0.9119681      a

out <- duncan.test(a,trt="transformation"); out$groups
#          correlation groups
# log2rpkm   0.9586887      a
# rld        0.8809711      b

out <- duncan.test(a,trt="netType"); out$groups
#                   correlation groups
# weighted.WGCNA.1    0.9837092      a
# weighted.WGCNA.24   0.9809889      a
# weighted.WGCNA.12   0.9800469      a
# binary.rank.0.99    0.9776646      a
# binary.rank.0.95    0.9774696      a
# binary.rank.0.995   0.9761518      a
# binary.Zscore.1.5   0.9723401      a
# binary.rank.0.999   0.9704955      a
# binary.Zscore.2     0.9346637      b
# binary.Zscore.2.5   0.7972762      c
# binary.Zscore.3     0.5673226      d

# Density A vs D vs A-D
library(ggpubr)
library(reshape2)
library(dplyr)
res=melt(data = sumCorr[,-1], id.vars = names(sumCorr)[8:11], measure.vars = names(sumCorr)[2:7],value.name ="density")
res$subnetwork = factor(gsub(".*[.]", "", res$variable))
res$dataset = factor(gsub("[.].*","", gsub("density[.]|", "", res$variable)))
res$group= paste(res$pipeline,res$transformation, res$netType, sep="_")
# ANOVA
fit<-lm(density~pipeline+transformation+netType+dataset+subnetwork,data= res)
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq   F value  Pr(>F)    
# pipeline          4   0.00   0.000 7.290e+00 7.46e-06 ***
# transformation    1   0.03   0.028 7.474e+02  < 2e-16 ***
# netType          10 144.16  14.416 3.882e+05  < 2e-16 ***
# dataset           1   0.00   0.000 3.640e-01    0.547    
# subnetwork        2   0.03   0.013 3.431e+02  < 2e-16 ***
# Residuals      6581   0.24   0.000                            
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.05826588      a
# rld      0.05416452      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                       auroc groups
# weighted.WGCNA.1  0.5212159183      a
# binary.rank.0.95  0.0526025594      b
# weighted.WGCNA.12 0.0203146053      c
# binary.rank.0.99  0.0111043032      d
# binary.rank.0.995 0.0057319997      e
# weighted.WGCNA.24 0.0041795932      f
# binary.Zscore.1.5 0.0014231346      g
# binary.rank.0.999 0.0013273155      g
# binary.Zscore.2   0.0001947385      h
# binary.Zscore.2.5 0.0001386831      h
# binary.Zscore.3   0.0001343776      h
out <- duncan.test(a,trt="subnetwork"); print(out$groups)
#            auroc groups
# A       0.05867968      a
# D       0.05609431      b
# interAD 0.05387162      c

pdf("s6.NC_density.pdf")
for(each in unique(res$group)){
    print(each)  # 11x 5 x 2 =110
    df=res[res$group==each,]
    title=paste0("Density: pipeline-",unique(df$pipeline),", transformation-",unique(df$transformation),", netType-",unique(df$netType))
    
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

## make figure
p=list()
res1=res[res$pipeline=="polycat"&res$transformation=="log2rpkm",]
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

pdf("s6.figure6.pdf", height=5,width=8)
ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],ncol = 4, nrow = 3,common.legend=TRUE) 
ggarrange(p[[2]],p[[6]],p[[10]],labels=c("A","B","C"),ncol = 3, nrow = 2,common.legend=TRUE) 
dev.off()


#############################
## Functional connectivity ##
#############################

load("s6.FC.090719.rdata") # sumCorr, 2 func x 5 pipeline x 2 norm x 11 x 7 permu = 1540 rows

fit<-lm(correlation~pipeline+transformation+netType+Func,data=sumCorr )
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq F value              Pr(>F)    
# pipeline          4   0.02   0.004   0.150              0.9631    
# transformation    1   0.09   0.089   3.302              0.0694 .  
# netType          10  32.18   3.218 118.913 <0.0000000000000002 ***
# Func              1   0.00   0.001   0.037              0.8484    
# Residuals      1519  41.11   0.027                                

# GO 
fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr[sumCorr$Func=="go",] )
print(summary(a<-aov(fit)))
#                 Df Sum Sq Mean Sq F value               Pr(>F)    
# pipeline         4  0.033  0.0082   0.444              0.77679    
# transformation   1  0.169  0.1688   9.153              0.00258 ** 
# netType         10 20.114  2.0114 109.072 < 0.0000000000000002 ***
# Residuals      644 11.876  0.0184                                 

# KEGG
fit<-lm(correlation~pipeline+transformation+netType,data=sumCorr[sumCorr$Func=="kegg",] )
print(summary(a<-aov(fit)))
# print(summary(a<-aov(fit)))
#                 Df Sum Sq Mean Sq F value              Pr(>F)    
# pipeline         4  0.141  0.0353   0.953               0.433    
# transformation   1  0.006  0.0061   0.166               0.684    
# netType         10 11.665  1.1665  31.497 <0.0000000000000002 ***
# Residuals      640 23.703  0.0370                                

pdf("s6.FC_corr.pdf")
ggplot(data=sumCorr[sumCorr$Func=="go",], aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Function Connectivity: go")
ggplot(data=sumCorr[sumCorr$Func=="kegg",], aes(x=netType, y=correlation, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Correlations of Function Connectivity: kegg")
dev.off()

# FC AUROC
pdf("s6.FC_auroc.pdf")
for(i in c("go","kegg"))
{
    p=ggplot(data=sumCorr[sumCorr$Func==i,], aes(x=netType, y=exp, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Expected Functional Connectivity: ",i))
    print(p)
    p=ggplot(data=sumCorr[sumCorr$Func==i,], aes(x=netType, y=obs, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(paste0("Observed Functional Connectivity: ",i))
    print(p)
}
dev.off()

fit<-lm(obs~pipeline+transformation+netType,data= sumCorr[sumCorr$Func=="go",])
print(summary(a<-aov(fit)))
#                 Df Sum Sq Mean Sq  F value   Pr(>F)    
# pipeline         4 0.0007 0.00018    6.903 1.85e-05 ***
# transformation   1 0.0045 0.00446  167.116  < 2e-16 ***
# netType         10 2.6125 0.26125 9792.012  < 2e-16 ***
# Residuals      754 0.0201 0.00003                      
out <- duncan.test(a,trt="pipeline"); print(out$groups)
#                obs groups
# rsem     0.5747569      a
# hylite   0.5744856     ab
# kallisto 0.5735097     bc
# salmon   0.5727364     cd
# polycat  0.5722126      d
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.5759466      a
# rld      0.5711339      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                         obs groups
# weighted.WGCNA.12 0.6619626      a
# weighted.WGCNA.24 0.6541643      b
# weighted.WGCNA.1  0.6377420      c
# binary.rank.0.95  0.6170664      d
# binary.rank.0.99  0.5815082      e
# binary.rank.0.995 0.5662655      f
# binary.Zscore.1.5 0.5391208      g
# binary.rank.0.999 0.5381444      g
# binary.Zscore.2   0.5105609      h
# binary.Zscore.2.5 0.5018858      i
# binary.Zscore.3   0.5005214      i

fit<-lm(obs~pipeline+transformation+netType,data= sumCorr[sumCorr$Func=="kegg",])
print(summary(a<-aov(fit)))
#                 Df Sum Sq Mean Sq  F value Pr(>F)    
# pipeline         4 0.0007 0.00019    2.188 0.0687 .  
# transformation   1 0.0096 0.00965  113.673 <2e-16 ***
# netType         10 2.6656 0.26656 3140.254 <2e-16 ***
# Residuals      754 0.0640 0.00008                    
out <- duncan.test(a,trt="pipeline"); print(out$groups)
#                obs groups
# rsem     0.5748750      a
# kallisto 0.5746074      a
# hylite   0.5739550     ab
# salmon   0.5728688     ab
# polycat  0.5723335      b
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.5772679      a
# rld      0.5701880      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                         obs groups
# weighted.WGCNA.24 0.6620320      a
# weighted.WGCNA.12 0.6590120      a
# binary.rank.0.95  0.6285811      b
# weighted.WGCNA.1  0.6270494      b
# binary.rank.0.99  0.5811237      c
# binary.rank.0.995 0.5643640      d
# binary.Zscore.1.5 0.5380434      e
# binary.rank.0.999 0.5368526      e
# binary.Zscore.2   0.5112027      f
# binary.Zscore.2.5 0.5022713      g
# binary.Zscore.3   0.5004750      g


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

# Density A vs D vs A-D
library(ggpubr)
library(reshape2)
library(dplyr)
res=melt(data = sumCorr[,c(3:5,7:14)], id.vars = names(sumCorr)[c(10:14)], measure.vars = names(sumCorr)[c(3:5,7:9)],value.name ="auroc")
res$subnetwork = factor(gsub(".*[.]", "", res$variable))
res$dataset = factor(gsub("[.].*","",res$variable))
res$group= paste(res$pipeline,res$transformation, res$netType, sep="_")
# aNOVA
fit<-lm(auroc~pipeline+transformation+netType+Func+dataset+subnetwork,data= res)
print(summary(a<-aov(fit)))
#                  Df Sum Sq Mean Sq   F value  Pr(>F)    
# pipeline          4   0.00    0.00     0.065 0.99233    
# transformation    1   0.06    0.06     9.530 0.00203 ** 
# netType          10  99.31    9.93  1511.992 < 2e-16 ***
# Func              1   0.03    0.03     5.207 0.02252 *  
# dataset           1   0.00    0.00     0.133 0.71533    
# subnetwork        2 248.94  124.47 18950.069 < 2e-16 ***
# Residuals      9220  60.56    0.01                      
out <- duncan.test(a,trt="transformation"); print(out$groups)
#                obs groups
# log2rpkm 0.4482661      a
# rld      0.4430606      b
out <- duncan.test(a,trt="netType"); print(out$groups)
#                       auroc groups
# weighted.WGCNA.12 0.5887069      a
# weighted.WGCNA.24 0.5875610      a
# weighted.WGCNA.1  0.5783241      b
# binary.rank.0.95  0.5449161      c
# binary.rank.0.99  0.4508022      d
# binary.rank.0.995 0.4141949      e
# binary.Zscore.1.5 0.3660475      f
# binary.rank.0.999 0.3639756      f
# binary.Zscore.2   0.3401754      g
# binary.Zscore.2.5 0.3343094      g
# binary.Zscore.3   0.3332842      g
out <- duncan.test(a,trt="subnetwork"); print(out$groups)
#            auroc groups
# D       0.5619642      a
# A       0.5614918      a
# interAD 0.2135342      b



pdf("s6.FC_subAD.pdf")
for(each in unique(res$group)){
    print(each)  # 11x 5 x 2 =110
    for(func in c("go","kegg"))
    {
       df=res[res$group==each&res$Func==func,]
       title=paste0(func," AUCROC: pipeline-",unique(df$pipeline),", transformation-",unique(df$transformation),", netType-",unique(df$netType))  
       # boxplot
       # p=ggplot(df, aes(x=subnetwork, y=auroc)) + geom_boxplot(aes(color = dataset)) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
       #
       # line plot with errorbar
       df.summary2 <- df[,7:9]  %>%
       group_by(subnetwork, dataset) %>%
       summarise( sd = sd(auroc), auroc = mean(auroc) )
       df.summary2
       p=ggplot(df, aes(x=subnetwork, y=auroc)) + geom_jitter( aes(color = dataset), position = position_jitter(0.2)) + geom_line( aes(group = dataset, color = dataset), data = df.summary2) + geom_errorbar( aes(ymin = auroc-sd, ymax = auroc+sd, color = dataset), data = df.summary2, width = 0.1) + scale_color_manual(values = c("#00AFBB", "#E7B800")) + theme_bw() + theme(panel.border=element_blank()) + stat_compare_means(aes(group = dataset),label = "p.signif") + ggtitle(title)
       print(p)
    }
}
dev.off()

# make a summary table of density mean and sd as supplemental table
res.summary <- res  %>%
group_by(subnetwork, dataset, group,Func) %>%
summarise( sd = sd(auroc), density = mean(auroc) )%>%
arrange(Func, group)
res.summary
write.table(res.summary,file = "s6.FC_subAD.txt", sep="\t",row.names=FALSE )



