# limit memory usage
# ulimit -v 1288490188


# limit memory usage
# ulimit -v 1288490188
system("mem")
system("ulimit -v")
# better use bigram for this analysis
library(WGCNA)
source("6.FUN.r")


################################
## Node Connectivity Analysis ##
################################
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


rdatafiles<-grep("R-05-dataInput",list.files(),value=TRUE)
rdatafiles
# "R-05-dataInput.hylite_log2rpkm.RData"
# "R-05-dataInput.hylite_rld.RData"
# "R-05-dataInput.kallisto_log2rpkm.RData"
# "R-05-dataInput.kallisto_rld.RData"
# "R-05-dataInput.polycat_log2rpkm.RData"
# "R-05-dataInput.polycat_rld.RData"
# "R-05-dataInput.rsem_log2rpkm.RData"
# "R-05-dataInput.rsem_rld.RData"
# "R-05-dataInput.salmon_log2rpkm.RData"
# "R-05-dataInput.salmon_rld.RData"
load(rdatafiles[1])
use=colnames(multiExpr[[1]]$data)
print(length(use))
for(i in rdatafiles[-1])
{
    load(i);print(i);
    ids =colnames(multiExpr[[1]]$data)
    print(length(ids))
    use=intersect(use, ids)
}
print(length(use)) #62656

for(file in rdatafiles)
{
    print(file)
    load(file)
    tag<-gsub(".*Input[.]|[.].*","",file)
    for(s in 1:nSets){
        multiExpr[[s]]$data = multiExpr[[s]]$data[,use]
    }
    
    message(paste0("Start dataset: ",tag))
    
    # get expression datasets of exp and obs
    exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    
    # Peason's correlation
    exp.r <- cor(t(exp),method = "pearson",nThread=8)
    obs.r <- cor(t(obs),method = "pearson",nThread=8)
    
    for(net in 1:nrow(netType))
    {
        start<-proc.time()
        # build networks
        if(netType$edge[net]=="binary")
        {
            message(paste0("----------------------------------\n...network type ",netType$net[net],": ", netType$desc[net]))
            exp.net <- build_binary_net(exp.r, metric = netType$metric[net], cutoff = netType$cutoff[net] ,subAD=FALSE )
            obs.net <- build_binary_net(obs.r, metric = netType$metric[net], cutoff = netType$cutoff[net], subAD=FALSE )
            
        } else if(netType$edge[net]=="weighted")
        {
            message(paste0("----------------------------------\n...network type ",netType$net[net],": ", netType$desc[net]))
            exp.net <- build_wgcna_net(exp.r, power=netType$cutoff[net], subAD=FALSE )
            obs.net <- build_wgcna_net(obs.r, power=netType$cutoff[net], subAD=FALSE )
            }
        }
        
        # node connectivity
        exp.k<-rowSums(exp.net)-1
        obs.k<-rowSums(obs.net)-1
        
        ks <- data.frame(exp=exp.k, obs=obs.k)
        adjs <- list (exp=exp.net, obs=obs.net)
        
        # density
        ## calculate density
        size <- nrow(exp.r)
        GeneID <-rownames(exp.r)
        at<-grep("a$",GeneID)
        dt<-grep("d$",GeneID)
        size.at<-length(at)
        size.dt<-length(dt)
        # Density = sum(Connectivity)/(Size * (Size - 1))
        density <- apply(ks,2,function(k){sum(k)/(size * (size - 1))} )
        at.density <- lapply(adjs,function(adj){(sum(adj[at,at])-size.at)/(size.at*(size.at-1))} )
        dt.density <- lapply(adjs,function(adj){(sum(adj[dt,dt])-size.dt)/(size.dt*(size.dt-1))} )
        ad.density <- lapply(adjs,function(adj){mean(adj[at,dt])})
        rr<- data.frame(net=names(ks),method=flag, netType= netType$desc[net], density=unlist(density), at.density=unlist(at.density), dt.density=unlist(dt.density), ad.density=unlist(ad.density))
        print(rr)
        if(!exists("sumDensity")){sumDensity<-rr}else{sumDensity<-rbind(sumDensity,rr)}
        
        ## corrlation of connectivity
        ids<-unique(gsub("a$|d$","",GeneID))
        a<-paste0(ids,"a")
        d<-paste0(ids,"d")
        corr <- cor(exp.k,obs.k)
        corr.expvsobs.A <-cor(exp.k[a],obs.k[a], use = 'pairwise.complete.obs')
        corr.expvsobs.D <-cor(exp.k[d],obs.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.exp <- cor(exp.k[a],exp.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.obs <- cor(obs.k[a],obs.k[d], use = 'pairwise.complete.obs')
        tt<-data.frame(correlation=c(corr,corr.expvsobs.A, corr.expvsobs.D, corr.AvsD.exp, corr.AvsD.obs), compr=c("expvsobs","expvsobs.A", "epvsobs.D", "AvsD.exp", "AvsD.obs"), mapping=tag, netType= netType$desc[net])
        print(tt)
        if(!exists("sumCorr")){sumCorr<-tt}else{sumCorr<-rbind(sumCorr,tt)}
    }
    save(sumCorr,sumDensity, file="s6.NC.rdata")
    gc()
}

--book
save(sumCorr,sumDensity, file="s6.NC.rdata")
        

###########
## PLOTS ##
###########
library(ggplot2)
theme_update(plot.title = element_text(hjust = 0.5),legend.position="bottom", axis.text.x=element_text(angle=90, hjust=1))

NC<-sumCorr
NC$homoeoAssign<-gsub("_.*","", NC$mapping)
NC$norm<-gsub(".*_","", NC$mapping)

pdf("s6.NC.plotCorr.pdf")
ggplot(data=NC, aes(x=homoeoAssign, y=correlation, color = norm)) + geom_boxplot() + facet_grid(.~compr)+ggtitle("Node connectivity: correlation of k")
ggplot(data=NC, aes(x=netType, y=correlation, color = norm)) + geom_boxplot() + facet_grid(.~compr)+ggtitle("Node connectivity: correlation of k")
ggplot(data=NC, aes(x=homoeoAssign, y=correlation, color = norm)) + geom_boxplot() + facet_grid(facet=compr~netType, labeller=label_both)
dev.off()

summary(aov(correlation~compr+homoeoAssign+norm+netType,data=NC))
# "
Df Sum Sq Mean Sq F value               Pr(>F)
compr          4 12.427  3.1069 483.995 < 0.0000000000000002 ***
homoeoAssign   2  0.008  0.0040   0.620              0.53958
norm           1  0.393  0.3932  61.255     0.00000000000118 ***
netType        4  0.127  0.0317   4.935              0.00095 ***
Residuals    138  0.886  0.0064
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# '


#############
## for NYC ##
#############
load("NC.rdata")
x<-NC[NC$compr=="expvsobs",]
x$homoeoAssign<-factor(x$homoeoAssign,levels=c("polycat","hylite","rsem" ))

fit<-lm(correlation~homoeoAssign+norm+netType,data=x)
summary(fit)
# correlation significanntly lower in Z network type, rpkm is better than rld
x$netType<-factor(x$netType, levels=c("binary.Zscore.2","binary.rank.0.99","weighted.WGCNA.12" ,"weighted.WGCNA.24", "weighted.WGCNA.1" ))

pdf("NC.obsvsexp.pdf")
ggplot(data=x, aes(x=netType, y=correlation, fill = homoeoAssign)) + geom_boxplot(lwd=0.2,position=position_dodge(0.8))+theme_bw()+theme(panel.border=element_blank(),axis.line.x=element_line(),legend.position="bottom") +scale_fill_manual(values=c( "#b4b2bc", "black","#22aad6"))

dev.off()



