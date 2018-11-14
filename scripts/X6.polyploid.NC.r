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

###################################
## Network Connectivity Analysis ##
###################################
system("mem")
system("ulimit -v")

library(WGCNA)
options(scipen=999) # disable scientifc format

rdatafiles<-grep("R-02-dataInput.*DP",list.files(),value=TRUE)
rdatafiles
# "R-02-dataInput.polycat_rld.DP.RData"
# "R-02-dataInput.polycat_rpkm.DP.RData"

netType<-data.frame( net=1:5, edge=c("binary","binary","weighted","weighted","weighted"),metric=c("rank","Zscore","WGCNA","WGCNA","WGCNA"),cutoff=c(0.99,2,1,24,12))
netType$desc <- paste(netType$edge,netType$metric, netType$cutoff,sep=".")
netType
# net     edge metric cutoff              desc
# 1   1   binary   rank   0.99  binary.rank.0.99
# 2   2   binary Zscore   2.00   binary.Zscore.2
# 3   3 weighted  WGCNA   1.00  weighted.WGCNA.1
# 4   4 weighted  WGCNA  24.00 weighted.WGCNA.24
# 5   5 weighted  WGCNA  12.00 weighted.WGCNA.12

rm(sumNC)
for(file in rdatafiles)
{
    # file = "R-02-dataInput.polycat_rld.DP.RData"
    print(file)
    flag <- gsub("R-02-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    message(paste0("Start dataset: ",flag))
    
    # get expression datasets of exp and obs
    A2D5<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    ADs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    # TM1<-as.data.frame(t(multiExpr[[which(shortLabels=="TM1")]]$data))
    Yuc<-as.data.frame(t(multiExpr[[which(shortLabels=="Yuc")]]$data))
    
    # functional connectivity analysis

        i=4
        message(paste0("----------------------------------\n...network type ",netType$net[i],": ", netType$desc[i]))
        
        A2D5.net <- adjacency(t(A2D5), power=netType$cutoff[i], type="signed" )
        A2D5.k<-rowSums(A2D5.net)-1
        
        ADs.net <- adjacency(t(ADs), power=netType$cutoff[i], type="signed" )
        ADs.k<-rowSums(ADs.net)-1
        
        #TM1.net <- adjacency(t(TM1), power=netType$cutoff[i], type="signed" )
        #TM1.k<-rowSums(TM1.net)-1
        
        Yuc.net <- adjacency(t(Yuc), power=netType$cutoff[i], type="signed" )
        Yuc.k<-rowSums(Yuc.net)-1

        ks <- data.frame(A2D5=A2D5.k, ADs=ADs.k, Yuc=Yuc.k)
        adjs <- list (A2D5=A2D5.net, ADs=ADs.net, Yuc=Yuc.net)
        
        ## calculate density
        size <- checkSets(multiExpr)$nGenes
        GeneID <-colnames(multiExpr[[1]]$data)
        at<-grep("a$",GeneID)
        dt<-grep("d$",GeneID)
        size.at<-length(at)
        size.dt<-length(dt)
        # Density = sum(Connectivity)/(Size * (Size - 1))
        density <- apply(ks,2,function(k){sum(k)/(size * (size - 1))} )
        at.density <- lapply(adjs,function(adj){(sum(adj[at,at])-size.at)/(size.at*(size.at-1))} )
        dt.density <- lapply(adjs,function(adj){(sum(adj[dt,dt])-size.dt)/(size.dt*(size.dt-1))} )
        ad.density <- lapply(adjs,function(adj){mean(adj[at,dt])})
        resSum<- data.frame(net=names(ks),density=unlist(density), at.density=unlist(at.density), dt.density=unlist(dt.density), ad.density=unlist(ad.density))
        print(resSum)
        save(ks,at,dt,resSum, file=paste0("s6.NC.",netType$desc[i],".",flag,".rdata"))
        
        ## plot node connectivity corrlations
        ids<-unique(gsub("a$|d$","",GeneID))
        a<-paste0(ids,"a")
        d<-paste0(ids,"d")
        pdf(paste0("s6.NC.",netType$desc[i],".",flag,".pdf"))
        # between species
        pairs(ks,lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))
        pairs(ks[a,],lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))
        pairs(ks[d,],lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))
        # At vs Dt
        resSum$AvsD.k.cor <- apply(ks,2,function(k)cor(k[a],k[d], use='pairwise.complete.obs'))
        dev.off()
}

---book
        # same as WGCNA function, which doesn't keep gene ID though
        # k<-softConnectivity(t(exp), power=netType$cutoff[i], type="signed")

pairs(k,lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))
pairs(k[a,],lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))
pairs(k[d,],lower.panel=panel.smooth, upper.panel=panel.cor, pch=".", col = rgb(0, 0, 0, 0.05))

cor(k[a,"A2D5"],k[d,"A2D5"], use = 'pairwise.complete.obs')

        
        # A vs D
        corr.AvsD.Yuc <- cor(Yuc.k[a],Yuc.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.TM1 <- cor(TM1.k[a],TM1.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.A2D5 <- cor(A2D5.k[a],A2D5.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.ADs <- cor(ADs.k[a],ADs.k[d], use = 'pairwise.complete.obs')
        
        # polyploid vs diploid
        corr.TM1vsA2D5 <- cor(TM1.k,A2D5.k, use = 'pairwise.complete.obs')
        corr.TM1vsADs  <- cor(TM1.k,ADs.k, use = 'pairwise.complete.obs')
        corr.YucvsA2D5 <- cor(Yuc.k,A2D5.k, use = 'pairwise.complete.obs')
        corr.YucvsADs  <- cor(Yuc.k,ADs.k, use = 'pairwise.complete.obs')
        # domestication
        corr.TM1vsYuc  <- cor(TM1.k,Yuc.k, use = 'pairwise.complete.obs')
        # diploid exp vs obs
        corr.ADsvsA2D5 <- cor(ADs.k,A2D5.k, use = 'pairwise.complete.obs')


# book

corr <- cor(exp.k,obs.k)
        corr.expvsobs.A <-cor(exp.k[a],obs.k[a], use = 'pairwise.complete.obs')
        corr.expvsobs.D <-cor(exp.k[d],obs.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.exp <- cor(exp.k[a],exp.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.obs <- cor(obs.k[a],obs.k[d], use = 'pairwise.complete.obs')
        
        tt<-data.frame(correlation=c(corr,corr.expvsobs.A, corr.expvsobs.D, corr.AvsD.exp, corr.AvsD.obs), compr=c("expvsobs","expvsobs.A", "epvsobs.D", "AvsD.exp", "AvsD.obs"), mapping=flag, netType= netType$desc[i])
        if(!exists("sumNC")){sumNC<-tt}else{sumNC<-rbind(sumNC,tt)}
    }
    gc()
}

NC_wgcna<-sumNC


###########
## wgcna ##
###########
library(WGCNA)

netType<-data.frame(net=c(1,4:5),edge="weighted",metric="WGCNA",cutoff=c(1, 24,12))
netType$desc <- paste(netType$edge,netType$metric, netType$cutoff,sep=".")
netType
# net     edge metric cutoff              desc
#   1 weighted  WGCNA      1  weighted.WGCNA.1
#   4 weighted  WGCNA     24 weighted.WGCNA.24
#   5 weighted  WGCNA     12 weighted.WGCNA.12


rm(sumNC)
for(file in rdatafiles)
{
    # file = "R-02-dataInput.polycat_rld.RData"
    print(file)
    flag <- gsub("R-02-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    message(paste0("Start dataset: ",flag))
    
    # get expression datasets of exp and obs
    exp<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    obs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    

    # functional connectivity analysis
    for(i in 1:nrow(netType))
    {
        message(paste0("----------------------------------\n...network type ",netType$net[i],": ", netType$desc[i]))
        
        exp.net <- adjacency(t(exp), power=netType$cutoff[i], type="signed" )
        exp.k<-rowSums(exp.net)-1
        
        obs.net <- adjacency(t(obs), power=netType$cutoff[i], type="signed" )
        obs.k<-rowSums(obs.net)-1
        
        # same as WGCNA function, which doesn't keep gene ID though
        # k<-softConnectivity(t(exp), power=netType$cutoff[i], type="signed")
        
        ids<-unique(gsub("a$|d$","",names(exp.k)))
        a<-paste0(ids,"a")
        d<-paste0(ids,"d")
        
        corr <- cor(exp.k,obs.k)
        corr.expvsobs.A <-cor(exp.k[a],obs.k[a], use = 'pairwise.complete.obs')
        corr.expvsobs.D <-cor(exp.k[d],obs.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.exp <- cor(exp.k[a],exp.k[d], use = 'pairwise.complete.obs')
        corr.AvsD.obs <- cor(obs.k[a],obs.k[d], use = 'pairwise.complete.obs')
        
        tt<-data.frame(correlation=c(corr,corr.expvsobs.A, corr.expvsobs.D, corr.AvsD.exp, corr.AvsD.obs), compr=c("expvsobs","expvsobs.A", "epvsobs.D", "AvsD.exp", "AvsD.obs"), mapping=flag, netType= netType$desc[i])
        if(!exists("sumNC")){sumNC<-tt}else{sumNC<-rbind(sumNC,tt)}
    }
    gc()
}

NC_wgcna<-sumNC

###########
## PLOTS ##
###########
library(ggplot2)
theme_update(plot.title = element_text(hjust = 0.5),legend.position="bottom", axis.text.x=element_text(angle=90, hjust=1))

NC<-rbind(NC_bi,NC_wgcna)
NC$homoeoAssign<-gsub("_.*","", NC$mapping)
NC$norm<-gsub(".*_","", NC$mapping)
save(NC, file="NC.rdata")

pdf("plotNC.pdf")
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



