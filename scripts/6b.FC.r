# limit memory usage
# ulimit -v 1288490188
# system("mem")
# system("ulimit -v")
# better use bigram for this analysis

library(EGAD)
library(WGCNA)

## generate functional annotation categories
#source("6.prepFunctionCategory.r")
load("GOnKEGGnGLs.Rdata")

## load functions
source("6.FUN.r")


######################################
## Functional Connectivity Analysis ##
######################################


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


rdatafiles<-list.files(pattern="R-05-dataInput")
rdatafiles
#  [1] "R-05-dataInput.bowtie_log2rpkm.RData"
#  [2] "R-05-dataInput.bowtie_rld.RData"
#  [3] "R-05-dataInput.eaglerc_log2rpkm.RData"
#  [4] "R-05-dataInput.eaglerc_rld.RData"
#  [5] "R-05-dataInput.hylite_log2rpkm.RData"
#  [6] "R-05-dataInput.hylite_rld.RData"
#  [7] "R-05-dataInput.kallisto_log2rpkm.RData"
#  [8] "R-05-dataInput.kallisto_rld.RData"
#  [9] "R-05-dataInput.polycat_log2rpkm.RData"
# [10] "R-05-dataInput.polycat_rld.RData"
# [11] "R-05-dataInput.rsem_log2rpkm.RData"
# [12] "R-05-dataInput.rsem_rld.RData"
# [13] "R-05-dataInput.salmon_log2rpkm.RData"
# [14] "R-05-dataInput.salmon_rld.RData"
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
print(length(use)) # 52000
# how many in pairs or solo
y=table(gsub("a$|d$","",use))
table(y==2) # 25907 pairs and 186 solo
idsP =names(y)[y==2]
netSize=5200

rm(sumCorr)
# calculating whole network r and adjacency require too much memory
# randomly choose 10000 genes to obstain corr, repeat 10 times

for(p in 1:10){
    message(paste0("Start permutation: ",p))
    sub=sample(1:length(idsP),netSize/2, replace=F)
    for(file in rdatafiles)
    {
        start<-proc.time()
        print(file)
        load(file)
        tag<-gsub(".*Input[.]|[.].*","",file)
        message(paste0("Start dataset: ",tag))
        
        # get expression datasets of exp and obs
        exp<-as.data.frame(multiExpr[[which(shortLabels=="A2D5")]]$data[,use])
        obs<-as.data.frame(multiExpr[[which(shortLabels=="ADs")]]$data[,use])
        
        # subsample
        indA = paste0(idsP[sub],"a")
        indD = paste0(idsP[sub],"d")
        id=c(indA, indD)
        
        # Peason's correlation
        exp.r <- cor(exp[,id],method = "pearson",nThread=8)
        obs.r <- cor(obs[,id],method = "pearson",nThread=8)
        
        for(net in 1:nrow(netType))
        {
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
            
            # functional connectivity analysis
            for(func in c("kegg","go"))
            {
                message(paste0("......",func))
                annotation <- get(paste0(func,".AD"))
                labels <- prep_annot(t(exp[,id]), annotation, min = 10, max = 500)
                # get fc
                message(paste0("......",length(indA)," At vs ",length(indD)," Dt genes."))
                exp.fc <- list(all=neighbor_voting(labels, exp.net, output = 'AUROC')[,"auc"], A=neighbor_voting(labels, exp.net[indA,indA], output = 'AUROC')[,"auc"], D=neighbor_voting(labels, exp.net[indD,indD], output = 'AUROC')[,"auc"], interAD=neighbor_voting(labels, exp.net[indA,indD], output = 'AUROC')[,"auc"])
                obs.fc <- list(all=neighbor_voting(labels, obs.net, output = 'AUROC')[,"auc"], A=neighbor_voting(labels, obs.net[indA,indA], output = 'AUROC')[,"auc"], D=neighbor_voting(labels, obs.net[indD,indD], output = 'AUROC')[,"auc"], interAD=neighbor_voting(labels, obs.net[indA,indD], output = 'AUROC')[,"auc"])
                # overall FC correlation
                corr =cor(exp.fc$all, obs.fc$all, use = "pairwise.complete.obs")
                tt<-data.frame(t(c(corr, unlist(lapply(obs.fc,function(x)mean(x,na.rm=TRUE))), unlist(lapply(exp.fc,function(x)mean(x,na.rm=TRUE))))))
                names(tt)=c("correlation", "obs", "obs.A", "obs.D", "obs.interAD", "exp","exp.A","exp.D","exp.interAD")
                tt$pipeline=gsub("_.*","",tag)
                tt$transformation=gsub(".*_","",tag)
                tt$netType= netType$desc[net]
                tt$permutation = p
                tt$Func=func
                print(tt)
                if(!exists("sumCorr")){sumCorr<-tt}else{sumCorr<-rbind(sumCorr,tt)}
            }
        }
        # report time used for each rdatafile at each iteration
        ptm<-proc.time()-start
        message(paste0("...running time: ", round(ptm[3]/60,1), " minutes."))
        gc()
        save(sumCorr,file="s6.FC.rdata")
    }
}
# about 30 min x 10 pipeline_norm x 10 interation x= 3000 min, about 50 hrs

q()
n
