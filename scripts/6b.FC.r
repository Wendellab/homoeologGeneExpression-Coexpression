# limit memory usage
# ulimit -v 1288490188
system("mem")
system("ulimit -v")
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
            exp.net <- build_binary_net(exp.r, metric = netType$metric[net], cutoff = netType$cutoff[net] ,subAD=TRUE )
            obs.net <- build_binary_net(obs.r, metric = netType$metric[net], cutoff = netType$cutoff[net], subAD=TRUE )
            
        } else if(netType$edge[net]=="weight")
        {
            message(paste0("----------------------------------\n...network type ",netType$net[net],": ", netType$desc[net]))
            exp.net <- build_wgcna_net(exp.r, power=netType$cutoff[net], subAD=TRUE )
            obs.net <- build_wgcna_net(obs.r, power=netType$cutoff[net], subAD=TRUE )
        }
        
        # functional connectivity analysis
        pdf(paste0("s6.FC/",tag,".",netType$desc[net],".pdf"))
        for(func in c("kegg","go","fa","ft")){
            message(paste0("......",func))
            annotation <- get(paste0(func,".AD"))
            labels <- prep_annot(exp, annotation, min = 10, max = 500)
            # get fc
            message(paste0("......",ncol(exp.net$A)," At vs ",ncol(exp.net$D)," Dt genes."))
            exp.fc <- list(A =neighbor_voting(labels, exp.net$A, output = 'AUROC'), D=neighbor_voting(labels, exp.net$D, output = 'AUROC'))
            obs.fc <- list(A =neighbor_voting(labels, obs.net$A, output = 'AUROC'), D=neighbor_voting(labels, obs.net$D, output = 'AUROC'))
            ## if network generated not for subgenomes
            # exp.fc <- run_nv_homoeo(exp.net, annotation)
            # obs.fc <- run_nv_homoeo(obs.net, annotation)
            
            # plot At vs Dt and exp vs obs
            compare_nvs(obs.fc, exp.fc, title=paste0("s6.FC/",tag,".",netType$desc[net],".",func), save.rdata=TRUE)
        }
        dev.off()
        
        # report time passed
        ptm<-proc.time()-start
        message(paste0("...running time: ", round(ptm[3]/60,1), " minutes."))
        
        gc()
    }
}
