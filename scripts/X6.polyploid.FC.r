######################################
## Functional connectivity analysis ##
######################################
system("mem")
system("ulimit -v")

library(EGAD)
library(WGCNA)
enableWGCNAThreads(nThreads = 8)
options(scipen=999) # disable scientifc format

source("6.FC.FUN.r")
load("GOnKEGGnGLs.Rdata")

rdatafiles<-grep("R-02-dataInput.*DP",list.files(),value=TRUE)
rdatafiles
# "R-02-dataInput.polycat_rld.DP.RData"
# "R-02-dataInput.polycat_rpkm.DP.RData"

options(scipen=999) # disable scientifc format
netType<-data.frame( net=1:5, edge=c("binary","binary","weighted","weighted","weighted"),metric=c("rank","Zscore","WGCNA","WGCNA","WGCNA"),cutoff=c(0.99,2,1,24,12))
netType$desc <- paste(netType$edge,netType$metric, netType$cutoff,sep=".")
netType
# net     edge metric cutoff               desc
#   1   binary     rank   0.99     binary.rank.0.99
#   2   binary   Zscore   2.00      binary.Zscore.2
#   3 weighted    WGCNA   1.00     weighted.WGCNA.1
#   4 weighted    WGCNA  24.00    weighted.WGCNA.24
#   5 weighted    WGCNA  12.00    weighted.WGCNA.12

for(file in rdatafiles)
{
    # file = "R-02-dataInput.polycat_rld.DP.RData"
    print(file)
    flag <- gsub("R-02-dataInput.|.RData","",file)
    fn<-load(file) # # multiExpr,nSets, setLabels, shortLabels
    message(paste0("Start dataset: ",flag))
    
    # get expression datasets
    A2D5<-as.data.frame(t(multiExpr[[which(shortLabels=="A2D5")]]$data))
    ADs<-as.data.frame(t(multiExpr[[which(shortLabels=="ADs")]]$data))
    Yuc<-as.data.frame(t(multiExpr[[which(shortLabels=="Yuc")]]$data))
    ids<-rownames(A2D5)
    
    net=4
    message(paste0("----------------------------------\n...network type ",netType$net[net],": ", netType$desc[net]))
    
    # WGCNA adjacency
    Yuc.A.net <- adjacency(t(Yuc[grep("a$",ids),]),power=netType$cutoff[net], type="signed")
    Yuc.D.net <- adjacency(t(Yuc[grep("d$",ids),]),power=netType$cutoff[net], type="signed")
    ADs.A.net <- adjacency(t(ADs[grep("a$",ids),]),power=netType$cutoff[net], type="signed")
    ADs.D.net <- adjacency(t(ADs[grep("d$",ids),]),power=netType$cutoff[net], type="signed")
    A2D5.A.net <- adjacency(t(A2D5[grep("a$",ids),]),power=netType$cutoff[net], type="signed")
    A2D5.D.net <- adjacency(t(A2D5[grep("d$",ids),]),power=netType$cutoff[net], type="signed")
    
    for(func in c("kegg","go","fa","ft")){
        message(paste0("......",func))
        annotation <- get(paste0(func,".AD"))
        labels <- prep_annot(A2D5, annotation, min = 9, max = 500)
        
        message(paste0("......",ncol(A2D5.A.net)," At vs ",ncol(A2D5.D.net)," Dt genes."))
        Yuc.fc <- list(A =neighbor_voting(labels, Yuc.A.net, output = 'AUROC'), D=neighbor_voting(labels, Yuc.D.net, output = 'AUROC'))
        ADs.fc <- list(A =neighbor_voting(labels, ADs.A.net, output = 'AUROC'), D=neighbor_voting(labels, ADs.D.net, output = 'AUROC'))
        A2D5.fc <- list(A =neighbor_voting(labels, A2D5.A.net, output = 'AUROC'), D=neighbor_voting(labels, A2D5.D.net, output = 'AUROC'))
        save(A2D5.fc,ADs.fc,Yuc.fc,file=paste("s6.FC",flag,netType$desc[net],func,"DP.rdata",sep="."))
    }
    gc()
}


## PLOTS ##
pdf("s6.FC.AvsD.DP.pdf")
for(func in c("kegg","go","fa","ft")){
    load(paste(flag,netType$desc[i],func,"rdata",sep="."))
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    plot_auc_compare(ADs.fc$A[,1], ADs.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="ADs",pch=".")
    plot_auc_compare(A2D5.fc$A[,1], A2D5.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="A2D5",pch=".")
    plot_auc_compare(Yuc.fc$A[,1], Yuc.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="Yuc",pch=".")
    plot_auc_compare(TM1.fc$A[,1], TM1.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="TM1",pch=".")
    mtext(func, outer = TRUE, cex = 2)
    
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    plot_auc_compare(ADs.fc$A[,1], TM1.fc$A[,1], xlab ="AUC - A2",ylab="AUC - At", title="ADs",pch=".")
    plot_auc_compare(ADs.fc$D[,1], TM1.fc$D[,1], xlab ="AUC -D5",ylab="AUC - Dt", title="ADs",pch=".")
    plot_auc_compare(A2D5.fc$A[,1], TM1.fc$A[,1], xlab ="AUC - A2",ylab="AUC - At", title="ADs",pch=".")
    plot_auc_compare(A2D5.fc$D[,1], TM1.fc$D[,1], xlab ="AUC -D5",ylab="AUC - Dt", title="A2D5",pch=".")
    
    
    plot_auc_compare(Yuc.fc$A[,1], Yuc.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="Yuc",pch=".")
    plot_auc_compare(TM1.fc$A[,1], TM1.fc$D[,1], xlab ="AUC - At",ylab="AUC - Dt", title="TM1",pch=".")
    mtext(func, outer = TRUE, cex = 2)
    
    diff<-data.frame(A2=ADs.fc$A[,1], At=TM1.fc$A[,1],D5=ADs.fc$D[,1], Dt=TM1.fc$D[,1] )
    diff$diffA = diff$At-diff$A2
    diff$diffD = diff$Dt-diff$D5
}
dev.off()

rm(sumFC)
for(func in c("kegg","go","fa")){
    load(paste(flag,netType$desc[i],func,"rdata",sep="."))
    for(genome in c("A2D5","ADs","Yuc","TM1"))
    {
        term<-paste0(genome, ".fc")
        df<-fc2df(get(term))
        df$species <-genome
        df$annotation<-func
        if(exists("sumFC"))
        { sumFC<- rbind(sumFC,df)
        }else{sumFC<-df}
    }
}
save(sumFC,file="NYC/sumFC.rdata")
fit<-lm(auc~genome+species+annotation,data=sumFC)

