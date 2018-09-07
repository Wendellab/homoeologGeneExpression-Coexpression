polycat <- load("R-01-polycatNetworkDatasets.RData")
raw.polycat <- networks
rld.polycat <- networks.rld
rpkm.polycat<- networks.rpkm

rsem    <- load("R-01-rsemNetworkDatasets.RData")
raw.rsem <- networks
rld.rsem <- networks.rld
rpkm.rsem<- networks.rpkm

hylite <- load("R-01-hyliteNetworkDatasets.RData")
raw.hylite <- networks
rld.hylite <- networks.rld
rpkm.hylite<- networks.rpkm

salmon <- load("R-01-salmonNetworkDatasets.RData")
raw.salmon <- networks
rld.salmon <- networks.rld
tpm.salmon<- networks.tpm

kallisto <- load("R-01-kallistoNetworkDatasets.RData")
raw.kallisto <- networks
rld.kallisto <- networks.rld
tpm.kallisto<- networks.tpm

# function to turn each network list into a big data frame with columns representing networks
library(ggplot2)
list2dataframe<-function(list)
{
    nn<- names(list);
    df<-cbind(as.numeric(as.matrix(list[[1]])), as.numeric(as.matrix(list[[2]])))
    if(length(nn)>2){
        for(i in 3:length(nn)){df<-cbind(df, as.numeric(as.matrix(list[[i]])))}
    }
    df<-as.data.frame(df)
    names(df)<-nn
    return(df)
}



# Plot grouping for each network list, visualizing distances among A2D5, ADs, A2D5.tech, ADs.ncorrect
input = c(ls(pattern=".polycat"),ls(pattern=".hylite"),ls(pattern=".rsem"),ls(pattern=".salmon"), ls(pattern=".kallisto") )
print(input)
pdf("s1.plotVariance.A2D5vsADs.pdf")
for(i in input){
    net.dat <- list2dataframe(get(i))
    # plot(hclust(net.dat))
    
    pca = prcomp(t(net.dat))
    dat = as.data.frame(pca$x)
    dat$group= rownames(dat)
    print( ggplot(aes(PC1, PC2, color=group),data=dat) + geom_point() + ggtitle(i) )
    
    # pairs(net.dat, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, main=paste0("Expression Profile Scatterplot Matrix",i))
}
# so empir corrected profiles are more similiar to uncorrected rpkm than theoretical eflen corrected profiles??
dev.off()
