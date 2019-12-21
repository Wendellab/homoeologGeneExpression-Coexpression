library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads(n=10);


##############################
## Prep multiExpr for WGCNA ##
##############################

filterMultiExpr<-function(multiExpr,nSets)
{
    gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
    print(gsg$allOK)
    if (!gsg$allOK)
    {
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
        for (set in 1:nSets)
        {
            if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
            # Remove the offending genes and samples
            multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
    }
    # Excluding genes from the calculation due to too many missing samples or zero variance.
    # Update exprSize
    print(checkSets(multiExpr))
    return(multiExpr)
}

rfiles<-grep("R-01-.*NetworkDatasets.RData",list.files(),value=TRUE)
rfiles
# [1] "R-01-hyliteNetworkDatasets.RData"
# [2] "R-01-kallistoNetworkDatasets.RData"
# [3] "R-01-polycatNetworkDatasets.RData"
# [4] "R-01-rsemNetworkDatasets.RData"
# [5] "R-01-salmonNetworkDatasets.RData"

pdf(file = "s5.multiExpr.clustering.pdf", width = 12, height = 12)
for(i in rfiles){
    # i<-rfiles[1]
    print(i)
    flag<-gsub("R-01-|NetworkDatasets.RData","",i)
    lnames<-load(i)
    print(lnames)
    
    ######## rlog ###########
    # For easier labeling of plots, create a vector holding descriptive names of the two sets.
    shortLabels = names(networks.rld) # "A2D5"         "A2D5.tech"    "ADs"
    nSets = length(shortLabels)
    multiExpr = vector(mode = "list", length = nSets)
    for(k in 1:nSets) {
        multiExpr[[k]] = list(data = t(networks.rld[[k]])) }
    # Check that the data has the correct format for many functions operating on multiple sets:
    print(checkSets(multiExpr))
    # Check that all genes and samples have sufficiently low numbers of missing values.
    multiExpr<-filterMultiExpr(multiExpr,nSets)
    # save datInput for wgcna
    print(paste0("Save file: R-05-dataInput.",flag,"_rld.RData"))
    save(multiExpr,nSets, shortLabels, file = paste0("R-05-dataInput.",flag,"_rld.RData"))
    # plot clustering
    par(mfrow=c(nSets,1))
    par(mar = c(0, 4, 2, 0))
    sampleTrees = list()
    for (set in 1:nSets)
    {
        sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
    }
    for (set in 1:nSets)
    plot(sampleTrees[[set]], main = paste0(flag,"_rld: ", shortLabels[set]),
    xlab="", sub="", cex = 0.7);
    
    ######## rpkm ###########
    # For easier labeling of plots, create a vector holding descriptive names of the two sets.
    if(flag%in%c("salmon","kallisto")){
        networks.rpkm = networks.tpm; message("tpm")
    }
    shortLabels = names(networks.rpkm) # "A2D5"         "A2D5.tech"    "ADs"
    nSets = length(shortLabels)
    multiExpr = vector(mode = "list", length = nSets)
    for(k in 1:nSets) {
        multiExpr[[k]] = list(data = t(log2(networks.rpkm[[k]]+1))) }
    # Check that the data has the correct format for many functions operating on multiple sets:
    print(checkSets(multiExpr))
    # Check that all genes and samples have sufficiently low numbers of missing values.
    multiExpr<-filterMultiExpr(multiExpr,nSets)
    # save datInput for wgcna
    print(paste0("Save file: R-05-dataInput.",flag,"_log2rpkm.RData"))
    save(multiExpr,nSets, shortLabels, file = paste0("R-05-dataInput.",flag,"_log2rpkm.RData"))
    # plot clustering
    par(mfrow=c(nSets,1))
    par(mar = c(0, 4, 2, 0))
    sampleTrees = list()
    for (set in 1:nSets)
    {
        sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
    }
    for (set in 1:nSets)
    plot(sampleTrees[[set]], main = paste0(flag,"_rpkm: ", shortLabels[set]),
    xlab="", sub="", cex = 0.7);
}
dev.off()


############################
## Choose soft thresholds ##
############################

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
    
    # Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
    # types<-c("unsigned", "signed", "signed hybrid")
    type <- "signed"
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type, blockSize=10000)[[2]])      }
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste0("s5.wgcna.choosePower_", tag,".pdf") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
        if (col==1)
        {
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
            labels=powers,cex=cex1,col=colors[set]);
        } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
        if (col==1)
        {
            legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
        } else
        legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    dev.off()
    assign(paste0("powerTables.",tag),powerTables)
}

# blockSize(65830, rectangularBlocks = TRUE) can calculates suitable block size for Biocrunch
# I just use 10000 here

#### save power tables
pts<-grep("powerTables.",ls(),value=TRUE)
save( list=pts, file = "R-05-chooseSoftThreshold.Rdata")

# Inspect "s2.ChooseSoftThresholdPower_???.pdf", 0.8 hard to reach, use 0.7 instead
load("R-05-chooseSoftThreshold.Rdata")->pts
for(pt in pts)
{
    print(pt)
    power<-unlist(lapply(get(pt), function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.7 & x$slope<0][1]} ) )
    print(power)
}
# "powerTables.hylite_log2rpkm"
# 28 24
# "powerTables.hylite_rld"
# 28 24
# "powerTables.kallisto_log2rpkm"
# 24 24 24
# "powerTables.kallisto_rld"
# 20 24 24
# "powerTables.polycat_log2rpkm"
# 22 22 22 22
# "powerTables.polycat_rld"
# 22 22 22 24
# "powerTables.rsem_log2rpkm"
# 20 20 22
# "powerTables.rsem_rld"
# 22 24 22
# "powerTables.salmon_log2rpkm"
# 22 22 22
# "powerTables.salmon_rld"
# 22 20 24

# ower=24 looks good for all
# check connectivity
for(pt in pts)
{
    print(pt)
    k<-unlist(lapply(get(pt), function(x){ x<-x$data; x[x$Power==24,5:7]} ) )
    print(k)
}
