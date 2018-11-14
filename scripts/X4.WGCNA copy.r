library(WGCNA);
library(RColorBrewer)
library(flashClust);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads();


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

rfiles<-grep("R-01-polycat.*Network",list.files(),value=TRUE)
rfiles
# "R-01-polycatNetworkDatasets.RData"
# "R-01-polycatNetworkDatasets.Yuc.RData"

lnames<-load("R-01-polycatNetworkDatasets.Yuc.RData")
print(lnames)
Yuc<-networks$Yuc
Yuc.rld<-networks.rld$Yuc
Yuc.rpkm<-networks.rpkm$Yuc

lnames<-load("R-01-polycatNetworkDatasets.RData")
print(lnames)
Yuc->networks$Yuc
Yuc.rld->networks.rld$Yuc
Yuc.rpkm->networks.rpkm$Yuc

use<-c("A2D5","ADs","Yuc")
networks<-networks[use]
networks.rld<-networks.rld[use]
networks.rpkm<-networks.rpkm[use]

pdf(file = "s4.wgcna.clustering.DP.pdf", width = 12, height = 12)

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
    print(paste0("Save file: R-02-dataInput.",flag,"_rld.DP.RData"))
    save(multiExpr,nSets, shortLabels, file = paste0("R-02-dataInput.",flag,"_rld.DP.RData"))
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
    print(paste0("Save file: R-02-dataInput.",flag,"_rpkm.DP.RData"))
    save(multiExpr,nSets, shortLabels, file = paste0("R-02-dataInput.",flag,"_rpkm.DP.RData"))
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

dev.off()

---book
############################
## Choose soft thresholds ##
############################

rdatafiles<-grep("R-02-dataInput",list.files(),value=TRUE)
rdatafiles
# "R-02-dataInput.hylite_rld.RData"   "R-02-dataInput.hylite_rpkm.RData"
# "R-02-dataInput.polycat_rld.RData"  "R-02-dataInput.polycat_rpkm.RData"
# "R-02-dataInput.rsem_rld.RData"     "R-02-dataInput.rsem_rpkm.RData"

for(file in rdatafiles)
{
    print(file)
    load(file)
    tag<-gsub(".*Input[.]|[.].*","",file)
    
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
    pdf(paste0("s4.wgcna.choosePower_", tag,".pdf") )
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

save( list=print(grep("powerTables.",ls(),value=TRUE)), file = "R-03-chooseSoftThreshold.Rdata")

# Inspect "s2.ChooseSoftThresholdPower_???.pdf", 0.8 hard to reach, use 0.7 instead
for(pt in pts)
{
    print(pt)
    power<-unlist(lapply(get(pt), function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.7 & x$slope<0][1]} ) )
    print(power)
}
# "powerTables.hylite_rld"   28 24
# "powerTables.hylite_rpkm"  28 24
# "powerTables.polycat_rld"  22 22 22 22
# "powerTables.polycat_rpkm" 22 22 22 22
# "powerTables.rsem_rld"     22 22 22
# "powerTables.rsem_rpkm"    20 20 20

# previous results
### polycat_rld: 24 28 24 24  =>28
### polycat_rpkm: 24 24 20 24 24 20  =>24
### rsem_rld: 30 28 26  => 30
### rsem_rpkm: 24,22,28 => 28

# I feel Power=24 looks be good for all
for(pt in pts)
{
    print(pt)
    k<-unlist(lapply(get(pt), function(x){ x<-x$data; x[x$Power==24,5:7]} ) )
    print(k)
}


##########################
## Network Construction ##
##########################
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

Powers <- rep(24,6)
names(Powers)<-c("polycat_rld","polycat_rpkm", "rsem_rld", "rsem_rpkm","hylite_rld","hylite_rpkm" )

rdatafiles<-grep("R-02-dataInput",list.files(),value=TRUE)
rdatafiles

ptm <- proc.time()

for(file in rdatafiles)
{
    print(file)
    ll<-load(file)
    compr<-gsub("R-02-dataInput.|.RData","",file)
    powerEach = Powers[compr]
    print(paste0(file,", construct network with b = ",powerEach))
    # print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    
    
    ###### calculate individual TOMs
    print("###### Calculate individual TOMs:")
    iTOMs = blockwiseIndividualTOMs(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    power = powerEach, networkType = "signed", TOMType = "signed",
    # Save individual TOMs?
    saveTOMs = TRUE,
    individualTOMFileNames = paste0(compr,".iTOM-%N-block.%b.RData")  )
    
    ###### calculate consensus modules
    print("###### Construct consensus networks:")
    cnet = blockwiseConsensusModules(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    power = powerEach, networkType = "signed", TOMType = "signed",
    # load previous TOMs
    individualTOMInfo = iTOMs,
    # Saving the consensus TOM
    saveConsensusTOMs = TRUE,
    consensusTOMFileNames = paste0(compr,".cTOM-Block%b.RData"),
    # Basic tree cut options
    deepSplit = 2,  #default, known to reasonable
    minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
    pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
    # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
    mergeCutHeight = 0.25,
    # others
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3)
    
    ###### calculate individual modules
    print("###### Construct individual networks:")
    # stupid blockwiseModules only load TOM rdata file with "TOM", not "tomDS"
    tomFiles<-grep(paste0(compr,".iTOM"),list.files(), value=TRUE)
    for(fl in tomFiles)
    {
        load(fl)
        TOM<-tomDS
        save(TOM, file=fl)
    }
    rm(TOM)
    collectGarbage()
    for(i in 1:nSets)
    {
        inet = blockwiseModules(
        # Input data
        multiExpr[[i]]$data,
        # Data checking options
        checkMissingData = TRUE,
        # Options for splitting data into blocks
        blocks =  iTOMs$blocks,
        #randomSeed = 12345,
        #maxBlockSize =  500,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
        # Network construction arguments: correlation options, use bicor instead of default pearson
        corType = "pearson",
        # Adjacency and topology overlap function options
        power = powerEach, networkType = "signed", TOMType = "signed",
        # load previous TOMs
        loadTOM = TRUE,
        saveTOMFileBase = paste0(compr,".iTOM-",i),
        # Basic tree cut options
        deepSplit = 2,  #default, known to reasonable
        minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
        pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
        # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
        mergeCutHeight = 0.25,
        # others
        reassignThreshold = 0,
        numericLabels = TRUE,
        verbose = 3)
        
        assign(paste0("inet",i), inet)
        
    }
    
    save(list=c( "iTOMs", "cnet", grep("inet.",ls(), value=TRUE)), file = paste0("R-04-buildNetwork.", compr,".RData"))
    collectGarbage()
    
}

proc.time() - ptm

--bookmark
# With maxBlockSize = 20000, it took roughly a week



############### Step 5.  General network topology analysis  ###############
## nohup R CMD BATCH s5.networkTopology.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

plotRefModulePres<-function(mp, file= "")
{
    pdf(file)
    for(set in 2:nSets ){
        # specify the reference and the test networks
        ref=1; test = set
        Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
        Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
        # Look at the observed preservation statistics
        # Obs.PreservationStats
        # Z statistics from the permutation test analysis
        # Z.PreservationStats
        
        # Let us now visualize the data.
        modIDs = rownames(Obs.PreservationStats)
        modColors=labels2colors(order(as.numeric(modIDs) )-1 )
        moduleSize = Obs.PreservationStats$moduleSize
        # we will omit the grey module (background genes)
        # and the gold module (random sample of genes)
        selectModules = !(modColors %in% c("grey", "gold"))
        # Text labels for points
        point.label = modIDs[selectModules]
        # Composite preservation statistics
        medianRank=Obs.PreservationStats$medianRank.pres
        Zsummary=Z.PreservationStats$Zsummary.pres
        
        par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
        # plot medianRank versus module size
        plot(moduleSize[selectModules],medianRank[selectModules],col=1, bg=modColors[selectModules],
        pch = 21,main=paste("medianRank -",shortLabels[ref], "vs",shortLabels[test]),
        cex = 2, ylab ="medianRank",xlab="Module size", log="x")
        labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1,offs=0.03)
        
        # plot Zsummary versus module size
        plot(moduleSize[selectModules],Zsummary[selectModules], col = 1, bg=modColors[selectModules],pch = 21,
        main=paste("Zsummary -",shortLabels[ref], "vs",shortLabels[test]),
        cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
        labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
        # Add threshold lines for Zsummary
        abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
    }
    dev.off()
}


ptm <- proc.time()

Powers <- c(28,24,30,28)
names(Powers)<-c("polycat_rld","polycat_rpkm", "rsem_rld", "rsem_rpkm")


rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file) # "multiExpr"   "nSets"       "setLabels"   "shortLabels"
    compr<-gsub("R-02-dataInput.|.RData","",file)
    powerEach = Powers[compr]
    print(paste0(file,", construct network with b = ",powerEach))
    print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    load(paste0("R-04-buildNetwork.", compr,".RData")) ->Names
    print(Names)  # "iTOMs" "cnet"  "inet1" "inet2" "inet3"
    
    ###  1.Comparing basic topological networks parameters
    ######################################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        gg<-net$goodGenes    #get the good genes descision
        
        ###  Get basic topological networks parameters
        print(paste("Start to extract network parameter for:",genome))
        # network concept
        adj <- adjacency(subDat, power = powerEach,type = "signed")   # calculate adjacency table
        # fundmentalNetworkConcept taking forever
        # net1C <- fundamentalNetworkConcepts(adj, GS = NULL)           # computes fundamental network concepts
        Size = dim(adj)[1]
        Connectivity = apply(adj, 2, sum)-1
        Density = sum(Connectivity)/(Size * (Size - 1))
        Centralization = Size * (max(Connectivity) - mean(Connectivity))/((Size -
        1) * (Size - 2))
        Heterogeneity = sqrt(Size * sum(Connectivity^2)/sum(Connectivity)^2 -
        1)
        # clustering coef takes forever to computate, skip
        # ClusterCoef = .ClusterCoef.fun(adj)
        fMAR = function(v) sum(v^2)/sum(v)
        MAR = apply(adj, 1, fMAR)
        ScaledConnectivity = Connectivity/max(Connectivity, na.rm = T)
        output = list(Connectivity = Connectivity, ScaledConnectivity = ScaledConnectivity,  MAR = MAR, Density = Density, Centralization = Centralization, Heterogeneity = Heterogeneity)
        param <- unlist(lapply(output,mean))
        
        ### Get module and dupllication related topology
        # how many modules detected??
        param["nModules"]<-length(unique(net$colors))
        # interal and crossing connections
        dt<-grep("d$",rownames(adj))
        size.dt<-length(dt)
        at<-grep("a$",rownames(adj))
        size.at<-length(at)
        param["at.density"] <- (sum(adj[at,at])-size.at)/(size.at*(size.at-1))
        param["dt.density"] <- (sum(adj[dt,dt])-size.dt)/(size.dt*(size.dt-1))
        param["ad.density"] <- mean(adj[at,dt])
        ### write results table
        if(!exists("Rtable"))
        {
            Rtable <- data.frame(param)
            names(Rtable)<-genome
        }else
        {Rtable[,genome] <- param }
        
        
        assign(paste0(genome,"Adj"),adj)
        assign(paste0(genome,"Concepts"),output)
        
    }
    #  Save basic topological networks parameters
    rownames(Rtable)[1:3]<- paste0("mean",rownames(Rtable)[1:3])
    write.table(Rtable, file = paste0("s5.parameters.",compr,".txt"),sep="\t")
    save( list=c( grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
    
    
    ###  2.Correlating node-specific network properties
    ###################################################
    pwset<-combn(nSets,2)
    
    # pdf(paste0("s5.correlatingParameters.",compr,".pdf") )
    # par(mfrow=c(2,2))
    # for( j in 1:3) # "Connectivity"       "ScaledConnectivity" "MAR"
    # {
    #    for(i in 1:ncol(pwset) )
    #    {   xnet <- get(paste0(shortLabels[ pwset[1,i]  ],"Concepts"))
    #        ynet <- get(paste0(shortLabels[ pwset[2,i]  ],"Concepts"))
    #        pp <- names(xnet)[j]
    #        corr <- cor.test(xnet[[j]], ynet[[j]])
    #        plot (xnet[[j]], ynet[[j]], main = paste(pp, ": cor=", round(corr$estimate,2), ", p=", corr$p.value, sep=""), xlab= shortLabels[pwset[1,i]], ylab = shortLabels[pwset[2,i]], type ="p", pch=16, col = rgb(0, 0, 0, 0.2) )
    #    }
    # }
    # dev.off()
    
    df_Connectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$Connectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$Connectivity)
    df_ScaledConnectivity <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$ScaledConnectivity, get(paste0(shortLabels[ 2  ],"Concepts"))$ScaledConnectivity)
    df_MAR <- cbind( get(paste0(shortLabels[ 1  ],"Concepts"))$MAR, get(paste0(shortLabels[ 2  ],"Concepts"))$MAR)
    df_expr <- cbind(  as.numeric(as.matrix(multiExpr[[1]]$data)),  as.numeric(as.matrix(multiExpr[[2]]$data))  )
    for(i in 3:nSets)
    {
        df_Connectivity <- cbind( df_Connectivity, get(paste0(shortLabels[ i  ],"Concepts"))$Connectivity)
        df_ScaledConnectivity <- cbind( df_ScaledConnectivity , get(paste0(shortLabels[ i  ],"Concepts"))$ScaledConnectivity)
        df_MAR <- cbind( df_MAR, get(paste0(shortLabels[ i  ],"Concepts"))$MAR)
        df_expr <- cbind(  df_expr,  as.numeric(as.matrix(multiExpr[[i]]$data))  )
    }
    colnames(df_Connectivity) <-shortLabels
    colnames(df_ScaledConnectivity) <-shortLabels
    colnames(df_MAR) <-shortLabels
    colnames(df_expr) <-shortLabels
    panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    
    pdf(paste0("s5.correlatingParameters.",compr,".k.pdf") )
    pairs(df_Connectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".sk.pdf") )
    pairs(df_ScaledConnectivity, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Scaled Connectivity Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".mar.pdf") )
    pairs(df_MAR, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene MAR Scatterplot Matrix")
    dev.off()
    pdf(paste0("s5.correlatingParameters.",compr,".expr.pdf") )
    pairs(df_expr, lower.panel=panel.smooth, upper.panel=panel.cor, pch=20, col = rgb(0, 0, 0, 0.05), main="Gene Expression Scatterplot Matrix")
    dev.off()
    
    
    
    ###  3.WGCNA provided methods for comparison
    ############################################
    # get all anova p and put together
    dpa=as.factor(c(10,10,10,20,20,30,30,30,40,40,40))
    anovaP<-list()
    for(set in 1:nSets)
    {
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        MEs<-net$MEs
        pval<-apply(MEs,2,function(x){round(anova(aov(x~dpa) )$"Pr(>F)"[1],4)})
        pval<-as.data.frame(pval)
        pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
        pval$numeric<-as.numeric(substring(rownames(pval),3) )
        pval<-pval[order(pval$numeric),]
        pval$symbol[1]<-" "  # ME0 always meaningless
        anovaP[[set]]<-pval
    }
    names(anovaP)<-shortLabels
    
    # Module correspondence test
    pdf(paste0("s5.moduleCorrepondence.",compr,".pdf"),width=10,height=7)
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(8, 10.4, 2.7, 1)+0.3);
    # loop pairwise comparison
    for(i in 1:ncol(pwset))
    {
        print(paste0("Plot for pairwise set ",i))
        colnet<-get(paste0("inet",pwset[1,i]) )
        rownet<-get(paste0("inet",pwset[2,i]) )
        coln<-ncol(colnet$MEs )  # number of MEs
        rown<-ncol(rownet$MEs )
        # color list of MEs in the color of decreasing numbers of memebers
        colModule  <- labels2colors(as.numeric(names(table(colnet$colors)) ))
        rowModule  <- labels2colors(as.numeric(names(table(rownet$colors)) ))
        # colors for each gene
        colColors  <- labels2colors(colnet$colors )
        rowColors  <- labels2colors(rownet$colors )   # colors for each gene
        # anova significance sybol
        colP  <- anovaP[[pwset[1,i]]]$symbol
        rowP  <- anovaP[[pwset[2,i]]]$symbol
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = rown, ncol = coln);
        CountTbl = matrix(0, nrow = rown, ncol = coln);
        # Execute all pairwaise comparisons
        for (rmod in 1:rown)
        for (cmod in 1:coln)
        {
            rMembers = (rowColors == rowModule[rmod] );
            cMembers = (colColors == colModule[cmod] );
            pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
            CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
        }
        
        # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        rModTotals = apply(CountTbl, 1, sum)
        cModTotals = apply(CountTbl, 2, sum)
        # Use function labeledHeatmap to produce the color-coded table with all the trimmings
        labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
        xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".",colModule, "_", cModTotals, colP),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".",rowModule, "_", rModTotals, rowP),
        textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
        
        # plotting only singficant modules
        labeledHeatmap( Matrix = pTable[rowP=="*",colP=="*"], colorLabels = TRUE,
        xLabels = paste(" ", colModule[colP=="*"]), yLabels = paste(" ", rowModule[rowP=="*"]),
        xSymbols = paste0(shortLabels[pwset[1,i]], ".", colModule[colP=="*"], ": ", cModTotals[colP=="*"], colP[colP=="*"]),
        ySymbols = paste0(shortLabels[pwset[2,i]], ".", rowModule[rowP=="*"], ": ", rModTotals[rowP=="*"], rowP[rowP=="*"]),
        textMatrix = CountTbl[rowP=="*",colP=="*"], colors = blueWhiteRed(100)[50:100],
        main = "Correspondence of modules, significant only",
        cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
    }
    dev.off()
    
    # Reference network based preservation test
    # set A2D5 as reference
    multiColor<-list()
    for(set in 1:nSets)
    {
        multiColor[[set]] <- get(paste0("inet",set) )$colors
    }
    names(multiColor) <-shortLabels
    # The number of permutations drives the computation time of the module preservation function. For a publication use 200 permutations.
    # But for brevity and testing, start with a small number
    nPermutations1=200
    # Set it to a low number (e.g. 3) if only the medianRank statistic and other observed statistics are needed.
    # Permutations are only needed for calculating Zsummary and other permutation test statistics.
    # set the random seed of the permutation test analysis
    set.seed(1)
    system.time({
        mp = modulePreservation(multiExpr, multiColor, networkType="signed", referenceNetworks = 1, testNetworks=NULL, nPermutations = nPermutations1,
        permutedStatisticsFile = paste0("s5.permutedStats-intrModule.", compr,".RData"),
        randomSeed = 1, quickCor = 0, verbose = 3)
    })
    # for 3 permutations, 4539.133s as roughly 1.26 hours
    # so it will probably takes 3.5 days to run 200 permutations
    save( list=c("anovaP", "mp", grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
    # plot preservation test results
    plotRefModulePres(mp, file= paste0("s5.refPreservation.", compr,".pdf"))
    
    
    ###  4.Homoeolog focused analysis
    ############################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        
        probes <- colnames(subDat )
        genes <- gsub(".$","",probes)
        subgenome <- gsub(".*00","",probes)
        homoeo <- c(length(intersect( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        colors <-net$colors
        inModuleNo<-aggregate(net$colors, by=list(genes), function(x)length(unique(x)))
        inSameModule <- length(which(inModuleNo$x==1))
        dominance<- as.data.frame(unclass(xtabs(~net$colors+subgenome) ) )
        dominance$chisq.p<-apply(dominance,1,function(x)chisq.test(x)$p.value)
        dominance$dominant <- apply(dominance,1,function(x)ifelse(x[3]>0.05,"-",ifelse(x[1]>x[2],"a","d")) )
        dominant <- factor(dominance$dominant, levels=c("a","b","-"))
        homoeo<-c(homoeo, tabulate(dominant) )
        names(homoeo)[4:6]<-levels(dominant)
        if(!exists("Rtableh"))
        { Rtableh<-data.frame(homoeo)
            names(Rtableh)<-genome
        }else {Rtableh[,genome]<-homoeo}
        
        assign(paste0(genome,"dominance"),dominance)
    }
    write.table(Rtableh, file = paste0("s5.homoeolog.",compr,".txt"),sep="\t")
    
    # save everything, AGAIN
    save( list=c("anovaP", "mp", grep(".dominance",ls(), value=TRUE), grep("Concepts",ls(), value=TRUE)), file = paste0("R-05-networkTopology.", compr,".RData"))
}

proc.time() - ptm



############### Step END.  Extra analysis  ###############
## nohup R CMD BATCH s.networkTopology.R &
################

# Plot grapgh parameters with PCA
library(ggplot2)
pdf("s6.graphParametersPCA.pdf")
for(i in c("polycat_rld","polycat_rpkm","rsem_rld","rsem_rpkm"))
{
    f<-read.table(file = paste0("s5.parameters.",i,".txt"), header=TRUE, sep="\t" )
    pca = prcomp(t(f), scale.=T)
    
    # scale = 0 ensures that arrows are scaled to represent the loadings.
    # focus on the extreme ends (top, bottom, left, right), e.g. first principal component corresponds to a measure of left and right arrows.
    # For exact measure of a variable in a component, pca$rotation
    biplot(pca, scale = 0 , main=i)
    
    dat = as.data.frame(pca$x)
    dat$group= rownames(dat)
    portion.var <- as.numeric( summary(pca)$importance[2,] )
    print(
    ggplot(aes(PC1, PC2, color=group),data=dat) +
    geom_point() + ggtitle(i) +
    xlab(paste0("PC1: ",portion.var[1]*100,"% variance")) +
    ylab(paste0("PC2: ",portion.var[2]*100,"% variance"))
    )
}
# so empir corrected profiles are more similiar to uncorrected rpkm than theoretical eflen corrected profiles??
dev.off()



# Comparing parameters directly, get Z scores
rdatafiles<-c("R-05-networkTopology.rsem_rpkm.RData", "R-05-networkTopology.rsem_rld.RData", "R-05-networkTopology.polycat_rpkm.RData", "R-05-networkTopology.polycat_rld.RData")
Zresult <- rbind(c("Method","Obs","parameter","Zscore"))
for(f in rdatafiles)
{
    print(f)
    print("-----------------")
    
    rm(list=grep("Concepts",ls(),value=TRUE))
    ll<-load(f)
    
    ref= "A2D5Concepts"
    testG<-grep("Concepts",ls(), value=TRUE)
    testG<-testG[testG!=ref]
    nameG<-names(get(ref))
    for(test in testG)
    {
        result <- c(gsub(".RData|R-05-networkTopology.","",f), test)
        print(paste0(test," is compared to reference ", ref))
        for(i in 1:3)
        {
            x<- get(test)[[i]] - get(ref)[[i]]
            print(paste0("     ", nameG[i], ": mean - ", mean(x), "; sd - ", sd(x)))
            Zresult <- rbind(Zresult, c(result, nameG[i], x/sd(x)))
            
        }
    }
}
Z<-data.frame(Zresult[-1,])
names(Z) <- Zresult[1,]
Z$Obs<-gsub("Concepts","", Z$Obs)
Z$Zscore <- as.numeric(Z$Zscore)
Z$p.value <- pnorm(Z$Zscore)

write.table(Z, file = "s6.parameterZscore.txt",sep="\t", row.names=FALSE)


# "R-05-networkTopology.rsem_rpkm.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - 0.781465672186864; sd - 98.6522941737412"
# "     ScaledConnectivity: mean - 0.000947613568534769; sd - 0.0217144740877903"
# "     MAR: mean - 0.00450310922680791; sd - 0.0363252496455325"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -36.6065926997929; sd - 269.453978935194"
# "     ScaledConnectivity: mean - -0.00689189844026479; sd - 0.0593641940042516"
# "     MAR: mean - -0.00665182524599488; sd - 0.0741680503434127"
[1] "R-05-networkTopology.rsem_rpkm.RData"
[1] "-----------------"
[1] "A2D5.techConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - 0.785287565474176; sd - 57.5508927554109"
[1] "     ScaledConnectivity: mean - 0.00123342597671302; sd - 0.0199828561365023"
[1] "     MAR: mean - 0.00495548095961635; sd - 0.0432059016086636"
[1] "ADsConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -13.8622253567885; sd - 164.335781159567"
[1] "     ScaledConnectivity: mean - -0.00747326262845879; sd - 0.0568774292673812"
[1] "     MAR: mean - -0.00551956296440461; sd - 0.0883282894623805"


# "R-05-networkTopology.rsem_rld.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - 6.03432521331369; sd - 47.9331601209624"
# "     ScaledConnectivity: mean - 0.00385747478587089; sd - 0.0311011472237645"
# "     MAR: mean - 0.0108177203452433; sd - 0.0665238512411918"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -15.8272119289457; sd - 85.875353209918"
# "     ScaledConnectivity: mean - -0.00152867570149998; sd - 0.0564456719757156"
# "     MAR: mean - 0.000917384941257472; sd - 0.0845454737261381"


# "R-05-networkTopology.polycat_rpkm.RData"
# "-----------------"
# "A2D5.tech.theoConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -21.9512000285965; sd - 240.137834569537"
# "     ScaledConnectivity: mean - -0.00280878884515633; sd - 0.0466000052518468"
# "     MAR: mean - 0.00740427188747321; sd - 0.0679811075400122"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -21.9512000285965; sd - 240.137834569537"
# "     ScaledConnectivity: mean - -0.00280878884515633; sd - 0.0466000052518468"
# "     MAR: mean - 0.00740427188747321; sd - 0.0679811075400122"
# "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -24.3006652983133; sd - 202.22935432233"
# "     ScaledConnectivity: mean - -0.00126881936965774; sd - 0.0391716183150736"
# "     MAR: mean - 0.00560018426316752; sd - 0.060594947256564"
# "ADs.theoConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -23.3773020911422; sd - 247.50336056002"
# "     ScaledConnectivity: mean - -0.00385318188826418; sd - 0.0479840016318954"
# "     MAR: mean - 0.00369977480503477; sd - 0.0628610552950785"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -23.3773020911422; sd - 247.50336056002"
# "     ScaledConnectivity: mean - -0.00385318188826418; sd - 0.0479840016318954"
# "     MAR: mean - 0.00369977480503477; sd - 0.0628610552950785"
[1] "R-05-networkTopology.polycat_rpkm.RData"
[1] "-----------------"
[1] "A2D5.techConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -22.0765434162768; sd - 163.648523219436"
[1] "     ScaledConnectivity: mean - -0.00357149736531711; sd - 0.044308874996079"
[1] "     MAR: mean - 0.00925051099948535; sd - 0.0796255066248835"
[1] "A2D5.tech.theoConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -23.4762894834495; sd - 163.815572860369"
[1] "     ScaledConnectivity: mean - -0.00353351048896584; sd - 0.0443603130168054"
[1] "     MAR: mean - 0.00908386695385239; sd - 0.0799041069808706"
[1] "ADsConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -22.7441966541389; sd - 166.634032186635"
[1] "     ScaledConnectivity: mean - -0.00430959600795669; sd - 0.045165704424851"
[1] "     MAR: mean - 0.00491785060911022; sd - 0.073303344500602"
[1] "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -21.0670290871299; sd - 137.304169417926"
[1] "     ScaledConnectivity: mean - -0.00211100698492692; sd - 0.0370154162797463"
[1] "     MAR: mean - 0.00673809793888078; sd - 0.0712921811660626"
[1] "ADs.theoConcepts is compared to reference A2D5Concepts"
[1] "     Connectivity: mean - -24.1557787577917; sd - 166.67064086685"
[1] "     ScaledConnectivity: mean - -0.00426954496386292; sd - 0.0451968808000639"
[1] "     MAR: mean - 0.00474241482081538; sd - 0.0735389334852459"

# "R-05-networkTopology.polycat_rld.RData"
# "-----------------"
# "A2D5.techConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -7.87451841852331; sd - 82.7567479726778"
# "     ScaledConnectivity: mean - -0.00128175364156708; sd - 0.0497757523553588"
# "     MAR: mean - 0.0126867937645758; sd - 0.0881149811226852"
# "ADs.ncorrectConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -4.52782204245409; sd - 73.9216718803422"
# "     ScaledConnectivity: mean - 0.00415615488438032; sd - 0.0450874147901161"
# "     MAR: mean - 0.00904897219457412; sd - 0.0806527881971959"
# "ADsConcepts is compared to reference A2D5Concepts"
# "     Connectivity: mean - -11.1175246601803; sd - 81.7896357500773"
# "     ScaledConnectivity: mean - -0.00145226483274745; sd - 0.0492663618513508"
# "     MAR: mean - 0.00810894907197572; sd - 0.0825466056186121"



rdatafiles<-c("R-05-networkTopology.rsem_rpkm.RData", "R-05-networkTopology.rsem_rld.RData", "R-05-networkTopology.polycat_rpkm.RData", "R-05-networkTopology.polycat_rld.RData")
# view A2D5 vs A2D5.tech median rank
# mp$preservation$observed[[1]][[2]] ->mr
# mr[order(mr$medianRank.pres,decreasing=TRUE),1:2]

pdf("s6.sumZsummary.pdf")
Zlist<- c("Zsummary.pres", "Zdensity.pres", "Zconnectivity.pres")
upper<-c(100,160,60)
for(i in 1:3)
{
    z<-Zlist[i]
    u<-upper[i]
    load("R-05-networkTopology.polycat_rpkm.RData")
    plot(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="blue", ylim=c(0,u), xlab="Module Size", ylab="Zsummary", main=z)
    points(mp$preservation$Z[[1]][[4]]$moduleSize, mp$preservation$Z[[1]][[4]][,z], col="blue", pch = 2) #ncorrect
    # points(mp$preservation$Z[[1]][[6]]$moduleSize, mp$preservation$Z[[1]][[6]]$Zsummary.pres, col="blue", pch=0)
    
    load("R-05-networkTopology.polycat_rld.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="brown")
    points(mp$preservation$Z[[1]][[4]]$moduleSize, mp$preservation$Z[[1]][[4]][,z], col="brown", pch = 2) #ncorrect
    
    load("R-05-networkTopology.rsem_rpkm.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="green")
    load("R-05-networkTopology.rsem_rld.RData")
    points(mp$preservation$Z[[1]][[3]]$moduleSize, mp$preservation$Z[[1]][[3]][,z], col="purple")
    
    abline(h=10,col="red")
    legend(100,u,legend=c("polycat_rpkm","polycat_rpkm: Npartition","polycat_rld","polycat_rld: Npartition","rsem_rpkm","rsem_rld"), col=c("blue","blue","brown","brown","green","purple"), pch=c(1,2,1,2,1,1))
    
}

#### Separability statistics performs poorly

load("R-05-networkTopology.polycat_rpkm.RData")
plot(mp$testSeparability$Z[[1]][[3]]$moduleSize,mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="blue", ylim=c(0,1500), xlab="Module Size", ylab="Zsummary", main="Z.separability.pres")
points(mp$testSeparability$Z[[1]][[4]]$moduleSize, mp$testSeparability$Z[[1]][[4]][,"Z.separability.pres"], col="blue", pch = 2) #ncorrect
# points(mp$preservation$Z[[1]][[6]]$moduleSize, mp$preservation$Z[[1]][[6]]$Zsummary.pres, col="blue", pch=0)

load("R-05-networkTopology.polycat_rld.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="brown")
points(mp$testSeparability$Z[[1]][[4]]$moduleSize, mp$testSeparability$Z[[1]][[4]][,"Z.separability.pres"], col="brown", pch = 2) #ncorrect

load("R-05-networkTopology.rsem_rpkm.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="green")
load("R-05-networkTopology.rsem_rld.RData")
points(mp$testSeparability$Z[[1]][[3]]$moduleSize, mp$testSeparability$Z[[1]][[3]][,"Z.separability.pres"], col="purple")

legend(100,1500,legend=c("polycat_rpkm","polycat_rpkm: Npartition","polycat_rld","polycat_rld: Npartition","rsem_rpkm","rsem_rld"), col=c("blue","blue","brown","brown","green","purple"), pch=c(1,2,1,2,1,1))
abline(h=10,col="red")
dev.off()


# how homoeolog pairs in the same module

rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.rsem_rld.RData", "R-02-dataInput.polycat_rpkm.RData", "R-02-dataInput.polycat_rld.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file) # "multiExpr"   "nSets"       "setLabels"   "shortLabels"
    compr<-gsub("R-02-dataInput.|.RData","",file)
    print(setLabels)
    print(checkSets(multiExpr)$nGenes)
    load(paste0("R-04-buildNetwork.", compr,".RData")) ->Names
    print(Names)  # "iTOMs" "cnet"  "inet1" "inet2" "inet3"
    ###  4.Homoeolog focused analysis
    ############################################
    for (set in 1:nSets )
    {
        # Extract total read counts for each genome
        subDat    <-  multiExpr[[set]]$data
        genome <-  shortLabels[set]
        net<-get(paste0("inet",set) )
        
        probes <- colnames(subDat )
        genes <- gsub(".$","",probes)
        subgenome <- gsub(".*00","",probes)
        inPairs<-intersect( genes[subgenome =="a"], genes[subgenome=="d"] )
        homoeo <- c(length( inPairs ), length(setdiff( genes[subgenome =="a"], genes[subgenome=="d"] ) ), length(setdiff( genes[subgenome =="d"], genes[subgenome=="a"] ) ))
        names(homoeo)<-c( "ADinNet", "AinNet", "DinNet" )
        colors <-net$colors
        inModuleNo<-aggregate(net$colors, by=list(genes), function(x)length(unique(x)))
        rownames(inModuleNo)<-inModuleNo$Group.1
        inModuleNo<-inModuleNo[inPairs,]
        inSameModule <- length(which(inModuleNo$x==1))
        homoeo["ADinSameModule"]<-inSameModule
        if(set==1)
        {
            true<-inModuleNo;
            homoeo["ADinSM:sensitivity"]<-"-";
            homoeo["ADinSM:specificity"]<-"-"
        }else
        {
            truePair <- true$Group.1[true$x==1]
            obsPair  <- inModuleNo$Group.1[inModuleNo$x==1]
            # sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
            homoeo["ADinSM:sensitivity"]<- length(intersect(truePair,obsPair) )/length(truePair)
            
            truePairN <- true$Group.1[true$x>1]
            obsPairN  <- inModuleNo$Group.1[inModuleNo$x>1]
            # specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)
            homoeo["ADinSM:specificity"]<- length(intersect(truePairN,obsPairN) )/length(truePairN)
        }
        dominance<- as.data.frame(unclass(xtabs(~net$colors+subgenome) ) )
        dominance$chisq.p<-apply(dominance,1,function(x)chisq.test(x)$p.value)
        dominance$dominant <- apply(dominance,1,function(x)ifelse(x[3]>0.05,"-",ifelse(x[1]>x[2],"a","d")) )
        homoeo<-c(homoeo, unclass(table(dominance$dominant) ) )
        if(!exists("Rtableh"))
        { Rtableh<-data.frame(homoeo)
            names(Rtableh)<-genome
        }else {Rtableh[,genome]<-homoeo}
        
        assign(paste0(genome,"dominance"),dominance)
    }
    write.table(Rtableh, file = paste0("s5.homoeologN.",compr,".txt"),sep="\t")
    rm(Rtableh)
}
# sensitivity, true-positive rate (TPR), both.sig / A2D5.sig
# specificity),true-negative rate (TNR)  both.not.sig / (37223-A2D5.not)



