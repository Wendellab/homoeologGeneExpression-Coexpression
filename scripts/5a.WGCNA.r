#FUN
plotRefModulePres<-function(mp,refIdx=1,testIdx=2, file= "")
{
    pdf(file)
    ref=refIdx
    for(test in testIdx ){
        # specify the reference and the test networks
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
##########################
## Network Construction ##
##########################
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=10)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

rdatafiles<-grep("R-05-dataInput",list.files(),value=TRUE)
rdatafiles
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

Powers <- c(1,12,24)

for(file in rdatafiles)
{
    # collect input dataset
    print(file)
    load(file) # "multiExpr"   "nSets"       "shortLabels"
    ori=multiExpr
    multiExpr = vector(mode = "list", length = 2)
    multiExpr[[1]] = list(data = ori[[which(shortLabels=="A2D5")]]$data[,use])
    multiExpr[[2]] = list(data = ori[[which(shortLabels=="ADs")]]$data[,use])
    print(checkSets(multiExpr)$nGenes)
    tag=gsub(".*Input[.]|[.].*","",file)
    
    # loop powers
    for(power in Powers){
        # construct A2D5 network
        refNet =  blockwiseModules(
            # Input data
            multiExpr[[1]]$data,
            # Data checking options
            checkMissingData = TRUE,
            # Options for splitting data into blocks
            maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
            #randomSeed = 12345,
            # Network construction arguments: correlation options, use bicor instead of default pearson
            corType = "pearson",
            # Adjacency and topology overlap function options
            power = power, networkType = "signed", TOMType = "signed",
            # load previous TOMs
            saveTOMs = FALSE,
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
        # preservation tests of ADs agains A2D5
        # The number of permutations drives the computation time of the module preservation function. For a publication use 200 permutations.
        # But for brevity and testing, start with a small number
        nP=200
        # Set it to a low number (e.g. 3) if only the medianRank statistic and other observed statistics are needed.
        # Permutations are only needed for calculating Zsummary and other permutation test statistics.
        # set the random seed of the permutation test analysis
        set.seed(1)
        multiColor=list(refNet$colors, refNet$colors)
        mp = modulePreservation(multiExpr, multiColor, networkType="signed", referenceNetworks = 1, testNetworks=2, nPermutations = nP, savePermutedStatistics = FALSE,
            # permutedStatisticsFile = paste0("s5.netPreservation.",tag,".B",power,".RData"),
            randomSeed = 1, quickCor = 0, verbose = 3)
        # for 3 permutations, roughly 20min
        save( refNet,mp, file = paste0("R-05-networkTopology.", tag,".B",power,".RData"))
        # plot preservation test results
        plotRefModulePres(mp, file= paste0("s5.refPreservation.", tag,".B",power,".pdf"))
        Zs =mp$preservation$Z[[1]][[2]]$Zsummary.pres
        # as long as above 10
        print(table(Zs>10))
    }
}
