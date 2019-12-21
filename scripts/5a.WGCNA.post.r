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

sessionInfo()

rdatafiles<-grep("R-05-networkTopology",list.files(),value=TRUE)
rdatafiles
rm(res)
for(file in rdatafiles)
{
    # collect input dataset
    print(file)
    load(file) #refNet,mp
    pipeline= gsub("R-05-networkTopology.|_.*","", file)
    transformation = gsub(".*_|.B.*","", file)
    netType=paste0("weighted.WGCNA.", gsub(".*B|.RData","", file))
    Zs =mp$preservation$Z[[1]][[2]]$Zsummary.pres # as long as above 10
    df=data.frame(pipeline, transformation, netType, moduleN = length(Zs), Zs)
    if(!exists("res")){res<-df}else{res<-rbind(res,df)}
   
}

library(ggplot2);
pdf("s5.Zsummary.pdf")
ggplot(data=res, aes(x=netType, y=Zs, fill=pipeline )) +  geom_boxplot(lwd=0.2,position=position_dodge(0.8)) + theme_bw() +theme(panel.border=element_blank(),axis.line.x=element_line(),axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Preservation of expected modules, Zsummary") + facet_grid(.~transformation)
dev.off()
