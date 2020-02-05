############### Step 1a. Prepare polycat network datasets  ###############
## unix command: R CMD BATCH s1.R
################
# 2020 revision
# Read classification - old polycat [v1.3], new polycat [v2?], old gsnap-polycat (original files)
# Counting - featureCounts, htseq, htseq primary alignments only

#FUN
readHTseqCounts<-function(dir){
    fileL<-grep("txt",list.files(dir),value=TRUE)
     # 3 species x 11 samples x 5(total, A, D, N, AD)  = 165
    for(file in fileL)
    {
        x <- read.table(paste0(dir,file), sep="\t")
        names(x) <- c("Geneid", file)
        if(!exists("allcount")) {allcount = x}
        else {allcount <- merge(allcount, x, by="Geneid")}
    }
    row.names(allcount) <- allcount$Geneid
    names(allcount)<-gsub("-",".",gsub(".txt", "", names(allcount) ))
    allcount<-allcount[grep("Gorai.0",rownames(allcount)),]
    return(allcount)
}


#### New GSNAP + New PolyCat + different counting
# gsnapD5.seed.newPolyCat.featureCounts = new gsnap against D5 ref using new polycat and feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.newPolyCat.featureCounts"
seed.NNF <- read.table(seed_file, sep="\t",header=T)
# htseq-count-new = new gsnap against D5 ref using new polycat and htseq (but not specifying --secondary-alignments=ignore)
seed_dir <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/htseq-count-new/"
seed.NNH=readHTseqCounts(seed_dir)
# htseq-count-new-primary-only = new gsnap against D5 ref using new polycat and htseq (specifying --secondary-alignments=ignore) * this should be more similar to featureCounts (edited)
seed_dir <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/htseq-count-new-primary-only/"
seed.NNHp=readHTseqCounts(seed_dir)


# gsnapD5.seed.oldPolyCat.featureCounts = new gsnap against D5 ref using old polycat version and feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.oldPolyCat.featureCounts"
seed.NOF <- read.table(seed_file, sep="\t",header=T)

# gsnapD5.seed.oldGSNAP.featureCounts = our original run; old gsnap against D5 with old polycat, but counted with feature counts
seed_file <- "/work/LAS/jfw-lab/hugj2006/eflen2020/gsnapD5/seed/differentCountTables/gsnapD5.seed.oldGSNAP.featureCounts"
# read files
seed.OOF <- read.table(seed_file, sep="\t",header=T)


##################
source("/work/LAS/jfw-lab/hugj2006/eflen2020/Ranalysis/2.1.FUN.r")
quickEval=function(count, withTotal =TRUE){
    names(count) =gsub(".bam|.sort","",names(count))
    names(count)[grep("[1-9]$",names(count))] = paste0(grep("[1-9]$",names(count),value=T),".T")
    A2.Total  <- count[,grep('A2.*[.]T$',names(count))]
    ADs.Total <- count[,grep('AD.*[.]T$',names(count))]
    D5.Total  <- count[,grep('D5.*[.]T$',names(count))]
    ##################
    A2.At  <- count[,grep('A2.*[.]A$',names(count))]
    ADs.At <- count[,grep('AD.*[.]A$',names(count))]
    D5.At  <- count[,grep('D5.*[.]A$',names(count))]
    ##################
    A2.Dt  <- count[,grep('A2.*[.]D$',names(count))]
    ADs.Dt <- count[,grep('AD.*[.]D$',names(count))]
    D5.Dt  <- count[,grep('D5.*[.]D$',names(count))]
    ##
    resBinA=binC(TP = A2.At, TN = D5.Dt, FP = D5.At ,FN = A2.Dt)
    resBinD=binC(TP = D5.Dt, TN = A2.At, FP = A2.Dt ,FN = D5.At)
    oa=unlist(resBinA$summary[1,])
    od=unlist(resBinD$summary[1,])
    ##
    if(withTotal){
        total <- D5.Total + A2.Total
        totalA <- A2.Total
        totalD <- D5.Total
        assigned <- ADs.At + ADs.Dt
        assignedA <- ADs.At
        assignedD <- ADs.Dt
        efficiency <- assigned/total
        efficiencyA <- assignedA/totalA
        efficiencyD <- assignedD/totalD
        e <- c(sum(assigned)/sum(total), sum(assignedA)/sum(totalA), sum(assignedD)/sum(totalD))
        ##
        different <- abs(ADs.At-A2.Total)+abs(ADs.Dt-D5.Total)
        differentA <- abs(ADs.At-A2.Total)
        differentD <- abs(ADs.Dt-D5.Total)
        discrepancy <- different/total
        discrepancyA <- differentA/totalA
        discrepancyD <- differentD/totalD
        d <- c(sum(different)/sum(total), sum(differentA)/sum(totalA), sum(differentD)/sum(totalD))
        ##
        me= c(e,d,oa[1],od[1],oa[2],od[2], oa[4],od[4],oa[3],oa[5])
    }else{
        me= c(rep("-",6),oa[1],od[1],oa[2],od[2], oa[4],od[4],oa[3],oa[5])
    }
    names(me) = c("Efficiency","Efficiency.At","Efficiency.Dt", "Discrepancy","Discrepancy.At","Discrepancy.Dt","Precision.At","Precision.Dt","Recall.At","Recall.Dt","F1.At","F1.Dt","Accuracy","MCC")
    return(me)
}

TotalLibSize =data.frame(OOF=colSums(seed.OOF[,grep("[1-9].sort.bam",names(seed.OOF))]), NOF=colSums(seed.NOF[,grep("[1-9]$",names(seed.NOF))]), NNF=colSums(seed.NOF[,grep("[1-9]$",names(seed.NNF))]), NNH=colSums(seed.NNH[,grep("T$",names(seed.NNH))]), NNHp=colSums(seed.NNHp[,grep("T$",names(seed.NNHp))]))

eval=data.frame(OOF=quickEval(seed.OOF), NOF=quickEval(seed.NOF), NNF= quickEval(seed.NNF, withTotal=FALSE), NNH=quickEval(seed.NNH),NNHp=quickEval(seed.NNHp))

eval
#                      OOF       NOF       NNH      NNHp               NNF
# Efficiency     0.7575119 0.7534309 0.7529318 0.7529318                 -
# Efficiency.At  0.7512763 0.7437990 0.7440454 0.7440454                 -
# Efficiency.Dt  0.7636301 0.7628783 0.7616178 0.7616178                 -
# Discrepancy    0.2434765 0.2481968 0.2487259 0.2487259                 -
# Discrepancy.At 0.2490796 0.2570309 0.2567937 0.2567937                 -
# Discrepancy.Dt 0.2379789 0.2395320 0.2408400 0.2408400                 -
# Precision.At   0.9855939 0.8822852 0.9327068 0.9327068 0.908715926502718
# Precision.Dt   0.9773479 0.8545057 0.9046178 0.9046178 0.878695323941279
# Recall.At      0.9807249 0.8938723 0.9164451 0.9164451 0.892485425645783
# Recall.Dt      0.9811259 0.8981645 0.9225111 0.9225111 0.897538837539669
# F1.At          0.9887703 0.9357274 0.9531590 0.9531590 0.932709079670616
# F1.Dt          0.9873936 0.9304554 0.9499935 0.9499935 0.927982823362997
# Accuracy       0.9846604 0.8984014 0.9200404 0.9200404 0.897802853049512
# MCC            0.9712633 0.7890756 0.8563179 0.8563179  0.78602794499771


