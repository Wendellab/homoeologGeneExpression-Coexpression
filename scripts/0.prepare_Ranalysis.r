## Analysis by Guanjing Hu from Nov 5th, 2016
## updated from previous version: Eflen_Networks021616.r, Eflen_Networks051916.r, Eflen_Networks090816.r

ssh hugj2006@bigram.ent.iastate.edu
cd /home/hugj2006/jfw-lab/Projects/Eflen_networks
ln -s /home/hugj2006/jfw-lab/Papers/GroverEflen2015/count_files count_files021616
screen -S eflen
module load R
# Module name: R                          Version: 3.3.1
R
# start R analysis



## The analysis was coducted in following steps
## 1. Basic data processing and cleaning - input DESeq2 rlog table and trait table for sample clustering and detecting outliers.
## 2. Choosing the soft-thresholding power - default or powers making good fit of scale free topology.
## 3. Network Construction - single block, corType = "pearson" (not "bicor", see discussion below), networkType = "signed"
## 4. General network topology analysis - produce a few plots for exploration


############## Install WGCNA ####################
source("http://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
biocLite("impute")
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("WGCNA")
install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") )


sessionInfo()