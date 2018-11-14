###################################
## Prepare Fucntional Annotation ##
###################################
# scricpts can be saved separately as "prep.annot.r"
# source("prep.annot.r")

go<-read.table("G.raimondii_JGI_221_v2.1_genes2GO.txt",sep="\t") #166616
go<-go[grepl("Gorai[.]0.*[.]1$",go$V1),1:2]
goterms <- unique(go[,2])
dim(go) # 75336     2
length(goterms) # 1710
go.AD <- data.frame(gene=c(gsub("[.]1$","a",go$V1),gsub("[.]1$","d",go$V1)), GO=rep(go$V2,2))

kegg<-read.table("G.raimondii_JGI_221_v2.1_KEGG.pathways.txt",sep="\t") #30596
kegg<-kegg[grepl("Gorai[.]0.*[.]1$",kegg$V1),]
kegg.AD <- data.frame(gene=c(gsub("[.]1$","a",kegg$V1),gsub("[.]1$","d",kegg$V1)), KEGG=rep(kegg$V2,2), desc=rep(kegg$V3,2))
keggterms <- unique(kegg.AD[,2])
length(keggterms) #305

# Oil related gene list derived from Hovav et al, TPG (2015) supplemental table S3
# https://dl.sciencesocieties.org/publications/tpg/abstracts/8/1/plantgenome2014.08.0041
FAs<-read.table("FAs.txt",header=TRUE,sep="\t")
dim(FAs) #657
length(FAgenes<-unique(FAs$nodeName)) # only 657 unique genes, oops
fa<-FAs[,1:2]
names(fa)<-c("gene","FA.bio.process")
fa.AD<- data.frame(gene=c(paste0(fa$gene,"a"),paste0(fa$gene,"d")), FA=rep(fa$FA.bio.process,2))

# Flowering time regulation related gene list derived from Grover et al, TPG (2015) supplemental table S1
# https://dl.sciencesocieties.org/publications/tpg/articles/8/2/plantgenome2014.12.0098
FTs<-read.table("FTs.txt",header=TRUE,sep="\t")
dim(FTs) #657
length(FTgenes<-unique(FTs$gene)) # only 657 unique genes, oops
ft<-FTs[,1:2]
names(ft)<-c("gene","FT.fun.category")
ft.AD<- data.frame(gene=c(paste0(ft$gene,"a"),paste0(ft$gene,"d")), FT=rep(ft$FT.fun.category,2))

save(go.AD, kegg.AD, ft.AD, fa.AD, file="GOnKEGGnGLs.Rdata")