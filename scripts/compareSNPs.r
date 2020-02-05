##################
## compare snps ##
##################
library(data.table)

### EAGLE-RC ###
system("head ../eagle/D.vs.A.raw.vcf")
# Gorai.001G000400.1	228	.	C	A	evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
# Gorai.001G000400.1	339	.	G	C	evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
# Gorai.001G000400.1	456	.	C	T	evm.model.Ga07G0005;rProp=0.7387;qProp=0.8388
system("wc -l ../eagle/D.vs.A.raw.vcf")
# 630,313 ../eagle/D.vs.A.raw.vcf 
system("wc -l ../eagle/A.vs.D.raw.vcf")
# 630309 ../eagle/A.vs.D.raw.vcf


### SNP index 4.0 ###
snpFile <- "D13.snp4.0.txt"
df = fread(snpFile,select=1:2)
snp=makeGRangesFromDataFrame(df, start.field="Pos", end.field="Pos")
length(snp) # 28540536
table(countOverlaps(query=snp, subject=gr))
#       0        1        2        3 
# 27288800  1245105     6614       17 
## a total of 1,251,736 SNPs

### HyLiTE
df = fread("../bowtie2hylite/resultsA2_011820/resultsA2_011820.snp.txt",select=5:7)
pa = paste(df$ADs,df$A2, df$D5,sep="@");table(pa)
# -1@-1@1   -1@0@1  -1@1@-1   -1@1@0   -1@1@1  0,0@0@1 0,0@1@-1  0,0@1@0 
#    43082    29626    44007    12231     9242     6355        1     3288 
#    0@0@1    0@1@0 1,0@-1@0 1,0@-1@1 1,0@0@-1  1,0@0@0  1,0@0@1 1,0@1@-1 
#    12482     2089    12753     4653    47264  1020631   692434    11407 
#  1,0@1@0  1,0@1@1   1@-1@1    1@0@1   1@1@-1    1@1@0    1@1@1 
#    62236    71953    12199     2928     4391      188    15964 
table(pa =="1,0@0@1" | pa =="1,0@1@0")
#   FALSE    TRUE 
# 1366734  754670
## a total of 754, 670 SNPs based on A reference

df = fread("/Volumes/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2hylite/results031417/results031417.snp.txt",select=5:7)
pa = paste(df$AD,df$A2, df$D5,sep="@");table(pa)
table(pa =="1,0@0@1" | pa =="1,0@1@0")
#   FALSE    TRUE 
# 2821354  871224
## a total of 871,224 SNPs based on D reference