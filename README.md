## Challenges and pitfall in the use of partitioned gene counts for homoeologous gene expression and co-expression network analyses
---

### RNA-seq mapping and homoeolog read estimation
* [gsnap2polycat_120116.sh](gsnap2polycat_120116.sh) - GNASP mapping followed by PolyCat homoeolog read participation.
* [bowtie2hylite.sh](bowtie2hylite.sh) - Bowtie2 mapping followed by HyLite SNP detection and read participation, using [this](sam2_protocol_file.txt) protocol file.
* [bowtie2rsem.sh](https://github.com/huguanjing/AD1_RNA-seq/blob/master/bowtie2rsem.sh)- Bowtie2 mapping with RSEM read estimation.


### Evalutation of homoeolog read estimation
* 

### Differential gene expression analysis
* [DE.r](DE033017.r) - Differential expression analysis using DESeq2 and EBSeq; the detection of 'true' DE genes according to reference dataset can be seen as a binary decision problem, where the performances of DE algorithms (DESeq2 vs EBSeq) in combination with homoeolog read estimation (polycat vs rsem vs hylite) were evaluated with ROC curves and AUC.  

### Differential gene-pair coexpression analysis
* [DC.all.r](DC.all.r) - Differential coexpression tests and classification for all gene pairs; scripts optimized for LARGE network (>40,000 genes) from DiffCorr and DGCA functions
* [DC.homoeoP.r](DC.homoeoP.r) - Differential coexpression tests and classification for At-Dt homoeo-pairs

### Coexpression network construction
* 

### Assessment of network topology and functional connectivity
* 
