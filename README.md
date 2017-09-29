## Challenges and pitfall in the use of partitioned gene counts for homoeologous gene expression and co-expression network analyses

---
### RNA-seq datasets
* Cotton seed development - 12 samples (4 time points x 3 biological replicates) per diploid or polyploid species, as deposited in NCBI BioProject [PRJNA179447](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA179447). SE reads.
* Flowering time regulation - 23 samples (8 tissue types x 3 biological replicates; only 2 reps for SDM) per diploid or polyploid species, as partially described [here](https://github.com/Wendellab/FloweringTimeDomestication/blob/master/sample.info). PE reads.


### RNA-seq mapping and homoeolog read estimation
* GNASP mapping followed by PolyCat homoeolog read participation: bash scritps for [SE](gsnap2polycat_120116.sh) and [PE](gsnap2polycat_PE.sh) reads.
* Bowtie2 mapping against reference transcripts followeed by RSEM read estimation: [SE](https://github.com/huguanjing/AD1_RNA-seq/blob/master/bowtie2rsem.sh) and [PE](bowtie2rsem.sh) bash scrpts.

* Bowtie2 mapping aginst D5 reference transcripts followed by HyLite SNP detection and read participation: [SE](bowtie2hylite.sh) bash stricpt with [this](sam2_protocol_file.txt) protocol file.

### Evalutation of homoeolog read estimation
* Detection of effective transcript regions that are diagnostic of homoeolog origins: [R script](detectEffectiveRegion.r)
* Calculating measures of ***Efficiency***, ***Accuracy***, and ***Discrepancy***, to evaluate the performance of homoeolog-specific read assignment: to be added
* Statistical modeling and prediction: to be added


### Differential gene expression analysis
* [DE.r](DE033017.r) - Differential expression analysis using DESeq2 and EBSeq; the detection of 'true' DE genes according to reference dataset can be seen as a binary decision problem, where the performances of DE algorithms (DESeq2 vs EBSeq) in combination with homoeolog read estimation (polycat vs rsem vs hylite) were evaluated with ROC curves and AUC.  

### Differential gene-pair coexpression analysis
* [DC.all.r](DC.all.r) - Differential coexpression tests and classification for all gene pairs; scripts optimized for LARGE network (>40,000 genes) from DiffCorr and DGCA functions
* [DC.homoeoP.r](DC.homoeoP.r) - Differential coexpression tests and classification for At-Dt homoeo-pairs

### Coexpression network construction
* [Eflen_Networks110516.r](Eflen_Networks110516.r)

### Assessment of network topology and functional connectivity
* 
