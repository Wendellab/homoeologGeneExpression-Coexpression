# Challenges and pitfall in the use of partitioned gene counts for homoeologous gene expression and co-expression network analyses

## Input RNA-seq datasets
* Cotton seed development - 12 samples (4 time points x 3 biological replicates) per diploid or polyploid species, as deposited in NCBI BioProject [PRJNA179447](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA179447). SE reads.
* Flowering time regulation - 23 samples (8 tissue types x 3 biological replicates; only 2 reps for SDM) per diploid or polyploid species, as partially described [here](https://github.com/Wendellab/FloweringTimeDomestication/blob/master/sample.info). PE reads.

## Data analysis workflow and scripts

### RNA-seq mapping and homoeolog read estimation
* GNASP mapping followed by PolyCat homoeolog read participation: bash scritps for SE ([gsnap2polycat_120116.sh](scripts/gsnap2polycat_120116.sh)) and PE ([gsnap2polycat_PE.sh](scripts/gsnap2polycat_PE.sh)) reads.
* Bowtie2 mapping against reference transcripts followeed by RSEM read estimation: SE ([here](https://github.com/huguanjing/AD1_RNA-seq/blob/master/bowtie2rsem.sh)) and PE ([bowtie2rsem_PE.sh](scripts/bowtie2rsem_PE.sh)) bash scrpts.
* Bowtie2 mapping against D5 reference transcripts followed by HyLite SNP detection and read participation: SE bash stricpt ([bowtie2hylite.sh](scripts/bowtie2hylite.sh) with protocol file [se_protocol_file.txt](scripts/se_protocol_file.txt).
* SALMON mapping and read estimation against reference transcripts: [salmon.flower.sh](scripts/salmon.flower.sh) and [salmon.seed.sh](scripts/salmon.seed.sh)
* kallisto [???]()


### Evalutation of homoeolog read estimation
* Detection of effective transcript regions that are diagnostic of homoeolog origins: [R script](effectiveRegion/detectEffectiveRegion.r)
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

## Result datasets

### Read count tables

* [polyCat raw counts](results/table.count.polycat.txt)
* [HyLiTE total read counts based D5 reference](results/table.count.hylite.total.txt)
* [HyLiTE partitioned read counts in polyploid](results/table.count.hylite.ADs.txt)* [RSEM estimated counts](results/table.count.rsem.txt)
* [RSEM estimated RPKMs](results/table.rpkm.rsem.txt)
* [Salmon raw read counts](results/table.count.salmon.txt)
* [Salmon tpm counts](results/table.tpm.salmon.txt)

### Explanation of other output files

Working dir: work/LAS/

Long term storage dir: lss/

**[method]** as "polycat", "hylite", "rsem", "salmon"

* pca.**[method]**.log2cpm.pdf
* R-01-**[method]**Datasets.RData
* R-01-**[method]**NetworkDatasets.RData

