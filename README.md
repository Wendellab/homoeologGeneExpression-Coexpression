# Challenges and pitfall in the use of partitioned gene counts for homoeologous gene expression and co-expression network analyses

## Input RNA-seq datasets
* Cotton seed development - 12 samples (4 time points x 3 biological replicates) per diploid or polyploid species, as deposited in NCBI BioProject [PRJNA179447](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA179447). SE reads.
* Flowering time regulation - 23 samples (8 tissue types x 3 biological replicates; only 2 reps for SDM) per diploid or polyploid species, as partially described [here](https://github.com/Wendellab/FloweringTimeDomestication/blob/master/sample.info). PE reads.

## RNA-seq mapping and homoeolog read estimation
* GNASP mapping followed by PolyCat homoeolog read participation: bash scritps for SE ([gsnap2polycat_120116.sh](scripts/gsnap2polycat_120116.sh)) and PE ([gsnap2polycat_PE.sh](scripts/gsnap2polycat_PE.sh)) reads.
* Bowtie2 mapping against reference transcripts followeed by RSEM read estimation: SE ([here](https://github.com/huguanjing/AD1_RNA-seq/blob/master/bowtie2rsem.sh)) and PE ([bowtie2rsem_PE.sh](scripts/bowtie2rsem_PE.sh)) bash scrpts.
* Bowtie2 mapping against D5 reference transcripts followed by HyLite SNP detection and read participation: SE bash stricpt ([bowtie2hylite.sh](scripts/bowtie2hylite.sh) with protocol file [se_protocol_file.txt](scripts/se_protocol_file.txt).
* SALMON mapping and read estimation against reference transcripts: [salmon.flower.sh](scripts/salmon.flower.sh) and [salmon.seed.sh](scripts/salmon.seed.sh)
* kallisto mapping and read estimation against reference transcripts: [kallisto.flower.sh](scripts/kallisto.flower.sh) and [kallisto.seed.sh](scripts/kallisto.seed.sh)

## Data analysis in R - workflow and scripts

Ongoing working dir: `/work/LAS/jfw-lab/hugj2006/eflen/output/`

LSS long term storage dir: `/lss/research/jfw-lab/Projects/Eflen/`

### Step 1. prepare read count tables

#### Scripts

* [0.prepare_Ranalysis.r](scripts/0.prepare_Ranalysis.r)
* [1a.read\_polycat\_datasets.r](scripts/1a.read_polycat_datasets.r)
* [1b.read\_RSEM\_datasets.r](scripts/1b.read_RSEM_datasets.r)
* [1c.read\_hylite\_datasets.r](scripts/1c.read_hylite_datasets.r)
* [1d.read\_salmon\_datasets.r](scripts/1d.read_salmon_datasets.r)
* [1e.read\_kallisto\_datasets.r](scripts/1e.read_kallisto_datasets.r)
* [1x.compare_datasets.r](scripts/1x.compare_datasets.r)

#### Output read count tables:

* PolyCat [raw counts](results/table.count.polycat.txt)
* HyLiTE [total read counts based D5 reference](results/table.count.hylite.total.txt) and [partitioned read counts in polyploid](results/table.count.hylite.ADs.txt)
* RSEM [estimated counts](results/table.count.rsem.txt) and [estimated RPKMs](results/table.rpkm.rsem.txt)
* Salmon [raw read counts](results/table.count.salmon.txt) and [tpm counts](results/table.tpm.salmon.txt)
* Kallisto [raw read counts](results/table.count.kallisto.txt) and [tpm counts](results/table.tpm.kallisto.txt)

#### Explanation of other output files

**[method]** as "polycat", "hylite", "rsem", "salmon", "kallisto".

* `pca.[method].log2cpm.pdf`............ PCA plot of 33 samples
* `R-01-[method]Datasets.RData`............ raw read counts of A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt
* `R-01-[method]NetworkDatasets.RData`............ dataset to build duplicated networks for A2D5 (expected from diploid true counts), A2D5.tech (technically expected from diploid counts accounting for program and reference errors), ADs (observed as estimated polyploid counts).
* `s1.plotVariance.A2D5vsADs.pdf`............ PCA plot of A2D5, A2D5.tech, ADs; not very useful.


### Step 2. Evalutation of homoeolog read estimation

Three metrics were calculated to evaluate the performance of homoeolog-specific read assignment - ***Efficiency***, ***Accuracy***, and ***Discrepancy***. 

For PolyCat and HyLiTE, ***Efficiency*** measures the proportion of total reads aligned to diploid reference that can be partitioned, while this measure for RSEM, Salmon and Kallisto approximates 1 given their different algorithms. 

***Accuracy*** measures the percentage of partitioned homoeolog read counts that are truely orginated from the predicted diploid genome. For PolyCat and HyLiTE, it is straight forward to look at alignments results `BAM` of partitioned reads and check whether their assignment agree with true origin. For RSEM, Salmon and Kallisto, no such "partitioned alignment" of polyploid can be inspected; therefore, we measured how accurate the _diploid_ reads were assigned to its corresponding subgenome based on alignment against the polyploid transcript reference.

***Discrepancy*** is like the percentage error of expected total diploid read counts, as abs(exp-obs)/exp.

Additionally from the perspective of [Precision and recall](https://en.wikipedia.org/wiki/Precision_and_recall), we obtained TP, TN, FP, FN from the diploid datasets for calculation. For example, At precision = A2.At/(A2.At + D5.At), recall = A2.At/(A2.At + A2.Dt); the recall is equivalent to previous ***Accuracy*** calculated for homoeologs. ***F measure*** combines precision and recall for evaluation: F= 2x(precision x recall) / (precision + recall).

#### Scripts

* [detectEffectiveRegion.r](effectiveRegion/detectEffectiveRegion.r): get effective transcript regions that are diagnostic of homoeolog origins, given certain RNA-seq fragment length - 100 bp for SE, 300 bp for PE here. 
* [2.0.get\_hylite_true.sh](scripts/2.0.get_hylite_true.sh): extract one-to-one mapping correspondence between ADs and diploid reads
* [2.1.evaluate\_read_assignment.r](scripts/2.evaluate_read_assignment.r) : metrics calculation
* [2.2.evaluation\_summary.r](scripts/2.evaluation_summary.r) : metrics comparison
* ####### **TO DO** ##########
* [2.3](): Statistical modeling and prediction by Meiling Liu; to be added

#### Explanation of output files

* `s2.eval.[method].pdf`............ histogram and pairs plot of metrics and etc.
* `s2.eval.[method].homoeolog.pdf`............ scatter plot of At vs Dt metrics
* `s2.eval.[method].summary.pdf`............ metric summary table
* `s2.assign_eval.[method].Rdata"
* `s2.assign_eval.summary.pdf`............ plots for comparing methods
* `2.2.evaluation_summary.r.txt`............ comparison analysis printout

**[method]** as "polycat", "hylite", "rsem", "salmon", "kallisto".




# bookmark, below NOT cleaned up yet

### Differential gene expression analysis
* [DE.r](DE033017.r) - Differential expression analysis using DESeq2 and EBSeq; the detection of 'true' DE genes according to reference dataset can be seen as a binary decision problem, where the performances of DE algorithms (DESeq2 vs EBSeq) in combination with homoeolog read estimation (polycat vs rsem vs hylite) were evaluated with ROC curves and AUC.  

#### Scripts

#### Output read count tables:

#### Explanation of other output files


### Differential gene-pair coexpression analysis
* [DC.all.r](DC.all.r) - Differential coexpression tests and classification for all gene pairs; scripts optimized for LARGE network (>40,000 genes) from DiffCorr and DGCA functions
* [DC.homoeoP.r](DC.homoeoP.r) - Differential coexpression tests and classification for At-Dt homoeo-pairs

#### Scripts

#### Output read count tables:

#### Explanation of other output files


### Coexpression network construction
* [Eflen_Networks110516.r](Eflen_Networks110516.r)

#### Scripts

#### Output read count tables:

#### Explanation of other output files


### Assessment of network topology and functional connectivity
* 

#### Scripts

#### Output read count tables:

#### Explanation of other output files

