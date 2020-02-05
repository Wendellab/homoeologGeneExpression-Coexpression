# Challenges and pitfall in the use of partitioned gene counts for homoeologous gene expression and co-expression network analyses

## Input RNA-seq datasets
* Cotton seed development - 12 samples (4 time points x 3 biological replicates) per diploid or polyploid species, as deposited in NCBI BioProject [PRJNA179447](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA179447). SE reads.
* Flowering time regulation - 23 samples (8 tissue types x 3 biological replicates; only 2 reps for SDM) per diploid or polyploid species, as partially described [here](https://github.com/Wendellab/FloweringTimeDomestication/blob/master/sample.info). PE reads.

### LSS location
* 33 SE seed RNA-seq files. `smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/Projects/Eflen/seed_for_eflen_paper/fastq/*fq.gz` **Exclude A2-20-R2.cut.fq.gz and D5-20-R2.cut.fq.gz.**
* 132 PE flowering time RNA-seq files. `smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/Projects/Eflen/flowerTimeDataset/fastq/*fq.gz`

* `smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/Projects/Eflen2020/flowerFastq`
* `smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/Projects/Eflen2020/seedFastq`


## RNA-seq mapping and homoeolog read estimation

### Specifically developed for polyploid systems

* [HyLite](https://hylite.sourceforge.io/index.html) automates Bowtie2 mapping followed by SNP detection and read participation [bowtie2hylite.sh](scripts/bowtie2hylite.sh)
    * A2 reference - slurm script [runBowtie2.A2TranscriptRef.slurm](scripts/runBowtie2.A2TranscriptRef.slurm) with protocol file [samA2\_protocol\_file.txt](scripts/samA2_protocol_file.txt)
    * D5 reference - slurm scripts [runBowtie2.D5TranscriptRef.slurm](scripts/runBowtie2.D5TranscriptRef.slurm) with protocol file [samD5\_protocol\_file.txt](scripts/samD5_protocol_file.txt).
* GNASP mapping followed by PolyCat homoeolog read participation: bash scritps for SE ([gsnap2polycat_120116.sh](scripts/gsnap2polycat_120116.sh)) and PE ([gsnap2polycat_PE.sh](scripts/gsnap2polycat_PE.sh)) reads.
* [EAGLE-rc](https://github.com/tony-kuo/eagle) classifies reads via calculating the likelihood of the read given Star alignments to A2 and D5 genome reference with scripts [eagle-rc.sh](scripts/eagle-rc.sh) and [runStarEagle-rc.slurm](scripts/runStarEagle-rc.slurm).

### Generic mapping tools

* [RSEM](http://deweylab.github.io/RSEM/) runs Bowtie2 mapping by default against polyploid reference transcripts followed by read estimation: SE ([here](https://github.com/huguanjing/AD1_RNA-seq/blob/master/bowtie2rsem.sh)) and PE ([bowtie2rsem_PE.sh](scripts/bowtie2rsem_PE.sh)) bash scrpts.
* [salmon](https://combine-lab.github.io/salmon/) mapping and read estimation against polyploid reference transcripts: [salmon.flower.sh](scripts/salmon.flower.sh) and [salmon.seed.sh](scripts/salmon.seed.sh)
* [kallisto](https://pachterlab.github.io/kallisto/) mapping and read estimation against reference transcripts: [kallisto.flower.sh](scripts/kallisto.flower.sh) and [kallisto.seed.sh](scripts/kallisto.seed.sh)
* We also applied the traditional approach using only Bowtie2 uniquely mapped reads against the polyploid transcriptome reference [runBowtie2only.slurm](scripts/runBowtie2only.slurm).

The A2D5 transcript sequences were used as reference, which were derived from the D5 gene models and A2-D5 SNP index (`smb://lss.its.iastate.edu/gluster-lss/research/jfw-lab/GenomicResources/pseudogenomes/A2D5.transcripts.fa`). 


## Data analysis in R - workflow and scripts

Ongoing working dir: `/work/LAS/jfw-lab/hugj2006/eflen/output/`

LSS long term storage dir: `/lss/research/jfw-lab/Projects/Eflen/`

### Step 1. prepare read count tables

#### Scripts

* [0.prepare_Ranalysis.r](scripts/0.prepare_Ranalysis.r)
* [1a.read\_polycat\_datasets.r](scripts/1a.read_polycat_datasets.r)
* [1b.read\_RSEM\_datasets.r](scripts/1b.read_RSEM_datasets.r)
* [1c.read\_hylite\_datasets.r](scripts/1c.read_hylite_datasets.r)
* [1c.read\_hyliteAref\_datasets.r](scripts/1c.read_hyliteAref_datasets.r)
* [1d.read\_salmon\_datasets.r](scripts/1d.read_salmon_datasets.r)
* [1e.read\_kallisto\_datasets.r](scripts/1e.read_kallisto_datasets.r)
* [1e.read\_eaglerc\_datasets.r](scripts/1e.read_eaglerc_datasets.r)
* [1e.read\_bowtie\_datasets.r](scripts/1e.read_bowtie_datasets.r)

#### Output read count tables:

* PolyCat [raw counts](results/table.count.polycat.txt)
* HyLiTE - D5 based [total read counts](results/table.count.hylite.total.txt) and [partitioned polyploid read counts](results/table.count.hylite.ADs.txt); A2 based [total read counts](results/table.count.hyliteAref.total.txt) and [partitioned polyploid read counts](results/table.count.hyliteAref.ADs.txt);
* RSEM [estimated counts](results/table.count.rsem.txt) and [estimated RPKMs](results/table.rpkm.rsem.txt)
* Salmon [raw read counts](results/table.count.salmon.txt) and [tpm counts](results/table.tpm.salmon.txt)
* Kallisto [raw read counts](results/table.count.kallisto.txt) and [tpm counts](results/table.tpm.kallisto.txt)
* Bowtie2 [raw read counts](results/table.count.bowtie.txt) 

#### Explanation of other output files

`[method]` as "polycat", "hylite", "rsem", "salmon", "kallisto".

* `pca.[method].log2cpm.pdf`............ PCA plot of 33 samples
* `R-01-[method]Datasets.RData`............ raw read counts of A2.Total, A2.At, A2.Dt, D5.Total, D5.At, D5.Dt, ADs.At, ADs.Dt
* `R-01-[method]NetworkDatasets.RData`............ dataset to build duplicated networks for A2D5 (expected from diploid true counts), A2D5.tech (technically expected from diploid counts accounting for program and reference errors), ADs (observed as estimated polyploid counts).


### Step 2. Evalutation of homoeolog read estimation

Knowing the true subgenome origin of each *in silico* polyploid ADs reads, we obtained the confusion matrix (TP, TN, FP, FN) to evaluate the homoeolog read classification. Several metrics including ***Precision/recall***, ***F1 score***, ***MCC***, ***Accuracy***metrics were calculated. 

In addition, custom measures of ***Efficiency*** and ***Discrepancy*** were used to examine how different methods deal with ambiguous read alignment (discard or statistical inference). For PolyCat and HyLiTE, ***Efficiency*** measures the proportion of total reads aligned to diploid reference that can be partitioned, while this measure for RSEM, Salmon and Kallisto approximates 1 given their different algorithms. 


#### Scripts

* [detectEffectiveRegion.r](effectiveRegion/detectEffectiveRegion.r): get effective transcript regions that are diagnostic of homoeolog origins, given certain RNA-seq fragment length - 100 bp for SE, 300 bp for PE here. ***Ambiguity* = 1 - %effective_region**
* [2.0.get\_hylite_true.sh](scripts/2.0.get_hylite_true.sh): extract one-to-one mapping correspondence between ADs and diploid reads
* [2.1.evaluate\_read_assignment.r](scripts/2.1.evaluate_read_assignment.r) : metrics calculation
* [2.2.evaluation\_summary.r](scripts/2.2.evaluation_summary.r): compare metrics and make summary table.

#### Explanation of output files

* `s2.eval.[method].pdf`............ histogram and pairs plot of metrics and etc.
* `s2.eval.[method].homoeolog.pdf`............ scatter plot of At vs Dt metrics
* `s2.eval.[method].summary.pdf`............ metric summary table
* `s2.assign_eval.[method].Rdata`............ metrics values
* `s2.evaluation_summary.txt`............ comparison analysis results

`[method]` as "polycat", "hylite", "rsem", "salmon", "kallisto".

### Step 3. Differential gene expression analysis

In comparison with the "expected" differentially expressed genes between homoeologs (A2 vs D5), we ask how homoeolog read estimation and DE methods tegoether affect the "observed" (ADs: At vs Dt) lists of DE genes. Two DE analysis algorithms - DESeq2 and EBSeq, in conjuction with each of the five homoeolog read estimation methods (polycat, hylite, rsem, salmon, kallisto) were tested. The detection of "Expected" DE genes can be seen as a binary decision problem, which were evaluated with Sensitivity (=recall), Specificity, Precision, F statistics, MCC, ROC curves and AUC. 

#### Scripts

* [3.differential_expression.r](scripts/3.differential_expression.r) 

#### Explanation of output files

* `s3.DE.summary.txt`............ summary table of DE gene numbers between A vs D.
* `s3.DE.summary.pdf`............ histogram and pairs plot of metrics and etc.
* `s3.DE.evals.txt`............ summary table of Sensitivity (=recall), Specificity, Precision, F statistics, and MCC.
* `s3.DE.evals.pdf`............ scatter plot of observed vs. expected log2FoldChanges.
* `s3.DE.ROC.txt`............ summary table of ROC curves and AUC.
* `s3.DE.ROC.pdf`............ ROC plots
* `s3.DE.performance.pdf`............ boxplot of DE number and binary classification metrics. **Important**
* `s3.anovaTests.rout.txt`............ ANOVA test results.

### Step 4. Differential gene-pair coexpression analysis

Coexpression of homoeologs and between all possible gene pairs were measured by Pearson's coefficients and then classified by contrasting *estimated* versus *true* patterns. Nine classes of DC patterns were resulted and tested for enrichment.

#### Scripts
* [4.DC.homoeoP.r](scripts/4.DC.homoeoP.r) 
* [4.DC.all.r](scripts/4.DC.all.r) 

#### Explanation of output files
* `s4.DC.homoeoPair.Rdata`............ result tables of homoeolog DC analysis for each pipeline normalized by rld and log2rpkm.
* `s4.DC.homoeoPair.pdf`............ enrichment analysis of DC categories between homoeolog gene pairs.
* `s4.DC.all.Rdata`............ summary result tables of all gene pairs DC analysis for each pipeline normalized by rld and log2rpkm.
* `s4.DC.all.pdf`............ enrichment analysis of DC categories for all gene pairs.

###  Step 5. Coexpression network construction

True and estimated homoeolog read count tables from 10 datasets (five mapping pipelines followed by rld or log2rpkm transformation) were subjected to weighted and unweighted coexpression network construction. The same set of genes were included in all networks for fair comparison.

#### Scripts

* [5a.WGCNA.prep.r](scripts/5a.WGCNA.prep.r) 
* [5a.WGCNA.r](scripts/5a.WGCNA.r) ([5a.slurm](scripts/5a.slurm))
* [5a.WGCNA.post.r](scripts/5a.WGCNA.post.r)
* [5b.BNA.r](scripts/5b.BNA.r) 

#### Explanation of output files

**WGCNA**: `[method]` as polycat, hylite, rsem, salmon or kallisto; `[transfomation]` as rld or log2rpkm.

* `R-05-dataInput.[method]_[transfomation].RData`............ network input multiExpr with rlog or log2rpkm transformation
* `s5.multiExpr.clustering.pdf`............ clustering analysis of network input multiExpr.
* `R-05-chooseSoftThreshold.Rdata`............ network input multiExpr with log2rpkm transformation
* `s5.wgcna.choosePower_[method]_[transfomation].pdf`............ plot of choose soft threshold
* `s5.refPresevation_[method]_[transfomation].B[power].pdf`............ preservation test results
* `s5.Zsummary.pdf`............ summary boxplot of Z preservation results


**Binary networks**

* `5.BNA.rout.txt`............ log history of binary network analysis.
* `s5.bna.AUC.Rdata`............ AUC results
* `s5.bna.plotAUC.sample6266.permutation10.pdf`............ plot of AUC results


### Step 6. Assessment of network topology and functional connectivity

#### Scripts

* [6.prepFunctionCategory.r](scripts/6.prepFunctionCategory.r)
* [6.FUN.r](scripts/6.FUN.r) 
* [6a.NC.r](scripts/6a.NC.r) 
* [6b.FC.r](scripts/6b.FC.r)
* [6.post.r](scripts/6.post.r)

#### Explanation of output files

* `GOnKEGGnGLs.Rdata`............ functional categories of GO, oil related gene families, and flowering time related gene families.
* -------------------- Node Connectivity --------------------
* `s6.NC.rdata`............ results of node connectivity correlation (exp vs obs) and A-/D-subnetwork density.
* `s6.NC_corr.pdf`............ boxplot of node connectivity correlation (exp vs obs) across pipeline, normalization and network Types.
* `s6.NC_denisty.pdf`............ line plot with error bars of subnetwork density for each combination of pipeline, normalization and network Types.
* `s6.NC_denisty.txt`............ result table of subnetwork density for each combination of pipeline, normalization and network Types.
* -------------------- Functional Connectivity --------------------
* `s6.FC.rdata`............ results of AUROC correlations (exp vs obs) and A-/D-subnetwork AUROCs.
* `s6.FC_corr.pdf`............ boxplot of functional connectivity correlation (exp vs obs) across pipeline, normalization and network Types.
* `s6.FC_auroc.pdf`............ boxplot of functional connectivity across pipeline, normalization and network Types.
* `s6.FC_subAD.pdf`............ ............ line plot with error bars of subnetwork auroc for each combination of pipeline, normalization and network Types.
* `s6.FC_subAD.txt`............ result table of subnetwork auroc for each combination of pipeline, normalization and network Types.



### Step 7. Examination of the impact of read ambiguity on performance

The metrics derived from read assignment, DE, DC and network analyses were correlayed with gene groups binned by ambiguity.

#### Scripts

* [7.readAmbiguity.r](scripts/7.readAmbiguity.r)

#### Explanation of output files

*  `s7.eflen.rdata`............ summary results of quantification, DE, DC, and network performance with respected to ambiguity levels.
