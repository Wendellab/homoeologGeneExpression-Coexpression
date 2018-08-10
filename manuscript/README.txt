# README

This analysis was conducted by Guanjing Hu. Files were last updated on Mar 18th, 2018.

* s5.rldExpr_all.txt - input file used
* differential_coexpression.r - r scripts
* s5.DC.all.Rdata - output rdata
* s5.DC.classes.pdf - output plot for DC gene pairs
* s5.DC.genes.txt - output table for DC genes

## Description

Briefly, differential co-expression (DC) analysis was performed by calculation Pearson correlation coefficients for all gene pairs followed by comparisons of corresponding gene pairs between wild and domesticated cotton fiber datasets. Differential correlation were tested based on Fisher's z-test using the R package DiffCorr. Categorization of differential correlation changes were conducted using the R package DGCA. 

Given that 3 wild and 3 domesticated accessions were used, I considered ALL 6 combinations of wild vs domesticated datasets to conduct DC analysis. That is, domesticated datasets in the order of CRB252, Maxxa and TM1 (each with 5, 10, 15, 20 dpa order), were compared to 6 order combinations of wild accessions as:
  
- [1,] "TX2095" "TX665"  "Yuc"   
- [2,] "TX2095" "Yuc"    "TX665" 
- [3,] "TX665"  "TX2095" "Yuc"   
- [4,] "TX665"  "Yuc"    "TX2095"
- [5,] "Yuc"    "TX2095" "TX665" 
- [6,] "Yuc"    "TX665"  "TX209

It turns out that all 6 combinations led to the same results: 6.24% of gene pairs (81091216 out of 1300270510) were differentially co-expressed between wild and domesticated cotton fibers.

Based on this average percentage of DC (p = 6.24%) among all possible gene pairs, differentially co-expressed genes were identified if its number of DC pairs was significantly higher than expected. That is, for a gene identified with k DC pairs among its all gene pairs n, the probability P of this gene to be significantly coexpression follows the binomial distribution model. The resulted P values were further corrected by the BH methods for multiple testing correction.

As shown in "s5.DC.genes.txt", 32.3% of genes (16503 out of 50996) were differentially co-expressed between wild and domesticated cotton fibers. We may next ask what are these genes? How many DC genes were found in homoeolog gene pairs, and how many were found as only At or Dt genes? Do they also exhibit different expression patterns between wild and domesticated fibers? etc.
 