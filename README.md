---
Title: "R2ROC"
Author: "Md Moksedul Momin and Hong Lee"
Date: "00/00/00" #to be updated
Output: pdf_document
---

# R2ROC
The ‘R2ROC’ package can be used to derive test statistics for AUC values from polygenic risk score (PRS) models (variance and covariance of AUC values, p-value, and 95% confidence intervals (CI)). For example, it can test if two sets of AUC values from two different PRS models are significantly different from each other and whether the two sets of PGS are independent or dependent. 

# INSTALLATION
To use R2ROC:
```
install.packages("devtools")
library(devtools)
devtools::install_github("mommy003/R2ROC")
library(R2ROC)
```
 or from CRAN
```
install.packages("R2ROC") 
library(R2ROC)
```
# QUICK START
We illustrate the usage of R2ROC using multiple sets of PRS estimated based on GWAS summary statistics from UK Biobank or Biobank Japan (reference datasets). In a target dataset, the phenotypes of target samples (y) can be predicted with PRS (a PRS model, e.g. y=PRS+e, where y and PRS are column-standardized 1 for pre adjusted phenotype. But for raw case-control data, phenotypes are assigned as 0 and 1. Note that the target individuals should be independent from reference individuals. We can test the significant differences in the AUC between a pair of PRS (see auc_var and auc_diff function and example in the manual).

# DATA PREPARATION
**a.	Statistical testing of significant difference between R2 values for p-value thresholds:** 
r2redux requires only phenotype and estimated PRS (from PLINK or any other software). Note that any missing value in the phenotypes and PRS tested in the model should be removed. If we want to test the significant difference of AUC values for two independent PRS, auc_diff function can be used with an input file that includes the following fields (please see dat1 and dat2 file embedded within the package and r2_diff function in the manual).
- Phenotype (y)
- PRS for discovery population 1 (x1)
- PRS for discovery population 2 (x2)
To get the test statistics for the difference between AUC(y=x[,v1]) and AUC(yx[,v2]). (here we define AUC= AUC(y=x[,v1])) and AUC=AUC(y=x[,v2])))

# References
1. Olkin, I. and  Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
2. DeLong, E.R., D.M. DeLong, and D.L. Clarke-Pearson, Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics, 1988: p. 837-845.
3. Heller, G., et al., Inference for the difference in the area under the ROC curve derived from nested binary regression models. Biostatistics, 2017. 18(2): p. 260-274.
4. Momin, M.M., Lee, S., Wray, N.R. and Lee S.h. 2022. Significance tests for R2 of out-of-sample prediction using polygenic scores. The American Journal of Human Genetics,110: p. 349-358. 

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Md Moksedul Momin (cvasu.momin@gmail.com) if you have any queries.
