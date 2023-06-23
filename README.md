---
Title: "GxEprs"
Author: "Dovini Jayasinghe, Md Moksedul Momin and Hong Lee"
Date: "00/00/00" #to be updated
Output: pdf_document
---


# R2ROC

The ‘R2ROC’ package can be used to derive test statistics for AUC values from polygenic risk score (PRS) models (variance and covariance of AUC values, p-value and 95% confidence intervals (CI)). For example, it can test if two sets of AUC values from two different PRS models are significantly different to each other whether the two sets of PGS are independent or dependent. used in the r2redux.  

# INSTALLATION
To use R2ROC:
```
install.packages("R2ROC") 
library(R2ROC)
```
 or
```
install.packages("devtools")
library(devtools)
devtools::install_github("mommy003/R2ROC")
library(R2ROC)
```

# References
1. Olkin, I. and  Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
2. DeLong, E.R., D.M. DeLong, and D.L. Clarke-Pearson, Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics, 1988: p. 837-845.
3. Heller, G., et al., Inference for the difference in the area under the ROC curve derived from nested binary regression models. Biostatistics, 2017. 18(2): p. 260-274.
4. Momin, M.M., Lee, S., Wray, N.R. and Lee S.h. 2022. Significance tests for R2 of out-of-sample prediction using polygenic scores. The American Journal of Human Genetics,110: p. 349-358. https://doi.org/10.1016/j.ajhg.2023.01.004 .

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Md Moksedul Momin (cvasu.momin@gmail.com) if you have any queries.
