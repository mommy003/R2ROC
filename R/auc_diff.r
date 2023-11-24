#' auc_diff function
#'
#' This function estimates var(AUC(y~x\[,v1]) - AUC(y~x\[,v2]))
#' where AUC is the Area Under ROC curve of the model,
#' y is N by 1 matrix having the dependent variable, and
#' x is N by M matrix having M explanatory variables.
#' v1 or v2 indicates the ith column in the x matrix 
#' (v1 or v2 can be multiple values between 1 - M, see Arguments below)
#' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
#' @param v1 This can be set as v1=c(1) or v1=c(1,2)
#' @param v2 This can be set as v2=c(2), v2=c(3), v2=c(1,3) or v2=c(3,4)
#' @param nv Sample size
#' @param kv Population prevalence
#' @keywords R2 variance information matrix
#' 
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will estimate significant difference between two PRS (either dependent or independent and joint or single). To get the test statistics for the difference between AUC(y~x\[,v1]) and AUC(y~x\[,v2]) (here we define AUC1=AUC(y~x\[,v1])) and AUC2=AUC(y~x\[,v2]))). The outputs are listed as follows.
#' \item{mean_diff}{AUC differences between AUC1 and AUC2}
#' \item{var}{Variances of AUC differences}
#' \item{upper_diff}{Upper value of the differences}
#' \item{lower_diff}{Upper value of the differences}
#' \item{p}{Two tailed P-value for significant difference between AUC1 and AUC2}
#' \item{p_one_tail}{One tailed P-value for significant difference}
#' \item{heller_p}{P-value based on Heller's test for significant difference}
#' \item{heller_upper_diff}{Upper limit of 95% CI for the difference basedon Heller's test}
#' \item{heller_lower_diff}{Lower limit of 95% CI for the difference basedon Heller's test}
#' @examples
#' #To get the test statistics for the difference between AUC(y=x[,1]) 
#' #and AUC(y=x[,2])
#' dat=dat1 #(this example embedded within the package)
#' nv=length(dat$V1)
#' kv=sum(dat$V1)/length(dat$V1)# pop. prevalence estimated from data
#' #R2ROC also allows users to estimate AUC using pre-adjusted phenotype
#' #In that case, users need to specify kv
#' #eg. kv=0.10 for dat2 (dat2 embedded within the package) 
#' v1=c(1)
#' v2=c(2)
#' output=auc_diff(dat,v1,v2,nv,kv)
#' 
#' #R2ROC output
#' #output$mean_diff (mean difference of AUC1 and AUC2)
#' #0.1756046
#' 
#' #output$var (variance of AUC difference)
#' #9.274356e-05
#' 
#' #output$upper_diff (upper limit of 95% CI for difference)
#' #0.1944801
#' 
#' #output$lower_diff (lower limit of 95% CI for difference)
#' #0.1567292
#' 
#' #output$p (two-tailed P-value for the differences is
#' #significantly different from zero)
#' #2.747031e-74
#' 
#' #output$p_one_tail (one-tailed P-value for the differences
#' #is significantly different from zero)
#' #1.373515e-74
#' 
#' 
#' #To get the test statistics for the difference between
#' #AUC(y=x[,1]+x[,2]) and AUC(y=x[,2])
#' dat=dat1 #(this example embedded within the package)
#' nv=length(dat$V1)
#' kv=sum(dat$V1)/length(dat$V1)# pop. prevalence estimated from data
#' #R2ROC also allows users to estimate AUC using pre-adjusted phenotype
#' #In that case, users need to specify kv
#' #eg. kv=0.10 for dat2 (dat2 embedded within the package) 
#' v1=c(1,2)
#' v2=c(2)
#' output=auc_diff(dat,v1,v2,nv,kv)
#' 
#' #R2ROC output
#' #output$mean_diff (mean difference of AUC1 and AUC2)
#' #0.1793682
#' 
#' #output$var (variance of AUC difference)
#' #0.0001190366
#' 
#' #output$upper_diff (upper limit of 95% CI for difference)
#' #0.2007526
#' 
#' #output$lower_diff (lower limit of 95% CI for difference)
#' #0.1579839
#' 
#' #output$p (two-tailed P-value for the differences is
#' #significantly different from zero)
#' #9.87014e-61
#' 
#' #output$p_one_tail (one-tailed P-value for the differences
#' #is significantly different from zero)
#' #4.93507e-61
#' 
#' #output$heller_p (two-tailed P-value based on Hellers test
#' #for the differences is significantly different from zero)
#' #4.2085e-237
#' 
#' #output$heller_upper_diff (upper limit of 95% CI for
#' #difference based on Hellers test)
#' #0.2013899
#' 
#' #output$heller_lower_diff (lower limit of 95% CI for
#' #difference based on Hellers test)
#' #0.1586212


auc_diff = function (dat,v1,v2,nv,kv) {
  
  omat=cor(dat)
  
  if (length(v1)==1 & length(v2)==1) {
    ord=c(1,(1+v1),(1+v2)) 
    
    aoa=olkin_auc1_2(omat[ord,ord],nv,kv)
    chi_dum=aoa$diff^2/aoa$var
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$diff+1.96*aoa$var^.5
    lci=aoa$diff-1.96*aoa$var^.5
    
    z=list(mean_diff=aoa$diff,var=aoa$var,upper_diff=uci,lower_diff=lci,p=p3,p_one_tail=p3/2)
    return(z)
    
  }
  
  if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==2) {
    if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]))}
    if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]))}
    
    aoa=olkin_auc12_1(omat[ord,ord],nv,kv)
    chi_dum=aoa$diff^2/aoa$var
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$diff+1.96*aoa$var^.5
    lci=aoa$diff-1.96*aoa$var^.5
    
    #Following Heller et al. (Biostatistics 2017, 18: 260)
    chi_dum=aoa$diff/(aoa$var/(4*aoa$diff))   #see p265 in Heller et al.
    heller_p=pchisq(chi_dum,1,lower.tail=F)
    heller_uci=(aoa$diff^.5+1.96*(aoa$var/(4*aoa$diff))^.5)^2
    heller_lci=(aoa$diff^.5-1.96*(aoa$var/(4*aoa$diff))^.5)^2
    
    z=list(mean_diff=aoa$diff,var=aoa$var,upper_diff=uci,lower_diff=lci,p=p3,p_one_tail=p3/2,heller_p=heller_p,heller_upper_diff=heller_uci,heller_lower_diff=heller_lci)
    return(z)
    
  }
  
  if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==3) {
    ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))
    
    aoa=olkin_auc12_3(omat[ord,ord],nv,kv)
    chi_dum=aoa$diff^2/aoa$var
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$diff+1.96*aoa$var^.5
    lci=aoa$diff-1.96*aoa$var^.5
    
    z=list(mean_diff=aoa$diff,var=aoa$var,upper_diff=uci,lower_diff=lci,p=p3,p_one_tail=p3/2)
    return(z)
    
  }
  
  if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==3) {
    if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[2]))}
    if (v1[1]==v2[2]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))}
    if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[2]))}
    if (v1[2]==v2[2]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[1]))}
    
    aoa=olkin_auc12_13(omat[ord,ord],nv,kv)
    chi_dum=aoa$diff^2/aoa$var
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$diff+1.96*aoa$var^.5
    lci=aoa$diff-1.96*aoa$var^.5
    
    z=list(mean_diff=aoa$diff,var=aoa$var,upper_diff=uci,lower_diff=lci,p=p3,p_one_tail=p3/2)
    return(z)
    
  }
  
  
  if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==4) {
    ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]),(1+v2[2]))
    
    aoa=olkin_auc12_34(omat[ord,ord],nv,kv)
    chi_dum=aoa$diff^2/aoa$var
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$diff+1.96*aoa$var^.5
    lci=aoa$diff-1.96*aoa$var^.5
    
    z=list(mean_diff=aoa$diff,var=aoa$var,upper_diff=uci,lower_diff=lci,p=p3,p_one_tail=p3/2)
    return(z)
    
  }
  
  
}
