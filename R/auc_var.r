  #' auc_var function
  #'
  #' This function estimates var(AUC(y~x\[,v1]))
  #' where AUC is the Area Under ROC curve of the model,
  #' y is N by 1 matrix having the dependent variable, and
  #' x is N by M matrix having M explanatory variables.
  #' v1 indicates the ith column in the x matrix
  #' (v1 can be multiple values between 1 - M, see Arguments below)
  #' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
  #' @param v1 This can be set as v1=c(1), v1=c(1,2) or possibly with more values
  #' @param nv Sample size
  #' @param kv Population prevalence
  #' @keywords R2 variance information matrix
  #' @export
  #' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
  #' @return  This function will test the null hypothesis for AUC. To get the test statistics for AUC(y~x\[,v1]). The outputs are listed as follows.
  #' \item{auc}{AUC}
  #' \item{var}{Variance of AUC}
  #' \item{upper_auc}{Upper limit of 95% CI for AUC}
  #' \item{lower_auc}{Lower limit of 95% CI for AUC}
  #' @examples
  #' #To get the AUC for AUC(y=x[,v1]) 
  #'
  #' dat=dat1 #(this example embedded within the package)
  #' nv=length(dat$V1)
  #' kv=sum(dat$V1)/length(dat$V1)# pop. prevalence estimated from data
  #' #R2ROC also allows users to estimate AUC using pre-adjusted phenotype
  #' #In that case, users need to specify kv
  #' #eg. kv=0.10 for dat2 (dat2 embedded within the package) 
  #' v1=c(1)
  #' output=auc_var(dat,v1,nv,kv)
  #' 
  #' #R2ROC output
  #' #output$auc (AUC)
  #' #0.7390354
  #' 
  #' #output$var (variance of AUC)
  #' #7.193337e-05
  #' 
  #' #output$upper_auc (upper limit of 95% CI for AUC)
  #' #0.7556589
  #' 
  #' #output$lower_auc (lower limit of 95% CI for AUC)
  #' #0.7224119



  auc_var = function (dat,v1,nv,kv) {
  
  dat=scale(dat);omat=cor(dat)
  
  ord=c(1,(1+v1)) 

  if (length(v1)==1) {

    aoa=olkin_auc1(omat[ord,ord],nv,kv)
    #chi_dum=aoa$auc^2/aoa$var
    #p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$auc+1.96*aoa$var^.5
    lci=aoa$auc-1.96*aoa$var^.5
    
    #z=list(auc=aoa$auc,var=aoa$var,upper_auc=uci,lower_auc=lci,p=p3,p_one_tail=p3/2)
    z=list(auc=aoa$auc,var=aoa$var,upper_auc=uci,lower_auc=lci)
    return(z)
    
  }

  if (length(v1)==2) {

    aoa=olkin_auc12(omat[ord,ord],nv,kv)
    #chi_dum=aoa$auc^2/aoa$var
    #p3=pchisq(chi_dum,1,lower.tail=F)
    uci=aoa$auc+1.96*aoa$var^.5
    lci=aoa$auc-1.96*aoa$var^.5
    
    #z=list(auc=aoa$auc,var=aoa$var,upper_auc=uci,lower_auc=lci,p=p3,p_one_tail=p3/2)
    z=list(auc=aoa$auc,var=aoa$var,upper_auc=uci,lower_auc=lci)
    return(z)
    
  }



  
}












