#' auc_trf function
#'
#' This function transforms the observed scale predictive ability (R2)
#' and its standard error (SE) to AUC with its SE
#' @references
#' Wray, Naomi R., et al. "The genetic interpretation of area under the ROC curve in genomic profiling." PLoS genetics 6.2 (2010): e1000864.
#' @references
#' Lee, Sang Hong, et al. "A better coefficient of determination for genetic profile analysis." Genetic epidemiology 36.3 (2012): 214-224.
#' @param R2 R2 or coefficient of determination on the observed scale
#' @param se Standard error of R2
#' @param kv Population prevalence
#' @keywords Transformation of observed R2 to AUC
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will transform the observed R2 and its s.e between to AUC. Output from the command is the lists of outcomes.
#' \item{auc}{Transformed AUC}
#' \item{se}{SE of transformed AUC}
#' @examples
#' #To get the transformed AUC
#'
#' output=auc_trf(0.04, 0.002, 0.05)
#' output
#'
#' #output$auc (transformed AUC)
#' #0.7522887
#'
#' #output$se (se of transformed AUC)
#' #0.005948364


auc_trf = function (R2,se,kv) {
  
  thd=-qnorm(kv)
  zv=dnorm(thd)
  iv=zv/kv
  iv2 = -iv*kv/(1-kv)
  cv=kv*(1-kv)/zv^2
  
  r1=sqrt(R2)
  f=expression(pnorm( R2*cv*(iv-iv2) / sqrt(R2*cv*(1-R2*cv*iv*(iv-thd)+(1-R2*cv*iv2*(iv2-thd))))) )
  
  #auc
  aova1=eval(f)
  
  #variance
  #aova2=eval(D(f,'r1'))^2*(1-R2)^2/nv
  aova2=eval(D(f,'R2'))^2*se^2
  
  ##R2 from AUC
  #qv=-qnorm(R2)
  #f2=expression( (2*qnorm(R2)^2 / ((iv2-iv)^2 + qnorm(R2)^2 * (iv*(iv-thd)+iv2*(iv2-thd))) ) /cv)
  
  ##R2 on the liability
  #R2l=eval(f2)
  ##variance
  #sel=sqrt( eval(D(f2,'R2'))^2*se^2)
  #z=list(auc=aova1,auc_se=sqrt(aova2),R2l=R2l,sel=sel)
  z=list(auc=aova1,auc_se=sqrt(aova2))
  return(z)
  
}

  