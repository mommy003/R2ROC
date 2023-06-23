
#' olkin_auc1 function
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm  
#' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
#' @param nv Sample size
#' @keywords source 
#' @return This function will be used as source code


  olkin_auc1 = function (omat,nv,kv) {

  thd=-qnorm(kv)
  zv=dnorm(thd)
  iv=zv/kv
  iv2 = -iv*kv/(1-kv)
  cv=kv*(1-kv)/zv^2

  #correlation matrix to estiamte AUC
  r1=omat[2,1]

  #aova in p158 in Olkin and Finn
  #f=expression(pnorm(sqrt ((r1^2*cv*(iv2-iv)^2 + r1^2*cv*iv2*(iv2-thd)) / (2-r1^2*cv*(iv*(iv-thd))))) )
  f=expression(pnorm( r1^2*cv*(iv-iv2) / sqrt(r1^2*cv*(1-r1^2*cv*iv*(iv-thd)+(1-r1^2*cv*iv2*(iv2-thd))))) )

  #auc
  aova1=eval(f)

  #variance
  aova2=eval(D(f,'r1'))^2*(1-omat[2,1]^2)^2/nv

  z=list(auc=aova1,var=aova2)
  return(z)

  }




