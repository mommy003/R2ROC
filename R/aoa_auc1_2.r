#' olkin_auc1_2 function
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm  
#' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
#' @param nv Sample size
#' @keywords source 
#' @return This function will be used as source code

  olkin_auc1_2 = function (omat,nv,kv) {

  thd=-qnorm(kv)
  zv=dnorm(thd)
  iv=zv/kv
  iv2 = -iv*kv/(1-kv)
  cv=kv*(1-kv)/zv^2

  r1=omat[2,1]
  r2=omat[3,1]

  #aova in p158 in Olkin and Finn

  f=expression(pnorm( r1^2*cv*(iv-iv2) / sqrt(r1^2*cv*(1-r1^2*cv*iv*(iv-thd)+(1-r1^2*cv*iv2*(iv2-thd))))) - pnorm( r2^2*cv*(iv-iv2) / sqrt(r2^2*cv*(1-r2^2*cv*iv*(iv-thd)+(1-r2^2*cv*iv2*(iv2-thd)))))  )


  av=array(0,3)
  av[1]=eval(D(f,'r1'))
  av[2]=eval(D(f,'r2'))
  av[3]=0

  ov=matrix(0,3,3)
  ov[1,1]=(1-omat[2,1]^2)^2/nv
  ov[2,2]=(1-omat[3,1]^2)^2/nv
  ov[3,3]=(1-omat[3,2]^2)^2/nv
  ov[2,1]=(0.5*(2*omat[3,2]-omat[2,1]*omat[3,1])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,2]^3)/nv
  ov[1,2]=ov[2,1]
  ov[3,1]=(0.5*(2*omat[3,1]-omat[2,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,1]^3)/nv
  ov[1,3]=ov[3,1]
  ov[3,2]=(0.5*(2*omat[2,1]-omat[3,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[2,1]^3)/nv
  ov[2,3]=ov[3,2]

  #auc difference
  aova1=eval(f)
  #variance of the difference
  aova2=t(av)%*%ov%*%(av)
  z=list(diff=aova1,var=aova2)
  return(z)

  }




