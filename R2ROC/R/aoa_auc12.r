#' olkin_auc12 function
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
#' @param nv Sample size
#' @param kv Population prevalance
#' @keywords source  
#' @return This function will be used as source code


  olkin_auc12 = function (omat,nv,kv) {

  thd=-qnorm(kv)
  zv=dnorm(thd)
  iv=zv/kv
  iv2 = -iv*kv/(1-kv)
  cv=kv*(1-kv)/zv^2

  #correlation matrix to estiamte multiple AUC 

  f=expression( pnorm((c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 -
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) *
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 *
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 -
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2) *cv*(iv-iv2) /

   sqrt((c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 -
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) *
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 *
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 -
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2) *cv* 

    (1 - (c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 -
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) *
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 *
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 -
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2) *cv*iv*(iv-thd) +

    (1 - (c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 -
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) *
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 *
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 -
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2) *cv*iv2*(iv2-thd))
    ))))


  c11=omat[1,1]
  c21=omat[2,1]
  c22=omat[2,2]
  c31=omat[3,1]
  c32=omat[3,2]
  c33=omat[3,3]

  #auc
  aova1=eval(f)

  av=array(0,3)
  av[1]=eval(D(f,'c21'))
  av[2]=eval(D(f,'c31'))
  av[3]=eval(D(f,'c32'))

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

  #variance of the difference
  aova2=t(av)%*%ov%*%(av)

  z=list(auc=aova1,var=aova2)
  return(z)


  }




