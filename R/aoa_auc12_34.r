#' olkin_auc12_34 function
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
#' @param nv Sample size
#' @keywords source  
#' @return This function will be used as source code


  olkin_auc12_34 = function (omat,nv,kv) {

  thd=-qnorm(kv)
  zv=dnorm(thd)
  iv=zv/kv
  iv2 = -iv*kv/(1-kv)
  cv=kv*(1-kv)/zv^2

  #aova in p158 in Olkin and Finn, but this is for var(r2(0,1,2) - r2(0,3,4))
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
    )))

    - 

    pnorm((c44*( (c55/(c44*c55-c54^2))*c41 + (c54/(c54^2-
    c44*c55))*c51 )^2 +  2*c54*(( (c55/(c44*c55-c54^2))*
    c41 + (c54/(c54^2-c44*c55))*c51 ) *  ( (c54/(c54^2-c44*
    c55))*c41 + (c44/(c44*c55-c54^2))*c51 )) +  c55*( (c54/(c54^2-
    c44*c55))*c41 + (c44/(c44*c55-c54^2))*c51 )^2)*cv*(iv-iv2) /

   sqrt((c44*( (c55/(c44*c55-c54^2))*c41 + (c54/(c54^2-
    c44*c55))*c51 )^2 +  2*c54*(( (c55/(c44*c55-c54^2))*
    c41 + (c54/(c54^2-c44*c55))*c51 ) *  ( (c54/(c54^2-c44*
    c55))*c41 + (c44/(c44*c55-c54^2))*c51 )) +  c55*( (c54/(c54^2-
    c44*c55))*c41 + (c44/(c44*c55-c54^2))*c51 )^2)*cv*

    (1 - (c44*( (c55/(c44*c55-c54^2))*c41 + (c54/(c54^2-
    c44*c55))*c51 )^2 +  2*c54*(( (c55/(c44*c55-c54^2))*
    c41 + (c54/(c54^2-c44*c55))*c51 ) *  ( (c54/(c54^2-c44*
    c55))*c41 + (c44/(c44*c55-c54^2))*c51 )) +  c55*( (c54/(c54^2-
    c44*c55))*c41 + (c44/(c44*c55-c54^2))*c51 )^2)*cv*iv*(iv-thd) +

    (1 - (c44*( (c55/(c44*c55-c54^2))*c41 + (c54/(c54^2-
    c44*c55))*c51 )^2 +  2*c54*(( (c55/(c44*c55-c54^2))*
    c41 + (c54/(c54^2-c44*c55))*c51 ) *  ( (c54/(c54^2-c44*
    c55))*c41 + (c44/(c44*c55-c54^2))*c51 )) +  c55*( (c54/(c54^2-
    c44*c55))*c41 + (c44/(c44*c55-c54^2))*c51 )^2)*cv*iv2*(iv2-thd))
    ))) )


  c11=omat[1,1]
  c21=omat[2,1]
  c22=omat[2,2]
  c31=omat[3,1]
  c32=omat[3,2]
  c33=omat[3,3]
  c41=omat[4,1]
  c42=omat[4,2]
  c43=omat[4,3]
  c44=omat[4,4]
  c51=omat[5,1]
  c52=omat[5,2]
  c53=omat[5,3]
  c54=omat[5,4]
  c55=omat[5,5]

  av=array(0,6)
  av[1]=eval(D(f,'c21'))
  av[2]=eval(D(f,'c31'))
  av[3]=eval(D(f,'c41'))
  av[4]=eval(D(f,'c51'))
  av[5]=eval(D(f,'c32'))
  av[6]=eval(D(f,'c54'))


  ov=matrix(0,6,6)
  ov[1,1]=(1-omat[2,1]^2)^2/nv
  ov[2,2]=(1-omat[3,1]^2)^2/nv
  ov[3,3]=(1-omat[4,1]^2)^2/nv
  ov[4,4]=(1-omat[5,1]^2)^2/nv
  ov[5,5]=(1-omat[3,2]^2)^2/nv
  ov[6,6]=(1-omat[5,4]^2)^2/nv

  ov[2,1]=(0.5*(2*omat[3,2]-omat[2,1]*omat[3,1])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,2]^3)/nv
  ov[1,2]=ov[2,1]

  ov[3,1]=(0.5*(2*omat[4,2]-omat[2,1]*omat[4,1])*(1-omat[4,2]^2-omat[2,1]^2-omat[4,1]^2)+omat[4,2]^3)/nv
  ov[1,3]=ov[3,1]
  ov[3,2]=(0.5*(2*omat[4,3]-omat[3,1]*omat[4,1])*(1-omat[4,3]^2-omat[3,1]^2-omat[4,1]^2)+omat[4,3]^3)/nv
  ov[2,3]=ov[3,2]

  ov[4,1]=(0.5*(2*omat[5,2]-omat[2,1]*omat[5,1])*(1-omat[5,2]^2-omat[2,1]^2-omat[5,1]^2)+omat[5,2]^3)/nv
  ov[1,4]=ov[4,1]
  ov[4,2]=(0.5*(2*omat[5,3]-omat[3,1]*omat[5,1])*(1-omat[5,3]^2-omat[3,1]^2-omat[5,1]^2)+omat[5,3]^3)/nv
  ov[2,4]=ov[4,2]
  ov[4,3]=(0.5*(2*omat[5,4]-omat[4,1]*omat[5,1])*(1-omat[5,4]^2-omat[4,1]^2-omat[5,1]^2)+omat[5,4]^3)/nv
  ov[3,4]=ov[4,3]

  ov[5,1]=(0.5*(2*omat[3,1]-omat[2,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,1]^3)/nv
  ov[1,5]=ov[5,1]
  ov[5,2]=(0.5*(2*omat[2,1]-omat[3,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[2,1]^3)/nv
  ov[2,5]=ov[5,2]

  ov[5,3]=omat[2,1]^2+omat[3,1]^2+omat[4,2]^2+omat[4,3]^2
  ov[5,3]=ov[5,3]*0.5*omat[4,1]*omat[3,2]
  ov[5,3]=ov[5,3]+omat[2,1]*omat[4,3]+omat[3,1]*omat[4,2]
  ov[5,3]=ov[5,3]-omat[4,1]*omat[2,1]*omat[3,1]
  ov[5,3]=ov[5,3]-omat[4,1]*omat[4,2]*omat[4,3]
  ov[5,3]=ov[5,3]-omat[2,1]*omat[4,2]*omat[3,2]
  ov[5,3]=ov[5,3]-omat[3,1]*omat[4,3]*omat[3,2]
  ov[5,3]=ov[5,3]/nv 
  ov[3,5]=ov[5,3] 

  ov[5,4]=omat[2,1]^2+omat[3,1]^2+omat[5,2]^2+omat[5,3]^2
  ov[5,4]=ov[5,4]*0.5*omat[5,1]*omat[3,2]
  ov[5,4]=ov[5,4]+omat[2,1]*omat[5,3]+omat[3,1]*omat[5,2]
  ov[5,4]=ov[5,4]-omat[5,1]*omat[2,1]*omat[3,1]
  ov[5,4]=ov[5,4]-omat[5,1]*omat[5,2]*omat[5,3]
  ov[5,4]=ov[5,4]-omat[2,1]*omat[5,2]*omat[3,2]
  ov[5,4]=ov[5,4]-omat[3,1]*omat[5,3]*omat[3,2]
  ov[5,4]=ov[5,4]/nv 
  ov[4,5]=ov[5,4] 

  ov[6,1]=omat[4,1]^2+omat[5,1]^2+omat[4,2]^2+omat[5,2]^2
  ov[6,1]=ov[6,1]*0.5*omat[2,1]*omat[5,4]
  ov[6,1]=ov[6,1]+omat[4,1]*omat[5,2]+omat[4,2]*omat[5,1]
  ov[6,1]=ov[6,1]-omat[2,1]*omat[4,1]*omat[5,1]
  ov[6,1]=ov[6,1]-omat[2,1]*omat[4,2]*omat[5,2]
  ov[6,1]=ov[6,1]-omat[4,1]*omat[4,2]*omat[5,4]
  ov[6,1]=ov[6,1]-omat[5,1]*omat[5,2]*omat[5,4]
  ov[6,1]=ov[6,1]/nv 
  ov[1,6]=ov[6,1] 

  ov[6,2]=omat[4,1]^2+omat[5,1]^2+omat[4,3]^2+omat[5,3]^2
  ov[6,2]=ov[6,2]*0.5*omat[3,1]*omat[5,4]
  ov[6,2]=ov[6,2]+omat[4,1]*omat[5,3]+omat[4,3]*omat[5,1]
  ov[6,2]=ov[6,2]-omat[3,1]*omat[4,1]*omat[5,1]
  ov[6,2]=ov[6,2]-omat[3,1]*omat[4,3]*omat[5,3]
  ov[6,2]=ov[6,2]-omat[4,1]*omat[4,3]*omat[5,4]
  ov[6,2]=ov[6,2]-omat[5,1]*omat[5,3]*omat[5,4]
  ov[6,2]=ov[6,2]/nv 
  ov[2,6]=ov[6,2] 

  ov[6,3]=(0.5*(2*omat[5,1]-omat[4,1]*omat[5,4])*(1-omat[5,1]^2-omat[5,4]^2-omat[4,1]^2)+omat[5,1]^3)/nv
  ov[3,6]=ov[6,3]
  ov[6,4]=(0.5*(2*omat[4,1]-omat[5,4]*omat[5,1])*(1-omat[4,1]^2-omat[5,4]^2-omat[5,1]^2)+omat[4,1]^3)/nv
  ov[4,6]=ov[6,4]

  ov[6,5]=omat[4,2]^2+omat[5,2]^2+omat[4,3]^2+omat[5,3]^2
  ov[6,5]=ov[6,5]*0.5*omat[3,2]*omat[5,4]
  ov[6,5]=ov[6,5]+omat[4,2]*omat[5,3]+omat[4,3]*omat[5,2]
  ov[6,5]=ov[6,5]-omat[3,2]*omat[4,2]*omat[5,2]
  ov[6,5]=ov[6,5]-omat[3,2]*omat[4,3]*omat[5,3]
  ov[6,5]=ov[6,5]-omat[4,2]*omat[4,3]*omat[5,4]
  ov[6,5]=ov[6,5]-omat[5,2]*omat[5,3]*omat[5,4]
  ov[6,5]=ov[6,5]/nv 
  ov[5,6]=ov[6,5] 

  #auc difference
  aova1=eval(f)
  #variance of the difference
  aova2=t(av)%*%ov%*%(av)
  z=list(diff=aova1,var=aova2)
  return(z)

  }

