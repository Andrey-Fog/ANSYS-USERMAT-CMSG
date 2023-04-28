*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5, 
     &                   var6, var7)
c
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMATPS"::usermatps
c
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(6),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
c
c***************** User defined part *************************************
c
      
      
      return
      end
       subroutine elasticstiffnes_ps(Ex, nu, mcomp, outmatrix)
c *** subroutine return the elastic stiffnes matrix for plane stress
c
c *** Input arguments 
c         Ex                          - Young's modulus
c         nu                          - Poisson's ratio
c         mcomp                       - array dimension     
c                                       (4 - for plane strain, 6 - for full 3d)
c *** Output arguments 
c         outmatrix(mcomp,mcomp)      - Output stiffnes matrix
c
c ***     local variables
c     eg2     (dp)        - Two shear modulus 
c
c ***     Additional variables      
c     tmpelam    (dp)        
c
      INTEGER icomp, jcomp, mcomp
      DOUBLE PRECISION 
     &                 Ex, nu, eg2, tmpelam
      DOUBLE PRECISION outmatrix(mcomp,mcomp)     
      
c
c***************** Body of the subroutine********************************   
c
          eg2=Ex/(1.d0+nu)
          tmpelam=(Ex/(1.d0-2.d0*nu)-eg2)/3.d0            

          outmatrix = 0
          do icomp=1,2
              do jcomp=1,2
                  outmatrix(jcomp,icomp)=tmpelam
              end do
              outmatrix(icomp,icomp)=eg2+tmpelam
          end do
          do icomp=3, mcomp
              outmatrix(icomp,icomp)=eg2/2.d0
          end do      
          return
      end
 


