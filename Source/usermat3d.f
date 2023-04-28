*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
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
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT3D"::usermat3d 
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
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
c
c***************** User defined part *************************************
c ***     local variables
c
c     Young   (dp)        - Young modulus
c     nue     (dp)        - Poisson's Ratio
c
c ***     Additional variables      
c    i, j, k  (int)       - counters      
c
      INTEGER i, j, k
      
      DOUBLE PRECISION 
     &                 Young, nue
c
c***************** Body of the subroutine********************************   
c
      cutFactor=0
      dsdePl = 0.d0
c     get material properties
      Young=prop(1)
      nue=prop(2)
c     create elastic matrix            
      call elasticstiffnes(Young, nue, ncomp, dsdePl)

c     Update stresses and strains           
      Strain = Strain + dStrain    
      Stress = matmul(dsdePl, Strain)
      
c     calculate elastic work
      call calcenergy(Strain, Stress, ncomp, sedEl)	
      return
      end


      
      subroutine elasticstiffnes(Ex, nu, mcomp, outmatrix)
c *** subroutine return the elastic stiffnes matrix
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
          outmatrix=0.d0
          do icomp=1,3
              do jcomp=1,3
                  outmatrix(jcomp,icomp)=tmpelam
              end do
              outmatrix(icomp,icomp)=eg2+tmpelam
          end do
          do icomp=4, mcomp
              outmatrix(icomp,icomp)=eg2/2.d0
          end do      
          return
      end

      subroutine calcenergy(StrainTens, StressTens, numcomp, sed)
c *** subroutine return the elastic work
c
c *** Input arguments 
c         StrainTens   (dp,ar(ncomp),i)       - Strain at beginning of time increment
c         StressTens   (dp,ar(ncomp),i)       - Stress
c         numcomp                             - array dimension     
c                                       (4 - for plane strain, 6 - for full 3d)
c *** Output arguments 
c         sed  - elastic work      
c
      INTEGER numcomp
      DOUBLE PRECISION StrainTens(numcomp), StressTens(numcomp), sed
c
c***************** Body of the subroutine********************************   
c
	  sed = 0
	  do i = 1, ncomp
              sed = sed + StrainTens(i)*StressTens(i)
	  end do
          return
      end