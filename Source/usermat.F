*deck,usermat      USERDISTRIB  parallel                                gal       
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, pVolDer, hrmflg, var3, var4,
     &                   var5, var6, var7)
c
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT"::usermat 
c
c*************************************************************************
c     *** primary function ***
c
c         The conventional theory of mechanism-based strain gradient plasticity model
c         If you use this code for research or industrial purposes, please cite:
c             V. Shlyannikov, E. Martinez-Paneda, A. Tumanov, A. Tartygasheva,
c             Crack tip fields and fracture resistance parameters based on strain gradient plasticity,
c             International Journal of Solids and Structures, Volumes 208–209, 2021, Pages 63-82      
c         https://doi.org/10.1016/j.ijsolstr.2020.10.015
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry or plane stress. 
c
c           A 3D material constitutive model present in usermat3d.F
c           Plane strain and axisymmetry cases can be found in usermatpd.F
c           For plane stress use usermatps.F
c           By default final library contained all this subroutines.       
c           
c           Be careful the state variable array is different in all subroutines!
c
c************************************************************************* 
c 
c
c*************************************************************************
c This part of comments is copied from ANSYS documentation
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c      hrmflg   (dp,sc,io)                flag to indicate harmonic analysis 
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      pVolDer  (dp,ar(3),o)              derivatives of volumetric potential wrt to J
c                                         pVolDer(1) = dU/dJ
c                                         pVolDer(2) = d^2U/dJ^2
c                                         pVolDer(3) = d^3U/dJ^3
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
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
     &                 dsdePl  (ncomp,ncomp),
     &                 pVolDer (3),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION hrmflg
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d
      EXTERNAL         usermat_harm

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      data             var1/0.0d0/
      data             var2/0.0d0/
c
c*************************************************************************
c

      IF(ncomp .eq. 6) THEN
c ***    3d
         call usermat3d (
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
         
      ELSE IF(ncomp .EQ. 4) THEN
c ***    plane strain and axisymmetric
         call usermatpd (
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

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
c ***    plane stress 
         call usermatps (
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

      ELSE IF(ncomp .EQ. 3) THEN
c ***    Following subroutines not realized here 
c        It's needed for correct ANSYS work
c
c ***    3d beam not realized
         call usermatbm (
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

      ELSE IF(ncomp .EQ. 1) THEN
c ***    1d beam not realized
         call usermat1d (
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

      END IF
      IF(nint(hrmflg)/=0)THEN
         call usermat_harm (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nProp,
     &                   Time,dTime,Temp,stress,dsdePl, tsstif,
     &                   Strain,prop,coords)
         return
       ENDIF
      return
      end      
*deck,usermat1d    USERDISTRIB  parallel                                gal
      subroutine usermat1d(
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
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT1D"::usermat1d       
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
c
      return
      end
*deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
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
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMATBM"::usermatbm 
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
     &                 sedEl,   sedPl,   epseq,   epsZZ,  cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
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
*deck,usermat_harm    USERDISTRIB  parallel                    jmgerken
      subroutine usermat_harm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nProp,
     &                   freq,dfreq,Temp,stress,jacobi,tsstif,
     &                   strain,prop,coords)
#include "impcom.inc"
c
c***************** Arguments *************************************
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nProp
      DOUBLE PRECISION freq,dfreq, Temp
      DOUBLE PRECISION 
     &                 stress  (ncomp,2),
     &                 jacobi  (ncomp,ncomp,2),
     &                 strain  (ncomp,2),
     &                 prop    (nProp),
     &                 coords  (3),
     &                 tsstif  (2)
c
      return
      end
