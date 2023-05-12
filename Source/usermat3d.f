*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
     &                   matId, noel,npt, layer, kspt,
     &                   jstep,kinc,keycut,
     &                   ndi,nshr,ntens,nstatv,nprop,
     &                   time,dTime,temp,dtemp,
     &                   stress,statev,ddsdde,sse,spd,epseq,
     &                   stran,dstran, epsPl, props, coords, 
     &                   var0, dfgrd0, dfgrd1,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c      
c     Custom user material subroutine for Conventional Gradient Plasticity Model.      
c     The procedure is created for SOLID186 element with Keyopt(2)=0.
c     SOLID186 have a 8 integration points. 
c     Coordinates of integration points in the coordinate system s-t-r = ± 0.577350269189626
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT3D"::usermat3d 
c
#include "impcom.inc"
c
      EXTERNAL         kdevia, keff,  ELPREV, REVERSE
      EXTERNAL         GET_ELMDATA, PUT_ELMDATA
              INTEGER 
     &                 matId, noel,
     &                 npt, Layer, jstep, kspt,
     &                 kdstep,kinc,keycut,
     &                 nDi,nShr,ntens,nStatv,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sse,   spd,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ntens  ), statev (nstatv),
     &                 ddsdde(ntens,ntens),
     &                 stran(ntens), dstran(ntens), 
     &                 epsPl(ntens), props(nprop), 
     &                 coords  (3),       
     &                 dfgrd0 (3,3), dfgrd1 (3,3),
     &                 tsstif  (2)
c
c     Integration Point Locations in local element CS
      DOUBLE PRECISION  s, t, r
c     Reserved by ANSYS variables
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8
c     state variables array for all integration points
      DOUBLE PRECISION SVAR(8,nstatv), SVARtmp(nstatv)
c     derivatives of shape functions at integration points
      DOUBLE PRECISION deriv(3,8)
c     auxiliary parameters for derivatives
      DOUBLE PRECISION P007, P057, P108
      PARAMETER       (P007       = 0.07216878364870325,
     &                 P057       = 0.577350269189626,
     &                 P108       = 1.d0/8.d0)
      
      DOUBLE PRECISION 
     1 ddsddt(ntens),drplde(ntens),
     2 predef(1),dpred(1),drot(3)

      parameter newton=1000, toler=1.0d-8
      DOUBLE PRECISION eelas(ntens)
c     matrices for Jacobian determination
      DOUBLE PRECISION xjacm(3,3), xjaci(3,3)
      DOUBLE PRECISION a1, a2, a3, a4, a5, a6, a7, a8,
     &  b1, b2, b3, b4, b5, b6, b7, b8,
     &  c1, c2, c3, c4, c5, c6, c7, c8, 
     &  eta(27)
      DOUBLE PRECISION xiden(3,3),dpstran(ntens),stra(3,3),
     1 destran(ntens),dstr(3,3),
     2 dpstrn(3,3),strad(3,3),dstre(6),
     3 dstra(3,3),dstrad(3,3),strain(3,3)
      
      DOUBLE PRECISION E, xnue, ebulk3, xk,eg2, eg, elam,
     1 syield, ele, sigmaf,sigmae,h, def, defi

      INTEGER k1, k2, i,kewton, k, kflag, j, prev,G(6)
      
      DOUBLE PRECISION  djacb, etat,
     1 rhs, dabs, deqpl, ep, q,  ene,
     2 xm ,b,gnd, ssd,td, sigi(ntens), sigElp(ntens)
      
c     temp parameters 
      DOUBLE PRECISION elast1, elast2
      
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0
     &                 )      
      
        DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
c     Some of the default names is changed here
c        
c     input arguments
c     ===============
c      cmname     (int,sc,i)               material #
c      noel       (int,sc,i)               element #
c      npt        (int,sc,i)               "k"th domain integration point
c      layer      (int,sc,i)               "k"th layer
c      kspty      (int,sc,i)               "k"th Section point
c      jstep      (int,sc,i)               load step number
c      kinc       (int,sc,i)               substep number
c      ndi        (int,sc,in)              # of direct components
c      nshr       (int,sc,in)              # of shear components
c      ntens      (int,sc,in)              nDirect + nShear
c      nstatv    (int,sc,i)               Number of state variables
c      nprops     (int,sc,i)               Number of material constants
c
c      temp       (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp      (dp,sc,in)               temperature increment 
c      time       (dp,sc,in)               time at beginning of increment (t)
c      dTime      (dp,sc,in)               current time increment (dt)
c
c      stran      (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStran     (dp,ar(ncomp),i)          Strain increment
c      props      (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords     (dp,ar(3),i)              current coordinates
c      dfgrd0     (dp,ar(3,3),i)            Deformation gradient at time t
c      dfgrd1     (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress     (dp,ar(ncomp),io)         stress
c      statev     (dp,ar(nstatev),io)      user state variables
c      sse        (dp,sc,io)                elastic work
c      spd        (dp,sc,io)                plastic work
c      epseq      (dp,sc,io)                equivalent plastic strain
c      epsPl      (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      ddsdde   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
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
      
      
c      ncomp   6   for 3D  (nshear=3)      
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D      
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |    
c
c     Structure of the ustatev() array
c     1-3 Integration point coordinates in global system
c     4-9 derivatives x,y,z,xy,yz,zx
c     10  gradient
c
c     obtain the gradient values from all integration points
      SVAR=0.d0
      do i=1,8
         call get_ElmData ('SVAR', noel,i, nstatv, SVARtmp)
          do j=1,nstatv
              SVAR(i,j)=SVARtmp(j)
         enddo
      enddo
c     elastic strains
      eelas=stran - epsPl
c *** Material properties      
c     Young modulus
      E=props(1)
c     Poisson ratio
      xnue=props(2)
c     Yeld stress
      syield=props(3)
c     characteristic length
      ele=props(4)
c     strain hardening exponent (0 < ene < 1)
      ene=props(5)
c     flag, statistically conserved dislocations. See Arsenlis and Parks (1998)
c     in most cases is eq to 1
      kflag=props(6)
       
      ebulk3=E/(1.d0-2.d0*xnue)
c     Bulk modulus
      xk=ebulk3/3.d0     
      eg2=E/(1.d0+xnue)
      eg=eg2/2.d0
      elam=(ebulk3-eg2)/3.d0 

c *** calculate elastic stiffness matrix (3d)
      ddsdde=0.d0   
       do k1=1,3
        do k2=1,3
         ddsdde(k2,k1)=elam
        enddo
        ddsdde(k1,k1)=eg2+elam
       enddo
       ddsdde(4,4)=eg 
       ddsdde(5,5)=eg 
       ddsdde(6,6)=eg 
c     it's needed for calculation of hourglass stiffness
      tsstif(1)=eg 
      
      if (all(dstran .eq. 0)) then 
           goto 100
          else
c *** gradients for all integration points
              
c     integration point coordinates in the isoparametric space              
          if (npt.eq.1) then   
             do k1=1,8  
                if (k1==1) then
                    s=-0.577350269189626
                    t=-0.577350269189626
                    r=-0.577350269189626
                elseif (k1==2) then 
                    s=0.577350269189626
                    t=-0.577350269189626
                    r=-0.577350269189626
                elseif (k1==3) then
                    s=0.577350269189626
                    t=0.577350269189626
                    r=-0.577350269189626
                elseif (k1==4) then      
                    s=-0.577350269189626
                    t=0.577350269189626
                    r=-0.577350269189626
                elseif (k1==5) then      
                    s=-0.577350269189626
                    t=-0.577350269189626 
                    r=0.577350269189626
                elseif (k1==6) then      
                    s=0.577350269189626
                    t=-0.577350269189626
                    r=0.577350269189626
                elseif (k1==7) then      
                    s=0.577350269189626
                    t=0.577350269189626
                    r=0.577350269189626
                elseif (k1==8) then      
                    s=-0.577350269189626
                    t=0.577350269189626
                    r=0.577350269189626
                end if
c     adopted linear shape functions     
         deriv(1,1)=-P007*(P057*r - 1)*(P057*t - 1)
         deriv(2,1)=-P057*(P057*r - 1)*(P007*s - P108)
         deriv(3,1)=-P057*(P057*t - 1)*(P007*s - P108)     
              
         deriv(1,2)=P007*(P057*r - 1)*(P057*t - 1)
         deriv(2,2)=P057*(P057*r - 1)*(P007*s + P108)
         deriv(3,2)=P057*(P057*t - 1)*(P007*s + P108)   
         
         deriv(1,3)=-P007*(P057*r - 1)*(P057*t + 1)
         deriv(2,3)=-P057*(P057*r - 1)*(P007*s + P108)
         deriv(3,3)=-P057*(P057*t + 1)*(P007*s + P108)  
         
         deriv(1,4)=P007*(P057*r - 1)*(P057*t + 1)
         deriv(2,4)=P057*(P057*r - 1)*(P007*s - P108)
         deriv(3,4)=P057*(P057*t + 1)*(P007*s - P108) 
         
         deriv(1,5)=P007*(P057*r + 1)*(P057*t - 1)
         deriv(2,5)=P057*(P057*r + 1)*(P007*s - P108)
         deriv(3,5)=P057*(P057*t - 1)*(P007*s - P108)
         
         deriv(1,6)=-P007*(P057*r + 1)*(P057*t - 1)
         deriv(2,6)=-P057*(P057*r + 1)*(P007*s + P108)
         deriv(3,6)=-P057*(P057*t - 1)*(P007*s + P108)
         
         deriv(1,7)=P007*(P057*r + 1)*(P057*t + 1)
         deriv(2,7)=P057*(P057*r + 1)*(P007*s + P108)
         deriv(3,7)=P057*(P057*t + 1)*(P007*s + P108)
         
         deriv(1,8)=-P007*(P057*r + 1)*(P057*t + 1)
         deriv(2,8)=-P057*(P057*r + 1)*(P007*s - P108)
         deriv(3,8)=-P057*(P057*t + 1)*(P007*s - P108)
          
c     xjacm array structure
c     columns are global x y z coordinates
c     lines are s t r
         xjacm(1,1)= deriv(1,1)*SVAR(1,1)+deriv(1,2)*SVAR(2,1)
     1 +deriv(1,3)*SVAR(3,1)+deriv(1,4)*SVAR(4,1) 
     2 +deriv(1,5)*SVAR(5,1)+deriv(1,6)*SVAR(6,1) 
     3 +deriv(1,7)*SVAR(7,1)+deriv(1,8)*SVAR(8,1) 
         
           xjacm(1,2)= deriv(2,1)*SVAR(1,1)+deriv(2,2)*SVAR(2,1)
     1 +deriv(2,3)*SVAR(3,1)+deriv(2,4)*SVAR(4,1) 
     2 +deriv(2,5)*SVAR(5,1)+deriv(2,6)*SVAR(6,1) 
     3 +deriv(2,7)*SVAR(7,1)+deriv(2,8)*SVAR(8,1) 
           
           xjacm(1,3)= deriv(3,1)*SVAR(1,1)+deriv(3,2)*SVAR(2,1)
     1 +deriv(3,3)*SVAR(3,1)+deriv(3,4)*SVAR(4,1) 
     2 +deriv(3,5)*SVAR(5,1)+deriv(3,6)*SVAR(6,1) 
     3 +deriv(3,7)*SVAR(7,1)+deriv(3,8)*SVAR(8,1) 
           
          xjacm(2,1)= deriv(1,1)*SVAR(1,2)+deriv(1,2)*SVAR(2,2)
     1 +deriv(1,3)*SVAR(3,2)+deriv(1,4)*SVAR(4,2) 
     2 +deriv(1,5)*SVAR(5,2)+deriv(1,6)*SVAR(6,2) 
     3 +deriv(1,7)*SVAR(7,2)+deriv(1,8)*SVAR(8,2) 
          
          xjacm(2,2)= deriv(2,1)*SVAR(1,2)+deriv(2,2)*SVAR(2,2)
     1 +deriv(2,3)*SVAR(3,2)+deriv(2,4)*SVAR(4,2) 
     2 +deriv(2,5)*SVAR(5,2)+deriv(2,6)*SVAR(6,2) 
     3 +deriv(2,7)*SVAR(7,2)+deriv(2,8)*SVAR(8,2) 
          
          xjacm(2,3)= deriv(3,1)*SVAR(1,2)+deriv(3,2)*SVAR(2,2)
     1 +deriv(3,3)*SVAR(3,2)+deriv(3,4)*SVAR(4,2) 
     2 +deriv(3,5)*SVAR(5,2)+deriv(3,6)*SVAR(6,2) 
     3 +deriv(3,7)*SVAR(7,2)+deriv(3,8)*SVAR(8,2) 
          
          
          xjacm(3,1)= deriv(1,1)*SVAR(1,3)+deriv(1,2)*SVAR(2,3)
     1 +deriv(1,3)*SVAR(3,3)+deriv(1,4)*SVAR(4,3) 
     2 +deriv(1,5)*SVAR(5,3)+deriv(1,6)*SVAR(6,3) 
     3 +deriv(1,7)*SVAR(7,3)+deriv(1,8)*SVAR(8,3) 
          
          xjacm(3,2)= deriv(2,1)*SVAR(1,3)+deriv(2,2)*SVAR(2,3)
     1 +deriv(2,3)*SVAR(3,3)+deriv(2,4)*SVAR(4,3) 
     2 +deriv(2,5)*SVAR(5,3)+deriv(2,6)*SVAR(6,3) 
     3 +deriv(2,7)*SVAR(7,3)+deriv(2,8)*SVAR(8,3) 
          
          xjacm(3,3)= deriv(3,1)*SVAR(1,3)+deriv(3,2)*SVAR(2,3)
     1 +deriv(3,3)*SVAR(3,3)+deriv(3,4)*SVAR(4,3) 
     2 +deriv(3,5)*SVAR(5,3)+deriv(3,6)*SVAR(6,3) 
     3 +deriv(3,7)*SVAR(7,3)+deriv(3,8)*SVAR(8,3) 
          
        
      if (all(xjacm .eq. 0)) then 
          xjaci=0
      else
          CALL REVERSE(3, xjacm, xjaci) 
      end if
      
      !dN/dx
         a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1)
     1                            + xjaci(1,3)*deriv(3,1) 
         a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2)
     1                            + xjaci(1,3)*deriv(3,2) 
         a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3)
     1                            + xjaci(1,3)*deriv(3,3) 
         a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4)
     1                            + xjaci(1,3)*deriv(3,4) 
         a5=xjaci(1,1)*deriv(1,5)+xjaci(1,2)*deriv(2,5)
     1                            + xjaci(1,3)*deriv(3,5)     
         a6=xjaci(1,1)*deriv(1,6)+xjaci(1,2)*deriv(2,6)
     1                            + xjaci(1,3)*deriv(3,6) 
         a7=xjaci(1,1)*deriv(1,7)+xjaci(1,2)*deriv(2,7)
     1                            + xjaci(1,3)*deriv(3,7)          
         a8=xjaci(1,1)*deriv(1,8)+xjaci(1,2)*deriv(2,8)
     1                            + xjaci(1,3)*deriv(3,8)          
      
       !dN/dy
         b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1)
     1                            + xjaci(2,3)*deriv(3,1) 
         b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2)
     1                            + xjaci(2,3)*deriv(3,2) 
         b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3)
     1                            + xjaci(2,3)*deriv(3,3) 
         b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4)
     1                            + xjaci(2,3)*deriv(3,4) 
         b5=xjaci(2,1)*deriv(1,5)+xjaci(2,2)*deriv(2,5)
     1                            + xjaci(2,3)*deriv(3,5)     
         b6=xjaci(2,1)*deriv(1,6)+xjaci(2,2)*deriv(2,6)
     1                            + xjaci(2,3)*deriv(3,6) 
         b7=xjaci(2,1)*deriv(1,7)+xjaci(2,2)*deriv(2,7)
     1                            + xjaci(2,3)*deriv(3,7)          
         b8=xjaci(2,1)*deriv(1,8)+xjaci(2,2)*deriv(2,8)
     1                            + xjaci(2,3)*deriv(3,8)     
         !dN/dz
         c1=xjaci(3,1)*deriv(1,1)+xjaci(3,2)*deriv(2,1)
     1                            + xjaci(3,3)*deriv(3,1) 
         c2=xjaci(3,1)*deriv(1,2)+xjaci(3,2)*deriv(2,2)
     1                            + xjaci(3,3)*deriv(3,2) 
         c3=xjaci(3,1)*deriv(1,3)+xjaci(3,2)*deriv(2,3)
     1                            + xjaci(3,3)*deriv(3,3) 
         c4=xjaci(3,1)*deriv(1,4)+xjaci(3,2)*deriv(2,4)
     1                            + xjaci(3,3)*deriv(3,4) 
         c5=xjaci(3,1)*deriv(1,5)+xjaci(3,2)*deriv(2,5)
     1                            + xjaci(3,3)*deriv(3,5)     
         c6=xjaci(3,1)*deriv(1,6)+xjaci(3,2)*deriv(2,6)
     1                            + xjaci(3,3)*deriv(3,6) 
         c7=xjaci(3,1)*deriv(1,7)+xjaci(3,2)*deriv(2,7)
     1                            + xjaci(3,3)*deriv(3,7)          
         c8=xjaci(3,1)*deriv(1,8)+xjaci(3,2)*deriv(2,8)
     1                            + xjaci(3,3)*deriv(3,8)  
         
         !dn111
         eta(1)=a1*SVAR(1,4) + a2*SVAR(2,4) + a3*SVAR(3,4)
     1        + a4*SVAR(4,4) + a5*SVAR(5,4) + a6*SVAR(6,4) 
     2        + a7*SVAR(7,4) + a8*SVAR(8,4) 
         !dn112
         eta(2)=2.d0*(a1*SVAR(1,7)+a2*SVAR(2,7)
     1    + a3*SVAR(3,7) + a4*SVAR(4,7)
     2    + a5*SVAR(5,7) + a6*SVAR(6,7)
     3    + a7*SVAR(7,7) + a8*SVAR(8,7))
     4    - b1*SVAR(1,4) - b2*SVAR(2,4) 
     5    - b3*SVAR(3,4) - b4*SVAR(4,4)
     6    - b5*SVAR(5,4) - b6*SVAR(6,4) 
     7    - b7*SVAR(7,4) - b8*SVAR(8,4)
         !dn113
         eta(3)=2.d0*(a1*SVAR(1,9)+a2*SVAR(2,9)
     1    + a3*SVAR(3,9) + a4*SVAR(4,9)
     2    + a5*SVAR(5,9) + a6*SVAR(6,9)
     3    + a7*SVAR(7,9) + a8*SVAR(8,9))
     4    - c1*SVAR(1,4) - c2*SVAR(2,4) 
     5    - c3*SVAR(3,4) - c4*SVAR(4,4)
     6    - c5*SVAR(5,4) - c6*SVAR(6,4) 
     7    - c7*SVAR(7,4) - c8*SVAR(8,4)
          !dn121
         eta(4)=b1*SVAR(1,4) + b2*SVAR(2,4) + b3*SVAR(3,4)
     1        + b4*SVAR(4,4) + b5*SVAR(5,4) + b6*SVAR(6,4) 
     2        + b7*SVAR(7,4) + b8*SVAR(8,4) 
         !dn122 
         eta(5)=a1*SVAR(1,5) + a2*SVAR(2,5) + a3*SVAR(3,5)
     1        + a4*SVAR(4,5) + a5*SVAR(5,5) + a6*SVAR(6,5) 
     2        + a7*SVAR(7,5) + a8*SVAR(8,5) 
         !dn123
         eta(6)=b1*SVAR(1,9) + b2*SVAR(2,9) + b3*SVAR(3,9)
     1        + b4*SVAR(4,9) + b5*SVAR(5,9) + b6*SVAR(6,9) 
     2        + b7*SVAR(7,9) + b8*SVAR(8,9) 
     3        + a1*SVAR(1,8) + a2*SVAR(2,8) + a3*SVAR(3,8)
     4        + a4*SVAR(4,8) + a5*SVAR(5,8) + a6*SVAR(6,8) 
     5        + a7*SVAR(7,8) + a8*SVAR(8,8) 
     4        - c1*SVAR(1,7) - c2*SVAR(2,7) 
     5        - c3*SVAR(3,7) - c4*SVAR(4,7)
     6        - c5*SVAR(5,7) - c6*SVAR(6,7) 
     7        - c7*SVAR(7,7) - c8*SVAR(8,7)
         !dn131
         eta(7)=c1*SVAR(1,4) + c2*SVAR(2,4) + c3*SVAR(3,4)
     1        + c4*SVAR(4,4) + c5*SVAR(5,4) + c6*SVAR(6,4) 
     2        + c7*SVAR(7,4) + c8*SVAR(8,4) 
         !dn132
         eta(8)=c1*SVAR(1,7) + c2*SVAR(2,7) + c3*SVAR(3,7)
     &        + c4*SVAR(4,7) + c5*SVAR(5,7) + c6*SVAR(6,7) 
     &        + c7*SVAR(7,7) + c8*SVAR(8,7)
     &        + a1*SVAR(1,8) + a2*SVAR(2,8)
     &        + a3*SVAR(3,8) + a4*SVAR(4,8)
     &        + a5*SVAR(5,8) + a6*SVAR(6,8)
     &        + a7*SVAR(7,8) + a8*SVAR(8,8)
     &        - b1*SVAR(1,9) - b2*SVAR(2,9) 
     &        - b3*SVAR(3,9) - b4*SVAR(4,9)
     &        - b5*SVAR(5,9) - b6*SVAR(6,9) 
     &        - b7*SVAR(7,9) - b8*SVAR(8,9)    
         !dn133
         eta(9)=a1*SVAR(1,8) + a2*SVAR(2,8)
     &        + a3*SVAR(3,8) + a4*SVAR(4,8)
     &        + a5*SVAR(5,8) + a6*SVAR(6,8)
         !dn211
        eta(10)=b1*SVAR(1,4) + b2*SVAR(2,4) + b3*SVAR(3,4)
     1        + b4*SVAR(4,4) + b5*SVAR(5,4) + b6*SVAR(6,4) 
     2        + b7*SVAR(7,4) + b8*SVAR(8,4) 
        !dn212 
        eta(11)=a1*SVAR(1,5) + a2*SVAR(2,5) + a3*SVAR(3,5)
     &        + a4*SVAR(4,5) + a5*SVAR(5,5) + a6*SVAR(6,5) 
     &        + a7*SVAR(7,5) + a8*SVAR(8,5)
        !dn213
        eta(12)=a1*SVAR(1,8) + a2*SVAR(2,8)
     &        + a3*SVAR(3,8) + a4*SVAR(4,8)
     &        + a5*SVAR(5,8) + a6*SVAR(6,8)
     &        + a7*SVAR(7,8) + a8*SVAR(8,8)
     &        + b1*SVAR(1,9) + b2*SVAR(2,9) + b3*SVAR(3,9)
     &        + b4*SVAR(4,9) + b5*SVAR(5,9) + b6*SVAR(6,9) 
     &        + b7*SVAR(7,9) + b8*SVAR(8,9) 
     &        - c1*SVAR(1,7) - c2*SVAR(2,7) 
     &        - c3*SVAR(3,7) - c4*SVAR(4,7)
     &        - c5*SVAR(5,7) - c6*SVAR(6,7) 
     &        - c7*SVAR(7,7) - c8*SVAR(8,7)       
        !dn221
        eta(13)=2.d0*(b1*SVAR(1,7) + b2*SVAR(2,7) + b3*SVAR(3,7)
     &        + b4*SVAR(4,7) + b5*SVAR(5,7) + b6*SVAR(6,7) 
     &        + b7*SVAR(7,7) + b8*SVAR(8,7)) 
     &        - a1*SVAR(1,5) - a2*SVAR(2,5) - a3*SVAR(3,5)
     &        - a4*SVAR(4,5) - a5*SVAR(5,5) - a6*SVAR(6,5) 
     &        - a7*SVAR(7,5) - a8*SVAR(8,5)
        !dn222
        eta(14)=b1*SVAR(1,5) + b2*SVAR(2,5) + b3*SVAR(3,5)
     1        + b4*SVAR(4,5) + b5*SVAR(5,5) + b6*SVAR(6,5) 
     2        + b7*SVAR(7,5) + b8*SVAR(8,5) 
        !dn223
        eta(15)=2.d0*(b1*SVAR(1,8) + b2*SVAR(2,8) + b3*SVAR(3,8)
     &        + b4*SVAR(4,8) + b5*SVAR(5,8) + b6*SVAR(6,8) 
     &        + b7*SVAR(7,8) + b8*SVAR(8,8))
     &        - c1*SVAR(1,5) - c2*SVAR(2,5) 
     &        - c3*SVAR(3,5) - c4*SVAR(4,5)
     &        - c5*SVAR(5,5) - c6*SVAR(6,5) 
     &        - c7*SVAR(7,5) - c8*SVAR(8,5) 
        !dn231
        eta(16)=c1*SVAR(1,7) + c2*SVAR(2,7) + c3*SVAR(3,7)
     &        + c4*SVAR(4,7) + c5*SVAR(5,7) + c6*SVAR(6,7) 
     &        + c7*SVAR(7,7) + c8*SVAR(8,7)
     &        + b1*SVAR(1,9) + b2*SVAR(2,9) + b3*SVAR(3,9)
     &        + b4*SVAR(4,9) + b5*SVAR(5,9) + b6*SVAR(6,9) 
     &        + b7*SVAR(7,9) + b8*SVAR(8,9)  
     &        - a1*SVAR(1,8) - a2*SVAR(2,8) - a3*SVAR(3,8)
     &        - a4*SVAR(4,8) - a5*SVAR(5,8) - a6*SVAR(6,8) 
     &        - a7*SVAR(7,8) - a8*SVAR(8,8)   
        !dn232
        eta(17)=c1*SVAR(1,5) + c2*SVAR(2,5) 
     &        + c3*SVAR(3,5) + c4*SVAR(4,5)
     &        + c5*SVAR(5,5) + c6*SVAR(6,5) 
     &        + c7*SVAR(7,5) + c8*SVAR(8,5)
        !dn233
        eta(18)=b1*SVAR(1,6) + b2*SVAR(2,6) + b3*SVAR(3,6)
     1        + b4*SVAR(4,6) + b5*SVAR(5,6) + b6*SVAR(6,6) 
     2        + b7*SVAR(7,6) + b8*SVAR(8,6)
        !dn311
        eta(19)=c1*SVAR(1,4) + c2*SVAR(2,4) + c3*SVAR(3,4)
     &        + c4*SVAR(4,4) + c5*SVAR(5,4) + c6*SVAR(6,4) 
     &        + c7*SVAR(7,4) + c8*SVAR(8,4)
        !dn312
        eta(20)=a1*SVAR(1,8) + a2*SVAR(2,8)
     &        + a3*SVAR(3,8) + a4*SVAR(4,8)
     &        + a5*SVAR(5,8) + a6*SVAR(6,8)
     &        + a7*SVAR(7,8) + a8*SVAR(8,8)
     &        + c1*SVAR(1,7) + c2*SVAR(2,7) 
     &        + c3*SVAR(3,7) + c4*SVAR(4,7)
     &        + c5*SVAR(5,7) + c6*SVAR(6,7) 
     &        + c7*SVAR(7,7) + c8*SVAR(8,7)    
     &        - b1*SVAR(1,9) - b2*SVAR(2,9) 
     &        - b3*SVAR(3,9) - b4*SVAR(4,9)
     &        - b5*SVAR(5,9) - b6*SVAR(6,9) 
     &        - b7*SVAR(7,9) - b8*SVAR(8,9)
        !dn313
        eta(21)=a1*SVAR(1,6) + a2*SVAR(2,6)
     &        + a3*SVAR(3,6) + a4*SVAR(4,6)
     &        + a5*SVAR(5,6) + a6*SVAR(6,6)
     &        + a7*SVAR(7,6) + a8*SVAR(8,6)
        !dn321
         eta(22)=b1*SVAR(1,9) + b2*SVAR(2,9) + b3*SVAR(3,9)
     &        + b4*SVAR(4,9) + b5*SVAR(5,9) + b6*SVAR(6,9) 
     &        + b7*SVAR(7,9) + b8*SVAR(8,9)
     &        + c1*SVAR(1,7) + c2*SVAR(2,7) + c3*SVAR(3,7)
     &        + c4*SVAR(4,7) + c5*SVAR(5,7) + c6*SVAR(6,7) 
     &        + c7*SVAR(7,7) + c8*SVAR(8,7)
     &        - a1*SVAR(1,8) - a2*SVAR(2,8) - a3*SVAR(3,8)
     &        - a4*SVAR(4,8) - a5*SVAR(5,8) - a6*SVAR(6,8) 
     &        - a7*SVAR(7,8) - a8*SVAR(8,8) 
        !dn322
         eta(23)=c1*SVAR(1,5) + c2*SVAR(2,5) 
     &        + c3*SVAR(3,5) + c4*SVAR(4,5)
     &        + c5*SVAR(5,5) + c6*SVAR(6,5) 
     &        + c7*SVAR(7,5) + c8*SVAR(8,5)
        !dn323
         eta(24)=b1*SVAR(1,6) + b2*SVAR(2,6) + b3*SVAR(3,6)
     1        + b4*SVAR(4,6) + b5*SVAR(5,6) + b6*SVAR(6,6) 
     2        + b7*SVAR(7,6) + b8*SVAR(8,6)
        !dn331
         eta(25)=2.d0*(c1*SVAR(1,9)+c2*SVAR(2,9)
     &        + c3*SVAR(3,9) + c4*SVAR(4,9)
     &        + c5*SVAR(5,9) + c6*SVAR(6,9)
     &        + c7*SVAR(7,9) + c8*SVAR(8,9))
     &        - a1*SVAR(1,6) - a2*SVAR(2,6) - a3*SVAR(3,6)
     &        - a4*SVAR(4,6) - a5*SVAR(5,6) - a6*SVAR(6,6) 
     &        - a7*SVAR(7,6) - a8*SVAR(8,6) 
         !dn332
         eta(26)=2.d0*(c1*SVAR(1,8)+c2*SVAR(2,8)
     &        + c3*SVAR(3,8) + c4*SVAR(4,8)
     &        + c5*SVAR(5,8) + c6*SVAR(6,8)
     &        + c7*SVAR(7,8) + c8*SVAR(8,8))
     &        - b1*SVAR(1,6) - b2*SVAR(2,6) 
     &        - b3*SVAR(3,6) - b4*SVAR(4,6)
     &        - b5*SVAR(5,6) - b6*SVAR(6,6) 
     &        - b7*SVAR(7,6) - b8*SVAR(8,6)
         !dn333
         eta(27)=c1*SVAR(1,6)+c2*SVAR(2,6)
     &        + c3*SVAR(3,6) + c4*SVAR(4,6)
     &        + c5*SVAR(5,6) + c6*SVAR(6,6)
     &        + c7*SVAR(7,6) + c8*SVAR(8,6)
         
         
        etat= eta(1)**2  + eta(2)**2  + eta(3)**2
     &      + eta(4)**2  + eta(5)**2  + eta(6)**2
     &      + eta(7)**2  + eta(8)**2  + eta(9)**2
     &      + eta(10)**2 + eta(11)**2 + eta(12)**2
     &      + eta(13)**2 + eta(14)**2 + eta(15)**2
     &      + eta(16)**2 + eta(17)**2 + eta(18)**2
     &      + eta(19)**2 + eta(20)**2 + eta(21)**2   
     &      + eta(22)**2 + eta(23)**2 + eta(24)**2
     &      + eta(25)**2 + eta(26)**2 + eta(27)**2
        
        !
        
        SVAR(k1,10) = SVAR(k1,10) + sqrt((1.d0/4.d0)*(etat))
    
        
          enddo
c     Save obtained gradients
          do k1=1,8
            do k2=1,nstatv
              SVARtmp(k2)=SVAR(k1,k2)
            enddo
           call put_ElmData ('SVAR', noel, k1, nstatv, SVARtmp)
          enddo
          
       endif
                            
       xiden=0.d0
       do i=1,3
        xiden(i,i)=1.d0   
       enddo    
       stra=0.d0   
c     Phisical strains
       do i=1,3
        stra(i,i)=eelas(i)
       enddo
       stra(1,2)=eelas(4)/2.d0
       stra(2,1)=eelas(4)/2.d0
       stra(1,3)=eelas(6)/2.d0
       stra(3,1)=eelas(6)/2.d0
       stra(2,3)=eelas(5)/2.d0
       stra(3,2)=eelas(5)/2.d0
c     deviatoric elastic strains 
       call kdevia(stra,xiden,strad)
       dstra=0.d0
c     elastic strains increment
       do i=1,3
        dstra(i,i)=dstran(i)
       enddo
       dstra(1,2)=dstran(4)/2.d0
       dstra(2,1)=dstran(4)/2.d0
       dstra(1,3)=dstran(6)/2.d0
       dstra(3,1)=dstran(6)/2.d0
       dstra(3,2)=dstran(5)/2.d0
       dstra(2,3)=dstran(5)/2.d0
c     deviatoric elastic strains increment
       call kdevia(dstra,xiden,dstrad)
       strain=strad+dstrad
c     equivalent strain
       call keff(strain,def)
c     equivalent strain increment
       call keff(dstrad,defi)
c *** Plastic      
c     flow stress
      sigmaf=syield*((E/syield)**ene)*sqrt((epseq+(syield/E))**(2*ene)
     1 + ele*SVAR(npt,10))
c          
       sigmae=syield
       h=0.d0
       do kewton=1, newton 
        rhs=3.d0*eg*(def-defi*((sigmae/sigmaf)**20))-sigmae
        sigmae=sigmae+rhs/(3.d0*eg*h+1.d0)
        h=20.d0*defi*((sigmae/sigmaf)**19)*(1.d0/sigmaf)
        if(dabs(rhs).lt.toler) goto 20
       end do
c       sigmae=0
c       write(7,19) newton
c   19 format(//,30x,'***warning - plasticity algorithm did not ',
c     1 'converged after',i3, 'iterations')
 20    continue   
          
       deqpl=def-sigmae/(3.d0*eg)
       epseq=epseq+deqpl
c
       dstr=strain*2.d0*eg/(1.d0+deqpl*3.d0*eg/sigmae)
c
       do i=1,3
        dstre(i)=dstr(i,i)
       enddo
       dstre(4)=dstr(1,2)
       dstre(5)=dstr(2,3)
       dstre(6)=dstr(3,1)
c
       dpstrn=dstr*(3.d0*deqpl)/(2.d0*sigmae)
c
       do i=1,3             
        dpstran(i)=dpstrn(i,i)
       enddo
       dpstran(4)=2.d0*dpstrn(1,2)
       dpstran(5)=2.d0*dpstrn(2,3)
       dpstran(6)=2.d0*dpstrn(3,1)
c       
       destran=dstran-dpstran
       epsPl=epsPl+dpstran
       eelas=eelas+destran
       ep=eelas(1)+eelas(2)+eelas(3)
       stran=eelas+epsPl
c      
       do j=1,3
        stress(j)=dstre(j)+xk*ep
       enddo
       stress(4)=dstre(4)
       stress(5)=dstre(5)
       stress(6)=dstre(6)
c       
       statev(4:9)= dpstran
c ***  material jacobian matrix  
       q=(2.d0/3.d0)*(sigmae/def)
       r=((h-(deqpl/sigmae))/(def*sigmae))*(3.d0*eg)/(1.d0+3.d0*eg*h)
       ddsdde=0
       do i=1,3
        do j=1,3
          ddsdde(i,j)=q*xiden(i,j)+(xk-q*1.d0/3.d0)-r*dstre(i)*dstre(j)
        end do
       end do
       do k=1,3
        ddsdde(k,4) = -r*dstre(k)*dstre(4)
        ddsdde(4,k) = ddsdde(k,4)
      end do
      do k=1,4
        ddsdde(k,5) = -r*dstre(k)*dstre(5)
        ddsdde(5,k) = ddsdde(k,5)
      end do
      do k=1,5
        ddsdde(k,6) = -r*dstre(k)*dstre(6)
        ddsdde(6,k) = ddsdde(k,6)
      end do
       ddsdde(4,4) = q/2.d0 - r*dstre(4)*dstre(4)
       ddsdde(5,5) = q/2.d0 - r*dstre(5)*dstre(5)
       ddsdde(6,6) = q/2.d0 - r*dstre(6)*dstre(6)

       
c *** OUTPUT
       kflag=props(6)
       if (kflag.eq.1) then
c     fcc
c     arsenlis and parks (1998)
              xm=3.06d0
              b=0.2555d-6
              r=1.9 
          else
c     bcc
              xm=2.9d0
              b=0.2725d-6
              r=1.9 
       endif
      
       if (ele.eq.0.d0) then
        gnd=0.d0
       else
        gnd=r*SVAR(npt,10)/b
       endif    
       ssd=((syield*(E/syield)**ene*(epseq+syield/E)**ene)/
     & (xm*0.5*eg*b))**2
       td=gnd+ssd
       
       
c     gradient 
       statev(10)= SVAR(npt,10) 
c     to plot ssd in m^-2 
       statev(11)=(1000000.0d0)*ssd
c     to plot gnd in m^-2 
       statev(12)=(1000.0d0)*gnd
c     total dislocations
       statev(13)=td 
  
      endif 
100   continue     
c      
      statev(1)=coords(1)
      statev(2)=coords(2)
      statev(3)=coords(3)
      return
      end
