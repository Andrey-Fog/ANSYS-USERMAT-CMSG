*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
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
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMATPS"::usermatps 
      use kplastic
#include "impcom.inc"   
      EXTERNAL         kdevia, keff, GET_ELMDATA, ELPREV ,EGEN
              INTEGER 
     &                 matId, noel,
     &                 npt, Layer, jstep, kspt,
     &                 kdstep,kinc,keycut,
     &                 nDi,nShr,ntens,nStatv,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sse,   spd,   epseq,   epsZZ, EGEN
      DOUBLE PRECISION 
     &                 stress  (ntens  ), statev (nstatv), 
     &                 ddsdde(ntens,ntens),
     &                 stran(ntens), dstran(ntens), 
     &                 epsPl   (ntens  ), props(nprop), 
     &                 coords  (3),       
     &                 dfgrd0 (3,3), dfgrd1 (3,3),
     &                 tsstif  (2), sigi(ntens)
c
     
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8, s, t
      
      
      DOUBLE PRECISION 
     1 ddsddt(ntens),drplde(ntens),
     2 predef(1),dpred(1),drot(3,3)
    
     

      parameter newton=1000, toler=1.0d-8

      DOUBLE PRECISION xiden(3,3),dpstran(4),stra(3,3),
     1 destran(4),dstr(3,3),strains(3,3),
     2 eelas(4),dpstrn(3,3),eplas(4),strad(3,3),dstre(4),xjaci(2,2),
     3 dstra(3,3),dstrad(3,3),strain(3,3),eta(14),deriv(2,4),xjacm(2,2)
      
      DOUBLE PRECISION E, xnue, ebulk3, xk,eg2, eg, elam,
     1 eqplas,syield, ele, sigmaf,sigmae,h, def, defi
     2 
      INTEGER k1, k2, i,kewton, k, kflag, j, prev
      
      DOUBLE PRECISION a1, a2, a3, a4, b1, b2, b3, b4, djacb, etat,
     1 rhs, dabs, deqpl, ep, q, r, ene,
     2 xm ,b,gnd, ssd,td
      
      DOUBLE PRECISION tmp,lamme
c
c *** Obtaining plastic strains from previous iteration
      eqplas=statev(1+2*ntens+2)
c      eplas =epsPl
      call get_ElmData ('EPTO', noel, npt, 4, eelas)
      call get_ElmData ('EPPL', noel, npt, 4, eplas)
      call get_ElmData ('SIG', noel,npt, 3, sigi)
       do i=1,4
          eelas(i) = eelas(i) - eplas (i)
       enddo 
c   
      E=props(1)
      xnue=props(2)
c      
c     it's needed for calculation of hourglass stiffness,
c     tsstif have no analog in ABAQUS
      eg2=E/(1.d0+xnue)
      eg=eg2/2.d0
      tsstif(1)=eg
      lamme=2*xnue*eg/(1-2*xnue)
c
      elam=E/(1.d0-xnue*xnue)
      ddsdde=0.d0 
      ddsdde(1,1)=elam
      ddsdde(1,2)=elam*xnue
      ddsdde(1,3)=0
      
      ddsdde(2,1)=elam*xnue
      ddsdde(2,2)=elam
      ddsdde(2,3)=0
      
      ddsdde(3,1)=0
      ddsdde(3,2)=0
      ddsdde(3,3)=elam*(1.d0-xnue)/2
      
       

  
      ebulk3=E/(1.d0-2.d0*xnue)
      xk=ebulk3/3.d0 !объемный модуль упругости
      
      
c     Rewrite  
      if (all(dstran .eq. 0)) then
          if (npt.eq.1) then
              strainp=0.d0    
              eqgrad=0.d0
          end if
c
      else
c     integration point coordinates in the isoparametric space
        if (npt.eq.1) then
        do k1=1,4
         if (k1==1) then
          s=-0.577350269189626
          t=-0.577350269189626
         elseif (k1==2) then 
          s=0.577350269189626
          t=-0.577350269189626
         elseif (k1==3) then
          s=0.577350269189626
          t=0.577350269189626  
         elseif (k1==4) then      
          s=-0.577350269189626
          t=0.577350269189626   
         end if
c     adopted linear shape functions 
         deriv(1,1)=0.083333333333333401314*t-0.1443375672974065 
         deriv(2,1)=0.083333333333333401314*s-0.1443375672974065
         
         deriv(1,2)=-0.083333333333333401314*t+0.1443375672974065
         deriv(2,2)=-0.083333333333333401314*s-0.1443375672974065
         
         
         deriv(1,3)=0.083333333333333401314*t+0.1443375672974065
         deriv(2,3)=0.083333333333333401314*s+0.1443375672974065
         
         
         deriv(1,4)=-0.083333333333333401314*t-0.1443375672974065
         deriv(2,4)=-0.083333333333333401314*s+0.1443375672974065

         
         xjacm(1,1)= deriv(1,1)*cordex(noel,1)+deriv(1,2)*cordex(noel,2)
     1 +deriv(1,3)*cordex(noel,3)+deriv(1,4)*cordex(noel,4)
    
         xjacm(1,2)= deriv(1,1)*cordey(noel,1)+deriv(1,2)*cordey(noel,2)
     1 +deriv(1,3)*cordey(noel,3)+deriv(1,4)*cordey(noel,4)
     
         xjacm(2,1)= deriv(2,1)*cordex(noel,1)+deriv(2,2)*cordex(noel,2)
     1 +deriv(2,3)*cordex(noel,3)+deriv(2,4)*cordex(noel,4)
      
         xjacm(2,2)= deriv(2,1)*cordey(noel,1)+deriv(2,2)*cordey(noel,2)
     1 +deriv(2,3)*cordey(noel,3)+deriv(2,4)*cordey(noel,4)

         djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1) 
      
         xjaci(1,1)=xjacm(2,2)/djacb 
         xjaci(1,2)=-xjacm(1,2)/djacb  
         xjaci(2,1)=-xjacm(2,1)/djacb   
         xjaci(2,2)=xjacm(1,1)/djacb

         a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1) 
         a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2) 
         a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3) 
         a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4) 
         b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1) 
         b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2)  
         b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3) 
         b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4) 
         !dn111
         eta(1)=a1*strainp(noel,1)+a2*strainp(noel,5)+a3*strainp(noel,9)
     1 +a4*strainp(noel,13)
         !dn112
         eta(2)=2.d0*(a1*strainp(noel,4)+a2*strainp(noel,8)
     1 +a3*strainp(noel,12)+a4*strainp(noel,16))
     2 -b1*strainp(noel,1)-b2*strainp(noel,5) 
     3 -b3*strainp(noel,9)-b4*strainp(noel,13)
         !dn121
         eta(3)=b1*strainp(noel,1)+b2*strainp(noel,5)+b3*strainp(noel,9)
     1 +b4*strainp(noel,13)
         !dn122    
        eta(4)=a1*strainp(noel,2)+a2*strainp(noel,6)+a3*strainp(noel,10)
     1 +a4*strainp(noel,14)
        !dn133   
        eta(5)=a1*strainp(noel,3)+a2*strainp(noel,7)+a3*strainp(noel,11)
     1 +a4*strainp(noel,15)
        !dn211 
         eta(6)=b1*strainp(noel,1)+b2*strainp(noel,5)+b3*strainp(noel,9)
     1 +b4*strainp(noel,13)
        !dn212
        eta(7)=a1*strainp(noel,2)+a2*strainp(noel,6)+a3*strainp(noel,10)
     1 +a4*strainp(noel,14)
        !dn221  
        eta(8)=2.d0*(b1*strainp(noel,4)+b2*strainp(noel,8)
     1 +b3*strainp(noel,12)+b4*strainp(noel,16))
     2 -a1*strainp(noel,2)-a2*strainp(noel,6)
     3 -a3*strainp(noel,10)-a4*strainp(noel,14)
        !dn222  
        eta(9)=b1*strainp(noel,2)+b2*strainp(noel,6)+b3*strainp(noel,10)
     1 +b4*strainp(noel,14)
        !dn233
         eta(10)=b1*strainp(noel,3)+b2*strainp(noel,7) 
     1 +b3*strainp(noel,11)+b4*strainp(noel,15)
        !dn313
         eta(11)=a1*strainp(noel,3)+a2*strainp(noel,7) 
     1 +a3*strainp(noel,11)+a4*strainp(noel,15)
        !dn323
         eta(12)=b1*strainp(noel,3) + b2*strainp(noel,7)
     1 +b3*strainp(noel,11)+b4*strainp(noel,15)
        !dn331
         eta(13)=-(a1*strainp(noel,3)+a2*strainp(noel,7) 
     1 +a3*strainp(noel,11)+a4*strainp(noel,15))
        !dn332 
         eta(14)=-(b1*strainp(noel,3)+b2*strainp(noel,7)
     1 +b3*strainp(noel,11)+b4*strainp(noel,15))

        etat=eta(1)**2+eta(2)**2+eta(3)**2+eta(4)**2+eta(5)**2+eta(6)**2
     1 +eta(7)**2+eta(8)**2+eta(9)**2+eta(10)**2+eta(11)**2+eta(12)**2
     2 +eta(13)**2+eta(14)**2

         eqgrad(noel,k1)=eqgrad(noel,k1)+sqrt((1.d0/4.d0)*(etat))
        enddo
       endif          
c              
       xiden=0.d0
       do i=1,3
        xiden(i,i)=1.d0   
       enddo    

       stra=0.d0
       
       strad=0.d0
       dstrad=0.d0         
       do i=1,3
          stra(i,i)=eelas(i)
       enddo
       stra(1,2)=eelas(4)/2.d0
       stra(2,1)=eelas(4)/2.d0
       stra(3,3)=-xnue/E*(sigi(1)+sigi(2))
       
c     Strain deviator
       call kdevia(stra,xiden,strad)
       dstra=0.d0
c     Strain increment deviator       
       do i=1,2
        dstra(i,i)=dstran(i)
       enddo
       dstra(1,2)=dstran(3)/2.d0
       dstra(2,1)=dstran(3)/2.d0  
c      new Z strain nd stress vector
c
       sigi(1)=elam*dstran(1)+xnue*elam*dstran(2)
       sigi(2)=xnue*elam*dstran(1)+elam*dstran(2)
       dstra(3,3)=-xnue/E*(sigi(1)+sigi(2))
       
       
c     Strain increment deviator
       call kdevia(dstra,xiden,dstrad)
       strain=strad+dstrad !ex,ey,ez,exy
c     Effective strains
       strains=stra+dstra
       
       call keff(strain,def)
       def=((strain(1,1))**2
     1 +(strain(2,2))**2
     2 +(strain(3,3))**2
     3 +2.d0*(strain(1,2))**2)
       def=sqrt(2.d0/3.d0*def)
c     Effective strains increment       
       call keff(dstrad,defi)
       defi=((dstrad(1,1))**2
     1 +(dstrad(2,2))**2
     2 +(dstrad(3,3))**2
     3 +2.d0*(dstrad(1,2))**2)
       defi=sqrt(2.d0/3.d0*defi)

c *** Plastic
       syield=props(3)
       ele=props(4)
       ene=props(5)
c     flow stress
      sigmaf=syield*((E/syield)**ene)*sqrt((eqplas+(syield/E))**(2*ene)
     1 + ele*eqgrad(noel,npt))
     
c          
       !eg=E/3
       sigmae=syield
       h=0.d0
       do kewton=1, newton 
        rhs=3.d0*eg*(def-defi*((sigmae/sigmaf)**20))-sigmae
        sigmae=sigmae+rhs/(3.d0*eg*h+1.d0)
        h=20.d0*defi*((sigmae/sigmaf)**19)*(1.d0/sigmaf)
        if(dabs(rhs).lt.toler) goto 20
       end do
!       write(7,19) newton
!   19 format(//,30x,'***warning - plasticity algorithm did not ',
!     1 'converge after',i3, 'iterations')
   20 continue

       deqpl=def-sigmae/(3.d0*eg) !3.15 
       eqplas=eqplas+deqpl
       
       
       dstr=strain*2.d0*eg/(1.d0+deqpl*3.d0*eg/sigmae) !3.10 
       dstre=0.d0
       do i=1,3
        dstre(i)=dstr(i,i)
       enddo
       dstre(4)=dstr(1,2)
        

       dpstrn=dstr*(3.d0*deqpl)/(2.d0*sigmae)     !3.5

       do i=1,3
        dpstran(i)=dpstrn(i,i)
       enddo
       dpstran(4)=2.d0*dpstrn(1,2)
       
       
       destran(1)=dstran(1)-dpstran(1)
       destran(2)=dstran(2)-dpstran(2) 
       destran(3)=dstra(3,3)-dpstran(3)
       destran(4)=dstran(3)-dpstran(4)
       
       eplas=eplas+dpstran
       eelas=eelas+destran
       ep=(eelas(1)+eelas(2)+eelas(3))
       
       stran(1)=eelas(1)+eplas(1)
       stran(2)=eelas(2)+eplas(2)
       stran(3)=eelas(4)+eplas(4)

       do j=1,3
        stress(j)=dstre(j)+xk*ep
       enddo
       stress(3)=dstre(4)

       if (npt==1) then
        strainp(noel,1:4)=dpstran
       elseif (npt==2) then
        strainp(noel,5:8)=dpstran 
       elseif (npt==3) then
        strainp(noel,9:12)=dpstran
       elseif (npt==4) then  
        strainp(noel,13:16)=dpstran      
       end if
       statev(1:4)=eelas                  !    Svar (1-4)
       statev(5:8)=eplas                  !    Svar (5-8)
c *** elastic strain and put into epsZZ
      epsZZ    = -xnue/E * (stress(1) + stress(2))
c *** add plastic strain to total strain epsZZ
      epsZZ    = epsZZ + (eplas(3))
      
      epspl=eplas
      epspl(3)=eplas(4)
c      
      q=(2.d0/3.d0)*(sigmae/def)
      r=((h-(deqpl/sigmae))/(def*sigmae))*(3.d0*eg)/(1.d0+3.d0*eg*h)
c
       do i=1,3
        do j=1,3
         ddsdde(i,j)=q*xiden(i,j)+(xk-q*1.d0/3.d0)
     1  -r*dstre(i)*dstre(j)
        end do
       end do
       do k=1,2
        ddsdde(k,3) = -r*dstre(k)*dstre(4)
        ddsdde(3,k) = ddsdde(k,3)
       end do
       ddsdde(3,3) = q/2.d0 - r*dstre(4)*dstre(4) 
       
       kflag=props(6)
       if (kflag.eq.1) then ! fcc
        xm=3.06d0
        b=0.2555d-6
        r=1.9 !arsenlis and parks (1998)
       else ! bcc
        xm=2.9d0
        b=0.2725d-6
        r=1.9 ! ?
       endif
      
       if (ele.eq.0.d0) then
        gnd=0.d0
       else
        gnd=r*eqgrad(noel,npt)/b
       endif    
       ssd=((syield*(E/syield)**ene*(eqplas+syield/E)**ene)/
     & (xm*0.5*eg*b))**2
       td=gnd+ssd
       statev(1+2*ntens+2)=eqplas                   !    Svar (9)
       epseq=eqplas
       
       statev(2+2*ntens+2)=eqgrad(noel,npt)         !    Svar (10) 
       
       statev(3+2*ntens+2)=(1000000.0d0)*ssd 
       ! to plot in m^-2 Svar (11)
       
       statev(4+2*ntens+2)=(1000.0d0)*gnd    
       ! to plot in m^-2 Svar (12)
       
       statev(7+2*ntens+2)=td-statev(5+2*ntens)     !    Svar (15)
       statev(5+2*ntens+2)=td                       !    Svar (13)
       statev(6+2*ntens+2)=deqpl                    !    Svar (14)
      end if

      cordex(noel,npt)=coords(1)
      cordey(noel,npt)=coords(2)    
      
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
      end subroutine
      subroutine kdevia(strss,xidenty,devito)
     
#include 'impcom.inc'

      DOUBLE PRECISION strss(3,3),xidenty(3,3),devito(3,3)
      DOUBLE PRECISION x
      INTEGER i,j
      x=0.0
      do 10 i=1,3
      do 10 j=1,3
      if(i .eq. j) then
      x=x+strss(i,j)
      else
      end if
10    continue


      do 30 i=1,3
      do 30 j=1,3
      if(i .eq. j) then
        devito(i,j)=strss(i,j)-((1./3.)*x*xidenty(i,j))
      else
        devito(i,j)=strss(i,j)
      end if
30    continue
      return
      end
      
      subroutine keff(eff1,val1)
      
#include 'impcom.inc'

      DOUBLE PRECISION eff1(3,3), val1
      DOUBLE PRECISION x
      INTEGER i,j
      
      x=0.0
      do 40 i=1,3
      do 40 j=1,3
       x=x+eff1(i,j)*eff1(i,j)
40    continue

      if(x .le. 0.0) go to 50      
      
      val1=sqrt((2.0/3.0)*x)
50    return
      end   
    
                  SUBROUTINE REVERSE(N,A,C)
            ! input ...
            ! a(n,n) - array of coefficients for matrix A
            ! n      - dimension
            ! output ...
            ! c(n,n) - inverse matrix of A
            ! comments ...
            ! the original matrix a(n,n) will be destroyed 
            ! during the calculation
            !===========================================================
            implicit none 
            integer n
            double precision A(n,n), C(n,n), OOO(n,2*n), lkl(n,2*n) 
            double precision ll,tt, gaus(n,2*n)
            integer i, j, pp, kk, kl
                  gaus = 0.0d0
                  do i = 1, N
                      do j = 1, N
                          gaus(i,j) = A(i,j)
                      end do
                  end do
                  do i=1, N
                      gaus(i,i+N) = 1.0d0
                  end do
                  OOO = gaus
                  do i = 1,N
                      lkl = OOO
                      kl = i+1
                      If(OOO(i,i).EQ.0.0d0) then
                          do while (OOO(i,i).EQ.0.0d0)
                              kl = kl
                          do j = 1, 2*N
                              OOO(i,j) = lkl(kl,j)
                              OOO(kl,j) = lkl(i,j)
                          end do
                          If(OOO(i,i).EQ.0.0d0) then
                              OOO = lkl
                          end if    
                          kl = kl +1
                          end do
                      end if    
                      tt = OOO(i,i)
                      do pp = 1, 2*N
                         OOO(i,pp) = OOO(i,pp)/tt 
                      end do   
                      do kk = 1,N-i
                          ll = OOO(kk+i,i)
                          do j = 1,2*N
                              OOO(kk+i,j) = OOO(kk+i,j)-OOO(i,j)*ll
                          end do
                      end do
                  end do   
c          обратный ход
                  do i = -N,-1
                      tt = OOO(-i,-i)
                      do pp = 1, 2*N
                         OOO(-i,pp) = OOO(-i,pp)/tt 
                      end do   
                      do kk = i+1, -1
                          ll = OOO(-kk,-i)
                          do j = 1,2*N
                              OOO(-kk,j) = OOO(-kk,j)-OOO(-i,j)*ll
                          end do
                      end do
                  end do   
                  do i=1,N
                      do j=1,N
                          C(i,j) = OOO(i,j+N)
                      end do
                  end do    
          end subroutine 



 


