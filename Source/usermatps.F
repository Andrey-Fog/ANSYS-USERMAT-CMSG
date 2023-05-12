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
#include "impcom.inc"   
      EXTERNAL         kdevia, keff
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
      
      DOUBLE PRECISION E, xnue, ebulk3, xk, eg2, eg, elam,
     1 syield, ele, sigmaf,sigmae,h, def, defi
      
      INTEGER k1, k2, i,kewton, k, kflag, j, prev
      
      DOUBLE PRECISION a1, a2, a3, a4, b1, b2, b3, b4, djacb, etat,
     1 rhs, dabs, deqpl, ep, q, r, ene,
     2 xm ,b,gnd, ssd,td
c     branch testps      
c     state variables array for all integration points
      DOUBLE PRECISION SVAR(4,nstatv), SVARtmp(nstatv)
      DOUBLE PRECISION tmp,lamme
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     obtain the gradient values from all integration points
      SVAR=0.d0
      do i=1,4
         call get_ElmData ('SVAR', noel,i, nstatv, SVARtmp)
          do j=1,nstatv
              SVAR(i,j)=SVARtmp(j)
         enddo
      enddo
      call get_ElmData ('EPTO', noel, npt, 4, eelas)
      call get_ElmData ('EPPL', noel, npt, 4, eplas)
      call get_ElmData ('SIG', noel,npt, 3, sigi)
       do i=1,4
          eelas(i) = eelas(i) - eplas (i)
       enddo 
      
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
c      
c     it's needed for calculation of hourglass stiffness,
c     tsstif have no analog in ABAQUS
      eg2=E/(1.d0+xnue)
      eg=eg2/2.d0
c     it's needed for calculation of hourglass stiffness
      tsstif(1)=eg
      lamme=2*xnue*eg/(1-2*xnue)
      elam=E/(1.d0-xnue*xnue)
c      
c     elastic stiffness matrix
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
      xk=ebulk3/3.d0 
   
c     
      if (all(dstran .eq. 0)) then
c     In first iteration we need only elastic stiffness matrix
          continue
      else
c     integration point coordinates in the isoparametric space
        if (npt.eq.1) then
c     Start integration point loop
        do k1=1,4
c     Integration point coordinates in the isoparametric space   
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
   
         xjacm(1,1)= deriv(1,1)*SVAR(1,1)+deriv(1,2)*SVAR(2,1)
     1 +deriv(1,3)*SVAR(3,1)+deriv(1,4)*SVAR(4,1)
    
         xjacm(1,2)= deriv(1,1)*SVAR(1,2)+deriv(1,2)*SVAR(2,2)
     1 +deriv(1,3)*SVAR(3,2)+deriv(1,4)*SVAR(4,2) 
     
         xjacm(2,1)= deriv(2,1)*SVAR(1,1)+deriv(2,2)*SVAR(2,1)
     1 +deriv(2,3)*SVAR(3,1)+deriv(2,4)*SVAR(4,1)
      
         xjacm(2,2)= deriv(2,1)*SVAR(1,2)+deriv(2,2)*SVAR(2,2)
     1 +deriv(2,3)*SVAR(3,2)+deriv(2,4)*SVAR(4,2) 
         
         if (all(xjacm .eq. 0)) then 
              xjaci=0
           else
             djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1) 
      
             xjaci(1,1)=xjacm(2,2)/djacb 
             xjaci(1,2)=-xjacm(1,2)/djacb  
             xjaci(2,1)=-xjacm(2,1)/djacb   
             xjaci(2,2)=xjacm(1,1)/djacb
         end if
         
         a1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1) 
         a2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2) 
         a3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3) 
         a4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4) 
         b1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1) 
         b2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2)  
         b3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3) 
         b4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4) 
         !dn111
         eta(1)=a1*SVAR(1,4) + a2*SVAR(2,4) + a3*SVAR(3,4)
     1        + a4*SVAR(4,4)
         !dn112
         eta(2)=2.d0*(a1*SVAR(1,7)+a2*SVAR(2,7)
     1    + a3*SVAR(3,7) + a4*SVAR(4,7))
     2    - b1*SVAR(1,4) - b2*SVAR(2,4) 
     3    - b3*SVAR(3,4) - b4*SVAR(4,4)
         !dn121
         eta(3)=b1*SVAR(1,4) + b2*SVAR(2,4) + b3*SVAR(3,4)
     1    + b4*SVAR(4,4)
         !dn122    
        eta(4)=a1*SVAR(1,5) + a2*SVAR(2,5) + a3*SVAR(3,5)
     1    + a4*SVAR(4,5)
        !dn133   
        eta(5)=a1*SVAR(1,6) + a2*SVAR(2,6)
     &    + a3*SVAR(3,6) + a4*SVAR(4,6)
        !dn211 
         eta(6)=b1*SVAR(1,4) + b2*SVAR(2,4) + b3*SVAR(3,4)
     1    + b4*SVAR(4,4)
        !dn212
        eta(7)=a1*SVAR(1,5) + a2*SVAR(2,5) + a3*SVAR(3,5)
     &    + a4*SVAR(4,5)
        !dn221  
        eta(8)=2.d0*(b1*SVAR(1,7) + b2*SVAR(2,7) + b3*SVAR(3,7)
     &    + b4*SVAR(4,7)) 
     &    - a1*SVAR(1,5) - a2*SVAR(2,5) - a3*SVAR(3,5)
     &    - a4*SVAR(4,5) 
        !dn222  
        eta(9)=b1*SVAR(1,5) + b2*SVAR(2,5) + b3*SVAR(3,5)
     1    + b4*SVAR(4,5)
        !dn233
         eta(10)=b1*SVAR(1,6) + b2*SVAR(2,6) + b3*SVAR(3,6)
     1    + b4*SVAR(4,6)
        !dn313
         eta(11)=a1*SVAR(1,6) + a2*SVAR(2,6)
     &    + a3*SVAR(3,6) + a4*SVAR(4,6)
        !dn323
         eta(12)= - b1*SVAR(1,6) - b2*SVAR(2,6) 
     &    - b3*SVAR(3,6) - b4*SVAR(4,6)
        !dn331
         eta(13)=
     &    - a1*SVAR(1,6) - a2*SVAR(2,6) - a3*SVAR(3,6)
     &    - a4*SVAR(4,6) 
        !dn332 
         eta(14)=- b1*SVAR(1,6) - b2*SVAR(2,6) 
     &    - b3*SVAR(3,6) - b4*SVAR(4,6)

        etat=eta(1)**2+eta(2)**2+eta(3)**2+eta(4)**2+eta(5)**2+eta(6)**2
     1 +eta(7)**2+eta(8)**2+eta(9)**2+eta(10)**2+eta(11)**2+eta(12)**2
     2 +eta(13)**2+eta(14)**2
         
        SVAR(k1,10) = SVAR(k1,10) + sqrt((1.d0/4.d0)*(etat))
!         eqgrad(noel,k1)=eqgrad(noel,k1)+sqrt((1.d0/4.d0)*(etat))
        enddo
c     Save obtained gradients
          do k1=1,4
            do k2=1,nstatv
              SVARtmp(k2)=SVAR(k1,k2)
            enddo
           call put_ElmData ('SVAR', noel, k1, nstatv, SVARtmp)
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

c     flow stress
      sigmaf=syield*((E/syield)**ene)*sqrt((epseq+(syield/E))**(2*ene)
     1 + ele*SVAR(npt,10))
     
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
      keycut=1
   20 continue

       deqpl=def-sigmae/(3.d0*eg) !3.15 
       epseq=epseq+deqpl
       
       
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

       statev(4:7)= dpstran      
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
        gnd=r*SVAR(npt,10)/b
       endif    
       ssd=((syield*(E/syield)**ene*(epseq+syield/E)**ene)/
     & (xm*0.5*eg*b))**2
       td=gnd+ssd
c *** Save user variable array     
c     gradient 
       statev(10)= SVAR(npt,10)    
c     to plot ssd in m^-2 
       statev(11)=(1000000.0d0)*ssd
c     to plot gnd in m^-2 
       statev(12)=(1000.0d0)*gnd
c     total dislocations
       statev(13)=td 
       
      end if
c      
      statev(1)=coords(1)
      statev(2)=coords(2)
      statev(3)=coords(3)
      
      return
      end
