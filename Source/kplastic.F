      module kplastic
      implicit none
c     The dimension of these arrays must be equivalent to
c     elements count or will be greater (eqgrad (noel,npts))
      real*8  cordex(100000,4),cordey(100000,4),strainp(100000,16),
     &        eqgrad(100000,4) 
      
      save
      end module
      
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
c
c     Main function: Return the inverse matrix
c     input ...
c     a(n,n) - array of coefficients for matrix A
c     n      - dimension
c     output ...
c     c(n,n) - inverse matrix of A
c     comments ...
c     the original matrix a(n,n) will be destroyed 
c     during the calculation
c *** ===========================================================
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
c          change way
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
