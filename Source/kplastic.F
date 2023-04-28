      module kplastic
      implicit none
c     Dimension of theese arrays must be equivalent to
c                     elements count or will be greater (eqgrad (noel,npts))
      real*8  cordex(100000,4),cordey(100000,4),strainp(100000,16),
     &        eqgrad(100000,4) 
      
      save
      end module