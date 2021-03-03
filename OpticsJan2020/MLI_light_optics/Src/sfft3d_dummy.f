! -*- Mode: Fortran; Modified: "Mon 17 Nov 2003 11:36:42 by dbs"; -*- 

! non-operative version of sfft3d() routine for building MLI without FFTW or ESSL

      subroutine sfft3d( f,N )
      implicit none
!Arguments
      integer N
      real*8 f(0:N,0:N,0:N)
      write(6,*) 'error: sfft3d: not implemented. '
     &    ,'Requires FFTW or ESSL.'
      stop 'NOSFFT3D'
      end
