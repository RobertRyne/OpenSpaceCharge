      subroutine SPCH3DAMR( c,ex,ey,ez,Msk,Np,Ntot                      &
     &                     ,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,phi )
      implicit none
!Arguments
      integer Np,Ntot,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj
      real*8,  dimension(6,Np)     :: c
      logical, dimension(Np)       :: Msk
      real*8,  dimension(Np)       :: ex,ey,ez
      real*8,  dimension(Nx,Ny,Nz) :: phi
      write(6,*) 'error: SPCH3DAMR: Chombo AMR is not included '
     &           ,'in this version.'
      stop 'NOAMR'
      end
