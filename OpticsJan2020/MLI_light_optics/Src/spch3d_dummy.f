      subroutine SPCH3D2( c,ex,ey,ez,Msk,Np,Ntot                        &
     &                   ,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,rhophi )
      implicit none
!Arguments
      integer Np,Ntot,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj
      real*8, dimension(6,Np) :: c
      logical, dimension(Np) :: Msk
      real*8, dimension(Np) :: ex,ey,ez
      real*8, dimension(Nx,Ny,Nz) :: rhophi
      write(6,*) 'error: SPCH3D2: ANAG Poisson is not included '
     &           ,'in this version.'
      stop 'NOANAG'
      return
      end

      subroutine SPCH3DBC0( c,ex,ey,ez,Msk,Np,Ntot                      &
     &                     ,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,rho )
!Arguments
      implicit none!   double precision(a-h,o-z)
      integer Np,Ntot,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj
      real*8, dimension(6,Np) :: c
      logical, dimension(Np) :: Msk
      real*8, dimension(Np) :: ex,ey,ez
      real*8, dimension(Nx,Ny,Nz) :: rho
      write(6,*) 'error: SPCH3D2: ANAG Poisson is not included '
     &           ,'in this version.'
      stop 'NOANAG'
      return
      end
