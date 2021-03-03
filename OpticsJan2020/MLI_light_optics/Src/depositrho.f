! IMPACT 3D space charge routines
! Copyright 2003 University of California
!
!      subroutine depositrho( C,Msk,Np,NpTot,Nx,Ny,Nz,Nadj,rho,rhotmp )
!
! Arguments:
!   C       :in  (1:5:2,Np) are x,y,z coordinates of particles
!                 Note: in the rest of MLI column 5 is time-of-flight (i.e. phase),
!                       but prior to calling the space charge routines it is
!                       converted to longitudinal position z.
!                 (2:6:2,Np) are momenta and not used here.
!   Msk      :in  flag indicating valid particles
!   Np       :in  number of particles (both good and bad) on this processor
!   Nptot    :in  number of (good) particles on all processors
!   Nx,Ny,Nz :in  dimensions of grid on the regular sized grid
!   Nadj     :in  =0 for open or Dirichlet BCs, >=1 for longitudinally periodic BCs
!   rho      :out array of charge density
!   rhotmp   :tmp scratch space
!
! Globals:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine depositrho( C,Msk,Np,NpTot,Nx,Ny,Nz,rho,rhotmp )
!Modules
      use parallel
      use ml_timer
      implicit none
!Arguments
      real*8  C(6,Np)
      logical Msk(Np)
      integer Np,Nptot,Nx,Ny,Nz,Nadj
      real*8  rho(Nx,Ny,Nz) ,rhotmp(Nx,Ny,Nz)
!Local variables
      integer ierror
!Globals
      integer       IVerbose
      common/SHOWME/IVerbose
      real*8          xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/GRIDSZ3D/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!Externals

!Execute
      if( IVerbose .GT. 5 ) write(6,*) 'in DEPOSITRHO'

! deposit charge on the grid:
      call increment_timer('rhoslo3d',0)
      call rhoslo3d(C,rho,Msk,Np,Nx,Ny,Nz,Nadj)
      if( NVP .GT. 1 )then
        call MPI_ALLREDUCE(rho,rhotmp,Nx*Ny*Nz,Mreal,Mpisum,Lworld,ierror)
        if( ierror .NE. 0 ) write(6,*)
     &      'depositrho: MPI_ALLREDUCE returned ',ierror
        rho=rhotmp
      endif
      call increment_timer('rhoslo3d',1)

!--------
!XXX!-- for iexactrho case
!XXX      glrhochk=sum(rho)
!XXX      if(idproc.eq.0)                                                 &
!XXX     &   write(6,*)'(exact rho): global sum of rho = ',glrhochk
!XXX      call getrms(c,xrms,yrms,zrms,np)
!XXX      if(idproc.eq.0)                                                 &
!XXX     &   write(6,*)'(exact rho): xrms,yrms,zrms=',xrms,yrms,zrms
!XXX      xmac=sqrt(5.d0)*xrms
!XXX      ymac=sqrt(5.d0)*yrms
!XXX      zmac=sqrt(5.d0)*zrms
!XXX      if(idproc.eq.0)                                                 &
!XXX     &   write(6,*)'(exact rho): xmac,ymac,zmac=',xmac,ymac,zmac
!XXX      rho=0.d0
!XXX      do k=1,nz
!XXX        z=zmin+(k-1)*hz
!XXX        do j=1,ny
!XXX          y=ymin+(j-1)*hy
!XXX          do i=1,nx
!XXX            x=xmin+(i-1)*hx
!XXX!     if((x/xbig)**2+(y/ybig)**2+(z/zbig)**2 .le.1.d0)rho(i,j,k)=1.d0
!XXX      if((x/xmac)**2+(y/ymac)**2+(z/zmac)**2 .le.1.d0)rho(i,j,k)=1.d0
!XXX          enddo
!XXX        enddo
!XXX      enddo
!XXX      glrhonew=sum(rho)
!XXX      if(idproc.eq.0)                                                 &
!XXX     &   write(6,*)'(exact rho): glrhonew=',glrhonew
!XXX      rho=rho*glrhochk/glrhonew
!XXX      glrhochk=sum(rho)
!XXX      if(idproc.eq.0)                                                 &
!XXX     &   write(6,*)'(exact rho): new global sum of rho = ',glrhochk
!XXX!-- end

!--------
cryne august 1, 2002
cryne moved normalization out of rhoslo to here:
      rho=rho/(hx*hy*hz*ntot)
      call system_clock(count=iticks1)
