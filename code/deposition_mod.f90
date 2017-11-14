module deposition_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine get_rho_and_mesh_spacing(ptcls,rho,dx,dy,dz,xmin,ymin,zmin,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi, &
     & ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,idecomp,npx,npy,npz,n1,nraysp,maxrayp)
      use mpi
      use data_movement_mod, only : localize
      implicit none
      integer :: ilo,ihi,jlo,jhi,klo,khi,ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,idecomp,npx,npy,npz,n1,nraysp,maxrayp
      real(dp) :: dx,dy,dz,xmin,ymin,zmin,chrgpermacro
      real(dp), dimension(n1,maxrayp) :: ptcls
      real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
!
      real(dp), dimension(maxrayp) :: xa,ya,za,lostflag !lostflag=1.0 if particle is "lost"
      real(dp) :: xmax,ymax,zmax !not needed
      integer :: ifail,n
      integer :: mprocs,myrank,ierr
      real(dp) :: xsml,xbig,ysml,ybig,zsml,zbig
      real(dp), parameter :: eps=1.d-15   !3.34d-16 is OK on my Mac
!     real(dp), parameter :: eps=0.5d0
      integer :: nx,ny,nz !temporaries
!
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!     if(myrank.eq.0)write(6,*)'hello from get_rho'
!     write(100+myrank,'(6(i5,1x))')ilo,ihi,jlo,jhi,klo,khi
!
      xa(1:nraysp)=ptcls(1,1:nraysp)
      ya(1:nraysp)=ptcls(3,1:nraysp)
      za(1:nraysp)=ptcls(5,1:nraysp)
!?    lostflag(1:nraysp)=ptcls(7,1:nraysp)
      lostflag(1:maxrayp)=ptcls(7,1:maxrayp)
! compute the bounding box and dx,dy,dz so particles can be localized to the correct proc:
      call getbeamboundingbox(xa,ya,za,lostflag,xsml,xbig,ysml,ybig,zsml,zbig,nraysp) !nraysp should be bigrayp?
      xbig=xbig*(1.d0+sign(1.d0,xbig)*eps)
      xsml=xsml*(1.d0-sign(1.d0,xsml)*eps)
      ybig=ybig*(1.d0+sign(1.d0,ybig)*eps)
      ysml=ysml*(1.d0-sign(1.d0,ysml)*eps)
      zbig=zbig*(1.d0+sign(1.d0,zbig)*eps)
      zsml=zsml*(1.d0-sign(1.d0,zsml)*eps)
      nx=ihi_gbl-ilo_gbl+1 !fixed (had forgotten gbl)
      ny=jhi_gbl-jlo_gbl+1 !fixed (had forgotten gbl)
      nz=khi_gbl-klo_gbl+1 !fixed (had forgotten gbl)
      dx=(xbig-xsml)/(nx-3)
      dy=(ybig-ysml)/(ny-3)
      dz=(zbig-zsml)/(nz-3)
      xmin=xsml-dx; xmax=xbig+dx
      ymin=ysml-dy; ymax=ybig+dy
      zmin=zsml-dz; zmax=zbig+dz
!-----
!     call getbeamboundingbox(xa,ya,za,lostflag,xmin,xmax,ymin,ymax,zmin,zmax,nraysp)
!!!!!!dx=(xmax-xmin+1.d-15)/((ihi_gbl-ilo_gbl+1)/2-3) !use when (ilo_gbl,ihi_gbl) denotes the doubled grid
!!!!!!dy=(ymax-ymin+1.d-15)/((jhi_gbl-jlo_gbl+1)/2-3)
!!!!!!dz=(zmax-zmin+1.d-15)/((khi_gbl-klo_gbl+1)/2-3)
!     dx=(xmax-xmin+1.d-15)/((ihi_gbl-ilo_gbl+1)-3)   !use when (ilo_gbl,ihi_gbl) denotes the physical grid
!     dy=(ymax-ymin+1.d-15)/((jhi_gbl-jlo_gbl+1)-3)
!     dz=(zmax-zmin+1.d-15)/((khi_gbl-klo_gbl+1)-3)
!     xmin=xmin-dx
!     ymin=ymin-dy
!     zmin=zmin-dz !last three statements, along with "-3" above, ensure no contribution to boundary grid points
!     xmax=xmin+(ihi_gbl-ilo_gbl+1-1)*dx
!     ymax=ymin+(jhi_gbl-jlo_gbl+1-1)*dy
!     zmax=zmin+(khi_gbl-klo_gbl+1-1)*dz
!     if(myrank.eq.0)then
!       write(6,*)'(get_rho_and_mesh_spacing) global bounding box:'
!       write(6,*)'xmin,xmax,dx=',xmin,xmax,dx
!       write(6,*)'ymin,ymax,dy=',ymin,ymax,dy
!       write(6,*)'zmin,zmax,dz=',zmin,zmax,dz
!       write(6,*)' '
!       write(6,*)'n1,nraysp,maxrayp=',n1,nraysp,maxrayp
!       write(6,*)'ilo_gbl,ihi_gbl=',ilo_gbl,ihi_gbl
!       write(6,*)'jlo_gbl,jhi_gbl=',jlo_gbl,jhi_gbl
!       write(6,*)'klo_gbl,khi_gbl=',klo_gbl,khi_gbl
!       write(6,*)'idecomp=',idecomp
!       write(6,*)'npx,npy,npz=',npx,npy,npz
!       write(6,*)'mprocs=',mprocs
!     endif
!
! move particles to the correct proc:
!xxx  ptcls(1,1:nraysp)=xa(1:nraysp)
!xxx  ptcls(2,1:nraysp)=ya(1:nraysp)
!xxx  ptcls(3,1:nraysp)=za(1:nraysp)
!     if(myrank.eq.0)write(6,*)'calling localize'
      call localize(ptcls,lostflag,xmin,ymin,zmin,dx,dy,dz,&
                    ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,n1,nraysp,maxrayp,idecomp,npx,npy,npz,mprocs)
!     if(myrank.eq.0)write(6,*)'back from localize'
!!deposit charge on the grid:
      xa(1:nraysp)=ptcls(1,1:nraysp)
      ya(1:nraysp)=ptcls(3,1:nraysp)
      za(1:nraysp)=ptcls(5,1:nraysp)
!
      lostflag(1:maxrayp)=1.d0
      lostflag(1:nraysp)=0.d0    !all particles are included in this example run
!!    call depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,dx,dy,dz,xmin,ymin,zmin,nraysp, &
!!   &                              ilo_gbl,jlo_gbl,klo_gbl,ifail)
!!    if(myrank.eq.0)write(6,*)'back from depose_rho_scalar'
      call newdepose_rho_scalar_NGP(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,dx,dy,dz,xmin,ymin,zmin,nraysp, &
     &                              ilo_gbl,jlo_gbl,klo_gbl,ifail)
      return
      end
!
      subroutine getbeamboundingbox(x,y,z,lostflag,xmin,xmax,ymin,ymax,zmin,zmax,nraysp)
      use mpi
      implicit none
      integer :: nraysp !# of particles per MPI process
      real(dp), dimension(*) :: x,y,z,lostflag !lostflag=1.0 if particle is "lost"
      real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax
      real(dp), dimension(6) :: veclcl,vecgbl
      real(dp) :: xminlcl,xmaxlcl,yminlcl,ymaxlcl,zminlcl,zmaxlcl
      real(dp) :: xwidthorig,ywidthorig,zwidthorig
      integer :: myrank,mpierr,ierror
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
!need this since, if nraysp=0, the next 6 statements are skipped
      veclcl(1:3)=-9999999.   !looks like a bug; a proc wiht nraysp=0 should not be included in ALLREDUCE below
      veclcl(4:6)=-9999999.
      if(nraysp.eq.0)goto 100
!     veclcl(1)=-minval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(2)=-minval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(3)=-minval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(4)=maxval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(5)=maxval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(6)=maxval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
      veclcl(1)=-minval(x(1:nraysp))  !deal with lostflag later; for testing, this should work
      veclcl(2)=-minval(y(1:nraysp))
      veclcl(3)=-minval(z(1:nraysp))
      veclcl(4)=maxval(x(1:nraysp))
      veclcl(5)=maxval(y(1:nraysp))
      veclcl(6)=maxval(z(1:nraysp))
  100 continue
      call MPI_ALLREDUCE(veclcl,vecgbl,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
      xmin=-vecgbl(1)
      ymin=-vecgbl(2)
      zmin=-vecgbl(3)
      xmax=vecgbl(4)
      ymax=vecgbl(5)
      zmax=vecgbl(6)
!hack
!make the grid bigger for testing purposes, e.g. see 1/r behavior of potential outside of a uniform beam
      xwidthorig=xmax-xmin
      ywidthorig=ymax-ymin
      zwidthorig=zmax-zmin
      xmin=xmin-0.5d0*xwidthorig
      xmax=xmax+0.5d0*xwidthorig
      ymin=ymin-0.5d0*ywidthorig
      ymax=ymax+0.5d0*ywidthorig
      zmin=zmin-0.5d0*zwidthorig
      zmax=zmax+0.5d0*zwidthorig
!end hack
      return
      end

!routines for charge deposition
      subroutine depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,      &
     &                                    dx,dy,dz,xmin,ymin,zmin,nraysp,ilogbl,jlogbl,klogbl,ifail)
      implicit none
      integer, intent(in) :: nraysp !# of particles per MPI process
      integer, intent(in) :: ilogbl,jlogbl,klogbl
      integer, intent(out) :: ifail
      real(dp), intent(in) :: chrgpermacro
      real(dp), intent(in), dimension(nraysp) :: xa,ya,za,lostflag
      integer :: ilo,jlo,klo,ihi,jhi,khi
      real(dp), intent(out), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
      real(dp) :: dx,dy,dz,xmin,ymin,zmin
      real(dp) :: dxi,dyi,dzi,ab,de,gh,sumrho
      integer :: n,ip,jp,kp
      integer :: myrank,ierr
      ifail=0
      dxi=1.d0/dx
      dyi=1.d0/dy
      dzi=1.d0/dz
      rho(ilo:ihi,jlo:jhi,klo:khi)=0.d0
      do n=1,nraysp
        if(lostflag(n).ne.0.d0)cycle
        ip=floor((xa(n)-xmin)*dxi+1) !this is 1-based; use ilogbl for the general case
        jp=floor((ya(n)-ymin)*dyi+1)
        kp=floor((za(n)-zmin)*dzi+1)
        ab=((xmin-xa(n))+ip*dx)*dxi
        de=((ymin-ya(n))+jp*dy)*dyi
        gh=((zmin-za(n))+kp*dz)*dzi
! this "if" statement slows things down, but I'm using it for debugging purposes
        if(ip.lt.ilo.or.jp.lt.jlo.or.kp.lt.klo.or.ip.gt.ihi.or.jp.gt.jhi.or.kp.gt.khi)then
          ifail=ifail+1
          cycle
        else
          rho(ip,jp,kp)=rho(ip,jp,kp)+ab*de*gh
          rho(ip,jp+1,kp)=rho(ip,jp+1,kp)+ab*(1.-de)*gh
          rho(ip,jp+1,kp+1)=rho(ip,jp+1,kp+1)+ab*(1.-de)*(1.-gh)
          rho(ip,jp,kp+1)=rho(ip,jp,kp+1)+ab*de*(1.-gh)
          rho(ip+1,jp,kp+1)=rho(ip+1,jp,kp+1)+(1.-ab)*de*(1.-gh)
          rho(ip+1,jp+1,kp+1)=rho(ip+1,jp+1,kp+1)+(1.-ab)*(1.-de)*(1.-gh)
          rho(ip+1,jp+1,kp)=rho(ip+1,jp+1,kp)+(1.-ab)*(1.-de)*gh
          rho(ip+1,jp,kp)=rho(ip+1,jp,kp)+(1.-ab)*de*gh
        endif
      enddo
      rho=rho*chrgpermacro
      if(ifail.ne.0)write(6,*)'(depose_rho_scalar) ifail=',ifail
      return
      end
!
      subroutine newdepose_rho_scalar_NGP(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi, &
     &                                    dx,dy,dz,xmin,ymin,zmin,nraysp,ilogbl,jlogbl,klogbl,ifail)
      use mpi
      implicit none
      integer, intent(in) :: nraysp !# of particles per MPI process
      integer, intent(in) :: ilogbl,jlogbl,klogbl
      integer, intent(out) :: ifail
      real(dp), intent(in) :: chrgpermacro
!!!!! real(dp), intent(in), dimension(*) :: xa,ya,za,lostflag !lostflag=1.0 if particle is "lost"
      real(dp), intent(in), dimension(nraysp) :: xa,ya,za,lostflag !lostflag=1.0 if particle is "lost"
      integer :: ilo,jlo,klo,ihi,jhi,khi
      real(dp), intent(out), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
      real(dp) :: dx,dy,dz,xmin,ymin,zmin
      real(dp) :: dxi,dyi,dzi,ab,de,gh,rhosum,rhosumgbl
      integer :: ip,i,j,k
      integer :: myrank,ierr
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      ifail=0
      dxi=1.d0/dx
      dyi=1.d0/dy
      dzi=1.d0/dz
      rho(ilo:ihi,jlo:jhi,klo:khi)=0.d0
      do ip=1,nraysp
        if(lostflag(ip).ne.0.d0)cycle
        ab=(xa(ip)-xmin)*dxi
        de=(ya(ip)-ymin)*dyi
        gh=(za(ip)-zmin)*dzi
        i=nint(ab)+ilogbl !floor(ab)+ilogbl
        j=nint(de)+jlogbl !floor(de)+jlogbl
        k=nint(gh)+klogbl !floor(gh)+klogbl
! this "if" statement slows things down, but I'm using it for debugging purposes
        if(i.lt.ilo.or.j.lt.jlo.or.k.lt.klo.or.i.gt.ihi.or.j.gt.jhi.or.k.gt.khi)then
          ifail=ifail+1
          cycle
        else
          rho(i,j,k)=rho(i,j,k)+1.d0
        endif
      enddo
      rho=rho*chrgpermacro
      rhosum=sum(rho)
      call MPI_ALLREDUCE(rhosum,rhosumgbl,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0)write(6,*)'global rhosum=',rhosumgbl
      if(ifail.ne.0)write(6,*)'(newdepose_rho_scalar_NGP) ifail=',ifail
      return
      end

end module
