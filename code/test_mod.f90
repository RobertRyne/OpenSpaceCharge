!This module contains extra routines (not related to OpenSC) used by the test code.
module test_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains





      subroutine gendist(ptcls,n1,nraysp,sigmat,gaussiancutoff,disttype,iseed,apipe,bpipe,rectpipe)
     ! use mpi
      implicit none
      integer :: n1,nraysp,disttype,iseed
      real(8) :: apipe,bpipe
      logical :: rectpipe
      real(8), dimension(nraysp,n1) :: ptcls
      real(8), dimension(6,6) :: sigmat
      real(8) :: gaussiancutoff
      integer :: n
      real(8) :: r1,r2,r3,r4,r5,r6,arg
      real(8), dimension(6) :: cent
      integer :: mprocs,myrank,ierr
      
     ! call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
     ! call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

      if(myrank.eq.0)print *, 'gendist'
      
      call initrandom(iseed)
      do n=1,nraysp
        if(disttype.eq.0)then
   10     continue
          call random_number(r1)
          call random_number(r3)
          call random_number(r5)
          r1=2.d0*r1-1.d0
          r3=2.d0*r3-1.d0
          r5=2.d0*r5-1.d0
          if(r1*r1+r3*r3+r5*r5.gt.1.d0)goto 10
          if(rectpipe)then
            if(r1**2*sigmat(1,1)/(0.5*apipe)**2+r3**2*sigmat(3,3)/(0.5*bpipe)**2.gt.0.9**2)goto 10
          endif
          r2=0.d0; r4=0.d0; r6=0.d0
        endif
        if(disttype.eq.1)then
   11     continue
          call normdv(r1,r2)
          call normdv(r3,r4)
          call normdv(r5,r6)
          if(r1*r1+r3*r3+r5*r5.gt.gaussiancutoff**2)goto 11
          if(rectpipe)then
            if(r1**2*sigmat(1,1)/(0.5*apipe)**2+r3**2*sigmat(3,3)/(0.5*bpipe)**2.gt.0.9**2)goto 11
          endif
        endif
        
! print *, 'gendist 1, n = ', n, r1*sqrt(sigmat(1,1))
! this version assumes that sigmat is 2x2 block diagonal
! this is a mess. don't know why I did it this way. Replace later. RDR
!x-gbx:
        ptcls(n,1)=r1*sqrt(sigmat(1,1))
       ! print *, 'gendist 2'
        ptcls(n,2)=r1*sigmat(1,2)/sqrt(sigmat(1,1))+r2*sqrt(sigmat(2,2)-sigmat(1,2)**2/sigmat(1,1))
!y-gby:
        ptcls(n,3)=r3*sqrt(sigmat(3,3))
        ptcls(n,4)=r3*sigmat(3,4)/sqrt(sigmat(3,3))+r4*sqrt(sigmat(4,4)-sigmat(3,4)**2/sigmat(3,3))
!z-gbz:
        ptcls(n,5)=r5*sqrt(sigmat(5,5))
        ptcls(n,6)=r5*sigmat(5,6)/sqrt(sigmat(5,5))+r6*sqrt(sigmat(6,6)-sigmat(5,6)**2/sigmat(5,5))
      enddo

      do n=1,6
        cent(n)=sum(ptcls(1:nraysp,n))/nraysp
      enddo
      do n=1,nraysp
        ptcls(n,1:6)=ptcls(n,1:6)-cent(1:6)
      enddo
      ptcls(1,1:5)=0.d0 !zero out the first particle to see how the origin gets kicked around by space charge
      write(6,*)'particle xmin,xmax=',minval(ptcls(1:nraysp,1)),maxval(ptcls(1:nraysp,1))
      write(6,*)'particle ymin,ymax=',minval(ptcls(1:nraysp,3)),maxval(ptcls(1:nraysp,3))
      write(6,*)'particle zmin,zmax=',minval(ptcls(1:nraysp,5)),maxval(ptcls(1:nraysp,5))
      return
      
       if(myrank.eq.0)print *, 'gendist end'
      end subroutine gendist
!
      subroutine normdv(d1,d2)
! routine to generate 2 normal deviates
      implicit none
      real(8) :: u1,u2
      real(8) twopi,d1,d2
      twopi=4.0d0*asin(1.0d0)
      call random_number(u1)
      call random_number(u2)
      d1=sqrt(-2.0d0*log(u1))*cos(twopi*u2)
      d2=sqrt(-2.0d0*log(u1))*sin(twopi*u2)
      return
      end subroutine normdv
!
      subroutine initrandom(inpseed)
    !  use mpi
      implicit none
      integer :: inpseed
      integer :: iseedsize,n,i
      integer, dimension(:), allocatable :: iputarray
      integer :: myrank,ierr
   !   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
! determine the seed size and allocate the seed array:
      call random_seed(size=iseedsize)
      allocate(iputarray(iseedsize))
! set a unique seed for each proc:
      n=inpseed
      do i=1,iseedsize*(myrank+1)
        n=mod(8121*n+28411,134456) !I don't remember why I did this
      enddo
      do i=1,iseedsize
        n=mod(8121*n+28411,134456)
        iputarray(i)=n
      enddo
      call random_seed(put=iputarray)
!?    deallocate(iputarray)
      return
      end subroutine initrandom
      
      
      
      !diagnostics routines
      subroutine prntall(nstep,n1,nraysp,nx,ny,nz,ptcls,efield,bfield,t,delta,umin,rectpipe)
      implicit none
      real(8), dimension(3) :: delta,umin
      integer :: nstep,n1,nraysp,nx,ny,nz
      logical :: rectpipe
      real(8) :: t,dx,dy,dz,xmin,ymin,zmin
      real(8), dimension(nx,ny,nz,3) :: efield,bfield
      real(8), dimension(nraysp,n1) :: ptcls
      character(32) :: pname,fname,fxname,fyname,fzname
      integer, parameter :: ndigits=4         !digits in filename suffix corr. to step#
      integer, parameter :: punit=2,funit=3,fxunit=7,fyunit=8,fzunit=9   !file unit numbers
      integer, parameter :: maxptclprnt=10000 !max# of particles printed
      
      dx=delta(1); dy=delta(2); dz=delta(3)
      xmin=umin(1); ymin=umin(2); zmin=umin(3)
      
      pname='ptcls' !prefix for particle filename
      if(rectpipe)then
        fname='fieldpipe' !prefix for field filename
        fxname='xlinepipe' !field vs x
        fyname='ylinepipe' !field vs y
        fzname='zlinepipe' !field vs z
      else
        fname='fieldfree'
        fxname='xlinefree'
        fyname='ylinefree'
        fzname='zlinefree'
      endif
!print particles:
      call openfile(pname,punit,nstep,ndigits)
      
      
      call pprnt(ptcls,n1,nraysp,maxptclprnt,punit)
      close(punit)
!print fields:
      call openfile(fname,funit,nstep,ndigits)
      call openfile(fxname,fxunit,nstep,ndigits)
      call openfile(fyname,fyunit,nstep,ndigits)
      call openfile(fzname,fzunit,nstep,ndigits)
      
      
      call ebprnt(efield(1,1,1,1),efield(1,1,1,2),efield(1,1,1,3),bfield(1,1,1,1),bfield(1,1,1,2),bfield(1,1,1,3), &
                  nx,ny,nz,funit,fxunit,fyunit,fzunit,dx,dy,dz,xmin,ymin,zmin)
      close(funit); close(fxunit); close(fyunit); close(fzunit)
      return
      end
      
      subroutine pprnt(ptcls,n1,nraysp,maxptclprnt,nfile)
      implicit none
      integer :: n1,nraysp,maxptclprnt,nfile
      real(8),  dimension(:,:) :: ptcls
      integer :: n
      do n=1,min(nraysp,maxptclprnt)
        write(nfile,'(20(1pg14.7,1x))')ptcls(n,1:n1)
      enddo
      return
      end
      
!old routine where the components of e and b are in separate 1d arrays
      subroutine ebprnt(ex,ey,ez,bx,by,bz,nx,ny,nz,funit,fxunit,fyunit,fzunit,dx,dy,dz,xmin,ymin,zmin)
      implicit none
      integer :: nx,ny,nz,funit,fxunit,fyunit,fzunit
      real(8) :: dx,dy,dz,xmin,ymin,zmin
      real(8), dimension(nx,ny,nz) :: ex,ey,ez,bx,by,bz
      integer :: i,j,k
! fields on the entire grid:
      do k=1,nz
      do j=1,ny
      do i=1,nx
        write(funit,'(9(1pg14.7,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    ex(i,j,k),ey(i,j,k),ez(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k),i,j,k
      enddo
      enddo
      enddo
!    
! fields vs x along a line: 
      j=ny/2+1   !   1-based, odd
      k=nz/2+1   !   1-based, odd
      do i=1,nx
        write(fxunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    ex(i,j,k),ey(i,j,k),ez(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k),i,j,k
      enddo
!
! fields vs y along a line:
      i=nx/2+1   !   1-based, odd
      k=nz/2+1   !   1-based, odd
      do j=1,ny
        write(fyunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    ex(i,j,k),ey(i,j,k),ez(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k),i,j,k
      enddo
!
! fields vs z along a line:
      i=nx/2+1   !   1-based, odd
      j=ny/2+1   !   1-based, odd
      do k=1,nz
        write(fzunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    ex(i,j,k),ey(i,j,k),ez(i,j,k),bx(i,j,k),by(i,j,k),bz(i,j,k),i,j,k
      enddo
      return
      end
!
      subroutine openfile(fname,nunit,numsuffix,ndigits)
      implicit none
      character(32) :: fname
      integer :: nunit,numsuffix,ndigits
!local variables:
      character(5) aseq
      integer :: j
      integer :: idproc
      idproc=0
!
      if(idproc.eq.0)then
        if(ndigits.gt.0)then
          call num2string(numsuffix,aseq,ndigits)
          j=len_trim(fname)
          fname=fname(1:j)//aseq(1:ndigits)
        endif
        open(unit=nunit,file=fname,status='unknown',err=110)
        goto 150
  110   continue
        write(66,*)'something wrong; error opening file ',fname
      endif
  150 continue
      return
      end
!
      subroutine num2string(num,a,ndigits)
! converts an integer "num" with at most "ndigits" digits
! into a character string "a"
! if ndigits exceeds the actual number of digits in num, then leading
! zeroes are inserted in the string
      integer num,ndigits
      character a*(*)
      integer n,m,k
      m=num
      do n=1,ndigits
      k=m/10**(ndigits-n)
      if(k.eq.0)a(n:n)='0'
      if(k.eq.1)a(n:n)='1'
      if(k.eq.2)a(n:n)='2'
      if(k.eq.3)a(n:n)='3'
      if(k.eq.4)a(n:n)='4'
      if(k.eq.5)a(n:n)='5'
      if(k.eq.6)a(n:n)='6'
      if(k.eq.7)a(n:n)='7'
      if(k.eq.8)a(n:n)='8'
      if(k.eq.9)a(n:n)='9'
      m=m-k*10**(ndigits-n)
      enddo
      return
      end   

subroutine get_mesh_quantities(xa,ya,za,lostflag, delta,umin, nlo,nhi,nlo_gbl,nhi_gbl,n1,nraysp,maxrayp,isetzmin)
!-!#ifdef MPIPARALLEL
!     USE mpi
!-!#endif
implicit none
integer, dimension(3) :: nlo,nhi,nlo_gbl,nhi_gbl
real(dp), dimension(3) :: umin,delta
integer :: n1,nraysp,maxrayp,isetzmin
!-! real(dp), dimension(maxrayp,n1) :: ptcls
!type (coord_struct) :: ptcls(maxrayp)
!
real(dp), dimension(maxrayp) :: xa,ya,za,lostflag !lostflag=1.0 if particle is "lost"
real(dp) :: xmax,ymax,zmax !not needed
integer :: ifail,n
integer :: mprocs,myrank,ierr
real(dp) :: xsml,xbig,ysml,ybig,zsml,zbig
      real(dp), parameter :: eps=1.d-13 ! 2.d-15   !3.34d-16 is OK on my Mac
!real(dp), parameter :: eps=0.5d0
integer :: nx,ny,nz !temporaries
!-!#ifdef MPIPARALLEL
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!-!#else
mprocs=1
myrank=0
!-!#endif
!     if(myrank.eq.0)write(6,*)'hello from get_rho'
      !xa(1:maxrayp)=ptcls(:)%vec(1) !-! ptcls(1:maxrayp,1)
      !ya(1:maxrayp)=ptcls(:)%vec(3) !-! ptcls(1:maxrayp,3)
      !za(1:maxrayp)=ptcls(:)%vec(5) !-! ptcls(1:maxrayp,5)
      !lostflag(1:maxrayp)= ptcls(:)%state  !-! ptcls(1:maxrayp,7)
! compute the bounding box and dx,dy,dz so particles can be localized to the correct proc:
call getbeamboundingbox(xa,ya,za,lostflag,xsml,xbig,ysml,ybig,zsml,zbig,nraysp) !nraysp should be bigrayp?
xbig=xbig*(1.d0+sign(1.d0,xbig)*eps)
xsml=xsml*(1.d0-sign(1.d0,xsml)*eps)
ybig=ybig*(1.d0+sign(1.d0,ybig)*eps)
ysml=ysml*(1.d0-sign(1.d0,ysml)*eps)
zbig=zbig*(1.d0+sign(1.d0,zbig)*eps)
zsml=zsml*(1.d0-sign(1.d0,zsml)*eps)
nx=nhi(1)-nlo(1)+1
ny=nhi(2)-nlo(2)+1
nz=nhi(3)-nlo(3)+1
delta(1)=(xbig-xsml)/(nx-3)
delta(2)=(ybig-ysml)/(ny-3)
delta(3)=(zbig-zsml)/(nz-3)
umin(1)=xsml-delta(1)  !why did I put xmax here? not needed.  ; xmax=xbig+delta(1)
umin(2)=ysml-delta(2)  !why did I put ymax here? not needed.  ; ymax=ybig+delta(2)
umin(3)=zsml-delta(3)  !why did I put zmax here? not needed.  ; zmax=zbig+delta(3)
if(isetzmin.eq.1)then !for cathode images, use this to see fields at the cathode surface
  zsml=1.d-9
  delta(3)=(zbig-zsml)/(nz-1)
  umin(3)=zsml
endif
!
!
!-!#ifdef MPIPARALLEL
! move particles to the correct proc:
!     ptcls(1,1:nraysp)=xa(1:nraysp)
!     ptcls(2,1:nraysp)=ya(1:nraysp)
!     ptcls(3,1:nraysp)=za(1:nraysp)
!     if(myrank.eq.0)write(6,*)'calling localize'
!     call localize(ptcls,lostflag,xmin,ymin,zmin,dx,dy,dz,&
!    &              ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,n1,nraysp,maxrayp,idecomp,npx,npy,npz,mprocs)
!     if(myrank.eq.0)write(6,*)'back from localize'
!!deposit charge on the grid:
!     xa(1:nraysp)=ptcls(1,1:nraysp)
!     ya(1:nraysp)=ptcls(2,1:nraysp)
!     za(1:nraysp)=ptcls(3,1:nraysp)
!-!#endif

end subroutine
!

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine getbeamboundingbox(x,y,z,lostflag,xmin,xmax,ymin,ymax,zmin,zmax,nraysp)
!-!#ifdef MPIPARALLEL
!     USE mpi
!-!#endif
implicit none
integer :: nraysp !!-!# of particles per MPI process
real(dp), dimension(*) :: x,y,z,lostflag !lostflag=1.0 if particle is "lost"
real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax
real(dp), dimension(6) :: veclcl,vecgbl
real(dp) :: xminlcl,xmaxlcl,yminlcl,ymaxlcl,zminlcl,zmaxlcl
real(dp) :: xwidthorig,ywidthorig,zwidthorig
integer :: ierror
!need this since, if nraysp=0, the next 6 statements are skipped
veclcl(1:3)=-9999999.
veclcl(4:6)=-9999999.
if(nraysp > 0) then
!     veclcl(1)=-minval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(2)=-minval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(3)=-minval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(4)=maxval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(5)=maxval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(6)=maxval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
  veclcl(1)=-minval(x(1:nraysp))
  veclcl(2)=-minval(y(1:nraysp))
  veclcl(3)=-minval(z(1:nraysp))
  veclcl(4)=maxval(x(1:nraysp))
  veclcl(5)=maxval(y(1:nraysp))
  veclcl(6)=maxval(z(1:nraysp))
endif
!-!#ifdef MPIPARALLEL
!     call MPI_ALLREDUCE(veclcl,vecgbl,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
!-!#else
vecgbl(1:6)=veclcl(1:6)
!-!#endif
xmin=-vecgbl(1)
ymin=-vecgbl(2)
zmin=-vecgbl(3)
xmax=vecgbl(4)
ymax=vecgbl(5)
zmax=vecgbl(6)

end subroutine

!routine for charge deposition
subroutine depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,nlo,delta,umin,nraysp,nlogbl,ifail)
!subroutine depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,      &
!                             delta,umin,nraysp,ilogbl,jlogbl,klogbl,ifail)
!use mpi
implicit none
integer, intent(in) :: nraysp !!-!# of particles per MPI process
integer, intent(out) :: ifail
real(dp), intent(in) :: chrgpermacro
real(dp), intent(in), dimension(:) :: xa,ya,za,lostflag
integer, intent(in), dimension(3) :: nlo,nlogbl
real(dp), intent(in), dimension(3) :: delta,umin
integer :: ilo,jlo,klo,ihi,jhi,khi
real(dp), intent(out), dimension(nlo(1):,nlo(2):,nlo(3):) :: rho
real(dp) :: dx,dy,dz,xmin,ymin,zmin
real(dp) :: dxi,dyi,dzi,ab,de,gh,sumrho
integer :: ilogbl,jlogbl,klogbl
integer :: n,ip,jp,kp
integer :: mprocs,myrank,ierr

!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

ilo=nlo(1); jlo=nlo(2); klo=nlo(3); ihi=size(rho,1); jhi=size(rho,2); khi=size(rho,3)
dx=delta(1); dy=delta(2); dz=delta(3)
xmin=umin(1); ymin=umin(2); zmin=umin(3)
ilogbl=nlogbl(1); jlogbl=nlogbl(2); klogbl=nlogbl(3)

ifail=0
dxi=1.d0/dx
dyi=1.d0/dy
dzi=1.d0/dz
rho(ilo:ihi,jlo:jhi,klo:khi)=0.d0
do n=1,nraysp
  if(lostflag(n).ne.0.d0)cycle
! ip=floor((xa(n)-xmin)*dxi+1) !this is 1-based; use ilogbl for the general case
! jp=floor((ya(n)-ymin)*dyi+1)
! kp=floor((za(n)-zmin)*dzi+1)
! ab=((xmin-xa(n))+ip*dx)*dxi
! de=((ymin-ya(n))+jp*dy)*dyi
! gh=((zmin-za(n))+kp*dz)*dzi
  ip=floor((xa(n)-xmin)*dxi+ilogbl)
  jp=floor((ya(n)-ymin)*dyi+jlogbl)
  kp=floor((za(n)-zmin)*dzi+klogbl)
  ab=((xmin-xa(n))+(ip-ilogbl+1)*dx)*dxi
  de=((ymin-ya(n))+(jp-jlogbl+1)*dy)*dyi
  gh=((zmin-za(n))+(kp-klogbl+1)*dz)*dzi
! this "if" statement slows things down, but I'm using it for debugging purposes
  if(ip<ilo .or. jp<jlo .or. kp<klo .or. ip>ihi .or. jp>jhi .or. kp>khi)then
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
if(ifail.ne.0)write(6,*)'(depose_rho_scalar) ifail=',ifail,' on process ',myrank

end subroutine
      
end module               
