module test_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains




      subroutine gendist(ptcls,n1,nraysp,sigmat,gaussiancutoff,disttype,iseed)
      implicit none
      integer :: n1,nraysp,disttype,iseed
   !   real(8), dimension(nraysp,n1) :: ptcls
      real(8), dimension(nraysp,n1) :: ptcls
      real(8), dimension(6,6) :: sigmat
      real(8) :: gaussiancutoff
      integer :: n
      real(8) :: r1,r2,r3,r4,r5,r6,arg
      real(8), dimension(6) :: cent
      
      print *, 'gendist'
      
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
        endif
        if(disttype.eq.1)then
   11     continue
          call normdv(r1,r2)
          call normdv(r3,r4)
          call normdv(r5,r6)
          if(r1*r1+r3*r3+r5*r5.gt.gaussiancutoff**2)goto 11
        endif
        
       ! print *, 'gendist 1, n = ', n, r1*sqrt(sigmat(1,1))
! this version assumes that sigmat is 2x2 block diagonal
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
      return
      
       print *, 'gendist end'
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
      implicit none
      integer :: inpseed
      integer :: iseedsize,n
      integer, dimension(:), allocatable :: iputarray
      integer :: idproc=0
! determine the seed size:
      call random_seed(size=iseedsize)
!     if(idproc.eq.0)write(6,*)'iseedsize=',iseedsize
      allocate(iputarray(iseedsize))
      iputarray(1)=inpseed
! initialize the remaining elements of iputarray:
      if(iseedsize.gt.1)then
        do n=2,iseedsize
          iputarray(n)=iputarray(n-1)+1
        enddo
      endif
! initialize random_number with the seed array:
      call random_seed(put=iputarray)
      return
      end subroutine initrandom
      
      
      
      !diagnostics routines
      subroutine prntall(nstep,n1,nraysp,nx,ny,nz,ptcls,hx,hy,hz,ex,ey,ez,t,dx,dy,dz,xmin,ymin,zmin)
      implicit none
      integer :: nstep,n1,nraysp,nx,ny,nz
      real(8) :: t,dx,dy,dz,xmin,ymin,zmin
      real(8), dimension(nx,ny,nz) :: hx,hy,hz,ex,ey,ez
      real(8), dimension(nraysp,n1) :: ptcls
      character(32) :: pname,fname,fxname,fyname,fzname
      integer, parameter :: ndigits=4         !digits in filename suffix corr. to step#
      integer, parameter :: punit=2,funit=3,fxunit=7,fyunit=8,fzunit=9   !file unit numbers
      integer, parameter :: maxptclprnt=10000 !max# of particles printed
      
      
      pname='ptcls' !prefix for particle filename
      fname='fields' !prefix for field filename
      fxname='xline' !field vs x
      fyname='yline' !field vs y
      fzname='zline' !field vs z
!print particles:
      call openfile(pname,punit,nstep,ndigits)
      
      
      call pprnt(ptcls,n1,nraysp,maxptclprnt,punit)
      close(punit)
!print fields:
      call openfile(fname,funit,nstep,ndigits)
      call openfile(fxname,fxunit,nstep,ndigits)
      call openfile(fyname,fyunit,nstep,ndigits)
      call openfile(fzname,fzunit,nstep,ndigits)
      
      
      call heprnt(hx,hy,hz,ex,ey,ez,nx,ny,nz,funit,fxunit,fyunit,fzunit,dx,dy,dz,xmin,ymin,zmin)
      close(funit); close(fxunit); close(fyunit); close(fzunit)
      return
      
      contains

!
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
      end
      
      subroutine heprnt(hx,hy,hz,ex,ey,ez,nx,ny,nz,funit,fxunit,fyunit,fzunit,dx,dy,dz,xmin,ymin,zmin)
      implicit none
      integer :: nx,ny,nz,funit,fxunit,fyunit,fzunit
      real(8) :: dx,dy,dz,xmin,ymin,zmin
      real(8), dimension(nx,ny,nz) :: hx,hy,hz,ex,ey,ez
      integer :: i,j,k
! fields on the entire grid:
      do k=1,nz
      do j=1,ny
      do i=1,nx
        write(funit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    hx(i,j,k),hy(i,j,k),hz(i,j,k),ex(i,j,k),ey(i,j,k),ez(i,j,k),i,j,k
      enddo
      enddo
      enddo
!    
! fields vs x along a line: 
      j=ny/2
      k=nz/2
      do i=1,nx
        write(fxunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    hx(i,j,k),hy(i,j,k),hz(i,j,k),ex(i,j,k),ey(i,j,k),ez(i,j,k),i,j,k
      enddo
!
! fields vs y along a line:
      i=nx/2
      k=nz/2
      do j=1,ny
        write(fyunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    hx(i,j,k),hy(i,j,k),hz(i,j,k),ex(i,j,k),ey(i,j,k),ez(i,j,k),i,j,k
      enddo
!
! fields vs z along a line:
      i=nx/2
      j=ny/2
      do k=1,nz
        write(fzunit,'(9(1pg12.5,1x),3i5)')xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,&
     &                                    hx(i,j,k),hy(i,j,k),hz(i,j,k),ex(i,j,k),ey(i,j,k),ez(i,j,k),i,j,k
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
      
end module               
