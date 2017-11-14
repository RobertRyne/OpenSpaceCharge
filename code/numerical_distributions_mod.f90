module numerical_distributions_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine gendist(ptcls,n1,maxrayp,nraysp,sigmat,gaussiancutoff,disttype,iseed)
      use mpi
      implicit none
      integer :: n1,nraysp,maxrayp,disttype,iseed
      real(dp), dimension(n1,maxrayp) :: ptcls
      real(dp), dimension(6,6) :: sigmat
      real(dp), dimension(6,6) :: sq
      real(dp), dimension(6) :: vec,newvec
      real(dp) :: gaussiancutoff
      integer :: n
      real(dp) :: r1,r2,r3,r4,r5,r6,arg
      real(dp), dimension(6) :: cent
      real(dp) :: epsx2,epsy2,epsz2
      integer :: myrank,mpierr,ierr

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      call initrandom(iseed)

      epsx2=sigmat(1,1)*sigmat(2,2)-sigmat(1,2)**2
      epsy2=sigmat(3,3)*sigmat(4,4)-sigmat(3,4)**2
      epsz2=sigmat(5,5)*sigmat(6,6)-sigmat(5,6)**2
      if(myrank.eq.0)write(6,"('epsx2,epsy2,epsz2=',3(1pe14.7,1x))")epsx2,epsy2,epsz2
      if(epsx2.lt.0.d0.or.epsy2.lt.0.d0.or.epsz2.lt.0.d0)stop

      call cholesky(sigmat,sq,6,6,ierr)
      if(ierr.ne.0)then
        if(myrank.eq.0)write(6,*)'error return from cholesky; ierr=',ierr
        stop
      endif

      do n=1,nraysp
        if(disttype.eq.0)then
   10     continue
          call random_number(r1); r1=2.d0*r1-1.d0
          call random_number(r2); r2=2.d0*r2-1.d0
          call random_number(r3); r3=2.d0*r3-1.d0
          call random_number(r4); r4=2.d0*r4-1.d0
          call random_number(r5); r5=2.d0*r5-1.d0
          call random_number(r6); r6=2.d0*r6-1.d0
          if(r1*r1+r3*r3+r5*r5.gt.1.d0)goto 10
!         if(r1*r1+r2*r2+r3*r3+r4*r4+r5*r5+r6*r6.gt.1.d0)goto 10
        endif
        if(disttype.eq.1)then
   11     continue
          call normdv(r1,r2); call normdv(r3,r4); call normdv(r5,r6)
          if(r1*r1+r3*r3+r5*r5.gt.gaussiancutoff**2)goto 11
!         if(r1*r1+r2*r2+r3*r3+r4*r4+r5*r5+r6*r6.gt.1.d0)goto 10
        endif

        vec(1)=r1;vec(2)=r2;vec(3)=r3;vec(4)=r4;vec(5)=r5;vec(6)=r6
        newvec=matmul(sq,vec)
        ptcls(1:6,n)=newvec(1:6)
      enddo

      ptcls(1:6,1)=0.d0 !zero out the first particle to see how the origin gets kicked around by space charge
      return
      end subroutine gendist
!
      subroutine normdv(d1,d2)
! routine to generate 2 normal deviates
      implicit none
      real(dp) :: u1,u2
      real(dp) twopi,d1,d2
      twopi=4.0d0*asin(1.0d0)
      call random_number(u1)
      call random_number(u2)
      d1=sqrt(-2.0d0*log(u1))*cos(twopi*u2)
      d2=sqrt(-2.0d0*log(u1))*sin(twopi*u2)
      return
      end subroutine normdv
!
      subroutine initrandom(inpseed)
      use mpi
      implicit none
      integer :: inpseed
      integer :: iseedsize,n,i
      integer, dimension(:), allocatable :: iputarray
      integer :: myrank,ierr
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
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
!
      subroutine serial_initrandom(inpseed)
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
      end subroutine serial_initrandom

!not sure where I got cholesky from
      subroutine cholesky(a,b,n,np,ierr)
      implicit none
      integer :: n,np,ierr
      real(dp), dimension(np,np) :: a,b
      real(dp), dimension(n) :: p
      integer :: i,j,k
      real(dp) :: sumx
      ierr=-1
      if(n.gt.np)return
      b(:,:)=a(:,:)
      do i=1,n
        do j=i,n
          sumx=b(i,j)
          do k=i-1,1,-1
            sumx=sumx-b(i,k)*b(j,k)
          enddo
          if(i.eq.j)then
            if(sumx.le.0.d0)ierr=-i
            if(sumx.le.0.d0)return
            p(i)=sqrt(sumx)
          else
            b(j,i)=sumx/p(i)
          endif
        enddo
      enddo
      do i=1,n
        do j=i,n
          b(i,j)=0.d0
          if(i.eq.j)b(i,i)=p(i)
        enddo
      enddo
      ierr=0
      return
      end

end module
