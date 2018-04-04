module data_movement_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+     
subroutine localize(ptcls,lostflag,xmin,ymin,zmin,dx,dy,dz,nxloarg,nxhiarg,nyloarg,nyhiarg,nzloarg,nzhiarg, &
                    n1,nraysp,maxrayp,idecomp,npx,npy,npz,nprocs)
use mpi
implicit none

integer :: nxloarg,nxhiarg,nyloarg,nyhiarg,nzloarg,nzhiarg,n1,nraysp,maxrayp,idecomp,npx,npy,npz,nprocs
real(dp) :: xmin,ymin,zmin,dx,dy,dz
real(dp), dimension(n1,maxrayp) :: ptcls,gtmp
real(dp), dimension(maxrayp) :: lostflag
integer :: ilo,ihi,jlo,jhi,klo,khi
!     integer, dimension(0:nprocs-1) :: iloa,ihia,jloa,jhia,kloa,khia
integer :: n,j,k,l,ip,ifail
real(dp) :: dxi,dyi,dzi,ab,de,gh
integer, dimension(1:maxrayp) :: iownerptcl,itmp,irnk !fixed was mistakenly a real(dp)
real(dp), dimension(1:maxrayp) :: rtmp
integer :: nxtot,nytot,nztot
real(dp) :: recipnx,recipny,recipnz,recipnxnpx,recipnynpy,recipnznpz
!     integer, external :: iowner
!     integer, external :: lowner
integer :: mprocs,myrank,ierr
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
if(nprocs.ne.mprocs)then
  write(6,*)'(error in localize) mprocs.ne.nprocs; mprocs,nprocs=',mprocs,nprocs 
!       call myexit
  stop
endif
!
nxtot=nxhiarg-nxloarg+1
nytot=nyhiarg-nyloarg+1
nztot=nzhiarg-nzloarg+1
recipnx=1.d0/nxtot
recipny=1.d0/nytot
recipnz=1.d0/nztot
recipnxnpx=recipnx*npx
recipnynpy=recipny*npy
recipnznpz=recipnz*npz
!
!     do n=0,nprocs-1
!       call decompose(n,nprocs,nxloarg,nxhiarg,nyloarg,nyhiarg,nzloarg,nzhiarg,idecomp,npx,npy,npz,ilo,ihi,jlo,jhi,klo,khi)
!       iloa(n)=ilo; ihia(n)=ihi
!       jloa(n)=jlo; jhia(n)=jhi
!       kloa(n)=klo; khia(n)=khi
!     enddo


!now compute the owner of each particle:
dxi=1.d0/dx
dyi=1.d0/dy
dzi=1.d0/dz
do ip=1,nraysp
  if(lostflag(ip).ne.0.d0)cycle !!!!!!!! looks like a bug here; iownerptcl(ip) is undefined if listflag(ip).ne.0 !!!!!
  ab=(ptcls(1,ip)-xmin)*dxi !NOTE WELL: this essentially hardwires the coordinate indices
  de=(ptcls(3,ip)-ymin)*dyi
  gh=(ptcls(5,ip)-zmin)*dzi
  j=nint(ab)+nxloarg !floor(ab)+ilo
  k=nint(de)+nyloarg !floor(de)+jlo
  l=nint(gh)+nzloarg !floor(gh)+klo
!       iownerptcl(ip)=iowner(iloa,ihia,jloa,jhia,kloa,khia,nprocs,nxloarg,nxhiarg,nyloarg,nyhiarg,&
!    &                        nzloarg,nzhiarg,j,k,l)
  iownerptcl(ip)=lowner(recipnxnpx,recipnynpy,recipnznpz,npx,npy,j,k,l)
!debugging:
  if( iownerptcl(ip).lt.0 .or. iownerptcl(ip).gt.mprocs-1 )then
    write(6,*)'iownerptcl error: rank,owner=',myrank,iownerptcl(ip)
  endif
enddo


!!!!!now that I know who owns it, sort the array:
itmp(1:nraysp)=iownerptcl(1:nraysp)
rtmp(1:nraysp)=iownerptcl(1:nraysp)
gtmp(1:n1,1:nraysp)=ptcls(1:n1,1:nraysp)
call dmrgrnk(rtmp,irnk,nraysp,nraysp)
do k=1,nraysp
  iownerptcl(k)=itmp(irnk(k))
  ptcls(1:n1,k)=gtmp(1:n1,irnk(k))
enddo
!!!!!
if(mprocs.gt.1)then
!     if(myrank.eq.0)write(6,*)'calling newmove'
  call pmove(ptcls,iownerptcl,n1,maxrayp,nraysp,2*nraysp,ifail)
!       call newmove(ptcls,iownerptcl,n1,maxrayp,nraysp) !,2*nraysp,ifail)
!     if(myrank.eq.0)write(6,*)'back from newmove'
endif
      
end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
function lowner(recipnxnpx,recipnynpy,recipnznpz,npx,npy,i,j,k)
use mpi
implicit none

integer :: lowner
real(dp) :: recipnxnpx,recipnynpy,recipnznpz
integer :: npx,npy,i,j,k
integer :: inew,jnew,knew
integer :: myrank,ierr
inew=(i-1)*recipnxnpx
jnew=(j-1)*recipnynpy
knew=(k-1)*recipnznpz
lowner=inew + jnew*npx + knew*npx*npy

end function


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+     
subroutine newmove(y,iownerptcl,n1,n2,nraysp)
use mpi
implicit none


integer :: n1,n2,nraysp
real(dp), dimension(n1,n2) :: y
integer, dimension(1:n2) :: iownerptcl

integer :: myrank,mprocs,ierr
integer, allocatable :: sdisp(:),scounts(:),rdisp(:),rcounts(:)
integer :: i,k,n
integer :: ssize,rsize
real(dp), allocatable :: s1d(:),r1d(:)

call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!     if(myrank.eq.0)write(6,*)'hello from newmove'
allocate(sdisp(0:mprocs-1),scounts(0:mprocs-1),rdisp(0:mprocs-1),rcounts(0:mprocs-1))
! create the send count array:
scounts(0:mprocs-1)=0
do k=1,nraysp
  scounts(iownerptcl(k))=scounts(iownerptcl(k))+1
enddo
!     if(myrank.eq.0)write(6,*)'done with first loop'
! tell the other procs how much data is coming:
call MPI_alltoall(scounts,1,MPI_INTEGER,rcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
!     if(myrank.eq.0)write(6,*)'done with alltoall'
! create the send- and receive-displacement arrays:
sdisp(0)=0
do i=1,mprocs-1
  sdisp(i)=scounts(i-1)+sdisp(i-1)
enddo
rdisp(0)=0
do i=1,mprocs-1
  rdisp(i)=rcounts(i-1)+rdisp(i-1)
enddo
! send and receive counts:
ssize=sum(scounts)
rsize=sum(rcounts)
!     if(myrank.eq.0)write(6,*)'(PE0) ssize,rsize=',ssize,rsize

if(ssize.gt.n2)then
  write(6,*)'ssize overrun'
!       call myexit
  stop
endif

if(rsize.gt.n2)then
  write(6,*)'rsize overrun'
!       call myexit
  stop
endif
!
allocate(s1d(ssize))
allocate(r1d(rsize))
!     if(myrank.eq.0)write(6,*)'starting final do loop'
!     write(100+myrank,*)'PE, old nraysp:',myrank,nraysp
!     flush(100+myrank)
do k=1,n1
  s1d(1:ssize)=y(k,1:ssize)
  call MPI_alltoallv(s1d,scounts,sdisp,MPI_DOUBLE,r1d,rcounts,rdisp,MPI_DOUBLE,MPI_COMM_WORLD,ierr)
  y(k,1:rsize)=r1d(1:rsize)
enddo
nraysp=rsize
!     write(100+myrank,*)'PE,new nraysp:',myrank,nraysp
!     flush(100+myrank)

end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+     
! Given a 2D array y(n1,n2), and a 1D array iownerptcl(n2) that describes who owns each row,
! move the rows to the proc that owns them.
! This routine could be greatly improved, but it appears to work.  R.D. Ryne
subroutine pmove(y,iownerptcl,n1,n2,nraysp,nbufsize,ifail)

use mpi
implicit none

integer :: n1,n2,nraysp,nbufsize,ifail,n,nsbuf,nrecv,nshift,ntarget,nsender,nraysgbl
real(dp), dimension(n1,n2) :: y
integer, dimension(1:n2) :: iownerptcl
real(dp), dimension(n1,nbufsize) :: sbuf,rbuf
integer, dimension(2) :: msid
real(dp), dimension(n1,n2) :: ynew
integer :: nrayspnew
integer, dimension(MPI_STATUS_SIZE) :: mpistat
integer :: mprocs,myrank,ierr
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!
call MPI_REDUCE(nraysp,nraysgbl,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
ifail=0
nrayspnew=0

!#    do nshift=1,mprocs-1 !the nshift'th nearest neighbor in order from 1 to mprocs-1
do nshift=0,mprocs-1 !if storing in a "new" array, then have to start with nshift=0
  ntarget=mod(myrank+nshift,mprocs)       !send to processor ntarget
  nsender=mod(myrank-nshift+mprocs,mprocs)   !receive from processor nsender
! buffer the particles to be sent:
        nsbuf=0;n=1 !# in the send buffer; counter that goes through particles on this proc
     4  continue
        if(n.gt.nraysp)goto 5
        if(iownerptcl(n).eq.ntarget)then
          nsbuf=nsbuf+1
          if(nsbuf.gt.nbufsize)then
            write(6,*)'send buffer overflow';ifail=-1;return
          endif
          sbuf(1:n1,nsbuf)=y(1:n1,n)
          y(1:n1,n)=y(1:n1,nraysp) !swap last good particle with location n
          iownerptcl(n)=iownerptcl(nraysp)
          nraysp=nraysp-1
          goto 4
        endif
        n=n+1
        goto 4
    5   continue
! send and receive data:
  call MPI_IRECV(rbuf,n1*nbufsize,MPI_DOUBLE_PRECISION,nsender,nshift,MPI_COMM_WORLD,msid(1),ierr) !tag value is nshift
  call MPI_ISEND(sbuf,n1*nsbuf,MPI_DOUBLE_PRECISION,ntarget,nshift,MPI_COMM_WORLD,msid(2),ierr)
  call MPI_WAIT(msid(1),mpistat,ierr)
  call MPI_GET_COUNT(mpistat,MPI_DOUBLE_PRECISION,nrecv,ierr)
  call MPI_WAIT(msid(2),mpistat,ierr)
!! for some reason, under some conditions the following causes the code to freeze:
!!      call MPI_SEND(sbuf,n1*nsbuf,MPI_DOUBLE_PRECISION,ntarget,nshift,MPI_COMM_WORLD,ierr)
!!      call MPI_RECV(rbuf,n1*nbufsize,MPI_DOUBLE_PRECISION,nsender,nshift,MPI_COMM_WORLD,mpistat,ierr) !tag value is nshift
!!      call MPI_GET_COUNT(mpistat,MPI_DOUBLE_PRECISION,nrecv,ierr)
  nrecv=nrecv/n1
  if(nrecv.eq.0)cycle
! store the received data:
  if(nrayspnew+nrecv.gt.n2)then
    write(6,'(20h(not enough storage),1x,5(i9,1x))')myrank,nrayspnew,nrecv,nrayspnew+nrecv,n2
  endif
  ynew(1:n1,nrayspnew+1:nrayspnew+nrecv)=rbuf(1:n1,1:nrecv)
  nrayspnew=nrayspnew+nrecv
enddo


nraysp=nrayspnew
y(1:n1,1:nraysp)=ynew(1:n1,1:nraysp)
call MPI_REDUCE(nraysp,nraysgbl,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

end subroutine 
!
!



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+     
subroutine dmrgrnk(XDONT, IRNGT, sizexdont,sizeirngt)
implicit none
integer :: sizexdont,sizeirngt
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
real(dp), Dimension (*), Intent (In) :: XDONT     !changed (:) to (*)
Integer, Dimension (*), Intent (Out) :: IRNGT   !changed (:) to (*)
! __________________________________________________________
real(dp) :: XVALA, XVALB
!
!!!!!!Integer, Dimension (SIZE(IRNGT)) :: JWRKT
Integer, Dimension (sizeirngt) :: JWRKT
Integer :: LMTNA, LMTNC, IRNG1, IRNG2
Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
!     write(6,*)'size(xdont)=',size(xdont)
!     write(6,*)'size(irngt)=',size(irngt)
!!!!!!NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
NVAL = Min (sizexdont, sizeirngt)
Select Case (NVAL)
Case (:0)
   Return
Case (1)
   IRNGT (1) = 1
   Return
Case Default
   Continue
End Select

!
!  Fill-in the index array, creating ordered couples
!
Do IIND = 2, NVAL, 2
   If (XDONT(IIND-1) <= XDONT(IIND)) Then
      IRNGT (IIND-1) = IIND - 1
      IRNGT (IIND) = IIND
   Else
      IRNGT (IIND-1) = IIND
      IRNGT (IIND) = IIND - 1
   End If
End Do
If (Modulo(NVAL, 2) /= 0) Then
   IRNGT (NVAL) = NVAL
End If

!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
LMTNA = 2
LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
Do
   If (NVAL <= 2) Exit
  !
  !   Loop on merges of A and B into C
  !
  Do IWRKD = 0, NVAL - 1, 4
    If ((IWRKD+4) > NVAL) Then
      If ((IWRKD+2) >= NVAL) Exit
      !
      !   1 2 3
      !
      If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
      !
      !   1 3 2
      !
      If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
        IRNG2 = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNG2
      !
      !   3 1 2
      !
      Else
        IRNG1 = IRNGT (IWRKD+1)
        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
        IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNG1
      End If
      Exit
    End If
    !
    !   1 2 3 4
    !
    If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
    !
    !   1 3 x x
    !
    If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
       IRNG2 = IRNGT (IWRKD+2)
       IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
       If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
         !   1 3 2 4
         IRNGT (IWRKD+3) = IRNG2
       Else
         !   1 3 4 2
         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
         IRNGT (IWRKD+4) = IRNG2
       End If
     !
     !   3 x x x
     !
    Else
      IRNG1 = IRNGT (IWRKD+1)
      IRNG2 = IRNGT (IWRKD+2)
      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
      If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
        IRNGT (IWRKD+2) = IRNG1
        If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
          !   3 1 2 4
          IRNGT (IWRKD+3) = IRNG2
        Else
        !   3 1 4 2
          IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
          IRNGT (IWRKD+4) = IRNG2
        End If
      Else
        !   3 4 1 2
        IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
        IRNGT (IWRKD+3) = IRNG1
        IRNGT (IWRKD+4) = IRNG2
      End If
    End If
  End Do
!
!  The Cs become As and Bs
!
  LMTNA = 4
  Exit
End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
Do
   If (LMTNA >= NVAL) Exit
   IWRKF = 0
   LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
  Do
    IWRK = IWRKF
    IWRKD = IWRKF + 1
    JINDA = IWRKF + LMTNA
    IWRKF = IWRKF + LMTNC
    If (IWRKF >= NVAL) Then
       If (JINDA >= NVAL) Exit
       IWRKF = NVAL
    End If
    IINDA = 1
    IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
    JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
    XVALA = XDONT (JWRKT(IINDA))
    XVALB = XDONT (IRNGT(IINDB))
!
    Do
      IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
      If (XVALA > XVALB) Then
         IRNGT (IWRK) = IRNGT (IINDB)
         IINDB = IINDB + 1
         If (IINDB > IWRKF) Then
          !  Only A still with unprocessed values
            IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
            Exit
         End If
         XVALB = XDONT (IRNGT(IINDB))
      Else
         IRNGT (IWRK) = JWRKT (IINDA)
         IINDA = IINDA + 1
         If (IINDA > LMTNA) Exit! Only B still with unprocessed values
         XVALA = XDONT (JWRKT(IINDA))
      End If
!
    End Do
  End Do
  !
  !  The Cs become As and Bs
  !
  LMTNA = 2 * LMTNA
End Do
!

!
end subroutine dmrgrnk

end module
