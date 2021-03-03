module beamdata
      implicit none
! beam parameters and scaling parameters
      real*8 :: brho,gamma,gamm1,beta,achg,pmass,bfreq,bcurr,c
      real*8 :: sl,p0sc,ts,omegascl,freqscl
      logical :: lflagmagu=.true.,lflagdynu=.false.
      save
end module beamdata

module acceldata
      implicit none

! lmnlbl(i) = name of ith element (numbered as they occur in #menu)
! p(k,i)    = kth parameter of ith element
! nt1(i)    = group index of the elements type (see /monics/)
! nt2(i)    = type index  .....
! mppr(i)    = menu parameters (real data)    for ith element
! note: for compatibility w/ existing MaryLie, I call this "mpp" not "mppr"
! mppc(i)    = menu parameters (char*16 data) for ith element
      character(len=16), dimension(:), allocatable :: lmnlbl,cmenu
! if there are, on average,
! 6 real parameters/element and 1 character parameter/element,
! then this is enough memory to store the parameters of mnumax elements
      real*8, dimension(:), allocatable :: pmenu
! atarray added by RDR 9/16/2004:
      real*8, dimension(:), allocatable :: atarray
!
      integer, dimension(:), allocatable :: mpp,mppc,nt1,nt2

! mnumax = maximum number of elements in #menu (default=5000)
! itmno  = maximum number of items (default=5000)
! itmnln = maximum number of elements per item (default=100)
! mlabor = maximum number of tasks in #labor (default=1000)
! maxcmt = maximum number of comment lines (default=250)
! nconmax = maximum number of defined constants (default=1000)
      integer :: mnumax=5000,itmno=5000,itmln=100
      integer :: mlabor=1000,maxcmt=250,nconmax=1000
!
      integer :: inmenu

! na    = actual number of elements
! nb    = actual number of items
! noble = actual number of tasks
! npcom = actual number of comment lines
! mpprpoi  = menu parameter pointer to "real" (i.e. numeric) data
! mppcpoi  = menu parameter pointer to character*16 data
      integer :: na=0,nb=0,noble=0,npcom=0,mpprpoi=0,mppcpoi=0


! lines,lumps and loops
! ilbl(i)   = name of ith item
! icon(k,i) = name of kth element/item in ith item
! irep(k,i) = repetition factor of icon(k,i)
! ityp(k,i) = type code = 2,3,4 if item=line,lump,loop resp.
! ilen(i)   = number of elements/items in ith item (=max k for ith item)
! lmade(i)  = flag, telling whether map of lump is in buffer and where.
!             preset to zero
      character(len=16), dimension(:), allocatable :: ilbl
      character(len=16), dimension(:,:), allocatable :: icon
      integer, dimension(:), allocatable :: ityp,ilen,lmade
      integer, dimension(:,:), allocatable :: irep

! latt(i) = the name of the ith task
! num(i) = the repetition factor of this task
      character(len=16), dimension(:), allocatable :: latt
      integer, dimension(:), allocatable :: num

      character(len=80), dimension(:), allocatable :: mline

!ryne 7/7/2002
!ryne this common block is used to store strings that have constant
!ryne values assigned to them.
!ryne nconst is equal to the number that have been assigned
!ryne This info should really be included in elments.inc, but I am
!ryne keeping it separate to avoid changes to the existing MaryLie
      character(len=16), dimension(:), allocatable :: constr
      real*8, dimension(:), allocatable :: conval
      integer :: nconst
!
!11/29/02 autoslicing info:
! defaults:
      character*16 :: slicetype='slices',sliceprecedence='local'
      real*8 :: slicevalue=1.0

      save

contains
      subroutine new_acceldata
      allocate(lmnlbl(mnumax))
      allocate(pmenu(6*mnumax),cmenu(1*mnumax),                         &
     &              mpp(mnumax),mppc(mnumax),nt1(mnumax),nt2(mnumax))
      allocate(atarray(1*mnumax))
!     write (6,*) 'pmenu size=',size(pmenu)
!     write (6,*) 'pmenu start,end=',loc(pmenu(1)),loc(pmenu(6*mnumax))

      allocate(ilbl(itmno),icon(itmln,itmno))
      allocate(irep(itmln,itmno),ityp(itmno),ilen(itmno),lmade(itmno))
      allocate(latt(mlabor),num(mlabor))
      allocate(mline(maxcmt))
      allocate(constr(nconmax),conval(nconmax))
        mpp(:)=0
        mppc(:)=0
        cmenu(:)=' '
        lmade(:)=0
        mline(:)=' '
        latt(:)=' '
        num(:)=1
        atarray(:)=-99999.d0
      end subroutine new_acceldata

      subroutine del_acceldata
      deallocate(lmnlbl)
      deallocate(pmenu,cmenu,mpp,mppc,nt1,nt2)
      deallocate(atarray)
      deallocate(ilbl,icon)
      deallocate(irep,ityp,ilen,lmade)
      deallocate(latt,num)
      deallocate(mline)
      deallocate(constr,conval)
      end subroutine del_acceldata


end module acceldata


module rays
      use parallel
      implicit none
! eventually make leading dimension "idimp" instead of 6
! zblock(k,i) kth coordinate of ray index i
! tblock(k,i) temporary array; might move out of this module later.
! zi (k)      initial coordinates
! zf (k)      final coordinates
! nrays       number of rays
! istat(i)    iturn, in which ray index i was lost, 0 otherwise
! ihist(1,j)  iturn, in which the jth lost particle was lost
! ihist(2,j)  index (i) of jth lost particle
! iturn       turn index (not used in cqlate)
! nlost       number of rays lost so far
      real*8, dimension(:,:), allocatable :: zblock,tblock
      real*8, dimension(:,:), allocatable :: pbh6t,uvect,tvect
      integer, dimension(:), allocatable :: istat
      integer, dimension(:,:), allocatable :: ihist
      real*8, dimension(6) :: zi,zf
! maxray is the initial size of the global particle array (default=30000)
      integer :: maxray=30000,nrays=0,nlost=0,iturn=0
! values per processor:
      integer, parameter :: nchunk=100
      integer :: maxrayp,nraysp=0
save

contains

      subroutine new_particledata
!cryne May 19, 2006      maxrayp=(maxray-1)/nvp + 1
      maxrayp=maxray/nvp + 1
!cryne May 19, 2006 commented out the following until particle mgr is installed:
!     maxrayp=2*maxrayp
!
      allocate(zblock(6,maxrayp),tblock(6,maxrayp),                        &
     &         pbh6t(6,923),uvect(923,nchunk),tvect(923,nchunk))
!     if(nvp.gt.1)allocate(tblock(6,maxrayp))
      allocate(istat(maxrayp),ihist(2,maxrayp))
      end subroutine new_particledata

      subroutine del_particledata
!     print *, 'zblock=',allocated(zblock),size(zblock)
!     print *, 'istat=',allocated(istat),size(istat)
!     print *, 'ihist=',allocated(ihist),size(ihist)
      integer status
      deallocate(zblock,tblock,pbh6t,uvect,tvect,stat=status)
      if(status.ne.0) print *,'del_particledata: deallocate1 returned ',status
!     if(allocated(tblock))deallocate(tblock)
      deallocate(istat,ihist,stat=status)
      if(status.ne.0) print *,'del_particledata: deallocate2 returned ',status
      end subroutine del_particledata
end module rays
