module multitrack
      use rays
      use lieaparam, only : monoms
      implicit none
      real*8, dimension(:), allocatable :: tlistin,ptlistin !initial t,pt
      real*8, dimension(:,:), allocatable :: rho56  !t-pt density (future use)
      integer,dimension(:),allocatable::inear,jnear !indices of nearest init t,pt
      real*8,dimension(:),allocatable::tdelt,ptdelt !distance to nearest init t,pt
      real*8, dimension(:,:,:,:), allocatable :: tmhlist !linear maps
      real*8, dimension(:,:,:), allocatable :: hlist     !nonlinear maps
      real*8, dimension(:,:,:,:),allocatable::pbhlist !Taylor tracking info
      real*8, dimension(:,:), allocatable :: tlistfin,ptlistfin !final t,pt
!if multitrac.ne.0, then use multiple maps when tracking
!cover the longitudinal phase space on a grid of size (imaps)x(jmaps)
!refentry5,refentry6 are values for the (single) ref particle at entrance
      integer :: multitrac,imaps,jmaps
      real*8 :: refentry5,refentry6
save

contains

      subroutine new_multiarrays
! create arrays "on the fly" based on imaps, jmaps, nraysp
      allocate(tlistin(imaps),ptlistin(jmaps),rho56(imaps,jmaps),          &
     &inear(nraysp),jnear(nraysp),tdelt(nraysp),ptdelt(nraysp),            &
     &tmhlist(6,6,imaps,jmaps),hlist(monoms,imaps,jmaps),                  &
     &tlistfin(imaps,jmaps),ptlistfin(imaps,jmaps))
      return
      end subroutine new_multiarrays

      subroutine del_multiarrays
! destroy arrays
      deallocate(tlistin,ptlistin,inear,jnear,tdelt,ptdelt,                &
     &           tmhlist,hlist)
      return
      end subroutine del_multiarrays

      subroutine multibrkts
! for nonlinear Taylor tracking, need to store the result of calling brkts
      include 'pbkh.inc'
      real*8 pbh
      integer i,j
      if(allocated(pbhlist))deallocate(pbhlist)
      if(idproc.eq.0)write(6,*)'(multibrkts): allocating pbhlist array'
      allocate(pbhlist(monoms,12,imaps,jmaps))
      do i=1,imaps
      do j=1,jmaps
        call brkts(hlist(1,i,j))
        pbhlist(:,:,i,j)=pbh(:,:)
      enddo
      enddo
      return
      end subroutine multibrkts

      subroutine multicanx
! for nonlinear symplectic tracking, need to store the result of calling canx
      if(idproc.eq.0)                                                     &
     &write(6,*)'error: multitrack symp tracking not yet implemented'
      call myexit
      end subroutine multicanx

end module multitrack
