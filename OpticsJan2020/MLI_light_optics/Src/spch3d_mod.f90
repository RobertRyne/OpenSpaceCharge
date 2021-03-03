module spchdata
      implicit none
!!!   real*8, dimension(:,:,:), allocatable :: rho
      complex*16, dimension(:,:,:), allocatable :: rho2
      complex*16, dimension(:,:,:), allocatable :: grnxtr,rho2xtr
!3/21/03 added three more arrays:
      complex*16, dimension(:,:,:), allocatable::xgrnxtr,ygrnxtr,zgrnxtr
save

contains

      subroutine new_spchdata(nx,ny,nz,n1,n2,n3a,isolve)
      implicit none
      integer :: nx,ny,nz,n1,n2,n3a,isolve
      integer :: ierror
!!!   allocate(rho(nx,ny,nz))
      if( isolve .ge. 10 )return  !ANAG solver doesnt need these arrays
      allocate(rho2(n1,n2,n3a),stat=ierror)
      if( ierror .NE. 0 )then
        print*,'error: new_spchdata: allocate(rho2) failed with code ',ierror
        stop 'NEWSPCH'
      endif
      allocate(grnxtr(n3a,n2,n1),rho2xtr(n3a,n2,n1),stat=ierror)
      if( ierror .NE. 0 )then
        print*,'error: new_spchdata: allocate({grn,rho2}xtr) failed with code ',ierror
        stop 'NEWSPCH'
      endif
!3/21/03 three more arrays:
      allocate(xgrnxtr(n3a,n2,n1),ygrnxtr(n3a,n2,n1),zgrnxtr(n3a,n2,n1),stat=ierror)
      if( ierror .NE. 0 )then
        print*,'error: new_spchdata: allocate({xyz}grnxtr) failed with code ',ierror
        stop 'NEWSPCH'
      endif
      end subroutine new_spchdata

      subroutine del_spchdata
!!!   deallocate(rho)
      deallocate(rho2)
      deallocate(grnxtr,rho2xtr)
!3/21/03 three more arrays:
      deallocate(xgrnxtr,ygrnxtr,zgrnxtr)
      end subroutine del_spchdata
end module spchdata
