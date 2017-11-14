module fft_mod



contains


! -----------------------------------------------------------------------
! perform FFTs

      subroutine fft_perform(data_in,data_out,idirection,nx,ny,nz,      &
     &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     &     ipermute,iscale,time)
      USE mpi !include "mpif.h"
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, parameter :: dp = REAL64
!     COMPLEX_DATA data_in(*)
      complex(dp) data_in(*)
      complex(dp) data_out(*)
      integer idirection,nx,ny,nz
      integer in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
      integer ipermute,iscale
      real(dp) time,time1,time2
      real(dp) plan1,plan2
      integer nbuf,ierror


!     if(idirection.ne.1 .and. idirection.ne.-1)then
!     if(myrank.eq.0)write(6,*)'(fft_perform) bad idirection argument'
!     call myexit
!     endif

      if(idirection.eq.1)then
      call fft_3d_create_plan_interface(mpi_comm_world,nx,ny,nz,        &
     &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     &     iscale,ipermute,nbuf,plan1)
      endif

      if(idirection.eq.-1)then
      if (ipermute.eq.0) then
        call fft_3d_create_plan_interface(mpi_comm_world,nx,ny,nz,      &
     &       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,           &
     &       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                 &
     &       iscale,0,nbuf,plan2)
      else if (ipermute.eq.1) then
        call fft_3d_create_plan_interface(mpi_comm_world,ny,nz,nx,      &
     &       out_jlo,out_jhi,out_klo,out_khi,out_ilo,out_ihi,           &
     &       in_jlo,in_jhi,in_klo,in_khi,in_ilo,in_ihi,                 &
     &       iscale,2,nbuf,plan2)   !typo??? ,2, -> ,1, ???
      else if (ipermute.eq.2) then
        call fft_3d_create_plan_interface(mpi_comm_world,nz,nx,ny,      &
     &       out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,           &
     &       in_klo,in_khi,in_ilo,in_ihi,in_jlo,in_jhi,                 &
     &       iscale,1,nbuf,plan2)   !typo??? ,1, -> ,2, ???
      endif
!     if(myrank.eq.0)write(6,*)
!    $ 'check calls to create_plan_interface re ipermute'
      endif

! LAMMPS-like calls

!      call fft_3d_create_plan_interface(mpi_comm_world,nx,ny,nz,
!     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
!     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
!     $     1,0,nbuf,plan1)
!
!      call fft_3d_create_plan_interface(mpi_comm_world,nx,ny,nz,
!     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
!     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
!     $     1,0,nbuf,plan2)

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      
!ryne do i = 1,iteration

!ryne   if (mod(i,2).eq.1) then
!ryne     call fft_3d_interface(data,data,1,plan1)
!ryne   else
!ryne     call fft_3d_interface(data,data,-1,plan2)
!ryne   endif

! LAMMPS-like calls

!        call fft_3d(data,data,1,plan1)
!        call fft_3d(data,data,-1,plan2)
!        call fft_3d(data,data,-1,plan2)
!        call fft_3d(data,data,-1,plan2)

!ryne enddo
      if(idirection.eq.1)then
        call fft_3d_interface(data_in,data_out,idirection,plan1)
      endif
      if(idirection.eq.-1)then
        call fft_3d_interface(data_in,data_out,idirection,plan2)
      endif
      
!ryne write(6,*)'exiting fft_perform'
!ryne call myexit
      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      
      time = time2 - time1

      if(idirection.eq.1)then
      call fft_3d_destroy_plan_interface(plan1)
      endif
      if(idirection.eq.-1)then
      call fft_3d_destroy_plan_interface(plan2)
      endif
      
      return
      end

! -----------------------------------------------------------------------
! perform remaps

      subroutine xpose_perform(data,buf,iteration,                      &
     &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     &     nqty,ipermute,memory,time)
      USE mpi !include "mpif.h"
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, parameter :: dp = REAL64
!     REAL_DATA data(*),buf(*)
      real(dp) data(*),buf(*)
      real(dp) time,time1,time2
      real(dp) plan1,plan2
      integer iteration,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
      integer nqty,ipermute,iscale,memory
      integer PREC,ierror,i

      call remap_3d_create_plan_interface(mpi_comm_world,               &
     &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     &     nqty,ipermute,memory,PREC,plan1)

      if (ipermute.eq.0) then
        call remap_3d_create_plan_interface(mpi_comm_world,             &
     &       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,           &
     &       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                 &
     &       nqty,0,memory,PREC,plan2)
      else if (ipermute.eq.1) then
        call remap_3d_create_plan_interface(mpi_comm_world,             &
     &       out_jlo,out_jhi,out_klo,out_khi,out_ilo,out_ihi,           &
     &       in_jlo,in_jhi,in_klo,in_khi,in_ilo,in_ihi,                 &
     &       nqty,2,memory,PREC,plan2)
      else if (ipermute.eq.2) then
        call remap_3d_create_plan_interface(mpi_comm_world,             &
     &       out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,           &
     &       in_klo,in_khi,in_ilo,in_ihi,in_jlo,in_jhi,                 &
     &       nqty,1,memory,PREC,plan2)
      endif

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      
      do i = 1,iteration
        if (mod(i,2).eq.1) then
          call remap_3d_interface(data,data,buf,plan1)
        else
          call remap_3d_interface(data,data,buf,plan2)
        endif
      enddo
      
      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      
      time = time2 - time1

      call remap_3d_destroy_plan_interface(plan1)
      call remap_3d_destroy_plan_interface(plan2)
      
      return
      end
! -----------------------------------------------------------------------


end module
