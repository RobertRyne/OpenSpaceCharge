module fft_interface_mod


implicit none

contains


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!
!subroutine gsl_fft(cdata, n_size, i_sign)
!implicit none
!integer :: n_size
!integer(fgsl_size_t):: n
!complex(fgsl_double) :: cdata(n_size)
!integer(fgsl_int) :: status, i_sign
!
!n = n_size
!status = fgsl_fft_complex_radix2_transform(cdata, 1_fgsl_size_t, n, i_sign)
!
!end subroutine


subroutine fftw_ccfft3d(a,b,idir,n1,n2,n3,iskiptrans)
use, intrinsic :: iso_c_binding
use omp_lib
implicit none
include 'fftw3.f03'

integer :: idir(3)
type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), dimension(:,:,:) ::a, b
integer :: dir, fdir, n1, n2, n3, iskiptrans, n_threads

if (idir(1) == 1) then
  fdir = FFTW_BACKWARD
else
  fdir = FFTW_FORWARD
endif

print *, 'fftw_execute_dft...'
!$ n_threads = omp_get_max_threads()
!$ if (n_threads > 1) then
!$   print *, 'n_threads: ', n_threads
!$   call fftw_plan_with_nthreads(n_threads)
!$ endif

plan = fftw_plan_dft_3d(n3,n2,n1, a,b, fdir,FFTW_ESTIMATE)
call fftw_execute_dft(plan,a, b)
call fftw_destroy_plan(plan)

!$ if (n_threads > 1) call fftw_cleanup_threads()

print *, '...done'

end subroutine



subroutine my_fft(data, n, dir)
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
integer :: dir, fdir, n
type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), dimension(n) :: data, out


if (dir == 1) then
  fdir = FFTW_BACKWARD
else
  fdir = FFTW_FORWARD
endif

plan = fftw_plan_dft_1d(n, data,data, fdir,FFTW_ESTIMATE)
call fftw_execute_dft(plan, data, data)
call fftw_destroy_plan(plan)

end subroutine





end module
