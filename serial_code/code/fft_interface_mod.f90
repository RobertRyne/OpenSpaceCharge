module fft_interface_mod

use fast_fourier_am

implicit none
! Fortran 2008
integer, parameter, private :: dp = REAL64

contains


#ifdef USE_FFTW
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


subroutine ccfft3d(a,b,idir,n1,n2,n3,iskiptrans)
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


#else



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine ccfft3d(a,b,idir,n1,n2,n3,iskiptrans)
implicit none
integer, dimension(3) :: idir
complex(dp), dimension(:,:,:) :: a,b
complex(dp), allocatable :: tmp1(:,:,:), tmp2(:,:,:)
integer :: n1,n2,n3
integer :: iskiptrans
integer :: ileft,iright,i,j,k
b(:,:,:)=a(:,:,:)
allocate(tmp1(n2,n1,n3))
allocate(tmp2(n3,n2,n1))
if(idir(1).eq.0.and.idir(2).eq.0.and.idir(3).eq.0)return
call mccfft1d(b,n2*n3,n1,idir(1))
forall(k=1:n3, j=1:n2, i=1:n1)
! there's some problem with the commented out statements
!       iright=(k-1)*n1*n2+(j-1)*n2+i
!       ileft =(k-1)*n2*n1+(i-1)*n1+j
!       tmp(ileft)=b(iright)
  tmp1(j,i,k)=b(i,j,k)
end forall

call mccfft1d(tmp1,n3*n1,n2,idir(2))
forall(k=1:n3, j=1:n2, i=1:n1)
!       iright=(k-1)*n2*n1+(i-1)*n1+j
!       ileft =(i-1)*n2*n3+(j-1)*n3+k
!       b(ileft)=tmp(iright)
  tmp2(k,j,i)=tmp1(j,i,k)
end forall

call mccfft1d(tmp2,n1*n2,n3,idir(3))
if(iskiptrans.eq.1)return
forall(k=1:n3, j=1:n2, i=1:n1)
!       ileft =(k-1)*n1*n2+(j-1)*n2+i
!       iright=(i-1)*n2*n3+(j-1)*n3+k
!       b(ileft)=tmp(iright)
        b(i,j,k)=tmp2(k,j,i)
end forall

end subroutine 
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!-------------------------------ccfftnr-----------------------------------------
!+
subroutine mccfft1d(a,ntot,lenfft,idir)
implicit none
complex(dp), dimension(*) :: a
integer :: ntot,lenfft,idir
integer :: n,ierr
do n=1,ntot*lenfft,lenfft
  call ccfftam(a(n),lenfft,idir,ierr)  ! Alan Miller version of FFT package
  !call my_fft(a(n),lenfft,idir) ! FFTW
  !ierr =  0
  
  if(ierr.ne.0)then
    write(6,*)'Error return from FFT package due to transform length.'
    write(6,*)'Try increasing the padding by 1 in each dimension and re-run'
    stop
  endif
! call ccfftnr(a(n),lenfft,idir) !numerical recipes power-of-2 routine
! call gsl_fft(a(n),lenfft,idir)
enddo

end



#endif



end module
