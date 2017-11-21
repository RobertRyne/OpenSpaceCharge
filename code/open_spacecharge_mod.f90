module open_spacecharge_mod

use, intrinsic :: iso_fortran_env

use fft_mod
use decomposition_mod

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128





type field_struct
  real(dp) :: E(3) = 0        ! electric field
  real(dp) :: B(3) = 0        ! magnetic field
end type

type mesh3d_struct
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma in the +z (3rd) direction
  real(dp), allocatable, dimension(:,:,:) :: rho ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: phi ! electric potential grid
  type(field_struct),  allocatable, dimension(:,:,:) :: field ! field grid
end type



contains

subroutine init_mesh3d(mesh3d, domain) 
type(mesh3d_struct) :: mesh3d
type (domain_decomposition_struct) :: domain



end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine getfields(rho,gam0,dx,dy,dz,phi,ex,ey,ez,hx,hy,hz,nx,ny,nz,idecomp,npx,npy,npz,igfflag,idirectfieldcalc, &
                ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl)
use mpi
implicit none
integer, intent(in) :: nx,ny,nz,idecomp,npx,npy,npz,igfflag,idirectfieldcalc
integer, intent(in) :: ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl
real(dp), intent(in) :: gam0,dx,dy,dz
!!!   real(dp), intent(in), dimension(nx,ny,nz) :: rho
!!!   real(dp), intent(out), dimension(nx,ny,nz) :: phi,hx,hy,hz,ex,ey,ez
real(dp), intent(in), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
real(dp), intent(out), dimension(ilo:ihi,jlo:jhi,klo:khi) :: phi,hx,hy,hz,ex,ey,ez
integer :: icomp,ierr
real(dp), parameter :: clight=299792458.d0
real(dp), parameter :: mu0=8.d0*asin(1.d0)*1.d-7
integer :: i,j,k
real(dp) :: gb0
integer :: mprocs,myrank,mpierr
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)


gb0=sqrt((gam0+1.d0)*(gam0-1.d0))

if(idirectfieldcalc.eq.0)then
  if(myrank.eq.0)write(6,*)'computing phi with OpenBC'
  icomp=0
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
         ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
!the following needs work due to incorrect location of the 1-sided diff and due to missing boundary value; fix later
  do k=klo,khi-1
    do j=jlo,jhi-1
      do i=ilo,ihi-1
        if(icomp.eq.0)then
          ex(i,j,k)=-(phi(i+1,j,k)-phi(i,j,k))/dx              !if desired, edit for centered diff and incr starting do index
          ey(i,j,k)=-(phi(i,j+1,k)-phi(i,j,k))/dy
          ez(i,j,k)=-(phi(i,j,k+1)-phi(i,j,k))/dz/gamma_ave**2 !Ez=-d/dz(phi(beta-ct) - beta/c d/dt(phi(beta-ct)) 
        endif
    enddo
  enddo
  enddo
endif
  
if(idirectfieldcalc.eq.1)then
  if(myrank.eq.0)write(6,*)'computing Ex,Ey,Ez with OpenBC'

  icomp=1 
  !if(myrank.eq.0)print *,'old'
  !call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
  !ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  !if(myrank.eq.0)print *,'new'
  call Zopenbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
  ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  do k=klo,khi
    do j=jlo,jhi
      do i=ilo,ihi
        ex(i,j,k)=phi(i,j,k)
      enddo
     enddo
   enddo
   
  icomp=2
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  do k=klo,khi
    do j=jlo,jhi
      do i=ilo,ihi
        ey(i,j,k)=phi(i,j,k)
      enddo
    enddo
  enddo
        
  icomp=3
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  do k=klo,khi
    do j=jlo,jhi
      do i=ilo,ihi
        ez(i,j,k)=phi(i,j,k)
      enddo
    enddo
  enddo
  
endif
      
! set the magnetic field:
do k=klo,khi
  do j=jlo,jhi
    do i=ilo,ihi
      hx(i,j,k)=-ey(i,j,k)/clight/mu0*gb0/gam0
      hy(i,j,k)= ex(i,j,k)/clight/mu0*gb0/gam0
      hz(i,j,k)=0.d0
    enddo
  enddo
enddo

end subroutine getfields







! Test only, uses same interface as openbcpotential
subroutine Zopenbcpotential(rho,phi,gam, &
    dx,dy,dz, &
     ilo,ihi,jlo,jhi,klo,khi, &
       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl, &
       idecomp,npx,npy,npz,icomp,igfflag,ierr)
 use mpi      
implicit none       
real(dp) :: gam,dx,dy,dz
integer :: ilo,ihi,jlo,jhi,klo,khi
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho,phi

type(mesh3d_struct) :: mesh3d
type (domain_decomposition_struct) :: domain

! Init domain
call MPI_COMM_SIZE(MPI_COMM_WORLD,domain%n_process,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,domain%rank,ierr)
domain%idecomp = idecomp
domain%global%lo = [ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl]
domain%global%hi = [ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl]
call init_domain_decomposition(domain, MPI_COMM_WORLD)

! Init mesh

mesh3d%delta = [dx, dy, dz]
mesh3d%gamma = gam
call allocate_real_3d(mesh3d%rho, domain%local)
call allocate_real_3d(mesh3d%phi, domain%local)
mesh3d%rho = rho

if (domain%local%lo(1) /= ilo) print *, 'error'

call Xopenbcpotential(mesh3d, domain, icomp, igfflag, ierr)

phi = mesh3d%phi

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine  Xopenbcpotential(mesh3d, domain, icomp, igfflag, ierr)

implicit none

type(mesh3d_struct) :: mesh3d
type (domain_decomposition_struct) :: domain, domain2, domainGF
real(dp) :: gam,dx,dy,dz
integer :: icomp,igfflag,ierr
complex(dp), allocatable, dimension(:,:,:) :: crho2,cphi2,cgrn1 !what is the "1" for in cgrn1?
real(dp) :: u,v,w,gval, time
integer :: idirection
integer :: i,j,k,ip,jp,kp, myrank


real(dp), parameter :: econst=299792458.d0**2*1.d-7
!
dx = mesh3d%delta(1)
dy = mesh3d%delta(2)
dz = mesh3d%delta(3)
gam =  mesh3d%gamma

!determine the lower and upper global indices of the single- and double-size arrays
!double-size charge density:
domain2 = domain
domain2%global%hi = domain2%global%lo + 2*domain2%global%n -1
call init_domain_decomposition(domain2)
call allocate_complex_3d(crho2, domain2%local)
call allocate_complex_3d(cphi2, domain2%local)


!
!green function (note: this has negative and positive indices, is in the convolution formula)
domainGF = domain
domainGF%global%lo = domain%global%lo - domain%global%hi
domainGF%global%hi = domain%global%hi - domain%global%lo +1 !+1 is padding
call init_domain_decomposition(domainGF)
call allocate_complex_3d(cgrn1, domainGF%local)



!store rho in a double-size complex array:
call movetodoublesizer2c(mesh3d%rho,crho2, &
domain%local%lo(1), domain%local%hi(1), &
domain%local%lo(2), domain%local%hi(2), &
domain%local%lo(3), domain%local%hi(3), &
domain2%local%lo(1), domain2%local%hi(1), &
domain2%local%lo(2), domain2%local%hi(2), &
domain2%local%lo(3), domain2%local%hi(3), &
domain2%global%lo(1), domain2%global%hi(1), &
domain2%global%lo(2), domain2%global%hi(2), &
domain2%global%lo(3), domain2%global%hi(3), &
domain2%idecomp, &
domain2%process%n(1), domain2%process%n(2), domain2%process%n(3))

    
! fft the charge density in place
idirection=1 
call Xfft_perform(crho2, domain2, idirection)



!     write(6,*)'igfflag,icomp=',igfflag,icomp
do k=domainGF%local%lo(3), domainGF%local%hi(3)
  kp=domainGF%global%lo(3)+mod(k-domainGF%global%lo(3)+domain2%global%n(3)/2-1,domain2%global%n(3) )
  w=kp*dz
  do j=domainGF%local%lo(2), domainGF%local%hi(2)
    jp=domainGF%global%lo(2)+mod(j-domainGF%global%lo(2)+domain2%global%n(2) /2-1,domain2%global%n(2) )
    v=jp*dy
   do i=domainGF%local%lo(1), domainGF%local%hi(1)
     ip=domainGF%global%lo(1)+mod(i-domainGF%global%lo(1)+domain2%global%n(1)/2-1,domain2%global%n(1) )
     u=ip*dx
!igfflag=0 for point-charge green function, =1 to use IGF i.e. to integrate point charge green function over a cell
!icomp=0 means calculate the potential; icomp=1,2,3 mean calculate Ex,Ey,Ez
!note: should modify this routine to return all field components simultaneously, not separately
!note: should provide missing options, e.g. ifflag=0 and icomp= 1 or 2 or 3
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn1(i,j,k)=cmplx(gval,0.d0)
    enddo
  enddo
enddo

! fft the Green function in place:
idirection=1
call Xfft_perform(cgrn1, domain2, idirection)

! multiply the fft'd charge density and green function:
cphi2(:,:,:)=crho2(:,:,:)*cgrn1(:,:,:)

! now do the inverse fft in place
idirection=-1
call Xfft_perform(cphi2, domain2, idirection)

!store the physical portion of the double-size complex array in a non-double-size real array:
call movetosinglesizec2r(cphi2, mesh3d%phi, &
domain2%local%lo(1), domain2%local%hi(1), &
domain2%local%lo(2), domain2%local%hi(2), &
domain2%local%lo(3), domain2%local%hi(3), &
domain%local%lo(1), domain%local%hi(1), &
domain%local%lo(2), domain%local%hi(2), &
domain%local%lo(3), domain%local%hi(3), &
domain%global%lo(1), domain%global%hi(1), &
domain%global%lo(2), domain%global%hi(2), &
domain%global%lo(3), domain%global%hi(3), &
domain2%idecomp, &
domain2%process%n(1), domain2%process%n(2), domain2%process%n(3))

! Why divide by double-sized n? -Chris
mesh3d%phi=mesh3d%phi*econst / &
  ((1.d0*domain2%global%n(1))*(1.d0*domain2%global%n(2))*(1.d0*domain2%global%n(3)))

end subroutine Xopenbcpotential



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine openbcpotential(rho,phi,gam,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
           ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
use mpi
use decomposition_mod
implicit none
real(dp) :: gam,dx,dy,dz
integer :: ilo,ihi,jlo,jhi,klo,khi
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho,phi
!
complex(dp), allocatable, dimension(:,:,:) :: crho2,cphi2,cgrn1 !what is the "1" for in cgrn1?
integer :: ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl
integer :: ilo_phi_gbl,ihi_phi_gbl,jlo_phi_gbl,jhi_phi_gbl,klo_phi_gbl,khi_phi_gbl
integer :: ilo_grn_gbl,ihi_grn_gbl,jlo_grn_gbl,jhi_grn_gbl,klo_grn_gbl,khi_grn_gbl
integer :: ilo2,ihi2,jlo2,jhi2,klo2,khi2
integer :: iloo,ihii,jloo,jhii,kloo,khii
integer :: iperiod,jperiod,kperiod,n
integer :: idecomp_in,idecomp_out,idirection,ipermute,iscale
integer :: in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
integer :: out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
real(dp) :: time
integer :: i,j,k,ip,jp,kp,ndx
real(dp) :: u,v,w,gval

integer :: mprocs,myrank
real(dp), parameter :: econst=299792458.d0**2*1.d-7
!

call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!determine the lower and upper global indices of the single- and double-size arrays
!double-size charge density:
ilo_rho2_gbl=ilo_rho_gbl
jlo_rho2_gbl=jlo_rho_gbl
klo_rho2_gbl=jlo_rho_gbl
ihi_rho2_gbl=ilo_rho2_gbl+2*(ihi_rho_gbl-ilo_rho_gbl+1)-1
jhi_rho2_gbl=jlo_rho2_gbl+2*(jhi_rho_gbl-jlo_rho_gbl+1)-1
khi_rho2_gbl=klo_rho2_gbl+2*(khi_rho_gbl-klo_rho_gbl+1)-1
!
!single-size potential:
ilo_phi_gbl=ilo_rho_gbl
ihi_phi_gbl=ihi_rho_gbl
jlo_phi_gbl=jlo_rho_gbl
jhi_phi_gbl=jhi_rho_gbl
klo_phi_gbl=klo_rho_gbl
khi_phi_gbl=khi_rho_gbl
!
!green function (note: this has negative and positive indices, is in the convolution formula)
ilo_grn_gbl=ilo_phi_gbl-ihi_rho_gbl
ihi_grn_gbl=ihi_phi_gbl-ilo_rho_gbl+1 !+1 is padding
jlo_grn_gbl=jlo_phi_gbl-jhi_rho_gbl
jhi_grn_gbl=jhi_phi_gbl-jlo_rho_gbl+1 !+1 is padding
klo_grn_gbl=klo_phi_gbl-khi_rho_gbl
khi_grn_gbl=khi_phi_gbl-klo_rho_gbl+1 !+1 is padding
!
!the period is the size of the doubled array:
iperiod=ihi_rho2_gbl-ilo_rho2_gbl+1
jperiod=jhi_rho2_gbl-jlo_rho2_gbl+1
kperiod=khi_rho2_gbl-klo_rho2_gbl+1
!
 
!allocate the double-size complex array crho2:
call decompose(myrank,mprocs,ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,&
               idecomp,npx,npy,npz,ilo2,ihi2,jlo2,jhi2,klo2,khi2)

allocate(crho2(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
allocate(cphi2(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
!

!store rho in a double-size complex array:
call movetodoublesizer2c(rho,crho2,ilo,ihi,jlo,jhi,klo,khi,ilo2,ihi2,jlo2,jhi2,klo2,khi2, &
                         ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,idecomp,npx,npy,npz)

!
!!prepare for fft:
idecomp_in=idecomp
idecomp_out=idecomp
idirection=1
ipermute=0
iscale=0
call decompose(myrank,mprocs, &
               ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl, &
               idecomp_in,npx,npy,npz,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
out_ilo = in_ilo
out_ihi = in_ihi
out_jlo = in_jlo
out_jhi = in_jhi
out_klo = in_klo
out_khi = in_khi               

! fft the charge density in place
call fft_perform(crho2,crho2,idirection,iperiod,jperiod,kperiod,  &
     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     ipermute,iscale,time)

! compute the Green function (called cgrn1 below):
call decompose(myrank,mprocs,& !why does the following line contain rho2 indices?????
ilo_grn_gbl,ihi_grn_gbl,jlo_grn_gbl,jhi_grn_gbl,klo_grn_gbl,khi_grn_gbl,idecomp,npx,npy,npz,iloo,ihii,jloo,jhii,kloo,khii)

allocate(cgrn1(iloo:ihii,jloo:jhii,kloo:khii))
!     write(6,*)'igfflag,icomp=',igfflag,icomp
do k=kloo,khii
  kp=klo_grn_gbl+mod(k-klo_grn_gbl+kperiod/2-1,kperiod)
  w=kp*dz
  do j=jloo,jhii
    jp=jlo_grn_gbl+mod(j-jlo_grn_gbl+jperiod/2-1,jperiod)
    v=jp*dy
   do i=iloo,ihii
     ip=ilo_grn_gbl+mod(i-ilo_grn_gbl+iperiod/2-1,iperiod)
     u=ip*dx
!igfflag=0 for point-charge green function, =1 to use IGF i.e. to integrate point charge green function over a cell
!icomp=0 means calculate the potential; icomp=1,2,3 mean calculate Ex,Ey,Ez
!note: should modify this routine to return all field components simultaneously, not separately
!note: should provide missing options, e.g. ifflag=0 and icomp= 1 or 2 or 3
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn1(i,j,k)=cmplx(gval,0.d0)
    enddo
  enddo
enddo           
      
! fft the Green function in place:
call fft_perform(cgrn1,cgrn1,idirection,iperiod,jperiod,kperiod,  &
     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     ipermute,iscale,time)

! multiply the fft'd charge density and green function:
cphi2(:,:,:)=crho2(:,:,:)*cgrn1(:,:,:)


! now do the inverse fft in place
idirection=-1
call fft_perform(cphi2,cphi2,idirection,iperiod,jperiod,kperiod,    &
     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
     ipermute,iscale,time)


!store the physical portion of the double-size complex array in a non-double-size real array:
call movetosinglesizec2r(cphi2,phi,ilo2,ihi2,jlo2,jhi2,klo2,khi2,ilo,ihi,jlo,jhi,klo,khi, &
                         ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz)

phi(:,:,:)=phi(:,:,:)*econst /((1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod))

end subroutine openbcpotential
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
 function coulombfun(u,v,w,gam) result(res)
!computes 1/r where r=sqrt(u**2+v**2+(gam*w)**2)
 implicit none
 real(dp) :: res
 real(dp) :: u,v,w,gam
 if(u.eq.0.d0 .and. v.eq.0.d0 .and. w.eq.0.d0)then
   res=0.d0
   return
 endif
 res=gam/sqrt(u**2+v**2+(gam*w)**2)  !coulomb
 return
 end function coulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfcoulombfun(u,v,w,gam,dx,dy,dz) result(res)
!computes the definite integral over a cell of 1/r. Cell vertices are combinations of (x1,x2),(y1,y2),(z1,z2)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
!-!         real(dp), external :: lafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
!     res=1.d0/sqrt(u**2+v**2+w**2)  !coulomb
!     res=u/(u**2+v**2+w**2)**1.5d0  !x-electric field
res=lafun(x2,y2,z2)-lafun(x1,y2,z2)-lafun(x2,y1,z2)-lafun(x2,y2,z1)-lafun(x1,y1,z1)+ &
    lafun(x1,y1,z2)+lafun(x1,y2,z1)+lafun(x2,y1,z1)
res=res/(dx*dy*dz)

end function igfcoulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function lafun(x,y,z) result(res)
!computes the indefinite integral of 1/r
! lafun is the function involving log and atan in the PRSTAB paper (I should find a better name for this function)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
res=-0.5d0*z**2*atan(x*y/(z*r))-0.5d0*y**2*atan(x*z/(y*r))-0.5d0*x**2*atan(y*z/(x*r)) &
     +y*z*log(x+r)+x*z*log(y+r)+x*y*log(z+r)

end function lafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!
!In the following, igfexfun and and xlafun are analogous to igfcoulombfun and lafun,
!but instead of computing the scalar potential, the quantity computed is Ex.
!Similary, there are functions further below for the computation of Ey, Ez.
!
function igfexfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !-13 causes NaN's
real(dp), parameter :: em=-1.d-10
!-!         real(dp), external :: xlafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=xlafun(em,em,em) &
     -xlafun(x1,em,em)-xlafun(em,y1,em)-xlafun(em,em,z1)-xlafun(x1,y1,z1) &
     +xlafun(x1,y1,em)+xlafun(x1,em,z1)+xlafun(em,y1,z1)+xlafun(x2,em,em) &
     -xlafun(ep,em,em)-xlafun(x2,y1,em)-xlafun(x2,em,z1)-xlafun(ep,y1,z1) &
     +xlafun(ep,y1,em)+xlafun(ep,em,z1)+xlafun(x2,y1,z1)+xlafun(ep,y2,ep) &
     -xlafun(x1,y2,ep)-xlafun(ep,ep,ep)-xlafun(ep,y2,z1)-xlafun(x1,ep,z1) &
     +xlafun(x1,ep,ep)+xlafun(x1,y2,z1)+xlafun(ep,ep,z1)+xlafun(ep,ep,z2) &
     -xlafun(x1,ep,z2)-xlafun(ep,y1,z2)-xlafun(ep,ep,ep)-xlafun(x1,y1,ep) &
     +xlafun(x1,y1,z2)+xlafun(x1,ep,ep)+xlafun(ep,y1,ep)+xlafun(x2,y2,ep) &
     -xlafun(ep,y2,ep)-xlafun(x2,ep,ep)-xlafun(x2,y2,z1)-xlafun(ep,ep,z1) &
     +xlafun(ep,ep,ep)+xlafun(ep,y2,z1)+xlafun(x2,ep,z1)+xlafun(x2,ep,z2) &
     -xlafun(ep,ep,z2)-xlafun(x2,y1,z2)-xlafun(x2,ep,ep)-xlafun(ep,y1,ep) &
     +xlafun(ep,y1,z2)+xlafun(ep,ep,ep)+xlafun(x2,y1,ep)+xlafun(ep,y2,z2) &
     -xlafun(x1,y2,z2)-xlafun(ep,ep,z2)-xlafun(ep,y2,ep)-xlafun(x1,ep,ep) &
     +xlafun(x1,ep,z2)+xlafun(x1,y2,ep)+xlafun(ep,ep,ep)+xlafun(x2,y2,z2) &
     -xlafun(ep,y2,z2)-xlafun(x2,ep,z2)-xlafun(x2,y2,ep)-xlafun(ep,ep,ep) &
     +xlafun(ep,ep,z2)+xlafun(ep,y2,ep)+xlafun(x2,ep,ep)
else
  res=xlafun(x2,y2,z2)-xlafun(x1,y2,z2)-xlafun(x2,y1,z2)-xlafun(x2,y2,z1) &
      -xlafun(x1,y1,z1)+xlafun(x1,y1,z2)+xlafun(x1,y2,z1)+xlafun(x2,y1,z1)
endif
res=res/(dx*dy*dz)
return
end function igfexfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfeyfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !should include em also, but probably doesn't matter
!-!         real(dp), external :: ylafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=ylafun(ep,ep,ep) &
     -ylafun(x1,ep,ep)-ylafun(ep,y1,ep)-ylafun(ep,ep,z1)-ylafun(x1,y1,z1) &
     +ylafun(x1,y1,ep)+ylafun(x1,ep,z1)+ylafun(ep,y1,z1)+ylafun(x2,ep,ep) &
     -ylafun(ep,ep,ep)-ylafun(x2,y1,ep)-ylafun(x2,ep,z1)-ylafun(ep,y1,z1) &
     +ylafun(ep,y1,ep)+ylafun(ep,ep,z1)+ylafun(x2,y1,z1)+ylafun(ep,y2,ep) &
     -ylafun(x1,y2,ep)-ylafun(ep,ep,ep)-ylafun(ep,y2,z1)-ylafun(x1,ep,z1) &
     +ylafun(x1,ep,ep)+ylafun(x1,y2,z1)+ylafun(ep,ep,z1)+ylafun(ep,ep,z2) &
     -ylafun(x1,ep,z2)-ylafun(ep,y1,z2)-ylafun(ep,ep,ep)-ylafun(x1,y1,ep) &
     +ylafun(x1,y1,z2)+ylafun(x1,ep,ep)+ylafun(ep,y1,ep)+ylafun(x2,y2,ep) &
     -ylafun(ep,y2,ep)-ylafun(x2,ep,ep)-ylafun(x2,y2,z1)-ylafun(ep,ep,z1) &
     +ylafun(ep,ep,ep)+ylafun(ep,y2,z1)+ylafun(x2,ep,z1)+ylafun(x2,ep,z2) &
     -ylafun(ep,ep,z2)-ylafun(x2,y1,z2)-ylafun(x2,ep,ep)-ylafun(ep,y1,ep) &
     +ylafun(ep,y1,z2)+ylafun(ep,ep,ep)+ylafun(x2,y1,ep)+ylafun(ep,y2,z2) &
     -ylafun(x1,y2,z2)-ylafun(ep,ep,z2)-ylafun(ep,y2,ep)-ylafun(x1,ep,ep) &
     +ylafun(x1,ep,z2)+ylafun(x1,y2,ep)+ylafun(ep,ep,ep)+ylafun(x2,y2,z2) &
     -ylafun(ep,y2,z2)-ylafun(x2,ep,z2)-ylafun(x2,y2,ep)-ylafun(ep,ep,ep) &
     +ylafun(ep,ep,z2)+ylafun(ep,y2,ep)+ylafun(x2,ep,ep)
else
  res=ylafun(x2,y2,z2)-ylafun(x1,y2,z2)-ylafun(x2,y1,z2)-ylafun(x2,y2,z1) &
     -ylafun(x1,y1,z1)+ylafun(x1,y1,z2)+ylafun(x1,y2,z1)+ylafun(x2,y1,z1)
endif

res=res/(dx*dy*dz)

end function igfeyfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
function igfezfun(u,v,w,gam,dx,dy,dz) result(res)
implicit none
real(dp) :: res
real(dp) :: u,v,w,gam,dx,dy,dz
real(dp) :: x1,x2,y1,y2,z1,z2
real(dp) :: x,y,z
real(dp), parameter :: ep=1.d-10 !should include em also, but probably doesn't matter
!-!         real(dp), external :: zlafun
x1=u-0.5d0*dx
x2=u+0.5d0*dx
y1=v-0.5d0*dy
y2=v+0.5d0*dy
z1=(w-0.5d0*dz)*gam
z2=(w+0.5d0*dz)*gam
if(x1<0.d0.and.x2>0.d0.and.y1<0.d0.and.y2>0.d0.and.z1<0.d0.and.z2>0.d0)then
  res=zlafun(ep,ep,ep) &
     -zlafun(x1,ep,ep)-zlafun(ep,y1,ep)-zlafun(ep,ep,z1)-zlafun(x1,y1,z1) &
     +zlafun(x1,y1,ep)+zlafun(x1,ep,z1)+zlafun(ep,y1,z1)+zlafun(x2,ep,ep) &
     -zlafun(ep,ep,ep)-zlafun(x2,y1,ep)-zlafun(x2,ep,z1)-zlafun(ep,y1,z1) &
     +zlafun(ep,y1,ep)+zlafun(ep,ep,z1)+zlafun(x2,y1,z1)+zlafun(ep,y2,ep) &
     -zlafun(x1,y2,ep)-zlafun(ep,ep,ep)-zlafun(ep,y2,z1)-zlafun(x1,ep,z1) &
     +zlafun(x1,ep,ep)+zlafun(x1,y2,z1)+zlafun(ep,ep,z1)+zlafun(ep,ep,z2) &
     -zlafun(x1,ep,z2)-zlafun(ep,y1,z2)-zlafun(ep,ep,ep)-zlafun(x1,y1,ep) &
     +zlafun(x1,y1,z2)+zlafun(x1,ep,ep)+zlafun(ep,y1,ep)+zlafun(x2,y2,ep) &
     -zlafun(ep,y2,ep)-zlafun(x2,ep,ep)-zlafun(x2,y2,z1)-zlafun(ep,ep,z1) &
     +zlafun(ep,ep,ep)+zlafun(ep,y2,z1)+zlafun(x2,ep,z1)+zlafun(x2,ep,z2) &
     -zlafun(ep,ep,z2)-zlafun(x2,y1,z2)-zlafun(x2,ep,ep)-zlafun(ep,y1,ep) &
     +zlafun(ep,y1,z2)+zlafun(ep,ep,ep)+zlafun(x2,y1,ep)+zlafun(ep,y2,z2) &
     -zlafun(x1,y2,z2)-zlafun(ep,ep,z2)-zlafun(ep,y2,ep)-zlafun(x1,ep,ep) &
     +zlafun(x1,ep,z2)+zlafun(x1,y2,ep)+zlafun(ep,ep,ep)+zlafun(x2,y2,z2) &
     -zlafun(ep,y2,z2)-zlafun(x2,ep,z2)-zlafun(x2,y2,ep)-zlafun(ep,ep,ep) &
     +zlafun(ep,ep,z2)+zlafun(ep,y2,ep)+zlafun(x2,ep,ep)
else
  res=zlafun(x2,y2,z2)-zlafun(x1,y2,z2)-zlafun(x2,y1,z2)-zlafun(x2,y2,z1) &
     -zlafun(x1,y1,z1)+zlafun(x1,y1,z2)+zlafun(x1,y2,z1)+zlafun(x2,y1,z1)
endif
res=res/(dx*dy*dz*gam) !note the factor of gam in the denominator here, as needed for Ez

end function igfezfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function xlafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
    !  if (z+r == 0 ) print *, 'ERROR:x, r, x+r', x, y, z,r, x+r
    !  if (y+r == 0 ) print *, 'ERROR:y, r, y+r', x, y, z,r, y+r
    !  if (z+r == 0 ) print *, 'ERROR:z, r, z+r', x, y, z,r, z+r
!     res=x*atan(y*z/(r*x))-z*atanh(r/y)-y*atanh(r/z)
!res=z-x*atan(z/x)+x*atan(y*z/(x*r))-z*log(y+r)-y*log(z+r)
res=z-x*atan(z/x)+x*atan(y*z/(x*r))
if (y+r /= 0) res = res -z*log(y+r)
if (z+r /= 0) res = res -y*log(z+r)
!     write(2,'(5(1pe12.5,1x))')x,y,z,r,res

end function xlafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function ylafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
!res=x-y*atan(x/y)+y*atan(z*x/(y*r))-x*log(z+r)-z*log(x+r)
res=x-y*atan(x/y)+y*atan(z*x/(y*r))
if (z+r /= 0) res = res -x*log(z+r)
if (x+r /= 0) res = res -z*log(x+r)   

end function ylafun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function zlafun(x,y,z) result(res)
implicit none
real(dp) :: res
real(dp) :: x,y,z,r
r=sqrt(x**2+y**2+z**2)
!res=y-z*atan(y/z)+z*atan(x*y/(z*r))-y*log(x+r)-x*log(y+r)
res=y-z*atan(y/z)+z*atan(x*y/(z*r))
if (x+r /= 0) res = res -y*log(x+r)
if (y+r /= 0) res = res -x*log(y+r)

end function zlafun


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine movetodoublesizer2c(rho,crho2,ilo,ihi,jlo,jhi,klo,khi,ilo2,ihi2,jlo2,jhi2,klo2,khi2, &
           ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,idecomp,npx,npy,npz)
!move the single-size real array rho to a double-size complex array crho2 for convolution over a double-size domain
!rho is stored in the lower left octant. the rest is zero.
use mpi
use data_movement_mod, only : lowner,dmrgrnk
implicit none
integer :: ilo,ihi,jlo,jhi,klo,khi,ilo2,ihi2,jlo2,jhi2,klo2,khi2
integer :: ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,idecomp,npx,npy,npz
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
complex(dp), dimension(ilo2:ihi2,jlo2:jhi2,klo2:khi2) :: crho2
!%%   integer, allocatable, dimension(:) :: iloa,ihia,jloa,jhia,kloa,khia
!#    real(dp), dimension(1:2,1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: gstuff,gtmp
integer, dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: iownermesh1d,itmp,irnk
real(dp), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: rownermesh1d
real(dp), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: rhotmp,rtmp
integer(INT64), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: indxtmp,jtmp
integer :: numvals,ifail
integer :: n,i,j,k,itot,jtot,iold,jold,kold
real(dp) :: aitot,ajtot
!     integer, external :: iowner    !I prefer not to have an external directive, alternative is to use a module
!     integer, external :: lowner    !I prefer not to have an external directive, alternative is to use a module
integer :: nxtot,nytot,nztot
real(dp) :: recipnx,recipny,recipnz,recipnxnpx,recipnynpy,recipnznpz
integer :: mprocs,myrank,ierr
integer :: iticks1,iticks2
integer, allocatable :: sdisp(:),scounts(:),rdisp(:),rcounts(:)
integer :: ssize,rsize
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
crho2(:,:,:)=(0.d0,0.d0) !this is important, since only data from the physical region will be moved
!
allocate(sdisp(0:mprocs-1),scounts(0:mprocs-1),rdisp(0:mprocs-1),rcounts(0:mprocs-1))
!%%   allocate(iloa(0:mprocs-1),ihia(0:mprocs-1),jloa(0:mprocs-1),jhia(0:mprocs-1),kloa(0:mprocs-1),khia(0:mprocs-1))
!%%   do n=0,mprocs-1
!%%    call decompose(n,mprocs,ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,&
!%%  &                 idecomp,npx,npy,npz,iloa(n),ihia(n),jloa(n),jhia(n),kloa(n),khia(n))
!%%   enddo
itot=ihi_rho2_gbl-ilo_rho2_gbl+1
jtot=jhi_rho2_gbl-jlo_rho2_gbl+1
aitot=1.d0*itot
ajtot=1.d0*jtot
!     if(myrank.eq.0)write(6,*)'(movetodouble) itot,jtot=',itot,jtot
n=0
nxtot=ihi_rho2_gbl-ilo_rho2_gbl+1
nytot=jhi_rho2_gbl-jlo_rho2_gbl+1
nztot=khi_rho2_gbl-klo_rho2_gbl+1
recipnx=1.d0/nxtot
recipny=1.d0/nytot
recipnz=1.d0/nztot
recipnxnpx=recipnx*npx
recipnynpy=recipny*npy
recipnznpz=recipnz*npz
!     call system_clock(count=iticks1)
do k=klo,khi
  do j=jlo,jhi
    do i=ilo,ihi
      n=n+1
!#    gstuff(1,n)=rho(i,j,k)
!#    gstuff(2,n)=(k-1)*aitot*ajtot + (j-1)*aitot + i !FIXED fixed
      rhotmp(n)=rho(i,j,k)
      indxtmp(n)=(k-1)*aitot*ajtot + (j-1)*aitot + i !FIXED fixed
!     iownermesh1d(n)=iowner(iloa,ihia,jloa,jhia,kloa,khia,mprocs,&
!    &       ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,i,j,k)
      iownermesh1d(n)=lowner(recipnxnpx,recipnynpy,recipnznpz,npx,npy,i,j,k)
    enddo
  enddo
enddo
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoubletripledo',iticks2-iticks1)
!!!!!now that I know who owns it, sort the array:
!     call system_clock(count=iticks1)
itmp(1:n)=iownermesh1d(1:n)
rownermesh1d(1:n)=iownermesh1d(1:n)
!#    gtmp(1:2,1:n)=gstuff(1:2,1:n)
jtmp(1:n)=indxtmp(1:n)
rtmp(1:n)=rhotmp(1:n)
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoublepredmrgrnk',iticks2-iticks1)
!     call system_clock(count=iticks1)
call dmrgrnk(rownermesh1d,irnk,n,n)
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoubledmrgrnk',iticks2-iticks1)
!     call system_clock(count=iticks1)
do k=1,n
  iownermesh1d(k)=itmp(irnk(k))
!#      gstuff(1:2,k)=gtmp(1:2,irnk(k))
  indxtmp(k)=jtmp(irnk(k))
 rhotmp(k)=rtmp(irnk(k))
enddo
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoubleafterdmrgrnk',iticks2-iticks1)
! create the send count array:
scounts(0:mprocs-1)=0
do k=1,n
  scounts(iownermesh1d(k))=scounts(iownermesh1d(k))+1
enddo
! tell the other procs how much data is coming:
call MPI_alltoall(scounts,1,MPI_INTEGER,rcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
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
!
!     call system_clock(count=iticks1)
call MPI_alltoallv(indxtmp,scounts,sdisp,MPI_INTEGER8,jtmp,rcounts,rdisp,MPI_INTEGER8,MPI_COMM_WORLD,ierr)
call MPI_alltoallv(rhotmp,scounts,sdisp,MPI_DOUBLE,rtmp,rcounts,rdisp,MPI_DOUBLE,MPI_COMM_WORLD,ierr)
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoublealltoallv',iticks2-iticks1)
!     numvals=size(rho,1)*size(rho,2)*size(rho,3) !the number of values to be moved
!     write(7000+myrank,*)ssize,rsize,numvals
!     flush(7000+myrank)
!     call system_clock(count=iticks1)
do n=1,rsize
  i=mod(jtmp(n)-1,int(itot, int64))+1
  j=mod((jtmp(n)-1)/itot,int(jtot, int64))+1
  k=(jtmp(n)-1)/(itot*jtot)+1
!       k=(jtmp(n)-1)/(itot*jtot)+1
!       j=(jtmp(n)-1-(k-1)*itot*jtot)/itot + 1
!       i=jtmp(n)-(k-1)*itot*jtot-(j-1)*itot
  crho2(i,j,k)=rtmp(n) !recall that crho2 is initialized to zero previously, as required.
enddo
!     call system_clock(count=iticks2)
!     call elapsedmma('movetodoublefinaldo',iticks2-iticks1)

end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine movetosinglesizec2r(crho2,rho,ilo2,ihi2,jlo2,jhi2,klo2,khi2,ilo,ihi,jlo,jhi,klo,khi, &
                         ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz)
!move the double-size complex array crho2 back to a single-size real array rho, dropping all but the lower left octant
use mpi
use data_movement_mod, only : lowner,dmrgrnk
implicit none
integer, parameter :: idebug=0
integer :: ilo,ihi,jlo,jhi,klo,khi,ilo2,ihi2,jlo2,jhi2,klo2,khi2
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
complex(dp), dimension(ilo2:ihi2,jlo2:jhi2,klo2:khi2) :: crho2
!%%   integer, allocatable, dimension(:) :: iloa,ihia,jloa,jhia,kloa,khia
!#    real(dp), dimension(1:2,1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: gstuff,gtmp
integer, dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: iownermesh1d,itmp,irnk
real(dp), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: rownermesh1d
real(dp), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: rhotmp,rtmp
integer(INT64), dimension(1:8*size(rho,1)*size(rho,2)*size(rho,3)) :: indxtmp,jtmp
integer :: numvals,ifail
integer :: n,i,j,k,ngbl,nsize,nsizegbl,itot,jtot,iold,jold,kold
real(dp) :: aitot,ajtot
!     integer, external :: iowner    !I prefer not to have an external directive, alternative is to use a module
!     integer, external :: lowner    !I prefer not to have an external directive, alternative is to use a module
integer :: nxtot,nytot,nztot
real(dp) :: recipnx,recipny,recipnz,recipnxnpx,recipnynpy,recipnznpz
integer :: mprocs,myrank,ierr
integer :: iticks1,iticks2
integer, allocatable :: sdisp(:),scounts(:),rdisp(:),rcounts(:)
integer :: ssize,rsize
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!
rho(:,:,:)=0.d0 !not necessary but it doesn't hurt, makes this routine look like the other
!
allocate(sdisp(0:mprocs-1),scounts(0:mprocs-1),rdisp(0:mprocs-1),rcounts(0:mprocs-1))
!%%   allocate(iloa(0:mprocs-1),ihia(0:mprocs-1),jloa(0:mprocs-1),jhia(0:mprocs-1),kloa(0:mprocs-1),khia(0:mprocs-1))
!%%   do n=0,mprocs-1
!%%    call decompose(n,mprocs,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,&
!%%  &                 idecomp,npx,npy,npz,iloa(n),ihia(n),jloa(n),jhia(n),kloa(n),khia(n))
!%%   enddo
itot=ihi_rho_gbl-ilo_rho_gbl+1
jtot=jhi_rho_gbl-jlo_rho_gbl+1
aitot=1.d0*itot
ajtot=1.d0*jtot
!     if(myrank.eq.0)write(6,*)'(movetosingle) itot,jtot=',itot,jtot
n=0
nxtot=ihi_rho_gbl-ilo_rho_gbl+1
nytot=jhi_rho_gbl-jlo_rho_gbl+1
nztot=khi_rho_gbl-klo_rho_gbl+1
recipnx=1.d0/nxtot
recipny=1.d0/nytot
recipnz=1.d0/nztot
recipnxnpx=recipnx*npx
recipnynpy=recipny*npy
recipnznpz=recipnz*npz
!     call system_clock(count=iticks1)
do k=klo2,khi2
  if(k.lt.klo_rho_gbl.or.k.gt.khi_rho_gbl)cycle !keep only the physical region
  do j=jlo2,jhi2
    if(j.lt.jlo_rho_gbl.or.j.gt.jhi_rho_gbl)cycle !keep only the physical region
    do i=ilo2,ihi2
      if(i.lt.ilo_rho_gbl.or.i.gt.ihi_rho_gbl)cycle !keep only the physical region
      n=n+1
!#      if(n.gt.size(gstuff,2))write(6,*)'trouble'
!#      gstuff(1,n)=real(crho2(i,j,k))
!#      gstuff(2,n)=(k-1)*aitot*ajtot + (j-1)*aitot + i !FIXED fixed
      rhotmp(n)=real(crho2(i,j,k))
      indxtmp(n)=(k-1)*aitot*ajtot + (j-1)*aitot + i !FIXED fixed
!true for rho but not for phi       if(gstuff(4,n).lt.0.d0)write(6,*)'(1) gstuff(4,n) error'
!       iownermesh1d(n)=iowner(iloa,ihia,jloa,jhia,kloa,khia,mprocs,&
!    &       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,i,j,k)
      iownermesh1d(n)=lowner(recipnxnpx,recipnynpy,recipnznpz,npx,npy,i,j,k)
    enddo
  enddo
enddo
!     call system_clock(count=iticks2)
!     call elapsedmma('movetosingletripledo',iticks2-iticks1)
!!!!!now that I know who owns it, sort the array:
itmp(1:n)=iownermesh1d(1:n)
rownermesh1d(1:n)=iownermesh1d(1:n)
!#    gtmp(1:2,1:n)=gstuff(1:2,1:n)
jtmp(1:n)=indxtmp(1:n)
rtmp(1:n)=rhotmp(1:n)
!     call system_clock(count=iticks1)
call dmrgrnk(rownermesh1d,irnk,n,n)
!     call system_clock(count=iticks2)
!     call elapsedmma('movetosingledmrgrnk',iticks2-iticks1)
!     call system_clock(count=iticks1)
do k=1,n
  iownermesh1d(k)=itmp(irnk(k))
!#      gstuff(1:2,k)=gtmp(1:2,irnk(k))
  indxtmp(k)=jtmp(irnk(k))
  rhotmp(k)=rtmp(irnk(k))
enddo
!     call system_clock(count=iticks2)
!     call elapsedmma('movetosingleafterdmrgrnk',iticks2-iticks1)
! create the send count array:
scounts(0:mprocs-1)=0
do k=1,n
  scounts(iownermesh1d(k))=scounts(iownermesh1d(k))+1
enddo
! tell the other procs how much data is coming:
call MPI_alltoall(scounts,1,MPI_INTEGER,rcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
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
!     call system_clock(count=iticks1)
!#    call MPI_alltoallv(indxtmp,scounts,sdisp,MPI_INTEGER8,jtmp,rcounts,rdisp,MPI_INTEGER8,MPI_COMM_WORLD,ierr)
!#    call MPI_alltoallv(rhotmp,scounts,sdisp,MPI_DOUBLE,rtmp,rcounts,rdisp,MPI_DOUBLE,MPI_COMM_WORLD,ierr)
call MPI_alltoallv(rhotmp,scounts,sdisp,MPI_DOUBLE,rho ,rcounts,rdisp,MPI_DOUBLE,MPI_COMM_WORLD,ierr)
!     call system_clock(count=iticks2)
!     call elapsedmma('movetosinglealltoallv',iticks2-iticks1)
!#    do n=1,rsize
!#      i=mod(jtmp(n)-1,itot)+1
!#      j=mod((jtmp(n)-1)/itot,jtot)+1
!#      k=(jtmp(n)-1)/(itot*jtot)+1
!#      rho(i,j,k)=rtmp(n)
!#    enddo

end subroutine

end module
