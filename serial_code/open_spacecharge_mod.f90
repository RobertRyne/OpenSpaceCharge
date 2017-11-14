module open_spacecharge_mod

use, intrinsic :: iso_fortran_env
use fft_interface_mod
!$ use omp_lib

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
  integer :: n(3) = [64, 64, 64]       ! Grid sizes in x, y, z (m)
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma
  real(dp), allocatable, dimension(:,:,:) :: charge            ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: phi               ! electric potential grid
  type(field_struct),  allocatable, dimension(:,:,:) :: field ! field grid
end type

contains


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine deposit_bunch_on_mesh(xa, ya, za, lostflag,  charge, mesh3d)
!use mpi
type(mesh3d_struct) :: mesh3d
real(dp) :: charge, chrgpermacro
real(dp), dimension(:) :: xa, ya, za, lostflag 
integer, parameter :: idecomp = 999, npx=999, npy=999, npz = 999, n1 = 6 
integer :: nraysp,maxrayp
integer :: error
integer :: mprocs,myrank,ierr

!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

if (.not. allocated(mesh3d%charge)) then
  allocate(mesh3d%charge(mesh3d%n(1),mesh3d%n(2),mesh3d%n(3)))
  allocate(mesh3d%phi(mesh3d%n(1),mesh3d%n(2),mesh3d%n(3)))
  allocate(mesh3d%field(mesh3d%n(1),mesh3d%n(2),mesh3d%n(3)))
endif

nraysp = size(xa)
maxrayp = nraysp
chrgpermacro =charge/ nraysp

call get_rho_and_mesh_spacing(xa, ya, za, lostflag, mesh3d%charge, mesh3d%delta(1),mesh3d%delta(2), mesh3d%delta(3), &
     mesh3d%min(1),mesh3d%min(2),mesh3d%min(3), chrgpermacro, &
     1, mesh3d%n(1), 1, mesh3d%n(2), 1, mesh3d%n(3), &
     1, mesh3d%n(1), 1, mesh3d%n(2), 1, mesh3d%n(3), &
     idecomp,npx,npy,npz,n1,nraysp,maxrayp)

mesh3d%max(1) = mesh3d%min(1)+(mesh3d%n(1)-1)*mesh3d%delta(1)
mesh3d%max(2) = mesh3d%min(2)+(mesh3d%n(2)-1)*mesh3d%delta(2)
mesh3d%max(3) = mesh3d%min(3)+(mesh3d%n(3)-1)*mesh3d%delta(3)

if(myrank.eq.0)call print_mesh3d(mesh3d)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine print_mesh3d(mesh3d)
type(mesh3d_struct) :: mesh3d
print *, '------------------------'
print *, 'Mesh: '
print *, 'n: ', mesh3d%n
print *, 'min: ', mesh3d%min
print *, 'max: ', mesh3d%max
print *, 'delta: ', mesh3d%delta
print *, 'gamma: ', mesh3d%gamma
print *, '------------------------'
end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine space_charge_field_calc(mesh3d, direct_field_calc, integrated_green_function)
type(mesh3d_struct) :: mesh3d
real(dp) :: gamma0
integer, parameter :: idecomp = 999, npx=999, npy=999, npz = 999
integer :: idirectfieldcalc=1 
integer :: igfflag=1 
logical, optional :: direct_field_calc, integrated_green_function

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0


call getfields(mesh3d%charge, mesh3d%gamma, &
  mesh3d%delta(1),mesh3d%delta(2), mesh3d%delta(3), &
  mesh3d%phi, &
  mesh3d%field%E(1), mesh3d%field%E(2), mesh3d%field%E(3), &
  mesh3d%field%B(1), mesh3d%field%B(2), mesh3d%field%B(3), &
  mesh3d%n(1), mesh3d%n(2), mesh3d%n(3), &
  idecomp,npx,npy,npz, &
  igfflag,idirectfieldcalc,  &
  1, mesh3d%n(1), 1, mesh3d%n(2), 1, mesh3d%n(3), &
  1, mesh3d%n(1), 1, mesh3d%n(2), 1, mesh3d%n(3) ) 

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine getfields(rho,gam0,dx,dy,dz,phi,ex,ey,ez,bx,by,bz,nx,ny,nz,idecomp,npx,npy,npz,igfflag,idirectfieldcalc, &
                     ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl)
!use mpi
implicit none
integer, intent(in) :: nx,ny,nz,idecomp,npx,npy,npz,igfflag,idirectfieldcalc
integer, intent(in) :: ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl
real(dp), intent(in) :: gam0,dx,dy,dz
real(dp), intent(in), dimension(:,:,:) :: rho
real(dp), intent(out), dimension(:,:,:) :: phi,bx,by,bz,ex,ey,ez
integer :: icomp
real(dp), parameter :: clight=299792458.d0
!real(dp), parameter :: mu0=8.d0*asin(1.d0)*1.d-7
integer :: i,j,k
real(dp) :: gb0
integer :: mprocs,myrank,ierr
!
!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)


gb0=sqrt((gam0+1.d0)*(gam0-1.d0))


if(idirectfieldcalc.eq.0)then
  if(myrank.eq.0)write(6,*)'computing phi with OpenBC'
  icomp=0
        call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
     &       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  if(myrank.eq.0)write(6,*)'...done'
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if(icomp.eq.0)then
          ex(i,j,k)=-(phi(i+1,j,k)-phi(i,j,k))/dx*gam0
          ey(i,j,k)=-(phi(i,j+1,k)-phi(i,j,k))/dy*gam0
          ez(i,j,k)=-(phi(i,j,k+1)-phi(i,j,k))/dz/gam0
        endif
      enddo
    enddo
  enddo
endif

if(idirectfieldcalc.eq.1)then
  if(myrank.eq.0)write(6,*)'computing Ex,Ey,Ez with OpenBC'
  icomp=1
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
     &       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  if(myrank.eq.0)print *, 'Ex test: ', phi(nx/2,ny/2,nz/2)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        ex(i,j,k)=phi(i,j,k)
      enddo
    enddo
  enddo
  icomp=2
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
     &       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  if(myrank.eq.0)print *, 'Ey test: ', phi(nx/2,ny/2,nz/2)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        ey(i,j,k)=phi(i,j,k)
      enddo
    enddo
  enddo
  icomp=3
  call openbcpotential(rho,phi,gam0,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
     &       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr)
  
  if(myrank.eq.0)print *, 'Ez test: ', phi(nx/2,ny/2,nz/2)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        ez(i,j,k)=phi(i,j,k)
      enddo
    enddo
  enddo
endif


! set the magnetic field:
do k=1,nz
  do j=1,ny
    do i=1,nx
      bx(i,j,k)=-ey(i,j,k)/clight/gb0/gam0
      by(i,j,k)= ex(i,j,k)/clight/gb0/gam0
      bz(i,j,k)=0.d0
    enddo
  enddo
enddo

end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

!routine for charge deposition
subroutine depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,      &
                             dx,dy,dz,xmin,ymin,zmin,nraysp,ilogbl,jlogbl,klogbl,ifail)
!use mpi
implicit none
integer, intent(in) :: nraysp !!-!# of particles per MPI process
integer, intent(in) :: ilogbl,jlogbl,klogbl
integer, intent(out) :: ifail
real(dp), intent(in) :: chrgpermacro
real(dp), intent(in), dimension(:) :: xa,ya,za,lostflag
integer :: ilo,jlo,klo,ihi,jhi,khi
real(dp), intent(out), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
real(dp) :: dx,dy,dz,xmin,ymin,zmin
real(dp) :: dxi,dyi,dzi,ab,de,gh,sumrho
integer :: n,ip,jp,kp
integer :: mprocs,myrank,ierr

!call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

ifail=0
dxi=1.d0/dx
dyi=1.d0/dy
dzi=1.d0/dz
rho(ilo:ihi,jlo:jhi,klo:khi)=0.d0
do n=1,nraysp
  if(lostflag(n).ne.0.d0)cycle
  ip=floor((xa(n)-xmin)*dxi+1) !this is 1-based; use ilogbl for the general case
  jp=floor((ya(n)-ymin)*dyi+1)
  kp=floor((za(n)-zmin)*dzi+1)
  ab=((xmin-xa(n))+ip*dx)*dxi
  de=((ymin-ya(n))+jp*dy)*dyi
  gh=((zmin-za(n))+kp*dz)*dzi
! this "if" statement slows things down, but I'm using it for debugging purposes
  if(ip<ilo .or. jp<jlo .or. kp<klo .or. ip>ihi .or. jp>jhi .or. kp>khi)then
    ifail=ifail+1
    cycle
  else
    rho(ip,jp,kp)=rho(ip,jp,kp)+ab*de*gh
    rho(ip,jp+1,kp)=rho(ip,jp+1,kp)+ab*(1.-de)*gh
    rho(ip,jp+1,kp+1)=rho(ip,jp+1,kp+1)+ab*(1.-de)*(1.-gh)
    rho(ip,jp,kp+1)=rho(ip,jp,kp+1)+ab*de*(1.-gh)
    rho(ip+1,jp,kp+1)=rho(ip+1,jp,kp+1)+(1.-ab)*de*(1.-gh)
    rho(ip+1,jp+1,kp+1)=rho(ip+1,jp+1,kp+1)+(1.-ab)*(1.-de)*(1.-gh)
    rho(ip+1,jp+1,kp)=rho(ip+1,jp+1,kp)+(1.-ab)*(1.-de)*gh
    rho(ip+1,jp,kp)=rho(ip+1,jp,kp)+(1.-ab)*de*gh
  endif
enddo
rho=rho*chrgpermacro
if(ifail.ne.0)write(6,*)'(depose_rho_scalar) ifail=',ifail,' on process ',myrank

end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Interpolation routine to calculate exyz for each particle

! Original:  subroutine ntrp3d(ptcls,msk,u,exyz,xmin,ymin,zmin,hx,hy,hz,n1,np,nx,ny,nz)
subroutine interpolate_mesh3d(x, y, z, E, mesh3d)
type(mesh3d_struct) ::  mesh3d
real(dp) :: x, y, z, E(3), H(3)
real(dp) :: hxi,hyi,hzi,ab,de,gh
integer :: ip,jp,kp,ip1,jp1,kp1
integer :: nflag
nflag=0
hxi=1.d0/mesh3d%delta(1); hyi=1.d0/mesh3d%delta(2); hzi=1.d0/mesh3d%delta(3)

ip=floor((x-mesh3d%min(1))*hxi+1)
jp=floor((y-mesh3d%min(2))*hyi+1)
kp=floor((z-mesh3d%min(3))*hzi+1)
if(ip<1 .or. ip>mesh3d%n(1)-1)then
  nflag=1
  write(6,*)'ierror: ip=', ip, (x-mesh3d%min(1)), (x-mesh3d%min(1))/mesh3d%delta(1)
  if(ip<1)then
    ip=1
  else
    ip=mesh3d%n(1)-1
  endif
endif

if(jp<1 .or. jp>mesh3d%n(2)-1)then
  nflag=1
  write(6,*)'jerror: jp=', jp
  write(6,*)ab,de,gh
  if(jp<1)then
    jp=1
  else
    jp=mesh3d%n(2)-1
  endif
endif

if(kp<1 .or. kp>mesh3d%n(3)-1)then
  nflag=1
  write(6,*)'kerror:  kp=',kp
!!!!!!!!!!write(6,*)ab,de,gh
  if(kp<1)then
    kp=1
  else
    kp=mesh3d%n(3)-1
  endif
endif
ab=((mesh3d%min(1)-x)+ip*mesh3d%delta(1))*hxi
de=((mesh3d%min(2)-y)+jp*mesh3d%delta(2))*hyi
gh=((mesh3d%min(3)-z)+kp*mesh3d%delta(3))*hzi
if(nflag.eq.1)then
  write(6,*)ab,de,gh
  nflag=0
endif

ip1=ip+1
jp1=jp+1
kp1=kp+1
E=mesh3d%field(ip,jp,kp)%E(:)*ab*de*gh+mesh3d%field(ip,jp1,kp)%E(:)*ab*(1.-de)*gh     &
  +mesh3d%field(ip,jp1,kp1)%E(:)*ab*(1.-de)*(1.-gh)+mesh3d%field(ip,jp,kp1)%E(:)*ab*de*(1.-gh)  &
  +mesh3d%field(ip1,jp,kp1)%E(:)*(1.-ab)*de*(1.-gh)                               &
  +mesh3d%field(ip1,jp1,kp1)%E(:)*(1.-ab)*(1.-de)*(1.-gh)                         &
  +mesh3d%field(ip1,jp1,kp)%E(:)*(1.-ab)*(1.-de)*gh+mesh3d%field(ip1,jp,kp)%E(:)*(1.-ab)*de*gh


end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine openbcpotential(rho,phi,gam,dx,dy,dz,ilo,ihi,jlo,jhi,klo,khi, &
                           ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl, &
                           idecomp,npx,npy,npz,icomp,igfflag,ierr)
!-!#ifdef MPIPARALLEL
!     USE mpi
!-!#endif
implicit none
real(dp) :: gam,dx,dy,dz
integer :: ilo,ihi,jlo,jhi,klo,khi
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,icomp,igfflag,ierr
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho,phi
!
complex(dp), allocatable, dimension(:,:,:) :: crho2,ctmp2,cphi2,cgrn1 !what is the "1" for in cgrn1?
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
!-!   real(dp), external :: coulombfun,igfcoulombfun !I prefer not to have an external directive, alternative is to use a module
!-!   real(dp), external :: igfexfun,igfeyfun,igfezfun
integer :: mprocs,myrank
real(dp), parameter :: econst=299792458.d0**2*1.d-7
!-!#ifdef MPIPARALLEL
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!-!#else
mprocs=1
myrank=0
!-!#endif
!
ilo_rho2_gbl=ilo_rho_gbl
jlo_rho2_gbl=jlo_rho_gbl
klo_rho2_gbl=jlo_rho_gbl
ihi_rho2_gbl=ilo_rho2_gbl+2*(ihi_rho_gbl-ilo_rho_gbl+1)-1
jhi_rho2_gbl=jlo_rho2_gbl+2*(jhi_rho_gbl-jlo_rho_gbl+1)-1
khi_rho2_gbl=klo_rho2_gbl+2*(khi_rho_gbl-klo_rho_gbl+1)-1
!
ilo_phi_gbl=ilo_rho_gbl
ihi_phi_gbl=ihi_rho_gbl
jlo_phi_gbl=jlo_rho_gbl
jhi_phi_gbl=jhi_rho_gbl
klo_phi_gbl=klo_rho_gbl
khi_phi_gbl=khi_rho_gbl
!
ilo_grn_gbl=ilo_phi_gbl-ihi_rho_gbl
ihi_grn_gbl=ihi_phi_gbl-ilo_rho_gbl+1 !+1 is padding
jlo_grn_gbl=jlo_phi_gbl-jhi_rho_gbl
jhi_grn_gbl=jhi_phi_gbl-jlo_rho_gbl+1 !+1 is padding
klo_grn_gbl=klo_phi_gbl-khi_rho_gbl
khi_grn_gbl=khi_phi_gbl-klo_rho_gbl+1 !+1 is padding
!
iperiod=ihi_rho2_gbl-ilo_rho2_gbl+1
jperiod=jhi_rho2_gbl-jlo_rho2_gbl+1
kperiod=khi_rho2_gbl-klo_rho2_gbl+1
!
!allocate the double-size complex array crho2:
!-!#ifdef MPIPARALLEL
!     call decompose(myrank,mprocs,ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,&
!    &               idecomp,npx,npy,npz,ilo2,ihi2,jlo2,jhi2,klo2,khi2)
!-!#else
      ilo2=ilo_rho2_gbl;ihi2=ihi_rho2_gbl; jlo2=jlo_rho2_gbl;jhi2=jhi_rho2_gbl; klo2=klo_rho2_gbl;khi2=khi_rho2_gbl
!-!#endif
allocate(crho2(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
allocate(ctmp2(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
allocate(cphi2(ilo2:ihi2,jlo2:jhi2,klo2:khi2)) !double-size array
!
!store rho in a double-size complex array:
!-!#ifdef MPIPARALLEL
!     call movetodoublesizer2c(rho,crho2,ilo,ihi,jlo,jhi,klo,khi,ilo2,ihi2,jlo2,jhi2,klo2,khi2, &
!    &                         ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl,idecomp,npx,npy,npz)
!-!#else
crho2(:,:,:)=0.d0
crho2(ilo:ihi,jlo:jhi,klo:khi)=rho(ilo:ihi,jlo:jhi,klo:khi)
!-!#endif
!
!-!#ifdef MPIPARALLEL
!!prepare for fft:
!     idecomp_in=idecomp
!     idecomp_out=idecomp
!     idirection=1
!     ipermute=0
!     iscale=0
!     call decompose(myrank,mprocs, &
!    & ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl, &
!    & idecomp_in,npx,npy,npz,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
!     call decompose(myrank,mprocs, &
!    & ilo_rho2_gbl,ihi_rho2_gbl,jlo_rho2_gbl,jhi_rho2_gbl,klo_rho2_gbl,khi_rho2_gbl, &
!    & idecomp_out,npx,npy,npz,out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi)
!-!#endif
!
! fft the charge density:
!-!#ifdef MPIPARALLEL
!     call fft_perform(crho2,ctmp2,idirection,iperiod,jperiod,kperiod,  &
!    &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
!    &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
!    &     ipermute,iscale,time)
!-!#else
!call ccfft3d(crho2,ctmp2,[1,1,1],iperiod,jperiod,kperiod,0)
call fftw_ccfft3d(crho2,ctmp2,[1,1,1],iperiod,jperiod,kperiod,0)
!-!#endif
crho2(:,:,:)=ctmp2(:,:,:)  !now ctmp2 can be reused for next fft
!
! compute the Green function (called cgrn1 below):
!-!#ifdef MPIPARALLEL
!     call decompose(myrank,mprocs,& !why does the following line contain rho2 indices?????
!    &ilo_grn_gbl,ihi_grn_gbl,jlo_grn_gbl,jhi_grn_gbl,klo_grn_gbl,khi_grn_gbl,idecomp,npx,npy,npz,iloo,ihii,jloo,jhii,kloo,khii)
!-!#else
iloo=ilo_grn_gbl;ihii=ihi_grn_gbl; jloo=jlo_grn_gbl;jhii=jhi_grn_gbl; kloo=klo_grn_gbl;khii=khi_grn_gbl
!-!#endif
allocate(cgrn1(iloo:ihii,jloo:jhii,kloo:khii))
!     write(6,*)'igfflag,icomp=',igfflag,icomp



!$ print *, 'omp_get_max_threads(): ', omp_get_max_threads()
!$ print *, 'Parallel Do'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn1)
do k=kloo,khii
  kp=klo_grn_gbl+mod(k-klo_grn_gbl+kperiod/2-1,kperiod)
  w=kp*dz
  do j=jloo,jhii
    jp=jlo_grn_gbl+mod(j-jlo_grn_gbl+jperiod/2-1,jperiod)
    v=jp*dy
   do i=iloo,ihii
     ip=ilo_grn_gbl+mod(i-ilo_grn_gbl+iperiod/2-1,iperiod)
     u=ip*dx
     if(igfflag.eq.0.and.icomp.eq.0)gval=coulombfun(u,v,w,gam)
     if(igfflag.eq.1.and.icomp.eq.0)gval=igfcoulombfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.1)gval=igfexfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.2)gval=igfeyfun(u,v,w,gam,dx,dy,dz)
     if(igfflag.eq.1.and.icomp.eq.3)gval=igfezfun(u,v,w,gam,dx,dy,dz)
     cgrn1(i,j,k)= cmplx(gval,0.d0, dp)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$ print *, 'End Parallel Do'      
      

      
! fft the Green function:
!-!#ifdef MPIPARALLEL
!     call fft_perform(cgrn1,ctmp2,idirection,iperiod,jperiod,kperiod,  &
!    &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
!    &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
!    &     ipermute,iscale,time)
!-!# else
!call ccfft3d(cgrn1,ctmp2,[1,1,1],iperiod,jperiod,kperiod,0)
call fftw_ccfft3d(cgrn1,ctmp2,[1,1,1],iperiod,jperiod,kperiod,0)
!-!#endif
! multiply the fft'd charge density and green function:
cphi2(:,:,:)=crho2(:,:,:)*ctmp2(:,:,:)
! now do the inverse fft:
      idirection=-1
!-!#ifdef MPIPARALLEL
!     call fft_perform(cphi2,ctmp2,idirection,iperiod,jperiod,kperiod,    &
!    &     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,                   &
!    &     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,             &
!    &     ipermute,iscale,time)
!-!#else
!call ccfft3d(cphi2,ctmp2,[-1,-1,-1],iperiod,jperiod,kperiod,0)
call fftw_ccfft3d(cphi2,ctmp2,[-1,-1,-1],iperiod,jperiod,kperiod,0)
!-!#endif
cphi2(:,:,:)=ctmp2(:,:,:)/((1.d0*iperiod)*(1.d0*jperiod)*(1.d0*kperiod))
!
!store the physical portion of the double-size complex array in a non-double-size real array:
!-!#ifdef MPIPARALLEL
!     call movetosinglesizec2r(cphi2,phi,ilo2,ihi2,jlo2,jhi2,klo2,khi2,ilo,ihi,jlo,jhi,klo,khi, &
!    &                         ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz)
!-!#else
phi(ilo:ihi,jlo:jhi,klo:khi)=real(cphi2(ilo:ihi,jlo:jhi,klo:khi), dp)
phi(ilo:ihi,jlo:jhi,klo:khi)=phi(ilo:ihi,jlo:jhi,klo:khi)*econst
!-!#endif

end subroutine
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
 function coulombfun(u,v,w,gam) result(res)
 implicit none
 real(dp) :: res
 real(dp) :: u,v,w,gam
 if(u.eq.0.d0 .and. v.eq.0.d0 .and. w.eq.0.d0)then
   res=0.d0
   return
 endif
 res=1.d0/sqrt(u**2+v**2+(gam*w)**2)  !coulomb
!     res=u/(u**2+v**2+(gam*w)**2)**1.5d0  !x-electric field
 return
 end function coulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function igfcoulombfun(u,v,w,gam,dx,dy,dz) result(res)
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
res=res/(dx*dy*dz*gam)

end function igfcoulombfun
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
function lafun(x,y,z) result(res)
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

subroutine get_rho_and_mesh_spacing(xa,ya,za,lostflag, rho,dx,dy,dz,xmin,ymin,zmin,chrgpermacro, &
                         ilo,ihi,jlo,jhi,klo,khi, &
                         ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl, &
                         idecomp,npx,npy,npz,n1,nraysp,maxrayp)
!-!#ifdef MPIPARALLEL
!     USE mpi
!-!#endif
implicit none

integer :: ilo,ihi,jlo,jhi,klo,khi,ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,idecomp,npx,npy,npz,n1,nraysp,maxrayp
real(dp) :: dx,dy,dz,xmin,ymin,zmin,chrgpermacro
!-! real(dp), dimension(maxrayp,n1) :: ptcls
!type (coord_struct) :: ptcls(maxrayp)
real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: rho
!
real(dp), dimension(maxrayp) :: xa,ya,za,lostflag !lostflag=1.0 if particle is "lost"
real(dp) :: xmax,ymax,zmax !not needed
integer :: ifail,n
integer :: mprocs,myrank,ierr
real(dp) :: xsml,xbig,ysml,ybig,zsml,zbig
!     real(dp), parameter :: eps=2.d-15   !3.34d-16 is OK on my Mac
real(dp), parameter :: eps=0.5d0
integer :: nx,ny,nz !temporaries
!-!#ifdef MPIPARALLEL
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!-!#else
mprocs=1
myrank=0
!-!#endif
!     if(myrank.eq.0)write(6,*)'hello from get_rho'
      !xa(1:maxrayp)=ptcls(:)%vec(1) !-! ptcls(1:maxrayp,1)
      !ya(1:maxrayp)=ptcls(:)%vec(3) !-! ptcls(1:maxrayp,3)
      !za(1:maxrayp)=ptcls(:)%vec(5) !-! ptcls(1:maxrayp,5)
      !lostflag(1:maxrayp)= ptcls(:)%state  !-! ptcls(1:maxrayp,7)
! compute the bounding box and dx,dy,dz so particles can be localized to the correct proc:
call getbeamboundingbox(xa,ya,za,lostflag,xsml,xbig,ysml,ybig,zsml,zbig,nraysp) !nraysp should be bigrayp?
xbig=xbig*(1.d0+sign(1.d0,xbig)*eps)
xsml=xsml*(1.d0-sign(1.d0,xsml)*eps)
ybig=ybig*(1.d0+sign(1.d0,ybig)*eps)
ysml=ysml*(1.d0-sign(1.d0,ysml)*eps)
zbig=zbig*(1.d0+sign(1.d0,zbig)*eps)
zsml=zsml*(1.d0-sign(1.d0,zsml)*eps)
nx=ihi-ilo+1
ny=jhi-jlo+1
nz=khi-klo+1
dx=(xbig-xsml)/(nx-3)
dy=(ybig-ysml)/(ny-3)
dz=(zbig-zsml)/(nz-3)
xmin=xsml-dx; xmax=xbig+dx
ymin=ysml-dy; ymax=ybig+dy
zmin=zsml-dz; zmax=zbig+dz
!
!
!-!#ifdef MPIPARALLEL
! move particles to the correct proc:
!     ptcls(1,1:nraysp)=xa(1:nraysp)
!     ptcls(2,1:nraysp)=ya(1:nraysp)
!     ptcls(3,1:nraysp)=za(1:nraysp)
!     if(myrank.eq.0)write(6,*)'calling localize'
!     call localize(ptcls,lostflag,xmin,ymin,zmin,dx,dy,dz,&
!    &              ilo_gbl,ihi_gbl,jlo_gbl,jhi_gbl,klo_gbl,khi_gbl,n1,nraysp,maxrayp,idecomp,npx,npy,npz,mprocs)
!     if(myrank.eq.0)write(6,*)'back from localize'
!!deposit charge on the grid:
!     xa(1:nraysp)=ptcls(1,1:nraysp)
!     ya(1:nraysp)=ptcls(2,1:nraysp)
!     za(1:nraysp)=ptcls(3,1:nraysp)
!-!#endif
lostflag(1:maxrayp)=1.d0
lostflag(1:nraysp)=0.d0    !all particles are included in this example run
call depose_rho_scalar(xa,ya,za,lostflag,rho,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi,dx,dy,dz,xmin,ymin,zmin,nraysp, &
                              ilo_gbl,jlo_gbl,klo_gbl,ifail)

end subroutine
!

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine getbeamboundingbox(x,y,z,lostflag,xmin,xmax,ymin,ymax,zmin,zmax,nraysp)
!-!#ifdef MPIPARALLEL
!     USE mpi
!-!#endif
implicit none
integer :: nraysp !!-!# of particles per MPI process
real(dp), dimension(*) :: x,y,z,lostflag !lostflag=1.0 if particle is "lost"
real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax
real(dp), dimension(6) :: veclcl,vecgbl
real(dp) :: xminlcl,xmaxlcl,yminlcl,ymaxlcl,zminlcl,zmaxlcl
real(dp) :: xwidthorig,ywidthorig,zwidthorig
integer :: ierror
!need this since, if nraysp=0, the next 6 statements are skipped
veclcl(1:3)=-9999999.
veclcl(4:6)=-9999999.
if(nraysp > 0) then
!     veclcl(1)=-minval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(2)=-minval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(3)=-minval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(4)=maxval(x(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(5)=maxval(y(1:nraysp),lostflag(1:nraysp).eq.0.d0)
!     veclcl(6)=maxval(z(1:nraysp),lostflag(1:nraysp).eq.0.d0)
  veclcl(1)=-minval(x(1:nraysp))
  veclcl(2)=-minval(y(1:nraysp))
  veclcl(3)=-minval(z(1:nraysp))
  veclcl(4)=maxval(x(1:nraysp))
  veclcl(5)=maxval(y(1:nraysp))
  veclcl(6)=maxval(z(1:nraysp))
endif
!-!#ifdef MPIPARALLEL
!     call MPI_ALLREDUCE(veclcl,vecgbl,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
!-!#else
vecgbl(1:6)=veclcl(1:6)
!-!#endif
xmin=-vecgbl(1)
ymin=-vecgbl(2)
zmin=-vecgbl(3)
xmax=vecgbl(4)
ymax=vecgbl(5)
zmax=vecgbl(6)

end subroutine

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
integer :: n
do n=1,ntot*lenfft,lenfft
 call ccfftnr(a(n),lenfft,idir)
 ! call gsl_fft(a(n),lenfft,idir)  
enddo
return
end

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine ccfftnr(cdata,nn,isign)
implicit none
integer nn,isign,i,j,n,m,istep,mmax
complex(dp) cdata
real(dp) data
dimension cdata(nn),data(2*nn)
real(dp) tempi,tempr,theta,wpr,wpi,wr,wi,wtemp,twopi
do i=1,nn
  data(2*i-1)=real(cdata(i))
  data(2*i) =aimag(cdata(i))
enddo
! bit reversal:
n=2*nn
j=1
do i=1,n,2
  if(j>i)then
    tempr=data(j)
    tempi=data(j+1)
    data(j)=data(i)
    data(j+1)=data(i+1)
    data(i)=tempr
    data(i+1)=tempi
   endif
  m=n/2
  do while ((m.ge.2).and.(j>m)) 
    j=j-m
    m=m/2
  enddo
  j=j+m
enddo
! Danielson-Lanczos:
twopi=4.0d0*asin(1.0d0)
mmax=2
do while (n>mmax)
  istep=2*mmax
  theta=twopi/(isign*mmax)
  wpr=-2.d0*sin(0.5d0*theta)**2
  wpi=sin(theta)
  wr=1.0
  wi=0.0
  do m=1,mmax,2
    do i=m,n,istep
      j=i+mmax
      tempr=wr*data(j)-wi*data(j+1)
      tempi=wr*data(j+1)+wi*data(j)
      data(j)=data(i)-tempr
      data(j+1)=data(i+1)-tempi
      data(i)=data(i)+tempr
      data(i+1)=data(i+1)+tempi
    enddo
    wtemp=wr
    wr=wr*wpr-wi*wpi+wr
    wi=wi*wpr+wtemp*wpi+wi
  enddo
  mmax=istep
enddo 
!     ezero=0.d0
!     eunit=1.d0
do i=1,nn
  cdata(i)=data(2*i-1)+(0.d0,1.d0)*data(2*i)
!       cdata(i)=data(2*i-1)+cmplx(ezero,eunit)*data(2*i)
enddo

end subroutine



end module
