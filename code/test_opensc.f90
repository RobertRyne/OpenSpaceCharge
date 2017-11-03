program openbc_test
! Code to test the OpenBC package

use open_spacecharge_mod
use test_mod

implicit none

type(mesh3d_struct) :: mesh3d
! This code has no input file. Most of the problem parameters are set here:
 integer :: nx, ny, nz ! # of grid points
!integer, parameter :: nx=64, ny=64, nz=64 ! # of grid points
integer :: n_particle
integer, parameter :: n1=7 ! (x,gbx,y,gby,z,gbz)(t) ; 7 is lost particle flag (not needed)
integer, parameter :: idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
integer, parameter :: igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
integer, parameter :: iseed=1234567 ! seed for random # generator
integer, parameter :: disttype=1 ! =0 for uniform, =1 for Gaussian
real(8), parameter :: gaussiancutoff=5 !cutoff if a Gaussian initial condition is used
!real(8), parameter :: ekin=100.e6 !250 MeV
real(8), parameter :: gam0=100e6 /0.510998910d6 !gamma for electrons with kinetic energy ekin
!real(8), parameter :: gam0=490.23783418637822
real(8), parameter :: chrgperbunch=0.25d-9 !charge per bunch; in this case 1 nC
real(8) :: chrgpermacro !charge per macroparticle
real(8), parameter :: tval=0.d0 !the time, an argument passed to diagnostic routines
!
real(8), allocatable, dimension(:,:) :: y                !particle array
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,ilo,ihi,jlo,jhi,klo,khi !global array sizes
integer :: idecomp,npx,npy,npz ! irrelevant to this serial code
!
real(8) :: dx,dy,dz,t !spatial and temporal quantities
real(8), dimension(6,6) :: sigmat !initial beam sigma (i.e. 2nd moment) matrix
real(8), dimension(6) :: cent     !initial beam centroid
real(8) :: gb0 !gamma*beta


namelist / openbc_test_params / &
    nx, ny, nz, n_particle



nx=64
ny=64
nz=64
n_particle = 1000000
print *, 'Begin'
chrgpermacro=chrgperbunch/n_particle

allocate(y(n_particle,n1))


! parameters needed to call openbc routines:
ilo_rho_gbl=1; ihi_rho_gbl=nx; jlo_rho_gbl=1; jhi_rho_gbl=ny; klo_rho_gbl=1; khi_rho_gbl=nz !assumes power-of-2
ilo=ilo_rho_gbl;ihi=ihi_rho_gbl; jlo=jlo_rho_gbl;jhi=jhi_rho_gbl; klo=klo_rho_gbl;khi=khi_rho_gbl !serial code
write(6,*)'grid sizes ihi,jhi,khi=',ihi,jhi,khi
idecomp=-1;npx=1;npy=1;npz=1 !none of these are relevant to this serial code
!
! initialize the particles:
write(6,*)'gam0=',gam0
sigmat(1:6,1:6)=0.d0
sigmat(1,1)=1.d-6; sigmat(3,3)=1.d-6; sigmat(5,5)=(1.d-4)**2 !2nd moment matrix of initial dist
cent(1:6)=0.d0
gb0=sqrt((gam0+1.d0)*(gam0-1.d0))
write(6,*)'beta0=',gb0/gam0
cent(6)=gb0  !if integrating in time, this is gamma*beta of the centroid (set to gam0 if integrating in z)

print *, 'Generating distribution...'
call gendist(y,n1,n_particle,sigmat,gaussiancutoff,disttype,iseed)
y(:,7)=0.d0 !lost particle flag (or use for some other purpose)
y(1:n_particle,6)=y(1:n_particle,6)+cent(6)
if(disttype.eq.0)write(6,*)'...done computing initial 3D uniform spatial distribution w/ cold velocity distribution'
if(disttype.eq.1)write(6,*)'...done computing initial 3D Gaussian spatial distribution w/ cold velocity distribution'
!

print *, 'depositing bunch on mesh...'
mesh3d%n = [nx, ny, nz]
call deposit_bunch_on_mesh(y(:,1),y(:,3), y(:,5), y(:,7), chrgperbunch, mesh3d, gam0)
print *, 'space charge field calc...'
call space_charge_field_calc(mesh3d)
print *, '...done'



! diagnostics:
call prntall(0,n1,n_particle,nx,ny,nz,y, &
  mesh3d%field%B(1), mesh3d%field%B(2), mesh3d%field%B(3), &
  mesh3d%field%E(1), mesh3d%field%E(2), mesh3d%field%E(3), &
  tval, &
  mesh3d%delta(1),  mesh3d%delta(2),  mesh3d%delta(3), &
  mesh3d%min(1), mesh3d%min(2), mesh3d%min(3)) 
  
!=========================================
!

end program

