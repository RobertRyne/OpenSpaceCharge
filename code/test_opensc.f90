program openbc_test
! Code to test the OpenBC package

use open_spacecharge_mod
use test_mod

implicit none
integer, parameter :: dp = REAL64

type(mesh3d_struct) :: mesh3d
integer :: nx, ny, nz ! # of grid points
integer :: n_particle
integer, parameter :: n1=7 ! (x,gbx,y,gby,z,gbz)(t) ; 7 is lost particle flag (not needed)
integer:: idirectfieldcalc 
integer :: igfflag
integer, parameter :: iseed=1234567 ! seed for random # generator
integer :: disttype
real(dp) :: gaussiancutoff !cutoff if a Gaussian initial condition is used
real(dp) :: gamma, mc2, e_tot, bunch_charge
real(dp) :: chrgpermacro !charge per macroparticle
real(dp), parameter :: tval=0.d0 !the time, an argument passed to diagnostic routines
real(dp), allocatable, dimension(:,:) :: y                !particle array

integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,ilo,ihi,jlo,jhi,klo,khi !global array sizes
integer :: idecomp,npx,npy,npz ! irrelevant to this serial code
!
real(dp), dimension(6,6) :: sigmat !initial beam sigma (i.e. 2nd moment) matrix
real(dp), dimension(6) :: cent     !initial beam centroid
real(dp) :: gb0 !gamma*beta
real(dp) :: sigma_x, sigma_y, sigma_z

integer :: open_status, namelist_file
character(40) :: in_file


namelist / opensc_test_params / &
    nx, ny, nz, n_particle, e_tot, bunch_charge, &
    sigma_x, sigma_y, sigma_z, gaussiancutoff, &
    idirectfieldcalc, igfflag

! Namelist defaults
nx=64
ny=64
nz=64 
e_tot = 100e6
n_particle = 1000000
mc2 = 0.510998910d6
bunch_charge=0.25d-9
sigma_x = 0.001
sigma_y = 0.001
sigma_z = 0.0001
gaussiancutoff=5
idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
disttype=1 ! =0 for uniform, =1 for Gaussian
 
!Read namelist
in_file = 'test_opensc.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)
open(newunit=namelist_file, file = in_file, status = 'old', iostat=open_status)
if (open_status /= 0) then
  print *, 'Input file missing: ', in_file
  print *, 'Using defaults'
else 
  read (namelist_file, nml = opensc_test_params)
  close (namelist_file)
endif

print *, '------------------------'
WRITE(*, opensc_test_params)
print *, '------------------------'

gamma = e_tot/mc2
chrgpermacro=bunch_charge/n_particle

allocate(y(n_particle,n1))


! parameters needed to call openbc routines:
ilo_rho_gbl=1; ihi_rho_gbl=nx; jlo_rho_gbl=1; jhi_rho_gbl=ny; klo_rho_gbl=1; khi_rho_gbl=nz !assumes power-of-2
ilo=ilo_rho_gbl;ihi=ihi_rho_gbl; jlo=jlo_rho_gbl;jhi=jhi_rho_gbl; klo=klo_rho_gbl;khi=khi_rho_gbl !serial code
write(6,*)'grid sizes ihi,jhi,khi=',ihi,jhi,khi
idecomp=-1;npx=1;npy=1;npz=1 !none of these are relevant to this serial code
!
! initialize the particles:
write(6,*)'gamma=',gamma
sigmat(1:6,1:6)=0.d0
sigmat(1,1)=sigma_x**2; sigmat(3,3)=sigma_y**2; sigmat(5,5)=sigma_z**2 !2nd moment matrix of initial dist
cent(1:6)=0.d0
gb0=sqrt((gamma+1.d0)*(gamma-1.d0))
write(6,*)'beta0=',gb0/gamma
cent(6)=gb0  !if integrating in time, this is gamma*beta of the centroid (set to gamma if integrating in z)

print *, 'Generating distribution...'
call gendist(y,n1,n_particle,sigmat,gaussiancutoff,disttype,iseed)
y(:,7)=0.d0 !lost particle flag (or use for some other purpose)
y(1:n_particle,6)=y(1:n_particle,6)+cent(6)
if(disttype.eq.0)write(6,*)'...done computing initial 3D uniform spatial distribution w/ cold velocity distribution'
if(disttype.eq.1)write(6,*)'...done computing initial 3D Gaussian spatial distribution w/ cold velocity distribution'
!

print *, 'depositing bunch on mesh...'
mesh3d%n = [nx, ny, nz]
mesh3d%gamma = gamma
call deposit_bunch_on_mesh(y(:,1),y(:,3), y(:,5), y(:,7), bunch_charge, mesh3d)
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

