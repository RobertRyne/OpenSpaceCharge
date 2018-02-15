! Test code for the OpenSC space-charge package
! R.D. Ryne and C. Mayes, February 2018
program opensc_test
use open_spacecharge_mod1
use open_spacecharge_mod, only : REAL64, osc_freespace_solver, osc_alloc_freespace_array, osc_rectpipe_solver, &
                                 osc_getgrnpipe, osc_alloc_rectpipe_arrays, osc_read_rectpipe_grn, osc_write_rectpipe_grn
use test_mod, only : gendist, get_mesh_quantities, depose_rho_scalar, prntall
implicit none
type(mesh3d_struct) :: mesh3d
integer, parameter :: dp = REAL64

real(dp), allocatable, dimension(:,:,:) :: rho,phi  !charge density, potential
real(dp), allocatable, dimension(:,:,:,:) :: efield,bfield  !fields
real(dp), allocatable, dimension(:,:) :: ptcl  !particle array
real(dp), allocatable, dimension(:) :: lostflag

integer, dimension(3) :: nlo,nhi !lower,upper indices of the charge density, potential, and field arrays
real(dp), dimension(3) :: umin,umax,delta !grid min, max, and spacing

!padding for FFT purposes. If the grid size is power-of-2, adding 1 to the domain of the free-space Green function
!leads to power-of-2 FFTs and a doubling of the charge density grid
integer, parameter :: ipad=1, jpad=1,kpad=1
integer, dimension(3) :: npad=(/ipad,jpad,kpad/)

integer :: n_particle ! # of particles used by the test program to produce the density array rho
integer, parameter :: n1=7 ! (x,gbx,y,gby,z,gbz)(t) ; 7 could be a lost particle flag (not used in this test code)

!particle distribution
integer :: disttype
real(dp) :: gaussiancutoff !cutoff if a Gaussian initial condition is used
real(dp) :: gamma, gambet, mc2, e_tot, bunch_charge
real(dp) :: sigma_x, sigma_y, sigma_z, chrgpermacro
real(dp), dimension(6,6) :: sigmat !initial beam sigma (i.e. 2nd moment) matrix
real(dp), dimension(6) :: cent     !initial beam centroid
integer :: iseed=1234567 ! seed for random # generator

!geometry
real(dp) :: apipe, bpipe !full width and height of rectangular pipe; ignored unless rectpipe is .true.

!Poisson solver parameters
logical :: direct_field_calc, integrated_green_function
logical :: rectpipe, read_rectpipe, write_rectpipe
integer :: igfflag,idirectfieldcalc

real(dp) :: tval=0.d0 !the time, an argument passed to diagnostic routines, irrelevant in this test code
integer :: ifail !return flag for charge deposition routine

integer :: open_status, namelist_file
character(40) :: in_file
integer :: nxlo,nxhi,nylo,nyhi,nzlo,nzhi
namelist / opensc_test_params / &
  nxlo,nxhi,nylo,nyhi,nzlo,nzhi, n_particle, e_tot, bunch_charge, sigma_x, sigma_y, sigma_z, gaussiancutoff, &
  direct_field_calc, integrated_green_function, &
  rectpipe, read_rectpipe, write_rectpipe, apipe, bpipe ! if rectangular pipe BC is being used, apipe=full width, bpipe=full height

!Namelist defaults
nxlo=1;  nylo=1;  nzlo=1
nxhi=64; nyhi=64; nzhi=64
e_tot = 5.d6    !5 MeV
n_particle = 30000000
mc2 = 0.510998910d6
bunch_charge=0.25d-9
sigma_x = 0.001d0; sigma_y = 0.001d0; sigma_z = 0.0001d0
gaussiancutoff=5
disttype=1 ! =0 for uniform, =1 for Gaussian
direct_field_calc = .false.
integrated_green_function = .true.
rectpipe=.false.
read_rectpipe=.false.
write_rectpipe=.false.
apipe=12.d-3 !the default beam has 5sigma=5mm, so make the rectangular pipe size +/- 6mm.
bpipe=12.d-3
 
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

!these control the size of the rho,phi,field grids and the length of the FFTs
nlo=(/nxlo,nylo,nzlo/)
nhi=(/nxhi,nyhi,nzhi/)
npad=(/ipad,jpad,kpad/)

gamma = e_tot/mc2
gambet=sqrt((gamma+1.d0)*(gamma-1.d0))
chrgpermacro=bunch_charge/n_particle
igfflag=1; if(.not.integrated_green_function)igfflag=0
idirectfieldcalc=0; if(direct_field_calc)idirectfieldcalc=1

allocate(ptcl(n_particle,n1))
allocate(lostflag(n_particle))

!initialize the particles:
write(6,*)'gamma=',gamma
sigmat(1:6,1:6)=0.d0
sigmat(1,1)=sigma_x**2; sigmat(3,3)=sigma_y**2; sigmat(5,5)=sigma_z**2 !2nd moment matrix of initial dist
cent(1:6)=0.d0
write(6,*)'beta0=',gambet/gamma
cent(6)=gambet  !if integrating in time, this is gamma*beta of the centroid (set to gamma if integrating in z)

call gendist(ptcl,n1,n_particle,sigmat,gaussiancutoff,disttype,iseed,apipe,bpipe,rectpipe)
ptcl(:,7)=0.d0 !lost particle flag (or use for some other purpose)
ptcl(1:n_particle,6)=ptcl(1:n_particle,6)+cent(6)
lostflag(1:n_particle)=0.d0
if(disttype.eq.0)write(6,*)'Done computing initial 3D uniform spatial distribution w/ cold velocity distribution'
if(disttype.eq.1)write(6,*)'Done computing initial 3D Gaussian spatial distribution w/ cold velocity distribution'

!allocate the charge density, potential, and field arrays
if(.not.allocated(rho))allocate(rho(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3)))
if(.not.allocated(phi))allocate(phi(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3)))
if(.not.allocated(efield))allocate(efield(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3))
if(.not.allocated(bfield))allocate(bfield(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3))

!determine mesh quantities xmin,xmax,dx and similarly for y,z
if(rectpipe)then
  umin(1)=-0.5d0*apipe; umax(1)=0.5d0*apipe; umin(2)=-0.5d0*bpipe; umax(2)=0.5d0*bpipe; umin(3)=-6.d0*sigma_z; umax(3)=6.d0*sigma_z
  delta(1)=apipe/(nhi(1)-nlo(1)); delta(2)=bpipe/(nhi(2)-nlo(2)); delta(3)=12.d0*sigma_z/(nhi(3)-nlo(3))
else
  call get_mesh_quantities(ptcl(:,1),ptcl(:,3),ptcl(:,5),lostflag,delta,umin,nlo,nhi,nlo,nhi,n1,n_particle,n_particle)
  umax(1)=umin(1)+(nhi(1)-nlo(1))*delta(1)
  umax(2)=umin(2)+(nhi(2)-nlo(2))*delta(2)
  umax(3)=umin(3)+(nhi(3)-nlo(3))*delta(3)
endif
write(6,*)'mesh xmin,xmax=',umin(1),umax(1)
write(6,*)'mesh ymin,ymax=',umin(2),umax(2)
write(6,*)'mesh zmin,zmax=',umin(3),umax(3)

call depose_rho_scalar(ptcl(:,1),ptcl(:,3),ptcl(:,5),lostflag,rho,chrgpermacro,nlo,delta,umin,n_particle,nlo,ifail)
write(6,*)'Done with charge deposition'

  mesh3d%nlo = [nxlo, nylo, nzlo]
  mesh3d%nhi = [nxhi, nyhi, nzhi]
  mesh3d%gamma = gamma
  mesh3d%delta = [delta(1), delta(2), delta(3)]
  mesh3d%min = [umin(1), umin(2), umin(3)]
  mesh3d%max = [umax(1), umax(2), umax(3)]
  mesh3d%npad = [npad(1), npad(2), npad(3)]
  mesh3d%rho = rho
! mesh3d%phi = phi
! mesh3d%efield = efield
! mesh3d%bfield = bfield
if(.not.allocated(mesh3d%phi))allocate(mesh3d%phi(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3)))
if(.not.allocated(mesh3d%efield))allocate(mesh3d%efield(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3))
if(.not.allocated(mesh3d%bfield))allocate(mesh3d%bfield(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3),1:3))

!======= Compute the space-charge fields =======
if(.not.rectpipe)then !FREE SPACE
  print *, 'Space charge field calc with free-space boundary condition...'
  call osc_alloc_freespace_array(nlo,nhi,npad)
!!call osc_freespace_solver(rho,gamma,delta,phi,efield,bfield,nlo,nhi,nlo,nhi,npad,idirectfieldcalc,igfflag)
  call space_charge_freespace(mesh3d, direct_field_calc, integrated_green_function)
endif
if(rectpipe)then      !RECTANGULAR PIPE
  print *, 'Space charge field calc with rectangular pipe boundary condition...'
!read in or compute the transformed rectpipe Green functions for future use
  call osc_alloc_rectpipe_arrays(nlo,nhi,npad)
  if(read_rectpipe)then
    call osc_read_rectpipe_grn
  else
    write(6,*)'Computing rectpipe Green functions and their transforms...'
    call osc_getgrnpipe(gamma,apipe,bpipe,delta,umin,npad)
    write(6,*)'...done computing Green functions and transforms'
    if(write_rectpipe)call osc_write_rectpipe_grn(apipe,bpipe,delta,umin,umax,nlo,nhi,gamma)
  endif
!!call osc_rectpipe_solver(rho,apipe,bpipe,gamma,delta,umin,phi,efield,bfield,nlo,nhi,nlo,nhi,idirectfieldcalc,igfflag)
  call space_charge_rectpipe(mesh3d,apipe,bpipe, direct_field_calc, integrated_green_function)
endif
print *, '...done'
!diagnostics:
call prntall(0,n1,n_particle,size(phi,1),size(phi,2),size(phi,3),ptcl,mesh3d%efield,mesh3d%bfield,tval,delta,umin,rectpipe)
end program
