! Test code for the OpenSC space-charge package
! R.D. Ryne and C. Mayes, February 2018
program opensc_test
use open_spacecharge_mod
use open_spacecharge_core_mod ! contains routines that begin with osc_ and some global variables that must persist for thread safety
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
real(dp) :: xcent,ycent,zcent

!geometry
real(dp) :: apipe, bpipe !full width and height of rectangular pipe; ignored unless rectpipe is .true.

!Poisson solver parameters
logical :: direct_field_calc, integrated_green_function
logical :: rectpipe, read_rectpipe, write_rectpipe
integer :: igfflag,idirectfieldcalc
logical :: cathode_images
integer :: image_method

real(dp) :: t1, t2 ! For timing
real(dp) :: tval=0.d0 !the time, an argument passed to diagnostic routines, irrelevant in this test code
integer :: ifail !return flag for charge deposition routine

integer :: open_status, namelist_file


character(40) :: in_file
integer :: nxlo,nxhi,nylo,nyhi,nzlo,nzhi
namelist / opensc_test_params / &
  nxlo,nxhi,nylo,nyhi,nzlo,nzhi, n_particle, e_tot, bunch_charge, sigma_x, sigma_y, sigma_z, gaussiancutoff, &
  direct_field_calc, integrated_green_function, cathode_images, image_method,&
  rectpipe, read_rectpipe, write_rectpipe, apipe, bpipe ! if rectangular pipe BC is being used, apipe=full width, bpipe=full height

!Namelist defaults
nxlo=1;  nylo=1;  nzlo=1
nxhi=64; nyhi=64; nzhi=64
e_tot = 5.d6 ! 250.d6    !5 MeV or 250 MeV
n_particle = 10000000
mc2 = 0.510998910d6
bunch_charge=0.25d-9
sigma_x = 0.001d0; sigma_y = 0.001d0; sigma_z = 0.0001d0
gaussiancutoff=4
disttype=0 ! =0 for uniform, =1 for Gaussian
direct_field_calc = .true. ! .false.
integrated_green_function = .true.
cathode_images=.true.
image_method=2 ! =1 for convolution/correlation, =2 for shifted Green function 
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
nlo=[nxlo,nylo,nzlo]
nhi=[nxhi,nyhi,nzhi]
npad=[ipad,jpad,kpad]

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
if(cathode_images .and. disttype.eq.0)cent(5)=1.2d-4 !image charge test with uniform bunch
if(cathode_images .and. disttype.eq.1)cent(5)=5.0d-4 !image charge test with Gaussian bunch
!if(cathode_images)write(6,*)'cent(5) set to ',cent(5)
write(6,*)'beta0=',gambet/gamma
cent(6)=gambet  !if integrating in time, this is gamma*beta of the centroid (set to gamma if integrating in z)

call gendist(ptcl,n1,n_particle,sigmat,gaussiancutoff,disttype,iseed,apipe,bpipe,rectpipe)
ptcl(:,7)=0.d0 !lost particle flag (or use for some other purpose)
ptcl(1:n_particle,5)=ptcl(1:n_particle,5)+cent(5)
write(6,*)'added zcentroid to z particle data, where zcentroid=',cent(5) 
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
  if(cathode_images)then
    call get_mesh_quantities(ptcl(:,1),ptcl(:,3),ptcl(:,5),lostflag,delta,umin,nlo,nhi,nlo,nhi,n1,n_particle,n_particle,1) !set umin(3)=1.d-9
  else
    call get_mesh_quantities(ptcl(:,1),ptcl(:,3),ptcl(:,5),lostflag,delta,umin,nlo,nhi,nlo,nhi,n1,n_particle,n_particle,0)
  endif
  umax(1)=umin(1)+(nhi(1)-nlo(1))*delta(1)
  umax(2)=umin(2)+(nhi(2)-nlo(2))*delta(2)
  umax(3)=umin(3)+(nhi(3)-nlo(3))*delta(3)
endif
write(6,*)'mesh xmin,xmax=',umin(1),umax(1)
write(6,*)'mesh ymin,ymax=',umin(2),umax(2)
write(6,*)'mesh zmin,zmax=',umin(3),umax(3)
write(6,*)'delta(1:3)=',delta(1:3)

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

call cpu_time(t1)
!======= Compute the space-charge fields =======
if(.not.rectpipe.and..not.cathode_images)then !FREE SPACE
  print *, 'Space charge field calc with free-space boundary condition...'
!!call osc_freespace_solver(rho,gamma,delta,phi,efield,bfield,nlo,nhi,nlo,nhi,npad,idirectfieldcalc,igfflag)
  call space_charge_freespace(mesh3d, direct_field_calc, integrated_green_function)
endif
if(.not.rectpipe.and.cathode_images)then !FREE SPACE BUT WITH CATHODE IMAGES
  print *, 'Space charge field calc with cathode images...'
  if(umin(3).lt.0.d0)write(6,*)'error: umin(3) is less than zcathode!'
  if(umin(3).lt.0.d0)stop
  call space_charge_cathodeimages(mesh3d, direct_field_calc, integrated_green_function, image_method)
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
call cpu_time(t2)
print *, 'Time for space charge calc (s): ', t2-t1


! New diagnostics
xcent=sum(ptcl(1:n_particle,1))/(1.d0*n_particle)
ycent=sum(ptcl(1:n_particle,3))/(1.d0*n_particle)
zcent=sum(ptcl(1:n_particle,5))/(1.d0*n_particle)
call write_lines(mesh3d,xcent,ycent,zcent)
call write_plane(mesh3d)

stop 

!old diagnostics:
call prntall(0,n1,n_particle,size(phi,1),size(phi,2),size(phi,3),ptcl,mesh3d%efield,mesh3d%bfield,tval,delta,umin,rectpipe)

contains

subroutine write_lines(mesh3d,xcent,ycent,zcent)
type(mesh3d_struct) :: mesh3d
real(dp) :: x, y, z, Evec(3), xcent,ycent,zcent 
integer :: i, outfile

open(newunit=outfile, file = 'x_lineout.dat')
z = zcent
y = ycent 
do i = mesh3d%nlo(1), mesh3d%nhi(1) -1 ! skip last point
  x = (i-1)*mesh3d%delta(1) + mesh3d%min(1) 
  call interpolate_field(x, y, z, mesh3d, E=Evec)
! write(outfile, *) x, Evec(1)
  write(outfile, '(6(1pe14.7,1x))') x,y,z,Evec(1:3)
enddo  
close(outfile) 

open(newunit=outfile, file = 'y_lineout.dat')
x = xcent
z = zcent 
do i = mesh3d%nlo(2), mesh3d%nhi(2) -1 ! skip last point
  y = (i-1)*mesh3d%delta(2) + mesh3d%min(2) 
  call interpolate_field(x, y, z, mesh3d, E=Evec)
! write(outfile, *) y, Evec(2)
  write(outfile, '(6(1pe14.7,1x))') x,y,z,Evec(1:3)
enddo  
close(outfile) 

open(newunit=outfile, file = 'z_lineout.dat')
x = xcent
y = ycent
do i = mesh3d%nlo(3), mesh3d%nhi(3) -1 ! skip last point
  z = (i-1)*mesh3d%delta(3) + mesh3d%min(3) 
  call interpolate_field(x, y, z, mesh3d, E=Evec)
! write(outfile, *) z, Evec(3)
  write(outfile, '(6(1pe14.7,1x))') x,y,z,Evec(1:3)
enddo  
close(outfile)  
end subroutine


 
subroutine write_plane(mesh3d)
type(mesh3d_struct) :: mesh3d
real(dp) :: x, y, z, Evec(3) 
integer :: i, k, outfile

open(newunit=outfile, file = 'x_z_Ex_Ez.dat')
y = 0 
do k = mesh3d%nlo(3), mesh3d%nhi(3) -1 ! skip last point
  z = (k-1)*mesh3d%delta(3) + mesh3d%min(3) 
  do i = mesh3d%nlo(1), mesh3d%nhi(1) -1 ! skip last point
    x = (i-1)*mesh3d%delta(1) + mesh3d%min(1) 
    call interpolate_field(x, y, z, mesh3d, E=Evec)
    write(outfile, *) x, z, Evec(1), Evec(3)
  enddo
enddo  
close(outfile)  
end subroutine


end program



