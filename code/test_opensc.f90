program openbc_test
! Code to test the OpenBC package

use mpi
use open_spacecharge_mod
use test_mod

implicit none
integer, parameter :: dp = REAL64

type(mesh3d_struct) :: mesh3d
integer :: nx, ny, nz ! # of grid points
integer :: n_particle
integer, parameter :: n1=7 ! (x,gbx,y,gby,z,gbz)(t) ; 7 is lost particle flag (not needed)
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

logical :: direct_field_calc, integrated_green_function

integer :: mprocs,myrank,ierr

namelist / opensc_test_params / &
    nx, ny, nz, n_particle, e_tot, bunch_charge, &
    sigma_x, sigma_y, sigma_z, gaussiancutoff, &
    direct_field_calc, integrated_green_function, &
    idecomp, npx, npy, npz
!
call MPI_INIT(ierr)      !initialize MPI
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

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
disttype=1 ! =0 for uniform, =1 for Gaussian
direct_field_calc = .true.  
integrated_green_function = .true.
idecomp=0
npx=1
npy=1
npz=1
 
!Read namelist
in_file = 'test_opensc.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)
open(newunit=namelist_file, file = in_file, status = 'old', iostat=open_status)
if (open_status /= 0) then
  if(myrank.eq.0)print *, 'Input file missing: ', in_file
  if(myrank.eq.0)print *, 'Using defaults'
else 
  read (namelist_file, nml = opensc_test_params)
  close (namelist_file)
endif

if(myrank.eq.0)print *, '------------------------'
if(myrank.eq.0)WRITE(*, opensc_test_params)
if(myrank.eq.0)print *, '------------------------'

gamma = e_tot/mc2
chrgpermacro=bunch_charge/n_particle

allocate(y(n_particle,n1))


! parameters needed to call openbc routines:
! global sizes:
ilo_rho_gbl=1; ihi_rho_gbl=nx; jlo_rho_gbl=1; jhi_rho_gbl=ny; klo_rho_gbl=1; khi_rho_gbl=nz !assumes power-of-2

! determine local array sizes:
if(idecomp.lt.-1.or.idecomp.gt.6)then
  if(myrank.eq.0)write(6,*)'error: idecomp=',idecomp
  stop
endif

! if idecomp=0,...,6 then idecomp defines the type of domain decomposition:  0=xyz, 1=xy, 2=yz, 3=xz, 4=x, 5=y, 6=z,
! and the code will determine npx,npy,npz. But if idecomp=-1 then the user has specified npx,npy,npz on input
if(idecomp.ge.0.and.idecomp.le.6)then
  call procgriddecomp(mprocs,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz)
  if(npx*npy*npz.ne.mprocs)then
    if(myrank.eq.0)write(6,*)'error: npx*npy*npz must equal mprocs but in this case npx,npy,npz,mprocs=',npx,npy,npz,mprocs
    stop
  endif
endif
if(myrank.eq.0)then
  write(6,*)'running with ',mprocs,' MPI processes'
  write(6,*)'global grid imin,imax=',ilo_rho_gbl,ihi_rho_gbl
  write(6,*)'global grid jmin,jmax=',jlo_rho_gbl,jhi_rho_gbl
  write(6,*)'global grid kmin,kmax=',klo_rho_gbl,khi_rho_gbl
  write(6,*)'grid domain decomposition:'
  write(6,*)'idecomp=',idecomp
  write(6,*)'# of processes in the x-direction =',npx
  write(6,*)'# of processes in the y-direction =',npy
  write(6,*)'# of processes in the z-direction =',npz
endif

! Important: the openBC routine assume that local arrays have indices that are subsets of the global indices.
! In other words, local arrays are NOT dimensioned (1:nxlcl,1:nylcl,1:nzlcl).
! The user must supply, as arguments to the routine, the starting and ending indices of local arrays in the GLOBAL system.
! The following call to routine decompose returns the required quantities, namely, ilo,ihi,jlo,jhi,klo,khi
call decompose(myrank,mprocs,&
     ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,ilo,ihi,jlo,jhi,klo,khi)
 
! initialize the particles:
if(myrank.eq.0)write(6,*)'gamma=',gamma
sigmat(1:6,1:6)=0.d0
sigmat(1,1)=sigma_x**2; sigmat(3,3)=sigma_y**2; sigmat(5,5)=sigma_z**2 !2nd moment matrix of initial dist
cent(1:6)=0.d0
gb0=sqrt((gamma+1.d0)*(gamma-1.d0))
if(myrank.eq.0)write(6,*)'beta0=',gb0/gamma
cent(6)=gb0  !if integrating in time, this is gamma*beta of the centroid (set to gamma if integrating in z)

if(myrank.eq.0)print *, 'Generating distribution...'
call gendist(y,n1,n_particle,sigmat,gaussiancutoff,disttype,iseed)
y(:,7)=0.d0 !lost particle flag (or use for some other purpose)
y(1:n_particle,6)=y(1:n_particle,6)+cent(6)
if(myrank.eq.0)then
  if(disttype.eq.0)write(6,*)'...done computing initial 3D uniform spatial distribution w/ cold velocity distribution'
  if(disttype.eq.1)write(6,*)'...done computing initial 3D Gaussian spatial distribution w/ cold velocity distribution'
endif
!diagnostics to make sure that the seed is different on each core:
if(myrank.eq.0)write(2,'(7(1pe12.5,1x))')y(1,1:7),y(2,1:7),y(3,1:7),y(4,1:7),y(5,1:7)
if(myrank.eq.1)write(3,'(7(1pe12.5,1x))')y(1,1:7),y(2,1:7),y(3,1:7),y(4,1:7),y(5,1:7)
if(myrank.eq.2)write(4,'(7(1pe12.5,1x))')y(1,1:7),y(2,1:7),y(3,1:7),y(4,1:7),y(5,1:7)
if(npx.ne.12345)stop !debuggin, so stop here
!

if(myrank.eq.0)print *, 'depositing bunch on mesh...'
mesh3d%n = [nx, ny, nz]
mesh3d%gamma = gamma
call deposit_bunch_on_mesh(y(:,1),y(:,3), y(:,5), y(:,7), bunch_charge, mesh3d)
if(myrank.eq.0)print *, 'space charge field calc...'
call space_charge_field_calc(mesh3d, direct_field_calc=direct_field_calc, integrated_green_function=integrated_green_function)
if(myrank.eq.0)print *, '...done'


! diagnostics:
if(myrank.eq.0)then
  call prntall(0,n1,n_particle,nx,ny,nz,y, &
  mesh3d%field%B(1), mesh3d%field%B(2), mesh3d%field%B(3), &
  mesh3d%field%E(1), mesh3d%field%E(2), mesh3d%field%E(3), &
  tval, &
  mesh3d%delta(1),  mesh3d%delta(2),  mesh3d%delta(3), &
  mesh3d%min(1), mesh3d%min(2), mesh3d%min(3)) 
endif
  
!=========================================
!
call MPI_FINALIZE(ierr)

end program

