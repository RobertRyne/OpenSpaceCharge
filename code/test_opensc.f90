program openbc_test
! Code to test the OpenBC package
use mpi
use decomposition_mod
use numerical_distributions_mod
use deposition_mod
use open_spacecharge_mod

use, intrinsic :: iso_fortran_env
implicit none
! Fortran 2008
integer, parameter  :: dp = REAL64


type (domain_decomposition_struct) :: domain
! This code has no input file. Most of the problem parameters are set here:
!integer, parameter :: nx_gbl=1024, ny_gbl=1024, nz_gbl=1024 ! global # of grid points
integer :: nx_gbl, ny_gbl, nz_gbl ! global # of grid points
!     integer, parameter :: n1=7,maxrayp=20000000 ! local particle arrey; (x,gbx,y,gby,z,gbz)(t) ; 7 is lost particle flag
integer :: n1=7,maxrayp ! local particle arrey; (x,gbx,y,gby,z,gbz)(t) ; 7 is lost particle flag
integer, parameter :: idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
integer, parameter :: igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
integer, parameter :: iseed=1234567 ! seed for random # generator
integer  :: disttype  ! =0 for uniform, =1 for Gaussian
real(dp) :: gaussiancutoff !cutoff if a Gaussian initial condition is used
real(dp) :: e_tot ! eV
real(dp) :: gamma0, mc2

real(dp) :: sigma_x, sigma_y, sigma_z

real(dp) :: bunch_charge !charge per bunch; in this case 0.25 nC
real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax !grid quantities
real(dp), parameter :: tval=0.d0 !the time, an argument passed to diagnostic routines
!
! in this version ptcls is (n1,maxrayp) to conform to what particle manager expects; can switch to (maxrayp,n1) in future versions.
real(dp), allocatable, dimension(:,:) :: ptcls            !particle array
real(dp), allocatable, dimension(:,:,:) :: hx,hy,hz,ex,ey,ez  !field arrays
real(dp), allocatable, dimension(:,:,:) :: rho,phi            !charge density, scalar potential
integer :: ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,ilo,ihi,jlo,jhi,klo,khi !global array sizes
integer :: idecomp,npx,npy,npz ! domain decomposition quantities used in the parallel code
!
real(dp) :: dx,dy,dz,t !spatial and temporal quantities
real(dp), dimension(6,6) :: sigmat !initial beam sigma (i.e. 2nd moment) matrix
real(dp), dimension(6) :: cent     !initial beam centroid
real(dp) :: gb0 !gamma*beta
!
integer :: nraysp,n
real(dp) :: chrgpermacro
integer :: mprocs,myrank,ierr
integer :: i,j,k !use for printout at the end of main

integer :: open_status, namelist_file
character(40) :: in_file

logical :: direct_field_calc, integrated_green_function

namelist / opensc_test_params / &
    nx_gbl, ny_gbl, nz_gbl, maxrayp, e_tot, bunch_charge, &
    sigma_x, sigma_y, sigma_z, gaussiancutoff 
   ! direct_field_calc, integrated_green_function

!
call MPI_INIT(ierr)      !initialize MPI
call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)



! Namelist defaults
nx_gbl=64
ny_gbl=64
nz_gbl=512
maxrayp=25000000
e_tot = 250e6

mc2 = 0.510998910d6
bunch_charge=0.25d-9
sigma_x = 0.001
sigma_y = 0.001
sigma_z = 0.0001
gaussiancutoff= 3.5
disttype=1 ! =0 for uniform, =1 for Gaussian
!direct_field_calc = .true.  
!integrated_green_function = .true.


if(myrank.eq.0) call print_title()
if(myrank.eq.0)write(6,*)'running with ',mprocs,' MPI processes'


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
if(myrank.eq.0)write(*, opensc_test_params)
if(myrank.eq.0)print *, '------------------------'

gamma0 = e_tot/mc2

allocate(ptcls(n1,maxrayp))


!
nraysp=1000000    !local # of rays (i.e. particles); note: nraysp must not be a parameter. It changes!
chrgpermacro=bunch_charge/(1.d0*nraysp)/(1.d0*mprocs)
!
! global grid indices:
 ilo_rho_gbl=1; ihi_rho_gbl=nx_gbl; jlo_rho_gbl=1; jhi_rho_gbl=ny_gbl; klo_rho_gbl=1; khi_rho_gbl=nz_gbl !assumes power-of-2
 if(myrank.eq.0)then
   write(6,*)'global grid imin,imax=',ilo_rho_gbl,ihi_rho_gbl
   write(6,*)'global grid jmin,jmax=',jlo_rho_gbl,jhi_rho_gbl
   write(6,*)'global grid kmin,kmax=',klo_rho_gbl,khi_rho_gbl
 endif
 
! idecomp for this test problem:
idecomp=3
if(myrank.eq.0)write(6,*)'grid domain decomposition: idecomp=',idecomp
! determine local grid indices:
if(idecomp.lt.-1.or.idecomp.gt.6)then
  if(myrank.eq.0)write(6,*)'error: idecomp=',idecomp
  stop
endif
!
! If idecomp=0,...,6 then idecomp defines the type of domain decomposition:  0=xyz, 1=xy, 2=yz, 3=xz, 4=x, 5=y, 6=z,
! and the code will determine npx,npy,npz.
! But if idecomp=-1 then the user has specified npx,npy,npz on input and the following should be skipped
if(idecomp.ge.0.and.idecomp.le.6)then
  call procgriddecomp(mprocs,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz)
  if(myrank.eq.0)then
    write(6,*)'returned from procgriddecomp'
    write(6,*)'# of processes in the x-direction =',npx
    write(6,*)'# of processes in the y-direction =',npy
    write(6,*)'# of processes in the z-direction =',npz
  endif
  if(npx*npy*npz.ne.mprocs)then
    if(myrank.eq.0)write(6,*)'error: npx*npy*npz must equal mprocs but in this case npx,npy,npz,mprocs=',npx,npy,npz,mprocs
    stop
  endif
  if(mod(nx_gbl,npx).ne.0)then
    if(myrank.eq.0)write(6,*)'error: global nx must be a multiple of npx'
    stop
  endif
  if(mod(ny_gbl,npy).ne.0)then
    if(myrank.eq.0)write(6,*)'error: global ny must be a multiple of npy'
    stop
  endif
  if(mod(nz_gbl,npz).ne.0)then
    if(myrank.eq.0)write(6,*)'error: global nz must be a multiple of npz'
    stop
  endif
endif
!
! Important: the openBC routine assume that local arrays have indices that are subsets of the global indices.
! In other words, local arrays are NOT dimensioned (1:nxlcl,1:nylcl,1:nzlcl).
! The user must supply, as arguments to the routine, the starting and ending indices of local arrays in the GLOBAL system.
! The following call to routine decompose returns the required quantities, namely, ilo,ihi,jlo,jhi,klo,khi
call decompose(myrank,mprocs,&
     ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,ilo,ihi,jlo,jhi,klo,khi)
!

! Test
domain%rank = myrank
domain%max_rank = mprocs
domain%idecomp = idecomp
domain%global%lo = [ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl]
domain%global%hi = [ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl]

call init_domain_decomposition(domain, MPI_COMM_WORLD)
do i=0, mprocs
  if(myrank.eq.0)print *, '------------------'
  domain%rank = i
  call init_domain_decomposition(domain, MPI_COMM_WORLD)
 ! if(myrank.eq.0)call print_domain_decomposition(domain)
end do

! Allocate the rho and phi and e and h arrays:
allocate(rho(ilo:ihi,jlo:jhi,klo:khi))
allocate(phi(ilo:ihi,jlo:jhi,klo:khi))
allocate(ex(ilo:ihi,jlo:jhi,klo:khi))
allocate(ey(ilo:ihi,jlo:jhi,klo:khi))
allocate(ez(ilo:ihi,jlo:jhi,klo:khi))
allocate(hx(ilo:ihi,jlo:jhi,klo:khi))
allocate(hy(ilo:ihi,jlo:jhi,klo:khi))
allocate(hz(ilo:ihi,jlo:jhi,klo:khi))
!
! initialize the particles:
if(myrank.eq.0)write(6,*)'generating numerical particle distribution'
if(myrank.eq.0)write(6,*)'particles/proc, total # of particles=',nraysp,nraysp*mprocs
sigmat(1:6,1:6)=0.d0
!!!!  sigmat(1,1)=1.d-6; sigmat(3,3)=1.d-6; sigmat(5,5)=(1.d-3/gamma0)**2 !2nd moment matrix of initial dist
!!!!  sigmat(1,1)=1.d-6; sigmat(3,3)=1.d-6; sigmat(5,5)=(1.d-3/490.237834186378d0)**2 !2nd moment matrix of initial dist
sigmat(1,1)=(sigma_x)**2; sigmat(3,3)=(sigma_y)**2; sigmat(5,5)=(sigma_z)**2 !2nd moment matrix of initial dist
sigmat(2,2)=(1.d-9)**2; sigmat(4,4)=(1.d-9)**2; sigmat(6,6)=(1.d-9)**2 !non-zero emittance for Cholesky
cent(1:6)=0.d0
gb0=sqrt((gamma0+1.d0)*(gamma0-1.d0))
if(myrank.eq.0)write(6,*)'gamma0,beta0=',gamma0,gb0/gamma0
cent(6)=gb0  !if integrating in time, this is gamma*beta of the centroid (set to gamma0 if integrating in z)
call gendist(ptcls,n1,maxrayp,nraysp,sigmat,gaussiancutoff,disttype,iseed)
ptcls(6,1:nraysp)=ptcls(6,1:nraysp)+cent(6) 
ptcls(7,1:maxrayp)=1.d0 !fill the empty array with lost particles
ptcls(7,1:nraysp)=0.d0  !lost particle flag (or use for some other purpose)
if(myrank.eq.0)then
  if(disttype.eq.0)write(6,*)'done computing initial 3D uniform spatial distribution w/ cold velocity dist'
  if(disttype.eq.1)write(6,*)'done computing initial 3D Gaussian spatial distribution w/ cold velocity dist'
endif
!     do n=1,10000
!       write(10000+myrank,'(3(1pe12.5,1x))')ptcls(1,n),ptcls(3,n),ptcls(5,n)
!     enddo
!     flush(10000+myrank)
!
! set grid min/max/dx and compute the charge density:
if(myrank.eq.0)write(6,*)'computing grid quantities and depositing charge...'
call get_rho_and_mesh_spacing(ptcls,rho,dx,dy,dz,xmin,ymin,zmin,chrgpermacro,ilo,ihi,jlo,jhi,klo,khi, &
       ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl,idecomp,npx,npy,npz,n1,nraysp,maxrayp)
!
if(myrank.eq.0)then
  write(6,*)'... done. grid quantities are:'
  write(6,'("xmin,xmax,dx=",3(1pe15.8,1x))')xmin,xmin+(nx_gbl-1)*dx,dx
  write(6,'("ymin,ymax,dy=",3(1pe15.8,1x))')ymin,ymin+(ny_gbl-1)*dy,dy
  write(6,'("zmin,zmax,dz=",3(1pe15.8,1x))')zmin,zmin+(nz_gbl-1)*dz,dz
endif
!
! call the OpenBC solver to compute the potential or fields:
call getfields(rho,gamma0,dx,dy,dz,phi,ex,ey,ez,hx,hy,hz,nx_gbl,ny_gbl,nz_gbl,idecomp,npx,npy,npz,igfflag,idirectfieldcalc, &
                ilo,ihi,jlo,jhi,klo,khi,ilo_rho_gbl,ihi_rho_gbl,jlo_rho_gbl,jhi_rho_gbl,klo_rho_gbl,khi_rho_gbl)
if(myrank.eq.0)write(6,*)'done with computation. writing results'
!
! diagnostics:
!     call prntall(0,n1,nraysp,nx,ny,nz,ptcls,hx,hy,hz,ex,ey,ez,tval,dx,dy,dz,xmin,ymin,zmin)
! As a simple diagnostic, write the potential along a line. Each proc write to a separate file. cat these together to plot.
! phi(i,jfixed,kfixed):
if(mprocs.gt.8192)stop !prevents writing a huge # of files
do k=lbound(phi,3),ubound(phi,3)
  do j=lbound(phi,2),ubound(phi,2)
    do i=lbound(phi,1),ubound(phi,1)
      if(j.ne.jlo_rho_gbl+(jhi_rho_gbl-jlo_rho_gbl+1)/2-1)cycle
      if(k.ne.klo_rho_gbl+(khi_rho_gbl-klo_rho_gbl+1)/2-1)cycle
      if(idirectfieldcalc.eq.0)then
        write(1000+myrank,'(5(1pe14.7,1x))')xmin+(i-ilo_rho_gbl)*dx,ymin+(j-jlo_rho_gbl)*dy,zmin+(k-klo_rho_gbl)*dz,phi(i,j,k)
      endif
      if(idirectfieldcalc.eq.1)then
        write(1000+myrank,'(5(1pe14.7,1x))')xmin+(i-ilo_rho_gbl)*dx,ymin+(j-jlo_rho_gbl)*dy,zmin+(k-klo_rho_gbl)*dz,ex(i,j,k)
      endif
    enddo
  enddo
enddo

    !!!  flush(1000+myrank)
! phi(ifixed,jfixed,k):
do k=lbound(phi,3),ubound(phi,3)
  do j=lbound(phi,2),ubound(phi,2)
    do i=lbound(phi,1),ubound(phi,1)
      if(i.ne.ilo_rho_gbl+(ihi_rho_gbl-ilo_rho_gbl+1)/2-1)cycle
      if(j.ne.jlo_rho_gbl+(jhi_rho_gbl-jlo_rho_gbl+1)/2-1)cycle
      if(idirectfieldcalc.eq.0)then
        write(6000+myrank,'(5(1pe14.7,1x))')xmin+(i-ilo_rho_gbl)*dx,ymin+(j-jlo_rho_gbl)*dy,zmin+(k-klo_rho_gbl)*dz,phi(i,j,k)
      endif
      if(idirectfieldcalc.eq.1)then
        write(6000+myrank,'(5(1pe14.7,1x))')xmin+(i-ilo_rho_gbl)*dx,ymin+(j-jlo_rho_gbl)*dy,zmin+(k-klo_rho_gbl)*dz,ez(i,j,k)
      endif
    enddo
  enddo
enddo
   !!!   flush(6000+myrank)
!

if(myrank.eq.0)write(6,*)'job finished'
call MPI_FINALIZE(ierr)



contains

subroutine print_title()
write(*,*)""
write(*,*)"   ___                     ____                          ____ _                           "
write(*,*)"  / _ \ _ __   ___ _ __   / ___| _ __   __ _  ___ ___   / ___| |__   __ _ _ __ __ _  ___  "
write(*,*)" | | | | '_ \ / _ \ '_ \  \___ \| '_ \ / _` |/ __/ _ \ | |   | '_ \ / _` | '__/ _` |/ _ \ "
write(*,*)" | |_| | |_) |  __/ | | |  ___) | |_) | (_| | (_|  __/ | |___| | | | (_| | | | (_| |  __/ "
write(*,*)"  \___/| .__/ \___|_| |_| |____/| .__/ \__,_|\___\___|  \____|_| |_|\__,_|_|  \__, |\___| "
write(*,*)"       |_|                      |_|                                           |___/       "
write(*,*)''
end subroutine

end
