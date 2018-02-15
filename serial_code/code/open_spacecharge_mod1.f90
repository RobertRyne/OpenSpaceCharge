module open_spacecharge_mod1
use, intrinsic :: iso_fortran_env
use open_spacecharge_mod
implicit none
! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

type mesh3d_struct
  integer :: nlo(3) = [ 1,  1,  1]       ! Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: nhi(3) = [64, 64, 64]       ! Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: npad(3) = [ 1,  1,  1]      ! Array padding for cyclic convolution
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma
  real(dp), allocatable, dimension(:,:,:) :: rho            ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: phi            ! electric potential grid
  real(dp), allocatable, dimension(:,:,:,:) :: efield
  real(dp), allocatable, dimension(:,:,:,:) :: bfield
end type

contains


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
subroutine print_mesh3d(mesh3d)
type(mesh3d_struct) :: mesh3d
print *, '------------------------'
print *, 'Mesh: '
print *, 'nlo: ', mesh3d%nlo
print *, 'nhi: ', mesh3d%nhi
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

subroutine space_charge_freespace(mesh3d, direct_field_calc, integrated_green_function)
type(mesh3d_struct) :: mesh3d
integer :: idirectfieldcalc=1
integer :: igfflag=1
logical, optional :: direct_field_calc, integrated_green_function

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0

call osc_freespace_solver(mesh3d%rho, mesh3d%gamma, &
  mesh3d%delta(1), mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, mesh3d%npad, idirectfieldcalc,igfflag)
end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine space_charge_rectpipe(mesh3d, apipe, bpipe, direct_field_calc, integrated_green_function)
type(mesh3d_struct) :: mesh3d
real(dp) :: apipe, bpipe
integer :: idirectfieldcalc=1
integer :: igfflag=1
logical, optional :: direct_field_calc, integrated_green_function

idirectfieldcalc=1 ! =0 to compute phi and use finite differences for E; =1 to compute E directly
if(present(direct_field_calc) .and. .not. direct_field_calc)idirectfieldcalc=0

igfflag=1 ! =0 for ordinary Green function; =1 for integrated Green function
if(present(integrated_green_function) .and. .not. integrated_green_function) igfflag = 0

call osc_rectpipe_solver(mesh3d%rho, apipe, bpipe, mesh3d%gamma, &
  mesh3d%delta(1), mesh3d%min, mesh3d%phi, mesh3d%efield, mesh3d%bfield, &
  mesh3d%nlo, mesh3d%nhi, mesh3d%nlo, mesh3d%nhi, idirectfieldcalc,igfflag)

end subroutine

end module
