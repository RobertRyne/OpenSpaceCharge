!***********************************************************************
!
! constants_mod: This file contains definitions of the modules
!                local, math_consts, and phys_consts
!
! Version: 0.1
! Author: D.T.Abell, Tech-X Corp., Jan.2005
!
! Comments
!   This should have been written a long time ago!
!
!   13.Jan.2005 (DTA): The _right_ way to do this would be to have
!   constants such as pi computed (eg., 2.d0*asin(1.d0)), so they're
!   guaranteed to be machine precision.  This, however, is a tough
!   problem for compilers (which would have to emulate an arbitrary
!   computer's floating point intrinsics).  But since we don't want a
!   user changing any of these constants (else why call them so?), we
!   want them declared parameters, which means the compiler must know
!   about them NOW.  The values were generated through 21 digits using
!   Mathematica.
!
!***********************************************************************
!***********************************************************************
!
! math_consts
!
! Description: This module defines some standard math constants.
!
!***********************************************************************
!
      module math_consts
        implicit none
!
! data
!
        double precision, parameter :: pi =     3.14159265358979323846d0
        double precision, parameter :: twopi =  6.28318530717958647693d0
        double precision, parameter :: fourpi = 12.5663706143591729539d0
        double precision, parameter :: euler_e =                        &
     &                                          2.71828182845904523536d0
        double precision, parameter :: golden_ratio =                   &
     &                                          1.61803398874989484820d0

!   conversion factors between radians and degrees
        double precision, parameter :: deg_per_rad =                    &
     &                                          57.2957795130823208768d0
        double precision, parameter :: rad_per_deg =                    &
     &                                        0.0174532925199432957692d0

      end module math_consts
!
!***********************************************************************
!
! phys_consts
!
! Description: This module defines the standard physical constants.  For
! each constant, the comment states the units and (in parentheses) the
! expected error in the last digits.  For more information, see the
! NIST web site: http://physics.nist.gov/cuu/index.html
! NB: The units are mostly SI, but some constants are given in alternate
! units.  Two examples: (1) 'electron_mass' holds the value in kg, but
! 'electron_restE' holds the value in eV; (2) 'planck_h' holds the value
! in J.s, but 'planck_h_ev' holds the value in eV.s.
!
!***********************************************************************
!
      module phys_consts
        implicit none
!
! data
!
! speed of light in vacuum /m.s^-1  (exact)
        double precision, parameter :: c_light = 299792458.0d0

! permittivity of free space /F.m^-1  (exact to digits shown)
        double precision, parameter :: epsilon_o =                      &
     &                                        8.85418781762038985054d-12

! 4.pi.{permittivity of free space} /F.m^-1  (exact to digits shown)
        double precision, parameter :: epsilon_o_4pi =                  &
     &                                        1.11265005605361843217d-10

! permeability of free space /H.m^-1  (exact to digits shown)
        double precision, parameter :: mu_o = 12.5663706143591729539d-7

! {permeability of free space}/4.pi /H.m^-1  (exact)
        double precision, parameter :: mu_o_ovr_4pi = 1.0d-7

! elementary charge /C  (63)
        double precision, parameter :: elem_charge = 1.602176462d-19

! elementary charge /esu  (19)
        double precision, parameter :: elem_charge_esu = 4.80320420d-10

! alpha particle mass /kg  (52)
        double precision, parameter :: alpha_mass = 6.64465598d-27

! alpha particle rest energy /eV  (15)
        double precision, parameter :: alpha_restE = 3727.37904d+6

! electron mass /kg  (72)
        double precision, parameter :: electron_mass = 9.10938188d-31

! electron rest energy /eV  (21)
        double precision, parameter :: electron_restE = 0.510998902d+6

! muon mass /kg  (16)
        double precision, parameter :: muon_mass = 1.88353109d-28

! muon rest energy /eV  (52)
        double precision, parameter :: muon_restE = 105.6583568d+6

! neutron mass /kg  (13)
        double precision, parameter :: neutron_mass = 1.67492716d-27

! neutron rest energy /eV  (38)
        double precision, parameter :: neutron_restE = 939.565330d+6

! proton mass /kg  (13)
        double precision, parameter :: proton_mass = 1.67262158d-27

! proton rest energy /eV  (38)
        double precision, parameter :: proton_restE = 938.271998d+6

! fine-structure constant  (27)
        double precision, parameter :: fine_structure_a = 7.297352533d-3

! fine-structure constant  (50)
        double precision, parameter :: inv_fine_structure_a =           &
     &                                                    137.03599976d0

! Planck constant /J.s  (52)
        double precision, parameter :: planck_h = 6.62606876d-34

! Planck constant /eV.s  (16)
        double precision, parameter :: planck_h_ev = 4.13566727d-15

! {Planck constant}/(2.pi) /J.s  (82)
        double precision, parameter :: planck_h_bar = 1.054571596d-34

! {Planck constant}/(2.pi) /eV.s  (26)
        double precision, parameter :: planck_hbar_ev = 6.58211889d-16

! Boltzmann constant /J.K^-1  (24)
        double precision, parameter :: boltzmann_k = 1.3806503d-23

! Boltzmann constant /eV.K^-1  (15)
        double precision, parameter :: boltzmann_k_ev = 8.617342d-5

      end module phys_consts
!
