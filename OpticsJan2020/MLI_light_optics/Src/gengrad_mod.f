!***********************************************************************
!
! gengrad: data module for generalized gradients
!
! Description: This module implements a derived type data container
! for generalized gradients for either electric or magnetic fields.
!
! Version: 0.1
! Author: D.T.Abell, Tech-X Corp., May.2005
!
! Comments
!   19.May.05 DTA: Implementing this so both electric and  magnetic
!     field generalized gradients can use the same data structure.
!
!***********************************************************************
!
      module gengrad_data
        implicit none
!
! derived types
!
      type gengrad
        integer :: nz_intrvl, maxj, maxm
        double precision :: zmin, zmax, dz
        double precision :: angfrq
        ! must allocate zvals(0:nz_intrvl+rkxtra)
        !               G1c(0:nz_intrvl+rkxtra,0/1:maxj,0:maxm)
        !               G1s(0:nz_intrvl+rkxtra,0/1:maxj,0:maxm)
        !               ...
        !               G3s(0:nz_intrvl+rkxtra,0/1:maxj,0:maxm)
        ! [ rkxtra = extra z values required by adam11 ]
        double precision, dimension(:), pointer :: zvals
        double precision, dimension(:,:,:), pointer :: G1c,G2c,G3c
        double precision, dimension(:,:,:), pointer :: G1s,G2s,G3s
      end type gengrad

      end module gengrad_data
