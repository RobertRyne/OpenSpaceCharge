# OpenSpaceCharge

The OpenSC package is an open-source software library written primarily in Fortran 2008 for calculating space charge fields [1]. It was originally developed as a reusable Poisson solver with free-space boundary conditions for use within the Warp framework [2], and will be incorporated into the Particle-In-Cell Scalable Application Resource, PICSAR [3].

OpenSC currently implements free-space and rectangular conducting pipe methods using integrated Green functions (IGFs) as described in [4] and [5], respectively. The package provides high-level routines to:

- Deposit weighted charged particles on a 3Drectangular grid.

- Calculate the space charge fields on this grid (various methods).

- Interpolate the field to an arbitrary point within its domain.

[1] C.E.Mayes, R.D.Ryne, D.C.Sagan, "3D Space Charge in Bmad", in Proceedings of IPAC2018, Vancouver, BC, Canada http://accelconf.web.cern.ch/AccelConf/ipac2018/papers/thpak085.pdf

[2] http://warp.lbl.gov/

[3] H. Vincenti, M. Lobet, R. Lehe, R. Sasanka and J-L Vay, "An efficient and portable SIMD algorithm for charge/current deposition in Particle-In-Cell codes," Comp. Phys. Comm. 210, 145-154 (2017) https://picsar.net/

[4] J.Qiang,S.Lidia,R.D.Ryne,andC.Limborg-Deprey,"Three-dimensional quasi-static model for high brightness beam Dynamics simulation," Phys. Rev. ST Accel. Beams, Vol. 9, 044204 (2006)

[5] J.Qiang,S.Lidia,R.D.Ryne,andC.Limborg-Deprey,"Three-dimensional quasi-static model for high brightness beam Dynamics simulation," Phys. Rev. ST Accel. Beams, Vol. 9, 044204 (2006)

## Compilation

To compile:

mkdir build

cd build

cmake ../

make

Run test program:

./test_opensc

or:

mpirun -n 8 ./test_opensc
