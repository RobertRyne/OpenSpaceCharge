      parameter (maxcoils = 100, maxshape = 6 )
      double precision acoil(maxcoils), alcoil(maxcoils)
      double precision zcoil(maxcoils), shape(maxcoils,maxshape)
      double precision glprod(maxcoils)
      integer itype(maxcoils)
      integer ncoil, mcoil(maxcoils)
      double precision ashield, permshield
      common/coils/mcoil, acoil, alcoil, zcoil, shape, glprod, itype
      common/coiln/ztotal, ncoil
      common/shield/ashield, permshield, shflag
      double precision xldr
      common/ldr/ xldr
      common/coilbl/lblcoil(maxcoils)
      character*8 lblcoil
