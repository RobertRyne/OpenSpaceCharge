c names of components of the Master Input File
c 9 types of input:
c 1=comments, 2=beam, 3=elements, 4=lines, 5=lumps,6=loops,
c 7=labor, 8=include 9=constants
c 10th component will be added shortly (wakes)
c 
c Also, kinds of data in the menu, i.e.  9 kinds of "menu elements"
c 1=simple elements (drift, quadrupole,...)
c 2=user-supplied elements (usr1, usr2,...)
c 3=parameter sets (ps1, ps2,...)
c 4=random elements
c 5=random user-supplied elements
c 6=random parameter sets
c 7=simple commands
c 8=advanced commands
c 9=procedures and fitting & optimization
c
c The following is somewhat wasteful of memory (e.g. 75 simple elements
c but far fewer other kinds of data), but I am leaving as-is for simplicity.
c Eventually, instead of ltc(9,nrpmax), there could be 9 1D arrays.
      parameter (nintypes=9)
      parameter(nrpmax=75)
      common/sharp/ling(nintypes)
      character*8 ling
c type-codes of menu elements (note: eventually ltc should be char*16) 
cryne Dec 1, 2002      character*8 ltc
      character*16 ltc
      common/monics/ltc(9,nrpmax)
c nrp=number of "real" (i.e. numeric) parameters
c ncp=number of character*16 parameters
c     that are associated with elements and commands in monics
      integer nrp,ncp,nrpold
      common/nparm/nrp(9,nrpmax),ncp(9,nrpmax),nrpold(9,nrpmax)
cryne 7/28/2002 added nrpold for backward compatability w/ original
c MaryLie input parameters
