c    storage for #include function in main Marylie input file
c    na1 = real menu count before the includes
c    na2 = real menu count after the includes
c    nb1 = real item count before the includes
c    nb2 = real item count after the includes
c    noble1 = real task count before the includes
c    noble2 = real task count after the includes
c    ninc = number of include files used
c    MIF fragment include file names for #include
      character*80 incfil
      common/mifrag/incfil(32)
      common/mainum/ninc,na1,na2,nb1,nb2,noble1,noble2
