c---------------------------------------------------------------------
c           variable definitions
c
      parameter (maxv=100)
      common/xvary/nva(3),nda(3),ivar,jdep,mm(maxv,3),idx(maxv,3),
     * idv(maxv,3),xbase(maxv,3),xslope(maxv,3)
      character*12 varnam
      common/vnames/varnam(maxv,3)
