c maxm  = maximum multipole order
c maxa  = maximum number of points around cylinder axis
c maxz  = maximum number of points along cylinder axis
c maxz2 = maximum number of points along cylinder axis,
c           after making field "symmetric"
c maxk  = maximum number of points along k (wave number) axis
c
c Notes:
c  1. The k axis will comprise nfilonK intervals---and thus a total of
c     nfilonK+1 points.  Hence, in order to use the subroutine filonarr,
c     we must make nfilonK even.
c  2. We ought to make maxm a user-definable parameter; or, rather, have
c     a variable 'maxnm' added to the harmonics common block to record
c     the maximum harmonic desired by the user.
c
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      parameter(maxm=0)
      parameter(maxa=601)
      parameter(maxz=2001)
      parameter(maxz2=2*maxz-1)
      parameter(maxk=2001)
      parameter(nfilonA=500)
      parameter(nfilonK=1000)
      parameter(nfilonZ=2001)
      common/zv/zvec(maxz)
      common/fittA/phin(maxa),aia(maxa),bia(maxa),cia(maxa),maxna
      common/fittZ/zn(maxz2),aiz(maxz2),biz(maxz2),ciz(maxz2),maxnz
      common/harmonics/frfourier(2,0:maxm,maxz),                        &
     &                 fphifourier(2,0:maxm,maxz),                      &
     &                 fzfourier(2,0:maxm,maxz)
      common/harmonics2/harm1(2,0:maxm,maxz2),harm2(2,0:maxm,maxz2)
      common/edata/Ez0j(0:3,maxz2),Er0j(0:2,maxz2)
