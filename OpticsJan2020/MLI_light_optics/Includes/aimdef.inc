c-------------------------------------------------------------------
c           aim definitions for inner and outer loops
c
      parameter (maxa=100)
      common/aimdef/maim,nsq(3),ksq(3),lsq(3,3),kyf(maxa,3),
     * nsf(maxa,3),idf(maxa,3),jdf(maxa,3),target(maxa,3),wts(maxa,3)
      character*8 qname
      common/qlabel/qname(maxa,3)
