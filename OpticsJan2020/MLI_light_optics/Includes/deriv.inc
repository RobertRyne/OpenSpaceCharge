c additional characteristics of the current map:
      double precision df,rjac,rdf,rrjac
      common/deriv/ df(6,monoms)
      common/rjacob/rjac(3,3,monoms)
c     common/rderiv/rdf(3,84)
c     common/rrjac/ rrjac(3,3,28)
      common/rderiv/rdf(3,monoms)
      common/rrjac/ rrjac(3,3,monoms)
