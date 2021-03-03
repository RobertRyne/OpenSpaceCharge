c maxlum = maximum number of lumps in /core/
      parameter (maxlum=20)
c storage for lumps
      common/core/thl(monoms,maxlum),tmhl(6,6,maxlum),
     #            dfl(6,monom1,maxlum),rdfl(3,monom1+1,maxlum),
     #            rrjacl(3,3,monom2+1,maxlum),inuse(maxlum)
