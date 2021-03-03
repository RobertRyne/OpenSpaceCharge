c maxlum = maximum number of lumps in /core/
      parameter (maxlum=20)
c storage for lumps
c     common/core/thl(monoms,maxlum),tmhl(6,6,maxlum),
c    #            dfl(6,monoms,maxlum),rdfl(3,monom1+1,maxlum),
c    #            rrjacl(3,3,monom2+1,maxlum),inuse(maxlum)
c storage for lumps
      common/core/thl(monoms,maxlum),tmhl(6,6,maxlum),
     #            dfl(6,monoms,maxlum),rdfl(3,monoms,maxlum),
     #            rrjacl(3,3,monoms,maxlum),inuse(maxlum)
