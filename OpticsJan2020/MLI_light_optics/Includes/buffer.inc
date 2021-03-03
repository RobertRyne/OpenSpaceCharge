      common/buffer/buf1a(monoms),buf2a(monoms),buf3a(monoms),
     #              buf4a(monoms),buf5a(monoms),
     #              buf1m(6,6),buf2m(6,6),buf3m(6,6),
     #              buf4m(6,6),buf5m(6,6)
      dimension bfa(monoms,5),bma(6,6,5)
      equivalence (bfa,buf1a), (bma,buf1m)
