      parameter ( maxdat = 1000)
      parameter (maxv = 10, maxq = (maxv+1)*(maxv+2)/2 )
      common /minvar/ jtyp, modu, ntot, nomin, minpt, maxpt, fmin, fmax,
     * regold, region, seed, oldmu, dmu, xdat(maxv,maxdat), fdat(maxdat)
      save /minvar/
