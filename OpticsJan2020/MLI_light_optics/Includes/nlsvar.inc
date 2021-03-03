      parameter (maxdat = 1000)
      parameter (maxv = 10, maxf = 15, maxq = maxv+1 )
      common /nlsvar/ modu, ntot, nomin, minpt, maxpt, nlsi, fmin,
     * fmax, regold, region, seed, oldmu, dmu, xdat(maxv,maxdat),
     * fdat(maxf,maxdat), fmdat(maxdat), fvmin(maxf), xcen(maxv)
      save /nlsvar/
