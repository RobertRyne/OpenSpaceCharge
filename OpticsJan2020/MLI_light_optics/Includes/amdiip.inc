c--------------------------------------------------------------
c   Adaptive Multidimensional Inverse Interpolation (AMDII) Pool
c    C. T. Mottershead  LANL AT-3   June 93
c
      parameter (maxv=40, mxvp = maxv + 1)
      common/amdiip/ic,maxcut,ncut,npts,level,mdii,errtol,antmin,antol,
     * delta,slo,fultgt(maxv),partgt(maxv),err(mxvp),xx(maxv,mxvp),
     * ff(maxv,mxvp),ga(maxv*mxvp), xbar(maxv), yy(maxv)
