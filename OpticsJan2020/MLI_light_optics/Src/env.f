      subroutine initenv(p,icorrel)
c R. Ryne 4/7/2004
c routine to initialize rms envelopes.
c this assumes that the units are those used in the current ML/I simulation.
c input parameters array elements:
c p(1)=x_rms,  p(2)=px_rms,  p(3)=xpx_rms,  p(4)=x_emittance_normalized_rms
c p(5)=y_rms,  p(6)=py_rms,  p(7)=ypy_rms,  p(8)=y_emittance_normalized_rms
c p(9)=t_rms,  p(10)=pt_rms, p(11)=tpt_rms, p(12)=t_emittance_normalized_rms
c p(13-15)=canonical momenta conjugate to p(1-3), respectively
c icorrel describes whether xpx, ypy, tpt are scaled (i.e. correlations) or not
c e.g. if icorrel=0, then <x px> is really <x px>
c      if icorrel=1, then <x px> is really <x px>/(xrms*pxrms)?
c      and similary for <y py> and <t pt>
c
c NOTE WELL: This routine assumes that, at the point the initial values
c            are set, the beam is not inside an rf cavity. In this case,
c            canonical mom. conjugate to x_rms = <x px>/xrms
c            canonical mom. conjugate to y_rms = <y py>/yrms
c            canonical mom. conjugate to t_rms = <t pt>/trms
c If the beam is initialized inside an rf cavity, this does not hold.
      use parallel, only : idproc
      include 'impli.inc'
      common/envdata/env(6),envold(6),emap(6,6)
      common/emitdata/emxn2,emyn2,emtn2
      dimension p(*)
      if(p(1).eq.0.d0)then
        if(idproc.eq.0)then
          write(6,*)'error(initenv): must specify rms beam sizes'
        endif
        call myexit
      endif
c
c X-PX:
      env(1)=p(1)
      if(p(13).ne.0.d0)then
c user is specifying cpx (canonical px) and emittance:
        env(2)=p(13)
        emxn2=p(4)**2
      else
c user is specifying px & xpx, px & emit, or xpx & emit:
        if(p(4).eq.0.d0)then   !user specifying px & xpx
          if(icorrel.eq.0)then
            env(2)=p(3)/p(1)
            emxn2=p(1)**2*p(2)**2-p(3)**2
          endif
          if(icorrel.eq.1)then
            env(2)=p(3)*p(2)
            emxn2=p(1)**2*p(2)**2*(1.d0-p(3)**2)
          endif
        else
          if(p(2).ne.0.d0)then !user specifying px & emit
            emxn2=p(4)**2
            xpx2=p(1)**2*p(2)**2-emxn2
            env(2)=sqrt(xpx2)/p(1)
          else                 !user specifying xpx & emit
            env(2)=p(3)/p(1)
            emxn2=p(4)**2
          endif
        endif
      endif
c
c Y-PY:
      env(3)=p(5)
      if(p(14).ne.0.d0)then
c user is specifying cpy (canonical py) and emittance:
        env(4)=p(14)
        emyn2=p(8)**2
      else
c user is specifying py & ypy, py & emit, or ypy & emit:
        if(p(8).eq.0.d0)then   !user specifying py & ypy
          if(icorrel.eq.0)then
            env(4)=p(7)/p(5)
            emyn2=p(5)**2*p(6)**2-p(7)**2
          endif
          if(icorrel.eq.1)then
            env(4)=p(7)*p(6)
            emyn2=p(5)**2*p(6)**2*(1.d0-p(7)**2)
          endif
        else
          if(p(6).ne.0.d0)then !user specifying py & emit
            emyn2=p(8)**2
            ypy2=p(5)**2*p(6)**2-emyn2
            env(4)=sqrt(ypy2)/p(5)
          else                 !user specifying ypy & emit
            env(4)=p(7)/p(5)
            emyn2=p(8)**2
          endif
        endif
      endif
c
c T-PT:
      env(5)=p(9)
      if(p(15).ne.0.d0)then
c user is specifying cpt (canonical pt) and emittance:
        env(6)=p(15)
        emtn2=p(12)**2
      else
c user is specifying pt & tpt, pt & emit, or tpt & emit:
        if(p(12).eq.0.d0)then   !user specifying pt & tpt
          if(icorrel.eq.0)then
            env(6)=p(11)/p(9)
            emtn2=p(9)**2*p(10)**2-p(11)**2
          endif
          if(icorrel.eq.1)then
            env(6)=p(11)*p(10)
            emtn2=p(9)**2*p(10)**2*(1.d0-p(11)**2)
          endif
        else
          if(p(10).ne.0.d0)then !user specifying pt & emit
            emtn2=p(12)**2
            tpt2=p(9)**2*p(10)**2-emtn2
            env(6)=sqrt(tpt2)/p(9)
          else                 !user specifying tpt & emit
            env(6)=p(11)/p(9)
            emtn2=p(12)**2
          endif
        endif
      endif
c initialize the matrix for the envelope map,
c in case the user wants to compute it:
      emap(1:6,1:6)=0.d0
      do i=1,6
        emap(i,i)=1.d0
      enddo
c store the initial values of the envelopes,
c in case the user wants to perform a contraction mapping later:
      envold(1:6)=env(1:6)
      return
      end
c
      subroutine contractenv(delta)
c perform contraction map on envelopes to find fixed point (rms matched beam).
c delta is the residual, i.e. a metric to describe how close the beam is
c to a match.
      use parallel, only : idproc
      use beamdata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
      include 'map.inc'
      include 'setref.inc'
      dimension d(6),v(6),ca(6),amat(6,6)
      common/envdata/env(6),envold(6),emap(6,6)
      dimension augm(6,7)
c    Apply the contraction mapping, C, to the vector a :
c              Ca = a + 1/(I-M) * (a-Na) 
c    where Na == Script m acting on a == numerical integration
c    using a for initial values.
c  
c     if(idproc.eq.0)then
c       write(6,*)'***********************entering contractenv'
c       write(6,*)'***********************entering contractenv, env='
c       write(6,51)env(1:6)
c     endif
c  51 format(6(1x,1pe12.5))
c     
      d(:)=envold(:)-env(:)
c       write(6,*)'***********************d(1:6)='
c       write(6,51)d(1:6)
c       write(6,*)'***********************emap(1:6,1:6)='
c       write(6,51)emap(1:6,1:6)
c  
        amat(:,:)=-emap(:,:)
      do i=1,6
        amat(i,i)=1.+amat(i,i)
      enddo
c  Compute 1/(I-M) * d  by solving (I-M)v=d  for v :
c     write(6,*)'calling leshs, v(1-6)=',v(1:6)
      call leshs(v,6,amat,d,augm,det)
c     write(6,*)'returned from leshs, v(1-6)=',v(1:6)
c  Now compute Ca:
      ca(:)=envold(:)-v(:)
c check for convergence:
      vnorm=sqrt(v(1)**2+v(2)**2+v(3)**2+v(4)**2+v(5)**2+v(6)**2)
      anorm=sqrt(envold(1)**2+envold(2)**2+envold(3)**2+                &
     &           envold(4)**2+envold(5)**2+envold(6)**2)
      delta=vnorm/anorm
c     if(idproc.eq.0)write(6,*)'delta=',delta
c     if(delta.lt.1.d-9)then
c       if(idproc.eq.0)then
c         write(6,*)'SEARCH CONVERGED'
c       endif
c     endif
c continue the contraction map procedure:
      envold(:)=ca(:)
      env(:)=envold(:)
      emap(:,:)=0.d0
      do i=1,6
        emap(i,i)=1.d0
      enddo
ccccc reftraj(1:6)=refsave(1,1:6)
ccccc arclen=arcsave(1)
ccccc brho=brhosav(1)
ccccc gamma=gamsav(1)
ccccc gamm1=gam1sav(1)
ccccc beta=betasav(1)
c
c diagnostic for debug:
c     if(idproc.eq.0)then
c       write(6,*)'matrix:'
c       write(6,51)(emap(1,j),j=1,6)
c       write(6,51)(emap(2,j),j=1,6)
c       write(6,51)(emap(3,j),j=1,6)
c       write(6,51)(emap(4,j),j=1,6)
c       write(6,51)(emap(5,j),j=1,6)
c       write(6,51)(emap(6,j),j=1,6)
c     endif
c       write(6,*)'***********************leaving contractenv'
c       write(6,*)'***********************leaving contractenv, env='
c       write(6,51)env(1:6)
      return
      end
c
      subroutine envtrace(mh)
c "Envelope tracking" routine applies the linear map to the envelopes.
c This is used along with the split-operator symplectic integrator
c to advance the envelope equations. (The space-charge and emittance
c terms are handled in mid-step in subroutine envkick)
c R. Ryne 4/7/2004
      include 'impli.inc'
      double precision mh(6,6)
      double precision vec(6)
      common/envdata/env(6),envold(6),emap(6,6)
c     write(6,*)'here I am in *envtrace*; env(1:6)='
c     write(6,51)env(1:6)
c     write(6,*)'mh(1:6,1:6)='
c     write(6,51)mh(1:6,1:6)
c  51 format(6(1x,1pe12.5))
c initialize zlm
      vec(1:6)=0.d0
c
      do 100 i=1,6
       do 90 j=1,6
        vec(i)=vec(i) + mh(i,j)*env(j)
   90  continue
  100 continue
      env(1:6)=vec(1:6)
c     write(6,*)'new values of env(1),env(2),env(5),env(6)='
c     write(6,*)env(1),env(2),env(5),env(6)
c code to advanced the momentum portion of the "envelope map," i.e. the
c linear map for envelope dynamics around the env reference envelope:
      call mmult(mh,emap,emap)
      return
      end
c
      subroutine envkick(tau)
      use parallel, only : idproc
      use beamdata
ccc   real*8 :: brho,gamma,gamm1,beta,achg,pmass,bfreq,bcurr,c
ccc   real*8 :: sl,p0sc,ts,omegascl,freqscl
      include 'impli.inc'
      dimension etmp(6,6)
      common/envdata/env(6),envold(6),emap(6,6)
      common/emitdata/emxn2,emyn2,emtn2
c     write(6,*)'here I am in envkick; env(1),env(3),env(5)='
c     write(6,*)env(1),env(3),env(5)
c formulas are coded in terms of xl,xp,xw; connect to parameters in acceldata:
      clite=299792458.d0
      pi=2.d0*asin(1.d0)
      fpei=(clite**2)*1.d-7
      q=1.d0
      xmc2=pmass
      xl=sl
      xp=p0sc
      xw=omegascl
c The commented out statement below is probably wrong, because
c bfreq should only affect the total charge.
c The statement makes sense only if w is the scale ang freq.
c     w=2.d0*pi*bfreq
      w=xw
      p0=gamma*beta*pmass/clite
      xmp0=1./(gamma*beta*clite)
c
      uu=env(1)**2
      vv=env(3)**2
      ww=(env(5)*gamma*beta*clite/(w*xl))**2
c     write(6,*)'new values of uu,vv,ww=',uu,vv,ww
cryne note: g311,g131,g113 all vary as ~1/[length**3]; g511 etc vary as 1/l**5
cryne In other words, if env(1),env(3),and env(5) are doubled, then
cryne g311,g131,g113 are reduced by a factor of 8 (i.e. 8 times smaller)
c     write(6,*)'calling scdrd'
      call scdrd(uu,vv,ww,g311,g131,g113,g511,g151,g115,g331,g313,g133)
c     write(6,*)'returned from scdrd'
cryne qtot=bcurr*2.*pi/w
      qtot=bcurr/bfreq
      xi0=xmc2*clite/(q*fpei)
      xlam3=1.d0/(5.d0*sqrt(5.d0))
      qcon=1.5d0*qtot*xlam3*clite/(xi0*(beta*gamma*xl)**2)*p0/xp
      tcon=qcon*(gamma*beta*clite/(w*xl))**2
c
      r11= qcon*g311*env(1)  +emxn2*(xp/(p0*xl))/env(1)**3
      r33= qcon*g131*env(3)  +emyn2*(xp/(p0*xl))/env(3)**3
      r55= tcon*g113*env(5)+emtn2/env(5)**3*(xl*xp*xw*xw*xmp0*xmp0/p0)
c
      env(2)=env(2)+r11*tau
      env(4)=env(4)+r33*tau
      env(6)=env(6)+r55*tau
c     write(6,*)'done updating env; preparing to update emap'
c
c Advance emap, which is transfer map around this envelope:
c Note: there is probably a bug in the original envelope code, fixsc3,
c which refers to the variable w0 (the cavity ang freq) in the formulas below.
c In fixsc3, w0 is set equal to the scale angular frequency.
c In other words, fixsc3 *only* works when the cav ang freq=scale ang freq.
c The two formulas are commented out below and replaced by the correct ones.
c
      xyfac=xp/(p0*xl)
c     xyfac=0.d0
c
      s11=3.*emxn2*xyfac/(uu**2)+qcon*(3.*uu*g511-g311)
      s33=3.*emyn2*xyfac/(vv**2)+qcon*(3.*vv*g151-g131)
c     s55=3.d0*(w0*xmp0*xl)**2*emtn2*xyfac/env(5)**4+                     &
      s55=3.d0*(xw*xmp0*xl)**2*emtn2*xyfac/env(5)**4+                     &
     &    tcon*(3.d0*ww*g115-g113)
      s13=qcon*env(1)*env(3)*g331
      s15=tcon*env(1)*env(5)*g313
      s35=tcon*env(3)*env(5)*g133
      s31=s13
      s51=s15
      s53=s35
      s22=xyfac
      s44=s22
c     s66=(w0*xmp0*xl)**2*s22
      s66=(xw*xmp0*xl)**2*s22
c
      etmp(1:6,1:6)=0.d0
      do i=1,6
        etmp(i,i)=1.d0
      enddo
      etmp(2,1)=-tau*s11
      etmp(4,3)=-tau*s33
      etmp(6,5)=-tau*s55
c
      etmp(2,3)=-tau*s13
      etmp(2,5)=-tau*s15
c
      etmp(4,1)=-tau*s13
      etmp(4,5)=-tau*s35
c
      etmp(6,1)=-tau*s15
      etmp(6,3)=-tau*s35
c
c     write(6,*)'calling mmult'
      call mmult(etmp,emap,emap)  !(left,right,out)
cxxx  call mmult(emap,etmp,emap)
c     write(6,*)'leaving envkick'
      return
      end 

c=====================================================================
      subroutine scdrd(uu,vv,ww,                                          &
     &h311,h131,h113,h511,h151,h115,h331,h313,h133)
      include 'impli.inc'
      alpha=3.d0/(uu+vv+ww)
      a=alpha*uu
      b=alpha*vv
      c=alpha*ww
      fac5=sqrt(alpha**3)
      fac7=sqrt(alpha**5)
      eps=1.d-4
c
      result1=drd(b,c,a,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      result1=2.d0/3.d0*result1
      h311=result1*fac5
c
      result2=drd(a,c,b,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      result2=2.d0/3.d0*result2
      h131=result2*fac5
c
      result3=drd(a,b,c,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      result3=2.d0/3.d0*result3
      h113=result3*fac5
c-----------------------------------------------------
      result1plus=drd(b*(1.d0+eps),c,a,ierr)
      result1minus=drd(b*(1.d0-eps),c,a,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*b)
      h331=-deriv*2.d0/3.d0*fac7*2.d0
c
      result1plus=drd(c*(1.d0+eps),b,a,ierr)
      result1minus=drd(c*(1.d0-eps),b,a,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*c)
      h313=-deriv*2.d0/3.d0*fac7*2.d0
c
      result1plus=drd(b*(1.d0+eps),a,c,ierr)
      result1minus=drd(b*(1.d0-eps),a,c,ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*b)
      h133=-deriv*2.d0/3.d0*fac7*2.d0
c-----------------------------------------------------
      result1plus=drd(b,c,a*(1.d0+eps),ierr)
      result1minus=drd(b,c,a*(1.d0-eps),ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*a)
      h511=-deriv*2.d0/3.d0*fac7*2.d0/3.d0
c
      result1plus=drd(a,c,b*(1.d0+eps),ierr)
      result1minus=drd(a,c,b*(1.d0-eps),ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*b)
      h151=-deriv*2.d0/3.d0*fac7*2.d0/3.d0
c
      result1plus=drd(a,b,c*(1.d0+eps),ierr)
      result1minus=drd(a,b,c*(1.d0-eps),ierr)
      if(ierr.gt.0)write(6,*)'trouble: ierr>0 in routine scdrd'
      deriv=(result1plus-result1minus)/(2.d0*eps*c)
      h115=-deriv*2.d0/3.d0*fac7*2.d0/3.d0
c-----------------------------------------------------
      return
      end
c
c====================================================================
c
      DOUBLE PRECISION FUNCTION DRD(X,Y,Z,IER)                          DRD    3
C***BEGIN PROLOGUE  DRD                                                 DRD    4
C***DATE WRITTEN   790801   (YYMMDD)                                    DRD    5
C***REVISION DATE  861211   (YYMMDD)                                    DRD    6
C***CATEGORY NO.  C14                                                   DRD    7
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(RD-S DRD-D),         DRD    8
C             COMPLETE ELLIPTIC INTEGRAL,DUPLICATION THEOREM,           DRD    9
C             INCOMPLETE ELLIPTIC INTEGRAL,INTEGRAL OF THE SECOND KIND, DRD   10
C             TAYLOR SERIES                                             DRD   11
C***AUTHOR  CARLSON, B.C., AMES LABORATORY-DOE                          DRD   12
C             IOWA STATE UNIVERSITY, AMES, IOWA  50011                  DRD   13
C           NOTIS, E.M., AMES LABORATORY-DOE                            DRD   14
C             IOWA STATE UNIVERSITY, AMES, IOWA  50011                  DRD   15
C           PEXTON, R.L., LAWRENCE LIVERMORE NATIONAL LABORATORY        DRD   16
C             LIVERMORE, CALIFORNIA  94550                              DRD   17
C***PURPOSE  Compute the INCOMPLETE or COMPLETE Elliptic integral of    DRD   18
C            the 2nd kind. For X and Y nonnegative, X+Y and Z positive, DRD   19
C            DRD(X,Y,Z) = Integral from ZERO to INFINITY of             DRD   20
C                                -1/2     -1/2     -3/2                 DRD   21
C                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.              DRD   22
C            If X or Y is zero, the integral is COMPLETE.               DRD   23
C***DESCRIPTION                                                         DRD   24
C                                                                       DRD   25
C   1.     DRD                                                          DRD   26
C          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL       DRD   27
C          of the second kind                                           DRD   28
C          Standard FORTRAN function routine                            DRD   29
C          Double precision version                                     DRD   30
C          The routine calculates an approximation result to            DRD   31
C          DRD(X,Y,Z) = Integral from zero to infinity of               DRD   32
C                              -1/2     -1/2     -3/2                   DRD   33
C                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,                DRD   34
C          where X and Y are nonnegative, X + Y is positive, and Z is   DRD   35
C          positive.  If X or Y is zero, the integral is COMPLETE.      DRD   36
C          The duplication theorem is iterated until the variables are  DRD   37
C          nearly equal, and the function is then expanded in Taylor    DRD   38
C          series to fifth order.                                       DRD   39
C                                                                       DRD   40
C   2.     Calling Sequence                                             DRD   41
C                                                                       DRD   42
C          DRD( X, Y, Z, IER )                                          DRD   43
C                                                                       DRD   44
C          Parameters On Entry                                          DRD   45
C          Values assigned by the calling routine                       DRD   46
C                                                                       DRD   47
C          X      - Double precision,nonnegative variable               DRD   48
C                                                                       DRD   49
C          Y      - Double precision,nonnegative variable               DRD   50
C                                                                       DRD   51
C                   X + Y is positive                                   DRD   52
C                                                                       DRD   53
C          Z      - Double precision,positive variable                  DRD   54
C                                                                       DRD   55
C                                                                       DRD   56
C                                                                       DRD   57
C          On Return    (values assigned by the DRD routine)            DRD   58
C                                                                       DRD   59
C          DRD     - Double precision approximation to the integral     DRD   60
C                                                                       DRD   61
C                                                                       DRD   62
C          IER    - Integer                                             DRD   63
C                                                                       DRD   64
C                   IER = 0 Normal and reliable termination of the      DRD   65
C                           routine. It is assumed that the requested   DRD   66
C                           accuracy has been achieved.                 DRD   67
C                                                                       DRD   68
C                   IER >  0 Abnormal termination of the routine        DRD   69
C                                                                       DRD   70
C                                                                       DRD   71
C          X, Y, Z are unaltered.                                       DRD   72
C                                                                       DRD   73
C   3.    Error Messages                                                DRD   74
C                                                                       DRD   75
C         Value of IER assigned by the DRD routine                      DRD   76
C                                                                       DRD   77
C                  Value assigned         Error message printed         DRD   78
C                  IER = 1                DMIN1(X,Y) .LT. 0.0D0         DRD   79
C                      = 2                DMIN1(X + Y, Z ) .LT. LOLIM   DRD   80
C                      = 3                DMAX1(X,Y,Z) .GT. UPLIM       DRD   81
C                                                                       DRD   82
C                                                                       DRD   83
C   4.     Control Parameters                                           DRD   84
C                                                                       DRD   85
C                  Values of LOLIM,UPLIM,and ERRTOL are set by the      DRD   86
C                  routine.                                             DRD   87
C                                                                       DRD   88
C          LOLIM and UPLIM determine the valid range of X, Y, and Z     DRD   89
C                                                                       DRD   90
C          LOLIM  - Lower limit of valid arguments                      DRD   91
C                                                                       DRD   92
C                    Not less  than 2 / (machine maximum) ** (2/3).     DRD   93
C                                                                       DRD   94
C          UPLIM  - Upper limit of valid arguments                      DRD   95
C                                                                       DRD   96
C                 Not greater than (0.1D0 * ERRTOL / machine            DRD   97
C                 minimum) ** (2/3), where ERRTOL is described below.   DRD   98
C                 In the following table it is assumed that ERRTOL will DRD   99
C                 never be chosen smaller than 1.0D-5.                  DRD  100
C                                                                       DRD  101
C                                                                       DRD  102
C                    Acceptable values for:   LOLIM      UPLIM          DRD  103
C                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48       DRD  104
C                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191      DRD  105
C                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201      DRD  106
C                    CRAY                 :   3.0D-1644   1.69D+1640    DRD  107
C                    VAX 11 SERIES        :   1.0D-25     4.5D+21       DRD  108
C                                                                       DRD  109
C                                                                       DRD  110
C          ERRTOL determines the accuracy of the answer                 DRD  111
C                                                                       DRD  112
C                 The value assigned by the routine will result         DRD  113
C                 in solution precision within 1-2 decimals of          DRD  114
C                 "machine precision".                                  DRD  115
C                                                                       DRD  116
C          ERRTOL    Relative error due to truncation is less than      DRD  117
C                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.               DRD  118
C                                                                       DRD  119
C                                                                       DRD  120
C                                                                       DRD  121
C        The accuracy of the computed approximation to the integral     DRD  122
C        can be controlled by choosing the value of ERRTOL.             DRD  123
C        Truncation of a Taylor series after terms of fifth order       DRD  124
C        introduces an error less than the amount shown in the          DRD  125
C        second column of the following table for each value of         DRD  126
C        ERRTOL in the first column.  In addition to the truncation     DRD  127
C        error there will be round-off error, but in practice the       DRD  128
C        total error from both sources is usually less than the         DRD  129
C        amount given in the table.                                     DRD  130
C                                                                       DRD  131
C                                                                       DRD  132
C                                                                       DRD  133
C                                                                       DRD  134
C          Sample choices:  ERRTOL   Relative truncation                DRD  135
C                                    error less than                    DRD  136
C                           1.0D-3    4.0D-18                           DRD  137
C                           3.0D-3    3.0D-15                           DRD  138
C                           1.0D-2    4.0D-12                           DRD  139
C                           3.0D-2    3.0D-9                            DRD  140
C                           1.0D-1    4.0D-6                            DRD  141
C                                                                       DRD  142
C                                                                       DRD  143
C                    Decreasing ERRTOL by a factor of 10 yields six moreDRD  144
C                    decimal digits of accuracy at the expense of one orDRD  145
C                    two more iterations of the duplication theorem.    DRD  146
C***LONG DESCRIPTION                                                    DRD  147
C                                                                       DRD  148
C   DRD Special Comments                                                DRD  149
C                                                                       DRD  150
C                                                                       DRD  151
C                                                                       DRD  152
C          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)                  DRD  153
C          = 3 / DSQRT(X * Y * Z), where X, Y, and Z are positive.      DRD  154
C                                                                       DRD  155
C                                                                       DRD  156
C          On Input:                                                    DRD  157
C                                                                       DRD  158
C          X, Y, and Z are the variables in the integral DRD(X,Y,Z).    DRD  159
C                                                                       DRD  160
C                                                                       DRD  161
C          On Output:                                                   DRD  162
C                                                                       DRD  163
C                                                                       DRD  164
C          X, Y, Z are unaltered.                                       DRD  165
C                                                                       DRD  166
C                                                                       DRD  167
C                                                                       DRD  168
C          ********************************************************     DRD  169
C                                                                       DRD  170
C          WARNING: Changes in the program may improve speed at the     DRD  171
C                   expense of robustness.                              DRD  172
C                                                                       DRD  173
C                                                                       DRD  174
C                                                                       DRD  175
C    -------------------------------------------------------------------DRD  176
C                                                                       DRD  177
C                                                                       DRD  178
C   Special double precision functions via DRD and DRF                  DRD  179
C                                                                       DRD  180
C                                                                       DRD  181
C                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind       DRD  182
C                                                                       DRD  183
C                  -----------------------------------------            DRD  184
C                                                                       DRD  185
C                                                                       DRD  186
C                                             2         2   2           DRD  187
C                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) - DRD  188
C                                                                       DRD  189
C                     2      3             2         2   2              DRD  190
C                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)      DRD  191
C                                                                       DRD  192
C                                                                       DRD  193
C                                  2        2            2              DRD  194
C                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)           DRD  195
C                                                                       DRD  196
C                         PI/2     2   2      1/2                       DRD  197
C                       = INT  (1-K SIN (PHI) )  D PHI                  DRD  198
C                          0                                            DRD  199
C                                                                       DRD  200
C                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind       DRD  201
C                                                                       DRD  202
C                  -----------------------------------------            DRD  203
C                                                                       DRD  204
C                                               2 2    2                DRD  205
C                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +             DRD  206
C                                                                       DRD  207
C                                              3          2 2    2      DRD  208
C                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )     DRD  209
C                                                                       DRD  210
C                                                                       DRD  211
C                                                                       DRD  212
C                                                                       DRD  213
C                  Legendre form of alternative ELLIPTIC INTEGRAL       DRD  214
C                  of 2nd kind                                          DRD  215
C                                                                       DRD  216
C                  -----------------------------------------            DRD  217
C                                                                       DRD  218
C                                                                       DRD  219
C                                                                       DRD  220
C                            Q     2       2   2  -1/2                  DRD  221
C                  D(Q,K) = INT SIN P  (1-K SIN P)     DP               DRD  222
C                            0                                          DRD  223
C                                                                       DRD  224
C                                                                       DRD  225
C                                                                       DRD  226
C                                     3          2     2   2            DRD  227
C                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)        DRD  228
C                                                                       DRD  229
C                                                                       DRD  230
C                                                                       DRD  231
C                                                                       DRD  232
C                  Lemniscate constant  B                               DRD  233
C                                                                       DRD  234
C                  -----------------------------------------            DRD  235
C                                                                       DRD  236
C                                                                       DRD  237
C                                                                       DRD  238
C                                                                       DRD  239
C                       1    2    4 -1/2                                DRD  240
C                  B = INT  S (1-S )    DS                              DRD  241
C                       0                                               DRD  242
C                                                                       DRD  243
C                                                                       DRD  244
C                  B = (1/3) DRD (0,2,1)                                DRD  245
C                                                                       DRD  246
C                                                                       DRD  247
C                  Heuman's LAMBDA function                             DRD  248
C                                                                       DRD  249
C                  -----------------------------------------            DRD  250
C                                                                       DRD  251
C                                                                       DRD  252
C                                                                       DRD  253
C                  (PI/2) LAMBDA0(A,B) =                                DRD  254
C                                                                       DRD  255
C                                    2                2                 DRD  256
C                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *            DRD  257
C                                                                       DRD  258
C                            2               2         2       2        DRD  259
C                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)  DRD  260
C                                                                       DRD  261
C                            2       3             2                    DRD  262
C                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *            DRD  263
C                                                                       DRD  264
C                           2         2       2                         DRD  265
C                   *DRD(COS (B),1-COS (A) SIN (B),1)                   DRD  266
C                                                                       DRD  267
C                                                                       DRD  268
C                                                                       DRD  269
C                  Jacobi ZETA function                                 DRD  270
C                                                                       DRD  271
C                  -----------------------------------------            DRD  272
C                                                                       DRD  273
C                             2                 2       2   2           DRD  274
C                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)     DRD  275
C                                                                       DRD  276
C                                                                       DRD  277
C                                       2             2                 DRD  278
C                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)              DRD  279
C                                                                       DRD  280
C                               2       3           2       2   2       DRD  281
C                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1) DRD  282
C                                                                       DRD  283
C                                                                       DRD  284
C --------------------------------------------------------------------- DRD  285
C          Subroutines or functions needed                              DRD  286
C              - XERROR                                                 DRD  287
C              - D1MACH                                                 DRD  288
C              - FORTRAN DABS, DMAX1,DMIN1, DSQRT                       DRD  289
C***REFERENCES  CARLSON, B.C. AND NOTIS,E .M.                           DRD  290
C                 ALGORITHMS FOR INCOMPLETE ELLIPTIC INTEGRALS          DRD  291
C                 ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,VOL.7,NO.3, DRD  292
C                 SEPT, 1981, PAGES 398-403                             DRD  293
C               CARLSON, B.C.                                           DRD  294
C                 COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION           DRD  295
C                 NUMER. MATH. 33, (1979), 1-16                         DRD  296
C               CARLSON, B.C.                                           DRD  297
C                 ELLIPTIC INTEGRALS OF THE FIRST KIND                  DRD  298
C                 SIAM J. MATH. ANAL. 8 (1977), 231-242                 DRD  299
C***ROUTINES CALLED  D1MACH,XERROR                                      DRD  300
C***END PROLOGUE  DRD                                                   DRD  301
      CHARACTER*176 MESSG                                               DRD  302
      INTEGER IER,ITODO                                                 DRD  303
      DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH             DRD  304
      DOUBLE PRECISION C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA        DRD  305
      DOUBLE PRECISION MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV          DRD  306
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV,      DRD  307
     * ZNROOT                                                           DRD  308
C                                                                       DRD  309
C                                                                       DRD  310
      SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,C4,ITODO                         DRD  311
C                                                                       DRD  312
C                                                                       DRD  313
      DATA ITODO/1/                                                     DRD  314
C                                                                       DRD  315
C                                                                       DRD  316
C                                                                       DRD  317
C***FIRST EXECUTABLE STATEMENT  DRD                                     DRD  318
      IF(ITODO.EQ.1)THEN                                                DRD  319
C                                                                       DRD  320
C                                                                       DRD  321
cryne      ERRTOL=(D1MACH(3)/3.0D0)**(1.0D0/6.0D0)                           DRD  322
c     errtol=1.d-3
      errtol=1.d-4
c     errtol=1.d-5
C                                                                       DRD  323
C                                                                       DRD  324
cryne      LOLIM = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)                          DRD  325
      lolim=1.d-25
C                                                                       DRD  326
cryne      UPLIM = D1MACH(1)**(1.0E0/3.0E0)                                  DRD  327
      uplim=1.d25
cryne      UPLIM = (0.10D0*ERRTOL)**(1.0E0/3.0E0)/UPLIM                      DRD  328
cryne      UPLIM = UPLIM**2.0D0                                              DRD  329
C                                                                       DRD  330
C                                                                       DRD  331
      C1 = 3.0D0/14.0D0                                                 DRD  332
      C2 = 1.0D0/6.0D0                                                  DRD  333
      C3 = 9.0D0/22.0D0                                                 DRD  334
      C4 = 3.0D0/26.0D0                                                 DRD  335
C                                                                       DRD  336
C                                                                       DRD  337
      ITODO=0                                                           DRD  338
C                                                                       DRD  339
      END IF                                                            DRD  340
C                                                                       DRD  341
C                                                                       DRD  342
C                                                                       DRD  343
C                                                                       DRD  344
C         CALL ERROR HANDLER IF NECESSARY.                              DRD  345
C                                                                       DRD  346
C                                                                       DRD  347
    5 DRD=0.0D0                                                         DRD  348
      IF( DMIN1(X,Y).LT.0.0D0) THEN                                     DRD  349
      IER=1                                                             DRD  350
      WRITE (MESSG,6) X, Y                                              DRD  351
    6 FORMAT('DRD - ERROR: DMIN1(X,Y).LT.0.0D0 WHERE X=', 1PD20.12,     DRD  352
     *   ' AND             Y=', D20.12)                                 DRD  353
cryne      CALL XERROR(MESSG(1:100),100,1,1)                                 DRD  354
      write(6,'(a100)')messg(1:100)
      RETURN                                                            DRD  355
      ENDIF                                                             DRD  356
      IF (DMAX1(X,Y,Z).GT.UPLIM) THEN                                   DRD  357
      IER=3                                                             DRD  358
      MESSG(1:43) = 'DRD - ERROR: DMAX1(X,Y,Z).GT.UPLIM WHERE X='       DRD  359
      WRITE (MESSG(44:176),7) X, Y, Z, UPLIM                            DRD  360
    7 FORMAT( 1PD20.12,                                                 DRD  361
     *  15X, 'Y=', D20.12, ' Z=', D20.12, ' AND', 23X, 'UPLIM=', D20.12)DRD  362
cryne      CALL XERROR(MESSG(1:176),176,3,1)                                 DRD  363
      write(6,'(a176)')messg(1:176)
      RETURN                                                            DRD  364
      ENDIF                                                             DRD  365
      IF (DMIN1(X+Y,Z).LT.LOLIM) THEN                                   DRD  366
      IER=2                                                             DRD  367
      MESSG(1:43)='DRD - ERROR: DMIN1(X+Y,Z).LT.LOLIM WHERE X='         DRD  368
      WRITE (MESSG(44:176),8) X, Y, Z, LOLIM                            DRD  369
    8 FORMAT( 1PD20.12,                                                 DRD  370
     *  15X, 'Y=', D20.12, ' Z=', D20.12, ' AND', 23X, 'LOLIM=', D20.12)DRD  371
cryne      CALL XERROR(MESSG(1:176),176,2,1)                                 DRD  372
      write(6,'(a176)')messg(1:176)
      RETURN                                                            DRD  373
      ENDIF                                                             DRD  374
C                                                                       DRD  375
C                                                                       DRD  376
C                                                                       DRD  377
C                                                                       DRD  378
   20 IER = 0                                                           DRD  379
      XN = X                                                            DRD  380
      YN = Y                                                            DRD  381
      ZN = Z                                                            DRD  382
      SIGMA = 0.0D0                                                     DRD  383
      POWER4 = 1.0D0                                                    DRD  384
C                                                                       DRD  385
C                                                                       DRD  386
C                                                                       DRD  387
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0                                      DRD  388
      XNDEV = (MU-XN)/MU                                                DRD  389
      YNDEV = (MU-YN)/MU                                                DRD  390
      ZNDEV = (MU-ZN)/MU                                                DRD  391
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))               DRD  392
      IF (EPSLON.LT.ERRTOL) GO TO 40                                    DRD  393
      XNROOT = DSQRT(XN)                                                DRD  394
      YNROOT = DSQRT(YN)                                                DRD  395
      ZNROOT = DSQRT(ZN)                                                DRD  396
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT                    DRD  397
      SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))                        DRD  398
      POWER4 = POWER4*0.250D0                                           DRD  399
      XN = (XN+LAMDA)*0.250D0                                           DRD  400
      YN = (YN+LAMDA)*0.250D0                                           DRD  401
      ZN = (ZN+LAMDA)*0.250D0                                           DRD  402
      GO TO 30                                                          DRD  403
C                                                                       DRD  404
C                                                                       DRD  405
C                                                                       DRD  406
C                                                                       DRD  407
   40 EA = XNDEV*YNDEV                                                  DRD  408
      EB = ZNDEV*ZNDEV                                                  DRD  409
      EC = EA - EB                                                      DRD  410
      ED = EA - 6.0D0*EB                                                DRD  411
      EF = ED + EC + EC                                                 DRD  412
      S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)                    DRD  413
      S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))                     DRD  414
      DRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*DSQRT(MU))           DRD  415
C                                                                       DRD  416
   50 RETURN                                                            DRD  417
C                                                                       DRD  418
C                                                                       DRD  419
C                                                                       DRD  420
C                                                                       DRD  421
      END                                                               DRD  422
