**********************************************************************
*  header:      COMBINED FUNCTION QUADRUPOLE ROUTINES (CFQD)         *
*    All routines for maps for combined function quadrupoles         *
**********************************************************************
c
      subroutine aarot(c, s, h, mh)
c
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision h(monoms),mh(6,6)
c
      call clear(h,mh)
c  Rotate coordinates
      mh(1,1)=c
      mh(1,3)=-s
      mh(3,1)=s
      mh(3,3)=c
c  Rotate momenta
      mh(2,2)=c
      mh(2,4)=-s
      mh(4,2)=s
      mh(4,4)=c
c  Don't touch flight time
      mh(5,5)=1.
      mh(6,6)=1.
c  Polynomials are zero (bless those linear maps).
      return
      end
c
**********************************************************************
c
      subroutine cfdrvr(pa,fa,fm)
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
cryne include 'parm.inc'
      include 'parset.inc'
c
      dimension fa(monoms)
      dimension fm(6,6)
      dimension pa(6),pb(6)
      dimension pc(6)
      dimension a(10), b(10)
c
      double precision gb0,k,k2,k3,k4,l,lk,lkm,lsc
      double precision sinv,cosv,ev,bedo
c
c  set up parameters and control indices
c
      al=pa(1)
      ilfrn=nint(pa(3))
      itfrn=nint(pa(4))
      ipset=nint(pa(2))
c
c  get multipole values from the parameter set ipset
c
      if (ipset.lt.1 .or. ipset.gt.maxpst) then
        do 50 i=1,6
   50   pb(i) = 0.0d0
      else
        do 60 i=1,6
   60   pb(i) = pst(i,ipset)
      endif
c
c compute multipole coefficients
c
      bsex=pb(3)
      asex=pb(4)
      boct=pb(5)
      aoct=pb(6)
      a(3) = asex
      a(4) = aoct
      b(3) = bsex
      b(4) = boct
      a(2) = pb(2)
      b(2) = pb(1)
c
      call gcfqd(al, b, a, fa, fm, ilfrn, itfrn)
c
      
      return
      end
c
**********************************************************************
c
      subroutine gcfqd(al, b, a, qa, qm, ilfr, itfr)
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
cryne include 'parm.inc'
c
      dimension qa(*), qm(6,6)
      dimension ha(monoms), hm(6,6)
c
      dimension a(*), b(*)
c
      absqd = b(2)**2 + a(2)**2
      if ( absqd .eq. 0.0d0 ) then
         asex = a(3)
         bsex = b(3)
         aoct = a(4)
         boct = b(4)
         call xcfdr(al, bsex, asex, boct, aoct, qa, qm)
      else
         if (a(2) .eq. 0.0d0 ) then
           if(b(2) .gt. 0 ) then
             cos1 = 1.d0
             cos2 = 1.d0
             cos3 = 1.d0
             cos4 = 1.d0
             sin1 = 0.d0
             sin2 = 0.d0
             sin3 = 0.d0
             sin4 = 0.d0
           else if ( b(2) .lt. 0.d0 ) then
             cos1 = 0.d0
             cos2 = -1.d0
             cos3 = 0.d0
             cos4 = 1.d0
             sin1 = 1.d0
             sin2 = 0.d0
             sin3 = -1.d0
             sin4 = 0.d0
           endif
         else
             phi2 = atan2(a(2), b(2))
             phi = phi2/2.d0
             phi3 = 3.d0*phi
             phi4 = 4.d0*phi
             cos1 = cos(phi)
             cos2 = cos(phi2)
             cos3 = cos(phi3)
             cos4 = cos(phi4)
             sin1 = sin(phi)
             sin2 = sin(phi2)
             sin3 = sin(phi3)
             sin4 = sin(phi4)
         endif
c
         bqd = sqrt(absqd)
c
         asex = a(3)*cos3 - b(3)*sin3
         aoct = a(4)*cos4 - b(4)*sin4
         bsex = b(3)*cos3 + a(3)*sin3
         boct = b(4)*cos4 + a(4)*sin4
c
         call aarot(cos1, sin1, qa, qm)
         if (ilfr .ne. 0 ) then
           call frquad(bqd,-1, ha, hm)
           call concat (qa, qm, ha, hm, qa, qm)
         endif
         call xcffqd(al, bqd, bsex, asex, boct, aoct, ha, hm)
         call concat (qa, qm, ha, hm, qa, qm)
         if (itfr .ne. 0 ) then
           call frquad(bqd, 1, ha, hm)
           call concat(qa, qm, ha, hm, qa, qm)
         endif
         call aarot(cos1, -sin1, ha, hm)
         call concat(qa, qm, ha, hm, qa, qm)
c
      endif
c
      return
      end
c
**********************************************************************
c
      subroutine xcfdr(al, bsex, asex, boct, aoct, fa, fm)
c
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
c
c  calling arrays
      dimension fa(monoms)
      dimension fm(6,6)
      dimension pa(6)
c
c  local arrays
      dimension ha(monoms)
      dimension hm(6,6)
      dimension pb(6)
c
      pb(1) = al
      pb(2) = bsex*al
      pb(3) = asex*al
      pb(4) = boct*al
      pb(5) = aoct*al
c
c  put map for cplm in fa,fm
      call cplm(pb,fa,fm)
c
c  put map for a drift of length al/2 in ha,hm
      alh=al/2.d0
      call drift(alh,ha,hm)
c
c  preceed and follow map for cplm with map of half-length
c  drift
      call concat(ha,hm,fa,fm,fa,fm)
      call concat(fa,fm,ha,hm,fa,fm)
c
      return
      end
c
**********************************************************************
c
      subroutine xcffqd(al, bqd, bsex, asex, boct, aoct, fa, fm)
c
c  This subroutine computes the map for a horizontally focusing
c  combined function quadrupole.
c  It uses analytic formulas produced by Johannes van Zeijts 
c  using the symbolic manipulation program Mathematica running
c  on a MacII-ci.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
cryne include 'parm.inc'
      include 'parset.inc'
c
      dimension fa(monoms)
      dimension fm(6,6)
      dimension pa(6),pb(6)
      dimension pc(6)
      dimension ha(monoms)
      dimension hm(6,6)
c
      double precision gb0,k,k2,k3,k4,l,lk,lkm,lsc
      double precision B3,B4,A3,A4,lambda
c
c     compute useful numbers
c
      sl2=sl*sl
      sl3=sl*sl2
      bet2=beta*beta
      bet3=beta*bet2
      gam2=gamma*gamma
c
      fqdnr=bqd/brho
      B3=bsex/brho
      A3=asex/brho
      B4=boct/brho
      A4=aoct/brho
c
      lsc=al/sl
	  lambda = lsc
c
c     evaluate k,lk
c
      arg=(bqd*(sl**2))/brho
      k=dsqrt(arg)
      lk=lsc*k
      lkm=(-1.0d0)*lk
c
c     set coefficients of fm
c
      chlk=(dexp(lk)+dexp(lkm))/(2.0d0)
      shlk=(dexp(lk)-dexp(lkm))/(2.0d0)
      clk=dcos(lk)
      slk=dsin(lk)
c
      call clear(fa,fm)
c
      fm(3,3)=+chlk
      fm(3,4)=+(shlk/k)
      fm(4,4)=+chlk
      fm(4,3)=+(k*shlk)
      fm(1,1)=+clk
      fm(1,2)=+(slk/k)
      fm(2,2)=+clk
      fm(2,1)=-(k*slk)
      fm(5,5)=+1.0d0
      fm(5,6)=+(lsc/((gamma**2)*(beta**2)))
      fm(6,6)=+1.0d0
c
      v = k*lsc
c  F3 and F4 terms ( computed by Johannes et al)
c
	  v = k*lambda
	  fa(28)=-(B3*(9*Sin(v) + Sin(3*v)))/(36*k)
      fa(29)=B3*(-4 + 3*Cos(v) + Cos(3*v))/(12*k**2)
      fa(30)=A3*(2*Cosh(v)*Sin(2*v) + 5*Sinh(v) + 
     # Cos(2*v)*Sinh(v))/(10*k)
      fa(31)=A3*(-6 + 5*Cosh(v) + Cos(2*v)*Cosh(v) + 
     # 2*Sin(2*v)*Sinh(v))/(10*k**2)
      fa(33)=k*(-2*v + Sin(2*v))/(8*beta)
      fa(34)=B3*(-3*Sin(v) + Sin(3*v))/(12*k**3)
      fa(35)=A3*(2 - 2*Cos(2*v)*Cosh(v) + Sin(2*v)*Sinh(v))/
     # (5*k**2)
      fa(36)=A3*(Cosh(v)*Sin(2*v) - 2*Cos(2*v)*Sinh(v))/
     # (5*k**3)
      fa(38)=(1 - Cos(2*v))/(4*beta)
      fa(39)=B3*(5*Sin(v) + Cosh(2*v)*Sin(v) + 
     # 2*Cos(v)*Sinh(2*v))/(10*k)
      fa(40)=B3*(-2 + 2*Cos(v)*Cosh(2*v) + Sin(v)*Sinh(2*v))/
     # (5*k**2)
      fa(43)=B3*(-5*Sin(v) + Cosh(2*v)*Sin(v) + 
     # 2*Cos(v)*Sinh(2*v))/(10*k**3)
      fa(49)=B3*(-8 + 9*Cos(v) - Cos(3*v))/(36*k**4)
      fa(50)=A3*(-2*Cosh(v)*Sin(2*v) + 5*Sinh(v) - 
     # Cos(2*v)*Sinh(v))/(10*k**3)
      fa(51)=A3*(-4 + 5*Cosh(v) - Cos(2*v)*Cosh(v) - 
     # 2*Sin(2*v)*Sinh(v))/(10*k**4)
      fa(53)=-(2*v + Sin(2*v))/(8*beta*k)
      fa(54)=B3*(6 - 5*Cos(v) - Cos(v)*Cosh(2*v) + 
     # 2*Sin(v)*Sinh(2*v))/(10*k**2)
      fa(55)=B3*(2*Cosh(2*v)*Sin(v) - Cos(v)*Sinh(2*v))/
     # (5*k**3)
      fa(58)=B3*(-4 + 5*Cos(v) - Cos(v)*Cosh(2*v) + 
     # 2*Sin(v)*Sinh(2*v))/(10*k**4)
      fa(64)=-(A3*(9*Sinh(v) + Sinh(3*v)))/(36*k)
      fa(65)=A3*(4 - 3*Cosh(v) - Cosh(3*v))/(12*k**2)
      fa(67)=k*(2*v - Sinh(2*v))/(8*beta)
      fa(68)=A3*(3*Sinh(v) - Sinh(3*v))/(12*k**3)
      fa(70)=(1 - Cosh(2*v))/(4*beta)
      fa(74)=A3*(-8 + 9*Cosh(v) - Cosh(3*v))/(36*k**4)
      fa(76)=-(2*v + Sinh(2*v))/(8*beta*k)
      fa(83)=-lambda/(2*beta**3*gamma**2)
      fa(84)=(-2700*k**7*lambda - 7920*v*A3**2 + 6000*v*B3**2 - 
     # 5400*k**3*lambda*B4 - 7200*B3**2*Sin(v) +1800*k**6*Sin(2*v)- 
     # 4320*A3**2*Sin(2*v) + 2400*B3**2*Sin(2*v) - 
     # 3600*k**2*B4*Sin(2*v) + 3456*A3**2*Cosh(v)*Sin(2*v) - 
     # 800*B3**2*Sin(3*v) - 225*k**6*Sin(4*v) - 
     # 180*A3**2*Sin(4*v) - 300*B3**2*Sin(4*v) - 
     # 450*k**2*B4*Sin(4*v) + 8640*A3**2*Sinh(v) + 
     # 1728*A3**2*Cos(2*v)*Sinh(v))/(57600*k**3)
      fa(85)=(45*k**6 - 156*A3**2 + 60*B3**2 -150*k**2*B4 -
     #  60*k**6*Cos(2*v) + 144*A3**2*Cos(2*v) - 80*B3**2*Cos(2*v) + 
     # 120*k**2*B4*Cos(2*v) + 15*k**6*Cos(4*v) + 
     # 12*A3**2*Cos(4*v) + 20*B3**2*Cos(4*v) + 
     # 30*k**2*B4*Cos(4*v) + 96*A3**2*Cosh(v) - 
     # 96*A3**2*Cos(2*v)*Cosh(v) + 
     # 96*A3**2*Sin(2*v)*Sinh(v))/(960*k**4)
      fa(86)=(270*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v) - 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 36*A3*B3*Cosh(2*v)*Sin(v) + 
     # 64*A3*B3*Cosh(v)*Sin(2*v) + 10*A3*B3*Sin(3*v) + 
     # 45*k**2*A4*Cosh(v)*Sin(3*v) - 
     # 12*A3*B3*Cosh(v)*Sin(3*v) + 160*A3*B3*Sinh(v) + 
     # 225*k**2*A4*Cos(v)*Sinh(v) - 
     # 380*A3*B3*Cos(v)*Sinh(v) + 
     # 32*A3*B3*Cos(2*v)*Sinh(v) + 
     # 15*k**2*A4*Cos(3*v)*Sinh(v) - 
     # 4*A3*B3*Cos(3*v)*Sinh(v) + 
     # 72*A3*B3*Cos(v)*Sinh(2*v))/(600*k**3)
      fa(87)=(-240*k**2*A4 + 384*A3*B3 - 90*A3*B3*Cos(v) - 
     # 30*A3*B3*Cos(3*v) + 40*A3*B3*Cosh(v) + 
     # 225*k**2*A4*Cos(v)*Cosh(v) - 
     # 380*A3*B3*Cos(v)*Cosh(v) + 
     # 8*A3*B3*Cos(2*v)*Cosh(v) + 
     # 15*k**2*A4*Cos(3*v)*Cosh(v) - 
     # 4*A3*B3*Cos(3*v)*Cosh(v) + 
     # 72*A3*B3*Cos(v)*Cosh(2*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) + 
     # 16*A3*B3*Sin(2*v)*Sinh(v) + 
     # 45*k**2*A4*Sin(3*v)*Sinh(v) - 
     # 12*A3*B3*Sin(3*v)*Sinh(v) + 
     # 36*A3*B3*Sin(v)*Sinh(2*v))/(600*k**4)
      fa(89)=B3*(-12*v - 9*v*Cos(v) - 3*v*Cos(3*v) - 
     # 3*Sin(v) + 6*Sin(2*v) + 5*Sin(3*v))/(144*beta*k)
      fa(90)=(-2700*k**7*lambda - 7920*v*A3**2 + 6000*v*B3**2 - 
     # 5400*k**3*lambda*B4 - 3600*B3**2*Sin(v) + 
     # 576*A3**2*Cosh(v)*Sin(2*v) - 2000*B3**2*Sin(3*v) + 
     # 675*k**6*Sin(4*v) + 540*A3**2*Sin(4*v) + 
     # 900*B3**2*Sin(4*v) + 1350*k**2*B4*Sin(4*v) + 
     # 7200*A3**2*Sinh(v) - 2592*A3**2*Cos(2*v)*Sinh(v))/
     # (28800*k**5)
      fa(91)=(360*k**2*A4 - 416*A3*B3 - 210*A3*B3*Cos(v) - 
     # 10*A3*B3*Cos(3*v) + 180*A3*B3*Cosh(v) - 
     # 225*k**2*A4*Cos(v)*Cosh(v) + 
     # 380*A3*B3*Cos(v)*Cosh(v) + 
     # 28*A3*B3*Cos(2*v)*Cosh(v) - 
     # 135*k**2*A4*Cos(3*v)*Cosh(v) + 
     # 36*A3*B3*Cos(3*v)*Cosh(v) + 
     # 12*A3*B3*Cos(v)*Cosh(2*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) + 
     # 76*A3*B3*Sin(2*v)*Sinh(v) + 
     # 45*k**2*A4*Sin(3*v)*Sinh(v) - 
     # 12*A3*B3*Sin(3*v)*Sinh(v) + 
     # 96*A3*B3*Sin(v)*Sinh(2*v))/(600*k**4)
      fa(92)=(-120*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v) - 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 96*A3*B3*Cosh(2*v)*Sin(v) + 
     # 4*A3*B3*Cosh(v)*Sin(2*v) - 80*A3*B3*Sin(3*v) + 
     # 45*k**2*A4*Cosh(v)*Sin(3*v) - 
     # 12*A3*B3*Cosh(v)*Sin(3*v) + 120*A3*B3*Sinh(v) - 
     # 225*k**2*A4*Cos(v)*Sinh(v) + 
     # 380*A3*B3*Cos(v)*Sinh(v) + 
     # 112*A3*B3*Cos(2*v)*Sinh(v) - 
     # 135*k**2*A4*Cos(3*v)*Sinh(v) + 
     # 36*A3*B3*Cos(3*v)*Sinh(v) + 
     # 12*A3*B3*Cos(v)*Sinh(2*v))/(600*k**5)
      fa(94)=B3*(-2 + 4*Cos(v) + 2*Cos(2*v) - 4*Cos(3*v) - 
     # 3*v*Sin(v) - 3*v*Sin(3*v))/(48*beta*k**2)
      fa(95)=(300*k**7*lambda + 1680*v*A3**2 - 1680*v*B3**2 + 
     # 1800*k**3*lambda*B4 + 2440*B3**2*Sin(v) + 
     # 272*B3**2*Cosh(2*v)*Sin(v) - 150*k**6*Sin(2*v) + 
     # 600*A3**2*Sin(2*v) - 440*B3**2*Sin(2*v) + 
     # 900*k**2*B4*Sin(2*v) - 544*A3**2*Cosh(v)*Sin(2*v) + 
     # 75*k**6*Cosh(2*v)*Sin(2*v) + 
     # 340*A3**2*Cosh(2*v)*Sin(2*v) + 
     # 140*B3**2*Cosh(2*v)*Sin(2*v) + 
     # 450*k**2*B4*Cosh(2*v)*Sin(2*v) + 120*B3**2*Sin(3*v) - 
     # 2440*A3**2*Sinh(v) - 272*A3**2*Cos(2*v)*Sinh(v) - 
     # 150*k**6*Sinh(2*v) + 440*A3**2*Sinh(2*v) - 
     # 600*B3**2*Sinh(2*v) + 900*k**2*B4*Sinh(2*v) + 
     # 544*B3**2*Cos(v)*Sinh(2*v) + 
     # 75*k**6*Cos(2*v)*Sinh(2*v) - 
     # 140*A3**2*Cos(2*v)*Sinh(2*v) - 
     # 340*B3**2*Cos(2*v)*Sinh(2*v) + 
     # 450*k**2*B4*Cos(2*v)*Sinh(2*v) - 120*A3**2*Sinh(3*v))/
     # (4800*k**3)
      fa(96)=(15*k**6 - 60*A3**2 + 188*B3**2 - 270*k**2*B4 - 
     # 24*B3**2*Cos(v) - 8*B3**2*Cos(3*v) - 88*A3**2*Cosh(v) + 
     # 112*A3**2*Cos(2*v)*Cosh(v) - 30*k**6*Cosh(2*v) + 
     # 88*A3**2*Cosh(2*v) - 120*B3**2*Cosh(2*v) + 
     # 180*k**2*B4*Cosh(2*v) + 32*B3**2*Cos(v)*Cosh(2*v) + 
     # 15*k**6*Cos(2*v)*Cosh(2*v) - 
     # 28*A3**2*Cos(2*v)*Cosh(2*v) - 
     # 68*B3**2*Cos(2*v)*Cosh(2*v) + 
     # 90*k**2*B4*Cos(2*v)*Cosh(2*v) - 24*A3**2*Cosh(3*v) - 
     # 64*A3**2*Sin(2*v)*Sinh(v) + 
     # 16*B3**2*Sin(v)*Sinh(2*v) + 
     # 15*k**6*Sin(2*v)*Sinh(2*v) + 
     # 68*A3**2*Sin(2*v)*Sinh(2*v) + 
     # 28*B3**2*Sin(2*v)*Sinh(2*v) + 
     # 90*k**2*B4*Sin(2*v)*Sinh(2*v))/(480*k**4)
      fa(98)=A3*(10*v + 5*v*Cosh(v) + 
     # 5*v*Cos(2*v)*Cosh(v) - 2*Sin(2*v) - 
     # 2*Cosh(v)*Sin(2*v) - Sinh(v) - 
     # 5*Cos(2*v)*Sinh(v) - 3*Sinh(2*v))/(40*beta*k)
      fa(99)=(-300*k**7*lambda - 1680*v*A3**2 + 1680*v*B3**2 - 
     # 1800*k**3*lambda*B4 - 160*B3**2*Sin(v) - 
     # 112*B3**2*Cosh(2*v)*Sin(v) + 150*k**6*Sin(2*v) - 
     # 600*A3**2*Sin(2*v) + 440*B3**2*Sin(2*v) - 
     # 900*k**2*B4*Sin(2*v) - 256*A3**2*Cosh(v)*Sin(2*v) + 
     # 75*k**6*Cosh(2*v)*Sin(2*v) + 
     # 340*A3**2*Cosh(2*v)*Sin(2*v) + 
     # 140*B3**2*Cosh(2*v)*Sin(2*v) + 
     # 450*k**2*B4*Cosh(2*v)*Sin(2*v) - 80*B3**2*Sin(3*v) + 
     # 1160*A3**2*Sinh(v) + 1312*A3**2*Cos(2*v)*Sinh(v) - 
     # 150*k**6*Sinh(2*v) + 440*A3**2*Sinh(2*v) - 
     # 600*B3**2*Sinh(2*v) + 900*k**2*B4*Sinh(2*v) - 
     # 224*B3**2*Cos(v)*Sinh(2*v) + 
     # 75*k**6*Cos(2*v)*Sinh(2*v) - 
     # 140*A3**2*Cos(2*v)*Sinh(2*v) - 
     # 340*B3**2*Cos(2*v)*Sinh(2*v) + 
     # 450*k**2*B4*Cos(2*v)*Sinh(2*v) - 120*A3**2*Sinh(3*v))/
     # (4800*k**5)
      fa(101)=A3*(-3 + 6*Cos(2*v) + 4*Cosh(v) - 
     # 4*Cos(2*v)*Cosh(v) - 3*Cosh(2*v) + 
     # 5*v*Sinh(v) + 5*v*Cos(2*v)*Sinh(v))/
     # (40*beta*k**2)
      fa(104)=k*(-8*v + 4*beta**2*v + 2*v*Cos(2*v) + 
     # 3*Sin(2*v) - 2*beta**2*Sin(2*v))/(32*beta**2)
      fa(105)=(75*k**6 - 132*A3**2 + 100*B3**2 - 90*k**2*B4 - 
     # 80*B3**2*Cos(v) - 60*k**6*Cos(2*v) + 144*A3**2*Cos(2*v) - 
     # 80*B3**2*Cos(2*v) + 120*k**2*B4*Cos(2*v) + 
     # 80*B3**2*Cos(3*v) - 15*k**6*Cos(4*v) - 
     # 12*A3**2*Cos(4*v) - 20*B3**2*Cos(4*v) - 
     # 30*k**2*B4*Cos(4*v) + 96*A3**2*Cosh(v) - 
     # 96*A3**2*Cos(2*v)*Cosh(v))/(960*k**6)
      fa(106)=(90*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v) - 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 72*A3*B3*Cosh(2*v)*Sin(v) + 
     # 208*A3*B3*Cosh(v)*Sin(2*v) + 10*A3*B3*Sin(3*v) - 
     # 135*k**2*A4*Cosh(v)*Sin(3*v) + 
     # 36*A3*B3*Cosh(v)*Sin(3*v) + 60*A3*B3*Sinh(v) + 
     # 225*k**2*A4*Cos(v)*Sinh(v) - 
     # 380*A3*B3*Cos(v)*Sinh(v) - 
     # 76*A3*B3*Cos(2*v)*Sinh(v) - 
     # 45*k**2*A4*Cos(3*v)*Sinh(v) + 
     # 12*A3*B3*Cos(3*v)*Sinh(v) + 
     # 24*A3*B3*Cos(v)*Sinh(2*v))/(600*k**5)
      fa(107)=(-180*k**2*A4 + 368*A3*B3 - 30*A3*B3*Cos(v)+ 
     # 70*A3*B3*Cos(3*v) - 60*A3*B3*Cosh(v) + 
     # 225*k**2*A4*Cos(v)*Cosh(v) - 
     # 380*A3*B3*Cos(v)*Cosh(v) - 
     # 4*A3*B3*Cos(2*v)*Cosh(v) - 
     # 45*k**2*A4*Cos(3*v)*Cosh(v) + 
     # 12*A3*B3*Cos(3*v)*Cosh(v) + 
     # 24*A3*B3*Cos(v)*Cosh(2*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) + 
     # 232*A3*B3*Sin(2*v)*Sinh(v) - 
     # 135*k**2*A4*Sin(3*v)*Sinh(v) + 
     # 36*A3*B3*Sin(3*v)*Sinh(v) + 
     # 72*A3*B3*Sin(v)*Sinh(2*v))/(600*k**6)
      fa(109)=B3*(-3*v*Cos(v) + 3*v*Cos(3*v) - 
     # 7*Sin(v) + 8*Sin(2*v) - 3*Sin(3*v))/(48*beta*k**3)
      fa(110)=(-15*k**6 + 188*A3**2 - 60*B3**2 + 270*k**2*B4- 
     # 88*B3**2*Cos(v) + 30*k**6*Cos(2*v) - 120*A3**2*Cos(2*v) + 
     # 88*B3**2*Cos(2*v) - 180*k**2*B4*Cos(2*v) - 
     # 24*B3**2*Cos(3*v) - 24*A3**2*Cosh(v) + 
     # 32*A3**2*Cos(2*v)*Cosh(v) + 
     # 112*B3**2*Cos(v)*Cosh(2*v) - 
     # 15*k**6*Cos(2*v)*Cosh(2*v) - 
     # 68*A3**2*Cos(2*v)*Cosh(2*v) - 
     # 28*B3**2*Cos(2*v)*Cosh(2*v) - 
     # 90*k**2*B4*Cos(2*v)*Cosh(2*v) - 8*A3**2*Cosh(3*v) - 
     # 16*A3**2*Sin(2*v)*Sinh(v) + 
     # 64*B3**2*Sin(v)*Sinh(2*v) + 
     # 15*k**6*Sin(2*v)*Sinh(2*v) - 
     # 28*A3**2*Sin(2*v)*Sinh(2*v) - 
     # 68*B3**2*Sin(2*v)*Sinh(2*v) + 
     # 90*k**2*B4*Sin(2*v)*Sinh(2*v))/(480*k**4)
      fa(111)=(-24*B3**2*Sin(v) + 16*B3**2*Cosh(2*v)*Sin(v) + 
     # 112*A3**2*Cosh(v)*Sin(2*v) + 
     # 15*k**6*Cosh(2*v)*Sin(2*v) - 
     # 28*A3**2*Cosh(2*v)*Sin(2*v) - 
     # 68*B3**2*Cosh(2*v)*Sin(2*v) + 
     # 90*k**2*B4*Cosh(2*v)*Sin(2*v) - 8*B3**2*Sin(3*v) - 
     # 24*A3**2*Sinh(v) + 16*A3**2*Cos(2*v)*Sinh(v) + 
     # 112*B3**2*Cos(v)*Sinh(2*v) - 
     # 15*k**6*Cos(2*v)*Sinh(2*v) - 
     # 68*A3**2*Cos(2*v)*Sinh(2*v) - 
     # 28*B3**2*Cos(2*v)*Sinh(2*v) - 
     # 90*k**2*B4*Cos(2*v)*Sinh(2*v) - 8*A3**2*Sinh(3*v))/
     # (240*k**5)
      fa(113)=A3*(1 - Cosh(2*v) + 5*v*Cosh(v)*Sin(2*v) - 
     # 4*Sin(2*v)*Sinh(v))/(20*beta*k**2)
      fa(114)=(45*k**6 - 52*A3**2 + 116*B3**2 - 90*k**2*B4 - 
     # 128*B3**2*Cos(v) - 30*k**6*Cos(2*v) + 
     # 120*A3**2*Cos(2*v) - 88*B3**2*Cos(2*v) + 
     # 180*k**2*B4*Cos(2*v) + 16*B3**2*Cos(3*v) - 
     # 24*A3**2*Cosh(v) + 32*A3**2*Cos(2*v)*Cosh(v) + 
     # 112*B3**2*Cos(v)*Cosh(2*v) - 
     # 15*k**6*Cos(2*v)*Cosh(2*v) - 
     # 68*A3**2*Cos(2*v)*Cosh(2*v) - 
     # 28*B3**2*Cos(2*v)*Cosh(2*v) - 
     # 90*k**2*B4*Cos(2*v)*Cosh(2*v) - 8*A3**2*Cosh(3*v) + 
     # 224*A3**2*Sin(2*v)*Sinh(v) - 
     # 32*B3**2*Sin(v)*Sinh(2*v) + 
     # 15*k**6*Sin(2*v)*Sinh(2*v) - 
     # 28*A3**2*Sin(2*v)*Sinh(2*v) - 
     # 68*B3**2*Sin(2*v)*Sinh(2*v) + 
     # 90*k**2*B4*Sin(2*v)*Sinh(2*v))/(480*k**6)
      fa(116)=A3*(5*Sin(2*v) - 3*Cosh(v)*Sin(2*v) - 
     # 2*Cos(2*v)*Sinh(v) + 5*v*Sin(2*v)*Sinh(v) - 
     # Sinh(2*v))/(20*beta*k**3)
      fa(119)=(2 - beta**2 - 2*Cos(2*v) + beta**2*Cos(2*v) + 
     # v*Sin(2*v))/(8*beta**2)
      fa(120)=(-160*A3*B3*Sin(v) - 225*k**2*A4*Cosh(v)*Sin(v)+ 
     # 380*A3*B3*Cosh(v)*Sin(v) - 
     # 32*A3*B3*Cosh(2*v)*Sin(v) - 
     # 15*k**2*A4*Cosh(3*v)*Sin(v) + 
     # 4*A3*B3*Cosh(3*v)*Sin(v) - 
     # 72*A3*B3*Cosh(v)*Sin(2*v) - 270*A3*B3*Sinh(v) - 
     # 225*k**2*A4*Cos(v)*Sinh(v) + 
     # 380*A3*B3*Cos(v)*Sinh(v) - 
     # 36*A3*B3*Cos(2*v)*Sinh(v) - 
     # 64*A3*B3*Cos(v)*Sinh(2*v) - 10*A3*B3*Sinh(3*v) - 
     # 45*k**2*A4*Cos(v)*Sinh(3*v) + 
     # 12*A3*B3*Cos(v)*Sinh(3*v))/(600*k**3)
      fa(121)=(360*k**2*A4 - 416*A3*B3 + 180*A3*B3*Cos(v)- 
     # 210*A3*B3*Cosh(v) - 225*k**2*A4*Cos(v)*Cosh(v) + 
     # 380*A3*B3*Cos(v)*Cosh(v) + 
     # 12*A3*B3*Cos(2*v)*Cosh(v) + 
     # 28*A3*B3*Cos(v)*Cosh(2*v) - 10*A3*B3*Cosh(3*v) - 
     # 135*k**2*A4*Cos(v)*Cosh(3*v) + 
     # 36*A3*B3*Cos(v)*Cosh(3*v) - 
     # 225*k**2*A4*Sin(v)*Sinh(v) + 
     # 380*A3*B3*Sin(v)*Sinh(v) - 
     # 96*A3*B3*Sin(2*v)*Sinh(v) - 
     # 76*A3*B3*Sin(v)*Sinh(2*v) - 
     # 45*k**2*A4*Sin(v)*Sinh(3*v) + 
     # 12*A3*B3*Sin(v)*Sinh(3*v))/(600*k**4)
      fa(123)=B3*(10*v + 5*v*Cos(v) + 
     # 5*v*Cos(v)*Cosh(2*v) - Sin(v) - 
     # 5*Cosh(2*v)*Sin(v) - 3*Sin(2*v) - 2*Sinh(2*v) - 
     # 2*Cos(v)*Sinh(2*v))/(40*beta*k)
      fa(124)=(60*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v) - 
     # 380*A3*B3*Cosh(v)*Sin(v) - 
     # 76*A3*B3*Cosh(2*v)*Sin(v) - 
     # 45*k**2*A4*Cosh(3*v)*Sin(v) + 
     # 12*A3*B3*Cosh(3*v)*Sin(v) + 
     # 24*A3*B3*Cosh(v)*Sin(2*v) + 90*A3*B3*Sinh(v) + 
     # 225*k**2*A4*Cos(v)*Sinh(v) - 
     # 380*A3*B3*Cos(v)*Sinh(v) + 
     # 72*A3*B3*Cos(2*v)*Sinh(v) + 
     # 208*A3*B3*Cos(v)*Sinh(2*v) + 10*A3*B3*Sinh(3*v) - 
     # 135*k**2*A4*Cos(v)*Sinh(3*v) + 
     # 36*A3*B3*Cos(v)*Sinh(3*v))/(600*k**5)
      fa(126)=B3*(-1 + Cos(2*v) + 5*v*Cos(v)*Sinh(2*v) - 
     # 4*Sin(v)*Sinh(2*v))/(20*beta*k**2)
      fa(130)=(-180*k**2*A4 + 368*A3*B3 - 180*A3*B3*Cos(v)+ 
     # 30*A3*B3*Cosh(v) + 225*k**2*A4*Cos(v)*Cosh(v) - 
     # 380*A3*B3*Cos(v)*Cosh(v) + 
     # 24*A3*B3*Cos(2*v)*Cosh(v) + 
     # 116*A3*B3*Cos(v)*Cosh(2*v) + 10*A3*B3*Cosh(3*v) - 
     # 45*k**2*A4*Cos(v)*Cosh(3*v) + 
     # 12*A3*B3*Cos(v)*Cosh(3*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) + 
     # 48*A3*B3*Sin(2*v)*Sinh(v) - 
     # 32*A3*B3*Sin(v)*Sinh(2*v) - 
     # 15*k**2*A4*Sin(v)*Sinh(3*v) + 
     # 4*A3*B3*Sin(v)*Sinh(3*v))/(600*k**6)
      fa(132)=B3*(-5*v*Cos(v) + 5*v*Cos(v)*Cosh(2*v) - 
     # 9*Sin(v) - 3*Cosh(2*v)*Sin(v) + 2*Sin(2*v) + 
     # 2*Sinh(2*v) + 2*Cos(v)*Sinh(2*v))/(40*beta*k**3)
      fa(140)=(-2700*k**7*lambda - 7920*v*A3**2 + 6000*v*B3**2 - 
     # 5400*k**3*lambda*B4 - 4800*B3**2*Sin(v) -1800*k**6*Sin(2*v)+ 
     # 4320*A3**2*Sin(2*v) - 2400*B3**2*Sin(2*v) + 
     # 3600*k**2*B4*Sin(2*v) - 2304*A3**2*Cosh(v)*Sin(2*v) + 
     # 1600*B3**2*Sin(3*v) - 225*k**6*Sin(4*v) - 
     # 180*A3**2*Sin(4*v) - 300*B3**2*Sin(4*v) - 
     # 450*k**2*B4*Sin(4*v) + 5760*A3**2*Sinh(v) - 
     # 1152*A3**2*Cos(2*v)*Sinh(v))/(57600*k**7)
      fa(141)=(180*k**2*A4 - 368*A3*B3 - 30*A3*B3*Cos(v) - 
     # 10*A3*B3*Cos(3*v) + 180*A3*B3*Cosh(v) - 
     # 225*k**2*A4*Cos(v)*Cosh(v) + 
     # 380*A3*B3*Cos(v)*Cosh(v) - 
     # 116*A3*B3*Cos(2*v)*Cosh(v) + 
     # 45*k**2*A4*Cos(3*v)*Cosh(v) - 
     # 12*A3*B3*Cos(3*v)*Cosh(v) - 
     # 24*A3*B3*Cos(v)*Cosh(2*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) - 
     # 32*A3*B3*Sin(2*v)*Sinh(v) - 
     # 15*k**2*A4*Sin(3*v)*Sinh(v) + 
     # 4*A3*B3*Sin(3*v)*Sinh(v) + 
     # 48*A3*B3*Sin(v)*Sinh(2*v))/(600*k**6)
      fa(142)=(-60*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v) - 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 48*A3*B3*Cosh(2*v)*Sin(v) - 
     # 8*A3*B3*Cosh(v)*Sin(2*v) + 20*A3*B3*Sin(3*v) - 
     # 15*k**2*A4*Cosh(v)*Sin(3*v) + 
     # 4*A3*B3*Cosh(v)*Sin(3*v) + 120*A3*B3*Sinh(v) - 
     # 225*k**2*A4*Cos(v)*Sinh(v) + 
     # 380*A3*B3*Cos(v)*Sinh(v) - 
     # 104*A3*B3*Cos(2*v)*Sinh(v) + 
     # 45*k**2*A4*Cos(3*v)*Sinh(v) - 
     # 12*A3*B3*Cos(3*v)*Sinh(v) - 
     # 24*A3*B3*Cos(v)*Sinh(2*v))/(600*k**7)
      fa(144)=B3*(-20 + 30*Cos(v) - 12*Cos(2*v) + 2*Cos(3*v) - 
     # 9*v*Sin(v) + 3*v*Sin(3*v))/(144*beta*k**4)
      fa(145)=(300*k**7*lambda + 1680*v*A3**2 - 1680*v*B3**2 + 
     # 1800*k**3*lambda*B4 + 1160*B3**2*Sin(v) + 
     # 1312*B3**2*Cosh(2*v)*Sin(v) + 150*k**6*Sin(2*v) - 
     # 600*A3**2*Sin(2*v) + 440*B3**2*Sin(2*v) - 
     # 900*k**2*B4*Sin(2*v) - 224*A3**2*Cosh(v)*Sin(2*v) - 
     # 75*k**6*Cosh(2*v)*Sin(2*v) - 
     # 340*A3**2*Cosh(2*v)*Sin(2*v) - 
     # 140*B3**2*Cosh(2*v)*Sin(2*v) - 
     # 450*k**2*B4*Cosh(2*v)*Sin(2*v) - 120*B3**2*Sin(3*v) - 
     # 160*A3**2*Sinh(v) - 112*A3**2*Cos(2*v)*Sinh(v) - 
     # 150*k**6*Sinh(2*v) + 440*A3**2*Sinh(2*v) - 
     # 600*B3**2*Sinh(2*v) + 900*k**2*B4*Sinh(2*v) - 
     # 256*B3**2*Cos(v)*Sinh(2*v) - 
     # 75*k**6*Cos(2*v)*Sinh(2*v) + 
     # 140*A3**2*Cos(2*v)*Sinh(2*v) + 
     # 340*B3**2*Cos(2*v)*Sinh(2*v) - 
     # 450*k**2*B4*Cos(2*v)*Sinh(2*v) - 80*A3**2*Sinh(3*v))/
     # (4800*k**5)
      fa(146)=(45*k**6 - 116*A3**2 + 52*B3**2 - 90*k**2*B4 + 
     # 24*B3**2*Cos(v) + 8*B3**2*Cos(3*v) + 
     # 128*A3**2*Cosh(v) - 112*A3**2*Cos(2*v)*Cosh(v) - 
     # 30*k**6*Cosh(2*v) + 88*A3**2*Cosh(2*v) - 
     # 120*B3**2*Cosh(2*v) + 180*k**2*B4*Cosh(2*v) - 
     # 32*B3**2*Cos(v)*Cosh(2*v) - 
     # 15*k**6*Cos(2*v)*Cosh(2*v) + 
     # 28*A3**2*Cos(2*v)*Cosh(2*v) + 
     # 68*B3**2*Cos(2*v)*Cosh(2*v) - 
     # 90*k**2*B4*Cos(2*v)*Cosh(2*v) - 16*A3**2*Cosh(3*v) - 
     # 32*A3**2*Sin(2*v)*Sinh(v) + 
     # 224*B3**2*Sin(v)*Sinh(2*v) - 
     # 15*k**6*Sin(2*v)*Sinh(2*v) - 
     # 68*A3**2*Sin(2*v)*Sinh(2*v) - 
     # 28*B3**2*Sin(2*v)*Sinh(2*v) - 
     # 90*k**2*B4*Sin(2*v)*Sinh(2*v))/(480*k**6)
      fa(148)=A3*(5*v*Cosh(v) - 5*v*Cos(2*v)*Cosh(v) - 
     # 2*Sin(2*v) - 2*Cosh(v)*Sin(2*v) + 9*Sinh(v) + 
     # 3*Cos(2*v)*Sinh(v) - 2*Sinh(2*v))/(40*beta*k**3)
      fa(149)=(-300*k**7*lambda - 1680*v*A3**2 + 1680*v*B3**2 - 
     # 1800*k**3*lambda*B4 - 1040*B3**2*Sin(v) + 
     # 928*B3**2*Cosh(2*v)*Sin(v) - 150*k**6*Sin(2*v) + 
     # 600*A3**2*Sin(2*v) - 440*B3**2*Sin(2*v) + 
     # 900*k**2*B4*Sin(2*v) + 64*A3**2*Cosh(v)*Sin(2*v) - 
     # 75*k**6*Cosh(2*v)*Sin(2*v) - 
     # 340*A3**2*Cosh(2*v)*Sin(2*v) - 
     # 140*B3**2*Cosh(2*v)*Sin(2*v) - 
     # 450*k**2*B4*Cosh(2*v)*Sin(2*v) + 80*B3**2*Sin(3*v) + 
     # 1040*A3**2*Sinh(v) - 928*A3**2*Cos(2*v)*Sinh(v) - 
     # 150*k**6*Sinh(2*v) + 440*A3**2*Sinh(2*v) - 
     # 600*B3**2*Sinh(2*v) + 900*k**2*B4*Sinh(2*v) - 
     # 64*B3**2*Cos(v)*Sinh(2*v) - 
     # 75*k**6*Cos(2*v)*Sinh(2*v) + 
     # 140*A3**2*Cos(2*v)*Sinh(2*v) + 
     # 340*B3**2*Cos(2*v)*Sinh(2*v) - 
     # 450*k**2*B4*Cos(2*v)*Sinh(2*v) - 80*A3**2*Sinh(3*v))/
     # (4800*k**7)
      fa(151)=A3*(-10 - 4*Cos(2*v) + 14*Cosh(v) + 
     # 2*Cos(2*v)*Cosh(v) - 2*Cosh(2*v) + 
     # 5*v*Sinh(v) - 5*v*Cos(2*v)*Sinh(v) - 
     # 4*Sin(2*v)*Sinh(v))/(40*beta*k**4)
      fa(154)=(-12*v + 4*beta**2*v - 2*v*Cos(2*v) - 
     # 5*Sin(2*v) + 2*beta**2*Sin(2*v))/(32*beta**2*k)
      fa(155)=(-240*k**2*A4 + 384*A3*B3 + 40*A3*B3*Cos(v) - 
     # 90*A3*B3*Cosh(v) + 225*k**2*A4*Cos(v)*Cosh(v) - 
     # 380*A3*B3*Cos(v)*Cosh(v) + 
     # 72*A3*B3*Cos(2*v)*Cosh(v) + 
     # 8*A3*B3*Cos(v)*Cosh(2*v) - 30*A3*B3*Cosh(3*v) + 
     # 15*k**2*A4*Cos(v)*Cosh(3*v) - 
     # 4*A3*B3*Cos(v)*Cosh(3*v) - 
     # 225*k**2*A4*Sin(v)*Sinh(v) + 
     # 380*A3*B3*Sin(v)*Sinh(v) - 
     # 36*A3*B3*Sin(2*v)*Sinh(v) - 
     # 16*A3*B3*Sin(v)*Sinh(2*v) - 
     # 45*k**2*A4*Sin(v)*Sinh(3*v) + 
     # 12*A3*B3*Sin(v)*Sinh(3*v))/(600*k**4)
      fa(156)=(120*A3*B3*Sin(v) - 225*k**2*A4*Cosh(v)*Sin(v) + 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 112*A3*B3*Cosh(2*v)*Sin(v) - 
     # 135*k**2*A4*Cosh(3*v)*Sin(v) + 
     # 36*A3*B3*Cosh(3*v)*Sin(v) + 
     # 12*A3*B3*Cosh(v)*Sin(2*v) - 120*A3*B3*Sinh(v) + 
     # 225*k**2*A4*Cos(v)*Sinh(v) - 
     # 380*A3*B3*Cos(v)*Sinh(v) + 
     # 96*A3*B3*Cos(2*v)*Sinh(v) + 
     # 4*A3*B3*Cos(v)*Sinh(2*v) - 80*A3*B3*Sinh(3*v) + 
     # 45*k**2*A4*Cos(v)*Sinh(3*v) - 
     # 12*A3*B3*Cos(v)*Sinh(3*v))/(600*k**5)
      fa(158)=B3*(3 - 4*Cos(v) + 3*Cos(2*v) - 6*Cosh(2*v) + 
     # 4*Cos(v)*Cosh(2*v) + 5*v*Sin(v) + 
     # 5*v*Cosh(2*v)*Sin(v))/(40*beta*k**2)
      fa(159)=(180*k**2*A4 - 368*A3*B3 + 60*A3*B3*Cos(v) + 
     # 30*A3*B3*Cosh(v) - 225*k**2*A4*Cos(v)*Cosh(v) + 
     # 380*A3*B3*Cos(v)*Cosh(v) - 
     # 24*A3*B3*Cos(2*v)*Cosh(v) + 
     # 4*A3*B3*Cos(v)*Cosh(2*v) - 70*A3*B3*Cosh(3*v) + 
     # 45*k**2*A4*Cos(v)*Cosh(3*v) - 
     # 12*A3*B3*Cos(v)*Cosh(3*v) + 
     # 225*k**2*A4*Sin(v)*Sinh(v) - 
     # 380*A3*B3*Sin(v)*Sinh(v) + 
     # 72*A3*B3*Sin(2*v)*Sinh(v) + 
     # 232*A3*B3*Sin(v)*Sinh(2*v) - 
     # 135*k**2*A4*Sin(v)*Sinh(3*v) + 
     # 36*A3*B3*Sin(v)*Sinh(3*v))/(600*k**6)
      fa(161)=B3*(2*Cosh(2*v)*Sin(v) + Sin(2*v) - 
     # 5*Sinh(2*v) + 3*Cos(v)*Sinh(2*v) + 
     # 5*v*Sin(v)*Sinh(2*v))/(20*beta*k**3)
      fa(165)=(-120*A3*B3*Sin(v) + 225*k**2*A4*Cosh(v)*Sin(v)- 
     # 380*A3*B3*Cosh(v)*Sin(v) + 
     # 104*A3*B3*Cosh(2*v)*Sin(v) - 
     # 45*k**2*A4*Cosh(3*v)*Sin(v) + 
     # 12*A3*B3*Cosh(3*v)*Sin(v) + 
     # 24*A3*B3*Cosh(v)*Sin(2*v) + 60*A3*B3*Sinh(v) - 
     # 225*k**2*A4*Cos(v)*Sinh(v) + 
     # 380*A3*B3*Cos(v)*Sinh(v) - 
     # 48*A3*B3*Cos(2*v)*Sinh(v) + 
     # 8*A3*B3*Cos(v)*Sinh(2*v) - 20*A3*B3*Sinh(3*v) + 
     # 15*k**2*A4*Cos(v)*Sinh(3*v) - 
     # 4*A3*B3*Cos(v)*Sinh(3*v))/(600*k**7)
      fa(167)=B3*(-10 + 14*Cos(v) - 2*Cos(2*v) - 4*Cosh(2*v) + 
     # 2*Cos(v)*Cosh(2*v) - 5*v*Sin(v) + 
     # 5*v*Cosh(2*v)*Sin(v) + 4*Sin(v)*Sinh(2*v))/
     # (40*beta*k**4)
      fa(175)=(-2700*k**7*lambda - 6000*v*A3**2 + 7920*v*B3**2 - 
     # 5400*k**3*lambda*B4 - 8640*B3**2*Sin(v) - 
     # 1728*B3**2*Cosh(2*v)*Sin(v) + 7200*A3**2*Sinh(v) + 
     # 1800*k**6*Sinh(2*v) - 2400*A3**2*Sinh(2*v) + 
     # 4320*B3**2*Sinh(2*v) - 3600*k**2*B4*Sinh(2*v) - 
     # 3456*B3**2*Cos(v)*Sinh(2*v) + 800*A3**2*Sinh(3*v) - 
     # 225*k**6*Sinh(4*v) + 300*A3**2*Sinh(4*v) + 
     # 180*B3**2*Sinh(4*v) - 450*k**2*B4*Sinh(4*v))/(57600*k**3)
      fa(176)=(-45*k**6 + 60*A3**2 - 156*B3**2 + 150*k**2*B4 + 
     # 96*B3**2*Cos(v) + 60*k**6*Cosh(2*v) - 
     # 80*A3**2*Cosh(2*v) + 144*B3**2*Cosh(2*v) - 
     # 120*k**2*B4*Cosh(2*v) - 96*B3**2*Cos(v)*Cosh(2*v) - 
     # 15*k**6*Cosh(4*v) + 20*A3**2*Cosh(4*v) + 
     # 12*B3**2*Cosh(4*v) - 30*k**2*B4*Cosh(4*v) - 
     # 96*B3**2*Sin(v)*Sinh(2*v))/(960*k**4)
      fa(178)=A3*(-12*v - 9*v*Cosh(v) - 3*v*Cosh(3*v) - 
     # 3*Sinh(v) + 6*Sinh(2*v) + 5*Sinh(3*v))/(144*beta*k)
      fa(179)=(2700*k**7*lambda + 6000*v*A3**2 - 7920*v*B3**2 + 
     # 5400*k**3*lambda*B4 + 7200*B3**2*Sin(v) - 
     # 2592*B3**2*Cosh(2*v)*Sin(v) - 3600*A3**2*Sinh(v) + 
     # 576*B3**2*Cos(v)*Sinh(2*v) - 2000*A3**2*Sinh(3*v) - 
     # 675*k**6*Sinh(4*v) + 900*A3**2*Sinh(4*v) + 
     # 540*B3**2*Sinh(4*v) - 1350*k**2*B4*Sinh(4*v))/(28800*k**5)
      fa(181)=A3*(2 - 4*Cosh(v) - 2*Cosh(2*v) + 4*Cosh(3*v) - 
     # 3*v*Sinh(v) - 3*v*Sinh(3*v))/(48*beta*k**2)
      fa(184)=k*(8*v - 4*beta**2*v - 2*v*Cosh(2*v) - 
     # 3*Sinh(2*v) + 2*beta**2*Sinh(2*v))/(32*beta**2)
      fa(185)=(75*k**6 - 100*A3**2 + 132*B3**2 - 90*k**2*B4 - 
     # 96*B3**2*Cos(v) + 80*A3**2*Cosh(v) - 60*k**6*Cosh(2*v) + 
     # 80*A3**2*Cosh(2*v) - 144*B3**2*Cosh(2*v) + 
     # 120*k**2*B4*Cosh(2*v) + 96*B3**2*Cos(v)*Cosh(2*v) - 
     # 80*A3**2*Cosh(3*v) - 15*k**6*Cosh(4*v) + 
     # 20*A3**2*Cosh(4*v) + 12*B3**2*Cosh(4*v) - 
     # 30*k**2*B4*Cosh(4*v))/(960*k**6)
      fa(187)=A3*(3*v*Cosh(v) - 3*v*Cosh(3*v) + 
     # 7*Sinh(v) - 8*Sinh(2*v) + 3*Sinh(3*v))/(48*beta*k**3)
      fa(190)=(2 - beta**2 - 2*Cosh(2*v) + beta**2*Cosh(2*v) - 
     # v*Sinh(2*v))/(8*beta**2)
      fa(195)=(-2700*k**7*lambda - 6000*v*A3**2 + 7920*v*B3**2 - 
     # 5400*k**3*lambda*B4 - 5760*B3**2*Sin(v) + 
     # 1152*B3**2*Cosh(2*v)*Sin(v) + 4800*A3**2*Sinh(v) - 
     # 1800*k**6*Sinh(2*v) + 2400*A3**2*Sinh(2*v) - 
     # 4320*B3**2*Sinh(2*v) + 3600*k**2*B4*Sinh(2*v) + 
     # 2304*B3**2*Cos(v)*Sinh(2*v) - 1600*A3**2*Sinh(3*v) - 
     # 225*k**6*Sinh(4*v) + 300*A3**2*Sinh(4*v) + 
     # 180*B3**2*Sinh(4*v) - 450*k**2*B4*Sinh(4*v))/(57600*k**7)
      fa(197)=A3*(-20 + 30*Cosh(v) - 12*Cosh(2*v) + 2*Cosh(3*v) + 
     # 9*v*Sinh(v) - 3*v*Sinh(3*v))/(144*beta*k**4)
      fa(200)=(-12*v + 4*beta**2*v - 2*v*Cosh(2*v) - 
     # 5*Sinh(2*v) + 2*beta**2*Sinh(2*v))/(32*beta**2*k)
      fa(209)=(-5 + beta**2)*lambda/(8*beta**4*gamma**2)
c
c  go from reverse to direct factorization:
c
      call revf(1,fa,fm)
c
c  fringe field effects are put on in the subroutine gcfqd
c
      return
      end
c
c end of file

