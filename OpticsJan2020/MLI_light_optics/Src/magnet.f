c  Minor changes made 6/04 by P. Walstrom to help it pass FTNCHEK
c    P Walstrom's current sheet magnet routines for all gtypes
c
      subroutine onaxgr(init,zfp,icoil,ndriv,dg)
c  modified 8-25 to use modified gtype6=gtipe6.
ctm   renamed back to gtype6 and gtype17 for MaryLie version
      implicit double precision(a-h,o-z)
c  Requires initialization calls for some values of icoil,according
c  to value of itype(icoil). See below.
c  Also, must call subroutines BINCOF,GAULEG and COEFFS once
c  before the above initialization calls
c   Applies to cylindrical current sheet coils or PMMs. Some types can h
c  an infinitely long, cylindrical infinite-mu shield surrrounding the
c  windings.
c
c  Written by P. L. Walstrom, Grumman Space Systems,6-90, LANL.
c  Modified by P. L. Walstrom 4/91 to include thick and thin Halbach PMM
c  Modified 8-96 to include coils with an iron cylinder or shield
c  around the windings. This is done by adding an equivalent cylindrcal
c  current sheet of radius b, the iron radius, that duplicates the
c  iron contribution for r<b.  Requires initialization call for each
c  such coil.
c
c  This routine calculates the on-axis gradient of a pure m-pole
c  surface winding.(Or PMM) That is, is calculates the lowest order
c   coefficient
c  F1(z) in the expansion of the scalar potential about r=0,
c
c   V(r,phi,z)=sin(m*phi)*(r**m)*(F1(z)+(r**2)*F2(z)+...)
c
c  andthe z derivatives of F1(z). F1(z) has units of Tesla*meters**(1-m
c  This definition of on-axis gradient differs by a factor of m from the
c  usual one.
c  If a surface-winding-type coil,the windings lie on a surface of radiu
c  areof infinitesimal radial thickness, and are continuous. Can also
c  have thick coils- this is done by Gaussian integration of thin
c  coils in radial depth.
c
c  Thesurface current densities are described by a stream function
c   psi(phi,z), where
c   psi(phi,z)=NIsin(m*phi)f(z)
c   NI=amp-turns/pole
c   Jz=mNIcos(m*phi)*f(z)/a
c   Jphi=-NIsin(m*phi)*df/dz
c  f(z) above is called the shape function.
c
c  Thegradient-length product of a coil can be zero or non-zero. If it
c  is non-zero, it can be used to determine NI. If not, NI must be defin
c  in some other way. For example, the higher harmonics of Lambertson
c  coils have shape functions that average to zero, and the correspondin
c  gradient-length products are also zero, for zero average shape fuctio
c  means zero integrated on-axis gradient.
c  Forsurface-current windings with non-zero average shape function ,
c
c  glprod=(mu0*NI/2)*Leff*a**(-m)*(1+(a/b)**2m)     Eq. (1)
c
c  where Leff=integral from -inf to +inf of f(z)dz
c  andb is the radius of a iron shield, if present.(None for this routi
c
c  Forcoils with zero average shape function, and user-specified shape
c  functions, THIS ROUTINE ASSIGNS VARIOUS MEANINGS TO GLPROD, as follow
c
c  Forhigher harmonics of Lambertsons (itype=4) NI is calculated with E
c  in the same way as for the fundamental, using Leff=(Lmax+Lmin)/2
c
c  Forthe user-specified shape function (itype=6), NI is calculated
c  with Leff=L (=alcoil(icoil)), irrespective of the real value of Leff,
c  using Eq. 1,i.e., glprod=(mu0*NI/2)*L*a**(-m)*(1+(a/b)**2m)
c
c  Coils with more than one m value and/or with both skew and non-skew
c  winding components are treated as separate coils from the point of
c  view of this subroutine. For example, a Lambertson m-pole is taken
c  to be two coincident coils-- one, a pure m-pole with a shape function
c  given by FZ from the subroutine FLAMB, and the second a pure 3m-pole
c  with a shape function given by FLAMB as F3Z.
c  Thehigher (2k-1)m shape functions, k=3,4,5... are set equal to zero
c  Lambertson coils.
c  Forgtype9 magnets (rings of z-oriented dipoles), the glprod is the
c  glprod that would arise if the dipoles were to be rotated 90 deg. aro
c  local tangent vector to the ring. See comments in subroutine gtype9 f
c  more details.
c
c   All units are MKS.
c
c  Input variables via subroutine call statement:
c   init=initialization code, =0 for initialization, 1 for gradients.
c   zfp=absolute field point z coordinate
c   icoil=index from 1 to maxcoil of coil.
c   ndriv=1 plus the number of times the on-axis gradient is to be
c   differentiated.
c
c  Input variables via common COILS:
c   mcoil(icoil)=multipole index
c   acoil(icoil)=coil winding radius
c   alcoil(icoil)=coil length
c   zcoil=z coordinate of center of coil
c   shape(icoil,i)-----has different meanings, depending upon itype(icoi
c   glprd(icoil)=For itype=1,2,3,5,7,8--=generalized gradient-length pro
c   For itype=4,6, and 9 ,see above comments.
c   itype(icoil)=integer specifying the type of coil winding:
c     itype=1 means quartic shape function in z.  shape(icoil,1)=eslope,
c     and the remaining parameters in shape(icoil,i) are not use
c     itype=2 means a thick Halbach permanent magnet multipole.
c     In this case, shape(icoil,1)=a2, the outer radius.
c     Also, acoil(icoil)=a1=inner radius.
c     itype=3 means a Lambertson coil fundamental shape function. In thi
c    case, shape(icoil,1)=xlmin,the length of the shortest straight segm
c     itype=4 means a Lambertson coil first harmonic shape function
c  In this case, mcoil has the value of 3m, where m is the mcoil value o
c  thefundamental.Shape(icoil,1) is again xlmin.
c     itype=5 means a truncated 2-D cosine coil-i.e. currents
c     are parallel to the z axis, except at the ends, where
c     they cross over in  filament loops.  shape is not used.
c     itype=6 means a user-defined shape function in the form of a data
c     file. The data file is called FORifile.dat, or defined to be
c     FORifile.dat by a command file during execution.ifile=shape(icoil
c     The no. of points in FORifile.dat is np=shape(icoil,2)
c  itype=7 means a shape function with a flat center and rounded ends.
c  Theends are given by parabolas tangent to the flattop.
c
c   itype=8 means a dipole sheet m-pole coil with a dipole density of
c  constant magnitude and varying direction, as in Halbach's permanent
c  magnet multipoles, ie.
c
c   M_r = -M_0 * sin( m*phi)
c
c  M_phi= M_0 * cos (m*phi)
c
c  Unlike Halbach's PMMs, these dipole sheets have no radial depth.
c  itype=2 magnets have the same angular dependence of the dipole moment
c  volume density  as the itype=8 magnets.
c
c  itype=9 means a ring with a z-oriented line dipole density varying as
c  sinmphi. This magnet type has zero integral strength-it is a sort of
c  all-fringe-field magnet that produces an on-axis gradient with a
c  positive and negative bump, with zero gradient at z=0. Since this
c  element has zero length, there is no integration over z in the gradie
c  evaluation.
c
c   itype=10 is a dipole sheet with sin mphi variation of a dipole densi
c   pointing in the z direction (like itype9, but finite length)
c
c   itype=3,4, and 6 require initialization calls- the others do not.
c
c  itype=11 is a current-sheet coil with a flat top and quartic+quadrati
c  ends. The ends can have different widths and end slopes.
c
c  itype=12 is a current-sheet coil with a flat top and cubic+quadratic
c  ends. The ends can have different widths and end slopes.
c
c  itype=14 is a current-sheet coil that has a shape function with a fla
c  topand rounded ends. The ends are formed by parabolae tangent to the
c  flat top. Thick analog of itype 7- uses Gaussian integration in radia
c  depth.
c
c  itype=15 is a current-sheet coil with a flat top and quartic+quadrati
c  ends. The ends can have different widths and end slopes.
c  Like gtype11, but has finite radial thickness- uses Gaussian integrat
c  of gtype11 gradients.
c
c  itype=16 is a current-sheet coil with a flat top and cubic+quadratic
c  ends. The ends can have different widths and end slopes.
c  Like gtype12, but has finite radial thickness- uses Gaussian integrat
c  of gtype12 gradients.
c
c  itype=17 is a coil with a flat top and quartic+quadratic ends and
c  an iron cylinder. The ends can have different widths and end slopes.
c  Theshape function is gtype11, but the coil has finite radial thickne
c  uses Gaussian integration in radial depth of gtype11 gradients.
c  Theiron contribution is computed from an equivalent surface winding
c  of radius b, the iron radius.  Requires an initialization call to
c  compute the iron shape function, which is represented by a piecewise
c  cubic.
c
c
c    Other itypes to be added later.
c
c  Output variables:
c     dg=vector with ndriv elements containing F1(z),dF1/dz,d2F1/dz2, e
c
c   The method of calculation of dg is to integrate a Green's
c   function g(z-zprime) and its z derivatives times f(zprime) from
c   zprime=-L/2 to +L/2.The integration is done by exact quadrature, if
c  theshape function can be represented by a quartic or lower piecewise
c  continuous polynomial. The user-defined shape function is treated
c  by a piecewise quadratic interpolation of the x and f(x) points suppl
c  TheLambertson coil fundamental and first harmonic are also treated
c  in this way. For current-sheet coils only (PMMs have similar expressi
c
c   g(z)=const.*(m/((z**2+a**2)**(m+1/2))-
c     (2m+1)*a**2/((z**2+a**2)**(m+3/2))
c    where
c  const.=mu0*NI*(a**m)*(1*3*5*...(2m-1))/(m!*2**(m+1))
c
c  Forshape functions that do not average to zero,
c
c  const.=glprod*(a**2m)*(1*3*5*...(2m-1))/(2**m*m!*xleff)
c
c
      parameter (maxcoils=100,maxshape=6)
      parameter(maxdrv=20)
      common/coils/ mcoil(maxcoils),acoil(maxcoils),alcoil(maxcoils),
     &zcoil(maxcoils),shape(maxcoils,maxshape),glprd(maxcoils),
     &itype(maxcoils)
      dimension dg(maxdrv)
c
c
c
      itipe=itype(icoil)
      a=acoil(icoil)
      xl=alcoil(icoil)
      z=zfp-zcoil(icoil)
      m=mcoil(icoil)
      glprod=glprd(icoil)
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),itipe
    1 eslope=shape(icoil,1)
c  Coil is type 1- symmetric quartic shape function of unit height.
c  dimensionless shape function end slope is shape(icoil,1)
      call gtype1(z,a,xl,eslope,glprod,m,ndriv,dg)
      return
c   Thick Halbach coil
    2 a1=a
      a2=shape(icoil,1)
      call gtype2(z,a1,a2,xl,glprod,m,ndriv,dg)
      return
c  Lambertson coil fundamental-k=0
    3 xlmin=shape(icoil,1)
      if(init.gt.0) go to 32
      call gtype3(0,0,icoil,xl,xlmin,0.d0,a,glprod,m,0,dg)
      return
   32 call gtype3(1,0,icoil,xl,xlmin,z,a,glprod,m,ndriv,dg)
      return
c   Lambertson coil first harmonic
c  Uses subroutine GTYPE3 with k=1 instead of k=0.
    4 xlmin=shape(icoil,1)
      if(init.gt.0) go to 42
      call gtype3(0,1,icoil,xl,xlmin,0.d0,a,glprod,m,0,dg)
      return
   42 call gtype3(1,1,icoil,xl,xlmin,z,a,glprod,m,ndriv,dg)
      return
c   Square shape function
    5 call gtype5(z,a,xl,glprod,m,ndriv,dg)
      return
c  User-defined shape function
    6 if(init.gt.0) go to 62
      ifile=shape(icoil,1)
c
c  Count points in GTIPE6- don't use SHAPE().
c     np=shape(icoil,2)
      call gtype6(0,icoil,xl,glprod,ifile,a,0.d0,m,0,dg)
      return
   62 call gtype6(1,icoil,xl,glprod,0,a,z,m,ndriv,dg)
      return
c  Flattop shape function with rounded ends.
    7 wflat=shape(icoil,1)
      call gtype7(z,a,xl,wflat,glprod,m,ndriv,dg)
      return
c  Dipole sheet coil
    8 call gtype8(z,a,xl,glprod,m,ndriv,dg)
      return
    9 call gtype9(z,a,glprod,m,ndriv,dg)
      return
   10 call gtype10(z,a,xl,glprod,m,ndriv,dg)
      return
c  Flattop with different quartic + quadratic ends
   11 w1=shape(icoil,1)
      w2=shape(icoil,2)
      s1=shape(icoil,3)
      s2=shape(icoil,4)
c     type *,alcoil(icoil)
c     type *,w1,w2
c     type *,s1,s2
      call gtype11(z,a,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      return
c  Flattop with different cubic +quadratic ends
   12 w1=shape(icoil,1)
      w2=shape(icoil,2)
      s1=shape(icoil,3)
      s2=shape(icoil,4)
c     type *,alcoil(icoil)
c     type *,w1,w2
c     type *,s1,s2
      call gtype12(z,a,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      return
c  Thick, symmetric quartic+quadratic with no flattop.
   13 s=shape(icoil,1)
      a2=shape(icoil,2)
      call gtype13(z,a,a2,xl,s,glprod,m,ndriv,dg)
      return
c  Thick symmetric flattop with parabolic ends.
   14 wflat=shape(icoil,1)
      a2=shape(icoil,2)
      call gtype14(z,a,a2,xl,wflat,glprod,m,ndriv,dg)
      return
c  Thick asymmetric flattop with different quartic + quadratic ends.
   15 w1=shape(icoil,1)
      w2=shape(icoil,2)
      s1=shape(icoil,3)
      s2=shape(icoil,4)
      a2=shape(icoil,5)
      call gtype15(z,a,a2,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      return
c  Thick asymmetric flattop with different cubic + quadratic ends.
   16 w1=shape(icoil,1)
      w2=shape(icoil,2)
      s1=shape(icoil,3)
      s2=shape(icoil,4)
      a2=shape(icoil,5)
      call gtype16(z,a,a2,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      return
c  Thick asymmetric flattop with different quartic + quadratic ends.
   17 continue
      w1=shape(icoil,1)
      w2=shape(icoil,2)
      s1=shape(icoil,3)
      s2=shape(icoil,4)
      a2=shape(icoil,5)
      b=shape(icoil,6)
      if(init.gt.0) go to 172
      z=0.d0
c     write(6,166) glprod,xL,a,a2,b
 166  format(' GL,xl,a,a2,b:',5f13.4)
c  write(6,167) w1,w2,s1,s2
 167  format(' w1,w2,s1,s2:',4f13.5)
      call gtype17(0,icoil,z,a,a2,b,xL,w1,w2,s1,s2,glprod,dg,
     & m,ndriv)
c write(6,*) 'NDRIV=',ndriv,'in ONAXGR after initialization'
      return
  172 continue
      call gtype17(1,icoil,z,a,a2,b,xL,w1,w2,s1,s2,glprod,dg,
     & m,ndriv)
      return
      end
c
c**************************************************
c
      subroutine bincof
c  calculates binomial coefficients
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /bicof/ bcoeff(maxcof,maxcof)
      dimension b(maxcof),bold(maxcof)
cryneneriwalstrom initialization of bcoeff moved to afro.f blockdata routine
cryneneriwalstrom Don't forget: if maxcof changes, also need to change in afro!
cryneneriwalstrom Eventually rewrite as a module
cryneneriwalstrom      data (bcoeff(k,1),k=1,2) /1.d0,1.d0/
cryneneriwalstrom      data (bcoeff(k,2),k=1,3) /1.d0,2.d0,1.d0/
      bold(1)=1.d0
      bold(2)=2.d0
      bold(3)=1.d0
      maxn=maxcof-1
      do 9 n=3,maxn
      b(1)=1.d0
      np1=n+1
      b(np1)=1.d0
      do 2 k=2,n
      km=k-1
    2 b(k)=bold(km)+bold(k)
      do 7 k=1,np1
      bcoeff(k,n)=b(k)
    7 bold(k)=b(k)
    9 continue
      return
      end
c
c**************************************************
c
      subroutine coeffs
      implicit double precision(a-h,o-z)
c calculates the coefficients
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)       Used in current-sheet magnets
c   and
c   cmp(m)=0.5*(2*m+1)*cm(m)                   Used for Halbach-type PMM
c   and
c   xmu0=mu0=4pi*10**-7                       Permeability of free space
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      cm(1)=0.5d0
      cmp(1)=0.75d0
      do 1 m=2,maxcof
      mm1=m-1
      cm(m)=cm(mm1)*0.5d0*dfloat(2*m-1)/dfloat(m)
      cmp(m)=0.5d0*cm(m)*dfloat(2*m+1)
    1 continue
      xmu0=8.d-7*dasin(1.d0)
      return
      end
c
c**************************************************
c
c   all active gtype routines
      subroutine gtype1(z,a,xl,eslope,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron,with an even quartic shape
c  function.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  eslope=dimensionless end slope of the shape function
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(Leff*mu0)
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  Leff=xl*(8+eslope)/15
c  mu0=4pi*1.d-7
      parameter (maxdrv=20,maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  typical value for eslope (for dipoles) is 1.
c  Shape function has form
c   f(x)=1+x**2*(eslope/2-2)*4/xl**2+x**4*(1-eslope/2)*16/xL**4
c     =1+bb*x**2+dd*x**4, -xl/2<x<+xl/2
c     =0, elsewhere.
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      m2=2*m
      xleff=xl*(8.d0+eslope)/15.d0
      c=glprod*cm(m)*a**m2/xleff
      hfl=0.5d0*xl
      hfl2=hfl**2
      hfl4=hfl2**2
      bb=(0.5d0*eslope-2.d0)/hfl2
      dd=(1.d0-eslope*0.5d0)/hfl4
      ai(1)=1.d0
      ai(2)=0.d0
      ai(3)=bb
      ai(4)=0.d0
      ai(5)=dd
      x1=-hfl
      x2=hfl
      call xngmin(m,5,ndriv,a,x1,x2,z,ai,xngi)
      do 1 n=1,ndriv
    1 dg(n)=xngi(n)*c
      return
      end
c
c**************************************************
c
      subroutine gtype2(z,a1,a2,xl,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  thick Halbach permanent multipole magnet without iron, of length xl.
c  Thedipole moment density depends on angle only and has components
c
c  M_r= -M_0 sin (m*phi)
c
c  M_phi= M_0 cos (m*phi)
c
c  M_0is given in terms of the integrated on-axis gradient glprod as
c  follows:
c   M_0 = -gLeff/(mu0*L*ln(a2/a1))   for m=1
c
c   M_0=-gLeff*(m-1)/(mu0*L*((1/a1)**(m-1)-(1/a2)**(m-1))   for m>1
c
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a1=inner magnet radius
c  a2=outer magnet radius
c  xl=coil length
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      parameter (maxdrv=20,maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      common /dnms/ denom(maxcof,2,2)
      dimension dg(maxdrv),dgdz(maxdrv)
c   cmp(m)=1*3*5*...(2m+1)/(m!*2**(m+1))
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Themagnetic moment per unit area has magnitude
c   M_0, -xl/2<x<+xl/2
c     =0, elsewhere.
c
c     The on-axis gradient for this type of magnet has the form
c
c     g_m(z)=mu0*M_0*(2m+1)!!/(m!*2**(m+1))       times
c  integral from a1 to a2 of the integral from -L/2-z to L/2-z  drho*ds
c
c     rho**(m+2)/((rho**2+s**2)**(m+3/2))
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  Must call BINCOF and COEFFS in an initialization call before using th
c    routine
c  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(m.gt.1) go to 44
c  m=1gradient constant
      c=glprod*cmp(1)/(xl*dlog(a2/a1))
      go to 45
c  m>1gradient constant
   44 mm1=m-1
      c=glprod*cmp(m)*dfloat(m-1)/(xl*(a1**(-mm1)-a2**(-mm1)))
c  Get0th derivative of g_m
   45 hfl=0.5d0*xl
      x1=-hfl-z
      x2=hfl-z
c  call g0hlb to get the above double integral
c
      call denoms(a1,a2,x1,x2,m,ndriv)
c
      call g0hlb(a1,a2,x1,x2,m,g0)
      dg(1)=-c*g0
      if(ndriv.lt.2) return
      s1=-x1
      s2=-x2
      ndrvm1=ndriv-1
c  Get1st and higher derivatives of g_m
      call dghlb(a1,a2,s1,s2,m,ndriv,dgdz)
      do 1 n=1,ndrvm1
      np1=n+1
    1 dg(np1)=c*dgdz(n)
      return
      end
c
c**************************************************
c
      subroutine gtype3(init,k,icoil,xl,xlmin,z,a,glprod,m,
     &ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calulates the on-axis gradient for both the fundament
c  andfirst harmonic of a Lambertson m-pole coil.
c  Called with init=0 to initialize, with init=1 for gradient.
c  k=0means fundamental
c  k=1means first harmonic
c  m=multipole index of harmonic in question.
c     i.e., if k=1, m=3*(m of fundamental)
      parameter(maxcoils=100,maxdrv=20,npoint=25)
      dimension aii(npoint,maxcoils),bii(npoint,maxcoils),
     &cii(npoint,maxcoils),amid(maxcoils),zn(npoint,maxcoils)
      dimension xi(npoint),yi(npoint),ai(npoint),bi(npoint)
      dimension ccoil(maxcoils),alfi(5),dg(maxdrv)
      dimension xngi(maxdrv),ci(npoint),points(npoint)
      data small/1.d-4/
      save aii,bii,cii,amid,zn,ccoil,npm1
      save small
      if(init.gt.0) go to 1
      dif=(xl-xlmin)/xl
      if(dabs(dif).lt.small) go to 1001
      m2=2*m
      hfl=xl*0.5d0
      hflmin=0.5d0*xlmin
      npm1=npoint-1
      dell=hfl-hflmin
      points(1)=0.d0
      points(npoint)=1.d0
      npm2=npoint-2
      dy=1.d0/dfloat(npm2)
      do 55 n=2,npm2
      y=dy*dfloat(n-1)
   55 points(n)=y**2*(2.d0-y)
      points(npm1)=0.5d0*points(npm2)+0.5d0
      xleff=0.5d0*(xl+xlmin)
      ccoil(icoil)=glprod*cm(m)*a**m2/xleff
      do 4 n=1,npoint
      x=dell*points(n)-hfl
      xi(n)=x
      zn(n,icoil)=x
      call flamb(xl,xlmin,x,fz,f3z)
      if(k.eq.0) yi(n)=fz
      if(k.eq.1) yi(n)=f3z
    4 continue
c  nowfit parabolae to xi,yi
      call parfit(npoint,xi,yi,ai,bi,ci)
      do 5 n=1,npm1
      aii(n,icoil)=ai(n)
      bii(n,icoil)=bi(n)
    5 cii(n,icoil)=ci(n)
      amid(icoil)=yi(npoint)
      return
c  Calculate gradients in subsequent calls.
    1 continue
      hflmin=0.5d0*xlmin
      c=ccoil(icoil)
      do 6 nd=1,ndriv
    6 dg(nd)=0.d0
      do 7 n=1,npm1
      np1=n+1
      x1=zn(n,icoil)
      x2=zn(np1,icoil)
      alfi(1)=aii(n,icoil)
      alfi(2)=bii(n,icoil)
      alfi(3)=cii(n,icoil)
      call xngmin(m,3,ndriv,a,x1,x2,z,alfi,xngi)
      do 8 nd=1,ndriv
    8 dg(nd)=dg(nd)+xngi(nd)
c  right hand side of coil is mirror image of left
      zm=-z
      call xngmin(m,3,ndriv,a,x1,x2,zm,alfi,xngi)
      sign=-1.d0
      do 9 nd=1,ndriv
      sign=-sign
    9 dg(nd)=dg(nd)+xngi(nd)*sign
    7 continue
c  getmiddle
      alfi(1)=amid(icoil)
      x1=-hflmin
      x2=hflmin
      call xngmin(m,1,ndriv,a,x1,x2,z,alfi,xngi)
      do 10 nd=1,ndriv
   10 dg(nd)=dg(nd)+xngi(nd)
      do 11 nd=1,ndriv
   11 dg(nd)=dg(nd)*c
      return
 1001 write(6,300)
  300 format(1x,'End zone of Lambertson too short-stopped in GTYPE3')
      stop
      end
c
c**************************************************
c
      subroutine gtype5(z,a,xl,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a flat shape function.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xl*mu0)
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      parameter (maxdrv=20)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c   f(x)=1, -xl/2<x<+xl/2
c     =0, elsewhere.
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      m2=2*m
      c=glprod*cm(m)*a**m2/xl
      hfl=0.5d0*xl
      ai(1)=1.d0
      x1=-hfl
      x2=hfl
      call xngmin(m,1,ndriv,a,x1,x2,z,ai,xngi)
      do 1 n=1,ndriv
    1 dg(n)=xngi(n)*c
      return
      end
c
c**************************************************
c
      subroutine gtype6(init,icoil,xl,glprod,ifile,a,z,m,
     &ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  Like gtype6, but with new fitting routine.
c  Also do not pass in np- np found in this routine by counting points.
c  This subroutine calulates the on-axis gradient for a user-defined
c  shape function.
c  Called with init=0 to initialize, with init=1 for gradient.
c   icoil=coil index
c   xl=coil length
c   glprod=is not the gradient-length product!
c   The user-defined shape function could have zero average value-
c  then glprod would be zero, but fields could be non-zero.
c  Forthis reason, we take glprod/xl instead of glprod/xleff in
c  determining the field constant c below.
c  ifile=integer specifying the data file containing the data for the
c  user-defined shape function
c  np=no. of lines in shape function data file.
c   a=coil radius
c   z=z coordinate of gradient evaluation point relative to coil center.
c   m=multipole index
c   ndriv=no. of derivatives of the gradient, plus 1.
c   dg=vector of ndriv numbers containing the on-axis gradient at z, plu
c     its ndriv-1 derivatives with respect to z.
c  Data for shape function are stored in disk file FORifile.dat, with
c   the format:
c  line 1:  x1,y1
c  line 2:  x2,y2
c   ...........
c  line np  xnp,ynp
c   Note: the points must be in order of increasing z.
c
      parameter(maxcoils=100,maxdrv=20,npoint=1001)
      dimension aii(npoint,maxcoils),bii(npoint,maxcoils),
     &cii(npoint,maxcoils),zn(npoint,maxcoils)
      dimension xi(npoint),yi(npoint),ci(npoint)
      dimension ccoil(maxcoils),alfi(5),dg(maxdrv)
      dimension xngi(maxdrv),npi(maxcoils),ai(npoint),bi(npoint)
      save aii,bii,cii,zn,npi
      if(init.gt.0) go to 1
      m2=2*m
      n=0
   89 continue
      read(ifile,*,end=90) x,y
      n=n+1
      xi(n)=x
      yi(n)=y
      zn(n,icoil)=xi(n)
      go to 89
   90 np=n
      write(6,*) 'np for icoil=',icoil,'=',np
      if(np.gt.1001) go to 1001
      npi(icoil)=np
c  Field constant.
      ccoil(icoil)=glprod*cm(m)*a**m2/xl
c  nowfit parabolae to xi,yi
c  newfitting routine
      call parfit1(np,xi,yi,ai,bi,ci)
      npm1=np-1
      do 5 n=1,npm1
      aii(n,icoil)=ai(n)
      bii(n,icoil)=bi(n)
    5 cii(n,icoil)=ci(n)
      return
c  Calculate gradients in subsequent calls.
    1 continue
      c=ccoil(icoil)
      npm1=npi(icoil)-1
      do 6 nd=1,ndriv
    6 dg(nd)=0.d0
      do 7 n=1,npm1
      np1=n+1
      h=zn(np1,icoil)-zn(n,icoil)
      zrel=z-zn(n,icoil)
      alfi(1)=aii(n,icoil)
      alfi(2)=bii(n,icoil)
      alfi(3)=cii(n,icoil)
      x1=0.d0
      x2=h
      call xngmin(m,3,ndriv,a,x1,x2,zrel,alfi,xngi)
      do 8 nd=1,ndriv
    8 dg(nd)=dg(nd)+xngi(nd)
    7 continue
      do 11 nd=1,ndriv
   11 dg(nd)=dg(nd)*c
      return
 1001 write(6,300)
  300 format(1x,'np>1001-stopped in GTYPE6')
      stop
      end
c
c**************************************************
c
      subroutine gtype7(z,a,xl,wflat,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arerepresented by a parabola tangent to the flat top.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  wflat=length of the flat section of the shape function.
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xleff*mu0)
c    where xleff=2/3*xl+1/3*wflat
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      parameter (maxdrv=20)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c   f(x)=1, -wflat/2<x<+wflat/2
c    f(x)=1-(wflat/2+x)**2/(xl/2-wflat/2)**2, -xl/2<x<-wflat/2
c    f(x)=1-(x-wflat/2)**2/(xl/2-wflat/2)**2, wflat/2<x<xl/2
c     =0, elsewhere.
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      data small/1.d-7/
      dif=(xl-wflat)/xl
      if(dabs(dif).lt.small) go to 1001
      o3rd=1.d0/3.d0
      t3rds=2.d0*o3rd
      m2=2*m
      xleff=t3rds*xl+o3rd*wflat
      c=glprod*cm(m)*a**m2/xleff
      hfw=0.5d0*wflat
      hfl=0.5d0*xl
      x1=-hfw
      x2=hfw
      ai(1)=1.d0
      call xngmin(m,1,ndriv,a,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=xngi(n)
      x1=hfw
      x2=hfl
      r=1.d0/(hfl-hfw)**2
      ai(1)=1.d0-hfw**2*r
      ai(2)=2.d0*hfw*r
      ai(3)=-1.d0*r
      call xngmin(m,3,ndriv,a,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+xngi(n)
      ai(2)=-ai(2)
      x1=-hfl
      x2=-hfw
      call xngmin(m,3,ndriv,a,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=c*(dg(n)+xngi(n))
      return
 1001 write(6,300)
  300 format(1x,'WFLAT in GTYPE7 was too close to xL-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype8(z,a,xl,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35,maxdrv=20)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface dipole sheet without iron, of length xl.
c  This is intended to model the large bore, weak, permanent magnet
c  multipoles.
c  Thedipole moment density has the components (as in Halbach PMMs)
c
c  M_r= -M_0 sin (m*phi)
c
c  M_phi= M_0 cos (m*phi)
c
c  M_0is given in terms of the integrated on-axis gradient glprod as
c  follows:
c   M_0 = -(gLeff * a**m )/ L
c
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      dimension dg(maxdrv),dipm1(maxdrv),dipm2(maxdrv)
c   cmp(m)=1*3*5*...(2m+1)/(m!*2**(m+1))
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Themagnetic moment per unit area has magnitude
c   M_0, -xl/2<x<+xl/2
c     =0, elsewhere.
c
c     The on-axis gradient for this type of magnet has the form
c
c     g_m(z)=((glprod * a**(m+2) * cmpm(m) )/ xl )*
c
c     Integral from -L/2-z to L/2-z of ds/((s**2+a**2)**(m+3/2))
c
c   cmp(m)=1*3*5*...(2m+1)/(m!*2**(m+1))
      m2p2=2*m+2
      c=glprod*cmp(m)*a**m2p2/xl
c  Get0th derivative of g_m
      hfl=0.5d0*xl
      x1=-hfl-z
      x2=hfl-z
c  call gm0int to get the integral from x1 to x2 of
c
c   1/((s**2+a**2)**(m+3/2))
c
      call pmmint(x1,x2,a,m,gm0int)
      dg(1)=-c*gm0int
      if(ndriv.lt.2) return
      ndrvm1=ndriv-1
c  Get1st and higher derivatives of g_m
      call dpmrv(m,ndrvm1,a,x1,dipm1)
      call dpmrv(m,ndrvm1,a,x2,dipm2)
      do 1 n=1,ndrvm1
      np1=n+1
      c=-c
    1 dg(np1)=c*(dipm1(n)-dipm2(n))
      return
      end
c
c**************************************************
c
      subroutine gtype9(z,a,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  Calculates the on-axis generalized gradient for a ZERO-integral-stren
c  ring multipole magnet.This magnet is composed of infintesimal dipoles
c  pointing in the z direction (plus or minus according to the angle) wi
c  a line density
c
c  M_z= M_0 sin (m*phi)
c
c
c  Since the integrated gradient is zero for this magnet type, we use a
c  BOGUS integrated strength to define M_0:
c
c  THEGLPROD IS DEFINED TO BE THE GLPROD THAT would ARISE IF ALL OF
c  THEDIPOLES WERE TO BE ROTATED 90 DEG. AROUND THE LOCAL PHI (TANGENT)
c  AXIS. This is 1/2 of the strength for a comparable Halbach ring, sinc
c  themagnitude of the dipole density is not constant, as in Halbach
c  magnets, but varies with angle as abs(sin m*phi)
c
c   M_0 = 2 * a**m * glprod / mu0
c
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  glprod=integral of on axis gradient
c  m=multipole index of ring density
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  glprod=A bogus integratedgradient-see above
c  mu0=4pi*1.d-7
      parameter (maxdrv=20,maxdp1=21,maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      dimension dg(maxdrv),di(maxdp1)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Thegradient has the form for this magnet type (with the ring at x=0)
c
c  g_m(z) =
c  glprod * (2m+1)!! * a**(2m+1) * z / [2**m * m! * (a**2+z**2)**(m+3/2)
c
      ndrvp1=ndriv+1
      m2p1=2*m+1
      n=m-1
      c=glprod*cm(m)*a**m2p1
      call dpmrv(n,ndrvp1,a,z,di)
      do 1 i=1,ndriv
    1 dg(i)=-c*di(i+1)
      return
      end
c
c**************************************************
c
      subroutine gtype10(z,a,xl,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  Calculates the on-axis generalized gradient for a ZERO-integral-stren
c  multipole magnet.This magnet is composed of infintesimal dipoles
c  pointing in the z direction (plus or minus according to the angle) wi
c  a constant surface density  (flat shape function)
c
c  M_z= M_0 sin (m*phi)
c
c
c  Since the integrated gradient is zero for this magnet type, we use a
c  BOGUS integrated strength to define M_0:
c
c  THEGLPROD IS DEFINED TO BE THE GLPROD THAT would ARISE IF ALL OF
c  THEDIPOLES WERE TO BE ROTATED 90 DEG. AROUND THE LOCAL PHI (TANGENT)
c  AXIS. This is 1/2 of the strength for a comparable Halbach sheet,sinc
c  themagnitude of the dipole density is not constant, as in Halbach
c  magnets, but varies with angle as abs(sin m*phi)
c
c   M_0 = 2 * a**m * glprod / (mu0 * L )
c
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=magnet length
c  glprod=integral of on axis gradient
c  m=multipole index of ring density
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  glprod=A bogus integrated gradient-see above
c  mu0=4pi*1.d-7
      parameter (maxdrv=20,maxdp1=21,maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      dimension dg(maxdrv),di(maxdp1)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Thegradient has the form for this magnet type (with ends at +- xl/2)
c
c  g_m(z) =
c     glprod * (2m-1)!! * a**(2m+1) * [ 1/D1 - 1/D2 ]  / [2**m * m! * xl
c   where
c
c    D1 = ( a**2 + (z-xl/2)**2 ) ** (m+1/2)
c    D2 = ( a**2 + (z+xl/2)**2 ) ** (m+1/2)
c
      m2p1=2*m+1
      c=glprod*cm(m)*a**m2p1/xl
      s=z-xl*0.5d0
      call rdrivs(m,ndriv,a,s,di)
      do 1 i=1,ndriv
    1 dg(i)=c*di(i)
      s=z+xl*0.5d0
      call rdrivs(m,ndriv,a,s,di)
      do 2 i=1,ndriv
    2 dg(i)=dg(i)-c*di(i)
      return
      end
c
c**************************************************
c
      subroutine gtype11(z,a,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arequartic in distance from flattop.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  w1=width of the curved section from -xl/2 to the flattop.
c  w2=width of the curved section from the flattop to xl/2.
c  glprod=integral of on-axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xleff*mu0)
c    where xleff=w1*(8+s1)/15+xl-w1-w2+w2*(8+s2)/15
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      parameter (maxdrv=20)
      parameter(fifteenth=1.d0/1.5d1,eight15th=8.d0/1.5d1)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c
c  f(x)=1+A1(x+xl/2-w1)**2 + B1*(x+xl/2-w1)**4  -xl/2<x<-xl/2+w1
c   where A1=(s1/2-2)/w1**2, B1=(1-s1/2)/w1**4  (Left-hand curved part)
c
c   f(x)=1, -xl/2+w1<x<xl/2-w2    (central flattop)
c
c  f(x)=1+A2(x-xl/2+w2)**2 + B2*(x-xl/2+w2)**4  xl/2-w2<x<xl/2
c   where A2=(s2/2-2)/w2**2, B2=(1-s2/2)/w2**4  (Right-hand curved part)
c
c   f(x)=0, elsewhere.
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      data small/1.d-9/
      xleff=0.d0
      m2=2*m
      hfl=0.5d0*xl
      dif=xl-w1-w2
      smal=-xl*small
      if(dif.lt.smal) go to 1001
      do 7 n=1,ndriv
    7 dg(n)=0.d0
      smal=-smal
      if(dabs(dif).le.smal) go to 70
c  Central flat section- skipped if w1+w2=xl.
      xleff=xl-w1-w2
      x1=-hfl+w1
      x2=hfl-w2
      ai(1)=1.d0
      call xngmin(m,1,ndriv,a,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=xngi(n)
c  Curved end sections.
   70 continue
c  Left hand curved section
      x1=-hfl
      x2=x1+w1
      w12=w1**2
      aa1=(0.5d0*s1-2.d0)/w12
      bb1=(1.d0-0.5d0*s1)/w12**2
      d1=hfl-w1
      d12=d1**2
      ai(1)=1.d0+d12*(aa1+d12*bb1)
      ai(2)=2.d0*d1*(aa1+2.d0*bb1*d12)
      ai(3)=aa1+6.d0*bb1*d12
      ai(4)=4.d0*bb1*d1
      ai(5)=bb1
      call xngmin(m,5,ndriv,a,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+xngi(n)
c  Right hand curved section
      x1=hfl-w2
      x2=hfl
      w22=w2**2
      aa2=(0.5d0*s2-2.d0)/w22
      bb2=(1.d0-0.5d0*s2)/w22**2
      d2=hfl-w2
      d22=d2**2
      ai(1)=1.d0+d22*(aa2+d22*bb2)
      ai(2)=-2.d0*d2*(aa2+2.d0*bb2*d22)
      ai(3)=aa2+6.d0*bb2*d22
      ai(4)=-4.d0*bb2*d2
      ai(5)=bb2
      call xngmin(m,5,ndriv,a,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+xngi(n)
      xleff=xleff+w1*(eight15th+s1*fifteenth)+
     &w2*(eight15th+s2*fifteenth)
      c=glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
 1001 write(6,300)
  300 format(1x,'w1+w2>xl in GTYPE11-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype12(z,a,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arequadratic + cubic in distance from flattop.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=coil radius
c  xl=coil length
c  w1=width of the curved section from -xl/2 to the flattop.
c  w2=width of the curved section from the flattop to xl/2.
c  glprod=integral of on-axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=glprod*a**mcoil/(xleff*mu0)
c    where xleff=w1*(8+s1)/15+xl-w1-w2+w2*(8+s2)/15
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
      parameter (maxdrv=20)
      parameter(twelveth=1.d0/1.2d1)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c
c  f(x)=1+A1(x+xl/2-w1)**2 + B1*(x+xl/2-w1)**3  -xl/2<x<-xl/2+w1
c   where A1=(s1-3)/w1**2, B1=(2-s1)/w1**3  (Left-hand curved part)
c
c   f(x)=1, -xl/2+w1<x<xl/2-w2    (central flattop)
c
c  f(x)=1+A2(x-xl/2+w2)**2 + B2*(x-xl/2+w2)**3  xl/2-w2<x<xl/2
c   where A2=(s2-3)/w2**2, B2=(2-s2)/w2**3  (Right-hand curved part)
c
c   f(x)=0, elsewhere.
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      data small/1.d-9/
      xleff=0.d0
      m2=2*m
      hfl=0.5d0*xl
      dif=xl-w1-w2
      smal=-xl*small
      if(dif.lt.smal) go to 1001
      do 7 n=1,ndriv
    7 dg(n)=0.d0
      smal=-smal
      if(dabs(dif).le.smal) go to 70
c  Central flat section- skipped if w1+w2=xl.
      xleff=xl-w1-w2
      x1=-hfl+w1
      x2=hfl-w2
      ai(1)=1.d0
      call xngmin(m,1,ndriv,a,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=xngi(n)
c  Curved end sections.
   70 continue
c  Left hand curved section
      x1=-hfl
      x2=x1+w1
      w12=w1**2
      aa1=(s1-3.d0)/w12
      bb1=(2.d0-s1)/(w12*w1)
      d1=hfl-w1
      d12=d1**2
      ai(1)=1.d0+d12*(aa1-d1*bb1)
      ai(2)=d1*(2.d0*aa1-3.d0*bb1*d1)
      ai(3)=aa1-3.d0*bb1*d1
      ai(4)=-bb1
      call xngmin(m,4,ndriv,a,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+xngi(n)
c  Right hand curved section
      x1=hfl-w2
      x2=hfl
      w22=w2**2
      aa2=(s2-3.d0)/w22
      bb2=(2.d0-s2)/(w22*w2)
      d2=hfl-w2
      d22=d2**2
      ai(1)=1.d0+d22*(aa2-d2*bb2)
      ai(2)=d2*(3.d0*bb2*d2-2.d0*aa2)
      ai(3)=aa2-3.d0*bb2*d2
      ai(4)=bb2
      call xngmin(m,4,ndriv,a,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+xngi(n)
      xleff=xleff+w1*(0.5d0+s1*twelveth)+
     &w2*(0.5d0+s2*twelveth)
      c=glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
 1001 write(6,300)
  300 format(1x,'w1+w2>xl in GTYPE12-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype13(z,a,a2,xl,eslope,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron,with an even quartic shape
c  function. Thick analog of gtype1.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=inner coil radius
c  a2 = outer coil radius.
c  xl=coil length
c  eslope=dimensionless end slope of the shape function
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, minus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(Leff*mu0)
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  typical value for eslope (for dipoles) is 1.
c  Shape function has form
c   f(x)=1+x**2*(eslope/2-2)*4/xl**2+x**4*(1-eslope/2)*16/xL**4
c     =1+bb*x**2+dd*x**4, -xl/2<x<+xl/2
c     =0, elsewhere.
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c  Leff=xl*(8+eslope)/15
c  mu0=4pi*1.d-7
      parameter (maxdrv=20,maxcof=35)
      parameter(maxgauss=100)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  Arrays for Gaussian quadrature. First index steps through Gauss point
c  second through orders.
      common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c  Determine no. of radial integration steps
      abar=0.5d0*(a+a2)
      adif=0.5d0*(a2-a)
      ratio=adif/a
      ngauss=4
      if(ratio.ge.0.1d0) ngauss=6
      if(ratio.gt.0.2d0) ngauss=8
      if(ratio.gt.0.3d0) ngauss=8
      if(ratio.gt.0.5d0) ngauss=12
      m2=2*m
      xleff=xl*(8.d0+eslope)/15.d0
      hfl=0.5d0*xl
      hfl2=hfl**2
      hfl4=hfl2**2
      bb=(0.5d0*eslope-2.d0)/hfl2
      dd=(1.d0-eslope*0.5d0)/hfl4
      ai(1)=1.d0
      ai(2)=0.d0
      ai(3)=bb
      ai(4)=0.d0
      ai(5)=dd
      x1=-hfl
      x2=hfl
      do 1 n=1,ndriv
    1 dg(n)=0.d0
      do 6 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  nowmultiply by field constant so that the integrated gradient
c  is the same as the specified gradient
      call thickf(a,a2,ratio,m)
      c=ratio*glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
      end
c
c**************************************************
c
      subroutine gtype14(z,a,a2,xl,wflat,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arerepresented by a parabola tangent to the flat top.
c  Thick analog of gtype7. Integration in radial depth by Gaussian
c  quadrature.
c  Assumes that NI/layer scales as ap, the radius of the layer,
c  while the shape function and Leff stay the same as the winding radius
c  increases. This is not quite right.
c  Since more turns are added to increase NI in the outer
c  layers, there is less room at the ends for crossover parts of the tur
c  Therefore, the shape function would be different, and Leff of the out
c  layers would be lower. For a better approximation of real windings,
c  multiple thick coils with increasing end zone widths w1 and w2 should
c  be used.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=inner coil radius
c  a2=outer coil radius
c  xl=coil length
c  wflat=length of the flat section of the shape function.
c  glprod=integral of on axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xleff*mu0)
c    where xleff=2/3*xl+1/3*wflat
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c   f(x)=1, -wflat/2<x<+wflat/2
c    f(x)=1-(wflat/2+x)**2/(xl/2-wflat/2)**2, -xl/2<x<-wflat/2
c    f(x)=1-(x-wflat/2)**2/(xl/2-wflat/2)**2, wflat/2<x<xl/2
c     =0, elsewhere.
c
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      parameter (maxdrv=20)
      parameter(maxcof=35,maxgauss=100,small=1.d-9)
      parameter(o3rd=1.d0/3.d0)
      parameter(t3rds=2.d0/3.d0)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  Arrays for Gaussian quadrature.  First index steps through Gauss poin
c  second through orders.
      common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
c
c  Determine no. of radial integration steps
      abar=0.5d0*(a+a2)
      adif=0.5d0*(a2-a)
      ratio=adif/a
      ngauss=4
      if(ratio.ge.0.1d0) ngauss=6
      if(ratio.gt.0.2d0) ngauss=8
      if(ratio.gt.0.3d0) ngauss=8
      if(ratio.gt.0.5d0) ngauss=12
      m2=2*m
      hfl=0.5d0*xl
      dif=xl-wflat
      smal=xl*small
      if(dabs(dif).lt.smal) go to 1001
      do 7 n=1,ndriv
    7 dg(n)=0.d0
      m2=2*m
      xleff=t3rds*xl+o3rd*wflat
      hfw=0.5d0*wflat
      hfl=0.5d0*xl
c  Central flat section.
      x1=-hfw
      x2=hfw
      ai(1)=1.d0
      do 4 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,1,ndriv,ap,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Curved end sections.
c  Right hand curved section
      x1=hfw
      x2=hfl
      r=1.d0/(hfl-hfw)**2
      ai(1)=1.d0-hfw**2*r
      ai(2)=2.d0*hfw*r
      ai(3)=-1.d0*r
      do 5 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,3,ndriv,ap,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Left hand curved section
      ai(2)=-ai(2)
      x1=-hfl
      x2=-hfw
      do 6 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,3,ndriv,ap,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  nowmultiply by field constant so that the integrated gradient
c  is the same as the specified gradient
      call thickf(a,a2,ratio,m)
      c=ratio*glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
 1001 write(6,300)
  300 format(1x,'wflat > or = xl in GTYPE14-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype15(z,a,a2,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arequadratic + quartic in distance from flattop.
c  Thick analog of gtype11. Integration in radial depth by Gaussian
c  quadrature.
c  Assumes that NI/layer scales as ap, the radius of the layer,
c  while the shape function and Leff stay the same as the winding radius
c  increases. This is not quite right.
c  Since more turns are added to increase NI in the outer
c  layers, there is less room at the ends for crossover parts of the tur
c  Therefore, the shape function would be different, and Leff of the out
c  layers would be lower. For a better approximation of real windings,
c  multiple thick coils with increasing end zone widths w1 and w2 should
c  be used.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=inner coil radius
c  a2=outer coil radius
c  xl=coil length
c  w1=width of the curved section from -xl/2 to the flattop.
c  w2=width of the curved section from the flattop to xl/2.
c   s1=Left hand end slope=df/dz(-xl/2)
c   s2=Right hand end slope=-df/dz(xl/2)
c  glprod=integral of on-axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xleff*mu0)
c    where xleff=w1*(8+s1)/15+xl-w1-w2+w2*(8+s2)/15
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c
c  f(x)=1+A1(x+xl/2-w1)**2 + B1*(x+xl/2-w1)**4  -xl/2<x<-xl/2+w1
c   where A1=(s1/2-2)/w1**2, B1=(1-s1/2)/w1**4  (Left-hand curved part)
c
c   f(x)=1, -xl/2+w1<x<xl/2-w2    (central flattop)
c
c  f(x)=1+A2(x-xl/2+w2)**2 + B2*(x-xl/2+w2)**4  xl/2-w2<x<xl/2
c   where A2=(s2/2-2)/w2**2, B2=(1-s2/2)/w2**4  (Right-hand curved part)
c
c   f(x)=0, elsewhere.
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      parameter (maxdrv=20)
      parameter(fifteenth=1.d0/1.5d1,eight15th=8.d0/1.5d1)
      parameter(third=1.d0/3.d0)
      parameter(maxcof=35,maxgauss=100,small=1.d-9)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  Arrays for Gaussian quadrature.  First index steps through Gauss poin
c  second through orders.
      common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
c
c Temporary++++++++
c  Determine no. of radial integration steps
      abar=0.5d0*(a+a2)
      adif=0.5d0*(a2-a)
      ratio=adif/a
      ngauss=4
      if(ratio.ge.0.1d0) ngauss=6
      if(ratio.gt.0.2d0) ngauss=8
      if(ratio.gt.0.3d0) ngauss=8
      if(ratio.gt.0.5d0) ngauss=12
      xleff=0.d0
      m2=2*m
      hfl=0.5d0*xl
      dif=xl-w1-w2
      smal=-xl*small
      if(dif.lt.smal) go to 1001
      do 7 n=1,ndriv
    7 dg(n)=0.d0
      smal=-smal
      if(dabs(dif).le.smal) go to 70
c  Central flat section- skipped if w1+w2=xl.
      xleff=xl-w1-w2
      x1=-hfl+w1
      x2=hfl-w2
      ai(1)=1.d0
      do 4 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,1,ndriv,ap,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Curved end sections.
   70 continue
c  Left hand curved section
      x1=-hfl
      x2=x1+w1
      w12=w1**2
      aa1=(0.5d0*s1-2.d0)/w12
      bb1=(1.d0-0.5d0*s1)/w12**2
      d1=hfl-w1
      d12=d1**2
      ai(1)=1.d0+d12*(aa1+d12*bb1)
      ai(2)=2.d0*d1*(aa1+2.d0*bb1*d12)
      ai(3)=aa1+6.d0*bb1*d12
      ai(4)=4.d0*bb1*d1
      ai(5)=bb1
      do 5 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Right hand curved section
      x1=hfl-w2
      x2=hfl
      w22=w2**2
      aa2=(0.5d0*s2-2.d0)/w22
      bb2=(1.d0-0.5d0*s2)/w22**2
      d2=hfl-w2
      d22=d2**2
      ai(1)=1.d0+d22*(aa2+d22*bb2)
      ai(2)=-2.d0*d2*(aa2+2.d0*bb2*d22)
      ai(3)=aa2+6.d0*bb2*d22
      ai(4)=-4.d0*bb2*d2
      ai(5)=bb2
      do 6 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
      xleff=xleff+w1*(eight15th+s1*fifteenth)+
     &w2*(eight15th+s2*fifteenth)
c  nowmultiply by field constant so that the integrated gradient
c  is the same as the specified gradient
      call thickf(a,a2,ratio,m)
      c=ratio*glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
 1001 write(6,300)
  300 format(1x,'w1+w2>xl in GTYPE15-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype16(z,a,a2,xl,w1,w2,s1,s2,glprod,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding without iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arequadratic + cubic in distance from flattop.
c  Thick analog of gtype12. Integration in radial depth by Gaussian
c  quadrature.
c  Assumes that NI/layer scales as ap, the radius of the layer,
c  while the shape function and Leff stay the same as the winding radius
c  increases. This is not quite right.
c  Since more turns are added to increase NI in the outer
c  layers, there is less room at the ends for crossover parts of the tur
c  Therefore, the shape function would be different, and Leff of the out
c  layers would be lower. For a better approximation of real windings,
c  multiple thick coils with increasing end zone widths w1 and w2 should
c  be used.
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=inner coil radius
c  a2=outer coil radius
c  xl=coil length
c  w1=width of the curved section from -xl/2 to the flattop.
c  w2=width of the curved section from the flattop to xl/2.
c   s1=Left hand end slope=df/dz(-xl/2)
c   s2=Right hand end slope=-df/dz(xl/2)
c  glprod=integral of on-axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  Forthis shape function, NI=2*glprod*a**mcoil/(xleff*mu0)
c    where xleff=w1*(8+s1)/15+xl-w1-w2+w2*(8+s2)/15
c  glprod=integral from z=-inf. to z=+inf. of g(z), g=on-axis generalize
c   gradient.
c  mu0=4pi*1.d-7
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c
c  f(x)=1+A1(x+xl/2-w1)**2 + B1*(x+xl/2-w1)**3  -xl/2<x<-xl/2+w1
c   where A1=(s1-3)/w1**2, B1=(2-s1)/w1**3  (Left-hand curved part)
c
c   f(x)=1, -xl/2+w1<x<xl/2-w2    (central flattop)
c
c  f(x)=1+A2(x-xl/2+w2)**2 + B2*(x-xl/2+w2)**3  xl/2-w2<x<xl/2
c   where A2=(s2-3)/w2**2, B2=(2-s2)/w2**3  (Right-hand curved part)
c
c   f(x)=0, elsewhere.
c
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      parameter (maxdrv=20)
      parameter(twelveth=1.d0/1.2d1)
      parameter(third=1.d0/3.d0)
      parameter(maxcof=35,maxgauss=100,small=1.d-9)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
c  Arrays for Gaussian quadrature.  First index steps through Gauss poin
c  second through orders.
      common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
c
c  Determine no. of radial integration steps
      abar=0.5d0*(a+a2)
      adif=0.5d0*(a2-a)
      ratio=adif/a
      ngauss=4
      if(ratio.ge.0.1d0) ngauss=6
      if(ratio.gt.0.2d0) ngauss=8
      if(ratio.gt.0.3d0) ngauss=8
      if(ratio.gt.0.5d0) ngauss=12
      xleff=0.d0
      m2=2*m
      hfl=0.5d0*xl
      dif=xl-w1-w2
      smal=-xl*small
      if(dif.lt.smal) go to 1001
      do 7 n=1,ndriv
    7 dg(n)=0.d0
      smal=-smal
      if(dabs(dif).le.smal) go to 70
c  Central flat section- skipped if w1+w2=xl.
      xleff=xl-w1-w2
      x1=-hfl+w1
      x2=hfl-w2
      ai(1)=1.d0
      do 4 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,1,ndriv,ap,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Curved end sections.
   70 continue
c  Left hand curved section
      x1=-hfl
      x2=x1+w1
      w12=w1**2
      aa1=(s1-3.d0)/w12
      bb1=(2.d0-s1)/(w12*w1)
      d1=hfl-w1
      d12=d1**2
      ai(1)=1.d0+d12*(aa1-d1*bb1)
      ai(2)=d1*(2.d0*aa1-3.d0*bb1*d1)
      ai(3)=aa1-3.d0*bb1*d1
      ai(4)=-bb1
      ai(5)=0.d0
      do 5 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Right hand curved section
      x1=hfl-w2
      x2=hfl
      w22=w2**2
      aa2=(s2-3.d0)/w22
      bb2=(2.d0-s2)/(w22*w2)
      d2=hfl-w2
      d22=d2**2
      ai(1)=1.d0+d22*(aa2-d2*bb2)
      ai(2)=d2*(3.d0*bb2*d2-2.d0*aa2)
      ai(3)=aa2-3.d0*bb2*d2
      ai(4)=bb2
      call xngmin(m,4,ndriv,a,x1,x2,z,ai,xngi)
      do 6 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 dg(n)=dg(n)+aratio*wgn(ng,ngauss)*xngi(n)
      xleff=xleff+w1*(0.5d0+s1*twelveth)+
     &w2*(0.5d0+s2*twelveth)
c  nowmultiply by field constant so that the integrated gradient
c  is the same as the specified gradient
      call thickf(a,a2,ratio,m)
      c=ratio*glprod*cm(m)*a**m2/xleff
      do 8 n=1,ndriv
    8 dg(n)=c*dg(n)
      return
 1001 write(6,300)
  300 format(1x,'w1+w2>xl in GTYPE16-stopped')
      stop
      end
c
c**************************************************
c
      subroutine gtype17(init,icoil,z,a,a2,b,xL,w1,w2,s1,s2,glprod,dg,
     & m,ndriv)
c  Like gtype17, except finer Gaussian quadrature.
c  This subroutine calculates the on-axis generalized gradient for a
c  cylindrical surface winding WITH iron, with a shape function
c  that is flat in the center and rounded at the ends. The ends
c  arequadratic + quartic in distance from flattop.
c
c  Coaxial mu=infinity iron cylinder surrounding a thick GTYPE15 coil.
c  Requires an initialization call to compute and store the equivalent i
c  shape function. This is the shape function for a fictitious
c  winding of radius b(=the iron radius) that mimics the incremental eff
c  of the iron for r<b.  The iron shape function is approximated
c  by a piece-wise-continuous cubic.
c  Thick analog of gtype11. Integration in radial depth by Gaussian
c  quadrature.
c  Assumes that NI/layer scales as ap, the radius of the layer,
c  while the shape function and Leff stay the same as the winding radius
c  increases. This is not quite right.
c  Since more turns are added to increase NI in the outer
c  layers, there is less room at the ends for crossover parts of the tur
c  Therefore, the shape function would be different, and Leff of the out
c  layers would be lower. For a better approximation of real windings,
c  multiple thick coils with increasing end zone widths w1 and w2 should
c  be used.
c
c  Usea Fourier integral method to find the iron shape function and its
c  z derivatives at the node points for the Piecewise-Continuous-Polynom
c  fit.
c
c  Need to link in the module ImKmpak (my portable modified Bessel-funct
c  package).
c
c  Subroutines called:
c
c  During initialization call: FBAR11,DIVIS,FILONIN,CUBINT,FITLENGTH,
c   FCONST,FITLENGTH
c  In gradient calls: XNGMIN.
c
c  Variables passed through common STUFF: ap,b,m
c   " "       "       "   GEOM : xL,w1,w2,s1,s2
c
c  Function names in EXTERNAL statement: FIRON,FDERIV.
c
c  Input  variables:
c
c  z=zcoordinate of field point with respect to coil center.
c  a=inner coil radius
c  a2=outer coil radius
c  b=iron radius.  To get an iron contribution, must have b>a2>a.  If b<
c     the iron contribution is omitted and the calculation
c     applies to a "bare" thick coil (just like GTYPE15).
c
c  xL=coil length
c  w1=width of the curved section from -xL/2 to the flattop.
c  w2=width of the curved section from the flattop to xL/2.
c   s1=Left hand end slope=df/dz(-xL/2)
c   s2=Right hand end slope=-df/dz(xL/2)
c  glprod=integral of on-axis gradient
c  m=multipole index of winding pattern-1=dipole, 2=quad.,etc.
c  ndriv=highest derivative of on-axis gradient, plus 1.
c  init=initialization code: 0 meams initialize, 1 means compute gradien
c
c
c  mu0=4pi*1.d-7
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
c
c  Output variables:
c  dg=vector with ndriv components containing the on-axis gradient, and
c     its ndriv-1 derivatives.
c
c  Shape function has form
c
c  f(x)=1+A1(x+xL/2-w1)**2 + B1*(x+xL/2-w1)**4  -xL/2<x<-xL/2+w1
c   where A1=(s1/2-2)/w1**2, B1=(1-s1/2)/w1**4  (Left-hand curved part)
c
c   f(x)=1, -xL/2+w1<x<xL/2-w2    (central flattop)
c
c  f(x)=1+A2(x-xL/2+w2)**2 + B2*(x-xL/2+w2)**4  xL/2-w2<x<xL/2
c   where A2=(s2/2-2)/w2**2, B2=(1-s2/2)/w2**4  (Right-hand curved part)
c
c   f(x)=0, elsewhere.
c
c   cm(m)=1*3*5*...(2m-1)/(m!*2**m)
      implicit double precision(a-h,o-z)
      double precision k1,k2,Lstar,Lstarn
      external firon,fderiv
      parameter (maxdrv=20)
      parameter (maxcoils=100)
      parameter(fifteenth=1.d0/1.5d1,eight15th=8.d0/1.5d1)
      parameter(third=1.d0/3.d0)
      parameter(mgauss=2)
      parameter(maxcof=35,maxgauss=100,small=1.d-9)
      dimension dg(maxdrv),ai(5),xngi(maxdrv)
      common /cmcof/ cm(maxcof),cmp(maxcof),xmu0
      common/stuff/ap,bp,mp
      common/geom1/xLp,w1p,w2p,s1p,s2p
c  Arrays for Gaussian quadrature.  First index steps through Gauss poin
c  second through orders.
      common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
      dimension yj(101),xi(103)
      dimension xj(101),ypj(101)
      dimension a1j(100),a2j(100),a3j(100),a4j(100)
c
      dimension a1i(103),a2i(103),a3i(103),a4i(103)
      dimension a1in(103,maxcoils), a2in(103,maxcoils),
     & a3in(103,maxcoils), a4in(103,maxcoils), xin(103,maxcoils),
     & nin(maxcoils), Lstarn(maxcoils)
      dimension ncoeff(103)
      dimension gbare(maxdrv),giron(maxdrv)
c  Save cubic coefficients, nodes, effective iron length, no. of nodes,
c  indexing them with coil index ICOIL.
      save a1in,a2in,a3in,a4in,xin,Lstarn,nin
c
      if(init.gt.0) go to 1
c
c  Initialization section
c  skip if b<a2
      if(b.le.a2)then
        return
      endif
c
c  Dummy variables for commons
      ap=a
      bp=b
      w1p=w1
      w2p=w2
      s1p=s1
      s2p=s2
      mp=m
      xLp=xL
c  Integration limits for the finite integral approx.
c  to the fourier transform of the shape function.
      xk1=40.d0/xL
      xk2=400.d0/xL
      if(dabs((b-a)/a).lt.1.e-3) go to 33
      k1=0.7d0*dfloat(m-1)/a+2.4d0/(b-a)
      k2=0.7d0*dfloat(m-1)/a+16.2d0/(b-a)
      if(k2.le.xk2) go to 34
   33 k2=xk2
      k1=xk1
   34 continue
c
c  Nowfind piecewise fit to iron shape function
c  First assign nodes.
c  First and last nodes found by linear extrapolation of function
c  to zero.
c  So fstar is represented by a cubic except for two end intervals
c  where it is linear, and equal zero at extreme L and R ends.
c xi are intervals for piecewise-continuous approx. to iorn
c  shape function
c xj =interior xi points, i.e. xj(1)=xi(2),xj(2)=xi(3), etc.
      call fbar11(xL,w1,w2,s1,s2,xbar)
      call divis(a,b,xL,w1,w2,s1,s2,xj,nj)
c     do 344 j=1,nj
c  344type *,'j=',j,' xj(j)=',xj(j)
      npoints=nj
      ni=nj+2
c  find fstar at knots.
      do 45 j=1,nj
c     x=xj(j)
c     type *,'j=',j,'xj(j)=',x
   45 xi(j+1)=xj(j)
c
c
c  Find iron shape function, and its derivatives, at nodes
c  by numerical Fourier integral, using Filon's quadrature method.
c  Need modified Bessel function object module ImKmpak here.
      abar=0.5d0*(a2+a)
      dela=0.5d0*(a2-a)
c  No.of quadrature steps/2
      nint=200
c
      ratio=dela/a
      ngauss=4
      if(ratio.ge.0.1d0) ngauss=6
      if(ratio.gt.0.2d0) ngauss=8
      if(ratio.gt.0.3d0) ngauss=8
      if(ratio.gt.0.5d0) ngauss=12
      do 44 j=1,nj
      z=xj(j)
      ffstar=0.d0
      dffstar=0.d0
      do 46 n=1,ngauss
      ap=abar+dela*xgn(n,ngauss)
      call filonin(0.d0,k1,z,nint,firon,sinint,cosint)
      ffstar=ffstar+ap*wgn(n,ngauss)*(sinint+cosint)
      call filonin(k1,k2,z,nint,firon,sinint,cosint)
      ffstar=ffstar+ap*wgn(n,ngauss)*(sinint+cosint)
c  Find first derivative of fstar at xj(j)
      call filonin(0.d0,k1,z,nint,fderiv,sinint,cosint)
      dffstar=dffstar+ap*wgn(n,ngauss)*(sinint+cosint)
      call filonin(k1,k2,z,nint,fderiv,sinint,cosint)
      dffstar=dffstar+ap*wgn(n,ngauss)*(sinint+cosint)
   46 continue
c  Gauss weights sum to 2 here- so divide by 2.
      ypj(j)=dffstar/(2.d0*a)
      yj(j)=ffstar/(2.d0*a)
   44 continue
c
c
      do 55 j=1,nj-1
      x1=xj(j)
      x2=xj(j+1)
      y1=yj(j)
      y2=yj(j+1)
      y1p=ypj(j)
      y2p=ypj(j+1)
      call cubint(x1,x2,y1,y2,y1p,y2p,aa1,aa2,aa3,aa4)
      a1j(j)=aa1
      a2j(j)=aa2
      a3j(j)=aa3
      a4j(j)=aa4
c     write(25,201) xj(j),yj(j)
   55 continue
      jcenter=nj/2
      ycenter=yj(jcenter)
      yc=yj(jcenter+1)
      if(yc.gt.ycenter) ycenter=yc
c     type *,'ycenter=',ycenter
c     write(25,201) xj(nj),yj(nj)
      if(ypj(1).le.0.d0) go to 15
      check=0.02d0*ycenter/a
      if(dabs(ypj(1)).lt.check) go to 15
      xi(1)=xi(2)-yj(1)/ypj(1)
c  LH linear function
      a1i(1)=-xi(1)*ypj(1)
      a2i(1)=ypj(1)
      a3i(1)=0.d0
      a4i(1)=0.d0
      go to 16
   15 continue
c     type *,'used alternate LH form'
c     type *,'ypj(1)=',ypj(1)
c  Alternate LH function
      xi(1)=xi(2)-2.d0*a
      a1i(1)=-xi(1)*yj(1)/(xi(2)-xi(1))
      a2i(1)=yj(1)/(xi(2)-xi(1))
      a3i(1)=0.d0
      a4i(1)=0.d0
   16 continue
      if(ypj(nj).ge.0.d0) go to 17
      check=0.02d0*ycenter/a
      if(dabs(ypj(nj)).lt.check) go to 17
      xi(ni)=xi(ni-1)-yj(nj)/ypj(nj)
c      RH end linear function A3=A4=0.
      a1i(ni-1)=-xi(ni)*ypj(nj)
      a2i(ni-1)=ypj(nj)
      a3i(ni-1)=0.d0
      a4i(ni-1)=0.d0
      go to 18
   17 continue
c     type *,'used alternate RH form'
c     type *,'ypj(nj)=',ypj(nj)
c  Alternate RH function
      xi(ni)=xi(ni-1)+2.d0*a
      a1i(ni-1)=xi(ni)*yj(nj)/(xi(ni)-xi(ni-1))
      a2i(ni-1)=-yj(nj)/(xi(ni)-xi(ni-1))
      a3i(1)=0.d0
      a4i(1)=0.d0
   18 continue
      ncoeff(1)=2
      ncoeff(ni-1)=2
c  Store nodes
      nin(icoil)=ni
c  interior intervals
c  Store cubic interpolant coefficients.
c
      do 47 j=1,nj-1
      a1i(j+1)=a1j(j)
      a2i(j+1)=a2j(j)
      a3i(j+1)=a3j(j)
      ncoeff(j+1)=4
   47 a4i(j+1)=a4j(j)
      do 48 i=1,ni-1
      a1in(i,icoil)=a1i(i)
      a2in(i,icoil)=a2i(i)
      a3in(i,icoil)=a3i(i)
      a4in(i,icoil)=a4i(i)
      xin(i,icoil)=xi(i)
   48 continue
      xin(ni,icoil)=xi(ni)
      call fitlength(Lstar,a1i,a2i,a3i,a4i,xi,ni)
      Lstarn(icoil)=Lstar
c     intv=ni-1
c     npts=ni
      write(6,*) 'Finished initializing GTYPE17:',icoil,'=icoil, using'
c     write(6,166) glprod,xL,a,a2,b
 166  format(' GL,xl,a,a2,b:',5f13.4)
c     write(6,167) w1,w2,s1,s2
 167  format(' w1,w2,s1,s2:',4f13.5)
      iniflg = 0
      return
c  Endof initialization section+++++++++++++++++++++++++++++++++++
c
c
c *************** Gradient calculation **************
c  Determine no. of radial integration steps
    1 continue
      abar=0.5d0*(a+a2)
      adif=0.5d0*(a2-a)
      ratio=adif/a
      ngauss=4*mgauss
      if(ratio.ge.0.1d0) ngauss=6*mgauss
      if(ratio.gt.0.2d0) ngauss=8*mgauss
      if(ratio.gt.0.3d0) ngauss=8*mgauss
      if(ratio.gt.0.5d0) ngauss=12*mgauss
      xLeff=0.d0
      m2=2*m
      hfl=0.5d0*xL
      dif=xL-w1-w2
      smal=-xL*small
      if(dif.lt.smal) go to 1001
      do 7 n=1,ndriv
    7 gbare(n)=0.d0
      smal=-smal
      if(dabs(dif).le.smal) go to 70
c  Central flat section- skipped if w1+w2=xL.
      xLeff=xL-w1-w2
      x1=-hfl+w1
      x2=hfl-w2
      ai(1)=1.d0
      do 4 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,1,ndriv,ap,x1,x2,z,ai,xngi)
      do 4 n=1,ndriv
    4 gbare(n)=gbare(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Curved end sections.
   70 continue
c  Left hand curved section
      x1=-hfl
      x2=x1+w1
      w12=w1**2
      aa1=(0.5d0*s1-2.d0)/w12
      bb1=(1.d0-0.5d0*s1)/w12**2
      d1=hfl-w1
      d12=d1**2
      ai(1)=1.d0+d12*(aa1+d12*bb1)
      ai(2)=2.d0*d1*(aa1+2.d0*bb1*d12)
      ai(3)=aa1+6.d0*bb1*d12
      ai(4)=4.d0*bb1*d1
      ai(5)=bb1
      do 5 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 5 n=1,ndriv
    5 gbare(n)=gbare(n)+aratio*wgn(ng,ngauss)*xngi(n)
c  Right hand curved section
      x1=hfl-w2
      x2=hfl
      w22=w2**2
      aa2=(0.5d0*s2-2.d0)/w22
      bb2=(1.d0-0.5d0*s2)/w22**2
      d2=hfl-w2
      d22=d2**2
      ai(1)=1.d0+d22*(aa2+d22*bb2)
      ai(2)=-2.d0*d2*(aa2+2.d0*bb2*d22)
      ai(3)=aa2+6.d0*bb2*d22
      ai(4)=-4.d0*bb2*d2
      ai(5)=bb2
      do 6 ng=1,ngauss
      ap=abar+xgn(ng,ngauss)*adif
      aratio=(ap/a)**(m+1)
      call xngmin(m,5,ndriv,ap,x1,x2,z,ai,xngi)
      do 6 n=1,ndriv
    6 gbare(n)=gbare(n)+aratio*wgn(ng,ngauss)*xngi(n)
      xLeff=xLeff+w1*(eight15th+s1*fifteenth)+
     &w2*(eight15th+s2*fifteenth)
c
c  Nowadd iron contribution (if b>a2)
      if(b.le.a2) go to 9
c  Find iron contribution to gradient
      do 10 n=1,ndriv
   10 giron(n)=0.d0
c  Step through intervals of piecewise-continuous cubic interpolant of
c  Fstar.
c  There are ni-1 intervals, ni nodes.
      ni=nin(icoil)
      do 11 i=1,ni-1
      ai(1)=a1in(i,icoil)
      ai(2)=a2in(i,icoil)
      ai(3)=a3in(i,icoil)
      ai(4)=a4in(i,icoil)
      ndeg=4
      if((i.eq.1).or.(i.eq.(ni-1))) ndeg=2
      x1=xin(i,icoil)
      x2=xin(i+1,icoil)
      call xngmin(m,ndeg,ndriv,b,x1,x2,z,ai,xngi)
      do 12 n=1,ndriv
   12 giron(n)=giron(n)+xngi(n)
   11 continue
    9 continue
c  calculate field constant C0=m*mu0*J0/2
c  This constant is chosen to make the integral gradient equal the
c  specified value, glprod.
c  Use"new" generalized gradient definition, including factor of m.
c  i.e, g_m(z)= m * lim as r goes to 0 of V_m(r,z)/r**m.
c  J0 is the winding current density at phi=0, z=0, and is assumed to
c  be constant with radial depth.
      denom=xLeff*aintgrl(a,a2,m)
      Lstar=Lstarn(icoil)
      if(b.gt.a2) denom=denom+adif*a*Lstar/b**m
      c0=glprod/denom
c  Nowadd up bare-coil and iron contributions to gradient with
c  appropriate weights.
      fact0=a**(m+1)*adif
      facti=b**m*a*adif
      c0cm=c0*cm(m)
      do 13 n=1,ndriv
      dg(n)=fact0*gbare(n)
      if(b.gt.a2) dg(n)=dg(n)+facti*giron(n)
   13 dg(n)=c0cm*dg(n)
      if(iniflg.eq.0) then
        iniflg = iniflg+1
c  write(6,297) fact0,facti,c0cm
c write(6,166) glprod,xL,a,a2,b
c rite(6,167) w1,w2,s1,s2
 297    format(' 1st grad, f0,fi,c)cm:',3f12.5)
      endif
      return
 1001 write(6,300)
  300 format(1x,'w1+w2>xL in GTYPE15-stopped')
      stop
      end
c
c**************************************************
c
c
c    rest of routines called by gtypes
c
      subroutine denoms(a1,a2,s1,s2,m,ndriv)
      implicit double precision(a-h,o-z)
c  Calculates the quantities
c
c   1/(rho**2+s**2)**(n+1/2))
c
c    for n=0,1,2,3,....m+ndriv-1
c     rho= a1,a2
c     s=   s1,s2
c
c   Results stored in common block DNMS
c
      parameter(maxcof=35)
      common /dnms/ denom(maxcof,2,2)
      dimension u(2,2)
c
c  First index in denoms is n+1, 2nd is 1 for a1, 2 for a2. 3rd index is
c   1 for s1, 2 for s2
c
c   ndriv is always equal to or greater than 1
      a12=a1**2
      a22=a2**2
      s12=s1**2
      s22=s2**2
      u(1,1)=1.d0/(a12+s12)
      u(1,2)=1.d0/(a12+s22)
      u(2,1)=1.d0/(a22+s12)
      u(2,2)=1.d0/(a22+s22)
      denom(1,1,1)=dsqrt(u(1,1))
      denom(1,1,2)=dsqrt(u(1,2))
      denom(1,2,1)=dsqrt(u(2,1))
      denom(1,2,2)=dsqrt(u(2,2))
      maxn=m+ndriv
      if(maxn.lt.2) return
      do 1 n=2,maxn
      nm1=n-1
      do 1 k=1,2
      do 1 l=1,2
    1 denom(n,k,l)=denom(nm1,k,l)*u(k,l)
c     do 5 k=1,maxn
c    5write(6,500) k,denom(k,1,2),denom(k,1,2),denom(k,2,1),
c    &denom(k,2,2)
c  500format(1x,'k,denom:',i2,4(1x,1pd14.7))
      return
      end
c
c**************************************************
c
      subroutine dghlb(a1,a2,s1,s2,m,ndriv,dg)
      implicit double precision(a-h,o-z)
c
c  Evaluates repeated derivatives with respect to z of the double
c  integral
c
c   s1to s2 ds   a1 to a2 drho of
c
c   rho**(m+2)/(rho**2+s**2)**(m+3/2)
c
c   s1=z+L/2
c   s2=z-L/2
c
c  Thefirst element of the output vector dg is the first derivative, et
c  dg has nd computed elements.
c
c   Note difference in definition from that of s1 and s2 in G0HLB- they
c   differ by a factor of -1.
c
c   Calls subroutine RHOINT
c
      parameter(maxdrv=20)
      parameter(small=1.d-6)
      dimension dg(maxdrv),ck(maxdrv),rhoin1(maxdrv),
     &rhoin2(maxdrv),ckold(maxdrv),pwr(maxdrv)
      dimension s1oddp(maxdrv),s2oddp(maxdrv),s1evnp(maxdrv),
     &s2evnp(maxdrv),xk(maxdrv),xkold(maxdrv)
      nd=ndriv-1
c     write(20,204) m,nd
c  204format(1x,'m=',i3,2x,'nd=',i3)
      call rhoint(a1,a2,s1,1,m,nd,rhoin1)
c     write(6,200) (rhoin1(k),k=1,nd)
  200 format(5(1x,1pd12.5))
      call rhoint(a1,a2,s2,2,m,nd,rhoin2)
c  First derivative
      dg(1)=rhoin2(1)-rhoin1(1)
      if(nd.lt.2) return
c  Second derivative
      nmin=m+1
      xkmin=dfloat(2*nmin+1)
      dg(2)=xkmin*(s1*rhoin1(2)-s2*rhoin2(2))
      if(nd.lt.3) return
c  Third derivative
      s12=s1**2
      s22=s2**2
      dg(3)=xkmin*(rhoin1(2)-rhoin2(2)+(xkmin+2.d0)*(s22*rhoin2(3)-
     &s12*rhoin1(3)))
      if(nd.lt.4) return
c  Fourth and higher derivatives.
c  Find out if nd is odd or even
      x=0.5d0*dfloat(nd)+small
      id=x
      ileft=0
      if((x-dfloat(id)).gt.0.1d0) ileft=1
c  ileft=0 if nd is even, =1 if nd is odd
c     write(20,205) id,ileft
c  205format(1x,'id=',i3,1x,'ileft=',i3)
      ckold(1)=-xkmin
      xkold(1)=xkmin
      xkold(2)=xkmin+2.d0
      ckold(2)=-ckold(1)*xkold(2)
      pwr(1)=0.d0
      pwr(2)=2.d0
      s1oddp(1)=s1
      s2oddp(1)=s2
      s1evnp(1)=1.d0
      s2evnp(1)=1.d0
      s1evnp(2)=s12
      s2evnp(2)=s22
      do 1 i=2,id
c     write(20,200) i
c     write(20,240)
  240 format(1x,'Even derivative Coefficients')
      im=i-1
      ip=i+1
c  first, even derivative in pair
      s1oddp(i)=s1oddp(im)*s12
      s2oddp(i)=s2oddp(im)*s22
      ng=2*i
      do 9 j=1,i
    9 xk(j)=xkold(j)+2.d0
      do 2 j=1,im
      jp=j+1
    2 ck(j)=-ckold(j)*xk(j)+ckold(jp)*pwr(jp)
      do 3 j=1,i
    3 pwr(j)=pwr(j)+1.d0
      xk(i)=xk(im)+2.d0
      ck(i)=-ckold(i)*xk(i)
      dg(ng)=0.d0
      do 4 j=1,i
      jcol=i+j
    4 dg(ng)=dg(ng)+(s2oddp(j)*rhoin2(jcol)-
     &s1oddp(j)*rhoin1(jcol))*ck(j)
c  Done with even derivative for this value of i
c  replace elements in ckold vector with elements in ck vector
      do 11 j=1,i
c     write(20,201) j,xk(j),pwr(j),ck(j),s1oddp(j),s2oddp(j)
   11 ckold(j)=ck(j)
c skipodd derivative calculation if not wanted- i.e. i=id & ileft=0
      if(i.lt.id) go to 5
      if(ileft.eq.0) go to 1
    5 continue
      s1evnp(ip)=s1evnp(i)*s12
      s2evnp(ip)=s2evnp(i)*s22
      ng=ng+1
      ck(1)=ckold(1)
      do 6 j=2,i
      jm=j-1
    6 ck(j)=-ckold(jm)*xk(j)+pwr(j)*ckold(j)
      xk(ip)=xk(i)+2.d0
      ck(ip)=-ckold(i)*xk(ip)
      do 7 j=1,i
    7 pwr(j)=pwr(j)-1.d0
      pwr(ip)=pwr(i)+2.d0
      dg(ng)=0.d0
      do 8 j=1,ip
      jcol=i+j
    8 dg(ng)=dg(ng)+(s2evnp(j)*rhoin2(jcol)-
     &s1evnp(j)*rhoin1(jcol))*ck(j)
      if(i.eq.id) go to 1
c     write(20,250)
c  250format(1x,'Odd Derivative Coefficients')
c  if not last value of i, load ck into ckold
      do 12 j=1,ip
c     write(20,201) j,xk(j),pwr(j),ck(j),s1evnp(j),s2evnp(j)
c  201format(1x,i3,5(1x,1pd14.7))
      xkold(j)=xk(j)
   12 ckold(j)=ck(j)
    1 continue
      return
      end
c
c**************************************************
c
      subroutine dpmrv(m,ndriv,a,s,dipm)
      implicit double precision(a-h,o-z)
c  Used in dipole sheet model of large bore PMMs.
c  This subroutine calculates repeated s derivatives of h_m=
c
c    1/(a**2+s**2)**(m+3/2))
c
c  dipm(1) = h_m, dipm(2)=d h_m /ds, dipm(3)= d**2 h_m /ds**2, etc.
c
c  Uses formula from Gradstein and Ryzhik, Table of Integrals,
c  Series, and Products, 1980, p.20
c
c  a=coil radius
c m=multipole index(m=1 for dipole,2 for quadrupole, 3 for sextupole,etc
c  s=x-z, where x=coil z coordinate, z= field point z coordinate.
c  ndriv= no. of entries in dipm (no. of times h_m is differentiated + 1
c
c  maxdrv is the maximum value of ndriv required.
      parameter (maxdrv=20,maxhlf=11)
c  maxhlf must be at least 1/2 of maxdrv if maxdrv is even,
c  1/2(maxdrv+1) if maxdrv is odd.
c
      dimension dipm(maxdrv)
      dimension cnk(maxhlf,maxdrv),
     &rnk2(maxhlf,maxdrv),snk(maxhlf,maxdrv),rkfac(maxhlf)
      aa=1.d0/a**2
      xm=dfloat(m)
      u=1.d0+aa*s**2
      ru=1.d0/u
      rru=dsqrt(ru)
      m2mm3=-2*m-3
      denom=ru**m*rru*a**m2mm3
      c2=denom*ru
      dipm(1)=c2
      if(ndriv.lt.2) return
      tas=2.d0*aa*s
      p2=-xm-1.5d0
      c2=c2*p2*ru
      dipm(2)=c2*tas
      if(ndriv.lt.3) return
      p2m=p2
      au=aa*u
      aunh=1.d0
      snk(1,1)=tas
c  n=order of derivative
c  find integers nhalf= nearest integer below or = (ndriv-1)/2,
c    and ileft=1+(ndriv-1)-2*nhalf
      xhalf=0.5d0*dfloat(ndriv-1)
      nhalf=xhalf
      dif=xhalf-dfloat(nhalf)
      ileft=1
      if(dif.gt.0.1d0) ileft=2
      n=1
      rkf=1.d0
      do 1 nh=1,nhalf
      rkf=rkf/dfloat(nh)
      rkfac(nh)=rkf
      nhp1=nh+1
      ilef=2
      if(nh.lt.nhalf) go to 6
      if(ileft.eq.1) ilef=1
    6 aunh=aunh*au
      nhp1=nh+1
      do 1 il=1,ilef
      nm1=n
      nm2=n-1
      n=n+1
      p2m=p2m-1.d0
      c2=c2*p2m*ru
      snk(nhp1,n)=aunh
      do 4 k=1,nh
    4 snk(k,n)=snk(k,nm1)*tas
      if(il.eq.2) snk(nhp1,n)=snk(nhp1,n)*tas
      cnk(2,n)=dfloat(n*(n-1))
      rnk2(2,n)=1.d0/p2m
      if(nhp1.lt.3) go to 8
      do 7 k=3,nhp1
      km1=k-1
      cnk(k,n)=cnk(km1,nm2)*cnk(2,n)
    7 rnk2(k,n)=rnk2(km1,nm1)*rnk2(2,n)
    8 sum2=snk(1,n)
      do 5 k=2,nhp1
      km1=k-1
    5 sum2=sum2+snk(k,n)*cnk(k,n)*rnk2(k,n)*rkfac(km1)
      np1=n+1
    1 dipm(np1)=sum2*c2
      return
      end
c
c**************************************************
c
      subroutine flamb(xl,xlmin,z,fz,f3z)
      implicit double precision(a-h,o-z)
c  This routine calculates the mth and 3mth(i.e., the two lowest) Fourie
c  coeffients of the stream function for an ideal Lambertson coil.
c  This type of Lambertson coil has equal z spacing of end crossover
c  turns. The stream function is dimensionless-i.e., the results
c  of this routine, fz, must be multiplied somewhere by NI,
c  where NI is the number of amp-turns/pole.
c
c  Input variables:
c  xl=coil length
c  xlmin=length of shortest turn (at center of pole)
c  z=axial coordinate measured from coil center
c
c Output
c  fz=mth Fourier coefficient of the dimensionless stream function
c  fora Lambertson m-pole coil. (All lower coefficients are zero).
c  f3z=3m th Fourier coefficient  "     "    "     "
c
c   Note that the shape functions f(z) and f3z(z) have no m dependence.
      data two/2.d0/,three/3.d0/,half/0.5d0/,zero/0.d0/
      data five/5.d0/,four/4.d0/
      data one/1.d0/,small/1.d-12/
      hfl=half*xl
      zet=dabs(z)
      fz=zero
      f3z=zero
      if(zet.gt.hfl) return
      dif=dabs(xl-xlmin)
      if(dif.gt.small) go to 2
      fz=one
      return
    2 hflmin=half*xlmin
      hfpi=dasin(one)
      xl2=xl**2
      qtrpi=half*hfpi
      xlmin2=xlmin**2
      third=one/three
      b=(xl2-xlmin2)/xl2
      b2=b**2
      c=xl/((xl-xlmin)*qtrpi)
      xkc2=(one-b)/(one+b)
      xkc=dsqrt(xkc2)
c  find  elliptic integrals E(pi/4,r),F(pi/4,r), where r=sqrt(2b/(1+b)).
      ee=el2(one,xkc,one,xkc2)
      ff=el2(one,xkc,one,one)
      rootab=dsqrt(one+b)
      sixth=half*third
      g2f=(one-b)*sixth/b
      g2e=(three*b-one)*sixth/b
      g2cos=-sixth/rootab
      denom=one/(30.d0*b2)
      g4f=-(five*b2-four*b-one)*denom
      g4e=(12.d0*b2-5.d0*b-one)*denom
      g4cos=-one/(60.d0*b*rootab)
      g6g2=-three*(one+b)/(14.d0*b)
      g6g4=(two+8.d0*b)/(7.d0*b)
      g6cos=one/(56.d0*b*rootab)
      dif=(zet-hflmin)/hflmin
      if(dif.gt.small) go to 1
c  z falls in center part of coil, ie. between crossover regions at ends
      dele=cel(xkc,one,one,xkc2)-ee
      delf=cel(xkc,one,one,one)-ff
      cm=third*(one+two*(rootab*dele-(one-b2)*delf/rootab)/b)
      fz=cm*c
      g0=dele
      g2=g2f*delf+g2e*dele-g2cos
      g4=g4f*delf+g4e*dele-g4cos*(10.d0*b-one)
      g6=g6g2*g2+g6g4*g4-g6cos
      c3m=third-two*rootab*(g0-18.d0*g2+48.d0*g4-32.d0*g6)
      f3z=c*c3m
      return
    1 continue
c  z falls in one of 2 crossover regions
      zr=zet/hfl
      zr2=zr**2
      sinfz=(one-zr2)/b
      cosfz=dsqrt(one-sinfz**2)
      x=dsqrt((one+sinfz)/(one-sinfz))
      dele=el2(x,xkc,one,xkc2)-ee
      delf=el2(x,xkc,one,one)-ff
      fz=third*(two*(rootab*dele-(one-b**2)*delf/rootab)/b+
     &one-cosfz*zr)
      fz=fz*c
      cos3fz=cosfz*(one-four*sinfz**2)
      g0=dele
      g2=g2cos*(cosfz*zr-one)+g2e*dele+g2f*delf
      g4=g4cos*((10.d0*b-one+three*b*sinfz)*cosfz*zr-
     &10.d0*b+one)+g4e*dele+g4f*delf
      g6=(cosfz*(one+sinfz)*zr*zr2-one)*g6cos+g6g2*g2+g6g4*g4
      f3z=third*(one-zr*cos3fz)-two*rootab*(g0-18.d0*g2+
     &48.d0*g4-32.d0*g6)
      f3z=f3z*c
      return
      end
c
c**************************************************
c
c
c
c     this is a function subprogram copied from numerical recipes
c     by w.h.press ,etl pp. 186-187.
c     evaluates the general incomplete elliptic integral of
c     the 2nd kind as the following function of four variables:
c
c     el2(x,kc,a,b)=
c     int[ (a+b*x**2)dx/{1+x**2)*sqrt((1+x**2)*(1+kc**2*x**2))} ]
c     from x=0 to x=x
c
c     the elliptic integral of the 1st kind is called by el2(x,kc,1,1)
c     the elliptic integral of the 2nd kind is called by el2(x,kc,1,kc**
c
      subroutine g0hlb(a1,a2,s1,s2,m,g0)
      implicit double precision(a-h,o-z)
      parameter(small=1.d-6)
      parameter(maxcof=35)
      dimension zinti(maxcof,2)
c   The array ZINTI contains the integrals from s1 to s2 of
c
c    ds/(a**2+s**2)**n+1/2, n=nmin,nmin+1,nmin+2,...m
c
c   zinti(i,1) is evaluated at a1, zinti(i,2) at a2.
c   The index i runs from 1 to m-nmin+1, with n=nmin for i=1, nmin+1 for
c   i=2, etc.
c
c
c
c   Used in calculating the gradient on axis of thick Halbach PMMs.
c   Evaluates the double integral, s=s1 to s2, rho=a1 to a2 of
c
c     (drho*ds*rho**(m+2))/(rho**2+s**2)**(m+3/2)
c
c   This is the gradient on axis , except for a constant factor
c
c   a1=inner radius
c   a2=outer radius
c   s1=-L/2-z
c   s2=L/2-z
c     where L is the length of the PMM, assumed to be rectangular in cro
c     section in the r-z plane, and z is the axial coordinate where the
c     gradient is evaluated.
c   m=multipole index:  m=1 for dipole, 2 for quad, etc.
c   g0=value of double integral
c   Uses repeated integration by parts to do the rho integral
c   Each term in resultant series is integrated in z, except if m is eve
c   Then the series ends in a double integral that is evaluated by the
c   special-case subroutine INTEND.
c  Check to see if m is even or odd.
      if(m.lt.2) go to 8
      mm1=m-1
      ra1=a1**(-mm1)
      ra2=a2**(-mm1)
      go to 9
    8 ra1=1.d0
      ra2=1.d0
    9 x=dfloat(m+3)*0.5d0+small
      iterms=x
      dif=x-dfloat(iterms)
      ileft=0
      if(dif.gt.0.1d0) ileft=1
c  ileft=0 if m is odd, =1 if m is even
      if(ileft.gt.0) go to 11
c   m is odd
      itop=0
      ibot=m-2
c  Exponent on denominator in innermost term=n+1/2
      nmin=(m-1)/2
      g0=0.d0
      go to 1
   11 continue
c     m is even- descending series ends in an integral evaluated by INTE
      itop=1
      ibot=m-1
      nmin=m/2
      call intend(a1,a2,s1,s2,m,xiend)
      g0=-xiend
    1 itrm1=iterms-1
      call zint(a1,a2,s1,s2,nmin,m,zinti)
      do 2 i=1,itrm1
      if(m.lt.2) go to 10
      g0=g0+zinti(i,2)*ra2-zinti(i,1)*ra1
      go to 14
   10 g0=g0+zinti(i,2)-zinti(i,1)
   14 itop=itop+2
      ibot=ibot+2
      g0=g0*dfloat(itop)/dfloat(ibot)
    2 continue
      ibot=ibot+2
      if(m.lt.2) go to 12
      g0=g0+zinti(iterms,2)*ra2-zinti(iterms,1)*ra1
      go to 13
   12 g0=g0+zinti(iterms,2)-zinti(iterms,1)
   13 g0=-g0/dfloat(ibot)
      return
      end
c
c**************************************************
c
      subroutine parfit(npoint,xi,yi,ai,bi,ci)
      implicit double precision(a-h,o-z)
      dimension xi(npoint),yi(npoint),ai(npoint),bi(npoint),
     &ci(npoint)
c fit parabola to 1st 3 points.
      x1=xi(1)
      y1=yi(1)
      x2=xi(2)
      y2=yi(2)
      x3=xi(3)
      y3=yi(3)
      h1=x2-x1
      h2=x3-x2
      d1=y1/(h1*(h1+h2))
      d2=y2/(h1*h2)
      d3=y3/(h2*(h1+h2))
c  parabola of form a1+b1*x+c1*x**2
      a1=x1*x2*d3+x2*x3*d1-x1*x3*d2
      b1=(x1+x3)*d2-(x2+x3)*d1-(x1+x2)*d3
      c1=d1-d2+d3
      ai(1)=a1
      bi(1)=b1
      ci(1)=c1
      npm2=npoint-2
c  except for ends, the parabola coefficients for a particular interval
c  theaverage of those calculated with the point on the right, with tho
c  calculated with the point on the left.
      do 1 i=2,npm2
      ip2=i+2
      x1=x2
      x2=x3
      y1=y2
      y2=y3
      h1=h2
      x3=xi(ip2)
      y3=yi(ip2)
      h2=x3-x2
      d1=y1/(h1*(h1+h2))
      d2=y2/(h1*h2)
      d3=y3/(h2*(h1+h2))
      a2=x1*x2*d3+x2*x3*d1-x1*x3*d2
      b2=(x1+x3)*d2-(x2+x3)*d1-(x1+x2)*d3
      c2=d1-d2+d3
      ai(i)=0.5d0*(a1+a2)
      bi(i)=0.5d0*(b1+b2)
      ci(i)=0.5d0*(c1+c2)
      a1=a2
      b1=b2
    1 c1=c2
      npm1=npoint-1
c  rightmost interval
      ai(npm1)=a2
      bi(npm1)=b2
      ci(npm1)=c2
      return
      end
c
c**************************************************
c
      subroutine parfit1(npoint,xi,yi,ai,bi,ci)
      implicit double precision(a-h,o-z)
      dimension xi(npoint),yi(npoint),ai(npoint),bi(npoint),
     &ci(npoint)
c  This subroutine finds coefficients for a piecewise-quadratic fit to a
c  of xi,yi points. There are npoint (xi,yi) points and npoint-1 interva
c  also npoint-1 ai values, npoint-1 bi values, npoint-1 ci values.
c  fitto ith interval is ai(i)+bi(i)*t+ci(i)*t**2, where t=x-xi(i).
c  Like PARFIT, but coefficients apply to x-xi(i), not x.  Roundoff erro
c  should be lower.
c fit parabola to 1st 3 points.
      y1=yi(1)
      y2=yi(2)
      y3=yi(3)
      h1=xi(2)-xi(1)
      h2=xi(3)-xi(2)
      denom=h1*h2*(h1+h2)
      aa=y1
      denom=h1*h2*(h1+h2)
      bb=((y2-y1)*(h1+h2)**2-(y3-y1)*h1**2)/denom
      cc=(h1*(y3-y2)-h2*(y2-y1))/denom
c  parabola of form a+b*t+c*t**2
c  where t=x-x1
      ai(1)=aa
      bi(1)=bb
      ci(1)=cc
c  except for ends, the parabola coefficients for a particular interval
c  theaverage of those calculated with the point on the right, with tho
c  calculated with the point on the left.
      if(npoint.lt.4) go to 2
      do 1 i=2,npoint-2
      h3=xi(i+2)-xi(i+1)
      y4=yi(i+2)
c  LH expression
      aa=y2
      denom=h2*h3*(h2+h3)
      bb=((y3-y2)*(h2+h3)**2-(y4-y2)*h2**2)/denom
      cc=(h2*(y4-y3)-h3*(y3-y2))/denom
c  RH expression
      denom=h1*h2*(h1+h2)
      a=y2
      b=(h2**2*(y2-y1)+h1**2*(y3-y2))/denom
      c=(h1*(y3-y2)-h2*(y2-y1))/denom
c  Average two expressions
      ai(i)=0.5d0*(a+aa)
      bi(i)=0.5d0*(b+bb)
      ci(i)=0.5d0*(c+cc)
c  Shift fit window, unless at RH end of array.
      if(i.ne.(npoint-2)) then
      h1=h2
      y1=y2
      h2=h3
      y2=y3
      y3=y4
      endif
    1 continue
    2 npm1=npoint-1
c  rightmost interval
      denom=h2*h3*(h2+h3)
      a=y3
      b=(h3**2*(y3-y2)+h2**2*(y4-y3))/denom
      c=(h2*(y4-y3)-h3*(y3-y2))/denom
      ai(npm1)=a
      bi(npm1)=b
      ci(npm1)=c
      ai(npoint)=0.d0
      bi(npoint)=0.d0
      ci(npoint)=0.d0
      return
      end
c
c**************************************************
c
      subroutine pmmint(x1,x2,a,m,gm0int)
      implicit double precision(a-h,o-z)
c  this subroutine calculates the integral from x1 to x2 of
c
c   (s**2+a**2)**(-m-3/2) ds
c
c  Uses formula for Gradstein and Ryzhik, p. 86
      parameter(maxcof=35)
      common /bicof/ bcoeff(maxcof,maxcof)
c     bcoeff contains the binomial expansion coefficients up to order
c   maxcof-1
      mp1=m+1
      gm0int=0.d0
      a2=a**2
      x=x2
      do 1 i=1,2
      xsq=x**2
      root2=xsq+a2
      root=dsqrt(root2)
      term=x/root
      temp=bcoeff(1,m)*term
      do 2 k=2,mp1
      term=-term*xsq/root2
      x2km1=dfloat(2*k-1)
    2 temp=temp+bcoeff(k,m)*term/x2km1
      gm0int=gm0int+temp
    1 x=-x1
      gm0int=gm0int*a2**(-mp1)
      return
      end
c
c**************************************************
c
      subroutine rdrivs(m,ndriv,a,s,di)
      implicit double precision(a-h,o-z)
c  This subroutine calculates repeated derivatives of g_m=
c
c    1/(a**2+s**2)**(m+1/2))
c
c  di(1) = g_m, di(2)=d g_m /ds, di(3)= d**2 g_m /ds**2, etc.
c
c  Uses formula from Gradstein and Ryzhik, Table of Integrals,
c  Series, and Products, 1980, p.20
c
c  a=coil radius
c  m=multipole index (m=1 for dipole,2 for quadrupole, 3 for sextupole,e
c  s=x-z, where x=coil z coordinate, z= field point z coordinate.
c  ndriv= no. of entries in di (no. of times g_m is differentiated + 1)
c
c  maxdrv is the maximum value of ndriv required.
      parameter (maxdrv=20,maxhlf=11)
c  maxhlf must be at least 1/2 of maxdrv if maxdrv is even,
c  1/2(maxdrv+1) if maxdrv is odd.
c
      dimension di(maxdrv)
      dimension cnk(maxhlf,maxdrv),rnk1(maxhlf,maxdrv),
     &snk(maxhlf,maxdrv),rkfac(maxhlf)
      aa=1.d0/a**2
      xm=dfloat(m)
      u=1.d0+aa*s**2
      ru=1.d0/u
      rru=dsqrt(ru)
      m2mm1=-2*m-1
      denom=ru**m*rru*a**m2mm1
      c1=xm*denom
      di(1)=c1
      if(ndriv.lt.2) return
      tas=2.d0*aa*s
      p1=-xm-0.5d0
      c1=c1*p1*ru
      di(2)=c1*tas
      if(ndriv.lt.3) return
      p1m=p1
      au=aa*u
      aunh=1.d0
      snk(1,1)=tas
c  n=order of derivative
c  find integers nhalf= nearest integer below or = (ndriv-1)/2,
c    and ileft=1+(ndriv-1)-2*nhalf
      xhalf=0.5d0*dfloat(ndriv-1)
      nhalf=xhalf
      dif=xhalf-dfloat(nhalf)
      ileft=1
      if(dif.gt.0.1d0) ileft=2
      n=1
      rkf=1.d0
      do 1 nh=1,nhalf
      rkf=rkf/dfloat(nh)
      rkfac(nh)=rkf
      nhp1=nh+1
      ilef=2
      if(nh.lt.nhalf) go to 6
      if(ileft.eq.1) ilef=1
    6 aunh=aunh*au
      nhp1=nh+1
      do 1 il=1,ilef
      nm1=n
      nm2=n-1
      n=n+1
      p1m=p1m-1.d0
      c1=c1*p1m*ru
      snk(nhp1,n)=aunh
      do 4 k=1,nh
    4 snk(k,n)=snk(k,nm1)*tas
      if(il.eq.2) snk(nhp1,n)=snk(nhp1,n)*tas
      cnk(2,n)=dfloat(n*(n-1))
      rnk1(2,n)=1.d0/p1m
      if(nhp1.lt.3) go to 8
      do 7 k=3,nhp1
      km1=k-1
      cnk(k,n)=cnk(km1,nm2)*cnk(2,n)
    7 rnk1(k,n)=rnk1(km1,nm1)*rnk1(2,n)
    8 sum1=snk(1,n)
      do 5 k=2,nhp1
      km1=k-1
    5 sum1=sum1+snk(k,n)*cnk(k,n)*rnk1(k,n)*rkfac(km1)
      np1=n+1
    1 di(np1)=sum1*c1
      return
      end
c
c**************************************************
c
      subroutine thickf(a,a2,ratio,m)
      implicit double precision(a-h,o-z)
c  calculates ratio of thin to thick coils for correction of field
c  constant to get desired glprod in gtype 13,14,15,16.
c  factor of 1/2 comes from the fact that sum of wgn is 2.
c  RATIO is  1/2 * (a2-a)/a**(m-1) / integral from a to a2 of dx*x**(1-m
c  this comes from the fact that if NI varies as x, and gL varies as
c   x**(-m), the product varies as x**(1-m) (x=dummy coil radius variabl
c  of integration).  For thin coils, RATIO approcahes 1/2.  For m=1, it
c  EXACTLY 1/2, for any thickness.
      parameter(c1=0.5d0,c2=0.5d0,c3=1.d0/6.d0,c4=1.d0/6.d0,
     &c5=1.9d1/9.0d1,c6=0.3d0,c7=8.63d2/1.89d3,c8=2.75d2/3.78d2)
      parameter(third=1.d0/3.d0)
      if(m.gt.1) go to 22
      ratio=0.5d0
      go to 50
   22 x=0.5d0*(a2-a)/a
      if(x.lt.0.002d0) go to 44
      if(m.gt.2) go to 23
c  m=2, not thin
      ratio=0.5d0*(a2-a)/(a*dlog(a2/a))
      go to 50
   23 if(m.gt.3) go to 24
c  m=3, not thin
      ratio=0.5d0*a2/a
      go to 50
c  m>3, not thin
   24 ratio=x*dfloat(m-2)/
     &(1.d0-(a/a2)**(m-2))
      go to 50
   44 continue
c  Thin-small value of x
      x=0.5d0*(a2-a)/a
      if(m.gt.2) go to 25
c  thin-m=2
      ratio=c1+x*(c2-x*(c3-x*(c4-x*(c5-x*c6))))
      go to 50
   25 if(m.gt.3) go to 26
c  Thin- m=3
      ratio=0.5d0*a2/a
      go to 50
c  Thin- m>3
   26 ratio=0.5d0*(1.d0+dfloat(m-1)*x*(1.d0+third*x*
     &dfloat(m-3)*(1.d0-x)))
   50 continue
      return
      end
c
c**************************************************
c
      subroutine xngmin(m,nterm,ndriv,a,x1,x2,z,ai,xngi)
      implicit double precision(a-h,o-z)
c  Evaluates analytical expressions for
c
c   d**n/ dz**n of Integral from x1 to x2 of f(x)g_m(x,z) dx,  where
c
c  f(x)=ai(1)+ai(2)*x+..ai(nterm)*x**(nterm-1)
c
c   g_m=a**(1-m) * d/da [a**m/(a**2+(x-z)**2)**(m+1/2)]
c
c for n=0 to ndriv-1.
c   Maximum value of nterm is 5.
      parameter(maxdrv=20)
      dimension c0(5,5),c1(4,4),c2(3,3),bi(5),ai(5)
c     dimension c3(2,2)
      dimension xgi1(5),xgi2(5),xngi(maxdrv),di1(maxdrv),di2(maxdrv)
      data (c0(k,1),k=1,5) /5*1.d0/
      data (c0(k,2),k=1,5) /0.d0,1.d0,2.d0,3.d0,4.d0/
      data (c0(k,3),k=1,5) /2*0.d0,1.d0,3.d0,6.d0/
      data (c0(k,4),k=1,5) /3*0.d0,1.d0,4.d0/
      data (c0(k,5),k=1,5) /4*0.d0,1.d0/
      data (c1(k,1),k=1,4) /1.d0,2.d0,3.d0,4.d0/
      data (c1(k,2),k=1,4) /0.d0,2.d0,6.d0,12.d0/
      data (c1(k,3),k=1,4) /2*0.d0,3.d0,12.0/
      data (c1(k,4),k=1,4) /3*0.d0,4.d0/
      data (c2(k,1),k=1,3) /2.d0,6.d0,12.d0/
      data (c2(k,2),k=1,3) /0.d0,6.d0,24.d0/
      data (c2(k,3),k=1,3) /2*0.d0,12.d0/
c     data (c3(k,1),k=1,2) /6.d0,24.d0/
c     data (c3(k,2),k=1,2) /0.d0,24.d0/
c     data c4/24.d0/
      save c0,c1,c2
      if(nterm.gt.5) go to 1001
      if(nterm.lt.1) go to 1002
      s1=x1-z
      s2=x2-z
      call xngint(m,nterm,a,s1,xgi1)
      call xngint(m,nterm,a,s2,xgi2)
c  Find coefficients for the expansion of f(x)=f(z+s) in powers of s.
      do 1 n=1,nterm
      kmax=nterm-n
      bi(n)=ai(nterm)*c0(nterm,n)
      if(kmax.lt.1) go to 1
      do 2 k=1,kmax
      l=nterm-k
    2 bi(n)=ai(l)*c0(l,n)+z*bi(n)
    1 continue
c  Find integral of f(s+z)*g_m(s)ds from s1 to s2, where s1=x1-z, x2=x2-
c   =zeroth derivative of integral.
      xngi(1)=0.d0
      do 3 n=1,nterm
    3 xngi(1)=xngi(1)+bi(n)*(xgi2(n)-xgi1(n))
      if(ndriv.lt.2) return
c  Getfirst derivative of integral
c  Getzeroth, first, and higher derivatives of g_m
c  Need only ndriv-1 terms in derivative vector, since everything is
c  integrated once.
      ndrm1=ndriv-1
      call drivs(m,ndrm1,a,s1,di1)
      call drivs(m,ndrm1,a,s2,di2)
      f1=ai(nterm)
      f2=f1
      if(nterm.lt.2) go to 6
      do 5 n=2,nterm
      k=nterm-n+1
      f1=ai(k)+x1*f1
    5 f2=ai(k)+x2*f2
    6 xngi(2)=f1*di1(1)-f2*di2(1)
      if(nterm.lt.2) go to 7
c calculate coefficients in expansion of df/dx(x)=df/dx(s+z) in powers o
      ntrm1=nterm-1
      do 8 n=2,nterm
      nm1=n-1
      kmax=nterm-n
      bi(n)=ai(nterm)*c1(ntrm1,nm1)
      if(kmax.lt.1) go to 8
      do 9 k=1,kmax
      l=ntrm1-k
      lp=nterm-k
    9 bi(n)=ai(lp)*c1(l,nm1)+z*bi(n)
    8 continue
      do 10 n=2,nterm
      nm1=n-1
   10 xngi(2)=xngi(2)+bi(n)*(xgi2(nm1)-xgi1(nm1))
    7 if(ndriv.lt.3) return
c    Now get second derivative of integral
      xngi(3)=f2*di2(2)-f1*di1(2)
      if(nterm.lt.2) go to 11
c  Evaluate df/dx at x1 and x2
      f1p=c1(ntrm1,1)*ai(nterm)
      f2p=f1p
      if(nterm.lt.3) go to 12
      do 13 n=3,nterm
      k=nterm-n+2
      km1=k-1
      f1p=ai(k)*c1(km1,1)+x1*f1p
   13 f2p=ai(k)*c1(km1,1)+x2*f2p
   12 xngi(3)=xngi(3)+f1p*di1(1)-f2p*di2(1)
      if(nterm.lt.3) go to 11
c  Find coefficients in expansion of d**2 f(x)/dx**2  in powers of s.
      ntrm2=nterm-2
      do 14 n=3,nterm
      nm2=n-2
      kmax=nterm-n
      bi(n)=ai(nterm)*c2(ntrm2,nm2)
      if(kmax.lt.1) go to 14
      do 15 k=1,kmax
      l=ntrm2-k
      lp=nterm-k
   15 bi(n)=ai(lp)*c2(l,nm2)+z*bi(n)
   14 continue
      do 16 n=3,nterm
      nm2=n-2
   16 xngi(3)=xngi(3)+bi(n)*(xgi2(nm2)-xgi1(nm2))
   11 if(ndriv.lt.4) return
c  Get3rd derivative of integral.
      xngi(4)=f1*di1(3)-f2*di2(3)
      if(nterm.lt.2) go to 17
      xngi(4)=xngi(4)+f2p*di2(2)-f1p*di1(2)
      if(nterm.lt.3) go to 17
c  evaluate d**2 f/ dz**2 at x1 and x2
      f1pp=c2(ntrm2,1)*ai(nterm)
      f2pp=f1pp
      if(nterm.lt.4) go to 18
      do 19 n=4,nterm
      k=nterm-n+3
      km2=k-2
      f1pp=ai(k)*c2(km2,1)+x1*f1pp
   19 f2pp=ai(k)*c2(km2,1)+x2*f2pp
   18 xngi(4)=xngi(4)+f1pp*di1(1)-f2pp*di2(1)
      if(nterm.lt.4) go to 17
      delxg=xgi2(1)-xgi1(1)
      xngi(4)=xngi(4)+6.d0*ai(4)*delxg
      if(nterm.lt.5) go to 17
      xngi(4)=xngi(4)+24.d0*ai(5)*(z*delxg+xgi2(2)-xgi1(2))
   17 if(ndriv.lt.5) return
c  Integral of 4th derivative
      xngi(5)=f2*di2(4)-f1*di1(4)
      if(nterm.lt.2) go to 20
      xngi(5)=xngi(5)-f2p*di2(3)+f1p*di1(3)
      if(nterm.lt.3) go to 20
      xngi(5)=xngi(5)+f2pp*di2(2)-f1pp*di1(2)
      if(nterm.lt.4) go to 20
      f1ppp=6.d0*ai(4)
      f2ppp=f1ppp
      if(nterm.lt.5) go to 21
      f1ppp=f1ppp+24.d0*ai(5)*x1
      f2ppp=f2ppp+24.d0*ai(5)*x2
   21 xngi(5)=xngi(5)+f1ppp*di1(1)-f2ppp*di2(1)
      if(nterm.lt.5) go to 20
      xngi(5)=xngi(5)+24.d0*ai(5)*(xgi2(1)-xgi1(1))
   20 if(ndriv.lt.6) return
c  ndriv=6 or higher- ie., 5th or higher derivatives
      sign=-1.d0
      do 22 n=6,ndriv
      sign=-sign
      nm1=n-1
      xngi(n)=f1*di1(nm1)-f2*di2(nm1)
      if(nterm.lt.2) go to 22
      nm2=n-2
      xngi(n)=xngi(n)+f2p*di2(nm2)-f1p*di1(nm2)
      if(nterm.lt.3) go to 22
      nm3=n-3
      xngi(n)=xngi(n)+f1pp*di1(nm3)-f2pp*di2(nm3)
      if(nterm.lt.4) go to 22
      nm4=n-4
      xngi(n)=xngi(n)+f2ppp*di2(nm4)-f1ppp*di1(nm4)
      if(nterm.lt.5) go to 22
      nm5=n-5
      xngi(n)=xngi(n)+24.d0*ai(5)*(di1(nm5)-di2(nm5))
   22 xngi(n)=xngi(n)*sign
      return
 1001 write(6,300)  nterm
  300 format(1x,'nterm=',i4,' was >5 in XNGMIN --stopped')
      stop
 1002 write(6,301)   nterm
  301 format(1x,'nterm=',i4,' was < 1 in XNGMIN --stopped')
      stop
      end
c
c**************************************************
c
c
c  next group of called routines
c
      function cel(qqc,pp,aa,bb)
      implicit double precision(a-h,o-z)
c
c     returns the general complete elliptic integral cel(kc,p,a,b) with
c     qqc=kc, pp=p, aa=a and bb=b.
c
      parameter (ca=.00001d0  , pio2=1.570796326795d0)
c
c     the desired accuracy is the square of ca
c
      if(qqc.eq.0.d0) go to 21
   22 format(1x,'xkc=0 in cel-stopped')
      qc=dabs(qqc)
      a=aa
      b=bb
      p=pp
      e=qc
      em=1.d0
      if(p.gt.0.d0)then
      p=dsqrt(p)
      b=b/p
      else
      f=qc*qc
      q=1.d0-f
      g=1.d0-p
      f=f-p
      q=q*(b-a*p)
      p=dsqrt(f/g)
      a=(a-b)/g
      b=-q/(g*g*p)+a*p
      endif
    1 f=a
      a=a+b/p
      g=e/p
      b=b+f*g
      b=b+b
      p=g+p
      g=em
      em=qc+em
      if(dabs(g-qc).gt.g*ca)then
      qc=dsqrt(e)
      qc=qc+qc
      e=qc*em
      go to 1
      endif
      cel=pio2*(b+a*em)/(em*(em+p))
      return
   21 write(6,22)
      stop
      end
c
c**************************************************
c
      subroutine drivs(m,ndriv,a,s,di)
      implicit double precision(a-h,o-z)
c  This subroutine calculates repeated derivatives of g_m=
c
c    a**(1-m)d/da(a**m/(a**2+s**2)**(m+1/2))
c
c  di(1) = g_m, di(2)=d g_m /ds, di(3)= d**2 g_m /ds**2, etc.
c
c  Uses formula from Gradstein and Ryzhik, Table of Integrals,
c  Series, and Products, 1980, p.20
c
c  a=coil radius
c  m=multipole index (m=1 for dipole,2 for quadrupole, 3 for sextupole,e
c  s=x-z, where x=coil z coordinate, z= field point z coordinate.
c  ndriv= no. of entries in di (no. of times g_m is differentiated + 1)
c
c  maxdrv is the maximum value of ndriv required.
      parameter (maxdrv=20,maxhlf=11)
c  maxhlf must be at least 1/2 of maxdrv if maxdrv is even,
c  1/2(maxdrv+1) if maxdrv is odd.
c
      dimension di(maxdrv)
      dimension cnk(maxhlf,maxdrv),rnk1(maxhlf,maxdrv),
     &rnk2(maxhlf,maxdrv),snk(maxhlf,maxdrv),rkfac(maxhlf)
      aa=1.d0/a**2
      xm=dfloat(m)
      u=1.d0+aa*s**2
      ru=1.d0/u
      rru=dsqrt(ru)
      m2mm1=-2*m-1
      denom=ru**m*rru*a**m2mm1
      c1=xm*denom
      c2=(2.d0*xm+1.d0)*denom*ru
      di(1)=c1-c2
      if(ndriv.lt.2) return
      tas=2.d0*aa*s
      p1=-xm-0.5d0
      p2=p1-1.d0
      c1=c1*p1*ru
      c2=c2*p2*ru
      di(2)=(c1-c2)*tas
      if(ndriv.lt.3) return
      p1m=p1
      p2m=p2
      au=aa*u
      aunh=1.d0
      snk(1,1)=tas
c  n=order of derivative
c  find integers nhalf= nearest integer below or = (ndriv-1)/2,
c    and ileft=1+(ndriv-1)-2*nhalf
      xhalf=0.5d0*dfloat(ndriv-1)
      nhalf=xhalf
      dif=xhalf-dfloat(nhalf)
      ileft=1
      if(dif.gt.0.1d0) ileft=2
      n=1
      rkf=1.d0
      do 1 nh=1,nhalf
      rkf=rkf/dfloat(nh)
      rkfac(nh)=rkf
      nhp1=nh+1
      ilef=2
      if(nh.lt.nhalf) go to 6
      if(ileft.eq.1) ilef=1
    6 aunh=aunh*au
      nhp1=nh+1
      do 1 il=1,ilef
      nm1=n
      nm2=n-1
      n=n+1
      p1m=p1m-1.d0
      p2m=p2m-1.d0
      c1=c1*p1m*ru
      c2=c2*p2m*ru
      snk(nhp1,n)=aunh
      do 4 k=1,nh
    4 snk(k,n)=snk(k,nm1)*tas
      if(il.eq.2) snk(nhp1,n)=snk(nhp1,n)*tas
      cnk(2,n)=dfloat(n*(n-1))
      rnk1(2,n)=1.d0/p1m
      rnk2(2,n)=1.d0/p2m
      if(nhp1.lt.3) go to 8
      do 7 k=3,nhp1
      km1=k-1
      cnk(k,n)=cnk(km1,nm2)*cnk(2,n)
      rnk1(k,n)=rnk1(km1,nm1)*rnk1(2,n)
    7 rnk2(k,n)=rnk2(km1,nm1)*rnk2(2,n)
    8 sum1=snk(1,n)
      sum2=sum1
      do 5 k=2,nhp1
      km1=k-1
      sum1=sum1+snk(k,n)*cnk(k,n)*rnk1(k,n)*rkfac(km1)
    5 sum2=sum2+snk(k,n)*cnk(k,n)*rnk2(k,n)*rkfac(km1)
      np1=n+1
    1 di(np1)=sum1*c1-sum2*c2
      return
      end
c
c**************************************************
c
      function el2(x,qqc,aa,bb)
      implicit double precision(a-h,o-z)
c
c     returns the general elliptic integral of the second kind,
c     el2(x,kc,a,b) with x.ge.0,qqc=kc,aa=a and bb=b.
c
      parameter (pi=3.14159265359d0  , ca=.00001d0,cb=1.d-11)
c
c     the desired accuracy is the square of ca, while cb should be
c     set to 0.01 times desired accuracy.
c
      if(x.eq.0.d0)then
      el2=0.d0
      else if(qqc.ne.0.d0) then
      qc=qqc
      a=aa
      b=bb
      c=x**2
      d=1.d0+c
      p=dsqrt((1.d0+qc**2*c)/d)
      d=x/d
      c=d/(2.d0*p)
      z=a-b
      eye=a
      a=0.5d0*(b+a)
      y=abs(1.d0/x)
      f=0.d0
      l=0
      em=1.d0
      qc=dabs(qc)
    1 b=eye*qc+b
      e=em*qc
      g=e/p
      d=f*g+d
      f=c
      eye=a
      p=g+p
      c=0.5d0*(d/p+c)
      g=em
      em=qc+em
      a=0.5d0*(b/em+a)
      y=-e/y+y
      if(y.eq.0)y=dsqrt(e)*cb
      if(dabs(g-qc).gt.ca*g)then
      qc=dsqrt(e)*2.d0
      l=l+l
      if(y.lt.0.) l=l+1
      go to 1
      endif
      if(y.lt.0.)l=l+1
      e=(datan(em/y)+pi*l)*a/em
      if(x.lt.0) e=-e
      el2=e+c*z
      else
      write(6,21)
   21 format(1x,'xkc=0 in el2-stopped')
      stop
c     argument qqc was zero
      endif
      return
      end
c
c**************************************************
c
c
c
c     this is a function subprogram copied from numerical recipes
c     by w.h.press ,etl pp. 186-187.
c     evaluates the general complete elliptic integral
c     as the following function of four variables:
c
c     cel(kc,p,a,b)=
c     int[ (a+b*x**2)dx/{1+p*x**2)*sqrt((1+x**2)*(1+kc**2*x**2))}
c     from x=0 to x=infinite
c
      subroutine intend(a1,a2,z1,z2,m,xiend)
c   Used for m even case, 0th derivative of on-axis gradient of a PMM.
      implicit double precision(a-h,o-z)
c  Evaluates the double integral from a1 to a2 and z1 to z2 of
c
c     drho*dz/((rho**2+z**2)**(m/2+1/2))
c
c    using a polar coordinate approach.
c    z1 or z2 can be zero; a1 and a2 are assumed to always be greater
c   than zero.
c
c    mMUST BE EVEN!
c
c    mgreater than or equal to 2
c
      parameter(mmax=26)
      dimension amk(mmax)
c     data small/1.d-1/
      xm=dfloat(m)
      xmm1=xm-1.d0
      mm1=m-1
      mhalf=m/2
      a12=a1**2
      a22=a2**2
      z12=z1**2
      z22=z2**2
      hyp12=a12+z22
      hyp22=a22+z22
      hyp32=a22+z12
      hyp42=a12+z12
      hyp1=dsqrt(hyp12)
      hyp2=dsqrt(hyp22)
      hyp3=dsqrt(hyp32)
      hyp4=dsqrt(hyp42)
      sin1=a1/hyp1
      sin2=a2/hyp2
      sin3=a2/hyp3
      sin4=a1/hyp4
      cos1=z2/hyp1
      cos2=z2/hyp2
      cos3=z1/hyp3
      cos4=z1/hyp4
c    Find coefficients for cosine**(m-1) and sine**(m-1) integrals--
c   they are the same for both--see Gradshtein&Ryzhik, p. 131.
c  There are a total of m/2-1 coefficients. 2l+1 in G&R is m-1 here.
c  Skip calculation of coeffficients if m=2
      if(m.lt.3) go to 2
      xnum=xm-2.d0
      dnom=xm-3.d0
      amk(1)=xnum/dnom
      mhm1=mhalf-1
      if(mhm1.lt.2) go to 23
      do 22 j=2,mhm1
      jm1=j-1
      xnum=xnum-2
      dnom=dnom-2
   22 amk(j)=amk(jm1)*xnum/dnom
   23 continue
c   m>2 case-- use subroutine INTCOS
c   RHS vertical line
      call intcos(cos1,cos2,sin1,sin2,hyp12,hyp22,z2,amk,m,cosint)
      xiend=cosint
c  LHSvertical line
      call intcos(cos3,cos4,sin3,sin4,hyp32,hyp42,z1,amk,m,cosint)
      xiend=xiend+cosint
c  Tophorizontal line
      call intcos(sin1,sin4,cos1,cos4,hyp12,hyp42,a1,amk,m,cosint)
      xiend=xiend+cosint
c  Bottom horizontal line
      call intcos(sin2,sin3,cos2,cos3,hyp22,hyp32,a2,amk,m,cosint)
      xiend=xiend-cosint
      xiend=xiend/dfloat(m-1)
      return
c  m=2case
    2 xiend=0.d0
      if(cos1.eq.0.d0) go to 3
      if(dabs(cos1).gt.0.1d0) go to 4
c  small angle formula used when z2 is small
c   Floating-point numbers in the expression below are the coefficients
c  theexpansion for (1+x)**(-1/2)
      e12=z22/a12
      e22=z22/a22
      t2=(0.5d0-e22*(0.375d0-e22*(0.3125d0-e22*(0.2734375d0-
     &e22*(0.24609375d0-e22*(0.2255859375d0-
     &e22*(0.20947265625d0-e22*(0.196380615234375d0-
     &e22*(0.1854705810546875d0-
     &e22*0.1761970520019531d0)))))))))/a22
      t1=(0.5d0-e12*(0.375d0-e12*(0.3125d0-e12*(0.2734375d0-
     &e12*(0.24609375d0-e12*(0.2255859375d0-
     &e12*(0.20947265625d0-e12*(0.196380615234375d0-
     &e12*(0.1854705810546875d0-
     &e12*0.1761970520019531d0)))))))))/a12
      xiend=z2*(t2-t1)
      go to 3
c   General (z2 not small) case
    4 xiend=(sin1-sin2)/z2
c  Left hand vertical line at z1
    3 if(cos4.eq.0.d0) go to 5
      if(dabs(cos4).gt.0.1d0) go to 6
c  Small angle formula for small z1
      e12=z12/a12
      e22=z12/a22
      t2=(0.5d0-e22*(0.375d0-e22*(0.3125d0-e22*(0.2734375d0-
     &e22*(0.24609375d0-e22*(0.2255859375d0-
     &e22*(0.20947265625d0-e22*(0.196380615234375d0-
     &e22*(0.1854705810546875d0-
     &e22*0.1761970520019531d0)))))))))/a22
      t1=(0.5d0-e12*(0.375d0-e12*(0.3125d0-e12*(0.2734375d0-
     &e12*(0.24609375d0-e12*(0.2255859375d0-
     &e12*(0.20947265625d0-e12*(0.196380615234375d0-
     &e12*(0.1854705810546875d0-
     &e12*0.1761970520019531d0)))))))))/a12
      xiend=xiend-z1*(t2-t1)
      go to 5
c   General (z1 not small) case
    6 xiend=xiend+(sin3-sin4)/z1
c   Upper and lower horizontal lines
    5 xiend=xiend+(cos1-cos4)/a1+(cos3-cos2)/a2
      return
      end
c
c**************************************************
c
      subroutine rhoint(a1,a2,s,is,m,nd,rhoin)
c  Checked by numerical integration
c  Used for gtype2 (thick Halbach)
      implicit double precision(a-h,o-z)
c  This subroutine evaluates the components of the vector rhoin, which
c  arethe integrals
c   integral from a1 to a2  dr of
c   (r**(m+2))/((r**2+s**2)*(m+k+3/2))  ,
c   k=0,1,2,3,...nd-1
c
c   m=multipole index.  m=1 for dipole, 2 for quad.,etc. m is > or = 1.
c   ndis number of derivatives. This subroutine called only if nd > or
c
      parameter(maxdrv=20,small=1.d-6,smal=5.d-2)
      parameter(nterms=15)
      parameter(maxcof=35)
      common /dnms/ denom(maxcof,2,2)
c  NTERMS is the number of terms used in the small-s expansion. Can be m
c  large, while SMAL is also made larger, for check.
      dimension rhoin(maxdrv),rhon1(maxdrv)
      dimension rhon2(maxdrv),xirend(maxdrv)
c
c  Calls subroutine RHOEND if m is even (2 or greater).
c
c  First check to see if s is small or zero. If so, use expansion in s/a
      x1=(s/a1)**2
      if(x1.lt.smal) go to 3
c
c  Next check to see if m is odd or even
c
      a12=a1**2
      a22=a2**2
      x=dfloat(m+3)*0.5d0+small
      iterms=x
      dif=x-dfloat(iterms)
      ileft=0
      if(dif.gt.0.1d0) ileft=1
c  ileft=0 if m is odd, =1 if m is even
      if(ileft.gt.0) go to 1
c   m is odd
      do 2 l=1,nd
    2 rhoin(l)=0.d0
c  Exponent on denominator in innermost term=n+k+1/2 ,k=0,1,..nd-1
      nmin=(m-1)/2
      a1prod=1.d0
      a2prod=1.d0
      itop=0
      ibot=m-2
      go to 5
    1 continue
c     m is even- descending series ends in an integral evaluated by RHOE
      nmin=m/2
      call rhoend(a1,a2,s,is,m,nd,xirend)
      do 6 l=1,nd
    6 rhoin(l)=-xirend(l)
      a1prod=a1
      a2prod=a2
      itop=1
      ibot=m-1
    5 itrm1=iterms-1
      np1=nmin+1
      do 8 i=1,itrm1
      itop=itop+2
      ibot=ibot+2
      do 10 l=1,nd
      k=l-1
      np1pk=np1+k
   10 rhoin(l)=(rhoin(l)+a2prod*denom(np1pk,2,is)-
     &a1prod*denom(np1pk,1,is))*dfloat(itop)/dfloat(ibot+2*k)
      a1prod=a1prod*a12
      a2prod=a2prod*a22
      np1=np1+1
    8 continue
      ibot=ibot+2
      do 7 l=1,nd
      k=l-1
      np1pk=np1+k
    7 rhoin(l)=-(rhoin(l)+a2prod*denom(np1pk,2,is)-
     &a1prod*denom(np1pk,1,is))/dfloat(ibot+2*k)
      return
    3 x2=(s/a2)**2
c  Small-s expansion algorithm
c  No.of terms should be increased if quadruple precision, etc. is used
      do 11 l=1,nd
      k=l-1
      rhon1(l)=1.d0/dfloat(m+2*k)
      rhon2(l)=rhon1(l)
      xnum=dfloat(m+k+1)-0.5d0
      fac=1.d0
      y1=1.d0
      y2=1.d0
      do 9 n=1,nterms
      xnum=xnum+1.d0
      y1=-y1*xnum*x1
      y2=-y2*xnum*x2
      fac=fac*dfloat(n)
      dnm=1.d0/(fac*dfloat(m+2*(k+n)))
      rhon1(l)=rhon1(l)+y1*dnm
      rhon2(l)=rhon2(l)+y2*dnm
    9 continue
      n=m+2*(l-1)
   11 rhoin(l)=rhon1(l)*a1**(-n)-rhon2(l)*a2**(-n)
      return
      end
c
c**************************************************
c
      subroutine xngint(m,nint,a,s,xgi)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      parameter(c83=8.d0/3.d0)
      parameter(t3rds=2.d0/3.d0)
      parameter(sv3rds=7.d0/3.d0)
      common /bicof/ bcoeff(maxcof,maxcof)
c  This subroutine calculates indefinite integrals to s of gm*x**(n-1),
c   n=1,nint, fixed m.
c
c  where gm=a**(1-m)*d/da((a**m)/(a**2+x**2)**(m+1/2))
c
c  maximum nint=5
c    "m  =9
c
c   m=multipole index
c   nint=no. of integrals, = no. of  components in xgi
c   xgi=vector of nint integrals
c   a=coil radius
      dimension xgi(5)
c   bcoeff contains the binomial expansion coefficients up to order maxc
c  input parameter checks
      if(m.gt.100) go to 1001
      if(m.lt.1) go to 1002
      if(nint.gt.5) go to 1003
      if(a.le.0.d0) go to 1004
c  m=1and m=2 are treated specially
      if(m.gt.1) go to 20
c  m=1case
      a2=a**2
      ra2=1.d0/a2
      s2=s**2
      u2=a2+s2
      ru2=1.d0/u2
      u=dsqrt(u2)
      ru=1.d0/u
      xgi(1)=-s*(ra2+ru2)*ru
      if(nint.lt.2) return
      xgi(2)=-ru+a2*ru*ru2
      if(nint.lt.3) return
      dl=dlog(s+u)
      xgi(3)=-s*ru+dl+a2*s*ru2*ru+a2*ru/(s+u)
      if(nint.lt.4) return
      xgi(4)=u+a2*ru*(4.d0-a2*ru2)
      if(nint.lt.5) return
      xgi(5)=s*(0.5d0*u+ru*a2*(4.d0+s2*ru2))-4.5d0*a2*dl
      return
c  m=2case
   20 if(m.gt.2) go to 30
      a2=a**2
      a4=a2**2
      s2=s**2
      u2=a2+s2
      u22=u2**2
      u=dsqrt(u2)
      ru=1.d0/u
      ru2=1.d0/u2
      xgi(1)=(s*ru*(c83*s2*ru2-3.d0-(s2*ru2)**2))/a4
      if(nint.lt.2) return
      xgi(2)=ru2*ru*(a2*ru2-t3rds)
      if(nint.lt.3) return
      xgi(3)=-s2*s/(u22*u)
      if(nint.lt.4) return
      xgi(4)=ru*(a2*ru2*(sv3rds-a2*ru2)-2.d0)
      if(nint.lt.5) return
      xgi(5)=2.d0*dlog(s+u)-
     &s*ru*(2.d0+s2*ru2*(t3rds+ru2*s2))
      return
   30 continue
c  m=3or greater
      xm=dfloat(m)
      mm1=m-1
      mm2=m-2
      mm3=m-3
      mp1=m+1
      xmm1=dfloat(mm1)
      xmm2=dfloat(mm2)
      a2=a**2
      a2m=a2**m
      a2mm2=a2m/a2
      a2mm4=a2mm2/a2
      x2mp1=dfloat(2*m+1)
      tmm1=dfloat(2*m-1)
      tmm3=dfloat(2*m-3)
      r1=xm/tmm3
      r2=dfloat(3*m+1)/tmm1
      z=s
      zsq=z**2
      u2=zsq+a2
      ru2=1.d0/u2
      ur2=a2/u2
      u=dsqrt(u2)
      rz2u2=zsq/u2
      rzu=z/u
      rend=rz2u2**m*rzu
      u2mm3=u2**mm2*u
      u2mm1=u2*u2mm3
c  n=1
      sign=-1.d0
      sum=0.d0
      do 2 k=1,m
      sign=-sign
      k2m1=2*k-1
      km1=k-1
      tkm1=dfloat(k2m1)
      if(km1.eq.0) go to 21
      term=rzu/tkm1*rz2u2**km1
      go to 22
   21 term=rzu/tkm1
   22 c=xm*bcoeff(k,mm1)-x2mp1*bcoeff(k,m)
    2 sum=sum+c*term*sign
      sum=sum+sign*rend
      sum=sum/a2m
      xgi(1)=sum
      if(nint.lt.2) return
c  n=2
      term=(ur2-xm/tmm1)/u2mm1
      xgi(2)=term
      if(nint.lt.3) return
c  n=3
      sign=-1.d0
      sum=0.d0
      do 3 k=1,mm1
      km1=k-1
      sign=-sign
      tkp1=dfloat(2*k+1)
      term=rzu/tkp1*rz2u2**k
      c=xm*bcoeff(k,mm2)-x2mp1*bcoeff(k,mm1)
    3 sum=sum+c*term*sign
      sum=sum+sign*rend
      sum=sum/a2mm2
      xgi(3)=sum
      if(nint.lt.4) return
c  n=4
      term=(r1-ur2*(r2-ur2))/u2mm3
      xgi(4)=-term
      if(nint.lt.5) return
c  n=5
      if(m.gt.3) go to 8
c  m=3is a special case
      term=(rend-0.8d0*rzu*rz2u2**2)/a2
      xgi(5)=term
      return
    8 sign=-1.d0
      sum=0.d0
      do 4 k=1,mm2
      tkp3=dfloat(2*k+3)
      kp1=k+1
      sign=-sign
      c=xm*bcoeff(k,mm3)-x2mp1*bcoeff(k,mm2)
      term=rzu/tkp3*rz2u2**kp1
    4 sum=sum+c*term*sign
      sum=sum+sign*rend
      sum=sum/a2mm4
      xgi(5)=sum
      return
 1001 write(6,3001)
 3001 format(1x,'m greater than 100 in xngint-stopped')
      stop
 1002 write(6,3002)
 3002 format(1x,'m less than 1 in xngint-stopped')
      stop
 1003 write(6,3003)
 3003 format(1x,'nint greater than 5 in xngint-stopped')
      stop
 1004 write(6,3004)
 3004 format(1x,'a less than or equal 0 in xngint-stopped')
      stop
      end
c
c**************************************************
c
      subroutine zint(a1,a2,s1,s2,nmin,m,zinti)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
c   The array ZINTI contains the integrals from s1 to s2 of
c
c    ds*a**2n/(a**2+s**2)**n+1/2, n=nmin,nmin+1,nmin+2,...m
c
c   zinti(i,1) is evaluated at a1, zinti(i,2) at a2.
c   The index i runs from 1 to m-nmin+1, with n=nmin for i=1, nmin+1 for
c   i=2, etc.
c
c
c   Uses analytical formula from Gradshtein and Ryzhik, p. 86.
c
c
      common /bicof/ bcoeff(maxcof,maxcof)
c  bcoeff(k,n) is an element of array containing the binomial coefficien
c  index k goes from 1 to maxcof; elements with k>n+1 are zero.
c  index n goes from 1 to maxcof
      common/dnms/denom(maxcof,2,2)
      dimension zinti(maxcof,2)
      if(nmin.gt.0) go to 1
c  nmin=0 --- First term is a logarithmic expression, since m=1
c  There are two terms.
      a12=a1**2
      a22=a2**2
      s12=s1**2
      s22=s2**2
c  fora1:
      zinti(1,1)=dlog((s2+dsqrt(s22+a12))/(s1+dsqrt(s12+a12)))
c  fora2:
      zinti(1,2)=dlog((s2+dsqrt(s22+a22))/(s1+dsqrt(s12+a22)))
c  fora1:
      zinti(2,1)=s2*denom(1,1,2)-s1*denom(1,1,1)
c  fora2:
      zinti(2,2)=s2*denom(1,2,2)-s1*denom(1,2,1)
      return
c  nmin=1 or greater
    1 nterm=m-nmin+1
c  nterm is the number of elements in each of zinti(n,1) and zinti(n,2),
c  separately.
      s12=s1**2
      s22=s2**2
      s1n=s1
      s2n=s2
      do 2 l=1,nterm
c  rho=a1
      zinti(l,1)=s2*denom(1,1,2)-s1*denom(1,1,1)
c  rho=a2
    2 zinti(l,2)=s2*denom(1,2,2)-s1*denom(1,2,1)
      do 3 k=2,m
      r2kp1=1.d0/dfloat(2*k-1)
      s1n=-s12*s1n
      s2n=-s22*s2n
      ck1=s1n*r2kp1
      ck2=s2n*r2kp1
      lmin=1
      if(k.gt.nmin) lmin=k-nmin+1
      do 3 l=lmin,nterm
      n=nmin+l-2
      ckl1=bcoeff(k,n)*ck1
      ckl2=bcoeff(k,n)*ck2
c  rho=a1
      zinti(l,1)=zinti(l,1)+ckl2*denom(k,1,2)-ckl1*denom(k,1,1)
c  rho=a2
    3 zinti(l,2)=zinti(l,2)+ckl2*denom(k,2,2)-ckl1*denom(k,2,1)
      return
      end
c
c**************************************************
c
      subroutine intcos(cos1,cos2,sin1,sin2,r12,r22,x,amk,m,cosint)
      implicit double precision(a-h,o-z)
      dimension amk(m)
c  This subroutine evaluates 1/x times the integral from phi1 to phi2 of
c
c    (cos phi)**(m-1) * dphi
c
c  used only for even m values
c
c  when x is small, by definition cos1 and cos2 are small, and the
c  result is well-behaved as x goes to zero. For small x, a small
c  x expansion is used.  Otherwise, the expression from Gradstein and
c  Ryzhik, p. 131, is used to evaluate the integral, and the result is
c  divided by x.
c   analogous sine integral obtained by switching sine for cosines in
c   call and changing the sign of COSINT.
      if(cos1.eq.0.d0) go to 1
      xm=dfloat(m)
      mhalf=m/2
      if(dabs(cos1).gt.0.1d0) go to 2
c  Usesmall angle Taylor series to evaluate integral
      cos12=cos1**2
      cos22=cos2**2
      t2=1.d0/xm+cos22*(0.5d0/(xm+2.d0)+cos22*(0.375d0/(xm+4.d0)+
     &cos22*(0.3125d0/(xm+6.d0)+cos22*(0.2734375d0/(xm+8.d0)+
     &cos22*(0.24609375d0/(xm+1.d1)+cos22*(0.2255859375d0/
     &(xm+1.2d1)+cos22*(0.20947265625d0/(xm+1.4d1)+
     &cos22*(0.196380615234375d0/(xm+1.6d1)+
     &cos22*(0.1854705810546875d0/(xm+1.8d1)+
     &cos22*0.1761970520019531d0/(xm+2.d1))))))))))
      t1=1.d0/xm+cos12*(0.5d0/(xm+2.d0)+cos12*(0.375d0/(xm+4.d0)+
     &cos12*(0.3125d0/(xm+6.d0)+cos12*(0.2734375d0/(xm+8.d0)+
     &cos12*(0.24609375d0/(xm+1.d1)+cos12*(0.2255859375d0/
     &(xm+1.2d1)+cos12*(0.20947265625d0/(xm+1.4d1)+
     &cos12*(0.196380615234375d0/(xm+1.6d1)+
     &cos12*(0.1854705810546875d0/(xm+1.8d1)+
     &cos12*0.1761970520019531d0/(xm+2.d1))))))))))
      r1m=r12**(-mhalf)
      r2m=r22**(-mhalf)
      cosint=x*(t2*r2m-t1*r1m)
      return
c   General expression for cosine integral
    2 cosi2=cos2**2
      sini=sin2
      mhm1=mhalf-1
      mm1=m-1
      xmm1=dfloat(mm1)
      temp=0.d0
      do 11 i=1,2
      cos2k=1.d0
      sum=0.d0
      do 12 j=1,mhm1
      k=mhm1-j+1
      sum=sum+amk(k)*cos2k
   12 cos2k=cos2k*cosi2
      sum=sum+cos2k
      sum=sum*sini
      if(i.eq.1) sum=-sum
      temp=temp+sum
      cosi2=cos1**2
   11 sini=sin1
      cosint=temp/(xmm1*x**mm1)
      return
    1 cosint=0.d0
      return
      end
c
c**************************************************
c
      subroutine rhoend(a1,a2,s,is,m,nd,xirend)
c  Checked by numerical integration
c  used for gtype2 (thick Halbach)
      implicit double precision(a-h,o-z)
      parameter(maxcof=35,maxdrv=20)
      common /dnms/ denom(maxcof,2,2)
      common /bicof/ bcoeff(maxcof,maxcof)
      dimension xirend(maxdrv)
c  This subroutine evaluates the integral from a1 to a2 of
c
c  drho/(rho**2+s**2)**(m/2+k+1/2) , k=0,2,3....nd-1
c
c  s is assumed not to be small or zero-- if so, the subroutine calling
c  one(i.e. RHOINT) uses a series expansion in s, and therefore does no
c  call this one.
c
c  THIS SUBROUTINE IS CALLED ONLY IF m IS EVEN (2 or GREATER).
c  Also nd must be 1 or greater.
c
c   Uses analytical formula from Gradshtein and Ryzhik, p. 86.
c
c  is =1 for s1, 2 for s2
c  xirend is a vector of length nd containing the inetgrals
      a12=a1**2
      a22=a2**2
      a1n=a1
      a2n=a2
      rs2=1.d0/s**2
      mo2=m/2
      kmax=mo2+nd-1
      do 2 l=1,nd
    2 xirend(l)=a2*denom(1,2,is)-a1*denom(1,1,is)
      do 3 k=2,kmax
      r2kp1=1.d0/dfloat(2*k-1)
      a1n=-a12*a1n
      a2n=-a22*a2n
      ck1=a1n*r2kp1
      ck2=a2n*r2kp1
      lmin=1
      if(k.gt.mo2) lmin=k-mo2+1
      do 3 l=lmin,nd
      n=mo2+l-2
      ckl1=bcoeff(k,n)*ck1
      ckl2=bcoeff(k,n)*ck2
    3 xirend(l)=xirend(l)+ckl2*denom(k,2,is)-ckl1*denom(k,1,is)
      if(m.gt.2) go to 6
      sfac=rs2
      go to 7
    6 sfac=rs2**mo2
    7 do 5 l=1,nd
      xirend(l)=xirend(l)*sfac
      if(l.eq.nd) go to 5
      sfac=sfac*rs2
    5 continue
      return
      end
c
c**************************************************
c
