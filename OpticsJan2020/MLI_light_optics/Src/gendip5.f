************************************************************************
* header              GENDIP (GENMAP for a "parallel" face dipole      *
*                     magnet with soft fringe fields                   *
*  All routines needed for this special GENMAP                         *
************************************************************************
c
      subroutine bderivs(z,y,b)
c  This routine computes b(z) = Int(-Inf to z) By(z') dz
c  and the second and fourth derivative b2(z) and b4(z) for
c  (a dipole magnet.)
      use beamdata
      include 'impli.inc'
      include 'dip.inc'
c----------------------------------------
      double precision z, y(*), b(0:6,0:6)
c----------------------------------------
      external amyf,f1,f2,f4,f6
c----------------------------------------
      epsz2 = gap/sl
      bb = (By*sl)
      do 10 i=0,6
        do 10 j=0,6
          b(i,j) = 0.0d0
 10   continue
c
c      b = (amyf(z-za,epsz2) - amyf(z-zb,epsz2)) * bb
      b(0,0)= (f1(z-za,epsz2) -f1(z-zb,epsz2)) * bb
      b(0,1)= (f2(z-za,epsz2) -f2(z-zb,epsz2)) * bb
      b(0,3)= (f4(z-za,epsz2) -f4(z-zb,epsz2)) * bb
c      b6= (f6(z-za,epsz2)-f6(z-zb,epsz2)) * bb
      return
      end
c
      function amyf(y,eps)
      double precision amyf,y,eps
      if (y/eps.lt.0.d0) then
        amyf = eps*Log(1 + Exp(2*y/eps))/2
      else
        amyf = y + eps*Log(Exp(-2*y/eps)+1)/2
      endif
      return
      end
c
      function f1(y,eps)
      double precision f1,y,eps
        f1 = tanh(y/eps)/2.d0 +.5d0
      return
      end
c
      function f2(y,eps)
      double precision f2,y,eps
      if (dabs(y/eps).lt.10.d0) then
        f2 = 2*Exp(2*y/eps)/(eps*(1 + Exp(2*y/eps))**2)
      else
        f2 = 0.d0
      endif
      return
      end
c
      function f4(y,eps)
      double precision f4,y,eps
      if (dabs(y/eps).lt.10.d0) then
        f4 = 8*Exp(2*y/eps)*
     &      (1 - 4*Exp(2*y/eps) + Exp(4*y/eps))/
     &   (eps**3*(1 + Exp(2*y/eps))**4)
      else
        f4 = 0.0d0
      endif
      return
      end
c
       function f6(y,eps)
       double precision f6,y,eps
      if (dabs(y/eps).lt.10.d0) then
         f6 = 32*Exp(2*y/eps)*
     &        (1 - 26*Exp(2*y/eps) + 
     &            66*Exp(4*y/eps) - 
     &            26*Exp(6*y/eps) + Exp(8*y/eps))/
     &   (eps**5*(1 + Exp(2*y/eps))**6)
      else
        f6 = 0.d0
      endif
      return
      end
c
************************************************************************
      subroutine gendip(p,fa,fm)
c This is a subroutine for computing the map for a soft edged dipole
c magnet ( F. Neri 5/16/89 ).
c The routine is based on Rob Ryne original gendip, but all the code has been
c rewritten.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'files.inc'
      include 'dip.inc'
      include 'pie.inc'
c
c  calling arrays
      dimension p(6)
      dimension fa(monoms), fm(6,6)
c
c  local arrays
      dimension pb(6)
      dimension y(monoms+15)
c
      real ttaa, ttbb
c
c use equivalence statement to make the various parameter sets pstj
c available as if they were in a two dimensional array
c
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c
c  get interval and number of steps from GENREC parameters
c
      za = p(1)
      zb = p(2)
      ns = 10000
      By = p(3)
      phi = p(4) * pi180
      s = sin(phi)
      gap = p(5)
      zi  = 0.0d0
      zf  = p(6)
c
      write(6,*) 'za=',za,'zb=',zb
      write(6,*) 'By=',By,'phi=',p(4)
      write(6,*) 'gap=',gap,'zf=',zf
c
      jtty = 1
      jdsk = 1
c
c      if((isend.eq.1).or.(isend.eq.3)) jtty = 1
c      if((isend.eq.2).or.(isend.eq.3)) jdsk = 1
c
      h=(zf-zi)/float(ns)
c
c     call VAX system routine for timing report
c
      ttaa = secnds(0.0)
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=0.d0
      y(2)= s
      y(3)=0.d0
      y(4)=0.d0
      y(5)=0.d0
      y(6)=-1.d0/beta
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c
c  initialize map to the identity map:
      ne=224
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c  do the computation:
      t=zi
      iflag = 3
      call adam11(h,ns,'start',t,y)
      call errchk(y(7))
      call putmap(y,fa,fm)
      s1 = -y(2)
      phi1 = dasin(s1)
      phi1deg = phi1/pi180
      write(jof , 991) phi1deg
      write(jodf, 991) phi1deg
991   format('  Final angle is ',f16.8,' degrees')
c
c     call VAX system routine for timing report
c
      ttbb = secnds(ttaa)
      if(jtty.eq.1) write( jof,567) ttbb
      if(jdsk.eq.1) write(jodf,567) ttbb
 567  format(' GENDIP integration time = ',f12.2,' sec.')
      end
c
**********************************************************************
c
      subroutine hmltn3(t,y,h)
c  this routine is used to specify h(z) for a dipole magnet.
c  Written by F. Neri, 5/16/89.
c  Modified 6/3/89 to handle different gauges.
c  The design orbit is assumed to be in the Y = 0 plane.
c  B Field derivatives on the design orbit are provided by the routine
c  BDERIVS as a function of z,and x = y(1), bderivs is called as
c      call bderivs(z,y,b)
c  The derivatives are stored in the array b(0:*,0:*), defined in the 
c  include file bfield.inc, and passed to other parts of the program 
c  in the common/bfield/b(0:*,0:*)
c  The array b is arranged so that B(0,0) is By, B(1,0) is d(By)/(dX),
c  b(1,1) is d2(By)/(dX dz), etc.
c
      use beamdata
      use lieaparam, only : monoms
      parameter (itop=209,iplus=224)
      include 'impli.inc'
      include 'dip.inc'
      include 'bfield.inc'
c
      dimension h(monoms),y(*)
c
      dimension X(0:itop),YY(0:itop),P1(0:itop),P2(0:itop)
      dimension A(0:12)
c  begin calculation
c
c  compute gradients
      call bderivs(t,y,b)
c
c  initialization
      do 10 i=1,monoms
   10 h(i)=0.d0
c
c cccc
c   Slow method to produce hamiltonian: use polynomial
c   expansion of square root: when this method works we will use it to
c   to check the "hardwired" version.
c cccc
c   Coefficients of expansion of -SQRT(1+X)+1
      A(0) = 0.d0
      A(1) = -1.d0/2.d0
      do 50 i=2, 6
        A(i) =  A(i-1) * (1.d0/2.d0 - (i-1.d0))/i
   50 continue
c
      do 60 i=0,209
        X(i) = 0.d0
   60 continue
      X(6) = -2.d0 / beta
c      X(13) = -1.d0
      X(22) = -1.d0
      X(27) =  1.d0
c
      maxord = 4
      nn = 4
      cos2 = 1 - y(2)**2
c P1 = px - q Ax
      do 101 i=0,itop
  101 P1(i) = 0.d0
c Note constant term in px - q Ax, coming from design orbit 
c y(2) = Px = sin(phi).
      P1(0) = y(2)
c
      P1(2) = 1.d0
      P1(18) = b(0,1)/(2.d0*brho)
      P1(39) = b(1,1)/(2.d0*brho)
      P1(95) = b(2,1)/(4.d0*brho)
      P1(175) = -(b(2,1)+b(0,3))/(24.d0*brho)
c P2 = P1**2
      call mypmult(P1,P1,P2,maxord)
c Zero order term subtracted ( Sum starts at 1 ):
c Divide by cos2
      do 102 i=1,itop
  102 X(i) = (X(i) - P2(i))/cos2
c X = pt**2 - 2/beta*pt - (px -q Ax)**2 - (py)**2
c YY = -Sqrt(1+X)
      call mypoly1(nn,A,X,YY,maxord)
c h = -Az
      h(7) = b(1,0)/(2.d0*brho)
      h(18) = -b(1,0)/(2.d0*brho)
      h(28) = b(2,0)/(6.d0*brho)
      h(39) = -b(2,0)/(6.d0*brho)
      h(84) = b(3,0)/(24.d0*brho)
      h(95) = -b(3,0)/(4.d0*brho)
      h(175) = (b(3,0)+b(1,2))/(24.d0*brho)
c h = -Sqrt(1+X) - Az
c Zero and first order terms subtracted ( Sum starts at 7 ):
c Scale by cos
      do 70 i=7, 209
        h(i) = h(i) + dsqrt(cos2)*YY(i)/sl
   70 continue
c
      return
      end
c
c ******************************************************************
c
c  Aux polynomial routines (really a poor man's DA package).
c  They don't really belong here.
c
c ******************************************************************
c
      subroutine mypoly1(N,A,X,Y,maxord)
      parameter (MN=6)
      include 'lims.inc'
      double precision X(0:209),Y(0:209)
      double precision A(0:N)
c
      double precision Vect(0:209,MN)
c
c   NO COMMENT
      do 100 i=0,top(maxord)
 100    Y(i) = 0.0d0
      do 200 i=0,top(maxord)
 200    Vect(i,1) = X(i)
      do 300 iord=2,N
        call mypmult(X,Vect(0,iord-1),Vect(0,iord),maxord)
 300  continue
      Y(0) = A(0)
      do 400 iord=1,N
        do 500 i=0,top(maxord)
          Y(i) = Y(i) + Vect(i,iord)*A(iord)
 500    continue
 400  continue
      return
      end
c
      subroutine mypmult(p1,p2,p3,maxord)
      double precision p1(0:209),p2(0:209),p3(0:209)
c
      include 'lims.inc'
      if (maxord.lt.0) return
c
      do 10 i=1,top(maxord)
   10 p3(i) = 0.0d0
      p3(0) = p1(0) * p2(0)
      do 100 mord=1,maxord
        call mypmadd(p1(1),mord,p2(0),p3(1))
        call mypmadd(p2(1),mord,p1(0),p3(1))
        do 200 nord1 = 1,mord-1
          nord2 = mord - nord1
          call myproduct(p1(1),nord1,p2(1),nord2,p3(1))
  200   continue
  100 continue
      return
      end
c
      subroutine mypmadd(f,n,coeff,h)
      implicit double precision (a-h,o-z)
      dimension f(209),h(209)
      include 'len.inc'
      include 'lims.inc'
      if(coeff.eq.1.d0) goto 20
      do 10 i=len(n-1)+1,len(n)
        h(i) = h(i) + f(i)*coeff
 10   continue
      return
 20   continue
      do 30 i = len(n-1)+1, len(n)
        h(i) = h(i) + f(i)
 30   continue
      return
      end
c
      subroutine myproduct(a,na,b,nb,c)
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'len.inc'
      include 'expon.inc'
      include 'vblist.inc'
      dimension a(209),b(209),c(209),l(6)
      if(na.eq.1) then
        ia1 = 1
      else
        ia1 = len(na-1)+1
      endif
      if(nb.eq.1) then
        ib1 = 1
      else
        ib1 = len(nb-1)+1
      endif
       do 200 ia=ia1,len(na)
           if(a(ia).eq.0.d0) goto 200
           do 20 ib = ib1,len(nb)
               if(b(ib).eq.0.d0) goto 20
               do 2 m=1,6
                   l(m) = expon(m,ia) +  expon(m,ib)
   2            continue
                n = ndex(l)
                c(n) = c(n) + a(ia)*b(ib)
  20        continue
 200   continue
       return
       end
c
***********************************************************************
c
      subroutine nndrift(l,phideg,h,mh)
c
c     Transverse entry drift.
c     generates linear matrix mh and
c     array h containing nonlinearities
c     for the transfer map describing
c     a drift section of length l meters,
c     with the design trajectory at an angle of phideg degrees
c     F. Neri 5/24/89
c     Actual code generated by Johannes Van Zeijts using
c     REDUCE.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      double precision l,h(monoms),mh(6,6)
c
      double precision lsc
      dimension j(6)
c
      DOUBLE PRECISION UL
      DOUBLE PRECISION UB
      DOUBLE PRECISION CO
      DOUBLE PRECISION Si
      DIMENSION UDERS(923)
      DOUBLE PRECISION UDERS
      DOUBLE PRECISION T1
c
      call clear(h,mh)
      lsc=l/sl
c
      phi = phideg*pi180
      SI = SIN(phi)
      CO = COS(phi)
c
c     add drift terms to mh
c
      do 40 k=1,6
      mh(k,k)=+1.0d0
   40 continue
      mh(1,2)=+lsc*(1.d0/CO + SI**2/CO**3)
      mh(1,6)=+lsc*SI/(beta*CO**3)
      mh(3,4)=+lsc*(1.d0/CO)
      mh(5,2)=+lsc*SI/(beta*CO**3)
      mh(5,6)=+lsc*(-1.d0/CO+1.d0/(beta**2*CO**3))
c     mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     add drift terms to h
      UL = lsc
      UB = beta
c
      do 1 i=1,923
        UDERS(i) = 0.0d0
    1 continue
c
c  From: CINCOM::JOHANNES     "Johannes van Zeijts" 23-MAY-1989 16:38
c
      UDERS(923) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7)) + (-2.1875D0
     &    ) * (UL / (UB ** 4 * CO ** 9)) + 1.3125D0 * (UL / (UB ** 6 
     &    * CO ** 11)) + (-6.25D-02) * (UL / CO ** 5)
      UDERS(910) = (-1.875D0) * (UL / (UB ** 2 * CO ** 7)) + 2.1875D0 
     &    * (UL / (UB ** 4 * CO ** 9)) + 0.1875D0 * (UL / CO ** 5)
      UDERS(901) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7)) + (-0.1875D0
     &    ) * (UL / CO ** 5)
      UDERS(896) = 6.25D-02 * (UL / CO ** 5)
      UDERS(839) = 1.875D0 * ((UL * SI) / (UB * CO ** 7)) + (
     &    -8.75D0) * ((UL * SI) / (UB ** 3 * CO ** 9)) + 7.875D0 * 
     &    ((UL * SI) / (UB ** 5 * CO ** 11))
      UDERS(828) = (-3.75D0) * ((UL * SI) / (UB * CO ** 7)) + 
     &    8.75D0 * ((UL * SI) / (UB ** 3 * CO ** 9))
      UDERS(821) = 1.875D0 * ((UL * SI) / (UB * CO ** 7))
      UDERS(783) = (-1.875D0) * (UL / (UB ** 2 * CO ** 7)) + 2.1875D0 
     &    * (UL / (UB ** 4 * CO ** 9)) + 0.1875D0 * (UL / CO ** 5) 
     &    + (-13.125D0) * ((UL * SI ** 2) / (UB ** 2 * CO ** 9)) + 
     &    19.6875D0 * ((UL * SI ** 2) / (UB ** 4 * CO ** 11)) + 
     &    0.9375D0 * ((UL * SI ** 2) / CO ** 7)
      UDERS(774) = 1.875D0 * (UL / (UB ** 2 * CO ** 7)) + (-0.375D0) 
     &    * (UL / CO ** 5) + 13.125D0 * ((UL * SI ** 2) / (UB ** 2 
     &    * CO ** 9)) + (-1.875D0) * ((UL * SI ** 2) / CO ** 7)
      UDERS(769) = 0.1875D0 * (UL / CO ** 5) + 0.9375D0 * ((UL * SI
     &     ** 2) / CO ** 7)
      UDERS(748) = (-3.75D0) * ((UL * SI) / (UB * CO ** 7)) + 
     &    8.75D0 * ((UL * SI) / (UB ** 3 * CO ** 9)) + (-8.75D0) * 
     &    ((UL * SI ** 3) / (UB * CO ** 9)) + 26.25D0 * ((UL * SI
     &     ** 3) / (UB ** 3 * CO ** 11))
      UDERS(741) = 3.75D0 * ((UL * SI) / (UB * CO ** 7)) + 8.75D0 * 
     &    ((UL * SI ** 3) / (UB * CO ** 9))
      UDERS(728) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7)) + (-0.1875D0
     &    ) * (UL / CO ** 5) + 13.125D0 * ((UL * SI ** 2) / (UB ** 
     &    2 * CO ** 9)) + (-1.875D0) * ((UL * SI ** 2) / CO ** 7) 
     &    + 19.6875D0 * ((UL * SI ** 4) / (UB ** 2 * CO ** 11)) + (
     &    -2.1875D0) * ((UL * SI ** 4) / CO ** 9)
      UDERS(723) = 0.1875D0 * (UL / CO ** 5) + 1.875D0 * ((UL * SI 
     &    ** 2) / CO ** 7) + 2.1875D0 * ((UL * SI ** 4) / CO ** 9
     &    )
      UDERS(718) = 1.875D0 * ((UL * SI) / (UB * CO ** 7)) + 8.75D0 
     &    * ((UL * SI ** 3) / (UB * CO ** 9)) + 7.875D0 * ((UL * 
     &    SI ** 5) / (UB * CO ** 11))
      UDERS(714) = 6.25D-02 * (UL / CO ** 5) + 0.9375D0 * ((UL * SI
     &     ** 2) / CO ** 7) + 2.1875D0 * ((UL * SI ** 4) / CO ** 
     &    9) + 1.3125D0 * ((UL * SI ** 6) / CO ** 11)
      UDERS(461) = 0.375D0 * (UL / (UB * CO ** 5)) + (-1.25D0) * (UL 
     &    / (UB ** 3 * CO ** 7)) + 0.875D0 * (UL / (UB ** 5 * CO **
     &     9))
      UDERS(450) = (-0.75D0) * (UL / (UB * CO ** 5)) + 1.25D0 * (UL / 
     &    (UB ** 3 * CO ** 7))
      UDERS(443) = 0.375D0 * (UL / (UB * CO ** 5))
      UDERS(405) = (-3.75D0) * ((UL * SI) / (UB ** 2 * CO ** 7)) + 
     &    4.375D0 * ((UL * SI) / (UB ** 4 * CO ** 9)) + 0.375D0 * (
     &    (UL * SI) / CO ** 5)
      UDERS(396) = 3.75D0 * ((UL * SI) / (UB ** 2 * CO ** 7)) + (
     &    -0.75D0) * ((UL * SI) / CO ** 5)
      UDERS(391) = 0.375D0 * ((UL * SI) / CO ** 5)
      UDERS(370) = (-0.75D0) * (UL / (UB * CO ** 5)) + 1.25D0 * (UL / 
     &    (UB ** 3 * CO ** 7)) + (-3.75D0) * ((UL * SI ** 2) / (UB 
     &    * CO ** 7)) + 8.75D0 * ((UL * SI ** 2) / (UB ** 3 * CO 
     &    ** 9))
      UDERS(363) = 0.75D0 * (UL / (UB * CO ** 5)) + 3.75D0 * ((UL * 
     &    SI ** 2) / (UB * CO ** 7))
      UDERS(350) = 3.75D0 * ((UL * SI) / (UB ** 2 * CO ** 7)) + (
     &    -0.75D0) * ((UL * SI) / CO ** 5) + 8.75D0 * ((UL * SI 
     &    ** 3) / (UB ** 2 * CO ** 9)) + (-1.25D0) * ((UL * SI ** 3
     &    ) / CO ** 7)
      UDERS(345) = 0.75D0 * ((UL * SI) / CO ** 5) + 1.25D0 * ((UL * 
     &    SI ** 3) / CO ** 7)
      UDERS(340) = 0.375D0 * (UL / (UB * CO ** 5)) + 3.75D0 * ((UL * 
     &    SI ** 2) / (UB * CO ** 7)) + 4.375D0 * ((UL * SI ** 4) 
     &    / (UB * CO ** 9))
      UDERS(336) = 0.375D0 * ((UL * SI) / CO ** 5) + 1.25D0 * ((UL 
     &    * SI ** 3) / CO ** 7) + 0.875D0 * ((UL * SI ** 5) / 
     &    CO ** 9)
      UDERS(209) = (-0.75D0) * (UL / (UB ** 2 * CO ** 5)) + 0.625D0 * 
     &    (UL / (UB ** 4 * CO ** 7)) + 0.125D0 * (UL / CO ** 3)
      UDERS(200) = 0.75D0 * (UL / (UB ** 2 * CO ** 5)) + (-0.25D0) * 
     &    (UL / CO ** 3)
      UDERS(195) = 0.125D0 * (UL / CO ** 3)
      UDERS(174) = (-1.5D0) * ((UL * SI) / (UB * CO ** 5)) + 2.5D0 
     &    * ((UL * SI) / (UB ** 3 * CO ** 7))
      UDERS(167) = 1.5D0 * ((UL * SI) / (UB * CO ** 5))
      UDERS(154) = 0.75D0 * (UL / (UB ** 2 * CO ** 5)) + (-0.25D0) * 
     &    (UL / CO ** 3) + 3.75D0 * ((UL * SI ** 2) / (UB ** 2 * 
     &    CO ** 7)) + (-0.75D0) * ((UL * SI ** 2) / CO ** 5)
      UDERS(149) = 0.25D0 * (UL / CO ** 3) + 0.75D0 * ((UL * SI ** 
     &    2) / CO ** 5)
      UDERS(144) = 1.5D0 * ((UL * SI) / (UB * CO ** 5)) + 2.5D0 * (
     &    (UL * SI ** 3) / (UB * CO ** 7))
      UDERS(140) = 0.125D0 * (UL / CO ** 3) + 0.75D0 * ((UL * SI **
     &     2) / CO ** 5) + 0.625D0 * ((UL * SI ** 4) / CO ** 7)
      UDERS(83) = (-0.5D0) * (UL / (UB * CO ** 3)) + 0.5D0 * (UL / (
     &    UB ** 3 * CO ** 5))
      UDERS(76) = 0.5D0 * (UL / (UB * CO ** 3))
      UDERS(63) = 1.5D0 * ((UL * SI) / (UB ** 2 * CO ** 5)) + (
     &    -0.5D0) * ((UL * SI) / CO ** 3)
      UDERS(58) = 0.5D0 * ((UL * SI) / CO ** 3)
      UDERS(53) = 0.5D0 * (UL / (UB * CO ** 3)) + 1.5D0 * ((UL * SI
     &     ** 2) / (UB * CO ** 5))
      UDERS(49) = 0.5D0 * ((UL * SI) / CO ** 3) + 0.5D0 * ((UL * 
     &    SI ** 3) / CO ** 5)
c
      do 2 i=1, monoms
        h(i) = -UDERS(i)
    2 continue
      return
      end
c
***********************************************************************
c
      subroutine myprot(phideg,h,mh)
c
c     High order PROT routine.
c     Actual code generated by Johannes van Zeijts using
c     REDUCE.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      double precision l,h(monoms),mh(6,6)
c
      dimension j(6)
c
      DOUBLE PRECISION B
      DOUBLE PRECISION CO
      DOUBLE PRECISION Si
c
      call clear(h,mh)
c
      B = beta
      phi = phideg*pi180
      SI = SIN(phi)
      CO = COS(phi)
c
      mh(1,1)=1.0d0/CO
      mh(2,2)=CO
      mh(2,6)=(-SI)/B
      mh(3,3)=1.0d0
      mh(4,4)=1.0d0
      mh(5,1)=SI/(B*CO)
      mh(5,5)=1.0d0
      mh(6,6)=1.0d0
c
      h(34) =(-SI)/(2.0d0*CO)
      h(43) =(-SI)/(2.0d0*CO)
      h(48) =(SI*(B**2-1))/(2.0d0*B**2*CO)
      h(105) =SI**2/(4.0d0*(SI**2-1))
      h(109) =(-SI)/(2.0d0*B*CO)
      h(114) =SI**2/(4.0d0*(SI**2-1))
      h(119) =(SI**2*(-B**2+1))/(4.0d0*B**2*(SI**2-1))
      h(132) =(-SI)/(2.0d0*B*CO)
      h(139) =(SI*(B**2-1))/(2.0d0*B**3*CO)
      h(266) =SI/(8.0d0*CO*(SI**2-1))
      h(270) =SI**2/(2.0d0*B*(SI**2-1))
      h(275) =(SI*(-2.0d0*SI**2+3))/(12.0d0*CO*(SI**2-1))
      h(280) =(SI*(2.0d0*B**2*SI**2-3.0d0*B**2-8.0d0*SI**2+9))/(
     & 12.0d0*B**2*CO*(SI**2-1))
      h(293) =SI**2/(2.0d0*B*(SI**2-1))
      h(300) =(SI**2*(-B**2+1))/(2.0d0*B**3*(SI**2-1))
      h(321) =(SI*(-4.0d0*SI**2+3))/(24.0d0*CO*(SI**2-1))
      h(326) =(SI*(4.0d0*B**2*SI**2-3.0d0*B**2-10.0d0*SI**2+9))/(
     & 12.0d0*B**2*CO*(SI**2-1))
      h(335) =(SI*(-4.0d0*B**4*SI**2+3.0d0*B**4+20.0d0*B**2*SI**2-
     & 18.0d0*B**2-16.0d0*SI**2+15))/(24.0d0*B**4*CO*(SI**2-1))
      h(588) =(SI**2*(SI**2-4))/(32.0d0*(SI**4-2.0d0*SI**2+1))
      h(592) =(SI*(SI**2+6))/(16.0d0*B*CO*(SI**2-1))
      h(597) =(SI**2*(3.0d0*SI**2-4))/(16.0d0*(SI**4-2.0d0*SI**2+
     & 1))
      h(602) =(SI**2*(-3.0d0*B**2*SI**2+4.0d0*B**2+15.0d0*SI**2-16.0d0
     & ))/(16.0d0*B**2*(SI**4-2.0d0*SI**2+1))
      h(615) =(3.0d0*SI*(-SI**2+2))/(8.0d0*B*CO*(SI**2-1))
      h(622) =(SI*(3.0d0*B**2*SI**2-6.0d0*B**2-7.0d0*SI**2+10))/(
     & 8.0d0*B**3*CO*(SI**2-1))
      h(643) =(SI**2*(5.0d0*SI**2-4))/(32.0d0*(SI**4-2.0d0*SI**2+
     & 1))
      h(648) =(SI**2*(-5.0d0*B**2*SI**2+4.0d0*B**2+17.0d0*SI**2-16.0d0
     & ))/(16.0d0*B**2*(SI**4-2.0d0*SI**2+1))
      h(657) =(SI**2*(5.0d0*B**4*SI**2-4.0d0*B**4-34.0d0*B**2*SI**2
     & +32.0d0*B**2+29.0d0*SI**2-28))/(32.0d0*B**4*(SI**4-2.0d0*SI**2+
     & 1))
      h(695) =(SI*(-7.0d0*SI**2+6))/(16.0d0*B*CO*(SI**2-1))
      h(702) =(SI*(7.0d0*B**2*SI**2-6.0d0*B**2-11.0d0*SI**2+10))/(
     & 8.0d0*B**3*CO*(SI**2-1))
      h(713) =(SI*(-7.0d0*B**4*SI**2+6.0d0*B**4+22.0d0*B**2*SI**2-
     & 20.0d0*B**2-15.0d0*SI**2+14))/(16.0d0*B**5*CO*(SI**2-1))
c
      call revf(1,h,mh)
c
      return
      end
c
c end of file
