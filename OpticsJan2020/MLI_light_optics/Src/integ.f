************************************************************************
* header              INTEGRATE  (GENMAP for a general string of       *
*                     magnets with soft fringe fields)                 *
*  All routines needed for this special GENMAP                         *
************************************************************************
      subroutine integ(p,fa,fm)
c This is a subroutine for computing the map for a soft edged dipole
c magnet ( F. Neri 5/16/89 ).
c Interface to subroutine BONAX by P. Walstrom. 
c Computes By (and derivatives on axis) for a sinPHI 
c (Walstrom) steering magnet. 2/15/90.
c Eventually use GRONAX for m diff. from 1.
c Done (F. Neri 3/13/90).
c Generalized to arbitrary m, etc. 3/23/90, F. Neri
c Added type 8 (REC dipole sheets).
c Corrected sign problem.
c 11-19-90.
c Interface to CFQD F. Neri 5/8/91.
c Spacers, etc. F. Neri 6/28/91.
c Thick Halbach magnets. 6/28/91.
c
c included in ML/I 
c
      use lieaparam, only : monoms
      use beamdata
cryneneriwalstrom      include 'param.inc'
      include 'impli.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
cryneneriwalstrom      include 'parm.inc'
      include 'files.inc'
      include 'dip.inc'
      include 'pie.inc'
      include 'gronax.inc'
      include 'multipole.inc'
      include 'zz.inc'
      include 'vecpot.inc'
      include 'nnprint.inc'
c
c  calling arrays
      dimension p(6)
      dimension fa(monoms), fm(6,6)
c      
c  local arrays      
      dimension ha(monoms), hm(6,6)
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
c
      zlength = xldr + ztotal
c
      nops = nint(p(2))
      shflag = nint(p(3))
      ashield = p(4)
      ns = nint(p(5))
      lun = nint(p(1))
      if (lun.gt.10) then
        nnprint = 2*nint(p(6))
      endif
      if ( nops .eq. 0 ) then
c          write (6,*) ' CFQD option not implemented!!!'
           call multcfq(zlength,fa,fm)
           return
      else if ( nops .lt. 0 ) then
           call mulplt(zlength)
c          write (6,*) ' PLOT option not implemented!!!'
           return
      endif
      zi = 0.d0
      zf = zlength
      h=(zf-zi)/float(ns)
c
c     call VAX system routine for timing report
c
      ttaa = secnds(0.0)
c
c  initial values for design orbit (in dimensionless units) :
c
      do 2 i = 1,6
         y(i) = 0.0d0
  2   continue
c
c I really don't understand this! FN:
c Actually it just means that the design trajectory is 
c ALWAYS on momentum.
c
      y(6)=-1.d0/beta
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c
c  initialize map to the identity map:
cryneneriwalstrom      ne=224
      ne=monoms+15
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c  do the computation:
      t=zi
      iflag = 3
cryne 1 August 2004 fix later:
cryne call adam11(h,ns,'start',t,y)
      nedummy=monoms+15
      call adam11(h,ns,'start',t,y,nedummy)
      call errchk(y(7))
      call putmap(y,fa,fm)
      s1 = -(y(2)-Ax(0))
      phi1 = dasin(s1)
      phi1deg = phi1/pi180
cx for test only the following lines are commented out:
c      call prot(phi1deg,2,ha,hm)
c      call concat(fa,fm,ha,hm,fa,fm)
c
      write(jof,*) ' Final conditions in fixed reference frame!'
      write(jodf,*) ' Final conditions in fixed reference frame!'
c
      do 22  i = 1,6
         zfinal(i) = y(i)
  22  continue   
c
      write(jof , 999) (zfinal(i),i=1,6)
      write(jodf, 999) (zfinal(i),i=1,6)
 999  format(' zf=',6e12.5)
c       
      write(jof , 991) phi1deg
      write(jodf, 991) phi1deg
991   format('  Final angle is ',f16.8,' degrees')
c
c     call VAX system routine for timing report
c
      ttbb = secnds(ttaa)
       write( jof,567) ttbb
       write(jodf,567) ttbb
 567  format(' Integration time = ',f12.2,' sec.')
c      call vecpot(Ax,Ay)
cryneneriwalstrom: added return statement
      return
      end
c
**********************************************************************
      subroutine vecpot(Ax,Ay)
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
cryneneriwalstrom      include 'parm.inc'
cryneneriwalstrom      include 'param.inc'
      include 'files.inc'
      double precision Ax(0:monoms), Ay(0:monoms)
      do 1 i=1,monoms
        if(Ax(i).ne. 0.d0) write(jodf,*) ' Ax(',i,')=',Ax(i)
        if(Ay(i).ne. 0.d0) write(jodf,*) ' Ay(',i,')=',Ay(i)
 1    continue
      return
      end
c
**********************************************************************
c
       subroutine gradients(z0)
c
c  F. Neri
c  Revised version 7/12/90.
c  Corrected bug    8/6/90.
c  Skew Multipoles  5/27/91.
c
       include 'impli.inc'
       include 'gronax.inc'
       include 'multipole.inc'
       include 'dip.inc'
       integer ifirst,nderiv
       logical initg
       data ifirst/0/
       double precision bb(14)
c
cryne this looks like old code that is not needed here,
cryne so I am commenting it out
cryne  data blprod/2.d0/,a1/0.75d0/,a2/1.d0/,xl/2.d0/,ss/1.d0/
cryne  data nmax/100/,kmax/100/
c
       save ifirst
c
c  Initialization, if necessary:
       if (ifirst .eq. 0 ) then
         do 10 j=1, ncoil
           initg = .false.
           jtyp = itype(j)
           if(jtyp.eq.3)  initg = .true.
           if(jtyp.eq.4)  initg = .true.
           if(jtyp.eq.6)  initg = .true.
           if(jtyp.eq.17) initg = .true.
cryneneriwalstrom      if(initg) call onaxgr(0,0.d0,j,bb)
           if(initg) call onaxgr(0,0.d0,j,nderiv,bb)
  10     continue
         ifirst = 1
       endif
c
ctm    plot net gradient in file 79
c
       gnet = 0.0d0
       do 1 j = 1, ncoil
c          zz = z0 - zcoil(j)
          zz = z0 - xldr
c          write(6,*) j,z0,zz
          call  onaxgr(1,zz,j,11,bb)
          gnet = gnet + bb(1)
c          write(79,*)  zz, bb(1)
c
          if ( mcoil(j) .gt. 0 ) then
            do 666 i= 0, 10
          gn(mcoil(j),i) = gn(mcoil(j),i) + bb(i+1)/float(mcoil(j))
 666        continue
          else
            do 777 i= 0, 10
          gs(-mcoil(j),i) = gs(-mcoil(j),i) + bb(i+1)/float(-mcoil(j))
 777        continue
          endif
c
c          write(6,*) ' z = ',zz, ' b =',bb(1)
 1     continue
       write(79,*)  zz,gnet
       return
       end
c   
**********************************************************************
c
      subroutine prenv(z,y)
      include 'impli.inc'
      include 'sigbuf.inc'
      include 'nnprint.inc'
      include 'gronax.inc'
c
      dimension y(*)      
c
      save ncount, nf
c
      if(ncount.eq.0) then
        nf = 1
        xx(1) = xxin
        ax(1) = axin
        px(1) = pxin
        yy(1) = yyin
        ay(1) = ayin
        py(1) = pyin
      endif
      ncount = ncount+1
      if (mod(ncount, nnprint).ne.0) return
c-------------------   x plane  ------------
      cfx = y(7)
      sfx = y(13)
      cpx = y(8)
      spx = y(14)
      txx = cfx**2
      tax = 2.0*cfx*sfx
      tpx = sfx**2
      txa = cfx*cpx
      taa = cfx*spx+sfx*cpx
      tpa = sfx*spx
      txp = cpx**2
      tap = 2.0*cpx*spx
      tpp = spx**2
c-------------------   y plane  ------------
      cfy = y(21)
      sfy = y(27)
      cpy = y(22)
      spy = y(28)
      ryy = cfy**2
      ray = 2.0*cfy*sfy
      rpy = sfy**2
      rya = cfy*cpy
      raa = cfy*spy+sfy*cpy
      rpa = sfy*spy
      ryp = cpy**2
      rap = 2.0*cpy*spy
      rpp = spy**2
c
c
c      compute x-plane final beam ellipse
c
      ni = 1
c
      nf = nf + 1
      zu = z
c
      xs = xx(ni)**2
      xa = -xx(ni)*ax(ni)*px(ni)/sqrt(1.0d0 + ax(ni)**2)
      xps = px(ni)**2
c      sigf1 = (cfx**2)*xs + 2.0*cfx*sfx*xa + (sfx**2)*xps
c      sigf2 = cfx*cpx*xs + (cfx*spx+sfx*cpx)*xa + sfx*spx*xps
c      sigf3 = (cpx**2)*xs + 2.0*cpx*spx*xa + (spx**2)*xps
      sigf1 = txx*xs + tax*xa + tpx*xps
      sigf2 = txa*xs + taa*xa + tpa*xps
      sigf3 = txp*xs + tap*xa + tpp*xps
      exsq = sigf1*sigf3 - sigf2**2
      if(exsq.le.0.0) exsq = 1.e-30
      ex = sqrt(exsq)
      xx(nf) = sqrt(sigf1)
      if(xx(nf).gt.xmax) then
         xmax = xx(nf)
         zxmax = zu
      endif
      if(xx(nf).gt.exmax) then
         exmax = xx(nf)
         zxm = zu
      endif
      if(xx(nf).lt.exmin) exmin = xx(nf)
      px(nf) = sqrt(sigf3)
      ax(nf) = -sigf2/ex
c
c     compute y-plane beam ellipse
c
      ys  = yy(ni)**2
      ya  = -yy(ni)*ay(ni)*py(ni)/sqrt(1.0d0 + ay(ni)**2)
      yps = py(ni)**2
c      sigf4 = (cfy**2)*ys + 2.0*cfy*sfy*ya + (sfy**2)*yps
c      sigf5 = cfy*cpy*ys + (cfy*spy+sfy*cpy)*ya + sfy*spy*yps
c      sigf6 = (cpy**2)*ys + 2.0*cpy*spy*ya + (spy**2)*yps
      sigf4 = ryy*ys + ray*ya + rpy*yps
      sigf5 = rya*ys + raa*ya + rpa*yps
      sigf6 = ryp*ys + rap*ya + rpp*yps
      eysq = sigf4*sigf6 - sigf5**2
      if(eysq.le.0.0) eysq = 1.e-30
      ey = sqrt(eysq)
      yy(nf) = sqrt(sigf4)
      if(yy(nf).gt.ymax) then
         ymax = yy(nf)
         zymax = zu
      endif
      if(yy(nf).gt.eymax) then
         eymax = yy(nf)
         zym = zu
      endif
      if(yy(nf).lt.eymin) eymin = yy(nf)
      py(nf) = sqrt(sigf6)
      ay(nf) = -sigf5/ey
c
      quadp = gn(2,0)*2.
      octp  = gn(4,0)*4.
c
      if(lun.gt.0) write(lun,27) zu,xx(nf),px(nf),yy(nf),py(nf),ax(nf)
     * ,ay(nf), quadp, octp
  27  format(f10.4,4(1pe11.4),4(1pe12.4))
      return
      end
c   
**********************************************************************
c
      subroutine hmltn3(t,y,h)
c  Originally written by F. Neri, 5/16/89.
c  Modified 6/3/89 to handle different gauges.
c  Modified 3/13/90 for arbitrary sin(m theta)+cos(m theta)
c  magnets, as given in common block GRONAX.
c  Interfaces to routine GRONAX by P. Walstrom.
c  Modified again 3/23/90 for new input format.
c
      use lieaparam, only : monoms
      use beamdata
cryneneriwalstrom It appears that iplus is not used. Comment out later.
cryneneriwalstrom Need to set itop to whatever makes sense for this version
c
cryneneriwalstrom      parameter (itop=209,iplus=224)
      parameter (itop=923,iplus=itop+15)
      include 'impli.inc'
cryneneriwalstrom      include 'param.inc'
cryneneriwalstrom      include 'parm.inc'
      include 'dip.inc'
      include 'bfield.inc'
      include 'gronax.inc'
      include 'multipole.inc'
      include 'vecpot.inc'
      include 'nnprint.inc'
c
      dimension h(monoms),y(*)
c
      dimension X(0:itop),YY(0:itop),P1(0:itop),P2(0:itop)
      dimension P12(0:itop),P22(0:itop)
      dimension A(0:12)
      double precision bb(14)
c
c  begin calculation
c
c  compute gradients on axis:
c
      do 10 i=0,10
        do 10 j=0,10
          gn(i,j) = 0.0d0
          gs(i,j) = 0.0d0
 10   continue
c
       zz = t
c       write(6,*) zz
c      call dbonax(1,a1,a2,xl,blprod,ss,nmax,kmax,ndriv,zz,bb)
c
       call gradients(zz)
c
c  print out envelopes, gradients, etc.
      if(lun.gt.10) then
        call prenv(zz,y)
      endif
c
c  initialization
      do 11 i=1,monoms
  11    h(i)=0.d0
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
        A(i) =  A(i-1) * (1.d0/2.d0 -(i-1.d0))/i
   50 continue
c
      do 60 i=0,monoms
        X(i) = 0.d0
   60 continue
      X(6) = -2.d0 / beta
c      X(13) = -1.d0
c      X(22) = -1.d0
      X(27) =  1.d0
c
cryneneriwalstrom      maxord = 4
cryneneriwalstrom      nn = 4
      maxord = 6
      nn = 6
c P1 = px -q Ax
      do 101 i=0,monoms
      Ax(i) = 0.d0
      Ay(i) = 0.d0
      Az(i) = 0.d0
      P2(i) = 0.d0
  101 P1(i) = 0.d0
c
      q = 1.d0/brho
      x0 = y(1)
      y0 = y(3)
c
c     write(6,*) ' x0=',x0,' y0=',y0
c 
c Mathematica generated block follows:
c      include 'a3.fop'
c NEED A SWITCH HERE THAT SAYS,  "IF(not an rf cavity)THEN
      call a3(q,x0,y0)
c ELSE
c     call rf cavity routines
c
      cos2 = 1.d0 -(y(2) -Ax(0))**2  -(y(4)-Ay(0))**2
c
c Note constant term in px -q Ax, coming from design orbit 
c y(2) = Px = sin(phi) (initially).
      P1(0) = y(2)
      P1(2) = 1.d0
c P12 = P1**2
      do 112 i=0,itop
  112 P1(i) = P1(i)-Ax(i)
      call pmult(P1,P1,P12,maxord)
c
      P2(0) = y(4)
      P2(4) = 1.d0
      do 122 i=0,itop
  122 P2(i) = P2(i)-Ay(i)
      call pmult(P2,P2,P22,maxord)
c
c Zeroth order term subtracted ( Sum starts at 1 ):
c Divide by cos2
      do 102 i=1,itop
  102 X(i) = (X(i)-P12(i)-P22(i))/cos2
c X = pt**2 -2/beta*pt -(px -q Ax)**2 -(py -q Ay)**2
c YY = -Sqrt(1+X)
      call poly1(nn,A,X,YY,maxord)
c h = -Az
      do 132 i=7, monoms
  132 h(i) = -Az(i)
c
c h = -Sqrt(1+X)-Az
c Zero and first order terms subtracted ( Sum starts at 7 ):
c Scale by cos
      do 70 i=7, monoms
        h(i) = h(i)+dsqrt(cos2)*YY(i)/sl
   70 continue
      return
      end
c
***********************************************************************
c
      subroutine poly1(N,A,X,Y,maxord)
      include 'lims.inc'
      parameter (MN=6)
cryne 7/23/2002      common /lims/bottom,top
cryne 7/23/3003      integer bottom(0:12),top(0:12)
      double precision X(0:923),Y(0:923)
      double precision A(0:N)

      double precision Vect(0:923,MN)

c   NO COMMENT
      do 100 i=0,top(maxord)
 100    Y(i) = 0.0d0
      do 200 i=0,top(maxord)
 200    Vect(i,1) = X(i)
      do 300 iord=2,N
        call pmult(X,Vect(0,iord-1),Vect(0,iord),maxord)
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
      subroutine pmult(p1,p2,p3,maxord)
      include 'lims.inc'
      double precision p1(0:923),p2(0:923),p3(0:923)
c
cryne 7/23/2002      common /lims/bottom,top
cryne 7/23/3003      integer bottom(0:12),top(0:12)
      if (maxord.lt.0) return
c
      do 10 i=1,top(maxord)
   10 p3(i) = 0.0d0
      p3(0) = p1(0) * p2(0)
      do 100 mord=1,maxord
        call pmadd(p1(1),mord,p2(0),p3(1))
        call pmadd(p2(1),mord,p1(0),p3(1))
        do 200 nord1 = 1,mord-1
          nord2 = mord - nord1
          call product(p1(1),nord1,p2(1),nord2,p3(1))
  200   continue
  100 continue
      return
      end
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
      use lieaparam, only : monoms
      use beamdata
c
      include 'impli.inc'
cryneneriwalstrom      include 'param.inc'
cryneneriwalstrom      include 'parm.inc'
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
      mh(1,2)=+lsc*(1.d0/CO+SI**2/CO**3)
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
c  From:    CINCOM::JOHANNES     "Johannes van Zeijts" 23-MAY-1989 16:38
c
      UDERS(923) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7))+(-2.1875D0
     &    ) * (UL / (UB ** 4 * CO ** 9))+1.3125D0 * (UL / (UB ** 6 
     &    * CO ** 11))+(-6.25D-02) * (UL / CO ** 5)
      UDERS(910) = (-1.875D0) * (UL / (UB ** 2 * CO ** 7))+2.1875D0 
     &    * (UL / (UB ** 4 * CO ** 9))+0.1875D0 * (UL / CO ** 5)
      UDERS(901) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7))+(-0.1875D0
     &    ) * (UL / CO ** 5)
      UDERS(896) = 6.25D-02 * (UL / CO ** 5)
      UDERS(839) = 1.875D0 * ((UL * SI) / (UB * CO ** 7))+(
     &    -8.75D0) * ((UL * SI) / (UB ** 3 * CO ** 9))+7.875D0 * 
     &    ((UL * SI) / (UB ** 5 * CO ** 11))
      UDERS(828) = (-3.75D0) * ((UL * SI) / (UB * CO ** 7))+
     &    8.75D0 * ((UL * SI) / (UB ** 3 * CO ** 9))
      UDERS(821) = 1.875D0 * ((UL * SI) / (UB * CO ** 7))
      UDERS(783) = (-1.875D0) * (UL / (UB ** 2 * CO ** 7))+2.1875D0 
     &    * (UL / (UB ** 4 * CO ** 9))+0.1875D0 * (UL / CO ** 5) 
     &   +(-13.125D0) * ((UL * SI ** 2) / (UB ** 2 * CO ** 9))+
     &    19.6875D0 * ((UL * SI ** 2) / (UB ** 4 * CO ** 11))+
     &    0.9375D0 * ((UL * SI ** 2) / CO ** 7)
      UDERS(774) = 1.875D0 * (UL / (UB ** 2 * CO ** 7))+(-0.375D0) 
     &    * (UL / CO ** 5)+13.125D0 * ((UL * SI ** 2) / (UB ** 2 
     &    * CO ** 9))+(-1.875D0) * ((UL * SI ** 2) / CO ** 7)
      UDERS(769) = 0.1875D0 * (UL / CO ** 5)+0.9375D0 * ((UL * SI
     &     ** 2) / CO ** 7)
      UDERS(748) = (-3.75D0) * ((UL * SI) / (UB * CO ** 7))+
     &    8.75D0 * ((UL * SI) / (UB ** 3 * CO ** 9))+(-8.75D0) * 
     &    ((UL * SI ** 3) / (UB * CO ** 9))+26.25D0 * ((UL * SI
     &     ** 3) / (UB ** 3 * CO ** 11))
      UDERS(741) = 3.75D0 * ((UL * SI) / (UB * CO ** 7))+8.75D0 * 
     &    ((UL * SI ** 3) / (UB * CO ** 9))
      UDERS(728) = 0.9375D0 * (UL / (UB ** 2 * CO ** 7))+(-0.1875D0
     &    ) * (UL / CO ** 5)+13.125D0 * ((UL * SI ** 2) / (UB ** 
     &    2 * CO ** 9))+(-1.875D0) * ((UL * SI ** 2) / CO ** 7) 
     &   +19.6875D0 * ((UL * SI ** 4) / (UB ** 2 * CO ** 11))+(
     &    -2.1875D0) * ((UL * SI ** 4) / CO ** 9)
      UDERS(723) = 0.1875D0 * (UL / CO ** 5)+1.875D0 * ((UL * SI 
     &    ** 2) / CO ** 7)+2.1875D0 * ((UL * SI ** 4) / CO ** 9
     &    )
      UDERS(718) = 1.875D0 * ((UL * SI) / (UB * CO ** 7))+8.75D0 
     &    * ((UL * SI ** 3) / (UB * CO ** 9))+7.875D0 * ((UL * 
     &    SI ** 5) / (UB * CO ** 11))
      UDERS(714) = 6.25D-02 * (UL / CO ** 5)+0.9375D0 * ((UL * SI
     &     ** 2) / CO ** 7)+2.1875D0 * ((UL * SI ** 4) / CO ** 
     &    9)+1.3125D0 * ((UL * SI ** 6) / CO ** 11)
      UDERS(461) = 0.375D0 * (UL / (UB * CO ** 5))+(-1.25D0) * (UL 
     &    / (UB ** 3 * CO ** 7))+0.875D0 * (UL / (UB ** 5 * CO **
     &     9))
      UDERS(450) = (-0.75D0) * (UL / (UB * CO ** 5))+1.25D0 * (UL / 
     &    (UB ** 3 * CO ** 7))
      UDERS(443) = 0.375D0 * (UL / (UB * CO ** 5))
      UDERS(405) = (-3.75D0) * ((UL * SI) / (UB ** 2 * CO ** 7))+
     &    4.375D0 * ((UL * SI) / (UB ** 4 * CO ** 9))+0.375D0 * (
     &    (UL * SI) / CO ** 5)
      UDERS(396) = 3.75D0 * ((UL * SI) / (UB ** 2 * CO ** 7))+(
     &    -0.75D0) * ((UL * SI) / CO ** 5)
      UDERS(391) = 0.375D0 * ((UL * SI) / CO ** 5)
      UDERS(370) = (-0.75D0) * (UL / (UB * CO ** 5))+1.25D0 * (UL / 
     &    (UB ** 3 * CO ** 7))+(-3.75D0) * ((UL * SI ** 2) / (UB 
     &    * CO ** 7))+8.75D0 * ((UL * SI ** 2) / (UB ** 3 * CO 
     &    ** 9))
      UDERS(363) = 0.75D0 * (UL / (UB * CO ** 5))+3.75D0 * ((UL * 
     &    SI ** 2) / (UB * CO ** 7))
      UDERS(350) = 3.75D0 * ((UL * SI) / (UB ** 2 * CO ** 7))+(
     &    -0.75D0) * ((UL * SI) / CO ** 5)+8.75D0 * ((UL * SI 
     &    ** 3) / (UB ** 2 * CO ** 9))+(-1.25D0) * ((UL * SI ** 3
     &    ) / CO ** 7)
      UDERS(345) = 0.75D0 * ((UL * SI) / CO ** 5)+1.25D0 * ((UL * 
     &    SI ** 3) / CO ** 7)
      UDERS(340) = 0.375D0 * (UL / (UB * CO ** 5))+3.75D0 * ((UL * 
     &    SI ** 2) / (UB * CO ** 7))+4.375D0 * ((UL * SI ** 4) 
     &    / (UB * CO ** 9))
      UDERS(336) = 0.375D0 * ((UL * SI) / CO ** 5)+1.25D0 * ((UL 
     &    * SI ** 3) / CO ** 7)+0.875D0 * ((UL * SI ** 5) / 
     &    CO ** 9)
      UDERS(209) = (-0.75D0) * (UL / (UB ** 2 * CO ** 5))+0.625D0 * 
     &    (UL / (UB ** 4 * CO ** 7))+0.125D0 * (UL / CO ** 3)
      UDERS(200) = 0.75D0 * (UL / (UB ** 2 * CO ** 5))+(-0.25D0) * 
     &    (UL / CO ** 3)
      UDERS(195) = 0.125D0 * (UL / CO ** 3)
      UDERS(174) = (-1.5D0) * ((UL * SI) / (UB * CO ** 5))+2.5D0 
     &    * ((UL * SI) / (UB ** 3 * CO ** 7))
      UDERS(167) = 1.5D0 * ((UL * SI) / (UB * CO ** 5))
      UDERS(154) = 0.75D0 * (UL / (UB ** 2 * CO ** 5))+(-0.25D0) * 
     &    (UL / CO ** 3)+3.75D0 * ((UL * SI ** 2) / (UB ** 2 * 
     &    CO ** 7))+(-0.75D0) * ((UL * SI ** 2) / CO ** 5)
      UDERS(149) = 0.25D0 * (UL / CO ** 3)+0.75D0 * ((UL * SI ** 
     &    2) / CO ** 5)
      UDERS(144) = 1.5D0 * ((UL * SI) / (UB * CO ** 5))+2.5D0 * (
     &    (UL * SI ** 3) / (UB * CO ** 7))
      UDERS(140) = 0.125D0 * (UL / CO ** 3)+0.75D0 * ((UL * SI **
     &     2) / CO ** 5)+0.625D0 * ((UL * SI ** 4) / CO ** 7)
      UDERS(83) = (-0.5D0) * (UL / (UB * CO ** 3))+0.5D0 * (UL / (
     &    UB ** 3 * CO ** 5))
      UDERS(76) = 0.5D0 * (UL / (UB * CO ** 3))
      UDERS(63) = 1.5D0 * ((UL * SI) / (UB ** 2 * CO ** 5))+(
     &    -0.5D0) * ((UL * SI) / CO ** 3)
      UDERS(58) = 0.5D0 * ((UL * SI) / CO ** 3)
      UDERS(53) = 0.5D0 * (UL / (UB * CO ** 3))+1.5D0 * ((UL * SI
     &     ** 2) / (UB * CO ** 5))
      UDERS(49) = 0.5D0 * ((UL * SI) / CO ** 3)+0.5D0 * ((UL * 
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
