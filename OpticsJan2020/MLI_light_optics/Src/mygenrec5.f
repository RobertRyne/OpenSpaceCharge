************************************************************************
* header              GENREC (GENMAP for a pattern of REC quads)       *
*  All routines needed for this special GENMAP                         *
************************************************************************
************************************************************************
      subroutine gnrec3(p,fa,fm)
c This is a subroutine for computing the map for patterns of REC quads.
c It is based on Rob Ryne's GENMAP as modified by Alex Dragt 12/18/86.
c
c Further modified 4 June 87 by Tom Mottershead to allow any number
c of cycles through a pattern of up to 5 REC quads, with arbitrary
c numbers of repeats. Actually probably more general than this. AJD
c
c Modified by Filippo Neri Jan. 16 1989 to include multipoles.
c Fifth order version by F. Neri Mar 13 1989.
c Modified to be like MaryLie 3.0 version by A. Dragt 4/24/95
c
c The input parameters come from recn, and the 
c auxilary parameter NPQUAD, MULTIPOLES, plus NQT pset defining the quads.

*******************************************************************************

c
c  The parameters of recm are:
c  1.    zi  Initial integration point.
c  2.    zf  Final integration point.
c  3.    NS  Number of integration steps.
c  4.    IFILE (profile file number).   
c  5.    MULTIPOLES index of pset containin multipole information
c                   for this integration region. If MULTIPOLES is 0, then
c                   the value of all multipoles is set to zero.
c  6.    NPQUAD     index of pset containing quad pattern information.
c---------------------------------------------------------------
c   The pset MULTIPOLES has the same format as used in the cfq element:
c  1.    Normal sextupole ( Tesla/meter^2 ).
c  2.    Skew   sextupole ( Tesla/meter^2 ).
c  3.    Normal octupole  ( Tesla/meter^3 ).
c  4.    Skew   octupole  ( Tesla/Meter^3 ).
c  5.    UNUSED ( must be there!)
c  6.    UNUSED ( nust be there!)
c---------------------------------------------------------------
c  The parameter set NPQUAD defines the set of Halbach quad units as follows:
c  1.    Di   Initial drift to first quad ( starting from Z = 0.)
c  2.    NQT  Number of psets used for quads.
c  3.    IPS  Index of inital pset used for quads. The following
c             psets in order define the successive quads.
c  4.    MAXQ Maximum number of quads actually used, including multiple
c             cycles thru pattern (but not repeats).
c  5.    NCYC Number of cycles thru pattern, except that the sequence
c             stops when MAXQ is reached.
c  6.    ISEND     output control:  0 = quiet running
c                                   1 = print on terminal (jof)
c                                   2 = print on std. output file (jodf)
c                                   3 = print on both
c---------------------------------------------------------------
c     The parameter set type codes are used to define the Halbach quad unit
c     cells in a manner parallel to the normal quad type codes:
c      ps(1) = ql, the quad length in meters.
c      ps(2) = fg, the field gradient in Tesla/meter (if the quad were
c                  infinitely long).
c      ps(3) = ra, the inner radius.   (These radii control the shape of the
c      ps(4) = rb, the outer radius.       fringe field)
c      ps(5) = td, the trailing drift to the front face of the next quad in
c                  the pattern (meters). Note: the trailing drift assigned to
c                  the last quad in the whole pattern is ignored; the
c                  integration proceeds to the final point zf anyway.
c      ps(6) = nr, the number of consecutive repetitions, or multiplicity,
c                  of this quad unit cell in the basic pattern.
c-------------------------------------------------------------------
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'quadpn.inc'
      include 'files.inc'
      include 'recmul.inc'
c
c  calling arrays
      dimension p(6)
      dimension fa(monoms), fm(6,6)
c
c  local arrays
      dimension pb(6)
      dimension y(monoms+15)
c
c timing variables
c      real ttaa, ttbb
c
c use equivalence statement to make the various parameter sets pstj
c available as if they were in a two dimensional array
c
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c  y(225-*) = f5
c  y(*-*) = f6
c
c  get interval and number of steps from GENREC parameters
c
      zi = p(1)
      zf = p(2)
      ns = nint(p(3))
      ifile = nint(p(4))
      mpole = nint(p(5))
      npquad = nint(p(6))
c
c  get multipole values from the parameter set mpole
c  Note that if mpole is zero, then the multipoles are set to zero.
c
      if (mpole.lt.1 .or. mpole.gt.maxpst) then
        do 500 i=1,6
  500   pb(i) = 0.0d0
      else
        do 600 i=1,6
  600   pb(i) = pst(i,mpole)
      endif
c
c compute multipole coefficients
c
      bsex=pb(1)
      asex=pb(2)
      boct=pb(3)
      aoct=pb(4)
      fsxnr=bsex/(3.d0*brho)
      fsxsk=asex/(3.d0*brho)
      focnr=boct/(4.d0*brho)
      focsk=aoct/(4.d0*brho)
c
c  get REC pattern information from pset npquad
c
      di = pst(1,npquad)
      nqt = nint(pst(2,npquad))
      if (nqt .gt. 6) then
      write(jof,*) ' input error: nqt=',nqt
      write(jof,*) ' this value is > 6 and therefore too large'
      call myexit
      endif
      ips = nint(pst(3,npquad))
      maxq = nint(pst(4,npquad))
      ncyc = nint(pst(5,npquad))
      isend = nint(pst(6,npquad))
      jtty = 0
      jdsk = 0
      if((isend.eq.1).or.(isend.eq.3)) jtty = 1
      if((isend.eq.2).or.(isend.eq.3)) jdsk = 1
c
      h=(zf-zi)/float(ns)
c
c     echo input parameters
c
      if(jtty.eq.1) write( jof,13) ns,zi,zf
      if(jdsk.eq.1) write(jodf,13) ns,zi,zf
  13  format(/' Integrating in ',i6,' steps from',f10.5,'=zi to',
     &f10.5,'=zf')
      if (mpole.lt.1 .or. mpole.gt.maxpst) then
      if(jtty.eq.1) write( jof,*) 
     &' mpole=0, all multipoles are zero'
      if(jdsk.eq.1) write( jodf,*) 
     &' mpole=0, all multipoles are zero'
      endif
      if (mpole.ge.1 .and. mpole.le.maxpst) then
      if(jtty.eq.1) write( jof,14) mpole,bsex,asex,boct,aoct
      if(jdsk.eq.1) write(jodf,14) mpole,bsex,asex,boct,aoct
  14  format(1x,'multipole strengths from pset',i2,':',/,
     &1x,' Sextupole:',2(1pg12.4),/,
     &1x,' Octupole: ',2(1pg12.4))
      endif
      if(jtty.eq.1) write( jof,15) npquad,ncyc,nqt,maxq,di
      if(jdsk.eq.1) write(jodf,15) npquad,ncyc,nqt,maxq,di
  15  format(' Pattern from pset',i2,':',/,
     &i4,' cycle(s) of',i3,' type(s) of section(s) with a maximum of'
     &,i3,' section(s)',/,'  di=',f10.5)
      if(jtty.eq.1) write( jof,*) ' types of section(s) are:'
      if(jdsk.eq.1) write(jodf,*) ' types of section(s) are:'
      if(jtty.eq.1) write( jof,16)
      if(jdsk.eq.1) write(jodf,16)
  16  format(1x,'pset',2x,'length',4x,'strength',3x,
     &'radii: inner',6x,'outer',6x,'tdrift',4x,'number')
c
c      inititialize the nqt quad types from the first nqt parameter sets
c
      do 20 n = ips, nqt+ips-1
      kn = n-ips+1
      wd(kn) = pst(1,n)
      ga(kn) = pst(2,n)
      ra(kn) = pst(3,n)
      rb(kn) = pst(4,n)
      dr(kn) = pst(5,n)
      nr(kn) = nint(pst(6,n))
      if(jtty.eq.1) write(jof,17)
     & n,wd(kn),ga(kn),ra(kn),rb(kn),dr(kn),nr(kn)
      if(jdsk.eq.1) write(jodf,17)
     & n,wd(kn),ga(kn),ra(kn),rb(kn),dr(kn),nr(kn)
  17  format(i4,5f12.6,i6)
  20  continue
c
c     call VAX system routine for timing report
c
c      ttaa = secnds(0.0)
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=0.d0
      y(2)=0.d0
      y(3)=0.d0
      y(4)=0.d0
      y(5)=0.d0
      y(6)=0.d0
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c      cmp2=(1.d0/(beta*gamma))**2
c
c  initialize map to the identity map:
      ne=monoms+15
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c  Set up multipoles:
c
c     compute useful numbers
c
      sl2=sl*sl
      sl3=sl*sl2
      bet2=beta*beta
      bet3=beta*bet2
      gam2=gamma*gamma
c
c  write out gradient and derivatives on file ifile:
      ipflag=0
      if (ifile .ne. 0) then
      if (ifile .lt.0 ) then
      ipflag=1
      ifile=-ifile
      endif
      do 100 i=0,ns
      s=zi + float(i)*h
      call g01234(s,g,gz,gzz,gzzz,gzzzz)
  100 write(ifile,101)s,g,gz,gzz,gzzz,gzzzz
  101 format(6(1x,1pg12.5))
      write(jof,*) ' '
      write(jof,*) ' profile written on file ',ifile
      endif
c
c  return identity map if ifile was < 0
      if (ipflag .eq. 1) then
      call ident(fa,fm)       
      return
      endif
c
c  do the computation:
      t=zi
      iflag = 2
cryne 1 August 2004 fix later:
cryne call adam11(h,ns,'start',t,y)
      nedummy=monoms+15
      call adam11(h,ns,'start',t,y,nedummy)
      call errchk(y(7))
      call putmap(y,fa,fm)
      call csym(1,fm,ans)
c
c     call VAX system routine for timing report
c
c      ttbb = secnds(ttaa)
c      if(jtty.eq.1) write( jof,567) ttbb
c      if(jdsk.eq.1) write(jodf,567) ttbb
c 567  format(' GENREC integration time = ',f12.2,' sec.')
c
      return
      end
c
**********************************************************************
c
      subroutine errchk(y)
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
cryne dimension y(monoms+15)
      dimension y(*)
      s1=1.d0 -( y(1)*y(8)-y(2)*y(7) )
      s2=1.d0 -( y(15)*y(22)-y(16)*y(21) )
      s3=1.d0 -( y(29)*y(36)-y(30)*y(35) )
      if (idproc.eq.0) then
        write(6,100) s1,s2,s3
      end if
  100 format(1x,'1. - 2 x 2 determinants = ',e14.7,1x,e14.7,1x,e14.7)
      return
      end
c
**********************************************************************
c
      subroutine hmltn2(t,y,h)
c  new version by R. Ryne 6/9/2002
c
c  this routine is used to specify h(z) for a REC quad triplet
c  This version ( f5 + f6 ) by F. Neri.
c  Jan 7 1988.
c  Multipoles added by F. Neri Mar 13 1989.
      use beamdata
      implicit double precision(a-h,o-z)
      parameter (itop=923,iplus=938)
      include 'combs.inc'
      include 'recmul.inc'
      dimension A(0:12),X(0:itop),YY(0:itop),P1(0:itop),P2(0:itop)
      dimension h(itop),y(*)
c
c  begin calculation
c
c  compute gradients
      call g01234(t,g,gz,gzz,gzzz,gzzzz)
      f=0.5d0*g
      fz=0.25d0*gz
      fzz=gzz/12.d0
      fzzz=gzzz/48.d0
      fzzzz=gzzzz/256.d0
      brho=1./qbyp
      beta=-1.d0/ptg
c
c  initialization
      do 10 i=7,itop
   10 h(i)=0.d0
c
c 2nd order
      h(27)=-(-1.d0 + beta**2)/(2.d0*beta**2*sl)
      h(13)=1.d0/(2.d0*sl)
      h(22)=1.d0/(2.d0*sl)
      h(7)=(f*sl)/brho
      h(18)=-((f*sl)/brho)
c 3rd order
      h(83)=-(-1.d0 + beta**2)/(2.d0*beta**3*sl)
      h(53)=1.d0/(2.d0*beta*sl)
      h(76)=1.d0/(2.d0*beta*sl)
c 4th order
      h(209)=(5.d0 - 6.d0*beta**2 + beta**4)/(8.d0*beta**4*sl)
      h(154)=-(-3.d0 + beta**2)/(4.d0*beta**2*sl)
      h(200)=-(-3.d0 + beta**2)/(4.d0*beta**2*sl)
      h(140)=1.d0/(8.d0*sl)
      h(149)=1.d0/(4.d0*sl)
      h(195)=1.d0/(8.d0*sl)
      h(85)=-((fz*sl**2)/brho)
      h(84)=-((fzz*sl**3)/brho)
      h(96)=-((fz*sl**2)/brho)
      h(110)=(fz*sl**2)/brho
      h(176)=(fz*sl**2)/brho
      h(175)=(fzz*sl**3)/brho
c 5th order
      h(370)=-(-5.d0 + 3.d0*beta**2)/(4.d0*beta**3*sl)
      h(450)=-(-5.d0 + 3.d0*beta**2)/(4.d0*beta**3*sl)
      h(461)=(7.d0 - 10.d0*beta**2 + 3.d0*beta**4)/(8.d0*beta**5*sl)
      h(340)=3.d0/(8.d0*beta*sl)
      h(363)=3.d0/(4.d0*beta*sl)
      h(443)=3.d0/(8.d0*beta*sl)
      h(220)=-((fz*sl**2)/(beta*brho))
      h(252)=-((fz*sl**2)/(beta*brho))
      h(284)=(fz*sl**2)/(beta*brho)
      h(412)=(fz*sl**2)/(beta*brho)
c 6th order
      h(783)=(35.d0 - 30.d0*beta**2 + 3.d0*beta**4)/(16.d0*beta**4*sl)
      h(910)=(35.d0 - 30.d0*beta**2 + 3.d0*beta**4)/(16.d0*beta**4*sl)
      h(728)=(-3.d0*(-5.d0 + beta**2))/(16.d0*beta**2*sl)
      h(901)=(-3.d0*(-5.d0 + beta**2))/(16.d0*beta**2*sl)
      h(923)=                                                           &
     &-(-21.d0+35.d0*beta**2-15.d0*beta**4+beta**6)/(16.d0*beta**6*sl)
      h(774)=(-3.d0*(-5.d0 + beta**2))/(8.d0*beta**2*sl)
      h(714)=1.d0/(16.d0*sl)
      h(723)=3.d0/(16.d0*sl)
      h(769)=3.d0/(16.d0*sl)
      h(896)=1.d0/(16.d0*sl)
      h(483)=-(fz*sl**2)/(2.d0*brho)
      h(492)=-(fz*sl**2)/(2.d0*brho)
      h(497)=((-3.d0 + beta**2)*fz*sl**2)/(2.d0*beta**2*brho)
      h(463)=(fzzz*sl**4)/brho
      h(462)=(sl**5*(fz**2 + 2.d0*brho*fzzzz))/(2.d0*brho**2)
      h(524)=-(fz*sl**2)/(2.d0*brho)
      h(563)=-(fz*sl**2)/(2.d0*brho)
      h(568)=((-3.d0 + beta**2)*fz*sl**2)/(2.d0*beta**2*brho)
      h(474)=(fzzz*sl**4)/brho
      h(593)=(fz*sl**2)/(2.d0*brho)
      h(627)=(fz*sl**2)/(2.d0*brho)
      h(632)=-((-3.d0 + beta**2)*fz*sl**2)/(2.d0*beta**2*brho)
      h(473)=(sl**5*(-fz**2 + 2.d0*brho*fzzzz))/(2.d0*brho**2)
      h(750)=(fz*sl**2)/(2.d0*brho)
      h(850)=(fz*sl**2)/(2.d0*brho)
      h(855)=-((-3.d0 + beta**2)*fz*sl**2)/(2.d0*beta**2*brho)
      h(623)=-((fzzz*sl**4)/brho)
      h(553)=-(sl**5*(fz**2 + 2.d0*brho*fzzzz))/(2.d0*brho**2)
      h(841)=-((fzzz*sl**4)/brho)
      h(840)=(sl**5*(fz**2 - 2.d0*brho*fzzzz))/(2.d0*brho**2)
c
c add sextupoles
      h(28)=fsxnr*sl2
      h(30)=-3.d0*fsxsk*sl2
      h(39)=-3.d0*fsxnr*sl2
      h(64)=fsxsk*sl2
c
c  add octupoles
c
      h(84)=h(84)+focnr*sl3
      h(86)=-4.0d0*focsk*sl3
      h(95)=-6.d0*focnr*sl3
      h(120)=4.d0*focsk*sl3
      h(175)=h(175)+focnr*sl3
c
      return
      end
c
************************************************************************
c
      subroutine g01234(z,g,gz,gzz,gzzz,gzzzz)
c  This routine computes g, dg/dz, and d/dz(dg/dz) for
c  (a REC Quadrupole Triplet.)
c  Written by Alex Dragt, Fall 1986, and based on work of
c  Rob Ryne and F. Neri
c  This version by T. Mottershead allows up to 9(?) different
c  REC Quadrupoles.
c  Extended to derivatives up to 4 by F. Neri Mar 13 1989. 
      include 'impli.inc'
      include 'quadpn.inc'
c----------------------------------------
      external f,fz,fzz,fzzz,fzzzz
c----------------------------------------
c
      g=0.d0
      gz=0.d0
      gzz=0.d0
      gzzz=0.d0
      gzzzz=0.d0
      kuad = 0.d0
      za = di
      do 60 i=1,ncyc
      do 50 k=1,nqt
      gk = ga(k)
      wk = wd(k)
      dk = dr(k)
      aa = ra(k)
      bb = rb(k)
      do 40 j=1,nr(k)
      zb = za + wk
      g  = g  +gk*(  f(z-zb,aa,bb)-f(z-za,aa,bb))
      gz = gz +gk*(  fz(z-zb,aa,bb)-fz(z-za,aa,bb))
      gzz= gzz+gk*(fzz(z-zb,aa,bb)-fzz(z-za,aa,bb))
      gzzz= gzzz+gk*(fzzz(z-zb,aa,bb)-fzzz(z-za,aa,bb))
      gzzzz= gzzzz+gk*(fzzzz(z-zb,aa,bb)-fzzzz(z-za,aa,bb))
      kuad = kuad + 1
      if(kuad.ge.maxq) return
      za = zb + dk
  40  continue
  50  continue
  60  continue
      return
      end
c
      double precision function f(z,r1,r2)
      implicit double precision (a-h,o-z)
      v(z,r)=1.d0/sqrt(1.d0+(z/r)**2)
      f = 0.5d0-.0625d0*z*
     &(1.d0/r1+1.d0/r2)*v(z,r1)**2*v(z,r2)**2/(v(z,r1)+v(z,r2))*(       &
     &v(z,r1)**2+v(z,r1)*v(z,r2)+v(z,r2)**2+4.d0+8.d0/(v(z,r1)*v(z,r2)))
      return
      end
c
      double precision function fz(z,r1,r2)
      implicit double precision (a-h,o-z)
      v(z,r)=1.d0/sqrt(1.d0+(z/r)**2)
      fz = -.1875d0*(1.d0/r1+1.d0/r2)*v(z,r1)**2*v(z,r2)**2*            &
     &(v(z,r1)**3+v(z,r2)**3 + v(z,r1)**2*v(z,r2)**2/(v(z,r1)+v(z,r2)))
      return
      end
c
      double precision function fzz(z,r1,r2)
      implicit double precision (a-h,o-z)
      v(z,r)=1.d0/sqrt(1.d0+(z/r)**2)
      fzz = .1875d0*(1.d0/r1+1.d0/r2)*z*v(z,r1)**2*v(z,r2)**2 *(        &
     &5.d0*(v(z,r1)**5/r1**2+v(z,r2)**5/r2**2)                          &
     &+v(z,r1)**2*v(z,r2)**2 *(                                         &
     & 2.d0*(v(z,r2)/r1**2+v(z,r1)/r2**2)                               &
     &+4.d0*(v(z,r1)**2/r1**2+v(z,r2)**2/r2**2)/(v(z,r1)+v(z,r2))       &
     &-(v(z,r1)**3/r1**2+v(z,r2)**3/r2**2)/(v(z,r1)+v(z,r2))**2 ))
      return
      end
c
c  Higher derivatives of the Halbach function by
c  F. Neri, Jan 1988.
c
      double precision function fzzz(z,r1,r2)
      implicit double precision (a-h,o-z)
      v(z,r)=1.d0/sqrt(1.d0+(z/r)**2)
      if(dabs(z) .ge. 0.0001d0*r1) then
      fzzz =.375d0*((3.5d0*(v(z,r1)**7/(r1**2*z**2))                    &
     &+3.d0*(v(z,r1)**7/z**4)+7.d0*(v(z,r1)**9/(r1**2*z**2))            &
     &+24.5d0*(v(z,r1)**9/r1**4)                                        &
     &                    -3.5d0*(v(z,r2)**7/(r2**2*z**2))              &
     &-3.d0*(v(z,r2)**7/z**4)-7.d0*(v(z,r2)**9/(r2**2*z**2))            &
     &-24.5d0*(v(z,r2)**9/r2**4))/(1.d0/r1-1.d0/r2))
      else
c--------------------------------------
c  Limiting expression for z ~ 0
      fzzz = 6.d0*(35.d0/128.d0)*(1./(r1*r2**2)+1./(r1**2*r2)           &
     &+ 1.d0/(r1**3)+1.d0/(r2**3))                                      &
     &+60.d0*(-63.d0/256.d0)*z**2*(1.d0/(r1*r2**4)+1.d0/(r1**2*r2**3)   &
     &+1.d0/(r1**3*r2**2)+1.d0/(r1**4*r2)+1.d0/(r1**5)+1.d0/(r2**5))    &
     &+210.d0*(495.d0/2048.d0)*z**4*(1.d0/(r1*r2**6)+1.d0/(r1**2*r2**5) &
     &+1.d0/(r1**3*r2**4)+1.d0/(r1**4*r2**3)+1.d0/(r1**5*r2**2)         &
     &+1.d0/(r1**6*r2)+1.d0/(r1**7)+1.d0/(r2**7))                       &
     &+504.d0*(-1001.d0/4096.d0)*z**6*(1./(r1*r2**8)+1.d0/(r1**2*r2**7) &
     &+1.d0/(r1**3*r2**6)+1.d0/(r1**4*r2**5)+1.d0/(r1**5*r2**4)         &
     &+1.d0/(r1**6*r2**3)+1.d0/(r1**7*r2**2)+1.d0/(r1**8*r2)            &
     &+1.d0/(r1**9)+1.d0/(r2**9))
      endif
      return
      end
c
      double precision function fzzzz(z,r1,r2)
      implicit double precision (a-h,o-z)
      v(z,r)=1.d0/sqrt(1.d0+(z/r)**2)
      if (dabs(z) .ge. 0.0001d0*r1) then
      fzzzz = 0.375d0*((                                                &
     &-220.5d0*((z*v(z,r1)**11)/r1**6)-7.d0*(v(z,r1)**7/(r1**2*z**3))   &
     &-12.d0*(v(z,r1)**7/z**5)-35.d0*(v(z,r1)**9/(r1**2*z**3))          &
     &-24.5d0*(v(z,r1)**9/(r1**4*z))-63.d0*(v(z,r1)**11/(r1**4*z))      &
     &+220.5d0*((z*v(z,r2)**11)/r2**6)+7.d0*(v(z,r2)**7/(r2**2*z**3))   &
     &+12.d0*(v(z,r2)**7/z**5)+35.d0*(v(z,r2)**9/(r2**2*z**3))          &
     &+24.5d0*(v(z,r2)**9/(r2**4*z))+63.d0*(v(z,r2)**11/(r2**4*z))      &
     & )/(1.d0/r1-1.d0/r2))
c----------------------------------------
      else
c  Limiting expression for z ~ 0
cryne 105.*18.?
cryne fzzzz = -105.*18*z*(1./r1**6-1./r2**6)/(64.*(1./r1-1./r2))
      fzzzz = -105.d0*18.d0*z*(1.d0/r1**6-1.d0/r2**6)/(64.d0*(1.d0/r1-  &
     &1.d0/r2))                                                         &
     &+840.d0*z**3*(495.d0/2048.d0)*( 1.d0/(r1*r2**6)                   &
     &+ 1.d0/(r1**2*r2**5)                                              &
     &+ 1.d0/(r1**3*r2**4) + 1.d0/(r1**4*r2**3) + 1.d0/(r1**5*r2**2)    &
     &+ 1.d0/(r1**6*r2) + 1.d0/r1**7 + 1.d0/r2**7 )                     &
     &-3024.d0*z**5*(1001.d0/4096.d0)*(1.d0/(r1*r2**8)                  &
     &+1.d0/(r1**2*r2**7)                                               &
     &+1.d0/(r1**3*r2**6)+1.d0/(r1**4*r2**5)+1.d0/(r1**5*r2**4)         &
     &+1.d0/(r1**6*r2**3)+1.d0/(r1**7*r2**2)+1.d0/(r1**8*r2)            &
     &+1.d0/r1**9 + 1.d0/r2**9 )
c
      endif
      return
      end
c
c end of file
