***********************************************************************
* header              ANALYSIS (advanced commands)                    *
*  Routines for advanced commands and advanced analysis               *
***********************************************************************
c
      subroutine amap(p,fa,fm)
c subroutine for applying a map to a function or a set of moments
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'usrdat.inc'
      character*3 kynd
c
      dimension p(6)
      dimension fa(monoms),fm(6,6)
c
      dimension ga(monoms),gm(6,6)
      dimension ha(monoms),hm(6,6)
      dimension t1a(monoms)
      dimension t2a(monoms)
c
c set up control indices
      mode=nint(p(1))
      isend=nint(p(2))
      ifile=nint(p(3))
      nopt=nint(p(4))
      nskip=nint(p(5))
      nmpo=nint(p(6))
c
c procedure for reading in function or moments when mode = 1, 2, or 3:
      if (mode.eq.1 .or. mode.eq.2 .or. mode.eq.3) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      kynd='gtm'
      call strget(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
c
c procedure for reading in moments when mode = 4 or 5:
      if (mode.eq.4 .or. mode.eq.5) then
c Filippo provide this procedure
      continue
      endif
c
c procedure for letting map act on a function:
      if(mode.eq.1) then
c let map characterized by fa,fm act on ga
c the result is the array ha
c this amounts to computing ha = (Dtranspose)*ga
c where D = D(fa,fm)
      call fxform(fa,fm,ga,ha)
      endif
c
c procedure for letting map act on moments:
c letting the map act on moments amounts to computing ha = D*ga
c
c procedure for mode = 2 or 3:
      if (mode.eq.2 .or. mode.eq.3) then
c
c clear arrays
      do 10 i=1,monoms
      t1a(i)=0.d0
   10 ha(i)=0.d0
c
c when mode = 3, compute all transformed moments through 4'th moments
      imax=209
c when mode = 2, compute transformed moments only through 2'nd moments
      if (mode.eq.2) imax=27
c
c perform calculation
      do 20 i=7,imax
      t1a(i)=1.d0
      call fxform(fa,fm,t1a,t2a)
      t1a(i)=0.d0
      do 30 j=1,209
   30 ha(i)=ha(i)+t2a(j)*ga(j)
   20 continue
      endif
c
c procedure for mode = 4 or 5:
      if (mode.eq.4 .or. mode.eq.5) then
c when mode = 5, compute all transformed moments through 4'th moments
      imax=209
c when mode = 4, compute transformed moments only through 2'nd moments
      if (mode.eq.4) imax=27
c Filippo put in code here for dealing with 6'th order moments:
c clear arrays
c perform calculation
      continue
      endif
c
c if map has been applied to moments, compute squares of rms emittances:
      if (mode.ne.1) then
      xemit2=ha(7)*ha(13)-ha(8)*ha(8)
      yemit2=ha(18)*ha(22)-ha(19)*ha(19)
      if (isend.eq.1.or.isend.eq.3) then
      write (jof,*) 'xemit2=',xemit2
      write (jof,*) 'yemit2=',yemit2
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write (jodf,*) 'xemit2=',xemit2
      write (jodf,*) 'yemit2=',yemit2
      endif
      endif
c
c write out result ha if nmpo > 0
      mpot=mpo
      mpo=nmpo
      if (nmpo.gt.0) call mapout(0,ha,hm)
      mpo=mpot
c
c write results into the array ucalc
c first clear the array
      do 40 i=1,250
   40 ucalc(i)=0.d0
c complete the task
      nuvar=211
      do 50 i=1,209
   50 ucalc(i)=ha(i)
      if (mode.ne.1) then
      ucalc(210)=xemit2
      ucalc(211)=yemit2
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine asni(p)
      use rays
      use lieaparam, only : monoms
c  This subroutine applies powers of script N inverse to phase space data.
c  It does this using analytic formulas.
c
      include 'impli.inc'
c
      character*3 kynd
c  calling array
      dimension p(6)
c
c  local array
      dimension fa(monoms), fr(monoms)
      dimension fm(6,6)
c
c  set up control parameters
      iopt=nint(p(1))
      nmap=nint(p(2))
      nfcf=nint(p(3))
      istart=nint(p(4))
      igroup=nint(p(5))
      nwrite=nint(p(6))
c
c  begin calculation
c
c  get script N from storage
      kynd='gtm'
      call strget(kynd,nmap,fa,fm)
c
c  procedure for a static map
c
      if( iopt.eq.1) then
c  compute linear phase advances and linear time of flight
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
      wt=fm(5,6)
c
c  transform nonlinear part of map to the static resonance basis
      call ctosr(fa,fr)
c
c  begin outer loop over the sets of particles
      nset=nint(float(nrays)/float(igroup))
      do 10 i=1,nset
c  begin inner loop over the particles within a set
      do 20 j=1,igroup
c
c  get phase-space coordinates
      iray=(i-1)*igroup+j
      do 30 k=1,6
   30 zi(k)=zblock(k,iray)
c
c  compute emittances
      ex2=zi(1)**2 + zi(2)**2
      ey2=zi(3)**2 + zi(4)**2
      pt=zi(6)
c
c  compute phase advances and time of flight terms
c  incorporate - sign needed for inverse in the definition of an
      an=-float(istart+(i-1)*nwrite)
c  compute x and y phase advances
      phix = wx - 2.*pt*fr(28) - 2.*pt*pt*fr(84)
     & - 4.*ex2*fr(87) - 2.*ey2*fr(89)
      phiy = wy - 2.*pt*fr(29) - 2.*pt*pt*fr(85)
     & - 4.*ey2*fr(88) - 2.*ex2*fr(89)
c  compute time-like drift terms
      drt = pt*wt - ex2*fr(28) - ey2*fr(29)
     & - 2.*pt*ex2*fr(84) -2.*pt*ey2*fr(85)
     & - 3.*pt*pt*fr(30) - 4.*pt*pt*pt*fr(86)
c
c  set up matrix quantities
      cx=cos(an*phix)
      sx=sin(an*phix)
      cy=cos(an*phiy)
      sy=sin(an*phiy)
      tof=an*drt
c
c  apply matrix to transverse coordinates
      zf(1)= cx*zi(1)+sx*zi(2)
      zf(2)=-sx*zi(1)+cx*zi(2)
      zf(3)= cy*zi(3)+sy*zi(4)
      zf(4)=-sy*zi(3)+cy*zi(4)
c  transform time deviation and energy deviation
      zf(5)= zi(5)+tof
      zf(6)= zi(6)
c
c  write out results
      write (nfcf,100) zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
  100 format(6(1x,1pe12.5))
c
   20 continue
   10 continue
c
      endif
c
c  procedure for a dynamic map
c
      if( iopt.eq.2) then
c  compute linear phase advances
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
      cwt=fm(5,5)
      swt=fm(5,6)
      wt=atan2(swt,cwt)
c
c  transform nonlinear part of map to the dynamic resonance basis
      call ctodr(fa,fr)
c
c  begin outer loop over the sets of particles
      nset=nint(float(nrays)/float(igroup))
      do 40 i=1,nset
c  begin inner loop over the particles within a set
      do 50 j=1,igroup
c
c  get phase-space coordinates
      iray=(i-1)*igroup+j
      do 60 k=1,6
   60 zi(k)=zblock(k,iray)
c
c  compute emittances
      ex2=zi(1)**2 + zi(2)**2
      ey2=zi(3)**2 + zi(4)**2
      et2=zi(5)**2 + zi(6)**2
c
c  compute phase advances
c  incorporate - sign needed for inverse in the definition of an
      an=-float(istart+(i-1)*nwrite)
c  this part of code not yet complete
      phix=wx
      phiy=wy
      phit=wt
c
c  set up matrix quantities
      cx=cos(an*phix)
      sx=sin(an*phix)
      cy=cos(an*phiy)
      sy=sin(an*phiy)
      ct=cos(an*phit)
      st=sin(an*phit)
c
c  apply matrix to coordinates
      zf(1)= cx*zi(1)+sx*zi(2)
      zf(2)=-sx*zi(1)+cx*zi(2)
      zf(3)= cy*zi(3)+sy*zi(4)
      zf(4)=-sy*zi(3)+cy*zi(4)
      zf(5)= ct*zi(5)+st*zi(6)
      zf(6)=-st*zi(5)+ct*zi(6)
c
c  write out results
      write (nfcf,100) zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
c
   50 continue
   40 continue
c
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine betmap(ana,anm,ba,bm)
c This is a subroutine for finding the betatron portion of a map.
c The map ana, anm is assumed to be a map about the fixed point.
c Written by Alex Dragt, 4 August 1986
      use lieaparam, only : monoms
      include 'impli.inc'
c
      dimension ana(monoms),anm(6,6)
      dimension ba(monoms),bm(6,6)
c
      dimension tempa(monoms),tempm(6,6)
      dimension ca(monoms),cm(6,6)
c
c Extraction of betatron term b.
c First pass: terms linear in pt.
      call clear(tempa,tempm)
      call matmat(anm,tempm)
      tempa(33)=ana(33)
      tempa(38)=ana(38)
      tempa(42)=ana(42)
      tempa(45)=ana(45)
      tempa(53)=ana(53)
      tempa(57)=ana(57)
      tempa(60)=ana(60)
      tempa(67)=ana(67)
      tempa(70)=ana(70)
      tempa(76)=ana(76)
      call mapmap(tempa,tempm,ba,bm)
c Second pass: terms quadratic in pt.
      call inv(ba,bm)
      call concat(ba,bm,ana,anm,ca,cm)
      tempa(104)=ca(104)
      tempa(119)=ca(119)
      tempa(129)=ca(129)
      tempa(135)=ca(135)
      tempa(154)=ca(154)
      tempa(164)=ca(164)
      tempa(170)=ca(170)
      tempa(184)=ca(184)
      tempa(190)=ca(190)
      tempa(200)=ca(200)
      call mapmap(tempa,tempm,ba,bm)
      return
      end
c
***********************************************************************
c
      subroutine cex(p,ga,gm,myorder)
c this routine computes t=exp(:f)
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'extalk.inc'
      include 'hmflag.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension p(6),ga(monoms),gm(6,6)
c
      dimension ta(monoms)
      dimension tm(6,6),em(6,6),fm(6,6)
cryne 3AM 7/11/2002      dimension y(224)
      dimension y(monoms+15)

c--------------------------------------------
c     write(6,*)'inside routine cex'
c     write(6,*)'input 6-vector p:'
c     do i=1,6
c     write(6,*)p(i)
c     enddo
c     write(6,*)'input matrix:'
c     do i=1,6
c     do j=1,6
c     if(gm(i,j).ne.0.d0)write(69,*)i,j,gm(i,j)
c     enddo
c     enddo
c
c     write(6,*)'input polynomial:'
c     do i=1,209
c     if(ga(i).ne.0.d0)write(69,*)i,ga(i)
c     enddo
c     if(p(1).ne.12345)stop
c--------------------------------------------
c
c set up power and control indices
      power=p(1)
      nmapf=nint(p(2))
      nmapg=nint(p(3))
c
c get map and clear arrays
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
      call clear(ta,tm)
c
c perform calculation
c
c set up exponent
      call csmul(power,fa,fa)
c
c compute a scaling factor to bring exponent within range of a
c taylor expansion or GENMAP
      call matify(em,fa)
      call mnorm(em,res)
      kmax=1
      scale=.5d0
   10 continue
      test=res*scale
      if (test.lt..1d0) goto 20
      kmax=kmax+1
      scale=scale/2.d0
      go to 10
   20 continue
c
c select procedure
      itest=1
      do 30 i=28,monoms
      if (fa(i).ne.0.d0) itest=2
      if (itest.ne.1) go to 40
   30 continue
   40 continue
      if (itest.eq.1) go to 50
      if (itest.eq.2) go to 80
c
c procedure using taylor series
   50 continue
      write (12,*) 'exp(:f:) computed using taylor series'
c rescale em
      call smmult(scale,em,em)
c compute taylor series result fm=exp(scale*em)
      call exptay(em,fm)
c raise the result to the 2**kmax (= 1/scale) power
      do 60 i=1,kmax
      call mmult(fm,fm,tm)
      call matmat(tm,fm)
   60 continue
      goto 200
c
c procedure using genmap
   80 continue
      write(12,*) 'exp(:f:) computed using GENMAP'
c rescale fa
      call csmul(scale,fa,fa)
c setup and initialize for GENMAP routines
      iflag=1
      t=0.d0
      ns=50.d0
      h=.02d0
cryne 3AM 7/11/2002      ne=224
cryne 1 August 2004      ne=monoms+15
cryne 1 August 2004 Initialize:
cryne do 90 i=1,ne
      do 90 i=1,monoms+15
   90 y(i)=0.d0
      do 100 i=1,6
      j=7*i
  100 y(j)=1.d0
c call GENMAP routines
cryne there is a better way to do this; fix later.
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c  y(225-476) = f5
c  y(477-938) = f6
      if(myorder.eq.1)ne=42    !36+6
      if(myorder.eq.2)ne=98    !83+15
      if(myorder.eq.3)ne=224    !209+15
      if(myorder.eq.4)ne=476    !461+15
      if(myorder.eq.5)ne=938    !923+15
c
      call adam11(h,ns,'start',t,y,ne)
      call putmap(y,fa,fm)
c
c raise the result to the 2**kmax (= 1/scale) power
      do 110 i=1,kmax
      call concat(fa,fm,fa,fm,ta,tm)
      call mapmap(ta,tm,fa,fm)
  110 continue
      go to 200
c
c decide where to put results
c
  200 continue
      if (nmapg.ge.1 .and. nmapg.le.5) then 
      kynd='stm'
      call strget(kynd,nmapg,ta,tm)      
      endif
c
      if (nmapg.eq.0) call mapmap(ta,tm,ga,gm)
c 
      if (nmapg.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
*******************************************************************************
c
      subroutine chrexp(iopt,delta,ta,tm,am1,am2,am3)
c
c This subroutine computes the chromatic expansion of the map ta,tm
c for the case in which the f3 and f4 parts of ta contain only terms
c linear and quadratic in pt, respectively.
c Written by Alex Dragt, Spring 1987
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
      dimension ta(monoms),tm(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
      dimension t1a(monoms),t1m(6,6),t2m(6,6)
c--------
c procedure when IOPT = 1 (delta=pt).
      if (iopt.eq.1) then
c Calculation of matrix associated with pt*f2 terms in ta.
c Set up f2a in t1a.
      call clear(t1a,am1)
      t1a(7)=ta(33)
      t1a(8)=ta(38)
      t1a(9)=ta(42)
      t1a(10)=ta(45)
      t1a(13)=ta(53)
      t1a(14)=ta(57)
      t1a(15)=ta(60)
      t1a(18)=ta(67)
      t1a(19)=ta(70)
      t1a(22)=ta(76)
c Compute matrix t1m corresponding to :f2a:.
      call matify(t1m,t1a)
c      write(6,*) 'result from chrexp'
c      call pcmap(1,0,0,0,t1a,t1m)
c Compute am1=t1m*tm
      call mmult(t1m,tm,am1)
c Calculation of matrix associated with (pt**2)*f2 terms in ta.
c Set up f2a in t1a.
      call clear(t1a,am2)
      t1a(7)=ta(104)
      t1a(8)=ta(119)
      t1a(9)=ta(129)
      t1a(10)=ta(135)
      t1a(13)=ta(154)
      t1a(14)=ta(164)
      t1a(15)=ta(170)
      t1a(18)=ta(184)
      t1a(19)=ta(190)
      t1a(22)=ta(200)
c Compute matrix t2m corresponding to :f2a:.
      call matify(t2m,t1a)
c      write(6,*) 'result from chrexp'
c      call pcmap(1,0,0,0,t1a,t2m)
c Compute am3=t2m+t1m*t1m/2.
      call mmult(t1m,t1m,am3)
      call smmult(.5d0,am3,am3)
      call madd(t2m,am3,am3)
c Compute am2=(t2m+t1m*t1m/2.)*tm
      call mmult(am3,tm,am2)
c Compute am3=tm+delta*am1+delta2*am2 for specific value of delta.
      delta2=delta*delta
      call matmat(tm,am3)
      call smmult(delta,am1,t1m)
      call madd(am3,t1m,am3)
      call smmult(delta2,am2,t1m)
      call madd(am3,t1m,am3)
      endif
c
c procedure when IOPT = 2 (delta=dp/p0).
      if (iopt.eq.2) then
      continue
      endif
c
      return
      end
c
      subroutine cod(p,th,tmh)
c This is a subroutine for computing closed orbit data.
c Written by Alex Dragt, 6 November 1985
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6)
      dimension th(monoms),tmh(6,6)
c
      dimension ana(monoms),anm(6,6)
      dimension ta(monoms),tm(6,6)
      dimension tempa(monoms),tempm(6,6)
      dimension ba(monoms),bm(6,6)
      dimension ca(monoms),cm(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
c Set up control indices:
      iopt=nint(p(1))
      delta=p(2)
      idata=nint(p(3))
      ipmaps=nint(p(4))
      isend=nint(p(5))
      iwmaps=nint(p(6))
c
c Write headings
      if (isend.eq.1 .or. isend.eq.3) write(jof,90)
      if (isend.eq.2 .or. isend.eq.3) write(jodf,90)
   90 format(/,1x,'closed orbit analysis for static map')
c
c Computation of fixed point and map about it.
      call fxpt(th,tmh,ana,anm,ta,tm)
c
c Procedure for output of closed orbit location.
      if(idata.eq.1 .or. idata.eq.3) then
c Procedure when IOPT = 1:
      if (iopt.eq.1) then
c Compute location of closed orbit for given delta value
      delta2=delta*delta
      delta3=delta*delta2
      xc=delta*tm(1,6)-delta2*ta(63)-delta3*ta(174)
      pxc=delta*tm(2,6)+delta2*ta(48)+delta3*ta(139)
      yc=delta*tm(3,6)-delta2*ta(79)-delta3*ta(204)
      pyc=delta*tm(4,6)+delta2*ta(73)+delta3*ta(194)
c Write out results
      do 5 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 5
      write(ifile,100)
  100 format(/,1x,'closed orbit data for delta defined in terms of',
     & 1x,'energy deviation:')
      write(ifile,120)
  120 format(1x,'location of closed orbit (x,px,y,py)')
      write(ifile,130)
  130 format(/,1x,'terms linear in delta')
      write(ifile,140) tm(1,6),tm(2,6),tm(3,6),tm(4,6)
  140 format(1x,4(d15.8,2x))
      write(ifile,150)
  150 format(/,1x,'terms quadratic in delta')
      write(ifile,140) -ta(63),ta(48),-ta(79),ta(73)
      write(ifile,160)
  160 format(/,1x,'terms cubic in delta')
      write(ifile,140) -ta(174),ta(139),-ta(204),ta(194)
      write(ifile,110) delta
  110 format(/,1x,'location of closed orbit when delta = ',d15.8)
      write(ifile,140) xc,pxc,yc,pyc
    5 continue
      endif
c Procedure when IOPT = 2
      if (iopt.eq.2) then
c Compute location of closed orbit for given delta value
c
c the code below needs to be modified
      delta2=delta*delta
      delta3=delta*delta2
      xc=delta*tm(1,6)-delta2*ta(63)-delta3*ta(174)
      pxc=delta*tm(2,6)+delta2*ta(48)+delta3*ta(139)
      yc=delta*tm(3,6)-delta2*ta(79)-delta3*ta(204)
      pyc=delta*tm(4,6)+delta2*ta(73)+delta3*ta(194)
c end of code to be modified
c
c Write out results
      do 7 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 7
      write(ifile,102)
  102 format(/,1x,'closed orbit data for delta defined in terms of',
     & 1x,'momentum deviation:')
c
c the code below needs to be modified
      write(ifile,120)
      write(ifile,130)
      write(ifile,140) tm(1,6),tm(2,6),tm(3,6),tm(4,6)
      write(ifile,150)
      write(ifile,140) -ta(63),ta(48),-ta(79),ta(73)
      write(ifile,160)
      write(ifile,140) -ta(174),ta(139),-ta(204),ta(194)
      write(ifile,110) delta
      write(ifile,140) xc,pxc,yc,pyc
c end of code to be modified
c
    7 continue
      write(6,*) 'IOPT = 2 case not yet installed completely'
      endif
      endif
c
c Factorization of map into betatron and remaining terms.
c Computation of betatron term script B.
      call betmap(ana,anm,ba,bm)
c Computation of remaining nonlinear correction map script C.
      call mapmap(ba,bm,tempa,tempm)
      call inv(tempa,tempm)
      call concat(tempa,tempm,ana,anm,ca,cm)
c Computation of B for specific value of delta.
      call chrexp(iopt,delta,ba,bm,am1,am2,am3)
c
c Procedure for output of twiss matrix and corrections.
      if(idata.eq.2 .or. idata.eq.3) then
c Procedure when IOPT =1
      if (iopt.eq.1) then
c Print out matrices bm, am1, and am2.
      do 9 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 9
      write(ifile,198)
  198 format(//,1x,'twiss matrix for delta defined in terms of energy',
     &1x,'deviation:')
      write(ifile,200)
  200 format(/,1x,'on energy twiss matrix')
      call pcmap(i,0,0,0,ba,bm)
      write(ifile,300)
  300 format(//,1x,'matrix for delta correction')
      call pcmap(i,0,0,0,ba,am1)
      write(ifile,400)
  400 format(//,1x,'matrix for delta**2 correction')
      call pcmap(i,0,0,0,ba,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
  402 format(//,1x,'value of twiss matrix when delta= ',d15.8)
      call pcmap(i,0,0,0,ba,am3)
    9 continue
      endif
c Procedure when IOPT = 2
      if (iopt.eq.2) then
c Print out matrices bm, am1, and am2.
      do 11 i=1,2
      if (i.eq.1) then
      iflag=1
      ifile=jof
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 11
      write(ifile,199)
  199 format(//,1x,'twiss matrix for delta defined in terms of',
     &1x,'momentum deviation:')
      write(ifile,202)
  202 format(/,1x,'on momentum twiss matrix')
      call pcmap(i,0,0,0,ba,bm)
      write(ifile,300)
      call pcmap(i,0,0,0,ba,am1)
      write(ifile,400)
      call pcmap(i,0,0,0,ba,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
      call pcmap(i,0,0,0,ba,am3)
   11 continue
      endif
      endif
c
c Procedure for output of the maps BC, B, and C.
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      do 13 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof 
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 13
      write(ifile,600)
  600 format(//,1x,'total transfer map about the closed orbit')
      call pcmap(i,i,0,0,ana,anm)
      write(ifile,700)
  700 format(//,1x,'betatron factor of transfer map')
      call pcmap(i,i,0,0,ba,bm)
      write(ifile,800)
  800 format(//,1x,'nonlinear factor of transfer map')
      call pcmap(i,i,0,0,ca,cm)
   13 continue
      endif
c
c Procedure for output of script T.
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      do 15 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof 
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 15
      write(ifile,900)
  900 format(//,1x,'transfer map script T to the closed orbit')
      call pcmap(i,i,0,0,ta,tm)
   15 continue
      endif
 
c
c Put maps in buffers
      call mapmap(ana,anm,buf1a,buf1m)
      call mapmap(ba,bm,buf2a,buf2m)
      call mapmap(ca,cm,buf3a,buf3m)
      call clear(buf4a,buf4m)
      call mapmap(buf4a,am3,buf4a,buf4m)
      call mapmap(ta,tm,buf5a,buf5m)
c
c Procedure for writing of maps.
      if(iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ana,anm)
      call mapout(0,ba,bm)
      call mapout(0,ca,cm)
      call mapout(0,buf4a,buf4m)
      call mapout(0,ta,tm)
      mpo=mpot
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine csym(isend,fm,ans)
c
c This subroutine checks the symplectic condition for the matrix fm
c Written by Alex Dragt, 4 October 1989
c
      include 'impli.inc'
      include 'files.inc'
c
c Calling arrays
      dimension fm(6,6)
c
c Temporary arrays
      dimension tm(6,6)
c
c-----Procedure
c
      call matmat(fm,tm)
      call minv(tm)
      call mmult(tm,fm,tm)
      do 10 i=1,6
   10 tm(i,i)=tm(i,i)-1.d0
      call mnorm(tm,ans)
c
c Write out results if desired
c
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' symplectic violation = ',ans
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' symplectic violation = ',ans
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine dia(p,fa,fm)
c this is a routine for computing invariants in the dynamic case
c Written by Alex Dragt, Spring 1987
      use rays
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c
c set up control indices
      iopt=nint(p(1))
      idata=nint(p(2))
      ipmaps=nint(p(3))
      isend=nint(p(4))
      iwmaps=nint(p(5))
      iwnum=nint(p(6))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'dynamic invariant analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'dynamic invariant analysis'
      endif
c
c begin calculation
c
c find the transforming (conjugating) map script A
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm)
c remove offensive terms from f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c remove offensive terms from f4 part of map:
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c put script A in buffer 1 and script N in buffer 2
      call mapmap(ta,tm,buf1a,buf1m)
      call mapmap(ga,gm,buf2a,buf2m)
c
c procedure for computing invariant
c invert script A
      call inv(ta,tm)
      call clear(t1a,t1m)
c computation of x invariant
      if(iopt.eq.1) then
      t1a(7)=1.d0
      t1a(13)=1.d0
      endif
c computation of y invariant
      if(iopt.eq.2) then
      t1a(18)=1.d0
      t1a(22)=1.d0
      endif
c computation of t invariant
      if(iopt.eq.3) then
      t1a(25)=1.d0
      t1a(27)=1.d0
      endif
c computation of mixed invariant
      if(iopt.ge.4) then
      read(iopt,*) anx,any,ant
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      t1a(7)=anx
      t1a(13)=anx
      t1a(18)=any
      t1a(22)=any
      t1a(25)=ant
      t1a(27)=ant
      endif
c put invariant in buffer 3
      call ident(buf3a,buf3m)
      call fxform(ta,tm,t1a,buf3a)
c
c procedure for putting out data
      if (idata.eq.1 .or. idata.eq.3) then
      do 10 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 10
      write(ifile,*)
      write(ifile,*) 'invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
   10 continue
      endif
      if (idata.eq.2 .or. idata.eq.3) then
      mpot=mpo
      mpo=18
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c procedure for printing out maps
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form map script N'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   20 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c procedure for computing values of invariant and writing them out
      if (iwnum.gt.0) then
      write(jof,*)
      write(jof,*) 'values of invariant written on file ',iwnum
      do 30 k=1,nraysp
      do 40 j=1,6
   40 zi(j)=zblock(j,k)
      call evalf(zi,buf3a,val2,val3,val4)
      write(iwnum,60) k,val2,val3,val4,0.,0.
   60 format(1x,i12,5(1x,1pe12.5))
   30 continue
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the invariant polynomial.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m) 
c
      return
      end
c
***********************************************************************
c
      subroutine dnor_old(p,fa,fm)
c this is a subroutine for normal form analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c
c set up control indices
      idata=nint(p(1))
      ipmaps=nint(p(2))
      isend=nint(p(3))
      iwmaps=nint(p(4))
c
c write headings
      if ((isend.eq.1 .or. isend.eq.3) .and. idproc.eq.0) then
      write (jof,*)
      write (jof,*) 'dynamic normal form analysis'
      endif
      if ((isend.eq.2 .or. isend.eq.3) .and. idproc.eq.0) then
      write (jodf,*)
      write (jodf,*) 'dynamic normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm,t2m)
c remove offensive terms from f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c remove offensive terms from f4 part of map:
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(ga,gm,buf2a,buf2m)  
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=gm(1,1)
      swx=gm(1,2)
      wx=atan2(swx,cwx)
      cwy=gm(3,3)
      swy=gm(3,4)
      wy=atan2(swy,cwy)
      cwt=gm(5,5)
      swt=gm(5,6)
      wt=atan2(swt,cwt)
c set up normal form for exponent
      do 10 i=1,27
   10 ga(i)=0.
      ga(7)=-wx/2.d0
      ga(13)=-wx/2.d0
      ga(18)=-wy/2.d0
      ga(22)=-wy/2.d0
      ga(25)=-wt/2.d0
      ga(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or. idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,ga,g1a)
c store results in buffer 3
      call mapmap(g1a,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,ga,gm)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)
     &write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m) 
c
      return
      end
c
*******************************************************************
c
      subroutine dnor(p,fa,fm)
c this is a subroutine for normal form analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
c
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6)
c
c set up control indices
      keep=  nint(p(1))
      idata= nint(p(2))
      ipmaps=nint(p(3))
      isend= nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      if(idproc.eq.0)write (jof,*)
      if(idproc.eq.0)write (jof,*) 'dynamic normal form analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      if(idproc.eq.0)write (jodf,*)
      if(idproc.eq.0)write (jodf,*) 'dynamic normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm,t1m)
c remove offensive terms from f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(ga,gm,buf2a,buf2m)
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=gm(1,1)
      swx=gm(1,2)
      wx=atan2(swx,cwx)
      cwy=gm(3,3)
      swy=gm(3,4)
      wy=atan2(swy,cwy)
      cwt=gm(5,5)
      swt=gm(5,6)
      wt=atan2(swt,cwt)
c set up normal form for exponent
      do 10 i=1,27
   10 ga(i)=0.
      ga(7)=-wx/2.d0
      ga(13)=-wx/2.d0
      ga(18)=-wy/2.d0
      ga(22)=-wy/2.d0
      ga(25)=-wt/2.d0
      ga(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or. idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,ga,g1a)
c store results in buffer 3
      call mapmap(g1a,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,ga,gm)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)
     &  write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
*******************************************************************
c
      subroutine fadm(p,fa,fm)
c  this subroutine fourier analyzes a dynamic transfer map
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'fadm not yet available'
      return
      end
c
*******************************************************************
c
      subroutine fasm(p,fa,fm)
c  this subroutine fourier analyzes a static transfer map
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'fasm not yet available'
      return
      end
c
*******************************************************************
c
      subroutine gbuf_old(p,fa,fm)
c  this subroutine gets a map from an auxiliary buffer and
c  concatenates it with the map in the main buffer.
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ta(monoms),tm(6,6)
c
c set up and test control index
      i=nint(p(1))
      if (i.lt.1 .or. i.gt.5) then
      write(6,*) 'trouble with index nmap in gbuf'
      return
      endif
c
c concatenate the map in bufi with the  existing map
      if(i.eq.1) call scncat(fa,fm,buf1a,buf1m,ta,tm)
      if(i.eq.2) call scncat(fa,fm,buf2a,buf2m,ta,tm)
      if(i.eq.3) call scncat(fa,fm,buf3a,buf3m,ta,tm)
      if(i.eq.4) call scncat(fa,fm,buf4a,buf4m,ta,tm)
      if(i.eq.5) call scncat(fa,fm,buf5a,buf5m,ta,tm)
      call mapmap(ta,tm,fa,fm)
c
      return
      end
c
c*****************************************************************************
c
      subroutine gbuf(p,fa,fm)
c  this subroutine gets a map from an auxiliary buffer and
c  either concatenates it with the map in the main buffer,
c  or uses it to replace the map in the main buffer.
c  Written by Alex Dragt, Spring 1987
c  Modified by Alex Dragt, 17 June 1988
c  Modified by Alex Dragt, 13 October 1988
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c  Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c set up and test control indices
      iopt=nint(p(1))
      i=nint(p(2))
      if (i.lt.1 .or. i.gt.5) then
      write(6,*) 'trouble with index nmap in gbuf'
      return
      endif
c
c if iopt=1, concatenate the existing map with the map in bufi
      if(iopt.eq.1) then
      if(i.eq.1) call concat(fa,fm,buf1a,buf1m,fa,fm)
      if(i.eq.2) call concat(fa,fm,buf2a,buf2m,fa,fm)
      if(i.eq.3) call concat(fa,fm,buf3a,buf3m,fa,fm)
      if(i.eq.4) call concat(fa,fm,buf4a,buf4m,fa,fm)
      if(i.eq.5) call concat(fa,fm,buf5a,buf5m,fa,fm)
      return
      endif
c
c if iopt=2, replace the existing map with the map in bufi
      if(iopt.eq.2) then
      if(i.eq.1) call mapmap(buf1a,buf1m,fa,fm)
      if(i.eq.2) call mapmap(buf2a,buf2m,fa,fm)
      if(i.eq.3) call mapmap(buf3a,buf3m,fa,fm)
      if(i.eq.4) call mapmap(buf4a,buf4m,fa,fm)
      if(i.eq.5) call mapmap(buf5a,buf5m,fa,fm)
      return
      endif
c
      write(6,*) 'trouble with index iopt in gbuf'
      return
      end
c
c*****************************************************************************
c
      subroutine geom(pp)
c
c   Routine to compute geometry of a loop.
c   Written by A. Dragt 8/27/92.
c   Modified 5/27/98 AJD.
c   Modified 1/8/99 AJD.
c   Based on the subroutines cqlate and pmif
c
      use beamdata
      use lieaparam
      use acceldata
      include 'impli.inc'
c
c common blocks
c
      include 'codes.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
      include 'pie.inc'
      include 'parset.inc'
      include 'usrdat.inc'
c
      dimension pp(6)
c
c local variables
c
      character*8 string(5),str
      dimension ex(3), ey(3), ez(3)
      dimension exl(3), eyl(3), ezl(3)
      dimension exlt(3), eylt(3), ezlt(3)
      dimension dims(6)
c
c  set up control indices
c
      iopt=nint(pp(1))
      si=pp(2)
      ti=pp(3)
      ipset1=nint(pp(4))
      ipset2=nint(pp(5))
      isend=nint(pp(6))
c
c get contents of psets
c
      if (ipset1.gt.0 .and. ipset1.le.maxpst) then
      xi=pst(1,ipset1)
      yi=pst(2,ipset1)
      zi=pst(3,ipset1)
      phid=pst(4,ipset1)
      thetad=pst(5,ipset1)
      psid=pst(6,ipset1)
      phir=phid*pi180
      thetar=thetad*pi180
      psir=psid*pi180
      endif
      if (ipset2.gt.0 .and. ipset2.le.maxpst) then
      ifile=nint(pst(1,ipset2))
      jfile=nint(pst(2,ipset2))
      kfile=nint(pst(3,ipset2))
      lfile=nint(pst(4,ipset2))
      mfile=nint(pst(5,ipset2))
      ndpts=nint(pst(6,ipset2))
      endif
c
c  set up variables and constants
c
      ss=si
      tt=ti
      vr=1.d0/(beta*c)
      x=xi
      y=yi
      z=zi
      do i=1,6
      dims(i)=0.d0
      end do
      do i=1,3
      ex(i)=0.d0
      ey(i)=0.d0
      ez(i)=0.d0
      end do
      ex(1)=1.d0
      ey(2)=1.d0
      ez(3)=1.d0
c      write(6,*) ' thetar =', thetar
      call erotv(phir,thetar,psir,ex,exl)
      call erotv(phir,thetar,psir,ey,eyl)
      call erotv(phir,thetar,psir,ez,ezl)
c
c  start routine
c
c write out header
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) '  '
      write(jof,*) ' Geometrical Analysis'
      write(jof,*) '  '
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) '  '
      write(jodf,*) ' Geometrical Analysis'
      write(jodf,*) '  '
      endif
c
c  see if a loop exists
      if(nloop.le.0) then
      write(jof ,*) ' error from geom: no loop has been specified'
      write(jodf,*) ' error from geom: no loop has been specified'
      return
      endif
c
c write out starting values
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' initial s,x,y,z,t;ex,ey,ez:'
      write(jof,237) si,xi,yi,zi,ti
      write(jof,*) exl(1),exl(2),exl(3)
      write(jof,*) eyl(1),eyl(2),eyl(3)
      write(jof,*) ezl(1),ezl(2),ezl(3)
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jof,*) ' initial s,x,y,z,t;ex,ey,ez:'
      write(jodf,237) si,xi,yi,zi,ti
      write(jodf,*) exl(1),exl(2),exl(3)
      write(jodf,*) eyl(1),eyl(2),eyl(3)
      write(jodf,*) ezl(1),ezl(2),ezl(3)
      endif
c
c scan the loop
c
      do 137 jk1=1,joy
c
c initialize various quantities
        icat=0
        grad = 0.d0
c
c record element category and set various quantities
c
c element
        if(mim(jk1).lt.0) then
          string(1)=lmnlbl(-mim(jk1))(1:8)
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(1)=lmnlbl(mim(jk1)-5000)(1:8)
c lump
        else
          string(1)=ilbl(inuse(mim(jk1)))(1:8)
        endif
      call lookup(string(1),itype,item)
c      write(6,513) string(1)
c  513 format(1x,a8)
c      write(6,*) 'itype and item are ',itype, item
c
c procedure for a menu item
c
      if(itype.eq.1) then
      k=item
      imax=nrp(nt1(k),nt2(k))
c
c see if item is a simple command
c
      if (nt1(k) .eq. 7) then
c dims
      if (nt2(k) .eq. 35) then
      do ii=1,6
      dims(ii)=pmenu(ii+mpp(k))
      end do
      endif
c
      endif
c
c see if item is a simple element
c
      if (nt1(k) .eq. 1) then
c
c drift
      if (nt2(k) .eq. 1) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c spce
      if (nt2(k) .eq. 25) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c arc
      if (nt2(k) .eq. 30) then
      icat=4
      d1 = pmenu(1+mpp(k))
      d2 = pmenu(2+mpp(k))
      aleng = pmenu(3+mpp(k))
      angd = pmenu(4+mpp(k))
      endif
c quad
      if (nt2(k) .eq. 9) then
      icat=1
      aleng = pmenu(1+mpp(k))
      grad = pmenu(2+mpp(k))
      endif
c cfqd
      if (nt2(k) .eq. 18) then
      icat=1
      aleng = pmenu(1+mpp(k))
      ipst = nint(pmenu(2+mpp(k)))
      if (ipst.gt.0 .and. ipst.le.maxpst) then
      grad = pst(1,ipst)
      endif
      endif
c sext
      if (nt2(k) .eq. 10) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c octm
      if (nt2(k) .eq. 11) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c recm
      if (nt2(k) .eq. 24) then
      icat=1
      aleng = pmenu(2+mpp(k)) - pmenu(1+mpp(k))
c
c computing the gradient for recm is complicated
c and so it is not done for the time being
c
      endif
c sol
      if (nt2(k) .eq. 20) then
      icat=1
      aleng = pmenu(2+mpp(k)) - pmenu(1+mpp(k))
      endif
c nbend
      if (nt2(k) .eq. 2) then
      icat=2
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c pbend
      if (nt2(k) .eq. 3) then
      icat=2
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c gbend
      if (nt2(k) .eq. 4) then
      icat=2
      angd=pmenu(1+mpp(k))
      b=pmenu(6+mpp(k))
      endif
c gbdy
      if (nt2(k) .eq. 6) then
      icat=2
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c cfbd
      if (nt2(k) .eq. 8) then
      icat=2
      angd=pmenu(1+mpp(k))
      b=pmenu(2+mpp(k))
      ipst = nint(pmenu(6+mpp(k)))
      if (ipst.gt.0 .and. ipst.le.maxpst) then
      grad = pst(1,ipst)
      endif
      endif
c arot
      if (nt2(k) .eq. 14) then
      icat=3
      angd=pmenu(1+mpp(k))
      endif
c prot
c      if (nt2(k) .eq. 5) then
c      angd=pmenu(1+mpp(k))
c      kind=nint(pmenu(2+mpp(k)))
c      endif
c
c carry out computations for items dims thru prot above
c
c compute and write out local geometry
      if( iopt .eq. 2 .or. iopt .eq. 3
     &  .or. iopt .eq. 12 .or. iopt .eq. 13) then
c
c procedure for a straight element
      if(icat .eq. 1) then
c
c write all descriptive information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
  605 format(1h ,1x,a8,1x,a8,1x,i5,1x,i5,1x,i5)
         if(imax.eq.0)goto 137
      write(mfile,607)(pmenu(i+mpp(k)),i=1,imax)
c Output using Mottershead's favorite pg format
  607 format((1h ,3(1x,1pg22.15)))
      write(mfile,611)  lmnlbl(k), aleng
  611 format(2x,a,f9.4)
      write(mfile,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
  608 format(1x,6f12.4)
      endif
c
c write all descriptive information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jof,607)(pmenu(i+mpp(k)),i=1,imax)
c Output using Mottershead's favorite pg format
      write(jof,609)  lmnlbl(k), aleng
  609 format(2x,a,'length =',f9.4,2x,'dimensions:')
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jof,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jodf,607)(pmenu(i+mpp(k)),i=1,imax)
      write(jodf,609)  lmnlbl(k), aleng
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jodf,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
c
c subdivide element
      ndptst=ndpts
      if(ndpts .lt. 1) ndptst=1
      daleng=aleng/float(ndptst)
      if( dabs(aleng) .lt. 1.0d-10) ndptst=0
      do i=0,ndptst
      sst = ss + daleng*float(i)
      ttt = tt + daleng*float(i)*vr
      xt = x + daleng*float(i)*ezl(1)
      yt = y + daleng*float(i)*ezl(2)
      zt = z + daleng*float(i)*ezl(3)
c
c write simple element descriptive information 
c and pathlength information to external file kfile
      if (kfile .gt. 0) then
      write(kfile,610) lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad
  610 format (1x,a8,2(i3),4(1x,1pe12.5))
      endif
c
c write simple element descriptive information and pathlength
c and coordinate information to external file lfile
      if (lfile .gt. 0) then
      write(lfile,710)
     & lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad,xt,yt,zt,ttt
  710 format (1x,a8,2(i3),8(1x,1pe12.5))
      endif
c
c write coordinate information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,236) sst,xt,yt,zt,ttt,0.0
  236 format(6(1x,1pe12.5))
      endif
c
c write coordinate information to terminal and/or drop file     
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,237) sst,xt,yt,zt,ttt
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,237) sst,xt,yt,zt,ttt
      endif
  237 format(5(1x,1pe12.5))
c
c write orientation information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,*) exl(1), exl(2), exl(3)
      write(mfile,*) eyl(1), eyl(2), eyl(3)
      write(mfile,*) ezl(1), ezl(2), ezl(3)
      endif
c
c write orientation information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) exl(1), exl(2), exl(3)
      write(jof,*) eyl(1), eyl(2), eyl(3)
      write(jof,*) ezl(1), ezl(2), ezl(3)
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) exl(1), exl(2), exl(3)
      write(jodf,*) eyl(1), eyl(2), eyl(3)
      write(jodf,*) ezl(1), ezl(2), ezl(3)
      endif
c
      end do
      endif
c
c procedure for a bending (dipole) element
      if(icat .eq. 2) then
      angr=angd*pi180
      rho=brho/b
      aleng=rho*angr
c
c write all descriptive information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(mfile,607)(pmenu(i+mpp(k)),i=1,imax)
c Output using Mottershead's favorite pg format
      write(mfile,611)  lmnlbl(k), aleng
      write(mfile,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      endif
c
c write all descriptive information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jof,607)(pmenu(i+mpp(k)),i=1,imax)
      write(jof,609) lmnlbl(k), aleng
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jof,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jodf,607)(pmenu(i+mpp(k)),i=1,imax)
      write(jodf,609) lmnlbl(k), aleng
      write(jodf,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jodf,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
c
c subdivide element
      ndptst=ndpts
      if(ndpts .lt. 1) ndptst=1
      daleng=aleng/float(ndptst)
      dangr=angr/float(ndptst)
      if( dabs(aleng) .lt. 1.0d-10) ndptst=0
      do i=0,ndptst
      sst = ss + daleng*float(i)
      ttt = tt + daleng*float(i)*vr
      angrt=dangr*float(i)
      sfact=rho*dsin(angrt)
      cfact=rho*(1.d0 - dcos(angrt))
c Note: signs have been adjusted to take into account that
c the rotation axis is - eyl.
      xt = x -cfact*exl(1) + sfact*ezl(1)
      yt = y -cfact*exl(2) + sfact*ezl(2)
      zt = z -cfact*exl(3) + sfact*ezl(3)
      angrt=-angrt
      call rotv(eyl,angrt,exl,exlt)
      call rotv(eyl,angrt,ezl,ezlt)
c
c write simple element descriptive information 
c and pathlength information to external file kfile
      if (kfile .gt. 0) then
      write(kfile,610) lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad
      endif
c
c
c write simple element descriptive information and pathlength
c and coordinate information to external file lfile
      if (lfile .gt. 0) then
      write(lfile,710)
     & lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad,xt,yt,zt,ttt
      endif
c
c write coordinate information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,236) sst,xt,yt,zt,ttt,0.0
      endif
c
c write coordinate information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,237) sst,xt,yt,zt,ttt
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,237) sst,xt,yt,zt,ttt
      endif
c
c write orientation information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,*) exlt(1), exlt(2), exlt(3)
      write(mfile,*) eyl(1), eyl(2), eyl(3)
      write(mfile,*) ezlt(1), ezlt(2), ezlt(3)
      endif
c
c write orientation information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) exlt(1), exlt(2), exlt(3)
      write(jof,*) eyl(1), eyl(2), eyl(3)
      write(jof,*) ezlt(1), ezlt(2), ezlt(3)
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) exlt(1), exlt(2), exlt(3)
      write(jodf,*) eyl(1), eyl(2), eyl(3)
      write(jodf,*) ezlt(1), ezlt(2), ezlt(3)
      endif
      end do
      endif
c
c procedure for an arc (arc)
      if(icat .eq. 4) then
c
c write all descriptive information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(mfile,607)(pmenu(i+mpp(k)),i=1,imax)
c Output using Mottershead's favorite pg format
      write(mfile,611)  lmnlbl(k), aleng
      write(mfile,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      endif
c
c write all descriptive information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jof,607)(pmenu(i+mpp(k)),i=1,imax)
      write(jof,609)  lmnlbl(k), aleng
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jof,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      write(jodf,607)(pmenu(i+mpp(k)),i=1,imax)
      write(jodf,609)  lmnlbl(k), aleng
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jodf,*) ' s,x,y,z,t;ex,ey,ez along element:'
      endif
c
c subdivide element
      ndptst=ndpts
      if(ndpts .lt. 1) ndptst=1
      daleng=aleng/float(ndptst)
      dd1=d1/float(ndptst)
      dd2=d2/float(ndptst)
      hypot=dsqrt(d1**2 + d2**2)
      angt=0.d0
      if (hypot .gt. 0.d0) then
      sin=d2/hypot
      cos=d1/hypot
      angt=atan2(sin,cos)
      endif
      call rotv(eyl,angt,exl,exlt)
      call rotv(eyl,angt,ezl,ezlt)
      if( dabs(aleng) .lt. 1.0d-10) ndptst=0
      do i=0,ndptst
      sst = ss + daleng*float(i)
      ttt = tt + daleng*float(i)*vr
      xt = x + (dd1*ezl(1)+dd2*exl(1))*float(i)
      yt = y + (dd1*ezl(2)+dd2*exl(2))*float(i)
      zt = z + (dd1*ezl(3)+dd2*exl(3))*float(i)
c
c write simple element descriptive information 
c and pathlength information to external file kfile
      if (kfile .gt. 0) then
      write(kfile,610) lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad
      endif
c
c write simple element descriptive information and pathlength
c and coordinate information to external file lfile
      if (lfile .gt. 0) then
      write(lfile,710)
     & lmnlbl(k),nt1(k),nt2(k),sst,dims(1),dims(2),grad,xt,yt,zt,ttt
      endif
c
c write coordinate information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,236) sst,xt,yt,zt,ttt,0.0
      endif
c
c write coordinate information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,237) sst,xt,yt,zt,ttt
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,237) sst,xt,yt,zt,ttt
      endif
c
c write orientation information to external file mfile
      if (mfile .gt. 0) then
      write(mfile,*) exlt(1), exlt(2), exlt(3)
      write(mfile,*) eyl(1), eyl(2), eyl(3)
      write(mfile,*) ezlt(1), ezlt(2), ezlt(3)
      endif
c
c write orientation information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) exlt(1), exlt(2), exlt(3)
      write(jof,*) eyl(1), eyl(2), eyl(3)
      write(jof,*) ezlt(1), ezlt(2), ezlt(3)
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) exlt(1), exlt(2), exlt(3)
      write(jodf,*) eyl(1), eyl(2), eyl(3)
      write(jodf,*) ezlt(1), ezlt(2), ezlt(3)
      endif
c
      end do
      endif
c
      endif
c
c update geometry for items dims thru prot above
c
c procedure for a straight element
      if(icat .eq. 1) then
      ss = ss + aleng
      tt = tt + aleng*vr
      x = x + aleng*ezl(1)
      y = y + aleng*ezl(2)
      z = z + aleng*ezl(3)
      endif
c
c procedure for a bending (dipole) element
      if(icat .eq. 2) then
      angr=angd*pi180
      rho=brho/b
      aleng=rho*angr
      ss = ss + aleng
      tt = tt + aleng*vr
      sfact=rho*dsin(angr)
      cfact=rho*(1.d0 - dcos(angr))
c Note: signs have been adjusted to take into account that
c the rotation axis is - eyl.
      x = x -cfact*exl(1) + sfact*ezl(1)
      y = y -cfact*exl(2) + sfact*ezl(2)
      z = z -cfact*exl(3) + sfact*ezl(3)
      angr=-angr
      call rotv(eyl,angr,exl,exl)
      call rotv(eyl,angr,ezl,ezl)
      endif
c
c procedure for an arot
      if(icat .eq. 3) then
      angr=angd*pi180
c Note: sign of angr has been adjusted in accord with Figure 6.14c
c of the MaryLie manual.
      angr=-angr
      call rotv(ezl,angr,exl,exl)
      call rotv(ezl,angr,eyl,eyl)
      endif
c
c procedure for an arc
      if(icat .eq. 4) then
      ss = ss + aleng
      tt = tt + aleng*vr
      x = x + d1*ezl(1) + d2*exl(1)
      y = y + d1*ezl(2) + d2*exl(2)
      z = z + d1*ezl(3) + d2*exl(3)
      angr=angd*pi180
      call rotv(eyl,angr,exl,exl)
      call rotv(eyl,angr,ezl,ezl)
      endif
c
c
c response to a data point
c
      if( iopt .eq. 1 .or. iopt .eq. 3
     &  .or. iopt .eq. 11 .or. iopt .eq. 13) then
c
      if (nt2(k) .eq. 23) then
c
c write all information to terminal and/or drop file
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
      write(jof,*) ' dimensions:'
      write(jof,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jof,*) ' s,x,y,z,t;ex,ey,ez:'
      write(jof,237) ss,x,y,z,tt
      write(jof,*) exl(1), exl(2), exl(3)
      write(jof,*) eyl(1), eyl(2), eyl(3)
      write(jof,*) ezl(1), ezl(2), ezl(3)
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
      write(jodf,*) ' dimensions:'
      write(jodf,608) dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)
      write(jodf,*) ' s,x,y,z,t;ex,ey,ez:'
      write(jodf,237) ss,x,y,z,tt
      write(jodf,*) exl(1), exl(2), exl(3)
      write(jodf,*) eyl(1), eyl(2), eyl(3)
      write(jodf,*) ezl(1), ezl(2), ezl(3)
      endif
c
c write to external files
c write coordinate information to file ifile
      if(ifile .gt. 0) then
      write(ifile,520) ss,x,y,z,tt,0.
  520 format(6(1x,1pe12.5))
      endif
c write coordinate and orientation information to file jfile
      if(jfile .gt. 0) then
      write(jfile,520) ss,x,y,z,tt,0.
      write(jfile,*) exl(1), exl(2), exl(3)
      write(jfile,*) eyl(1), eyl(2), eyl(3)
      write(jfile,*) ezl(1), ezl(2), ezl(3)
      endif
c
c put result in ucalc
      if( iopt .eq. 11 .or. iopt .eq. 13) then
      ucalc(1)=ss
      ucalc(2)=x
      ucalc(3)=y
      ucalc(4)=z
      ucalc(5)=tt
      ucalc(6)= exl(1)
      ucalc(7)= exl(2)
      ucalc(8)= exl(3)
      ucalc(9)= eyl(1)
      ucalc(10)=eyl(2)
      ucalc(11)=eyl(3)
      ucalc(12)=ezl(1)
      ucalc(13)=ezl(2)
      ucalc(14)=ezl(3)
      endif
c
      endif
c
      endif
c
      endif
c
      endif
c
c procedure for a lump
c
      if(itype .eq. 3) then
c      write(ifile,515) string(1),mim(jk1)
  515 format(1x,1x,a8,1x,'lump',9x,'0',4x,'-1',2x,i4)
      endif
c
  137 continue
c
      return
      end
c
************************************************************************
c
      subroutine hmltn1(h)
c  this subroutine is used to specify a constant hamiltonian
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'extalk.inc'
      include 'hmflag.inc'
      dimension h(monoms)
c
c  begin computation
      iflag=0
      do 10 i=1,monoms
   10 h(i)=-fa(i)
c
      return
      end
c
********************************************************************
c
      subroutine lnf(p,fa,fm)
c this is a subroutine that takes the log of a map
c in normal form
c written by Alex Dragt 1/30/99
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
c
c set up parameters
c
      job=nint(p(1))
      iopt=nint(p(2))
      mult=nint(p(3))
c
c compute exponent in the static case
c
      if(job .eq. 1) then
      call ident(buf1a,buf1m)
c horizontal plane f2
      cos = fm(1,1)
      sin = fm(1,2)
      ang = datan2(sin,cos)
      buf1a(7) = -ang/2.d0
      buf1a(13) = -ang/2.d0
c vertical plane f2
      cos = fm(3,3)
      sin = fm(3,4)
      ang = datan2(sin,cos)
      buf1a(18) = -ang/2.d0
      buf1a(22) = -ang/2.d0
c temporal plane f2
      endif
c
c compute exponent in the dynamic case
c
      if(job .eq. 2) then
      call ident(buf1a,buf1m)
c horizontal plane f2
      cos = fm(1,1)
      sin = fm(1,2)
      ang = datan2(sin,cos)
      buf1a(7) = -ang/2.d0
      buf1a(13) = -ang/2.d0
c vertical plane f2
      cos = fm(3,3)
      sin = fm(3,4)
      ang = datan2(sin,cos)
      buf1a(18) = -ang/2.d0
      buf1a(22) = -ang/2.d0
c temporal plane f2
      endif
c
c higher-order f's
      if(iopt .eq. 2) then
      do i=28,monoms
      buf1a(i)=fa(i)
      enddo
      endif
c
c scale by (-1/pi)
      if(mult .eq. 2) then
      pi=4.d0*datan(1.d0)
      scale=-1.d0/pi
      call csmul(scale,buf1a,buf1a) 
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine moma(p)
c this subroutine is for moment analysis
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
      character*3 kynd
c Calling arrays
      dimension p(6),fm(6,6),gm(6,6),hm(6,6)
cryne 8/9/2002      dimension ga(monoms5)
cryne 8/9/2002      dimension ha(monoms5)
cryne 8/9/2002      dimension t1a(monoms5)
cryne 8/9/2002      dimension t2a(monoms5)
      dimension ga(monoms)
      dimension ha(monoms)
      dimension t1a(monoms)
      dimension t2a(monoms)
c
c set up control indices
      job=nint(p(1))
      isend=nint(p(2))
      ifile=nint(p(3))
      nopt=nint(p(4))
      nskip=nint(p(5))
      nmpo=nint(p(6))
c
c procedure for reading in function or moments when job = 11,21, or 31:
      if(job.eq.11 .or. job.eq.21 .or. job.eq.31) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      call ident(ga,gm)
      kynd='gtm'
      call strget(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
c
c procedure for reading in moments when job = 4,5, or 6:
c      if (job.eq.4 .or. job.eq.5 .or. job.eq.6) then
c test for file read or internal map fetch
c      if(ifile.lt.0) then
c      nmap=-ifile
c      kynd='gtm'
c	call strget5(kynd,nmap,ga,gm)
c      else
c      mpit=mpi
c      mpi=ifile
c      call mapin5(nopt,nskip,ga,gm)
c      mpi=mpit
c      endif
c      endif
c      continue
c
c compute eigen emittances 
      if(job .eq. 11) call eigemt(2,ga)
      if(job. eq. 21) call eigemt(4,ga)
      if(job .eq. 31) call eigemt(6,ga)
c compute mean square eigen-emittances and put results in fitbuf
      wex=buf1a(7)**2
      wey=buf1a(18)**2
      wet=buf1a(25)**2
c
c write out results if desired
c code needs to be modified to write out eigen emittances
      if (isend.eq.1.or.isend.eq.3) then
      write (jof,*) 'xee2=',wex
      write (jof,*) 'yee2=',wey
      write (jof,*) 'tee2=',wet
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write (jodf,*) 'xee2=',wex
      write (jodf,*) 'yee2=',wey
      write (jodf,*) 'tee2=',wet
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine padd(p,ha,hm)
c this subroutine adds two polynomials
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension p(6),ha(monoms),hm(6,6)
c
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
c
c set up control indices
      nmapf=nint(p(1))
      nmapg=nint(p(2))
      nmaph=nint(p(3))
c
c get maps and clear arrays
      kynd='gtm'
      if (nmapf.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) call strget(kynd,nmapf,fa,fm)
      if (nmapg.eq.0) call mapmap(ha,hm,ga,gm)
      if (nmapg.ge.1 .and. nmapg.le.5) call strget(kynd,nmapg,ga,gm)
      call ident(ta,tm)
c
c perform calculation
      call cpadd(fa,ga,ta)
c
c decide where to put results
c
      if (nmaph.ge.1 .and. nmaph.le.5) then 
      kynd='stm'
      call strget(kynd,nmaph,ta,tm)      
      endif
c
      if (nmaph.eq.0) call mapmap(ta,tm,ha,hm)
c 
      if (nmaph.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmaph.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmaph.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmaph.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmaph.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
********************************************************************
c
      subroutine pbpol(p,fa,fm)
c this subroutine poisson brackets two polynomials
c written 5/22/02 AJD
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      dimension p(6),fa(monoms),fm(6,6)
      dimension ta(monoms),tm(6,6)
      dimension t1a(monoms)
      dimension t2a(monoms)
      character*3 kynd
c
      write (6,*) 'in subroutine pbpol'
c
c set up and test parameters
      map1in=nint(p(1))
      map2in=nint(p(2))
      mapout=nint(p(3))
      if((map1in .lt. 0) .or. (map1in .gt. 9)) go to 100
      if((map2in .lt. 0) .or. (map2in .gt. 9)) go to 100
      if((mapout .lt. -6) .or. (mapout .gt. 9)) go to 100
c
c get maps
      if(map1in .eq. 0) then
      call mapmap(fa,fm,t1a,tm)
      else
      kynd='gtm'
      call strget(kynd,map1in,t1a,tm)
      endif
      if(map2in .eq. 0) then
      call mapmap(fa,fm,t2a,tm)
      else
      kynd='gtm'
      call strget(kynd,map2in,t2a,tm)
      endif
c
c compute Poisson bracket
      call mclear(tm)
      call cppb(t1a,t2a,ta)
c
c deposit result and return
      if(mapout .lt. 0) then
      kbuf=-mapout
      if(kbuf .eq. 1) call mapmap(ta,tm,buf1a,buf1m)
      if(kbuf .eq. 2) call mapmap(ta,tm,buf2a,buf2m)
      if(kbuf .eq. 3) call mapmap(ta,tm,buf3a,buf3m)
      if(kbuf .eq. 4) call mapmap(ta,tm,buf4a,buf4m)
      if(kbuf .eq. 5) call mapmap(ta,tm,buf5a,buf5m)
      if(kbuf .eq. 6) call mapmap(ta,tm,buf6a,buf6m)
      endif
      if(mapout .eq. 0) call mapmap(ta,tm,fa,fm)
      if(mapout .gt. 0) then
      kynd='stm'
      call strget(kynd,mapout,ta,tm)
      endif
      return
c
 100  continue
c error return
      write(6,*) 'control parameters out of bounds in pbpol'
      write(6,*)'map1in=',map1in
      write(6,*)'map2in=',map2in
      write(6,*)'mapout=',mapout
      return
      end
c
**********************************************************************
c
      subroutine pdnf(p,ha,hm)
c this subroutine computes powers of a dynamic normal form
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      character*3 kynd
c
      dimension p(6),ha(monoms),hm(6,6)
c
      dimension fa(monoms),ta(monoms)
      dimension fm(6,6),tm(6,6)
c
c set up control indices
      jinopt=nint(p(1))
      if (jinopt.eq.1) pow=p(2)
      if (jinopt.eq.2) npowf=nint(p(2))
      nmapi=nint(p(3))
      joutop=nint(p(4))
      nmapo=nint(p(5))
c
c get map and clear arrays
      kynd='gtm'
      if (nmapi.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapi.ge.1 .and. nmapi.le.5) call strget(kynd,nmapi,fa,fm)
      call ident(ta,tm)
c
c perform calculation
c
c procedure when jinopt=1
      if (jinopt.eq.1) then
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
c
c procedure when jinopt=2
      if (jinopt.eq.2) then
c
c procedure when npowf .lt. 0
      if (npowf.lt.0) then
      do 10 k=1,6
      pow=0.d0
c      if (npowf.eq.-1) pow=pst1(k)
c      if (npowf.eq.-2) pow=pst2(k)
c      if (npowf.eq.-3) pow=pst3(k)
c      if (npowf.eq.-4) pow=pst4(k)
c      if (npowf.eq.-5) pow=pst5(k)
      pow=pst(-npowf,k)
      if (pow.ne.0.d0) then
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
   10 continue
      endif
c procedure when npowf .gt.0
      if (npowf.gt.0) then
   20 continue
      read(npowf,*,end=30) pow
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      goto 20
   30 continue
      endif
c
      endif
c
      return
      end
c
********************************************************************
c
      subroutine pmul(p,ha,hm)
c this subroutine multiplies two polynomials
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension p(6),ha(monoms),hm(6,6)
c
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
c
c set up control indices
      nmapf=nint(p(1))
      nmapg=nint(p(2))
      nmaph=nint(p(3))
c
c get maps and clear arrays
      kynd='gtm'
      if (nmapf.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) call strget(kynd,nmapf,fa,fm)
      if (nmapg.eq.0) call mapmap(ha,hm,ga,gm)
      if (nmapg.ge.1 .and. nmapg.le.5) call strget(kynd,nmapg,ga,gm)
      call ident(ta,tm)
c
c perform calculation
      call cpmul(fa,ga,ta)
c
c decide where to put results
c
      if (nmaph.ge.1 .and. nmaph.le.5) then 
      kynd='stm'
      call strget(kynd,nmaph,ta,tm)      
      endif
c
      if (nmaph.eq.0) call mapmap(ta,tm,ha,hm)
c 
      if (nmaph.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmaph.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmaph.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmaph.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmaph.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine pnlp(p,ga,gm)
c this subroutine raises the nonlinear part of a map to a power
c Written by Alex Dragt, Fall 1988
c Fifth order by F. Neri (1989).
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      character*3 kynd
c
c Calling arrays
      dimension p(6),ga(monoms),gm(6,6)
c
c Local arrays
      dimension fa(monoms)
      dimension fm(6,6)
      dimension t5(461),t6(923),ff(monoms)
c
c set up scalar and control indices
      iopt=nint(p(1))
      power=p(2)
      nmapf=nint(p(3))
      nmapg=nint(p(4))
c
c get map
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
c
c perform calculation
      if(iopt.eq.0) then
      call mclear(fm)
      do 10 i=1,6
  10  fm(i,i)=1.d0
      endif
      call csmul(power,fa,ff)
c
c f5 commutator:
      call pbkt1(fa,3,fa,4,t5)
      call pmadd(t5,5,(power*(1.d0-power)/2.d0),ff)
c f6 commutators:
      call pbkt1(fa,3,t5,5,t6)
      call pmadd(t6,6,(power-3.d0*power**2+2.d0*power**3)/12.d0,ff)
      call pbkt1(fa,3,fa,5,t6)
      call pmadd(t6,6,(power*(1.d0-power)/2.d0),ff)
c
c copy result back to fa:
      call mapmap(ff,fm,fa,fm)
c
c decide where to put results
c
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='stm'
      call strget(kynd,nmapg,fa,fm)
      endif
c
      if (nmapg.eq.0) call mapmap(fa,fm,ga,gm)
c
      if (nmapg.eq.-1) call mapmap(fa,fm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(fa,fm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(fa,fm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(fa,fm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(fa,fm,buf5a,buf5m)
c
      return
      end
c
********************************************************************************
c
      subroutine pold(p,fa,fm)
c this is a subroutine for polar decomposition
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ta(monoms)
      dimension tm(6,6),rm(6,6),pdsm(6,6),revec(6,6)
      dimension reval(6)
c
c set up control indices
      mapin=nint(p(1))
      isend=nint(p(2))
      idata=nint(p(3))
      ipmaps=nint(p(4))
      iwmaps=nint(p(5))
c
c write heading(s)
c
c begin calculation
      call polr(fa,fm,rm,pdsm,reval,revec)
c put out data
      call ident(ta,tm)
c
c put out maps
c
c put maps in buffers
      call ident(ta,tm)
      call mapmap(ta,pdsm,buf1a,buf1m)
      call mapmap(ta,rm,buf2a,buf2m)
      call mapmap(fa,tm,buf3a,buf3m)
      do 10 i=1,6
      tm(i,i)=reval(i)
   10 continue
      call mapmap(ta,tm,buf4a,buf4m)
      call mapmap(ta,revec,buf5a,buf5m)
c
c write out maps
c
      return
      end
c
**********************************************************************
c
      subroutine ppa(p,fa,fm)
c
c     This routine computes focal lengths and principal planes for the
c     current map.
c         C. T. Mottershead  LANL AT-3 / 2 Oct 89
c---------------------------------------------------------------------
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'fitdat.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6)
c
      isend = int(p(1))
      if(iquiet.eq.1) isend = 0
      eps = 1.e-9
c
c       x-plane
c
      finv = fm(2,1)
      if(abs(finv).lt.eps) finv = eps
      f = -1.0/finv
      z2 = f*(1.0 - fm(1,1))
      z1 = f*(1.0 - fm(2,2))
      fx = f
      xb = z1
      xa = z2
      xu = f - z1
      xd = f - z2
c
c       y-plane
c
      finv = fm(4,3)
      if(abs(finv).lt.eps) finv = eps
      f = -1.0/finv
      z2 = f*(1.0 - fm(3,3))
      z1 = f*(1.0 - fm(4,4))
      fy = f
      yb = z1
      ya = z2
      yu = f - z1
      yd = f - z2
c
c     print the matrix and focal lengths
c
      if(isend.lt.2) go to 400
      lun = jodf
 300  continue
      write(lun,17) fx,fy
  17  format(' * PPA  Focal lengths : fx =',1pg15.7,'  fy =',1pg15.7)
      write(lun,33)
  33  format(' * Principal Planes (before and after):')
      write(lun,37) xb,xa,yb,ya
  37  format(1x,' xb=',1pg15.7,' xa=',1pg15.7,' yb=',1pg15.7,' ya=',
     & 1pg15.7)
      write(lun,27) xu,yu
  27  format(' * Focal points (upstream): xu =',1pg15.7,
     & '  yu =',1pg15.7)
      write(lun,29) xd,yd
  29  format('              (downstream): xd =',1pg15.7,
     & '  yd =',1pg15.7)
 400  if((isend.eq.1).or.(isend.eq.3)) then
         lun = jof
         isend = 0
         go to 300
      endif
c
      return
      end
c
**********************************************************************
c
      subroutine psnf(p,ha,hm)
c this subroutine computes powers of a static normal form
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      character*3 kynd
c
      dimension p(6),ha(monoms),hm(6,6)
c
      dimension fa(monoms),ta(monoms)
      dimension fm(6,6),tm(6,6)
c
c set up control indices
      jinopt=nint(p(1))
      if (jinopt.eq.1) pow=p(2)
      if (jinopt.eq.2) npowf=nint(p(2))
      nmapi=nint(p(3))
      joutop=nint(p(4))
      nmapo=nint(p(5))
c
c get map and clear arrays
      kynd='gtm'
      if (nmapi.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapi.ge.1 .and. nmapi.le.5) call strget(kynd,nmapi,fa,fm)
      call ident(ta,tm)
c
c perform calculation
c
c procedure when jinopt=1
      if (jinopt.eq.1) then
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
c
c procedure when jinopt=2
      if (jinopt.eq.2) then
c
c procedure when npowf .lt. 0
      if (npowf.lt.0) then
      do 10 k=1,6
      pow=0.d0
c      if (npowf.eq.-1) pow=pst1(k)
c      if (npowf.eq.-2) pow=pst2(k)
c      if (npowf.eq.-3) pow=pst3(k)
c      if (npowf.eq.-4) pow=pst4(k)
c      if (npowf.eq.-5) pow=pst5(k)
      pow = pst(-npowf,k)
      if (pow.ne.0.d0) then
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
   10 continue
      endif
c procedure when npowf .gt.0
      if (npowf.gt.0) then
   20 continue
      read(npowf,*,end=30) pow
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      goto 20
   30 continue
      endif
c
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine psp(p,ha,hm)
c this subroutine computes the scalar product of two polynomials
c Written by Alex Dragt, 10/23/89
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      include 'usrdat.inc'
c
      character*3 kynd
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ga(monoms)
      dimension tm(6,6)
c
c set up control indices
      job=nint(p(1))
      nmapf=nint(p(2))
      nmapg=nint(p(3))
      isend=nint(p(4))
c
c  clear arrays and get maps
      if (nmapf.eq.0) call mapmap(ha,hm,fa,tm)
      if (nmapf.ge.1 .and. nmapf.le.9) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,tm)
      endif
      if (nmapg.eq.0) call mapmap(ha,hm,ga,tm)
      if (nmapg.ge.1 .and. nmapg.le.9) then
      kynd='gtm'
      call strget(kynd,nmapg,ga,tm)
      endif
c
c perform calculation
      call cpsp(job,fa,ga,ans1,ans2,ans3,ans4)
ccccc call old_cpsp(fa,ga,ans1,ans2,ans3,ans4)
c
c decide where to send and put results
c
c     write(6,*) ' ans1=',ans1,' ans2=',ans2
c     write(6,*) ' ans3=',ans3,' ans4=',ans4
      write(6,"('ans1234=',4(1pe15.8,1x))")ans1,ans2,ans3,ans4
      ucalc(141)=ans1
      ucalc(142)=ans2
      ucalc(143)=ans3
      ucalc(144)=ans4
      ucalc(145)=ans1+ans2+ans3+ans4
c
      return
      end
c
********************************************************************
c
      subroutine pval(p,ga,gm)
c this is a routine for evaluating a polynomial
c Written by Alex Dragt, Spring 1987
      use rays
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'parset.inc'
      character*3 kynd
c
      dimension p(6),ga(monoms),gm(6,6)
c
      dimension fa(monoms),fm(6,6)
      dimension zit(6)
c
c set up control indices
      mapin=nint(p(1))
      idata=nint(p(2))
      iwnum=nint(p(3))
c
c get polynomial
c
      if (mapin.eq.0) call mapmap(ga,gm,fa,fm)
      if (mapin.ge.1 .and. mapin.le.5) then
      kynd='gtm'
      call strget(kynd,mapin,fa,fm)
      endif
c
c compute value(s) of polynomial and write them out
c
      if (iwnum.gt.0) then
      write(jof,*)
      write(jof,*) 'value(s) of polymomial written on file ',iwnum
c
c procedure when idata > 0
      if (idata.gt.0) then
      ipset=idata
c get phase space data from the parameter set ipset
c
c      do 5 i=1,6
c    5 zit(i)=0.d0
c      goto (10,20,30,40,50),ipset
c      goto 60
c   10 do 11 i=1,6
c   11 zit(i)=pst1(i)
c      goto 60
c   20 do 21 i=1,6
c   21 zit(i)=pst2(i)
c      goto 60
c   30 do 31 i=1,6
c   31 zit(i)=pst3(i)
c      goto 60
c   40 do 41 i=1,6
c   41 zit(i)=pst4(i)
c      goto 60
c   50 do 51 i=1,6
c   51 zit(i)=pst5(i)
c   60 continue
      if(ipset.lt.1 .or. ipset.gt.maxpst) then
       do 50 i=1,6
  50   zit(i) = 0.0d0
      else
       do 60 i=1,6
  60   zit(i) = pst(ipset,i)
      endif
c carry out computation
      call evalf(zit,fa,val2,val3,val4)
      k=1
      write(iwnum,100) k,val2,val3,val4,0.,0.
  100 format(1x,i12,5(1x,1pe12.5))
      endif
c
c procedure when idata = 0
      if (idata.eq.0) then
      do 70 k=1,nraysp
c check on status of k'th ray; skip this ray if its status is 'lost'
      if (istat(k).ne.0) goto 70
      do 80 j=1,6
   80 zit(j)=zblock(j,k)
      call evalf(zit,fa,val2,val3,val4)
      write(iwnum,100) k,val2,val3,val4,0.,0.
   70 continue
      endif
c
      endif
c
      return
      end
c
**********************************************************************
c
      subroutine radm(p,fa,fm)
c this is a subroutine for resonance analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension tm(6,6),t1m(6,6),t2m(6,6)
      dimension look(3),pmask(6)
c
c set up control indices
      iopt=nint(p(1))
      i2=nint(p(2))
      i3=nint(p(3))
      i4=nint(p(4))
      iwmaps=nint(p(5))
c
c compute isend
      do 10 j=1,3
      look(j)=0
      if (i2.eq.j .or. i3.eq.j .or. i4.eq.j) look(j)=1
   10 continue
      isend=0
      if (look(1).eq.1) isend=1
      if (look(2).eq.1) isend=2
      if (look(1).eq.1 .and. look(2).eq.1) isend=3
      if (look(3).eq.1) isend=3
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*)
      write(jof,*) 'resonance analysis of dynamic map'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'resonance analysis of dynamic map'
      endif
c
c beginning of calculation
c
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,buf2a,buf2m,ta,tm)
c
c procedure for removing third order terms
      if(iopt.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'third order terms removed'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'third order terms removed'
      endif
      call dpur3(buf2a,buf2m,buf3a,buf3m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c put maps in proper places
      call mapmap(buf3a,buf3m,buf2a,buf2m)
      call mapmap(t2a,t2m,ta,tm)
      endif
c
c resonance decompose purified map:
      call matmat(buf2m,buf3m)
      call ctodr(buf2a,buf3a)
c
c procedure for writing resonance driving terms at terminal (file jof)
      if (isend.eq.1 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.1 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.1 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.1 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jof,*)
      write(jof,*) 'requested resonance driving terms written as a map'
      call pdrmap(0,1,buf4a,buf4m)
      endif
c
c procedure for writing resonance driving terms on external file
c (file jodf)
      if (isend.eq.2 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.2 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.2 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.2 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jodf,*)
      write(jodf,*) 'requested resonance driving terms written as a map'
      call pdrmap(0,2,buf4a,buf4m)
      endif
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ta,tm)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c put the transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c buffers 2 and 3 already contain the purified map in the cartesian
c and dynamic resonance bases, respectively
c clear the remaining buffers
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine rasm(p,fa,fm)
c this is a subroutine for resonance analysis of static maps
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension tm(6,6),t1m(6,6),t2m(6,6)
      dimension look(3),pmask(6)
 
c
c set up control indices
      iopt=nint(p(1))
      i2=nint(p(2))
      i3=nint(p(3))
      i4=nint(p(4))
      iwmaps=nint(p(5))
c
c compute isend
      do 10 j=1,3
      look(j)=0
      if (i2.eq.j .or. i3.eq.j .or. i4.eq.j) look(j)=1
   10 continue
      isend=0
      if (look(1).eq.1) isend=1
      if (look(2).eq.1) isend=2
      if (look(1).eq.1 .and. look(2).eq.1) isend=3
      if (look(3).eq.1) isend=3
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*)
      write(jof,*) 'resonance analysis of static map'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'resonance analysis of static map'
      endif
c
c beginning of calculation
c
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,buf2a,buf2m,ta,tm,t2m)
c
c procedure for removing third order terms
      if(iopt.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'third order terms removed'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'third order terms removed'
      endif
c remove offensive chromatic terms from f3 part of map:
      call scpur3(buf2a,buf2m,buf3a,buf3m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(buf3a,buf3m,buf2a,buf2m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,t2a,t2m,ta,tm)
      endif
c
c resonance decompose purified map:
      call matmat(buf2m,buf3m)
      call ctosr(buf2a,buf3a)
c
c procedure for writing resonance driving terms at terminal (file jof)
      if (isend.eq.1 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.1 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.1 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.1 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jof,*)
      write(jof,*) 'requested resonance driving terms written as a map'
      call psrmap(0,1,buf4a,buf4m)
      endif
c
c procedure for writing resonance driving terms on external file (file jodf)
      if (isend.eq.2 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.2 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.2 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.2 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jodf,*)
      write(jodf,*) 'requested resonance driving terms written as a map'
      call psrmap(0,2,buf4a,buf4m)
      endif
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ta,tm)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c put the transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c buffers 2 and 3 already contain the purified map in the cartesian
c and static resonance bases, respectively
c clear the remaining buffers
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine sia(p,fa,fm)
c this is a routine for computing invariants in the static case
c Written by Alex Dragt, Spring 1987
      use rays
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c
c set up control indices
      iopt=nint(p(1))
      idata=nint(p(2))
      ipmaps=nint(p(3))
      isend=nint(p(4))
      iwmaps=nint(p(5))
      iwnum=nint(p(6))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'static invariant analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'static invariant analysis'
      endif
c
c begin calculation
c
c find the transforming (conjugating) map script A
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,ga,gm,t1a,t1m,t2m)
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t2a,t2m)
c accumulate transforming map:
      call concat(t2a,t2m,t1a,t1m,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c remove offensive terms from f4 part of map:
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c put script A in buffer 1 and script N in buffer 2
      call mapmap(ta,tm,buf1a,buf1m)
      call mapmap(g1a,g1m,buf2a,buf2m)
c
c procedure for computing invariant
c invert script A
      call inv(ta,tm)
      call clear(t1a,t1m)
c computation of x invariant
      if(iopt.eq.1) then
      t1a(7)=1.d0
      t1a(13)=1.d0
      endif
c computation of y invariant
      if(iopt.eq.2) then
      t1a(18)=1.d0
      t1a(22)=1.d0
      endif
c computation of mixed invariant
      if(iopt.ge.3) then
      read(iopt,*) anx,any
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      t1a(7)=anx
      t1a(13)=anx
      t1a(18)=any
      t1a(22)=any
      endif
c put invariant in buffer 3
      call ident(buf3a,buf3m)
      call fxform(ta,tm,t1a,buf3a)
c
c procedure for putting out data
      if (idata.eq.1 .or. idata.eq.3) then
      do 10 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 10
      write(ifile,*)
      write(ifile,*) 'invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
   10 continue
      endif
      if (idata.eq.2 .or. idata.eq.3) then
      mpot=mpo
      mpo=18
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c procedure for printing out maps
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form map script N'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   20 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c procedure for computing values of invariant and writing them out
      if (iwnum.gt.0) then
      write(jof,*)
      write(jof,*) 'values of invariant written on file ',iwnum
      do 30 k=1,nraysp
      do 40 j=1,6
   40 zi(j)=zblock(j,k)
      call evalf(zi,buf3a,val2,val3,val4)
      write(iwnum,60) k,val2,val3,val4,0.,0.
   60 format(1x,i12,5(1x,1pe12.5))
   30 continue
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the invariant polynomial.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m) 
c
      return
      end
c
***********************************************************************
c
      subroutine smul(p,ga,gm)
c this subroutine multiplies a polynomial by a scalar
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension p(6),ga(monoms),gm(6,6)
c
      dimension fa(monoms),ta(monoms)
      dimension fm(6,6),tm(6,6)
c
c set up scalar and control indices
      scalar=p(1)
      nmapf=nint(p(2))
      nmapg=nint(p(3))
c
c get map and clear arrays
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
      call clear(ta,tm)
c
c perform calculation
      call csmul(scalar,fa,ta)
c
c decide where to put results
c
      if (nmapg.ge.1 .and. nmapg.le.5) then 
      kynd='stm'
      call strget(kynd,nmapg,ta,tm)      
      endif
c
      if (nmapg.eq.0) call mapmap(ta,tm,ga,gm)
c 
      if (nmapg.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine snor_old(p,fa,fm)
c this is a subroutine for normal form analysis of static maps
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c
c set up control indices
      idata=nint(p(1))
      ipmaps=nint(p(2))
      isend=nint(p(3))
      iwmaps=nint(p(4))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'static normal form analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'static normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,ga,gm,t1a,t1m,t2m)
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t2a,t2m)
c accumulate transforming map:
      call concat(t2a,t2m,t1a,t1m,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c remove offensive terms from f4 part of map:
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(g1a,g1m,buf2a,buf2m)  
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=g1m(1,1)
      swx=g1m(1,2)
      wx=atan2(swx,cwx)
      cwy=g1m(3,3)
      swy=g1m(3,4)
      wy=atan2(swy,cwy)
c compute momentum compaction
      wt=g1m(5,6)
c set up normal form for exponent
      do 10 i=1,27
   10 g1a(i)=0.
      g1a(7)=-wx/2.d0
      g1a(13)=-wx/2.d0
      g1a(18)=-wy/2.d0
      g1a(22)=-wy/2.d0
      g1a(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or.idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,g1a,ga)
c store results in buffer 3
      call mapmap(ga,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,g1a,g1m)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m) 
c
      return
      end
c
***********************************************************************
c
      subroutine snor(p,fa,fm)
c this is a subroutine for normal form analysis of static maps
c Written by Alex Dragt, Spring 1987
c Modified 20 June 1988
c
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6)
      dimension tm(6,6),t1m(6,6)
c
c set up control indices
      keep=  nint(p(1))
      idata= nint(p(2))
      ipmaps=nint(p(3))
      isend= nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if ((isend.eq.1 .or. isend.eq.3) .and. idproc.eq.0) then
      write (jof,*)
      write (jof,*) 'static normal form analysis'
      endif
      if ((isend.eq.2 .or. isend.eq.3) .and. idproc.eq.0) then
      write (jodf,*)
      write (jodf,*) 'static normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,ga,gm,ta,tm,t1m)
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(g1a,g1m,buf2a,buf2m)
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=g1m(1,1)
      swx=g1m(1,2)
      wx=atan2(swx,cwx)
      cwy=g1m(3,3)
      swy=g1m(3,4)
      wy=atan2(swy,cwy)
c compute momentum compaction
      wt=g1m(5,6)
c set up normal form for exponent
      do 10 i=1,27
   10 g1a(i)=0.
      g1a(7)=-wx/2.d0
      g1a(13)=-wx/2.d0
      g1a(18)=-wy/2.d0
      g1a(22)=-wy/2.d0
      g1a(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or.idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,g1a,ga)
c store results in buffer 3
      call mapmap(ga,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,g1a,g1m)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)
     &write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
c
***********************************************************************
c
      subroutine submn(p,ha,hm)
c This subroutine computes the norm of a matrix
c Written by Alex Dragt, 10/20/90
c
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension tm(6,6)
c
c set up control indices
      iopt=nint(p(1))
      isend=nint(p(2))
c
c copy matrix
      call matmat(hm,tm)
c
c perform calculation
c
c modify matrix if required
      if (iopt .eq. 1) then
      do 10 i=1,6
 10   tm(i,i) = tm(i,i) -1.d0
      endif
c
c compute norm
      call mnorm(tm,ans)
c
c decide where to send and put results
c
      if(idproc.eq.0)write(6,*) ' matrix norm = ',ans
c
      return
      end
c
***********************************************************************
c
      subroutine tadm(p,fa,fm)
c this is a subroutine for twiss analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
c this program will eventually have to be rewritten to improve the
c output data and its format
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'buffer.inc'
c
      dimension p(6),fa(monoms),fm(6,6)
c
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c
c set up control indices
      idata=nint(p(1))
      ipmaps=nint(p(2))
      isend=nint(p(3))
      iwmaps=nint(p(4))
c
c write headings
      if ((isend.eq.1 .or. isend.eq.3) .and. idproc.eq.0) then
      write(jof,*)
      write(jof,*) 'twiss analysis of dynamic map'
      endif
      if ((isend.eq.2 .or. isend.eq.3) .and. idproc.eq.0) then
      write(jodf,*)
      write(jodf,*) 'twiss analysis of dynamic map'
      endif
c
c beginning of calculation
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm)
c store purifying map script A2 in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c
c compute tunes:
      cwx=gm(1,1)
      swx=gm(1,2)
      cwy=gm(3,3)
      swy=gm(3,4)
      cwt=gm(5,5)
      swt=gm(5,6)
      pi=4.*atan(1.d0)
      wx=atan2(swx,cwx)
      if (wx.lt.0.) wx=wx+2.*pi
      wy=atan2(swy,cwy)
      if (wy.lt.0.) wy=wy+2.*pi
      wt=atan2(swt,cwt)
c     if (wt.lt.0.) wt=wt+2.*pi
      tnx=wx/(2.*pi)
      tny=wy/(2.*pi)
      tnt=wt/(2.*pi)
c
c remove f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,t2a,t2m)
c resonance decompose purified map:
      call ctodr(g1a,ga)
c
c compute dependence of tune on betatron amplitude:
      hh=ga(84)
      vv=ga(85)
      tt=ga(86)
      hv=ga(87)
      ht=ga(88)
      vt=ga(89)
      hhn=-2.d0*hh/pi
      vvn=-2.d0*vv/pi
      ttn=-2.d0*tt/pi
      hvn=-hv/pi
      htn=-ht/pi
      vtn=-vt/pi
c
c remove f4 part of map
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c
c compute envelopes and twiss parameters 
c use the map in buffer 1
      call mapmap(buf1a,buf1m,t1a,t1m)
c compute envelopes
      exh2=t1m(1,1)**2+t1m(1,2)**2
      exv2=t1m(1,3)**2+t1m(1,4)**2
      ext2=t1m(1,5)**2+t1m(1,6)**2
      epxh2=t1m(2,1)**2+t1m(2,2)**2
      epxv2=t1m(2,3)**2+t1m(2,4)**2
      epxt2=t1m(2,5)**2+t1m(2,6)**2
      eyh2=t1m(3,1)**2+t1m(3,2)**2
      eyv2=t1m(3,3)**2+t1m(3,4)**2
      eyt2=t1m(3,5)**2+t1m(3,6)**2
      epyh2=t1m(4,1)**2+t1m(4,2)**2
      epyv2=t1m(4,3)**2+t1m(4,4)**2
      epyt2=t1m(4,5)**2+t1m(4,6)**2
      eth2=t1m(5,1)**2+t1m(5,2)**2
      etv2=t1m(5,3)**2+t1m(5,4)**2
      ett2=t1m(5,5)**2+t1m(5,6)**2
      epth2=t1m(6,1)**2+t1m(6,2)**2
      eptv2=t1m(6,3)**2+t1m(6,4)**2
      eptt2=t1m(6,5)**2+t1m(6,6)**2
      exh=sqrt(exh2)
      exv=sqrt(exv2)
      ext=sqrt(ext2)
      epxh=sqrt(epxh2)
      epxv=sqrt(epxv2)
      epxt=sqrt(epxt2)
      eyh=sqrt(eyh2)
      eyv=sqrt(eyv2)
      eyt=sqrt(eyt2)
      epyh=sqrt(epyh2)
      epyv=sqrt(epyv2)
      epyt=sqrt(epyt2)
      eth=sqrt(eth2)
      etv=sqrt(etv2)
      ett=sqrt(ett2)
      epth=sqrt(epth2)
      eptv=sqrt(eptv2)
      eptt=sqrt(eptt2)
c compute twiss parameters
      call inv(t1a,t1m)
c compute invariants using buffers 2 thru 5
c computation of x invariant
      call clear(buf2a,buf2m)
      buf2a(7)=1.d0
      buf2a(13)=1.d0
      call fxform(t1a,t1m,buf2a,buf3a)
c computation of y invariant
      call clear(buf2a,buf2m)
      buf2a(18)=1.d0
      buf2a(22)=1.d0
      call fxform(t1a,t1m,buf2a,buf4a)
c computation of t invariant
      call clear(buf2a,buf2m)
      buf2a(25)=1.d0
      buf2a(27)=1.d0
      call fxform(t1a,t1m,buf2a,buf5a)
c get twiss parameters from the invariants
c terms for horizontal (x) plane
c 'diagonal' terms
      ax=buf3a(8)/2.d0
      bx=buf3a(13)
      gx=buf3a(7)
c skew terms
c remove diagonal terms, and later print buf3a as a map
      buf3a(8)=0.
      buf3a(13)=0.
      buf3a(7)=0.
c terms for vertical (y) plane
c 'diagonal' terms
      ay=buf4a(19)/2.d0
      by=buf4a(22)
      gy=buf4a(18)
c skew terms
c remove diagonal terms, and later print buf4a as a map
      buf4a(19)=0.
      buf4a(22)=0.
      buf4a(18)=0.
c terms for temporal (t) plane
c 'diagonal' terms
      at=buf5a(26)/2.d0
      bt=buf5a(27)
      gt=buf5a(25)
c skew terms
c remove diagonal terms, and later print buf5a as a map
      buf5a(26)=0.
      buf5a(27)=0.
      buf5a(25)=0.
c
c procedure for writing out tunes and anharmonicities
      if (idata.eq.1 .or. idata.eq.12 .or.
     & idata.eq.13 .or. idata.eq.123) then
      do 10 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 10
      if (iflag.eq.0) goto 10
c write out tunes:
      if(idproc.eq.0)then
      write (ifile,*)
      write (ifile,*) 'horizontal tune =',tnx
      write (ifile,*) 'vertical tune =',tny
      write (ifile,*) 'temporal tune =',tnt
c write out normalized anharmonicities
      write(ifile,*)
      write(ifile,*) 'normalized anharmonicities'
      write(ifile,*) ' hhn=',hhn
      write(ifile,*) ' vvn=',vvn
      write(ifile,*) ' ttn=',ttn
      write(ifile,*) ' hvn=',hvn
      write(ifile,*) ' htn=',htn
      write(ifile,*) ' vtn=',vtn
      endif
   10 continue
      endif
c
c procedure for printing out twiss parameters and envelopes
      if (idata.eq.2 .or. idata.eq.12 .or.
     & idata.eq.23 .or. idata.eq.123) then
      do 20 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 20
      if (iflag.eq.0) goto 20
      if(idproc.eq.0)then
      write(ifile,*)
      write(ifile,*) 'horizontal twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) ax,bx,gx
      write(ifile,*) 'skew terms written as a map'
      endif
      call pcmap(0,i,0,0,buf3a,buf3m)
      if(idproc.eq.0)then
      write(ifile,*)
      write(ifile,*) 'vertical twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) ay,by,gy
      write(ifile,*) 'skew terms written as a map'
      endif
      call pcmap(0,i,0,0,buf4a,buf4m)
      if(idproc.eq.0)then
      write(ifile,*)
      write(ifile,*) 'temporal twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) at,bt,gt
      write(ifile,*) 'skew terms written as a map'
      endif
      call pcmap(0,i,0,0,buf5a,buf5m)
      if(idproc.eq.0)then
      write(ifile,*)
      write(ifile,*) 'horizontal envelopes (exh,exv,ext;epxh,epxv,epxt)'
      write(ifile,*) exh,exv,ext
      write(ifile,*) epxh,epxv,epxt
      write(ifile,*)
      write(ifile,*) 'vertical envelopes (eyh,eyv,eyt;epyh,epyv,epyt)'
      write(ifile,*) eyh,eyv,eyt
      write(ifile,*) epyh,epyv,epyt
      write(ifile,*)
      write(ifile,*) 'temporal envelopes (eth,etv,ett;epth,eptv,eptt)'
      write(ifile,*) eth,etv,ett
      write(ifile,*) epth,eptv,eptt
      endif
   20 continue
      endif
c
c procedure for printing out eigenvectors
      if (idata.eq.3 .or. idata.eq.13 .or.
     & idata.eq.23 .or. idata.eq.123) then
      do 30 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 30
      if (iflag.eq.0) goto 30
c write out the matrix buf1m
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'matrix of eigenvectors'
      call pcmap (i,0,0,0,buf1a,buf1m)
   30 continue
      endif
c
c procedure for printing of maps
c
c put out script A
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      do 40 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 40
      if (iflag.eq.0) goto 40
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'transforming map script A'
      call pcmap(i,i,0,0,ta,tm)
   40 continue
      endif
c
c put out script N
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      do 50 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 50
      if (iflag.eq.0) goto 50
      if(idproc.eq.0)write(ifile,*)
      if(idproc.eq.0)write(ifile,*) 'normal form map script N'
      call pcmap(i,i,0,0,ga,gm)
   50 continue
      endif
c
c procedure for writing of maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,ta,tm)
      call mapout(0,ga,gm)
      mpo=mpot
      endif
c
c put maps in buffers
c buffer 1 already contains script A2
      call mapmap(ta,tm,buf2a,buf2m)
      call mapmap(ga,gm,buf3a,buf3m)
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
      subroutine tasm(p,fa,fm)
c this is a subroutine for twiss analysis of static maps
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 20 June 1988
c Again modified by Alex Dragt, 19 August 1988
c
      use parallel, only : idproc
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
c
      include 'files.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6)
      dimension tm(6,6),t1m(6,6),t2m(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
c temporary local arrays
      dimension temp1a(monoms),temp2a(monoms),temp3a(monoms)
      dimension temp1m(6,6),temp2m(6,6),temp3m(6,6)
      dimension bm1(6,6),bm2(6,6),bm3(6,6)
      dimension rm1(6,6),rm2(6,6),rm3(6,6)
      dimension em3(6,6)
      dimension pm1(6,6),pm2(6,6),pm3(6,6)
      dimension um1(6,6)
      dimension dm1(6,6)
c
c set up control indices
      iopt=nint(p(1))
      delta=p(2)
      idata=nint(p(3))
      ipmaps=nint(p(4))
      isend=nint(p(5))
      iwmaps=nint(p(6))
c
c write headings
      if (isend.eq.1.or.isend.eq.3) then
      write(jof,*)
      write(jof,*) 'twiss analysis of static map'
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'twiss analysis of static map'
      endif
c
c first compute closed orbit to get dispersion functions
c this routine does not put out dispersions (subroutine cod does)
c but does put them in common /fitdat/ array for possible fitting
c or plotting
      call fxpt(fa,fm,temp1a,temp1m,temp3a,temp3m)
      dz(1)=temp3m(1,6)
      dz(2)=temp3m(2,6)
      dz(3)=temp3m(3,6)
      dz(4)=temp3m(4,6)
c
c preparatory steps for starting main calculation
c one objective of this calculation is to find script Ac,
c the transforming map with respect to the closed orbit
c remove offensive terms from matrix part of map:
      call clear(buf5a,buf5m)
      call spur2(fa,fm,ga,gm,t1a,t1m,buf5m)
c temporarily save the transforming map associated with sa2 in
c buffer 5 for later use
c
c compute tunes:
      cwx=gm(1,1)
      swx=gm(1,2)
      cwy=gm(3,3)
      swy=gm(3,4)
      pi=4.*atan(1.d0)
      wx=atan2(swx,cwx)
      if (wx.lt.0.) wx=wx+2.*pi
      wy=atan2(swy,cwy)
      if (wy.lt.0.) wy=wy+2.*pi
      tx=wx/(2.*pi)
      ty=wy/(2.*pi)
c put results in commom/fitdat/ array
      tux=tx
      tuy=ty
c
c preparatory steps for continuing calculation
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,buf5a,buf5m,ta,tm)
c resonance decompose purified map:
      call ctosr(g1a,t2a)
c
c compute chromaticities
c
c procedure when IOPT = 1
      if (iopt.eq.1) then
c compute first order chromaticities:
      chrox1=-(1.d0/pi)*t2a(28)
      chroy1=-(1.d0/pi)*t2a(29)
c compute second order chromaticities:
      chrox2=-(2.d0/pi)*t2a(84)
      chroy2=-(2.d0/pi)*t2a(85)
c put results in commom/fitdat/ array
      cx=chrox1
      cy=chroy1
      qx=chrox2
      qy=chroy2
c
c compute tunes about closed orbit
      delta2=delta*delta
      txc=tx+delta*chrox1+delta2*chrox2
      tyc=ty+delta*chroy1+delta2*chroy2
      tsc=txc-tyc
c
      endif
c
c procedure when IOPT = 2
      if (iopt.eq.2) then
c compute first order chromaticities:
      chrox1=(beta/pi)*t2a(28)
      chroy1=(beta/pi)*t2a(29)
c compute second order chromaticities:
      beta2=beta*beta
      beta3=beta*beta2
      chrox2=(t2a(28)*(beta-beta3)-2.d0*t2a(84)*beta2)/pi
      chroy2=(t2a(29)*(beta-beta3)-2.d0*t2a(85)*beta2)/pi
c put results in commom/fitdat/ array
      cx=chrox1
      cy=chroy1
      qx=chrox2
      qy=chroy2
c
c compute tune about closed orbit
      delta2=delta*delta
      txc=tx+delta*chrox1+delta2*chrox2
      tyc=ty+delta*chroy1+delta2*chroy2
      tsc=txc-tyc
c
      endif
c
c write out tunes and chromaticities
      if (idata.eq.1 .or. idata.eq.12 .or.
     # idata.eq.13 .or. idata.eq.123) then
      if(isend.eq.0) goto 11
      do 10 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 10
      write (ifile,*)
      if(iopt.eq.1) then
      write (ifile,*) 'tunes and chromaticities for delta defined in',
     #' terms of P sub tau:'
      endif
      if(iopt.eq.2) then
      write (ifile,*) 'tunes and chromaticities for delta defined in',
     #' terms of momentum deviation:'
      endif
      write (ifile,*)
      write (ifile,*) 'horizontal tune =',tx
      write (ifile,*) 'first order horizontal chromaticity =',chrox1
      write (ifile,*) 'second order horizontal chromaticity =',chrox2
      write (ifile,*) 'horizontal tune when delta =',delta
      write (ifile,*) txc
      write (ifile,*)
      write (ifile,*) 'vertical tune =',ty
      write (ifile,*) 'first order vertical chromaticity =',chroy1
      write (ifile,*) 'second order vertical chromaticity =',chroy2
      write (ifile,*) 'vertical tune when delta =',delta
      write (ifile,*) tyc
      write (ifile,*)
      write (ifile,*) 'tune separation when delta=',delta
      write (ifile,*) tsc
   10 continue
   11 continue
      endif
c
c preparatory steps for continuing calculation
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t2a,t2m)
c accumulate transforming map
      call concat(t2a,t2m,ta,tm,t1a,t1m)
c resonance decompose purified map:
      call ctosr(ga,t2a)
c
c compute dependence of tune on betatron amplitude:
      hhi=t2a(87)
      vvi=t2a(88)
      hvi=t2a(89)
      hhn=-2.d0*hhi/pi
      vvn=-2.d0*vvi/pi
      hvn=-hvi/pi
c put results in commom/fitdat/ array
      hh=hhn
      vv=vvn
      hv=hvn
c
c write out normalized anharmonicities
      if (idata.eq.1 .or. idata.eq.12 .or.
     # idata.eq.13 .or. idata.eq.123) then
      if(isend.eq.0) goto 21
      do 20 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 20
      write(ifile,*)
      write(ifile,*) 'normalized anharmonicities'
      write(ifile,*) ' hhn=',hhn
      write(ifile,*) ' vvn=',vvn
      write(ifile,*) ' hvn=',hvn
   20 continue
   21 continue
      endif
c
c complete computation of script Ac and script N
c store script Ac in buffer 1 and script N in buffer 2
      call spur4(ga,gm,buf2a,buf2m,t2a,t2m)
c accumulate transforming map to get script Ac
      call concat(t2a,t2m,t1a,t1m,buf1a,buf1m)
c
c compute twiss parameter expansions (invariants) and envelopes
c
c preparatory steps for continuing calculation
c compute fixed point and map around it, and put the transforming
c map to the closed orbit in buffer 3
      call fxpt(fa,fm,g1a,g1m,buf3a,buf3m)
c extract the betatron factor of the map and store it in buffer 4
      call betmap(g1a,g1m,buf4a,buf4m)
c find the transforming (conjugating) map script Ab for the betatron factor
c use the map in buffer 5 to purify the f2 part of the betatron factor
      call sndwch(buf5a,buf5m,buf4a,buf4m,g1a,g1m)
c remove offensive chromatic terms from f3 part of betatron factor
      call scpur3(g1a,g1m,ga,gm,ta,tm)
c accumulate transforming map:
      call concat(ta,tm,buf5a,buf5m,t2a,t2m)
c remove offensive terms from f4 part of betatron factor
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map to get script Ab:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c store script Ab in buffer 5
      call mapmap(ta,tm,buf5a,buf5m)
c invert script Ab
      call inv(ta,tm)
c
c compute invariants
c computation of x invariant
      call clear(ga,gm)
      ga(7)=1.d0
      ga(13)=1.d0
      call fxform(ta,tm,ga,t1a)
c computation of y invariant
      call clear(ga,gm)
      ga(18)=1.d0
      ga(22)=1.d0
      call fxform(ta,tm,ga,t2a)
c
c preliminary calculations required for envelopes and eigenvalues
      if (idata.eq.2 .or. idata.eq.3 .or.
     # idata.eq.12 .or. idata.eq.13 .or.
     # idata.eq.23 .or. idata.eq.123) then
c put script Ab in ta,tm
      call mapmap(buf5a,buf5m,ta,tm)
c make chromatic expansion of ta,tm
c the result will be used to compute both envelopes and eigenvectors
      call chrexp(iopt,delta,ta,tm,am1,am2,am3)
      endif
c
c see if output of twiss parameters and envelopes is desired
      if (idata.eq.2 .or. idata.eq.12 .or.
     # idata.eq.23 .or. idata.eq.123) then
c
c continue with calculation
c
c terms for horizontal (x) plane
c 'diagonal terms'
      ax0=t1a(8)/2.d0
      bx0=t1a(13)
      gx0=t1a(7)
c all terms: later print t1a as a map
c terms for vertical (y) plane
c 'diagonal'terms
      ay0=t2a(19)/2.d0
      by0=t2a(22)
      gy0=t2a(18)
c all terms: later print t2a as a map
c
c put horizontal and vertical results in commom/fitdat/ array
      ax=ax0
      bx=bx0
      gx=gx0
      ay=ay0
      by=by0
      gy=gy0
c
c compute envelopes
      exhc2=am3(1,1)**2+am3(1,2)**2
      exvc2=am3(1,3)**2+am3(1,4)**2
      epxhc2=am3(2,1)**2+am3(2,2)**2
      epxvc2=am3(2,3)**2+am3(2,4)**2
      eyhc2=am3(3,1)**2+am3(3,2)**2
      eyvc2=am3(3,3)**2+am3(3,4)**2
      epyhc2=am3(4,1)**2+am3(4,2)**2
      epyvc2=am3(4,3)**2+am3(4,4)**2
      exhc=sqrt(exhc2)
      exvc=sqrt(exvc2)
      epxhc=sqrt(epxhc2)
      epxvc=sqrt(epxvc2)
      eyhc=sqrt(eyhc2)
      eyvc=sqrt(eyvc2)
      epyhc=sqrt(epyhc2)
      epyvc=sqrt(epyvc2)
c
c write out twiss functions invariants, and envelopes
      if(isend.eq.0) goto 31
      do 30 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 30
      write (ifile,*)
      write (ifile,*) 'twiss parameters, invariants, and envelopes'
c
c write twiss functions and invariants
      write (ifile,*)
      write (ifile,*) 'horizontal parameters'
      write (ifile,*) 'on energy diagonal terms (alpha,beta,gamma)'
      write (ifile,*) ax0,bx0,gx0
      write (ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,t1a,t1m)
      write (ifile,*)
      write (ifile,*) 'vertical parameters'
      write (ifile,*) 'on energy diagonal terms (alpha,beta,gamma)'
      write (ifile,*) ay0,by0,gy0
      write (ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,t2a,t2m)
c
c write out envelopes
      if(iopt.eq.1) then
      write (ifile,*)
      write (ifile,*) 'envelopes for delta defined in terms of',
     #' P sub tau with delta =',delta
      endif
      if(iopt.eq.2) then
      write (ifile,*) 'envelopes for delta defined in terms of',
     #' momentum deviation with delta =',delta
      endif
      write (ifile,*)
      write (ifile,*) 'normalized horizontal envelope coefficients'
     #,' (exhc,exvc;epxhc,epxvc)'
      write (ifile,*) exhc,exvc
      write (ifile,*) epxhc,epxvc
      write (ifile,*)
      write (ifile,*) 'normalized vertical envelope coefficients'
     #,' (eyhc,eyvc;epyhc,epyvc)'
      write (ifile,*) eyhc,eyvc
      write (ifile,*) epyhc,epyvc
c
   30 continue
   31 continue
      endif
c
c Procedure for output of eigenvectors.
      if (idata.eq.3 .or. idata.eq.13 .or.
     # idata.eq.23 .or. idata.eq.123) then
c Print out matrices tm, am1, and am2.
      if(isend.eq.0) goto 41
      do 40 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 40
c
c procedure when IOPT = 1
      if (iopt.eq.1) then
      write(ifile,198)
  198 format(//,1x,'eigenvector expansion for delta defined',
     #1x,'in terms of P sub tau:')
      write(ifile,200)
  200 format(/,1x,'on energy matrix of eigenvectors')
      endif
c
c procedure when IOPT = 2
      if (iopt.eq.2) then
      write(ifile,199)
  199 format(//,1x,'eigenvector expansion for delta defined',
     #1x,'in terms of momentum deviation:')
      write(ifile,201)
  201 format(/,1x,'on momentum matrix of eigenvectors')
      endif
c
      call pcmap(i,0,0,0,ta,tm)
      write(ifile,300)
  300 format(//,1x,'delta correction')
      call pcmap(i,0,0,0,ta,am1)
      write(ifile,400)
  400 format(//,1x,'delta**2 correction')
      call pcmap(i,0,0,0,ta,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
  402 format(//,1x,'matrix of eigenvectors when delta= ',d15.8)
      call pcmap(i,0,0,0,ta,am3)
c
c test results
c
      write(ifile,*) 'test results'
**************************************************************************
***************************************************************************
c compute matrix for betatron portion of map
      call chrexp(iopt,delta,buf4a,buf4m,bm1,bm2,bm3)
c compute tune matrix
      call spur2(fa,fm,ga,gm,t1a,t1m,t2m)
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
      call clear(temp1a,temp1m)
      call clear(temp3a,temp3m)
      call matmat(gm,temp1m)
      call ctosr(g1a,temp2a)
      temp3a(28)=temp2a(28)
      temp3a(29)=temp2a(29)
      temp3a(84)=temp2a(84)
      temp3a(85)=temp2a(85)
      call srtoc(temp3a,temp1a)
      call chrexp(iopt,delta,temp1a,temp1m,rm1,rm2,rm3)
c set up matrix of eigenvectors
      call matmat(am3,em3)
c form the product pm1=em3*rm3
      call mmult(em3,rm3,pm1)
c invert pm1
      call inv(t1a,pm1)
c form the product pm2=bm3*em3
      call mmult(bm3,em3,pm2)
c form the product pm3=(pm1 inverse)*pm2
      call mmult(pm1,pm2,pm3)
c form the negative identity matrix
      call ident(temp1a,um1)
      call smmult(-1.d0,um1,um1)
c form the difference dm1=pm3-um1
      call madd(pm3,um1,dm1)
c print the results bm3 and dm1
      call pcmap(i,0,0,0,ta,bm3)
      call pcmap(i,0,0,0,ta,dm1)
**************************************************************
**************************************************************
c
   40 continue
   41 continue
      endif
c
c Procedure for printing of maps.
c
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
c print out script Ac and script N
      if(isend.eq.0) goto 51
      do 50 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 50
      write(ifile,600)
  600 format(//,1x,'transforming map with respect to
     # the closed orbit')
      call pcmap(i,i,0,0,buf1a,buf1m)
      write(ifile,700)
  700 format(//,1x,'normal form for transfer map')
      call pcmap(i,i,0,0,buf2a,buf2m)
   50 continue
   51 continue
      endif
c
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
c print out betatron portion of map and script Ab
      if(isend.eq.0) goto 61
      do 60 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 60
      write(ifile,610)
  610 format(//,1x,'betatron factor of transfer map')
      call pcmap(i,i,0,0,buf4a,buf4m)
      write(ifile,710)
  710 format(//,1x,'transforming map for betatron factor')
      call pcmap(i,i,0,0,buf5a,buf5m)
   60 continue
   61 continue
      endif
c
c Procedure for writing of maps.
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      call mapout(0,buf4a,buf4m)
      call mapout(0,buf5a,buf5m)
      mpo=mpot
      endif
c
      return
      end
***********************************************************************
c
      subroutine tbas(p,fa,fm)
c this routine translates bases
c Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),fa(monoms),fm(6,6)
      dimension ga(monoms),gm(6,6)
c
      call mapmap(fa,fm,ga,gm)
      iopt=nint(p(1))
      if (iopt.eq.1) call ctosr(ga,fa)
      if (iopt.eq.2) call ctodr(ga,fa)
      if (iopt.eq.3) call srtoc(ga,fa)
      if (iopt.eq.4) call drtoc(ga,fa)
      return
      end
c
*******************************************************************
c
      subroutine trda(p,fa,fm)
c  this subroutine transports a dynamic script A
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'trda not yet available'
      return
      end
c
*******************************************************************
c
      subroutine trsa(p,fa,fm)
c  this subroutine transports a static script A
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'trsa not yet available'
      return
      end
c
c  end of file
