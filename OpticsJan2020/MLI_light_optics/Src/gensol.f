************************************************************************
* header                    GENSOL                                     *
*         (GENMAP for a solenoid magnet with soft fringe fields)       *
*  All routines needed for this special GENMAP                         *
************************************************************************
c
      subroutine gensol(p,fa,fm,jsl,nsl,slfr)
c
c This routine computes the map for a solenoid, by numerical integration.
c F. Neri, 8/18/89; A. Dragt, 10/4/89
c The routine is based on Rob Ryne original solnsc, but all the code
c has been rewritten.
c Modified by D.T. Abell (15.Jan.07) to include 5th and 6th-order terms
c and enable slicing of solenoids.
c
      use beamdata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'files.inc'
      include 'sol.inc'
c
c  calling arrays
      dimension p(*)
      dimension fa(monoms), fm(6,6)
c
c  local arrays
      dimension y(monoms+15)
c
c use equivalence statement to make the various parameter sets pstj
c available as if they were in a two dimensional array
c
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c  y(225-476) = f5
c  y(477-938) = f6
c
c  get interval and number of steps from GENSOL parameters
c
      zi = p(1)
      zf = p(2)
      ns = nint(p(3))
      ifile = nint(p(4))
      ips = nint(p(5))
      mpole = nint(p(6))  ! not used
c
      if(ips.ge.1 .and. ips.le.9)then
c  get other parameters from pset
        di =  pst(1,ips)
        tl =  pst(2,ips)
        cl =  pst(3,ips)
        bz0 = pst(4,ips)
        iecho = nint(pst(5,ips))
        ioptr = nint(pst(6,ips))
c
      elseif(ips.eq.-1)then
c  If using sif-style input (type code 'solenoid' instead of 'sol')
c  then the code will have set p(5)=-1.
c  In that case, the other parameters are in array elements p(7)-p(12)
        di =  p(7)
        tl =  p(8)
        cl =  p(9)
        bz0 = p(10)
        iecho = nint(p(11))
        ioptr = nint(p(12))
      else
        write(6,*)'error in solenoid specification'
        write(6,*)'problem with element 5 of array passed to gensol'
        stop
      endif
c
c dabell Mon Jan 15 09:48:01 PST 2007
c MaryLie manual defines di as length from zi to start of solenoid body
c but code treats di as z at start of solenoid body; so here we redefine
c di as the code expects.
      di = zi + di
c
c complain if cl equals zero
      if (cl.eq.0.d0) then
        if (idproc.eq.0) then
          write(6,*) ' <*** ERROR ***>  solenoid specified with'
          write(6,*) '  characteristic length of zero!'
        end if
        call myexit()
      end if
c
c set parameters for this slice
      if (nsl.gt.1) then
        dz=slfr*(zf-zi)
        z1=zi+dz*(jsl-1)
        z2=zi+dz*jsl
        nstep=nint(slfr*ns)
        dz=dz/real(nstep)
      else
        dz = (zf - zi) / real(ns)
        z1 = zi
        z2 = zf
        nstep = ns
      end if
c
c
c  write out parameters, if desired
c
      if (iecho.eq.1.or.iecho.eq.3) then
        if (jsl.eq.1) then
          write(jof,*)
          write(jof,*) ' zi=',zi,' zf=',zf,' nsl=',nsl
          write(jof,*) ' ns=',ns
          write(jof,*) ' di=',di,' length=',tl
          write(jof,*) ' cl=',cl,' B=',bz0
          write(jof,*) ' iopt=',ioptr
          write(jof,*)
        end if
        if (nsl.ne.1) then
          write(jof,*) ' z1=',z1,' z2=',z2,' dz=',dz
          write(jof,*) ' jsl=',jsl,' nstep=',nstep
          write(jof,*)
        end if
      end if
      if (iecho.eq.2.or.iecho.eq.3) then
        if (jsl.eq.1) then
          write(jodf,*)
          write(jodf,*) ' zi=',zi,' zf=',zf,' nsl=',nsl
          write(jodf,*) ' ns=',ns
          write(jodf,*) ' di=',di,' length=',tl
          write(jodf,*) ' cl=',cl,' B=',bz0
          write(jodf,*) ' iopt=',ioptr
          write(jodf,*)
        end if
        if (nsl.ne.1) then
          write(jodf,*) ' z1=',z1,' z2=',z2,' dz=',dz
          write(jodf,*) ' jsl=',jsl,' nstep=',nstep
          write(jodf,*)
        end if
      end if
c
c
c  write out gradient and derivatives on file ifile:
      ipflag=0
      if (ifile.ne.0) then
        if (ifile.lt.0) then
          ipflag=1
          ifile=-ifile
        endif
        zz = zi
        h = (zf - zi) / real(ns)
        do 999 ii = 0, ns
          call bz04(zz,b0,b2,b4)
          write(ifile,137) zz, b0, b2, b4, 0., 0.
 137      format(6(1x,1pg12.5))
          zz = zz + h
 999    continue
        write(jof,*) ' profile written on file ',ifile
      endif
c
c  return identity map if ifile was < 0
c
      if (ipflag .eq. 1) then
        call ident(fa,fm)
        return
      endif
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=  0.d0
      y(2)=  0.d0
      y(3)=  0.d0
      y(4)=  0.d0
      y(5)=  0.d0       ! updated after integration (in afro.f)
      y(6)= -1.d0/beta  ! for static units (change to dynamic in afro)
c  set constants
      qbyp=  1.d0/brho
      ptg=  -1.d0/beta
c
c  initialize map to the identity:
c
      ne=monoms+15
      do 40 i=7,ne
   40   y(i)=0.d0
      do 50 i=1,6
        j=7*i
   50   y(j)=1.d0
c
c  do the computation:
      t=z1
      iflag = 4
cryne 1 August 2004 fix later:
cryne call adam11(h,ns,'start',t,y)
      call adam11(dz,nstep,'start',t,y,ne)
      call putmap(y,fa,fm)
      call csym(1,fm,ans)
c
      return
      end
c
*************************************************************************
c
      subroutine bz02(z,b0,b2)
c  This routine computes b0(z) = Bz(z)
c  and the second derivative b2(z) on axis
c  Alex Dragt 10/4/89
c
      include 'impli.inc'
      include 'sol.inc'
c
      zz = z - di
      call bump0(zz,cl,tl,ans0)
      call bump2(zz,cl,tl,ans2)
      b0 = bz0*ans0
      b2 = bz0*ans2
c
      return
      end
c
*************************************************************************
      subroutine bz04(z,b0,b2,b4)
c
c Return the zeroth, second, and fourth derivatives of the on-axis
c solenoidal field described by the parameters in 'sol.inc'.
c Dan Abell, January 2007
      implicit none
      double precision, intent(in) :: z
      double precision, intent(out) :: b0,b2,b4
c-----!----------------------------------------------------------------!
      include 'sol.inc'
      double precision :: ta,ta2,tb,tb2,zz
c
      zz = z - di
      ta = tanh(zz/cl)
      tb = tanh((zz-tl)/cl)
      ta2=ta**2
      tb2=tb**2
      b0 =  bz0 * 0.5d0 * (ta - tb)
      b2 = -bz0 * (ta * (1.d0 - ta2) - tb * (1.d0 - tb2)) / (cl**2)
      b4 =  bz0 * 4.d0 * (ta * (1.d0 - ta2) * (2.d0 - 3.d0 * ta2)       &
     &                    - tb * (1.d0 - tb2) * (2.d0 - 3.d0 * tb2))    &
     &                   / (cl**4)
c
      return
      end
c
**********************************************************************
c
      subroutine bump0(z,cl,tl,ans0)
c
c This routine computes the soft-edge bump function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
c-----------------------------------------------------------
      sgn0(z,cl)=tanh(z/cl)
c-----------------------------------------------------------
      ans0=( sgn0(z,cl) - sgn0(z-tl,cl) )/2.
c
      return
      end
c
************************************************************************
c
      subroutine bump2(z,cl,tl,ans2)
c
c This routine computes the second derivative of the soft-edge bump function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
      zl=z
      zr=z-tl
      call sgn2(zl,cl,ansl)
      call sgn2(zr,cl,ansr)
      ans2=(ansl - ansr)/2.
c
      return
      end
c
*********************************************************************
c
      subroutine sgn2(z,cl,ans)
c
c This subroutine computes the second derivative of the approximating
c signum function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
      ans=0.
      if( abs(z/cl) .lt. 30.) then
      ans=-(2./(cl**2))*(tanh(z/cl))/((cosh(z/cl))**2)
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine hmltn4(t,y,h)
c
c  This routine is used to specify the Hamiltonian h for a solenoid.
c  Written by A. Dragt 9/28/89
c  Modified by D.T. Abell to include 5th and 6th-order terms (15.Jan.07)
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'sol.inc'
c
c calling arrays
      dimension h(monoms)
      dimension y(monoms+15)
c
c  begin calculation
c
c  compute gradients
      call bz04(t,b0,b2,b4)
c
c scale gradients
c
      b0=sl*b0/brho
      b2=(sl**3)*b2/brho
c
c compute terms in hamiltonian
c
      h=0.d0
c
c terms of degree 2
c
      h( 7)= b0**2/(8.*sl)
      h(10)= -b0/(2.*sl)
      h(13)= 1./(2.*sl)
      h(14)= b0/(2.*sl)
      h(18)= b0**2/(8.*sl)
      h(22)= 1./(2.*sl)
      h(27)= 1./(2.*beta**2*gamma**2*sl)
      if (ioptr .ne. 0) then
        h(10)= 0.d0
        h(14)= 0.d0
      endif
c
c terms of degree 3
c
      h(33)= b0**2/(8.*beta*sl)
      h(45)= -b0/(2.*beta*sl)
      h(53)= 1./(2.*beta*sl)
      h(57)= b0/(2.*beta*sl)
      h(67)= b0**2/(8.*beta*sl)
      h(76)= 1./(2.*beta*sl)
      h(83)= 1./(2.*beta**3*gamma**2*sl)
c
c terms of degree 4
c
      h( 84)= b0*(b0**3 - 4.*b2)/(128.*sl)
      h( 87)= -(b0**3 - b2)/(16.*sl)
      h( 90)= b0**2/(16.*sl)
      h( 91)= (b0**3 - b2)/(16.*sl)
      h( 95)= b0*(b0**3 - 4.*b2)/(64.*sl)
      h( 99)= 3.*b0**2/(16.*sl)
      h(104)= b0**2*(3. - beta**2)/(16.*beta**2*sl)
      h(107)= -b0/(4.*sl)
      h(111)= -(b0**2)/(4.*sl)
      h(121)= -(b0**3 - b2)/(16.*sl)
      h(130)= -b0/(4.*sl)
      h(135)= -b0*(3. - beta**2)/(4.*beta**2*sl)
      h(140)= 1./(8.*sl)
      h(141)= b0/(4.*sl)
      h(145)= 3.*b0**2/(16.*sl)
      h(149)= 1./(4.*sl)
      h(154)= (3. - beta**2)/(4.*beta**2*sl)
      h(155)= (b0**3 - b2)/(16.*sl)
      h(159)= b0/(4.*sl)
      h(164)= b0*(3. - beta**2)/(4.*beta**2*sl)
      h(175)= b0*(b0**3 - 4.*b2)/(128.*sl)
      h(179)= b0**2/(16.*sl)
      h(184)= b0**2*(3. - beta**2)/(16.*beta**2*sl)
      h(195)= 1./(8.*sl)
      h(200)= (3. - beta**2)/(4.*beta**2*sl)
      h(209)= (5. - beta**2)/(8.*beta**4*gamma**2*sl)
c
c terms of degree 5
c
      h(215)= b0*(3.*b0**3 - 4.*b2)/(128.*beta*sl)
      h(227)= -(3.*b0**3 - b2)/(16.*beta*sl)
      h(235)= 3.*b0**2/(16.*beta*sl)
      h(239)= (3.*b0**3 - b2)/(16.*beta*sl)
      h(249)= b0*(3.*b0**3 - 4.*b2)/(64.*beta*sl)
      h(258)= 9.*b0**2/(16.*beta*sl)
      h(265)= b0**2*(5. - 3.*beta**2)/(16.*beta**3*sl)
      h(277)= -3.*b0/(4.*beta*sl)
      h(287)= -3.*b0**2/(4.*beta*sl)
      h(307)= -(3.*b0**3 - b2)/(16.*beta*sl)
      h(323)= -3.*b0/(4.*beta*sl)
      h(330)= -b0*(5. - 3.*beta**2)/(4.*beta**3*sl)
      h(340)= 3./(8.*beta*sl)
      h(344)= 3.*b0/(4.*beta*sl)
      h(354)= 9.*b0**2/(16.*beta*sl)
      h(363)= 3./(4.*beta*sl)
      h(370)= (5. - 3.*beta**2)/(4.*beta**3*sl)
      h(374)= (3.*b0**3 - b2)/(16.*beta*sl)
      h(383)= 3.*b0/(4.*beta*sl)
      h(390)= b0*(5. - 3.*beta**2)/(4.*beta**3*sl)
      h(409)= b0*(3.*b0**3 - 4.*b2)/(128.*beta*sl)
      h(418)= 3.*b0**2/(16.*beta*sl)
      h(425)= b0**2*(5. - 3.*beta**2)/(16.*beta**3*sl)
      h(443)= 3./(8.*beta*sl)
      h(450)= (5. - 3.*beta**2)/(4.*beta**3*sl)
      h(461)= (7. - 3.*beta**2)/(8.*beta**5*gamma**2*sl)
c
c terms of degree 6
c
      h(462)= (3.*b0**6 - 12.*b0**3*b2 + 6.*b2**2 + 4.*b0*b4)/(3072.*sl)
      h(465)= (-9.*b0**5 + 18.*b0**2*b2 - 2.*b4)/(768.*sl)
      h(468)= b0*(3.*b0**3 - 4.*b2)/(256.*sl)
      h(469)= (9.*b0**5 - 18.*b0**2*b2 + 2.*b4)/(768.*sl)
      h(473)= (3.*b0**6 - 12.*b0**3*b2 + 6.*b2**2 + 4.*b0*b4)/(1024.*sl)
      h(477)= b0*(15.*b0**3 - 12.*b2)/(256.*sl)
      h(482)= b0*(b0**3*(15. - 3.*beta**2) - b2*(12. - 4.*beta**2))     &
     &        /(256.*beta**2*sl)
      h(485)= -(3.*b0**3 - b2)/(32.*sl)
      h(489)= -b0*(3.*b0**3 - 2.*b2)/(32.*sl)
      h(499)= -(9.*b0**5 - 18.*b0**2*b2 + 2.*b4)/(384.*sl)
      h(508)= -(5.*b0**3 - b2)/(32.*sl)
      h(513)= -(b0**3*(15. - 3.*beta**2) - b2*(3. - beta**2))           &
     &        /(32.*beta**2*sl)
      h(518)= 3.*b0**2/(64.*sl)
      h(519)= (3.*b0**3 - b2)/(32.*sl)
      h(523)= b0*(9.*b0**3 - 8.*b2)/(128.*sl)
      h(527)= 9.*b0**2/(32.*sl)
      h(532)= 3.*b0**2*(5. - beta**2)/(32.*beta**2*sl)
      h(533)= (9.*b0**5 - 18.*b0**2*b2 + 2.*b4)/(384.*sl)
      h(537)= (9.*b0**3 - b2)/(32.*sl)
      h(542)= (b0**3*(15. - 3.*beta**2) - b2*(3. - beta**2))            &
     &        /(32.*beta**2*sl)
      h(553)= (3.*b0**6 - 12.*b0**3*b2 + 6.*b2**2 + 4.*b0*b4)/(1024.*sl)
      h(557)= b0*(9.*b0**3 - 8.*b2)/(128.*sl)
      h(562)= b0*(b0**3*(15. - 3.*beta**2) - b2*(12. - 4.*beta**2))     &
     &        /(128.*beta**2*sl)
      h(573)= 15.*b0**2/(64.*sl)
      h(578)= b0**2*(45. - 9.*beta**2)/(32.*beta**2*sl)
      h(587)= b0**2*(35. - 30.*beta**2 + 3.*beta**4)/(64.*beta**4*sl)
      h(590)= -3.*b0/(16.*sl)
      h(594)= -3.*b0**2/(8.*sl)
      h(604)= -(9.*b0**3 - b2)/(32.*sl)
      h(613)= -3.*b0/(8.*sl)
      h(618)= -b0*(15 - 3.*beta**2)/(8.*beta**2*sl)
      h(624)= -b0*(3.*b0**3 - 2.*b2)/(32.*sl)
      h(633)= -3.*b0**2/(8.*sl)
      h(638)= -b0**2*(15. - 3.*beta**2)/(8.*beta**2*sl)
      h(659)= -(9.*b0**5 - 18.*b0**2*b2 + 2.*b4)/(768.*sl)
      h(668)= -(3.*b0**3 - b2)/(32.*sl)
      h(673)= -(b0**3*(15. - 3.*beta**2) - b2*(3. - beta**2))           &
     &        /(32.*beta**2*sl)
      h(693)= -3.*b0/(16.*sl)
      h(698)= -b0*(15. - 3.*beta**2)/(8.*beta**2*sl)
      h(707)= -b0*(35. - 30.*beta**2 + 3.*beta**4)/(16.*beta**4*sl)
      h(714)= 1./(16.*sl)
      h(715)= 3.*b0/(16.*sl)
      h(719)= 15.*b0**2/(64.*sl)
      h(723)= 3./(16.*sl)
      h(728)= (15. - 3.*beta**2)/(16.*beta**2*sl)
      h(729)= (5.*b0**3 - b2)/(32.*sl)
      h(733)= 3.*b0/(8.*sl)
      h(738)= b0*(15. - 3.*beta**2)/(8.*beta**2*sl)
      h(749)= b0*(15.*b0**3 - 12.*b2)/(256.*sl)
      h(753)= 9.*b0**2/(32.*sl)
      h(758)= b0**2*(45. - 9.*beta**2)/(32.*beta**2*sl)
      h(769)= 3./(16.*sl)
      h(774)= (15. - 3.*beta**2)/(8.*beta**2*sl)
      h(783)= (35. - 30.*beta**2 + 3.*beta**4)/(16.*beta**4*sl)
      h(784)= (9.*b0**5 - 18.*b0**2*b2 + 2.*b4)/(768.*sl)
      h(788)= (3.*b0**3 - b2)/(32.*sl)
      h(793)= (b0**3*(15. - 3.*beta**2) - b2*(3. - beta**2))            &
     &        /(32.*beta**2*sl)
      h(804)= 3.*b0/(16.*sl)
      h(809)= b0*(15. -3.* beta**2)/(8.*beta**2*sl)
      h(818)= b0*(35. - 30.*beta**2 + 3.*beta**4)/(16.*beta**4*sl)
      h(840)= (3.*b0**6 - 12.*b0**3*b2 + 6.*b2**2 + 4.*b0*b4)/(3072.*sl)
      h(844)= b0*(3.*b0**3 - 4.*b2)/(256.*sl)
      h(849)= b0*(b0**3*(15. - 3.*beta**2) - b2*(12. - 4.*beta**2))     &
     &        /(256.*beta**2*sl)
      h(860)= 3.*b0**2/(64.*sl)
      h(865)= b0**2*(15. - 3.*beta**2)/(32.*beta**2*sl)
      h(874)= b0**2*(35. - 30.*beta**2 + 3.*beta**4)/(64.*beta**4*sl)
      h(896)= 1./(16.*sl)
      h(901)= (15. - 3.*beta**2)/(16.*beta**2*sl)
      h(910)= (35. - 30.*beta**2 + 3.*beta**4)/(16.*beta**4*sl)
      h(923)= (21. - 14.*beta**2 + beta**4)/(16.*beta**6*gamma**2*sl)
c
c add sextupoles
c      h(28)=fsxnr*sl2
c      h(30)=-3.d0*fsxsk*sl2
c      h(39)=-3.d0*fsxnr*sl2
c      h(64)=fsxsk*sl2
c
c  add octupoles
c
c      h(84)=h(84)+focnr*sl3
c      h(86)=-4.0d0*focsk*sl3
c      h(95)=-6.d0*focnr*sl3
c      h(120)=4.d0*focsk*sl3
c      h(175)=h(175)+focnr*sl3
c
      return
      end
c
c end of file
