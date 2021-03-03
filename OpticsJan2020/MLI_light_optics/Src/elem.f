***********************************************************************
* header:            ELEMENT LIBRARY                                  *
*  Common beamline elements                                           *
***********************************************************************
c
      subroutine arot(ang,h,mh)
c
c  Rotates axes counterclockwise in x-y plane by angle 'ang',
c  looking in the direction of the beam.
c  In order to get the map for an element, e.g., a quad,
c  rotated on its axis by theta clockwise looking in the direction
c  of the beam, the element map should be preceded a
c  by arot(theta) and followed by arot(-theta).
c  Written by Liam Healy, June 12, 1984.
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision h(monoms),mh(6,6)
c
      call clear(h,mh)
c  Rotate coordinates
      mh(1,1)=cos(ang)
      mh(1,3)=-sin(ang)
      mh(3,1)=sin(ang)
      mh(3,3)=cos(ang)
c  Rotate momenta
      mh(2,2)=cos(ang)
      mh(2,4)=-sin(ang)
      mh(4,2)=sin(ang)
      mh(4,4)=cos(ang)
c  Don't touch flight time
      mh(5,5)=1.
      mh(6,6)=1.
c  Polynomials are zero (bless those linear maps).
      return
      end
c
***********************************************************************
c
      subroutine jmap(h,mh)
c  Creates a map consisting of the matrix J (used in the definition of
c  symplectic matrices), with polynomials = zero.
c  Written by Liam Healy, April 16, 1985.
c
c----Variables----
      include 'impli.inc'
      include 'symp.inc'
      double precision h(*),mh(6,6)
c
c----Routine----
      call clear(h,mh)
      do 100 i=1,6
        do 100 j=1,6
          mh(i,j)=jm(i,j)
  100 continue
      return
      end
c
**********************************************************************
c
      subroutine dquad(qlength,qgrad,h,hm)
c   
c subroutine for a horizontally defocussing quad
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include  'files.inc'
      dimension h(monoms),hm(6,6)
c
c change sign of gradient
      grad=-qgrad
      call sssquad(qlength,grad,h,hm)
c
      return
      end
c
***********************************************************************
c
      subroutine sssquad(qlength,qgrad,h,hm)
c
c     Computes matrix hm and polynomial array
c     h of a quadrupole of length qlength (meters) and with field
c     gradient qgrad (tesla/meter) using the SSS algorithm.
c     If qgrad is positive the quad is focusing in the horizontal
c     plane.
c
c     MARYLIE5.0 upgrade.
c     Written by M.Venturini 5 Aug 1997.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include  'files.inc'
      dimension h(monoms),hm(6,6)
c
      dimension p(6)
c
      call clear(h,hm)
c
      p(1) = qlength/sl
      p(2) =0.d0
      p(3) =0.d0
c      p(4) =.1d0
      p(4) =.05d0
c      p(5) =1.d0
      p(5) =0.d0
c      p(6) =12
      p(6) = 0.d0
c
      call hamdrift(h)
cryne 6/21/2002 modified to multiply by sl**2:
      Qbrho=qgrad/2.d0/brho*sl**2
c
      h(7)=Qbrho
      h(18)=-Qbrho
c
       write(jodf,*)'sssquad has been used'
       write(jodf,*)'eps=',p(4)
c
c        call unixtime(ti1)
       call sss(p,h,hm)
c        call unixtime(ti2)
c       t3=ti2-ti1
c       write(jof,*) 'Execution time in seconds =',t3
c       write(jof,*) '  '
c       write(jodf,*) 'Execution time in seconds =',t3
c       write(jodf,*) '  '
c
      return
      end
c
***********************************************************************
c
      subroutine drift(el,h,hm)
c
c     Generates linear matrix hm and
c     array h containing nonlinearities
c     for the transfer map describing
c     a drift section of length el meters
c     Written by A. Dragt 3 January 1995.
c     Modified by A. Dragt 12 April 1997 to add degree 5 and 6 terms.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision hm(6,6),h(monoms)
c
      call clear(h,hm)
c
      elsc=el/sl
c
c     compute matrix part
c
      do 40 k=1,6
      hm(k,k)=+1.0d0
   40 continue
      hm(1,2)=+elsc
      hm(3,4)=+elsc
      hm(5,6)=+(elsc/((gamma**2)*(beta**2)))
c
c     compute array part h
c
c     degree 3
c
      h(53)=-(elsc/(2.0d0*beta))
      h(76)=-(elsc/(2.0d0*beta))
      h(83)=-(elsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      h(140)=-elsc/8.0d0
      h(149)=-elsc/4.0d0
      h(154)=+(elsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(195)=-elsc/8.0d0
      h(200)=+(elsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(209)=+elsc*(1.0d0-(5.0d0/beta**2))/(8.0d0*gamma**2*beta**2)
c
c     degree 5
c
      h(340)=-3.d0*elsc/(8.d0*beta)
      h(363)=-3.d0*elsc/(4.d0*beta)
      h(370)=elsc*(-5.d0+3.d0*beta**2)/(4.d0*beta**3)
      h(443)=-3.d0*elsc/(8.d0*beta)
      h(450)=elsc*(-5.d0+3.d0*beta**2)/(4.d0*beta**3)
      h(461)=elsc*(-7.d0+3.d0*beta**2)/
     &(8.d0*beta**5*gamma**2)
c
c     degree 6
c
      h(714)=-elsc/(16.d0)
      h(723)=-3.d0*elsc/(16.d0)
      h(728)=3.d0*elsc*(-5.d0+beta**2)/(16.d0*beta**2)
      h(769)=-3.d0*elsc/(16.d0)
      h(774)=3.d0*elsc*(-5.d0+beta**2)/(8.d0*beta**2)
      h(783)=-elsc*(35.d0-30.d0*beta**2+3.d0*beta**4)/
     &(16.d0*beta**4)
      h(896)=-elsc/(16.d0)
      h(901)=3.d0*elsc*(-5.d0+beta**2)/(16.d0*beta**2)
      h(910)=-elsc*(35.d0-30.d0*beta**2+3.d0*beta**4)/
     &(16.d0*beta**4)
      h(923)=-elsc*(21.d0-14.d0*beta**2+beta**4)/
     &(16.d0*beta**6*gamma**2)
c
      return
      end
c
***********************************************************************
c
      subroutine gfrngg(psideg,rho,iedge,ha,hm,gap,xk1)
c
c     subroutine to generate lie transformation
c     for fringe fields of a general bending magnet.
c  "Gap" correction a la TRANSPORT added by F. Neri (5/7/89).
c
c     psi is the angle between the design
c     orbit and the normal to the pole face,
c     rho is the magnet design orbit radius in meters.
c
c     iedge=1 for leading edge transformation
c     iedge=2 for trailing edge transformation
c
c     the fringe field map is symplectic through
c     all orders, but is guaranteed to match the physical
c     map only through order 2.   that is, the
c     fringe map is a symplectic approximation
c     to the exact map, accurate through order 2.
c     Written by Liam Healy, ca 1985
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
c
      dimension ha(monoms)
      dimension hm(6,6)
c     write(6,*)'inside gfrngg with gap=',gap
c
      call clear(ha,hm)
c  set up needed quantities
      psi=psideg*pi180
c  AAARGH......................
      cpsi=dsin(psi)
      spsi=dcos(psi)
c  AAARGH......................
      cot=cpsi/spsi
      cot2=cot*cot
      csc=1.0d0/spsi
      csc2=csc*csc
      csc3=csc2*csc
      csc4=csc2*csc2
      csc5=csc2*csc3
      rhosc=rho/sl
      gsc = gap/sl
c
c     choice of leading or trailing edge
      if (iedge.gt.1) go to 1100
c
c     leading edge fringe field map
c
c     matrix array (containing linear effects)
c
      do 160 i=1,6
      hm(i,i)=1.0d0
  160 continue
      bet0 = datan(cot)
      hm(4,3)= -tan(bet0-xk1*(1.d0+cpsi**2)*csc*gsc/rhosc)/rhosc
c
c     arrays contaning generators of nonlinearities
c
c     degree 3
c
      ha(54)=-csc3/(rhosc*2.0d0)
      ha(67)=-cot*csc2/(rhosc*beta*2.0d0)
c
c     degree 4
c
      ha(145)=-3.d0*cot*csc4/(4.d0*rhosc)
      ha(158)=-cot2*csc3/(rhosc*beta)-csc5/(2.d0*rhosc*beta)
      ha(184)=-3.d0*cot*csc4/(4.d0*beta**2*rhosc)
     &+cot*csc2/(4.d0*rhosc)
c
      return
c
 1100 continue
c
c     trailing edge fringe field map
c
c     matrix array (containing linear effects)
c
      do 190 i=1,6
      hm(i,i)=+1.0d0
  190 continue
      bet0 = datan(cot)
      hm(4,3)= -tan(bet0-xk1*(1.d0+cpsi**2)*csc*gsc/rhosc)/rhosc
c
c     arrays containing nonlinearities
c
c     degree 3
c
      ha(54)=+csc3/(rhosc*2.0d0)
      ha(67)=-cot*csc2/(rhosc*beta*2.0d0)
c
c     degree 4
c
      ha(145)=-3.d0*cot*csc4/(4.d0*rhosc)
      ha(158)=+cot2*csc3/(beta*rhosc)+csc5/(2.d0*beta*rhosc)
      ha(184)=-3.d0*cot*csc4/(4.d0*beta**2*rhosc)
     &+cot*csc2/(4.d0*rhosc)
c
      return
      end
c
***********************************************************************
c
      subroutine nfrng(rho,iedge,h,mh)
c
c     generates lie transformation for fringe fields
c     of a normal entry bend with a design orbit
c     radius of rho meters
c
c     iedge=1 for leading edge
c
c     iedge=2 for trailing edge
c
c     the map here is the psi=pi/2 case
c     of the map employed in gfrng
c     Written by Alex Dragt, Fall 1986
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision mh(6,6),h(monoms),rho
c
      call clear(h,mh)
c
c     set coefficients of arrays
c
c     set mh equal to the identity
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
c
c     set coefficients of h
c
c     choose leading or trailing edge
c
      if (iedge.eq.2) go to 100
c
c     leading edge coefficients
c
c     degree 3
c
      h(54)=-(sl/(2.0d0*rho))
c
c     degree 4
c
      h(158)=-sl/(2.d0*beta*rho)
      return
c
  100 continue
c     trailing edge coefficients
c
c     degree 3
c
      h(54)=+(sl/(2.0d0*rho))
c
c     degree 4
c
      h(158)=+sl/(2.d0*beta*rho)
      return
      end
c
***********************************************************************
c
      subroutine gfrng(psideg,rho,iedge,ha,hm)
c
c     subroutine to generate lie transformation
c     for fringe fields of a general bending magnet.
c
c     psi is the angle between the design
c     orbit and the normal to the pole face,
c     rho is the magnet design orbit radius in meters.
c
c     iedge=1 for leading edge transformation
c     iedge=2 for trailing edge transformation
c
c     the fringe field map is symplectic through
c     all orders, but is guaranteed to match the physical
c     map only through order 2.   that is, the
c     fringe map is a symplectic approximation
c     to the exact map, accurate through order 2.
c     Written by Liam Healy, ca 1985
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
c
      dimension ha(monoms)
      dimension hm(6,6)
c
      call clear(ha,hm)
c  set up needed quantities
      psi=psideg*pi180
      cpsi=dsin(psi)
      spsi=dcos(psi)
      cot=cpsi/spsi
      cot2=cot*cot
      csc=1.0d0/spsi
      csc2=csc*csc
      csc3=csc2*csc
      csc4=csc2*csc2
      csc5=csc2*csc3
      rhosc=rho/sl
c
c     choice of leading or trailing edge
      if (iedge.gt.1) go to 1100
c
c     leading edge fringe field map
c
c     matrix array (containing linear effects)
c
      do 160 i=1,6
      hm(i,i)=1.0d0
  160 continue
      hm(4,3)= -cot/rhosc
c
c     arrays contaning generators of nonlinearities
c
c     degree 3
c
      ha(54)=-csc3/(rhosc*2.0d0)
      ha(67)=-cot*csc2/(rhosc*beta*2.0d0)
c
c     degree 4
c
      ha(145)=-3.d0*cot*csc4/(4.d0*rhosc)
      ha(158)=-cot2*csc3/(rhosc*beta)-csc5/(2.d0*rhosc*beta)
      ha(184)=-3.d0*cot*csc4/(4.d0*beta**2*rhosc)
     &+cot*csc2/(4.d0*rhosc)
c
      return
c
 1100 continue
c
c     trailing edge fringe field map
c
c     matrix array (containing linear effects)
c
      do 190 i=1,6
      hm(i,i)=+1.0d0
  190 continue
      hm(4,3)=-cot/rhosc
c
c     arrays containing nonlinearities
c
c     degree 3
c
      ha(54)=+csc3/(rhosc*2.0d0)
      ha(67)=-cot*csc2/(rhosc*beta*2.0d0)
c
c     degree 4
c
      ha(145)=-3.d0*cot*csc4/(4.d0*rhosc)
      ha(158)=+cot2*csc3/(beta*rhosc)+csc5/(2.d0*beta*rhosc)
      ha(184)=-3.d0*cot*csc4/(4.d0*beta**2*rhosc)
     &+cot*csc2/(4.d0*rhosc)
      return
      end
c
***********************************************************************
c
      subroutine octe(l,phi0,h,mh)
c
c     computes matrix mh and polynomial h for
c     electric octupole with length l meters and scalar
c     potential of the form phi0*(x**4-6*x**2*y**2+y**4)
c     Written by D. Douglas, ca 1982
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision l,h(monoms),mh(6,6)
      double precision lsc
c
      call clear(h,mh)
      lsc=l/sl
c
c     enter matrix elements
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
      mh(1,2)=+lsc
      mh(3,4)=+lsc
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     set coefficients of h
c
c     degree 3
c
      h(53)=-(lsc/(2.0d0*beta))
      h(76)=-(lsc/(2.0d0*beta))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      ap=((sl**4)*phi0)/(c*beta*brho)
      h(84)=+lsc*ap
      h(85)=-2.0d0*(lsc**2)*ap
      h(90)=+2.0d0*(lsc**3)*ap
      h(95)=-6.0d0*lsc*ap
      h(96)=+6.d0*lsc**2*ap
      h(99)=-(2.0d0*(lsc**3)*ap)
      h(105)=-(lsc**4)*ap
      h(110)=+6.0d0*(lsc**2)*ap
      h(111)=-(8.0d0*(lsc**3)*ap)
      h(114)=+3.0d0*(lsc**4)*ap
      h(140)=-lsc/8.0d0+(lsc**5)*ap/5.0d0
      h(145)=-(2.0d0*(lsc**3)*ap)
      h(146)=+3.0d0*(lsc**4)*ap
      h(149)=-lsc/4.0d0-(6.0d0*(lsc**5)*ap)/5.0d0
      h(154)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(175)=+lsc*ap
      h(176)=-2.0d0*(lsc**2)*ap
      h(179)=+2.0d0*(lsc**3)*ap
      h(185)=-(lsc**4)*ap
      h(195)=-lsc/8.0d0+(lsc**5)*ap/5.0d0
      h(200)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(209)=+(lsc*(1.d0-(5.d0/beta**2)))/(8.d0*gamma**2*beta**2)
      return
      end
c
***********************************************************************
c
      subroutine octm(el,gb0,h,hm)
c
c     computes matrix hm and polynomial h for
c     magnetic octupole with length l meters and vector
c     potential of the form -(gb0/4.)*(+x**4-6*x**2*y**2+y**4)
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension h(monoms),hm(6,6)
c
      call sssoctm(el,gb0,h,hm)
c
      return
      end
c
***********************************************************************
c
      subroutine sssoctm(olength,oct,h,hm)
c
c     computes matrix mh and polynomial array
c     h of a octupole of length olength (meters) and with strength
c     oct (tesla/meter^3) using the SSS algorithm.
c     The vector potential is A_z=-(oct/4)(x^4 - 6x^2*y^2 + y^4).
c
c     MARYLIE5.0 upgrade.
c     Written by M.Venturini 5 Aug 1997.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      dimension h(monoms),hm(6,6)
c
      dimension p(6)
c
      call clear(h,hm)
c
      p(1) = olength/sl
      p(2) =0.d0
      p(3) =0.d0
cryne 6/21/2002      p(4) =.1d0
      p(4) =1.d-4
      p(5) =0.d0
      p(6) =0.d0
c
      call hamdrift(h)
c
cryne 6/21/2002 modified to multiply by sl**4:
      oct4brho=oct/4.d0/brho*sl**4
c
      h(84)=oct4brho
      h(95)=-6*oct4brho
      h(175)=oct4brho
c
       write(jodf,*)'sssoctm has been used'
       write(jodf,*)'eps=',p(4)
       call sss(p,h,hm)
c
      return
      end
c
***********************************************************************
c
      subroutine pbend(rho,phideg,h,mh)
c
c     computes generators for the lie transformation
c     describing a parallel-faced bending magnet
c     subtending an angle of phi radians
c     with a design orbit radius of rho meters
c     Written by D. Douglas, ca 1982
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      double precision h(monoms),mh(6,6)
c
      call clear(h,mh)
      phi=phideg*pi180
      alpha=phi/2.0d0
      rhosc=rho/sl
      cal=dcos(alpha)
      sal=dsin(alpha)
      tan=sal/cal
c
c     enter matrix elements into mh
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
      mh(1,2)=+rhosc*tan*2.0d0
      mh(3,4)=+phi*rhosc
      mh(5,6)=-phi*rhosc+(2.0d0*rhosc*tan)/(beta**2)
c
c     add coefficients of generators of nonlinearities to h
c
c     degree 3
c
      h(53)=-(rhosc*tan)/((cal**2)*beta)
      h(76)=-(rhosc*tan)/beta
      h(83)=+(rhosc*tan)/beta-(tan+tan**3/3.d0)*rhosc/beta**3
c
c     degree 4
c
      sec=1.d0/cal
      sec2=sec*sec
      sec4=sec2*sec2
      h(140)=-rhosc*tan*sec4/4.d0
      h(149)=-rhosc*tan*sec2/2.d0
      h(154)=-3.d0*rhosc*tan*sec4/(2.d0*beta**2)
     &+rhosc*tan*sec2/2.d0
      h(195)=-rhosc*tan/4.d0
      h(200)=-rhosc*tan*sec2/(2.d0*beta**2)
     &-rhosc*tan/(2.d0*gamma**2*beta**2)
     &-rhosc*tan/(2.d0*beta**2)
      h(209)=-rhosc*tan*tan*tan*tan*tan/(4.d0*beta**4)
     &-5.d0*rhosc*tan*tan*tan/(6.d0*beta**4)
     &-5.d0*rhosc*tan/(4.d0*beta**4)
     &+rhosc*tan*tan*tan/(2.d0*beta**2)
     &+3.d0*rhosc*tan/(2.d0*beta**2)
     &-rhosc*tan/4.d0
      return
      end
c
***********************************************************************
c
      subroutine nbend(rho,phideg,h,mh)
c
c     computes generators for normal entry dipole bending
c     magnet subtending angle phi radians and having
c     design orbit radius rho meters
c     Written by D. Douglas, ca 1982
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      double precision h(monoms),mh(6,6)
c
      call clear(h,mh)
c
c     compute useful numbers
c
      rhosc=rho/sl
      phi=phideg*pi180
      cphi=dcos(phi)
      sphi=dsin(phi)
c
c     set coefficients of mh
c
      mh(1,1)=+cphi
      mh(1,2)=+(rhosc*sphi)
      mh(1,6)=-(((1.0d0-cphi)*rhosc)/beta)
      mh(2,1)=-(sphi/rhosc)
      mh(2,2)=+cphi
      mh(2,6)=-(sphi/beta)
      mh(3,3)=+1.d0
      mh(3,4)=+rhosc*phi
      mh(4,4)=+1.d0
      mh(5,1)=+sphi/beta
      mh(5,2)=+(((1.d0-cphi)*rhosc)/beta)
      mh(5,5)=+1.d0
      mh(5,6)=-(rhosc*phi)+((rhosc*sphi)/(beta**2))
      mh(6,6)=+1.d0
c
c     set coefficients of array h
c
c     degree 3
c
      sphi2=sphi*sphi
      sphi3=sphi2*sphi
      cphi2=cphi*cphi
      cphi3=cphi2*cphi
      h(28)=-(sphi3/(6.0d0*(rhosc**2)))
      h(29)=-((cphi*sphi2)/((2.0d0)*rhosc))
      h(33)=-(sphi3/(2.0d0*rhosc*beta))
      h(34)=-((sphi*cphi2)/2.0d0)
      h(38)=-((sphi2*cphi)/beta)
      h(43)=-(sphi/2.0d0)
      h(48)=-(sphi/(2.0d0*(gamma**2)*(beta**2)))
     &-(sphi3/(2.0d0*(beta**2)))
      h(49)=+(((1.0d0-cphi3)*rhosc)/6.d0)
      h(53)=-((sphi*cphi2*rhosc)/(beta*2.d0))
      h(58)=+(((1.0d0-cphi)*rhosc)/2.0d0)
      h(63)=+(((1.0d0-cphi)*rhosc)/(2.0d0*(gamma**2)*(beta**2)))
     &-((sphi2*cphi*rhosc)/(2.0d0*(beta**2)))
      h(76)=-((rhosc*sphi)/(2.0d0*beta))
      h(83)=-((rhosc*sphi)/(2.0d0*(gamma**2)*(beta**3)))
     &-((rhosc*sphi3)/(6.0d0*(beta**3)))
c
c     degree 4
c
      h(89)=-sphi3/(6.d0*beta*rhosc**2)
      h(90)=-sphi3/(8.d0*rhosc)
      h(94)=-cphi*sphi2/(2.d0*beta*rhosc)
      h(99)=-sphi3/(8.d0*rhosc)
      h(104)=-sphi3*(5.d0-beta**2)/(8.d0*beta**2*rhosc)
      h(105)=-cphi*sphi2/4.d0
      h(109)=+sphi3/(4.d0*beta)-sphi/(2.d0*beta)
      h(114)=-cphi*sphi2/4.d0
      h(119)=-cphi*sphi2/(4.d0*beta**2*gamma**2)
     &-cphi*sphi2/beta**2
      h(132)=-sphi3/(4.d0*beta)-sphi/(2.d0*beta)
      h(139)=-sphi3/(4.d0*beta**3*gamma**2)
     &-sphi3/(2.d0*beta**3)-sphi/(2.d0*beta**3*gamma**2)
      h(140)=-rhosc*sphi*cphi2/8.d0
      h(144)=-rhosc*cphi*sphi2/(1.2d1*beta)
     &+rhosc*(1.d0-cphi)/(6.d0*beta)
      h(149)=+rhosc*sphi3/8.d0-rhosc*sphi/4.d0
      h(154)=+rhosc*sphi3*(4.d0-beta**2)/(8.d0*beta**2)
     &-rhosc*sphi/(4.d0*beta**2*gamma**2)
     &-rhosc*sphi/(2.d0*beta**2)
      h(167)=-rhosc*cphi*sphi2/(4.d0*beta)
     &+rhosc*(1.d0-cphi)/(2.d0*beta)
      h(174)=-rhosc*cphi*sphi2/(4.d0*beta**3*gamma**2)
     &-rhosc*cphi*sphi2/(2.d0*beta**3)
     &+rhosc*(1.d0-cphi)/(2.d0*beta**3*gamma**2)
      h(195)=-rhosc*sphi/8.d0
      h(200)=-rhosc*sphi3/(8.d0*beta**2)
     &-rhosc*sphi/(4.d0*beta**2*gamma**2)
     &-rhosc*sphi/(2.d0*beta**2)
      h(209)=-rhosc*sphi3/(8.d0*beta**4*gamma**2)
     &-rhosc*sphi3/(6.d0*beta**4)-rhosc*sphi/(2.d0*beta**4*gamma**2)
     &-rhosc*sphi/(8.d0*beta**4*gamma**4)
      return
      end
c
***********************************************************************
c
cryne Aug 6, 2003: routine name changed to prot3
cryne The fifth order version by Johannes van Zeijts is now called 'prot'
      subroutine prot3(psideg,kind,ha,hm)
c
c     subroutine to generate lie transformation
c     for rotation of reference plane.
c     used primarily  in connection with computing
c     map for dipoles with rotated pole faces.
c
c     psideg is the rotation angle angle in degrees.
c
c     kind=1 for transition from the normal reference plane
c       to a rotated reference plane
c     kind=2 for transition from a rotated reference plane
c       to a normal reference plane
c     Written by Alex Dragt, Fall 1986
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      double precision ha(monoms),hm(6,6)
c
      dimension ha1(monoms),ha2(monoms)
      dimension hm1(6,6),hm2(6,6)
c
      call clear(ha,hm)
      call clear(ha1,hm1)
      call clear(ha2,hm2)
c
c  set up needed quantities
c
      psi=psideg*pi180
      cpsi=dsin(psi)
      spsi=dcos(psi)
      cot=cpsi/spsi
      cot2=cot*cot
      cot3=cot2*cot
      cot4=cot2*cot2
      cot5=cot2*cot3
      csc=1.0d0/spsi
      csc2=csc*csc
c
c     choice of kind
      if (kind.gt.1) go to 100
c
c     transition from normal to rotated reference plane
c
c     matrix arrays (containing linear effects)
c
      do 40 i=1,6
      hm1(i,i)=+1.0d0
   40 continue
      hm1(2,6)=-cot/beta
      hm1(5,1)=+cot/beta
      do 50 i=3,6
      hm2(i,i)=+1.0d0
   50 continue
      hm2(1,1)=hm2(1,1)+1.0d0/spsi
      hm2(2,2)=hm2(2,2)+spsi
c
c     arrays containing generators of nonlinearities
c
c     degree 3
c
      ha1(34)=-cot/2.0d0
      ha1(38)=-cot2/beta
      ha1(43)=-cot/2.0d0
      ha1(48)=-cot/(gamma**2*beta**2*2.0d0)
     &-cot3/(beta**2*2.0d0)
c
c     degree 4
c
      ha1(105)=-cot2/4.d0
      ha1(109)=-cot/(2.d0*beta)-3.d0*cot3/(4.d0*beta)
      ha1(114)=-cot2/4.d0
      ha1(119)=-cot2/beta**2-3.d0*cot4/(4.d0*beta**2)
     &-cot2/(4.d0*gamma**2*beta**2)
      ha1(132)=-cot/(2.d0*beta)-cot3/(4.d0*beta)
      ha1(139)=-cot3/(2.d0*beta**3)-cot5/(4.d0*beta**3)
     &-cot/(2.d0*gamma**2*beta**3)
     &-cot3/(4.d0*gamma**2*beta**3)
c
c     compute map
      call concat(ha1,hm1,ha2,hm2,ha,hm)
      return
c
  100 continue
c
c     transition from a rotated to a normal reference plane
c
c     matrix arrays (containing linear effects)
c
      do 60 i=1,6
      hm1(i,i)=+1.0d0
   60 continue
      hm1(2,6)=-cpsi/beta
      hm1(5,1)=+cpsi/beta
      do 70 i=3,6
      hm2(i,i)=+1.0d0
   70 continue
      hm2(1,1)=hm2(1,1)+spsi
      hm2(2,2)=hm2(2,2)+csc
c
c     arrays containing generators of nonlinearties
c
c     degree 3
c
      ha1(34)=-cot*csc/2.0d0
      ha1(43)=-cpsi/2.0d0
      ha1(48)=-cpsi/(gamma**2*beta**2*2.0d0)
c
c     degree 4
c
      ha1(105)=+cot2*csc2/4.d0
      ha1(109)=-cot*csc/(2.d0*beta)
      ha1(114)=+cot2/4.d0
      ha1(119)=+cot2/(4.d0*gamma**2*beta**2)
      ha1(132)=-cpsi/(2.d0*beta)
      ha1(139)=-cpsi/(2.d0*gamma**2*beta**3)
c
c     compute map
      call concat(ha1,hm1,ha2,hm2,ha,hm)
      return
      end
c
***********************************************************************
c
cryne Aug 6, 2003: routine name changed from myprot5 to prot
      subroutine prot(angdeg,ijkind,h,mh)
c
c     High order PROT routine.
c     Actual code generated by Johannes van Zeijts using
c     REDUCE.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
c     include 'param.inc'
c     include 'parm.inc'
      double precision l,h(monoms),mh(6,6)
c
      dimension j(6)
c
      DOUBLE PRECISION B
      DOUBLE PRECISION CO
      DOUBLE PRECISION Si
c
cryne mods to allow for leading/trailing option:
cryne ijkind=1 for normal-to-rotated, =2 for rotated-to-normal
      phideg=angdeg
      if(ijkind.eq.2)phideg=-angdeg

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
cryne this line added Aug 6, 2003:
      if(ijkind.eq.2)call inv(h,mh)
      return
      end
c
c end of file
c
***********************************************************************
c
      subroutine fquad(qlength,qgrad,h,hm)
c
c subroutine for a horizontally focussing quad
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include  'files.inc'
      dimension h(monoms),hm(6,6)
c
      grad=qgrad
      call sssquad(qlength,grad,h,hm)
c
      return
      end
c
***********************************************************************
c
      subroutine frquad(gb0,ifr,h,mh)
c
c  computes hard edge fringe field map for quads
c  Written by E. Forest, ca 1984
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision gb0,h(monoms),mh(6,6)
c
c  set up linear part as identity matrix
      do 1 i=1,6
      do 2 j=1,6
      mh(i,j)=0.d0
      if(i.ne.j) goto 2
      mh(i,j)=1.d0
 2    continue
 1    continue
c  clear polynomial array
      do 3 i=1,monoms
      h(i)=0.d0
 3    continue
      gb=gb0
c  see if fringe field is leading or trailing
      if(ifr.lt.0)gb=-gb0
c  compute nonlinear part of map
      arg=(gb*(sl**2))/brho
      h(85)=arg/12.d0
      h(176)=-arg/12.d0
      h(110)=arg/4.d0
      h(96)=-arg/4.d0
      return
      end
c
***********************************************************************
c
      subroutine gbend (rho,bndang,psideg,phideg,h,mh)
c  a dipole bending magnet with arbitrary entrance and
c  exit angles (psi and phi respectively) and arbitrary
c  bend andgle (bndang).
c  written by L. Healy, April 19, 1984.
c  modified by A. Dragt, 5 Oct 1986.
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision mh(6,6),mht1(6,6),mht2(6,6),mht3(6,6)
      double precision h(monoms),ht1(monoms),ht2(monoms),ht3(monoms)
c
c  calculate each of the 3 individual pieces that make
c  up the general bending magnet, and concatenate them:
c  gbend=hpf1*nbend*hpf2
c
c  compute hpf1 and nbend
      call hpf(rho,1,psideg,ht1,mht1)
      aldeg=bndang-psideg-phideg
      call nbend(rho,aldeg,ht2,mht2)
c  form the product hpf1*nbend
      call concat(ht1,mht1,ht2,mht2,ht3,mht3)
c  compute hpf2
      call hpf(rho,-1,phideg,ht1,mht1)
c  form the complete product hpf1*nbend*hpf2
      call concat(ht3,mht3,ht1,mht1,h,mh)
      return
      end
c
***********************************************************************
c
      subroutine hpf(rho,which,phideg,h,mh)
c  This generates matrix elements and monomial coeffs
c  for half of a parallel face bending magnet.
c  The bending angle is phi, and 'which' indicates
c  whether it is the leading half (1) or trailing half (-1).
c  Written by Liam Healy, April 19, 1984.
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      integer which
      double precision mh(6,6),h(monoms)
      include 'pie.inc'
      call clear(h,mh)
c  set trig functions and scale rho
      phi=phideg*pi180
      rhosc=rho/sl
      tn=tan(phi)
      tn2=tn*tn
      tn3=tn*tn2
      tn5=tn2*tn3
      sec=1./cos(phi)
      sec2=sec*sec
      sec3=sec*sec2
      sec4=sec2*sec2
      sec5=sec2*sec3
c  set matrix
      do 120 i=1,6
  120 mh(i,i)=1.
      mh(1,2)=rhosc*tn
      mh(1,6)=-which*rhosc*(1.-sec)/beta
      mh(3,4)=rhosc*phi
      mh(5,2)=-which*rhosc*(1.-sec)/beta
      mh(5,6)=rhosc*(tn/beta**2-phi)
c  set monomial coefficients
      h(49)=which*rhosc*(1.-sec3)/6.
      h(53)=-rhosc*tn*sec2/(2.*beta)
      h(58)=which*rhosc*(1.-sec)/2.
      h(63)=which*rhosc*(1./beta**2-1.+sec*(1.-sec2/beta**2))/2.
      h(76)=-rhosc*tn/(2.*beta)
      h(83)=-rhosc*tn*(1.-beta**2+tn2/3.)/(2.*beta**3)
      h(140)=-rhosc*tn*sec4/8.
      h(144)=which*rhosc*(sec3/3.-sec5/2.+1.d0/6.d0)/beta
      h(149)=-rhosc*tn*sec2/4.
      h(154)=rhosc*tn*sec2*(1.-3.*sec2/beta**2)/4.
      h(167)=which*rhosc*(1.-sec3)/(2.*beta)
      h(174)=-which*rhosc*((1.-sec3)/(2.*beta)-(1.-sec5)/(2.*beta**3))
      h(195)=-rhosc*tn/8.
      h(200)=-tn*rhosc*((sec2+2.)/beta**2-1.)/4.
      h(209)=rhosc*(-(2.5*tn+5.*tn3/3.+tn5/2.)/beta**4
     &      +tn*(3.+tn2)/beta**2-tn/2.)/4.
c
      return
      end
c
***********************************************************************
c
      subroutine cfbend(pa,pb,fa,fm)
c  This subroutine computes the map for a combined function bend
c  numerically by use of the GENMAP exponentiation routine.
c  It should eventually be replaced by a collection of analytic formulas.
c  Written by Alex Dragt, 30 March 1987.
c  Corrected by Alex Dragt, 9 Sept 1989, based on calculations in the
c  paper "Third Order Transfer Map for Combined Function Dipole" by
c  Dragt et al. (1989).
cryne modified by Rob Ryne 7/13/2002 to get multipole coeffcients from
cryne an array passed in the parameter list. (previously this argument
cryne was an integer that pointed to a pset)
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      include 'parset.inc'
      include 'files.inc'
c
      dimension fa(monoms)
      dimension fm(6,6)
      dimension pa(6),pb(6)
      dimension pc(6)
      dimension ha(monoms)
      dimension hm(6,6)
c
c  set up parameters and control indices
c
      phideg=pa(1)
      b=pa(2)
      ilfrn=nint(pa(3))
      itfrn=nint(pa(4))
      ijopt=nint(pa(5))
cryne 7/13/2002      ipset=nint(pa(6))
c
cryne 1 August 2004:
      myorder=nint(pa(6))
c
c  compute iopt and iecho
c
      iopt=mod(ijopt,10)
c     write(6,*)'inside cfbend; ijopt,iopt=',ijopt,iopt
c     write(6,*)'pb(1-6)=',pb(1),pb(2),pb(3),pb(4),pb(5),pb(6)
      iecho=(ijopt-iopt)/10
c      write(6,*) ' iopt=',iopt,' iecho=',iecho
c
c     compute useful numbers
c
      rho=brho/b
      rhosc=rho/sl
      phi=phideg*pi180
c     write(6,*)'phi,phideg=',phi,phideg
c
cryne 7/13/2002  get multipole values from the parameter set ipset
cryne get multipole values from the array in the argument list
c
c      if (ipset.lt.1 .or. ipset.gt.maxpst) then
c        do 50 i=1,6
c  50    pb(i) = 0.0d0
c      else
c        do 60 i=1,6
c  60    pb(i) = pst(i,ipset)
c      endif
c
c compute multipole coefficients
c procedure when iopt=1
      if (iopt.eq.1) then
c     write(6,*)'(cfbend) iopt is 1'
      bqd=pb(1)
      aqd=pb(2)
      bsex=pb(3)
      asex=pb(4)
      boct=pb(5)
      aoct=pb(6)
      endif
c procedure when iopt=2
      if (iopt.eq.2) then
c     write(6,*)'(cfbend) iopt is 2'
      tay1=pb(1)
      aqd=pb(2)
      tay2=pb(3)
      asex=pb(4)
      tay3=pb(5)
      aoct=pb(6)
      bqd=tay1
      bsex=tay2+(5.d0*tay1)/(8.d0*rho)
      boct=tay3-tay1/(48.d0*rho*rho)+(2.d0*tay2)/(3.d0*rho)
      endif
c procedure when iopt=3
      if (iopt.eq.3) then
c     write(6,*)'(cfbend) iopt is 3; brho,rho,sl=',brho,rho,sl
      tay1=brho*pb(1)
      aqd=brho*pb(2)
      tay2=brho*pb(3)
      asex=brho*pb(4)
      tay3=brho*pb(5)
      aoct=brho*pb(6)
      bqd=tay1
c     write(6,*)'pb(1),tay1,bqd=',pb(1),tay1,bqd
      bsex=tay2+(5.d0*tay1)/(8.d0*rho)
      boct=tay3-tay1/(48.d0*rho*rho)+(2.d0*tay2)/(3.d0*rho)
      endif
c
c  write out multipole and taylor coefficients if desired
c
      if (iecho .ne. 0) then
      ttay1=bqd
      ttay2=bsex-5.d0*bqd/(8.d0*rho)
      ttay3=boct+7.d0*bqd/(16.d0*rho*rho)-2.d0*bsex/(3.d0*rho)
      if (iecho .eq. 1 .or. iecho .eq. 3) then
      write (jof,137)
  137 format(/,1x,'multipole strengths bqd,aqd,bsex,asex,boct,aoct')
      write (jof,*) bqd,aqd,bsex,asex,boct,aoct
      write (jof,138)
  138 format(1x,'taylor coefficients tay1,tay2,tay3')
      write (jof,*) ttay1,ttay2,ttay3
      write (jof,*)
      endif
      if (iecho .eq. 2 .or. iecho .eq. 3) then
      write (jodf,137)
      write (jodf,*) bqd,aqd,bsex,asex,boct,aoct
      write (jodf,138)
      write (jodf,*) ttay1,ttay2,ttay3
      write (jodf,*)
      endif
      endif
c
c  scale multipole strengths and compute scaled curvature feed up terms
c
      sbqdf=1.d0/(2.d0*rhosc)
      saqd=aqd*sl*rho/brho
      sbqd=bqd*sl*rho/(2.d0*brho)
c     write(6,*)'saqd,sbqd=',saqd,sbqd
c
      sasexf=aqd*(sl**2)/(8.d0*brho)
      sbsexf=bqd*(sl**2)/(8.d0*brho)
      sasex=asex*(sl**2)*rho/(3.d0*brho)
      sbsex=bsex*(sl**2)*rho/(3.d0*brho)
c
      sscalf=-bqd*(sl**3)/(64.d0*rho*brho)
      saoctf=(asex/6.d0-aqd/(16.d0*rho))*(sl**3)/brho
      sboctf=(bsex/12.d0-bqd/(32.d0*rho))*(sl**3)/brho
      saoct=aoct*(sl**3)*rho/brho
      sboct=boct*(sl**3)*rho/(4.d0*brho)
c
c  compute the Hamiltonian
c
      call clear(ha,hm)
c
c  degree 2
c
      ha(7)=sbqdf + sbqd
      ha(9)=-saqd
      ha(12)=1.d0/beta
      ha(13)=rhosc/2.d0
      ha(18)=-sbqd
      ha(22)=rhosc/2.d0
      ha(27)=rhosc/(2.d0*(beta*gamma)**2)
c
c  degree 3
c
      ha(28)=sbsexf + sbsex
      ha(30)=-sasexf - 3.d0*sasex
      ha(34)=1.d0/2.0d0
      ha(39)=sbsexf - 3.d0*sbsex
      ha(43)=1.d0/2.d0
      ha(48)=1.d0/(2.0d0*((gamma*beta)**2))
      ha(53)=rhosc/(2.d0*beta)
      ha(64)=-sasexf+sasex
      ha(76)=rhosc/(2.d0*beta)
      ha(83)=rhosc/(2.d0*(gamma**2)*(beta**3))
c
c  degree 4
c
      ha(84)=sscalf + sboctf + sboct
      ha(86)=-saoctf - saoct
      ha(95)=2.d0*sscalf -6.d0*sboct
      ha(109)=1.d0/(2.d0*beta)
      ha(120)=-saoctf + saoct
      ha(132)=1.d0/(2.d0*beta)
      ha(139)=1.d0/(2.d0*(beta**3)*(gamma**2))
      ha(140)=rhosc/8.d0
      ha(149)=rhosc/4.d0
      ha(154)=rhosc/(4.d0*(beta**2)*(gamma**2))+rhosc/(2.d0*(beta**2))
      ha(175)=sscalf - sboctf + sboct
      ha(195)=rhosc/8.d0
      ha(200)=rhosc/(4.d0*(beta**2)*(gamma**2))+rhosc/(2.d0*(beta**2))
      ha(209)=rhosc/(2.d0*(beta**4)*(gamma**2))
     &+rhosc/(8.d0*(beta**4)*(gamma**4))
c
c  call GENMAP exponentiation routine
      pc(1)=-phi
      pc(2)=0.d0
      pc(3)=0.d0
c---------------------
c     write(6,*)'pc(1)=',pc(1)
c     write(6,*)'pc(2)=',pc(2)
c     write(6,*)'pc(3)=',pc(3)
c     write(6,*)'calling cex from routine cfbend'
c---------------------
      call cex(pc,ha,hm,myorder)
c     write(6,*)'returned from cex'
c     write(6,*)'hm='
c     do i=1,6
c     write(6,1232)(hm(i,j),j=1,6)
c1232 format(6(1pe12.5,1x))
c     enddo
      call mapmap(ha,hm,fa,fm)
c
c  fringe field effects (both dipole and quadrupole) are put on 
c  in the subroutine lmnt
c
      return
      end
c
***********************************************************************
c
      subroutine cfbend_old(pa,fa,fm)
c  This subroutine computes the map for a combined function bend
c  numerically by use of the GENMAP exponentiation routine.
c  It should eventually be replaced by a collection of analytic formulas
c  Written by Alex Dragt, 30 March 1987.
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension pa(6),pb(6)
      dimension pc(6)
      dimension ha(monoms)
      dimension hm(6,6)
c
      include 'pie.inc'
      include 'parset.inc'
c
c  set up parameters and control indices
c
      phideg=pa(1)
      b=pa(2)
      ilfrn=nint(pa(3))
      itfrn=nint(pa(4))
      iopt=nint(pa(5))
      ipset=nint(pa(6))
c
c     compute useful numbers
c
      rho=brho/b
      rhosc=rho/sl
      phi=phideg*pi180
c
c  get multipole values from the parameter set ipset
c
c      do 5 i=1,6
c    5 pb(i)=0.d0
c      goto (10,20,30,40,50),ipset
c      goto 60
c   10 do 11 i=1,6
c   11 pb(i)=pst1(i)
c      goto 60
c   20 do 21 i=1,6
c   21 pb(i)=pst2(i)
c      goto 60
c   30 do 31 i=1,6
c   31 pb(i)=pst3(i)
c      goto 60
c   40 do 41 i=1,6
c   41 pb(i)=pst4(i)
c      goto 60
c   50 do 51 i=1,6
c   51 pb(i)=pst5(i)
c   60 continue
       if (ipset.lt.1 .or. ipset.gt.maxpst) then
         do 50 i=1,6
   50     pb(i) = 0.0d0
       else
         do 60 i=1,6
   60     pb(i) = pst(ipset,i)
       endif
c
c compute multipole coefficients
c procedure when iopt=1
      if (iopt.eq.1) then
c      pb(3)=-pb(3)
c      pb(4)=-pb(4)
      bqd=pb(1)
      aqd=pb(2)
      bsex=pb(3)
      asex=pb(4)
      boct=pb(5)
      aoct=pb(6)
      endif
c procedure when iopt=2
      if (iopt.eq.2) then
      tay1=pb(1)
      aqd=pb(2)
      tay2=pb(3)
      asex=pb(4)
      tay3=pb(5)
      aoct=pb(6)
      bqd=tay1
      bsex=tay2-(5.d0*tay1)/(8.d0*rho)
      boct=tay3-(41.d0*tay1)/(48.d0*rho*rho)+(2.d0*tay2)/(3.d0*rho)
      endif
c procedure when iopt=3
      if (iopt.eq.3) then
      tay1=brho*pb(1)
      aqd=brho*pb(2)
      tay2=brho*pb(3)
      asex=brho*pb(4)
      tay3=brho*pb(5)
      aoct=brho*pb(6)
      bqd=tay1
      bsex=tay2-(5.d0*tay1)/(8.d0*rho)
      boct=tay3-(41.d0*tay1)/(48.d0*rho*rho)+(2.d0*tay2)/(3.d0*rho)
      endif
c
c  scale multipole strengths and compute scaled curvature feed up terms
c
      sbqdf=1.d0/(2.d0*rhosc)
      saqd=aqd*sl*rho/brho
      sbqd=bqd*sl*rho/(2.d0*brho)
c
      sasexf=aqd*(sl**2)/(8.d0*brho)
      sbsexf=bqd*(sl**2)/(8.d0*brho)
      sasex=asex*(sl**2)*rho/(3.d0*brho)
      sbsex=bsex*(sl**2)*rho/(3.d0*brho)
c
      sscalf=-bqd*(sl**3)/(64.d0*rho*brho)
      saoctf=(asex/6.d0-aqd/(16.d0*rho))*(sl**3)/brho
      sboctf=(bsex/12.d0-bqd/(32.d0*rho))*(sl**3)/brho
      saoct=aoct*(sl**3)*rho/brho
      sboct=boct*(sl**3)*rho/(4.d0*brho)
c
c  compute the Hamiltonian
c
      call clear(ha,hm)
c
c  degree 2
c
      ha(7)=sbqdf + sbqd
      ha(8)=saqd
      ha(12)=1.d0/beta
      ha(13)=rhosc/2.d0
      ha(18)=-sbqd
      ha(22)=rhosc/2.d0
      ha(27)=rhosc/(2.d0*(beta*gamma)**2)
c
c  degree 3
c
      ha(28)=sbsexf + sbsex
      ha(30)=sasexf + 3.d0*sasex
      ha(34)=1.d0/2.0d0
      ha(39)=sbsexf - 3.d0*sbsex
      ha(43)=1.d0/2.d0
      ha(44)=sasexf-sasex
      ha(48)=1.d0/(2.0d0*((gamma*beta)**2))
      ha(53)=rhosc/(2.d0*beta)
      ha(76)=rhosc/(2.d0*beta)
      ha(83)=rhosc/(2.d0*(gamma**2)*(beta**3))
c
c  degree 4
c
      ha(84)=sscalf + sboctf + sboct
      ha(86)=saoctf + saoct
      ha(95)=2.d0*sscalf -6.d0*sboct
      ha(109)=1.d0/(2.d0*beta)
      ha(120)=saoctf - saoct
      ha(132)=1.d0/(2.d0*beta)
      ha(139)=1.d0/(2.d0*(beta**3)*(gamma**2))
      ha(140)=rhosc/8.d0
      ha(149)=rhosc/4.d0
      ha(154)=rhosc/(4.d0*(beta**2)*(gamma**2))+rhosc/(2.d0*(beta**2))
      ha(175)=sscalf - sboctf + sboct
      ha(195)=rhosc/8.d0
      ha(200)=rhosc/(4.d0*(beta**2)*(gamma**2))+rhosc/(2.d0*(beta**2))
      ha(209)=rhosc/(2.d0*(beta**4)*(gamma**2))
     &+rhosc/(8.d0*(beta**4)*(gamma**4))
c
c  call GENMAP exponentiation routine
      pc(1)=-phi
      pc(2)=0.d0
      pc(3)=0.d0
      call cex(pc,ha,hm)
      call mapmap(ha,hm,fa,fm)
c
c  add fringe field effects
c  not yet implemented
c
      return
      end
c
***********************************************************************
c
      subroutine srfc(phi0,w,h,mh)
c
c     computes matrix representation mh of linear
c     portion of transfer map and an
c     array h containing coefficients of polynomial
c     generators of nonlinearities
c     for transfer map of a short rf cavity with
c     max potential drop phi0 and
c     frequency w
c     Written by Alex Dragt, ca 1983
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision mh
      dimension h(monoms)
      dimension mh(6,6)
c
      call clear(h,mh)
c
c     set coefficients of mh
c
      do 40 k=1,6
      mh(k,k)=+1.0d0
   40 continue
      mh(6,5)=-(phi0*w*ts)/(brho*c)
c
c     degree 4 only in short buncher limit
c
      h(205)=+(phi0*(ts**3)*(w**3))/(24.d0*c*brho)
      return
      end
c
***********************************************************************
c
      subroutine ssssext(slength,sex,h,hm)
c
c     computes matrix hm and polynomial array
c     h of a sextupole of length slength (meters) and with strength
c     sex (tesla/meter^2) using the SSS algorithm.
c     The vector potential is A_z=-(sex/3)(x**3-3xy**2).
c
c     MARYLIE5.0 upgrade.
c     Written by M.Venturini 5 Aug 1997.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include  'files.inc'
      dimension h(monoms),hm(6,6)
c
      dimension p(6)
c
      call clear(h,hm)
c
      p(1) = slength/sl
      p(2) =0.d0
      p(3) =0.d0
cryne 6/21/2002      p(4) =.1d0
      p(4) =1.d-4
      p(5) =0.d0
      p(6) =0.d0
cryne      p(6) =12
c
      call hamdrift(h)
c
cryne 6/21/2002 modified to multiply by sl**3:
      sex3brho=sex/3.d0/brho*sl**3
c
      h(28)=sex3brho
      h(39)=-3*sex3brho
c
       write(jodf,*)'ssssext has been used'
       write(jodf,*)'eps=',p(4)
c
c        call unixtime(ti1)
        call sss(p,h,hm)
c        call unixtime(ti2)
c       t3=ti2-ti1
c       write(jof,*) 'Execution time in seconds =',t3
c       write(jof,*) '  '
c       write(jodf,*) 'Execution time in seconds =',t3
c       write(jodf,*) '  '
c
      return
      end
c
************************************************************************
c
      subroutine sext(el,gb0,h,hm)
c
c     computes matrix mh and polynomial h
c     for sextupole of length l with strength
c     gb0 (in tesla/meter**2).  routine assumes
c     a vector potential of the form
c     (gb0/3)(x**3-3xy**2)
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension h(monoms)
      dimension hm(6,6)
c
      call ssssext(el,gb0,h,hm)
c
      return
      end
c
***********************************************************************
c
      subroutine thnl (lsqdnr,lsqdsk,lssxnr,lssxsk,lsocnr,lsocsk,h,mh)
c  Does the thin lens map up through octupole.
c  The length-strength product may be given for normal and skew
c  Quads, sextupoles, and octupoles.  The skew elements are
c  rotated clockwise looking in the direction of the beam
c  by half the pole symmetry angle
c  (45, 30, 22.5 degrees respectively).
c  Written by Liam Healy, February 28, 1985.
c
      use beamdata
      use lieaparam, only : monoms
c----Variables----
c  h, mh = output array and matrix
      double precision h(monoms),mh(6,6)
c  ls, f = length*strength (l*gb0), factor (used in calculation)
c  qd, sx, oc : quad, sextupole, octupole
c  nr, sk : normal, skew
      double precision lsqdnr,lsqdsk,lssxnr,lssxsk,lsocnr,lsocsk
      double precision fqdnr,fqdsk,fsxnr,fsxsk,focnr,focsk
c  common quantities
c      double precision brho,c,gamma,gamm1,beta,achg,sl,ts
c
c----Routine----
      fqdnr=lsqdnr*sl/brho
      fqdsk=lsqdsk*sl/brho
      fsxnr=-lssxnr*sl**2/(3.*brho)
      fsxsk=-lssxsk*sl**2/(3.*brho)
      focnr=-lsocnr*sl**3/(4.*brho)
      focsk=-lsocsk*sl**3/(4.*brho)
      call ident(h,mh)
c  Matrix (quadrupole)
      mh(2,1)=-fqdnr
      mh(2,3)=fqdsk
      mh(4,3)=fqdnr
      mh(4,1)=fqdsk
c  Polynomials (sextupole)
      h(28)=fsxnr
      h(30)=-3.*fsxsk
      h(39)=-3.*fsxnr
      h(64)=fsxsk
c  Polynomials (octupole)
      h(84)=focnr
      h(86)=-4.*focsk
      h(95)=-6.*focnr
      h(120)=4.*focsk
      h(175)=focnr
c  That's all folks
      return
      end
c
***********************************************************************
c
      subroutine twsm(iplane,phad,alpha,beta,fa,fm)
c this is a subroutine for computing a linear transfer map
c described in terms of twiss parameters.
c Written b y Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      dimension fa(monoms),fm(6,6)
c-----
c set map to the identity
      call ident(fa,fm)
c compute entries
c      w=phad*pi/(180.d0)
      w=phad*pi180
      cw=cos(w)
      sw=sin(w)
      gam=(1.+alpha*alpha)/beta
c set up subscripts
      k=2*(iplane-1)
      isub1=1+k
      isub2=2+k
c set up matrix
      fm(isub1,isub1)=cw+alpha*sw
      fm(isub1,isub2)=beta*sw
      fm(isub2,isub1)=-gam*sw
      fm(isub2,isub2)=cw-alpha*sw
      return
      end
c
***********************************************************************
c
      subroutine cplm(p,h,mh)
c compressed low order multipole
c Written by Alex Dragt and E. Forest, Fall 1986
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision l,lsc,l1,l2,l3,l4,l5,l6,l7,mh
      dimension h(monoms),p(6)
      dimension mh(6,6)
c
      l=p(1)
c
c set up multipole strengths and normalizations
c eventually the division by l should be removed
c and correspondingly multiplications by l should be removed elsewhere
c in the code
      gs=-p(2)/l
      sgs=-p(3)/l
      go=-p(4)/(4.d0*l)
      sgo=-p(5)/(4.d0*l)
c
      call clear(h,mh)
      lsc=l/sl
c
c     evaluate coefficients of h
c
c
c     set coefficients of mh
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
c
c     set coefficients of h
c
c     degree 3
c
      srs=(sgs*(sl**3))/brho
      sro=(sgo*sl**4)/brho
      rs=(gs*(sl**3))/brho
      ro=(go*sl**4)/brho
      l1=lsc
      l2=lsc**2
      l3=lsc**3
      l4=lsc**4
      l5=lsc**5
      l6=lsc**6
      l7=lsc**7
      h(28)=l1*rs/3.d0
      h(34)=l3*rs/12.d0
      h(39)=-l1*rs
      h(43)=-l3*rs/12.d0
      h(55)=-l3*rs/6.d0
      h(64)=l1*srs/3.d0
      h(68)=l3*srs/12.d0
      h(30)=-l1*srs
      h(50)=-l3*srs/12.d0
      h(36)=-l3*srs/6.d0
c
c     degree 4
c
      sr2=srs**2
      r2=rs**2
      h(84)=l1*ro+l3*r2/12.d0+l3*sr2/12.d0
      h(90)=l3*ro/2.d0
      h(95)=l3*r2/6.d0+l3*sr2/6.d0-6.d0*l1*ro
      h(99)=-l5*r2/30.d0-l5*sr2/30.d0-l3*ro/2.d0
      h(109)=l3*rs/6.d0/beta
      h(111)=l5*r2/15.d0+l5*sr2/15.d0-2.d0*l3*ro
      h(132)=-l3*rs/6.d0/beta
      h(140)=l7*r2/1344.d0+l7*sr2/1344.d0+l5*ro/80.d0
      h(145)=-l5*r2/30.d0-l5*sr2/30.d0-l3*ro/2.d0
      h(149)=l7*r2/672.d0+l7*sr2/672.d0-3.d0*l5*ro/40.d0
      h(161)=-l3*rs/3.d0/beta
      h(175)=l3*r2/12.d0+l3*sr2/12.d0+l1*ro
      h(179)=l3*ro/2.d0
      h(195)=l7*r2/1344.d0+l7*sr2/1344.d0+l5*ro/80.d0
      h(142)=-l5*sro/20.d0
      h(120)=4.d0*l1*sro
      h(156)=l3*sro
      h(86)=-4.d0*l1*sro
      h(124)=l3*sro
      h(106)=-l3*sro
      h(187)=l3*srs/beta/6.d0
      h(148)=-l3*srs/6.d0/beta
      h(92)=-l3*sro
      h(116)=-l3*srs/3.d0/beta
      h(165)=l5*sro/20.d0
      return
      end
c
***********************************************************************
c
      subroutine iftm(p,fa,fm)
c  subroutine for a linear matrix via initial
c  and final twiss parameters
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6)
      write(6,*) 'iftm not yet available'
      return
      end
c
***********************************************************************
c
      subroutine sol(p,fa,fm)
c  subroutine for a solenoid
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6)
      write(6,*) 'sol not yet available'
      return
      end
c
***********************************************************************
c
      subroutine cfqd(pa,ha,hm)
c
c  This subroutine computes the map for a combined function quadrupole.
c  Written by A. Dragt on 9 January 1988.
c  Modified by A. Dragt 8/28/97 to use SSS routine.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      dimension pa(6),pb(6)
      dimension pc(6)
      dimension ha(monoms)
      dimension hm(6,6)
c
      include 'parset.inc'
c----
c  set up parameters and control indices
c
      al=pa(1)
      ipset=nint(pa(2))
      ilfrn=nint(pa(3))
      itfrn=nint(pa(4))
c
c     compute useful numbers
c
      sl2=sl*sl
      sl3=sl*sl2
c
c  get multipole values from the parameter set ipset
c
      if (ipset.lt.1 .or. ipset.gt.maxpst) then
        do 50 i=1,6
   50    pb(i) = 0.0d0
      else
        do 60 i=1,6
   60    pb(i) = pst(i,ipset)
      endif
c
c compute multipole coefficients
c
      bquad=pb(1)
      aquad=pb(2)
      bsex=pb(3)
      asex=pb(4)
      boct=pb(5)
      aoct=pb(6)
c
      fqdnr=bquad/brho
      fqdsk=aquad/brho
      fsxnr=bsex/(3.d0*brho)
      fsxsk=asex/(3.d0*brho)
      focnr=boct/(4.d0*brho)
      focsk=aoct/(4.d0*brho)
c
c  compute the Hamiltonian
c
      call clear(ha,hm)
c
c drift part
c
      call hamdrift(ha)
c
c  quad terms
c
      ha(7)=fqdnr*sl/2.d0
      ha(18)=-ha(7)
      ha(9)=-fqdsk*sl
c
c  sext terms
c
      ha(28)=fsxnr*sl2
      ha(30)=-3.d0*fsxsk*sl2
      ha(39)=-3.d0*fsxnr*sl2
      ha(64)=fsxsk*sl2
c
c  oct terms
c
      ha(84)=focnr*sl3
      ha(86)=-4.0d0*focsk*sl3
      ha(95)=-6.d0*focnr*sl3
      ha(120)=4.d0*focsk*sl3
      ha(175)=focnr*sl3
c
c set up and call sss routine
c
      pc(1)= al/sl
      pc(2)= 0.d0
      pc(3)= 0.d0
      pc(4)= .05d0
      pc(5)= 0.d0
      pc(6)= 0.d0
c
      write(jodf,*) 'ssscfqd has been used'
      write(jodf,*) 'eps=',pc(4)
c
      call sss(pc,ha,hm)
c
c  fringe field effects are put on in the subroutine lmnt
c
      return
      end
c
c*******************************************************************
c
      subroutine thnh (p,h,mh)
c  Does the thin lens map for decapoles and duodecapoles.
c  The length-strength product may be given for normal and skew
c  Decapoles and Duodecapoles.  The skew elements are
c  rotated clockwise looking in the direction of the beam
c  by half the pole symmetry angle ( 18 degs. for decapoles,
c  15 degs. for duodecapoles).
c  Written by F. Neri, July 29, 1988.
c
      use beamdata
      use lieaparam, only : monoms
      implicit none
c----Variables----
c  h, mh = output array and matrix
      double precision h(monoms),mh(6,6),p(6)
c  ls, f = length*strength (l*gb0), factor (used in calculation)
c  qd, sx, oc : quad, sextupole, octupole
c  nr, sk : normal, skew
      double precision lsdenr,lsdesk,lsdunr,lsdusk
      double precision fdenr,fdesk,fdunr,fdusk
c  common quantities
c      double precision brho,c,gamma,gamm1,beta,achg,sl,ts
      external ident
c
c----Routine----
      lsdenr = p(1)
      lsdesk = p(2)
      lsdunr = p(3)
      lsdusk = p(4)
      fdenr=-lsdenr*sl**4/(5.*brho)
      fdesk=-lsdesk*sl**4/(5.*brho)
      fdunr=-lsdunr*sl**5/(6.*brho)
      fdusk=-lsdusk*sl**5/(6.*brho)
      call ident(h,mh)
c  Matrix is identity
c  Polynomials (decapole)
      h(210)=fdenr
      h(212)=-5.*fdesk
      h(221)=-10.*fdenr
      h(246)=10.*fdesk
      h(301)=5.*fdenr
      h(406)=-fdesk
c  Polynomials (duodecapole)
      h(462)=fdunr
      h(464)=-6.*fdusk
      h(473)=-15.*fdunr
      h(498)=20.*fdusk
      h(553)=15.*fdunr
      h(658)=-6.*fdusk
      h(840)=-fdunr
c  That's all folks
      return
      end
c======================================================================
c                       THIRD ORDER ROUTINES
c======================================================================
c
      subroutine drift3(l,h,mh)
c
c     generates linear matrix mh and
c     array h containing nonlinearities
c     for the transfer map describing
c     a drift section of length l meters
c     Written by D. Douglas, ca 1982
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision l,h(monoms),mh(6,6)
c
      double precision lsc
c      dimension j(6)
c
      call clear(h,mh)
      lsc=l/sl
c
c     add drift terms to mh
c
      do 40 k=1,6
      mh(k,k)=+1.0d0
   40 continue
      mh(1,2)=+lsc
      mh(3,4)=+lsc
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     add drift terms to h
c
c     degree 3
c
      h(53)=-(lsc/(2.0d0*beta))
      h(76)=-(lsc/(2.0d0*beta))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      h(140)=-lsc/8.0d0
      h(149)=-lsc/4.0d0
      h(154)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(195)=-lsc/8.0d0
      h(200)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(209)=+lsc*(1.0d0-(5.0d0/beta**2))/(8.0d0*gamma**2*beta**2)
c
      if(monoms.gt.209)h(210:monoms)=0.d0
      return
      end
c
c======================================================================
c
      subroutine dquad3(l,gb0,h,mh)
c
c     computes matrix mh and polynomial array
c     h for horizontally defocussing quadrupole
c     of length l meters and with field
c     gradient of gb0 tesla/meter
c     Written by D. Douglas, ca 1982
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision h(monoms),mh(6,6)
      double precision gb0,k,k2,k3,k4,l,lk,lkm,lsc
c
c     write(6,*)'USING DQUAD3'
      call clear(h,mh)
      lsc=l/sl
c
c     evaluate k,lk
c
      arg=(gb0*(sl**2))/brho
      k=dsqrt(arg)
      lk=lsc*k
      lkm=(-1.0d0)*lk
c
c     set coefficients of mh
c
      chlk=(dexp(lk)+dexp(lkm))/(2.0d0)
      shlk=(dexp(lk)-dexp(lkm))/(2.0d0)
      clk=dcos(lk)
      slk=dsin(lk)
      mh(1,1)=+chlk
      mh(1,2)=+(shlk/k)
      mh(2,2)=+chlk
      mh(2,1)=+(k*shlk)
      mh(3,3)=+clk
      mh(3,4)=+(slk/k)
      mh(4,4)=+clk
      mh(4,3)=-(k*slk)
      mh(5,5)=+1.0d0
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
      mh(6,6)=+1.0d0
c
c     set coefficients of h
c
      tlk=(2.0d0)*lk
      tlkm=(-1.0d0)*tlk
      flk=(4.0d0)*lk
      flkm=(-1.0d0)*flk
      chtlk=(dexp(tlk)+dexp(tlkm))/2.0d0
      shtlk=(dexp(tlk)-dexp(tlkm))/2.0d0
      stlk=dsin(tlk)
      ctlk=dcos(tlk)
      chflk=(dexp(flk)+dexp(flkm))/2.0d0
      shflk=(dexp(flk)-dexp(flkm))/2.0d0
      sflk=dsin(flk)
      cflk=dcos(flk)
c
c      degree 3
c
      h(33)=+((k*(tlk-shtlk))/(8.0d0*beta))
      h(38)=-((1.0d0-chtlk)/(4.0d0*beta))
      h(53)=-((tlk+shtlk)/(8.0d0*k*beta))
      h(67)=-((k*(tlk-stlk))/(8.0d0*beta))
      h(70)=-((1.0d0-ctlk)/(4.0d0*beta))
      h(76)=-((tlk+stlk)/(8.0d0*beta*k))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      k2=k*k
      k3=k*k2
      k4=k2*k2
      c2=clk*clk
      s2=slk*slk
      c4=c2*c2
      s4=s2*s2
      ch2=chlk*chlk
      sh2=shlk*shlk
      ch4=ch2*ch2
      sh4=sh2*sh2
      h(84)=-k3*shflk/2.56d2+k3*shtlk/3.2d1-3.d0*k4*lsc/6.4d1
      h(85)=+k2*sh4/8.d0
      h(90)=-3.d0*k*shflk/1.28d2+3.d0*k2*lsc/3.2d1
      h(95)=+k3*(chtlk-2.d0)*stlk/6.4d1+k4*lsc/1.6d1
     &+k3*(ctlk-2.d0)*shtlk/6.4d1
      h(96)=-k2*shtlk*stlk/3.2d1+k2*(chtlk-2.d0)*ctlk/3.2d1
     &+k2/3.2d1
      h(99)=-k*(chtlk-2.d0)*stlk/6.4d1
     &-k*(ctlk+2.d0)*shtlk/6.4d1
     &+k2*lsc/1.6d1
      h(104)=3.d0*k*shtlk/(3.2d1*beta**2)
     &-k2*lsc*(chtlk+2.d0)/(1.6d1*beta**2)
     &+k*(1.d0-3.d0/beta**2)*(shtlk-tlk)/1.6d1
      h(105)=+(ch4-1.d0)/8.d0
      h(110)=-k2*shtlk*stlk/3.2d1-k2*chtlk*(ctlk-2.d0)/3.2d1
     &-k2/3.2d1
      h(111)=+k*chtlk*stlk/1.6d1-k*shtlk*ctlk/1.6d1
      h(114)=+shtlk*stlk/3.2d1-3.d0/3.2d1
     &+chtlk*(ctlk+2.d0)/3.2d1
      h(119)=-chflk/(6.4d1*beta**2)
     &+lk*shtlk/(8.d0*beta**2)
     &-chtlk/(1.6d1*beta**2)+ch4/(8.d0*beta**2)-ch2/(4.d0*beta**2)
     &+1.3d1/(6.4d1*beta**2)
     &+(1.d0-3.d0/beta**2)*(1.d0-ch2)/4.d0
      h(140)=-shflk/(2.56d2*k)-shtlk/(3.2d1*k)-3.d0*lsc/6.4d1
      h(145)=+k*(chtlk+2.d0)*stlk/6.4d1
     &+k*(ctlk-2.d0)*shtlk/6.4d1-k2*lsc/1.6d1
      h(146)=-shtlk*stlk/3.2d1-3.d0/3.2d1
     &+(chtlk+2.d0)*ctlk/3.2d1
      h(149)=-(chtlk+2.d0)*stlk/(6.4d1*k)
     &-(ctlk+2.d0)*shtlk/(6.4d1*k)-lsc/1.6d1
      h(154)=+shtlk/(k*3.2d1*beta**2)-lsc*chtlk/(1.6d1*beta**2)
     &+(1.d0-3.d0/beta**2)*(shtlk+tlk)/(1.6d1*k)
      h(175)=-k3*sflk/2.56d2+k3*stlk/3.2d1-3.d0*k4*lsc/6.4d1
      h(176)=-k2*s4/8.d0
      h(179)=+3.d0*k*sflk/1.28d2-3.d0*k2*lsc/3.2d1
      h(184)=-k*(1.d0-3.d0/beta**2)*(stlk-tlk)/1.6d1
     &-3.d0*k*stlk/(3.2d1*beta**2)
     &+k2*lsc*(ctlk+2.d0)/(1.6d1*beta**2)
      h(185)=+(c4-1.d0)/8.d0
      h(190)=+(1.d0-3.d0/beta**2)*(1.d0-c2)/4.d0
     &-cflk/(6.4d1*beta**2)-lk*stlk/(8.d0*beta**2)
     &-ctlk/(1.6d1*beta**2)+c4/(8.d0*beta**2)
     &-c2/(4.d0*beta**2)+1.3d1/(6.4d1*beta**2)
      h(195)=-sflk/(2.56d2*k)-stlk/(3.2d1*k)-3.d0*lsc/6.4d1
      h(200)=+(1.d0-3.d0/beta**2)*(tlk+stlk)/(1.6d1*k)
     &+stlk/(k*3.2d1*beta**2)-lsc*ctlk/(1.6d1*beta**2)
      h(209)=+lsc*(1.d0-5.d0/beta**2)/(8.d0*gamma**2*beta**2)
c
      if(monoms.gt.209)h(210:monoms)=0.d0
      return
      end
c
c======================================================================
c
      subroutine fquad3(l,gb0,h,mh)
c
c     computes matrix mh and polynomial array
c     h for horizontally focussing quadrupole
c     of length l meters and with field
c     gradient of gb0 tesla/meter
c     Written by D. Douglas, ca 1982
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision h(monoms),mh(6,6)
c
      double precision gb0,k,k2,k3,k4,l,lk,lkm,lsc
c
c     write(6,*)'USING FQUAD3'
      call clear(h,mh)
      lsc=l/sl
c
c     evaluate k,lk
c
      arg=(gb0*(sl**2))/brho
      k=dsqrt(arg)
      lk=lsc*k
      lkm=(-1.0d0)*lk
c
c     set coefficients of mh
c
      chlk=(dexp(lk)+dexp(lkm))/(2.0d0)
      shlk=(dexp(lk)-dexp(lkm))/(2.0d0)
      clk=dcos(lk)
      slk=dsin(lk)
      mh(3,3)=+chlk
      mh(3,4)=+(shlk/k)
      mh(4,4)=+chlk
      mh(4,3)=+(k*shlk)
      mh(1,1)=+clk
      mh(1,2)=+(slk/k)
      mh(2,2)=+clk
      mh(2,1)=-(k*slk)
      mh(5,5)=+1.0d0
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
      mh(6,6)=+1.0d0
c
c     set coefficients of h
c
      tlk=(2.0d0)*lk
      tlkm=(-1.0d0)*tlk
      flk=(4.0d0)*lk
      flkm=(-1.0d0)*flk
      chtlk=(dexp(tlk)+dexp(tlkm))/2.0d0
      shtlk=(dexp(tlk)-dexp(tlkm))/2.0d0
      stlk=dsin(tlk)
      ctlk=dcos(tlk)
      chflk=(dexp(flk)+dexp(flkm))/2.0d0
      shflk=(dexp(flk)-dexp(flkm))/2.0d0
      sflk=dsin(flk)
      cflk=dcos(flk)
c
c      degree 3
c
      h(67)=+((k*(tlk-shtlk))/(8.0d0*beta))
      h(70)=-((1.0d0-chtlk)/(4.0d0*beta))
      h(76)=-((tlk+shtlk)/(8.0d0*k*beta))
      h(33)=-((k*(tlk-stlk))/(8.0d0*beta))
      h(38)=-((1.0d0-ctlk)/(4.0d0*beta))
      h(53)=-((tlk+stlk)/(8.0d0*beta*k))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      k2=k*k
      k3=k*k2
      k4=k2*k2
      c2=clk*clk
      s2=slk*slk
      c4=c2*c2
      s4=s2*s2
      ch2=chlk*chlk
      sh2=shlk*shlk
      ch4=ch2*ch2
      sh4=sh2*sh2
      h(175)=-k3*shflk/2.56d2+k3*shtlk/3.2d1-3.d0*k4*lsc/6.4d1
      h(176)=+k2*sh4/8.d0
      h(179)=-3.d0*k*shflk/1.28d2+3.d0*k2*lsc/3.2d1
      h(95)=+k3*(chtlk-2.d0)*stlk/6.4d1+k4*lsc/1.6d1
     &+k3*(ctlk-2.d0)*shtlk/6.4d1
      h(110)=-k2*shtlk*stlk/3.2d1+k2*(chtlk-2.d0)*ctlk/3.2d1
     &+k2/3.2d1
      h(145)=-k*(chtlk-2.d0)*stlk/6.4d1
     &-k*(ctlk+2.d0)*shtlk/6.4d1
     &+k2*lsc/1.6d1
      h(184)=+3.d0*k*shtlk/(3.2d1*beta**2)
     &-k2*lsc*(chtlk+2.d0)/(1.6d1*beta**2)
     &+k*(1.d0-3.d0/beta**2)*(shtlk-tlk)/1.6d1
      h(185)=+(ch4-1.d0)/8.d0
      h(96)=-k2*shtlk*stlk/3.2d1-k2*chtlk*(ctlk-2.d0)/3.2d1
     &-k2/3.2d1
      h(111)=+k*chtlk*stlk/1.6d1-k*shtlk*ctlk/1.6d1
      h(146)=+shtlk*stlk/3.2d1-3.d0/3.2d1
     &+chtlk*(ctlk+2.d0)/3.2d1
      h(190)=-chflk/(6.4d1*beta**2)
     &+lk*shtlk/(8.d0*beta**2)
     &-chtlk/(1.6d1*beta**2)+ch4/(8.d0*beta**2)-ch2/(4.d0*beta**2)
     &+1.3d1/(6.4d1*beta**2)
     &+(1.d0-3.d0/beta**2)*(1.d0-ch2)/4.d0
      h(195)=-shflk/(2.56d2*k)-shtlk/(3.2d1*k)-3.d0*lsc/6.4d1
      h(99)=+k*(chtlk+2.d0)*stlk/6.4d1
     &+k*(ctlk-2.d0)*shtlk/6.4d1-k2*lsc/1.6d1
      h(114)=-shtlk*stlk/3.2d1-3.d0/3.2d1
     &+(chtlk+2.d0)*ctlk/3.2d1
      h(149)=-(chtlk+2.d0)*stlk/(6.4d1*k)
     &-(ctlk+2.d0)*shtlk/(6.4d1*k)-lsc/1.6d1
      h(200)=+shtlk/(k*3.2d1*beta**2)-lsc*chtlk/(1.6d1*beta**2)
     &+(1.d0-3.d0/beta**2)*(shtlk+tlk)/(1.6d1*k)
      h(84)=-k3*sflk/2.56d2+k3*stlk/3.2d1-3.d0*k4*lsc/6.4d1
      h(85)=-k2*s4/8.d0
      h(90)=+3.d0*k*sflk/1.28d2-3.d0*k2*lsc/3.2d1
      h(104)=-k*(1.d0-3.d0/beta**2)*(stlk-tlk)/1.6d1
     &-3.d0*k*stlk/(3.2d1*beta**2)
     &+k2*lsc*(ctlk+2.d0)/(1.6d1*beta**2)
      h(105)=+(c4-1.d0)/8.d0
      h(119)=+(1.d0-3.d0/beta**2)*(1.d0-c2)/4.d0
     &-cflk/(6.4d1*beta**2)-lk*stlk/(8.d0*beta**2)
     &-ctlk/(1.6d1*beta**2)+c4/(8.d0*beta**2)
     &-c2/(4.d0*beta**2)+1.3d1/(6.4d1*beta**2)
      h(140)=-sflk/(2.56d2*k)-stlk/(3.2d1*k)-3.d0*lsc/6.4d1
      h(154)=+(1.d0-3.d0/beta**2)*(tlk+stlk)/(1.6d1*k)
     &+stlk/(k*3.2d1*beta**2)-lsc*ctlk/(1.6d1*beta**2)
      h(209)=+lsc*(1.d0-5.d0/beta**2)/(8.d0*gamma**2*beta**2)
c
      if(monoms.gt.209)h(210:monoms)=0.d0
      return
      end
c
***********************************************************************
c
      subroutine sext3(l,gb0,h,mh)
c
c     computes matrix mh and polynomial h
c     for sextupole of length l with strength
c     gb0 (in tesla/meter**2).  routine assumes
c     a vector potential of the form
c     (gb0/3)(x**3-3xy**2)
c     Written by D. Douglas, ca 1982
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision l,lsc,lsc2,lsc3,lsc4,lsc5,lsc6,lsc7,mh
      dimension h(monoms)
      dimension mh(6,6)
c
      call clear(h,mh)
      lsc=l/sl
c
c     evaluate coefficients of h
c
c
c     set coefficients of mh
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
      mh(1,2)=+lsc
      mh(3,4)=+lsc
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     set coefficients of h
c
c     degree 3
c
      rk=-(gb0*(sl**3))/(brho*3.0d0)
      h(28)=+(lsc*rk)
      h(29)=-((3.0d0*(lsc**2)*rk)/2.0d0)
      h(34)=+((lsc**3)*rk)
      h(39)=-(3.0d0*lsc*rk)
      h(40)=+(3.0d0*(lsc**2)*rk)
      h(43)=-((lsc**3)*rk)
      h(49)=-(((lsc**4)*rk)/4.0d0)
      h(53)=-(lsc/(2.0d0*beta))
      h(54)=+((3.0d0*(lsc**2)*rk)/2.0d0)
      h(55)=-(2.0d0*(lsc**3)*rk)
      h(58)=+((3.0d0*(lsc**4)*rk)/4.0d0)
      h(76)=-(lsc/(2.0d0*beta))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      lsc2=lsc**2
      lsc3=lsc**3
      lsc4=lsc**4
      lsc5=lsc**5
      lsc6=lsc**6
      lsc7=lsc**7
      rk2=rk**2
      h(84)=+7.5d-1*rk2*lsc3
      h(85)=-1.5d0*rk2*lsc4
      h(90)=+1.125d0*rk2*lsc5
      h(95)=+1.5d0*rk2*lsc3
      h(96)=-1.5d0*rk2*lsc4
      h(99)=+3.d0*rk2*lsc5/4.d1
      h(105)=-3.75d-1*rk2*lsc6
      h(109)=+rk*lsc3/(2.0d0*beta)
      h(110)=-1.5d0*rk2*lsc4
      h(111)=+2.1d0*rk2*lsc5
      h(114)=-3.75d-1*rk2*lsc6
      h(132)=-rk*lsc3/(2.0d0*beta)
      h(140)=+3.0d0*rk2*lsc7/5.6d1-lsc/8.0d0
      h(144)=-rk*lsc4/(4.0d0*beta)
      h(145)=+3.0d0*rk2*lsc5/4.d1
      h(146)=-3.75d-1*rk2*lsc6
      h(149)=+3.0d0*rk2*lsc7/2.8d1-lsc/4.0d0
      h(154)=+lsc*(1.0d0-(3.0d0/(beta**2)))/4.0d0
      h(161)=-rk*lsc3/beta
      h(167)=+7.5d-1*rk*lsc4/beta
      h(175)=+7.5d-1*rk2*lsc3
      h(176)=-1.5d0*rk2*lsc4
      h(179)=+1.125d0*rk2*lsc5
      h(185)=-3.75d-1*rk2*lsc6
      h(195)=+3.0d0*rk2*lsc7/5.6d1-lsc/8.0d0
      h(200)=+lsc*(1.d0-(3.d0/(beta**2)))/4.d0
      h(209)=+(lsc*(1.d0-(5.0d0/beta**2)))/(8.0d0*gamma**2*beta**2)
c
      if(monoms.gt.209)h(210:monoms)=0.d0
      return
      end
c
***********************************************************************
c
      subroutine octm3(l,gb0,h,mh)
c
c     computes matrix mh and polynomial h for
c     magnetic octupole with length l meters and vector
c     potential of the form (gb0/4.)*(+x**4-6*x**2*y**2+y**4)
c     Written by D. Douglas, ca 1982
c
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision l,h(monoms),mh(6,6)
      double precision lsc
c
      call clear(h,mh)
      lsc=l/sl
c
c     enter matrix elements
c
      do 40 i=1,6
      mh(i,i)=+1.0d0
   40 continue
      mh(1,2)=+lsc
      mh(3,4)=+lsc
      mh(5,6)=+(lsc/((gamma**2)*(beta**2)))
c
c     set coefficients of h
c
c     degree 3
c
      h(53)=-(lsc/(2.0d0*beta))
      h(76)=-(lsc/(2.0d0*beta))
      h(83)=-(lsc/(2.0d0*(gamma**2)*(beta**3)))
c
c     degree 4
c
      ap=-(sl**4)*gb0/(4.d0*brho)
      h(84)=+lsc*ap
      h(85)=-2.0d0*(lsc**2)*ap
      h(90)=2.0d0*(lsc**3)*ap
      h(95)=-6.0d0*lsc*ap
      h(96)=6.d0*lsc**2*ap
      h(99)=-(2.0d0*(lsc**3)*ap)
      h(105)=-(lsc**4)*ap
      h(110)=+6.0d0*(lsc**2)*ap
      h(111)=-(8.0d0*(lsc**3)*ap)
      h(114)=+3.0d0*(lsc**4)*ap
      h(140)=-lsc/8.0d0+(lsc**5)*ap/5.0d0
      h(145)=-(2.0d0*(lsc**3)*ap)
      h(146)=+3.0d0*(lsc**4)*ap
      h(149)=-lsc/4.0d0-(6.0d0*(lsc**5)*ap)/5.0d0
      h(154)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(175)=+lsc*ap
      h(176)=-2.0d0*(lsc**2)*ap
      h(179)=+2.0d0*(lsc**3)*ap
      h(185)=-(lsc**4)*ap
      h(195)=-lsc/8.0d0+(lsc**5)*ap/5.0d0
      h(200)=+(lsc*(1.0d0-(3.0d0/(beta**2))))/4.0d0
      h(209)=+(lsc*(1.d0-(5.d0/beta**2)))/(8.d0*gamma**2*beta**2)
c
      if(monoms.gt.209)h(210:monoms)=0.d0
      return
      end
c
***********************************************************************
