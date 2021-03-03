      subroutine transit(el,en,h,hm)
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision hm(6,6),h(monoms)
      write(6,*)'inside transit; el,en=',el,en
c
      call clear(h,hm)
c
c matrix:
      do k=1,6
      hm(k,k)=+1.0d0
      enddo
      hm(1,2)=el/en
      hm(3,4)=el/en
c f4:
      h(140)=-0.125d0*el/en**3
      h(149)=-2.d0*0.125d0*el/en**3
      h(195)=-0.125d0*el/en**3
c f6:
      h(714)=-0.0625d0*el/en**5
      h(723)=-3.d0*0.0625d0*el/en**5
      h(769)=-3.d0*0.0625d0*el/en**5
      h(896)=-0.0625d0*el/en**5
      return
      end
c
c
c
      subroutine interface(enm,enp,b2,b4,b6,h,hm)
      use lieaparam, only : monoms
ccccc include 'impli.inc'
      implicit none
      integer k
      double precision enm,enp,b2,b4,b6
      double precision endiff,term1,term2,term3,term4,term5
      double precision term6,term7,term8,term9
      double precision hm(6,6),h(monoms)
      write(6,*)'inside interface; enm,enp=',enm,enp
      write(6,*)'b2,b4,b6=',b2,b4,b6
c
      endiff=enm-enp
      call clear(h,hm)
c
c matrix:
      do k=1,6
      hm(k,k)=+1.0d0
      enddo
      hm(2,1)=2.d0*b2*endiff
      hm(4,3)=2.d0*b2*endiff
c f4:
      term1=-endiff/enm*(enm*(2.d0*b2**3-b4) - 2.d0*b2**3*enp)
      term2=2.d0*b2**2*endiff/enm
      term3=0.5d0*b2*endiff/(enm*enp)
c (q^2)^2
      h(84)=term1
      h(95)=2.d0*term1
      h(175)=term1
c (q^2)p.q
      h(85)=term2
      h(96)=term2
      h(110)=term2
      h(176)=term2
c (q^2)(p^2)
      h(90)=term3
      h(99)=term3
      h(145)=term3
      h(179)=term3
c f6:
      term4=(1.d0/enm**3)*( (b6-6.d0*b2**2*b4+2.d0*b2**5)*enm**4           &
     &     + (-b6+12.d0*b2**2*b4-4.d0*b2**5)*enm**3*enp                    &
     &     -6.d0*b2**2*b4*enm**2*enp**2+4.d0*b2**5*enm*enp**3              &
     &     -2.d0*b2**5*enp**4 )
      term5=(1.d0/(enm**3*enp))*(                                          &
     &(2.d0*b2*b4-4.d0*b2**4)*enm**4+(2.d0*b2*b4+6.d0*b2**4)*enm**3*enp    &
     &-(4.d0*b2*b4+4.d0*b2**4)*enm**2*enp**2                               &
     &+6.d0*b2**4*enm*enp**3-4.d0*b2**4*enp**4 )
      term6=0.5d0/(enm**3*enp)*(                                           &
     &b4*enm**3-b4*enm**2*enp+2.d0*b2**3*enm*enp**2-2.d0*b2**3*enp**3)
      term9=2.d0*b2**3/(enm**3*enp)*                                       &
     &(enm**3-enm**2*enp+enm*enp**2-enp**3)
      term7=0.5d0*b2**2/(enm**3*enp**2)*                                   &
     &(enm**3+enm*enp**2-2.d0*enp**3)
      term8=0.125d0*b2/(enm**3*enp**3)*(enm**3-enp**3)
c f6:
c (q^2)^3
      h(462)=term4
      h(473)=3.d0*term4
      h(553)=3.d0*term4
      h(840)=term4
c (q^2)^2 q.p
      h(463)=term5
      h(474)=term5
      h(488)=2.d0*term5
      h(554)=2.d0*term5
      h(623)=term5
      h(841)=term5
c (q^2)^2 p^2
      h(468)=term6
      h(477)=term6
      h(523)=2.d0*term6
      h(557)=2.d0*term6
      h(749)=term6
      h(844)=term6
c q^2 p^2 q.p
      h(483)=term7
      h(492)=term7
      h(524)=term7
      h(563)=term7
      h(593)=term7
      h(750)=term7
      h(627)=term7
      h(850)=term7
c q^2 (p^2)^2
      h(518)=term8
      h(527)=2.d0*term8
      h(573)=term8
      h(719)=term8
      h(753)=2.d0*term8
      h(860)=term8
c q^2 (p.q)^2
      h(468)=h(468)+term9
      h(489)=2.d0*term9
      h(523)=h(523)+term9
      h(557)=h(557)+term9
      h(624)=2.d0*term9
      h(844)=h(844)+term9
      return
      end
c
c
c
      subroutine rootmap(en,b2,b4,b6,h,hm)
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision hm(6,6),h(monoms)
c
      call clear(h,hm)
c
c matrix:
      do k=1,6
      hm(k,k)=+1.0d0
      enddo
      hm(2,1)=2.d0*b2*en
      hm(4,3)=2.d0*b2*en
c f4:
      term1=en*(b4-2.d0*b2**3)
      term2=2.d0*b2**2
      term3=-0.5d0*b2/en
c (q^2)^2
      h(84)=term1
      h(95)=2.d0*term1
      h(175)=term1
c (q^2)p.q
      h(85)=term2
      h(96)=term2
      h(110)=term2
      h(176)=term2
c (q^2)(p^2)
      h(90)=term3
      h(99)=term3
      h(145)=term3
      h(179)=term3
c f6:
      term4=en*(b6-6.d0*b2**2*b4+2.d0*b2**5)
      term5=4.d0*b2*b4-2.d0*b2**4
      term6=-0.5d0*b4/en
      term7=0.5d0*b2**2/en**2
      term8=-0.125d0*b2/en**3
c (q^2)^3
      h(462)=term4
      h(473)=3.d0*term4
      h(553)=3.d0*term4
      h(840)=term4
c (q^2)^2 q.p
      h(463)=term5
      h(474)=term5
      h(488)=2.d0*term5
      h(554)=2.d0*term5
      h(623)=term5
      h(841)=term5
c (q^2)^2 p^2
      h(468)=term6
      h(477)=term6
      h(523)=2.d0*term6
      h(557)=2.d0*term6
      h(749)=term6
      h(844)=term6
c q^2 p^2 q.p
      h(483)=term7
      h(492)=term7
      h(524)=term7
      h(563)=term7
      h(593)=term7
      h(750)=term7
      h(627)=term7
      h(850)=term7
c q^2 (p^2)^2
      h(518)=term8
      h(527)=2.d0*term8
      h(573)=term8
      h(719)=term8
      h(753)=2.d0*term8
      h(860)=term8
      return
      end
c
      subroutine optirot(angdeg,ijkind,h,mh)
c
c     High order PROT routine.
c     Actual code generated by Johannes van Zeijts using
c     REDUCE.
c
cryne use beamdata
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
cryne B = beta
      B = 1.0d0
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
