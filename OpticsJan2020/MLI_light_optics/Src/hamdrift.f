c
***************************************************
c
      subroutine hamdrift(h)
c
c     computes the Hamiltonian for a drift.
c     MARYLIE5.0 upgrade.
c     Written by M.Venturini 5 Aug 1997.
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension h(monoms)
c
      do i=7,monoms
       h(i)=0.d0
      enddo
c
      ptg=-1/beta
      ptg2=ptg*ptg
      ptg3=ptg*ptg2
      ptg4=ptg*ptg3
      ptg5=ptg*ptg4
      ptg6=ptg*ptg5
c
      h(13)=1/2.d0
      h(22)=1/2.d0
      h(27)=(-1 + ptg2)/2.d0
      h(53)=-ptg/2.d0
      h(76)=-ptg/2.d0
      h(83)=(3*ptg - 3*ptg3)/6.d0
      h(140)=1/8.d0
      h(149)=1/4.d0
      h(154)=(-1 + 3*ptg2)/4.d0
      h(195)=1/8.d0
      h(200)=(-1 + 3*ptg2)/4.d0
      h(209)=(3 - 18*ptg2 + 15*ptg4)/24.d0
      h(340)=-3*ptg/8.d0
      h(363)=-3*ptg/4.d0
      h(370)=(9*ptg - 15*ptg3)/12.d0
      h(443)=-3*ptg/8.d0
      h(450)=(9*ptg - 15*ptg3)/12.d0
      h(461)=(-45*ptg + 150*ptg3 - 105*ptg5)/120.d0
      h(714)=1/16.d0
      h(723)=3/16.d0
      h(728)=(-9 + 45*ptg2)/48.d0
      h(769)=3/16.d0
      h(774)=(-3 + 15*ptg2)/8.d0
      h(783)=(9 - 90*ptg2 + 105*ptg4)/48.d0
      h(896)=1/16.d0
      h(901)=(-9 + 45*ptg2)/48.d0
      h(910)=(9 - 90*ptg2 + 105*ptg4)/48.d0
      h(923)=(-45 + 675*ptg2 - 1575*ptg4 + 945*ptg6)/720.d0
c
      return
      end
c end of file

