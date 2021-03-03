c
      subroutine cgbend(rho,bndang,psideg,phideg,ijopt,coefpset,h,mh)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision mh(6,6),mht1(6,6),mht2(6,6),mht3(6,6)
      double precision h(monoms),ht1(monoms),ht2(monoms),ht3(monoms)
      double precision ptmp(6),coefpset(6)
c
c this routine is modeled on the Healy/Dragt subroutine gbend
c
c  calculate each of the 3 individual pieces that make
c  up the general bending magnet, and concatenate them:
c  cgbend=hpf1*cfbend*hpf2
c****NOV 7 "call hpf" replaced w/ calls to gbend+cfbend since routine hpf
c    does not know about combined function magnets.  RDR
c    (note: "gbend" is the subroutine name for gbdy)
c    (note: "cfbend" is the subroutine name for normal entry comb fxn magnet)
c
c 5 parameters and a 6-vector are passed to cfbend:
c     bend angle, b field, ilfrn, itfrn, ijopt, coefpset(1-6)
c but note that, in the current implementation of cfbend, ilfrn and itfrn
c are not used.
      aldeg=bndang-psideg-phideg
c     write(6,*)'brho=',brho
c     write(6,*)'rho=',rho
      ptmp(2)=brho/rho
      ptmp(3)=0.
      ptmp(4)=0.
      ptmp(5)=ijopt
c   ( ptmp(6) is not used )
      ptmp(6)=0.
c  compute hpf1 and nbend
c          call hpf(rho,1,psideg,ht1,mht1)
      azero=0.d0
cryne 11/9/02      ptmp(1)=0.5d0*bndang
      ptmp(1)=psideg
      call gbend(rho,azero,psideg,azero,ht1,mht1)
      call cfbend(ptmp,coefpset,ht2,mht2)
      call concat(ht1,mht1,ht2,mht2,ht1,mht1)
      ptmp(1)=aldeg
      call cfbend(ptmp,coefpset,ht2,mht2)
c  form the product hpf1*nbend
      call concat(ht1,mht1,ht2,mht2,ht3,mht3)
c  compute hpf2
c          call hpf(rho,-1,phideg,ht1,mht1)
c     write(6,'(6(1pe12.5,1x))')((mht1(i,j),j=1,6),i=1,6)
c     if(1.gt.0)call myexit
cryne 11/9/02      ptmp(1)=0.5d0*bndang
      ptmp(1)=phideg
      call cfbend(ptmp,coefpset,ht1,mht1)
      call gbend(rho,azero,azero,phideg,ht2,mht2)
      call concat(ht1,mht1,ht2,mht2,ht1,mht1)
c  form the complete product hpf1*nbend*hpf2
      call concat(ht3,mht3,ht1,mht1,h,mh)
      return
      end
c
c------------------------------------------------------------------
c
      subroutine cgfrngg(iw,eangle,rho,kfrn,gap,fint,iopt,coefpset,h,mh)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      include 'pie.inc'
cryne 7/13/2002      include 'frnt.inc'
      double precision mh(6,6),mht1(6,6),mht2(6,6),mht3(6,6)
      double precision h(monoms),ht1(monoms),ht2(monoms),ht3(monoms)
      integer lfrn,tfrn,kfrn
      dimension coefpset(6)
c     write(6,*)'inside cgfrngg w/ coefpset='
c     write(6,*)coefpset
c     write(6,*)'and kfrn=',kfrn
c kfrn = 1,2,or 3  for dipole and/or quad fringe fields
c iw denotes which fringe:  =1 for leading,  = 2 for trailing
c
c this is not the most efficient routine, but it is easy to follow.
c
      call ident(h,mh)
c
      if(iw.eq.2)goto 200
      lfrn=kfrn
c leading fringe (iw=1):
c put on leading dipole and quad fringe fields
      if(lfrn.eq.0)return
c compute and put on leading dipole fringe field
cryne 7/13/2002      gap=cfbgap
cryne 7/13/2002      xk1=cfblk1
      xk1=fint
      if((lfrn.eq.1).or.(lfrn.eq.3))then
        call gfrngg(eangle,rho,1,ht1,mht1,gap,xk1)
        call concat(ht1,mht1,h,mh,h,mh)
      endif
c compute and put on leading quad fringe fields if lfrn=2 or 3
      if(lfrn.eq.1)return
      if(iopt.eq.1) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.2) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.3) then
      bqd=brho*coefpset(1)
      aqd=brho*coefpset(2)
      endif
c compute and put on normal quad fringe field
c     write(6,*)'about to call frquad (1st time) w/ bqd=',bqd
      call frquad(bqd,-1,ht1,mht1)
c     write(6,*)'returned from 1st call to frqad w/ mht1='
c     do i=1,6
c     write(6,'(6(1pe12.5,1x))')(mht1(i,j),j=1,6)
c     enddo
      call concat(ht1,mht1,h,mh,h,mh)
c compute and put on skew quad fringe field
cryne 12/21/2004      write(6,*)'cgfrngg: check this IF test'
      if(aqd.ne.0.d0)then
      angr=-pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c     write(6,*)'about to call frquad (2nd time) w/ aqd=',aqd
      call frquad(aqd,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      endif
      return
c
  200 continue
      tfrn=kfrn
c trailing fringe (iw=2):
c put on trailing dipole and quad fringe fields
      if(tfrn.eq.0)return
c compute and put on trailing dipole fringe field
cryne 7/13/2002      gap=cfbgap
cryne 7/13/2002      xk1=cfbtk1
      xk1=fint
      if((tfrn.eq.1).or.(tfrn.eq.3))then
        call gfrngg(eangle,rho,2,ht1,mht1,gap,xk1)
        call concat(h,mh,ht1,mht1,h,mh)
      endif
c compute and put on trailing quad fringe fields if tfrn=2 or 3
      if(tfrn.eq.1)return
      if(iopt.eq.1) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.2) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.3) then
      bqd=brho*coefpset(1)
      aqd=brho*coefpset(2)
      endif
c compute and put on normal quad fringe field
c     write(6,*)'about to call frquad (3rd time) w/ bqd=',bqd
      call frquad(bqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c compute and put on skew quad fringe field
cryne 12/21/2004      write(6,*)'cgfrngg: also check this IF test'
      if(aqd.ne.0.d0)then
      angr=pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c     write(6,*)'about to call frquad (4th time) w/ aqd=',aqd
      call frquad(aqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      endif
      return
      end
