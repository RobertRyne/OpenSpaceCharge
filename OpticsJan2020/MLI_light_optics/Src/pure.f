***********************************************************************
* header              PURIFYING (PURE) ROUTINES                       *
*  Routines for computing conjugacy classes                           *
***********************************************************************
*
      subroutine da2(fm,a2a,a2m)
c this is a subroutine that generates the matrix a2m
c that brings fm to block form in the dynamic case
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fm(6,6),a2a(monoms),a2m(6,6)
      dimension reval(6),aieval(6)
      dimension revec(6,6),aievec(6,6)
c clear the array a2a
      do 10 i=1,monoms
   10 a2a(i)=0.
c compute the eigenvectors of fm
      call eig6(fm,reval,aieval,revec,aievec)
c sort these eigenvectors
      call dvsort(revec,aievec)
c compute a2m from these eigenvectors
      do 20 i=1,6
      a2m(i,1)=revec(1,i)
      a2m(i,2)=aievec(1,i)
      a2m(i,3)=revec(3,i)
      a2m(i,4)=aievec(3,i)
      a2m(i,5)=revec(5,i)
      a2m(i,6)=aievec(5,i)
   20 continue
c rephase the result
      call drphse(a2m)
      return
      end
c
***********************************************************************
c
      subroutine dpur2(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the f2 (matrix)
c part of a dynamic map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map.
c ta,tm is the transforming map. that is g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
c find the map a2:
      call da2(fm,ta,tm)
c go to floquet variables:
      call sndwch(ta,tm,fa,fm,ga,gm)
      return
      end
c
********************************************************************
c
      subroutine dpur3(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the f3
c part of a dynamic map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map, and ta,tm is the purifying transformation.
c that is, g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      ibrief = 1
      detmin = 1.d-12
      min=28
      max=83
      call gdpur(fa,fm,ga,gm,ta,tm,min,max,ibrief,detmin)
      return
      end
c
******************************************************
c
      subroutine dpur4(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the f4
c part of a dynamic map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map, and ta,tm is the purifying transformation.
c that is, g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      ibrief = 1
      detmin = 1.d-12
      min=90
      max=208
      call gdpur(fa,fm,ga,gm,ta,tm,min,max,ibrief,detmin)
      return
      end
c
************************************************************
c
      subroutine fxpt(fa,fm,ana,anm,ta,tm)
c This is a subroutine for finding the fixed point of a map.
c The initial map is given by fa and fm. It is unchanged by the
c subroutine. The map about the closed orbit corresponding to
c the fixed point is given by ana and anm. The transformation
c to the fixed point is given by the map ta,tm.
c Written by Alex Dragt, 11 October 1985.
      use parallel, only : idproc
      implicit double precision (a-h,o-z)
      dimension fa(923),fm(6,6)
      dimension ana(923),anm(6,6)
      dimension ta(923),tm(6,6)
      dimension tta(923),ttm(6,6)
      dimension tempa(923),tempm(6,6)
      dimension am(4,4),av(4),alphv(4),augm(4,5)
c Computation of an1:
c Set up am = 4x4 block of fm-i:
      do 10 i=1,4
      do 10 j=1,4
      am(i,j)=fm(i,j)
      if(i.eq.j) am(i,j)=am(i,j)-1.d0
   10 continue
c Set up av:
      do 20 i=1,4
   20 av(i)=fm(i,6)
c Compute alphv:
      call leshs(alphv,4,am,av,augm,det)
c Write out message about determinant:
cryne Jan 3, 2005 added ".lt. 1.d-3" to only print det if it is small
      if(idproc.eq.0 .and. det.lt.1.d-3)write(6,600) det
  600 format(1x,'det in fxpt is',e15.4)
c Compute t1:
      call ident(ta,tm)
      do 30 i=1,4
   30 tm(i,6)=-alphv(i)
      tm(5,1)=alphv(2)
      tm(5,2)=-alphv(1)
      tm(5,3)=alphv(4)
      tm(5,4)=-alphv(3)
c Store t1
      call mapmap(ta,tm,tta,ttm)
c Compute t1*m*t1inv where m is initial map:
      call sndwch(ta,tm,fa,fm,ana,anm)
c      call mycat(6,ta,tm,fa,fm,tempa,tempm)
c      call inv(ta,tm)
c      call mycat(6,tempa,tempm,ta,tm,ana,anm)
c Computation of an2:
c Set up am:
      call mapmap(ana,anm,tempa,tempm)
      call inv(tempa,tempm)
      do 40 i=1,4
      do 40 j=1,4
      am(i,j)=tempm(j,i)
      if(i.eq.j) am(i,j)=am(i,j)-1.d0
   40 continue
c Set up av:
      av(1)=-ana(48)
      av(2)=-ana(63)
      av(3)=-ana(73)
      av(4)=-ana(79)
c Compute alphv:
      call leshs(alphv,4,am,av,augm,det)
c Compute t2:
      call ident(ta,tm)
      ta(48)=alphv(1)
      ta(63)=alphv(2)
      ta(73)=alphv(3)
      ta(79)=alphv(4)
c Compute and store t2*t1
c      call mycat(6,ta,tm,tta,ttm,tempa,tempm)
cryne 8/16/02      call concat(6,ta,tm,tta,ttm,tempa,tempm)
      call concat(ta,tm,tta,ttm,tempa,tempm)
      call mapmap(tempa,tempm,tta,ttm)
c Compute t2*an1*t2inv:
      call sndwch(ta,tm,ana,anm,ana,anm)
c      call mycat(6,ta,tm,ana,anm,tempa,tempm)
c      call inv(ta,tm)
c      call mycat(6,tempa,tempm,ta,tm,ana,anm)
c Computation of an3:
c The matrix am is unchanged from the previous step.
c Set up av:
      av(1)=-ana(139)
      av(2)=-ana(174)
      av(3)=-ana(194)
      av(4)=-ana(204)
c Compute alphv:
      call leshs(alphv,4,am,av,augm,det)
c Compute t3:
      call ident(ta,tm)
      ta(139)=alphv(1)
      ta(174)=alphv(2)
      ta(194)=alphv(3)
      ta(204)=alphv(4)
c Compute and store t3*t2*t1:
c      call mycat(6,ta,tm,tta,ttm,tempa,tempm)
cryne 8/16/02      call concat(6,ta,tm,tta,ttm,tempa,tempm)
      call concat(ta,tm,tta,ttm,tempa,tempm)
      call mapmap(tempa,tempm,tta,ttm)
c Compute t3*an2*t3inv:
      call sndwch(ta,tm,ana,anm,ana,anm)
c      call mycat(6,ta,tm,ana,anm,tempa,tempm)
c      call inv(ta,tm)
c      call mycat(6,tempa,tempm,ta,tm,ana,anm)
c Computation of an4:
c The matrix am is the same for all orders.
c Set up av:
      av(1)=-ana(335)
      av(2)=-ana(405)
      av(3)=-ana(440)
      av(4)=-ana(455)
c Compute alphv:
      call leshs(alphv,4,am,av,augm,det)
c Compute t4:
      call ident(ta,tm)
      ta(335)=alphv(1)
      ta(405)=alphv(2)
      ta(440)=alphv(3)
      ta(455)=alphv(4)
c Compute and store t4*t3*t2*t1:
c      call mycat(6,ta,tm,tta,ttm,tempa,tempm)
cryne 8/16/02      call concat(6,ta,tm,tta,ttm,tempa,tempm)
      call concat(ta,tm,tta,ttm,tempa,tempm)
      call mapmap(tempa,tempm,tta,ttm)
c Compute t4*an3*t4inv:
      call sndwch(ta,tm,ana,anm,ana,anm)
c Computation of an5:
c Set up av:
      av(1)=-ana(713)
      av(2)=-ana(839)
      av(3)=-ana(895)
      av(4)=-ana(916)
c Compute alphv:
      call leshs(alphv,4,am,av,augm,det)
c Compute t5:
      call ident(ta,tm)
      ta(713)=alphv(1)
      ta(839)=alphv(2)
      ta(895)=alphv(3)
      ta(916)=alphv(4)
c Compute and store t5*t4*t3*t2*t1:
c      call mycat(6,ta,tm,tta,ttm,tempa,tempm)
cryne 8/16/02      call concat(6,ta,tm,tta,ttm,tempa,tempm)
      call concat(ta,tm,tta,ttm,tempa,tempm)
      call mapmap(tempa,tempm,tta,ttm)
c Compute an5 = t5*an4*t5inv:
      call sndwch(ta,tm,ana,anm,ana,anm)
c Store t5*t4*t3*t2*t1 in t:
      call mapmap(tta,ttm,ta,tm)
      return
      end
c
***********************************************************************
c
      subroutine gdpur(fa,fm,ga,gm,ta,tm,min,max,
     &                 ibrief,detmin)
c  This is a generic purifying routine for the
c  dynamic resonance case
c  fa,fm is the original map, which is left unchanged.
c  ga,gm is the purified map, and ta,tm is the purifying
c  transformation, that is g=t*f*(t inverse)
c  min and max are the minimum and maximum index
c  of the monomials to be purified
c Written by Alex Dragt and F. Neri, Spring 1986
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(*),ga(*),ta(*),t1a(monoms),t2a(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6),t1m(6,6),t2m(6,6)
      dimension ax(-4:4),bx(-4:4),ay(-4:4),by(-4:4),at(-4:4),bt(-4:4)
      include 'dr.inc'
c
c  Set up multiple angle arrays
      cwx = fm(1,1)
      swx = fm(1,2)
      cwy = fm(3,3)
      swy = fm(3,4)
      cwt = fm(5,5)
      swt = fm(5,6)
      ax(0)=1.d0
      bx(0)=0.d0
      ay(0)=1.d0
      by(0)=0.d0
      at(0)=1.d0
      bt(0)=0.d0
      ax(1)=cwx
      bx(1)=swx
      ay(1)=cwy
      by(1)=swy
      at(1)=cwt
      bt(1)=swt
c
      do 10 n=2,4
        ax(n)=ax(1)*ax(n-1)-bx(1)*bx(n-1)
        bx(n)=bx(1)*ax(n-1)+ax(1)*bx(n-1)
   10 continue
c
      do 20 n=2,4
        ay(n)=ay(1)*ay(n-1)-by(1)*by(n-1)
        by(n)=by(1)*ay(n-1)+ay(1)*by(n-1)
   20 continue
c
      do 30 n=2,4
        at(n)=at(1)*at(n-1)-bt(1)*bt(n-1)
        bt(n)=bt(1)*at(n-1)+at(1)*bt(n-1)
   30 continue
c
      do 40 n=1,4
        ax(-n)=ax(n)
        bx(-n)=-bx(n)
        ay(-n)=ay(n)
        by(-n)=-by(n)
        at(-n)=at(n)
        bt(-n)=-bt(n)
   40 continue
c
c  Resonance decompose map:
      call ctodr(fa,t1a)
c
c  Set up map to remove offensive terms of index min->max:
      call ident(t2a,t2m)
      do 100 k=min,max
        if (drexp(0,k).ne.0) then
c  compute cth=cos(theta) and sth=sin(theta) for
c  theta = drexp(1,k)*wx + drexp(2,k)*wy + drexp(3,k)*wt
          nx=drexp(1,k)
          ny=drexp(2,k)
          nt=drexp(3,k)
          axnx=ax(nx)
          bxnx=bx(nx)
          ayny=ay(ny)
          byny=by(ny)
          atnt=at(nt)
          btnt=bt(nt)
          cth=(axnx*ayny-bxnx*byny)*atnt-(axnx*byny+bxnx*ayny)*btnt
          sth=(axnx*ayny-bxnx*byny)*btnt+(axnx*byny+bxnx*ayny)*atnt
c
c  Carry out rest of calculation
          det = 2.d0 * (1.d0 - cth)
          k1 = k
          k2 = k+1
          if ( ibrief.ne.1 ) then
            if(idproc.eq.0)
     &      write(6,*) ' Det in subspace',k1,',',k2,' is',det
          endif
          if ( dabs(det) .gt. detmin) then
            t2a(k1) = ((1.-cth)*t1a(k1) + sth*t1a(k2))/det
            t2a(k2) = ((1.-cth)*t1a(k2) - sth*t1a(k1))/det
          else
            if(idproc.eq.0)
     &      write(6,*) ' Det(',k1,',',k2,') =',det,' not removed'
          endif
        endif
  100 continue
c transfom map t2 to cartesian basis; the result is the map t:
      call ident(ta,tm)
      call drtoc(t2a,ta)
c remove offensive terms ( min->max )
      call sndwch(ta,tm,fa,fm,ga,gm)
      return
      end
c
************************************************************************
c
      subroutine gspur(fa,fm,ga,gm,ta,tm,min,max,
     &                 ibrief,detmin)
c  This is a generic purifying routine for the
c  static resonance case
c  fa,fm is the original map, which is left unchanged.
c  ga,gm is the purified map, and ta,tm is the purifying
c  transformation, that is g=t*f*(t inverse)
c  min and max are the minimum and maximum index
c  of the monomials to be purified
c Written by Alex Dragt and F. Neri, Spring 1986
      use parallel, only : idproc
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(*),ga(*),ta(*),t1a(monoms),t2a(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6),t1m(6,6),t2m(6,6)
      dimension ax(-4:4),bx(-4:4),ay(-4:4),by(-4:4)
      include 'sr.inc'
c
c  Set up multiple angle arrays
      cwx = fm(1,1)
      swx = fm(1,2)
      cwy = fm(3,3)
      swy = fm(3,4)
      ax(0)=1.d0
      bx(0)=0.d0
      ay(0)=1.d0
      by(0)=0.d0
      ax(1)=cwx
      bx(1)=swx
      ay(1)=cwy
      by(1)=swy
c
      do 10 n=2,4
        ax(n)=ax(1)*ax(n-1)-bx(1)*bx(n-1)
        bx(n)=bx(1)*ax(n-1)+ax(1)*bx(n-1)
   10 continue
c
      do 20 n=2,4
        ay(n)=ay(1)*ay(n-1)-by(1)*by(n-1)
        by(n)=by(1)*ay(n-1)+ay(1)*by(n-1)
   20 continue
c
      do 30 n=1,4
        ax(-n)=ax(n)
        bx(-n)=-bx(n)
        ay(-n)=ay(n)
        by(-n)=-by(n)
   30 continue
c
c  Resonance decompose map:
      call ctosr(fa,t1a)
c
c  Set up map to remove offensive terms of index min->max:
      call ident(t2a,t2m)
      do 100 k=min,max
        if (srexp(0,k).ne.0) then
c  compute cth =cos(theta) and sth=sin(theta) for
c  theta = srexp(1,k)*wx + srexp(2,k)*wy
          nx=srexp(1,k)
          ny=srexp(2,k)
          axnx=ax(nx)
          bxnx=bx(nx)
          ayny=ay(ny)
          byny=by(ny)
          cth=axnx*ayny-bxnx*byny
          sth=bxnx*ayny+axnx*byny
c
c  Carry out rest of calculation
          det = 2.d0 * (1.d0 - cth)
          k1 = k
          k2 = k+1
          if ( ibrief.ne.1 ) then
            if(idproc.eq.0)
     &      write(6,*) ' Det in subspace',k1,',',k2,' is',det
          endif
          if ( dabs(det) .gt. detmin) then
            t2a(k1) = ((1.-cth)*t1a(k1) + sth*t1a(k2))/det
            t2a(k2) = ((1.-cth)*t1a(k2) - sth*t1a(k1))/det
          else
            if(idproc.eq.0)
     &      write(6,*) ' Det(',k1,',',k2,') =',det,' not removed'
          endif
        endif
  100 continue
c  transform map to cartesian basis; the result is the map t:
      call ident(ta,tm)
      call srtoc(t2a,ta)
c  remove offensive terms ( min->max )
      call sndwch(ta,tm,fa,fm,ga,gm)
      return
      end
c
*************************************************************
c
      subroutine sa2(fm,a2a,a2m)
c this is a subroutine that generates the matrix a2m
c that brings fm to block form in the static case
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fm(6,6),a2a(monoms),a2m(6,6)
      dimension reval(6),aieval(6)
      dimension revec(6,6),aievec(6,6)
c clear the array a2a
      do 10 i=1,monoms
   10 a2a(i)=0.
c compute the eigenvectors of fm
      call eig4(fm,reval,aieval,revec,aievec)
c sort these eigenvectors
       call svsort(revec,aievec)
c compute a2m from these eigenvectors
      do 20 i=1,6
      a2m(i,1)=revec(1,i)
      a2m(i,2)=aievec(1,i)
      a2m(i,3)=revec(3,i)
      a2m(i,4)=aievec(3,i)
      a2m(i,5)=revec(5,i)
      a2m(i,6)=revec(6,i)
   20 continue
c rephase the result
      call srphse(a2m)
      return
      end
c
***********************************************************************
c
      subroutine scpur3(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the chromatic f3
c part of a static map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map, and ta,tm is the purifying transformation.
c that is, g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      ibrief = 1
      detmin = 1.d-12
      min=31
      max=42
      call gspur(fa,fm,ga,gm,ta,tm,min,max,ibrief,detmin)
      return
      end
c
**************************************************
c
      subroutine sgpur3(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the geometric f3
c part of a static map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map, and ta,tm is the purifying transformation.
c that is, g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      ibrief = 1
      detmin = 1.d-12
      min=43
      max=62
      call gspur(fa,fm,ga,gm,ta,tm,min,max,ibrief,detmin)
      return
      end
c
******************************************************
c
      subroutine spur2(fa,fm,ga,gm,ta,tm,t2m)
c this is a subroutine for purifying the f2 (matrix)
c part of a static map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map.
c ta,tm is the transforming map. that is g=t*f*(t inverse).
c the matrix t2m associated with a2 is also returned.
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms),t1a(monoms),t2a(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6),t1m(6,6),t2m(6,6)
c go to the off momentum closed orbit:
      call fxpt(fa,fm,ta,tm,t1a,t1m)
c find the map a2:
      call sa2(tm,t2a,t2m)
c go to floquet variables:
      call sndwch(t2a,t2m,ta,tm,ga,gm)
c accumulate transforming map:
      call concat(t2a,t2m,t1a,t1m,ta,tm)
      return
      end
c
********************************************************************
c
c
      subroutine spur4(fa,fm,ga,gm,ta,tm)
c this is a subroutine for purifying the f4
c part of a static map.
c fa,fm is the original map, and it is left unchanged.
c ga,gm is the purified map, and ta,tm is the purifying transformation.
c that is, g=t*f*(t inverse).
c Written by Alex Dragt, Spring 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      ibrief = 1
      detmin = 1.d-12
      min=90
      max=153
      call gspur(fa,fm,ga,gm,ta,tm,min,max,ibrief,detmin)
      return
      end
c
c  end of file
