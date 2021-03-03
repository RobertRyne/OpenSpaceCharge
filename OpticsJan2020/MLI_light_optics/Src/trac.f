************************************************************************
* header:                 TRACK                                        *
*  Particle tracking, both symplectic and non-symplectic.              *
************************************************************************
c
      subroutine eval(amh,norder,ntrace,nwrite,nunit,idmin,idmax,         &
     &                nprecision)
c
c     evaluates the image zf of initial
c     data zi under the lie transformation
c     with linear part having matrix representation
c     mh and nonlinear part having generators
c     with coefficients stored in the array h
c
      use multitrack   !includes 'use rays'  May 9, 2006
      use lieaparam, only : monoms
      include 'impli.inc'
cryne May 22, 2006 need reftraj from map.inc
      include 'map.inc'
      include 'ind.inc'
      include 'pbkh.inc'
      include 'lims.inc'
      include 'vblist.inc'
      include 'files.inc'
c
c xmh is a local array, amh is in the argument list:
      dimension xmh(6,6), amh(6,6)
      dimension tempz(6),vect(923)
c
c     write(6,*)'HELLO FROM TRACE'
      write(6,*)'inside eval with norder=',norder
      if(.not.allocated(tblock))then
        write(6,*)'error:tblock not allocated'
        stop
      endif
c
c     generate tempz=amh*zi (linear contribution)
c
      ktrace=ntrace
      if(ktrace.eq.0)ktrace=1
c
      do 9999 idum=1,ktrace
      do 5678 nnnn=1,nraysp
c
c store 6-vector in zi.
cNote that zi(5), zi(6) get changed below if multitrac.ne.0
      zi(1:6)=zblock(1:6,nnnn)
c
c Default (multitrac=0) is to use the current transfer map mh:
      if(multitrac.eq.0)xmh(1:6,1:6)=amh(1:6,1:6)
c
      if(multitrac.ne.0)then
        imin=inear(nnnn)
        jmin=jnear(nnnn)
        xmh(1:6,1:6)=tmhlist(1:6,1:6,imin,jmin)
        zi(5)=tdelt(nnnn)
        zi(6)=ptdelt(nnnn)
      endif
c
c
c perform matrix-vector multiply:
      tempz(1:6)=0.d0
      do 30 i=1,6
      do 20 k=1,6
        tempz(i)=tempz(i)+xmh(i,k)*zi(k)
   20 continue
   30 continue
      tblock(1:6,nnnn)=tempz(1:6)

      if(norder.eq.1)goto 5001
c
c     loop to compute final values zf(i), including nonlinearities
        vect(1:923)=0.d0
        vect(1:6) = tblock(1:6,nnnn)
        do i=7,27
          vect(i) = vect(index1(i)) * vect(index2(i))
        enddo
c       if(norder.ge.3)then
          do i=28,83
            vect(i) = vect(index1(i)) * vect(index2(i))
          enddo
c       endif
c       if(norder.ge.4)then
          do i=84,209
            vect(i) = vect(index1(i)) * vect(index2(i))
          enddo
c       endif
c       if(norder.ge.5)then
          do i=210,461
            vect(i) = vect(index1(i)) * vect(index2(i))
          enddo
c       endif
c       if(norder.ge.6)then
          do i=462,923
            vect(i) = vect(index1(i)) * vect(index2(i))
          enddo
c       endif

        klast=923
        if(norder.eq.2)klast=27
        if(norder.eq.3)klast=83
        if(norder.eq.4)klast=209
        if(norder.eq.5)klast=461
        if(norder.eq.6)klast=923
c       if(multitrac.eq.0)then
c       <nothing to do because brkts(h) was called from routine trace>
c       endif
        if(multitrac.ne.0)then
          pbh(:,:)=pbhlist(:,:,imin,jmin)
c         call brkts(hlist(1,imin,jmin))
        endif
        do i=1,6
          tmp = 0.d0
          do k=7,klast
            tmp = tmp + pbh(k,i)*vect(k)
          enddo
          tblock(i,nnnn)=tblock(i,nnnn) + tmp
        enddo
 5001 continue
c
c if using multiple maps, then shift the result back to with respect
c to the main reference trajectory:
      if(multitrac.ne.0)then
        tblock(5,nnnn)=tblock(5,nnnn)+ tlistfin(imin,jmin)-reftraj(5)
        tblock(6,nnnn)=tblock(6,nnnn)+ptlistfin(imin,jmin)-reftraj(6)
      endif

 5678 continue

cryne 3/25/04 : if ntrace.ge.0, copy the tblock data to zblock, 
cryne           then print what is in zblock.
cryne      Otherwise, print what is tblock. NOTE WELL: this will
cryne      ruin the tblock array, but that is OK since it is no
cryne      longer needed after that data have been printed.
c
c Case with ntrace.ge.1:
      if(ntrace.ge.1)then
!       zblock(:,1:nraysp)=tblock(:,1:nraysp)
        do j=1,nraysp
          if(istat(j).eq.0)zblock(1:6,j)=tblock(1:6,j)
        enddo
          if(nwrite.ne.0)then
            if( (mod(idum,nwrite).eq.0) .and. (jfcf.gt.0) )then
            call pwritez(nunit,idmin,idmax,nprecision,0,0)
            endif
          endif
          if(nwrite.ne.0)then
            if( (mod(idum,nwrite).eq.0) .and. (jfcf.lt.0) )then
            jfcf=-jfcf
            call pwritezfp
            jfcf=-jfcf
            endif
          endif
      endif
c Case with ntrace.eq.0 (print rays trace results if requested,
c but do not do tracking, i.e. leave zblock unchanged):
      if(ntrace.eq.0)then
        if(nwrite.ne.0)then
          if( (mod(idum,nwrite).eq.0) .and. (jfcf.gt.0) )then
          call pwritet(nunit,idmin,idmax,nprecision)
          endif
        endif
        if(nwrite.ne.0)then
          if( (mod(idum,nwrite).eq.0) .and. (jfcf.lt.0) )then
          jfcf=-jfcf
          call pwritetfp
          jfcf=-jfcf
          endif
        endif
      endif
 9999 continue
c     write(6,*)'LEAVING EVAL'
      return
      end
c
***********************************************************************
c
      subroutine evalsr_old(mh,zi,zf,df,rdf,rrjac)
cryne 8/13/02      sbroutne evalsr_old(mh,zi,zf)
c  subroutine to evaluate the image zf of initial data zi under the
c  lie transformation with linear part mh and nonlinear part
c  represented by a standard representation stored in common/deriv/df
c  this subroutine modified on 8/13/86 to change the
c  order in which exp(:f2:) acts.
c  canx was also modified accordingly.  this
c  whole subroutine needs rewriting badly AJD
c
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
cryne 8/13/02      include 'deriv.inc'
      include 'files.inc'
      include 'ind.inc'
      include 'len.inc'
      double precision mh
      dimension mh(6,6)
!!!!! dimension df(6,monoms),rdf(3,monom1+1),rrjac(3,3,monom2+1)
      dimension df(6,monoms),rdf(3,monoms),rrjac(3,3,monoms)
      dimension ztmp1(6),ztmp2(6),ztmp3(6)
      dimension vect(monoms)
c zlm added for ordering modification 8/13/86
      dimension zlm(6)
      dimension zi(6),zf(6)
c
c  initialize arrays
         zf(1:6)=0.d0
         ztmp1(1:6)=0.d0
         ztmp2(1:6)=0.d0
         ztmp3(1:6)=0.d0
c also initialize zlm
         zlm(1:6)=0.d0
c  the following section is added due to reordering:
c
c  compute effect of linear part of the map on zi
c
      do 40 i=1,6
         do 4 j=1,6
            zlm(i)=zlm(i) + mh(i,j)*zi(j)
    4    continue
   40 continue
c
c  call newton search routine to find value of new momentum
c
c  the following statement is commented out due to reordering
c     call newt(zi,ztmp1,rdf,rrjac)
c  the following statement is added due to reordering
      call newt(zlm,ztmp1,rdf,rrjac)
cryne 8/15/02      call newt(zlm,ztmp1)
c
c  compute vector containing values of basis monomials
c
      do 2 i=1,6
          vect(i) = ztmp1(i)
    2 continue
c
c  (only terms through order imaxi-1 are required)
c
      do 20 i = 7,len(imaxi-1)
      vect(i) = vect(index1(i))*vect(index2(i))
 20   continue
c  compute values of the new coords and old momenta using the standard
c  rep of the transfer map which is stored in df
c
      do 3000 i=1,5,2
         ivalue=i+1
         do 300 j=1,len(imaxi-1)
c
c  new coordinate values
            ztmp2(i)=ztmp2(i)+df(ivalue,j)*vect(j)
c      ztmp2(i)= ddot(len(imaxi-1),vect,1,df(ivalue,1),6)
c  old momentum values (done as a check)
            ztmp3(ivalue)=ztmp3(ivalue)+df(i,j)*vect(j)
  300 continue
c      ztmp3(ivalue)=ddot(len(imaxi-1),vect,1,df(i,1),6)
c
c  transfer new momentum values to ztmp2 (these were returned from
c  newt in ztmp1)
c
         ztmp2(ivalue)=ztmp2(ivalue)+ztmp1(ivalue)
 3000 continue
c
c  at this point, the nonlinearities are fixed.  the image of
c  the initial data under the nonlinear portion of the
c  transformation is stored in ztmp2.  before applying the linear
c  part of the map to this, check to see if the old momentum values
c  which were read in match those computed using the standard
c  representation (done immediately above.)
c
c  nonlinearities check
c
      delpx=zi(2)-ztmp3(2)
      delpy=zi(4)-ztmp3(4)
      delpz=zi(6)-ztmp3(6)
      if(ibrief.ne.2)goto 302
      write(jof,9250)
 9250 format(1h ,' **momentum deviations** ')
      write(jof,9300)delpx,delpy,delpz
 9300 format(1h ,3(d22.15,2x))
  302 continue
c  the following statements are commented out due to reordering
c
c  compute effect of linear part of the map on ztmp2
c
c     do 40 i=1,6
c        do 4 j=1,6
c           zf(i)=zf(i) + mh(i,j)*ztmp2(j)
c   4    continue
c  40 continue
c  add the following lines due to reordering:
      do 137 i=1,6
  137 zf(i)=ztmp2(i)
      return
      end
c
***********************************************************************
c
      subroutine evalsr(mh,df,rdf,rrjac,norder,ntrace,nwrite,nunit,          &
     &                  idmin,idmax,nprecision)
c  subroutine to evaluate the image zf of initial data zi under the
c  lie transformation with linear part mh and nonlinear part
c  represented by a standard representation stored in common/deriv/df
c  this subroutine modified on 8/13/86 to change the
c  order in which exp(:f2:) acts.
c  canx was also modified accordingly.  this
c  whole subroutine needs rewriting badly AJD
c
      use rays
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
cryne 8/13/02      include 'deriv.inc'
      include 'files.inc'
      include 'ind.inc'
      include 'len.inc'
      double precision mh
      dimension mh(6,6)
!!!!! dimension df(6,monoms),rdf(3,monom1+1),rrjac(3,3,monom2+1)
      dimension df(6,monoms),rdf(3,monoms),rrjac(3,3,monoms)
      dimension ztmp1(6),ztmp2(6),ztmp3(6)
      dimension vect(monoms)
c zlm added for ordering modification 8/13/86
      dimension zlm(6)
!!!!! dimension zi(6),zf(6)
c
      ktrace=ntrace
      if(ktrace.eq.0)ktrace=1
c
      do 9999 idum=1,ktrace
c     write(6,*)'idum=',idum
c
      lda=6
      ldb=6
      ldc=6
      call DGEMM('N','N',6,nraysp,6,1.0d0,mh,lda,
     &zblock,ldb,0.0d0,tblock,ldc)
c
      if(norder.eq.1)then
        goto 5001
      else
        write(6,*) '(evalsr) computing nonlinear part'
      endif
c
      do 8500 nnn=1,nraysp
         zi(1:6)=zblock(1:6,nnn)
         zf(1:6)=0.d0
         ztmp1(1:6)=0.d0
         ztmp2(1:6)=0.d0
         ztmp3(1:6)=0.d0
c
c  the following section is added due to reordering:
c
c  compute effect of linear part of the map on zi
c
!     do 40 i=1,6
!        do 4 j=1,6
!           zlm(i)=zlm(i) + mh(i,j)*zi(j)
!   4    continue
!  40 continue
      zlm(1:6)=tblock(1:6,nnn)
c
c  call newton search routine to find value of new momentum
c
c  the following statement is commented out due to reordering
c     call newt(zi,ztmp1,rdf,rrjac)
c  the following statement is added due to reordering
      call newt(zlm,ztmp1,rdf,rrjac)
cryne 8/15/02      call newt(zlm,ztmp1)
c
c  compute vector containing values of basis monomials
c
      do 2 i=1,6
          vect(i) = ztmp1(i)
    2 continue
c
c  (only terms through order imaxi-1 are required)
c
      do 20 i = 7,len(imaxi-1)
      vect(i) = vect(index1(i))*vect(index2(i))
 20   continue
c  compute values of the new coords and old momenta using the standard
c  rep of the transfer map which is stored in df
c
      do 3000 i=1,5,2
         ivalue=i+1
         do 300 j=1,len(imaxi-1)
c
c  new coordinate values
            ztmp2(i)=ztmp2(i)+df(ivalue,j)*vect(j)
c      ztmp2(i)= ddot(len(imaxi-1),vect,1,df(ivalue,1),6)
c  old momentum values (done as a check)
            ztmp3(ivalue)=ztmp3(ivalue)+df(i,j)*vect(j)
  300 continue
c      ztmp3(ivalue)=ddot(len(imaxi-1),vect,1,df(i,1),6)
c
c  transfer new momentum values to ztmp2 (these were returned from
c  newt in ztmp1)
c
         ztmp2(ivalue)=ztmp2(ivalue)+ztmp1(ivalue)
 3000 continue
c
c  at this point, the nonlinearities are fixed.  the image of
c  the initial data under the nonlinear portion of the
c  transformation is stored in ztmp2.  before applying the linear
c  part of the map to this, check to see if the old momentum values
c  which were read in match those computed using the standard
c  representation (done immediately above.)
c
c  nonlinearities check
c
      delpx=zi(2)-ztmp3(2)
      delpy=zi(4)-ztmp3(4)
      delpz=zi(6)-ztmp3(6)
      if(ibrief.ne.2)goto 302
      write(jof,9250)
 9250 format(1h ,' **momentum deviations** ')
      write(jof,9300)delpx,delpy,delpz
 9300 format(1h ,3(d22.15,2x))
  302 continue
c  the following statements are commented out due to reordering
c
c  compute effect of linear part of the map on ztmp2
c
c     do 40 i=1,6
c        do 4 j=1,6
c           zf(i)=zf(i) + mh(i,j)*ztmp2(j)
c   4    continue
c  40 continue
c  add the following lines due to reordering:
      do 137 i=1,6
  137 zf(i)=ztmp2(i)
c
cryne May 22, 2006      if(ntrace.ge.1)then
cryne May 22, 2006        zblock(1:6,nnn)=zf(1:6)
cryne May 22, 2006      endif
c As in routine eval, the result should temporarily be in tblock:
        tblock(1:6,nnn)=zf(1:6)
cryne May 22, 2006
 8500 continue
c
 5001 continue
c
cryne May 22, 2006 -- Same code as from routine eval:
c
c Case with ntrace.ge.1:
      if(ntrace.ge.1)then
!       zblock(:,1:nraysp)=tblock(:,1:nraysp)
        do j=1,nraysp
          if(istat(j).eq.0)zblock(1:6,j)=tblock(1:6,j)
        enddo
          if(nwrite.ne.0)then
            if( (mod(idum,nwrite).eq.0) .and. (jfcf.gt.0) )then
            call pwritez(nunit,idmin,idmax,nprecision,0,0)
            endif
          endif
          if(nwrite.ne.0)then
            if( (mod(idum,nwrite).eq.0) .and. (jfcf.lt.0) )then
            jfcf=-jfcf
            call pwritezfp
            jfcf=-jfcf
            endif
          endif
      endif
c Case with ntrace.eq.0 (print rays trace results if requested,
c but do not do tracking, i.e. leave zblock unchanged):
      if(ntrace.eq.0)then
        if(nwrite.ne.0)then
          if( (mod(idum,nwrite).eq.0) .and. (jfcf.gt.0) )then
          call pwritet(nunit,idmin,idmax,nprecision)
          endif
        endif
        if(nwrite.ne.0)then
          if( (mod(idum,nwrite).eq.0) .and. (jfcf.lt.0) )then
          jfcf=-jfcf
          call pwritetfp
          jfcf=-jfcf
          endif
        endif
      endif
 9999 continue
c     write(6,*)'LEAVING EVALSR'
      return
      end
c
***********************************************************************
c
c  subroutine canx
c  the following line is commented out due to reordering;
c  see the note below.
c     subroutine canx(mh,h)
c  the following line is inserted due to reordering.
      subroutine canx(mh,h,norder)
c  subroutine to establish standard rep of transfer map for
c  use in evaluating ray traces.
c
cryne 5/22/02      implicit none
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'deriv.inc'
      include 'len.inc'
      include 'ind.inc'
      external pbkt1,pbkt2,pmadd,product
      integer i,j,k,l,m,n,nn,ior
      integer ivalue,jvalue
      integer index1,index2,jv,len,imaxi
      integer iq,ip,jq,jp,kq,kp
cryne 5/22/02      double precision df,rjac
      double precision deriv3
      double precision mh(6,6)
      double precision dg(923,6)
      double precision h(923),g(923),gtemp(923)
      double precision dgtmp(923),crstrm(923),f(923)
      double precision temp(923),dum(923),tmp1(923),tmp2(923)
      double precision t1(6),t2(27),t3(83),t4(209),t5(461),t6(923)
cryne      common/deriv/df(6,923)
cryne      common/rjacob/rjac(3,3,923)
cryne 7/23/2002      common /len/len(16)
cryne 7/23/2002      common/ind/imaxi,jv(923),index1(923),index2(923)
c
c  initialize arrays
c
      do 1000 n=1,len(imaxi)
        g(n)=0.d0
        gtemp(n)=0.d0
        dgtmp(n)=0.d0
        crstrm(n)=0.d0
        f(n)=0.d0
        t6(n) = 0.0d0
        temp(n)=0.d0
        tmp1(n)=0.d0
        tmp2(n)=0.d0
        dum(n)=0.d0
        do 1 m=1,6
          dg(n,m)=0.d0
          df(m,n)=0.d0
    1   continue
        do 100 l=1,3
          do 10 m=1,3
            rjac(l,m,n)=0.d0
   10     continue
  100   continue
 1000 continue
c
c  transfor nonlinear arrays by linear piece
c
c      write(6,*) '  The Transfer Map is: '
c      do 7777 k=1,len(imaxi)
c        if(dabs(h(k)).gt.1.d-10) write(6,*) k,h(k)
c 7777 continue
c      do 200 k=3,imaxi
c         do 2 l=1,len(imaxi)
c            gtemp(l)=0.d0
c    2    continue
c         call xform5(h,k,mh,gtemp)
c         do 20 l=1,len(imaxi)
c            g(l)=g(l)+gtemp(l)
c   20    continue
c  200 continue
c
c   The nonlinear part of h is not transformed, just copied to g:
       do 20 l=1,len(imaxi)
         g(l) = h(l)
  20   continue
c
c  determine derivs of array g
c
      do 303 ior = 3,imaxi-1
      do 300 i=1,6
         do 3 j=1,len(imaxi)
            dgtmp(j)=0.d0
    3    continue
         call pbkt2(g,ior,i,dgtmp)
         do 30 j=1,len(imaxi)
            dg(j,i) = dg(j,i) + dgtmp(j)
   30    continue
  300 continue
  303 continue
c
       ivalue = 0
       do 7 i= 1,5,2
           ivalue = ivalue+2
          call product(dg(1,i),2,dg(1,ivalue),2,crstrm)
          if(imaxi.gt.4) then
            call product(dg(1,i),2,dg(1,ivalue),3,crstrm)
            call product(dg(1,i),3,dg(1,ivalue),2,crstrm)
          endif
          if (imaxi.gt.5) then
c  The following implements the 2 derivative part of the contribution to F6.
            call product(dg(1,i),3,dg(1,ivalue),3,crstrm)
            call product(dg(1,i),2,dg(1,ivalue),4,crstrm)
            call product(dg(1,i),4,dg(1,ivalue),2,crstrm)
          endif
   7   continue
c
c  sum all the contributions to F excluding terms
c  with more than 2 derivs and commutators.
c
      do 4 n=1,len(imaxi)
         f(n)=f(n)-g(n)- crstrm(n)/2.0d0
    4 continue
c
c
c   Terms with 4 derivatives
      if(imaxi.gt.4) then
          do 808 nn = 1,461
  808       t5(nn) = 0.0d0
        ip = 0
        do 800 iq = 1,5,2
          ip = ip+2
          call pbkt2(crstrm,4,iq,t3)
          call product(t3,3,dg(1,ip),2,t5)
  800   continue
        ip = 0
        do 900 iq = 1,5,2
          ip = ip+2
          jp = 0
          do 9000 jq = 1,5,2
            jp = jp+2
            call pbkt2(dg(1,ip),2,jp,t1)
            do 999 nn=1,83
  999       t3(nn) = 0.0d0
            call product(dg(1,iq),2,t1,1,t3)
            call product(dg(1,jq),2,t3,3,t5)
 9000     continue
  900   continue
c
        call pmadd(t5,5,-1.0d0/6.0d0,f)
c  Commutator term
        call pbkt1(g,3,g,4,temp)
        call pmadd(temp,5,-1.0d0/2.0d0,f)
      endif
      if ( imaxi .gt. 5 ) then
c  four derivative terms ( 1f4 x 2f3 )
c  ( See Dave's Notes )
c  term #1 and #2:
        do 66 iq = 1,5,2
          ip = iq+1
            do 666 jq = 1,5,2
              jp = jq+1
              do 6666 nn=1,209
 6666           t4(nn) = 0.0d0
              call pbkt2(dg(1,jq),2,ip,t1)
              call product(t1,1,dg(1,jp),3,t4)
              call pbkt2(dg(1,jp),3,ip,t2)
              call product(t2,2,dg(1,jq),2,t4)
c  add up
              call product(dg(1,iq),2,t4,4,tmp1)
  666     continue
   66   continue
c  term #3
        do 36 iq = 1,5,2
          ip=iq+1
          do 366 jq = 1,5,2
            jp = jq+1
            do 3666 nn = 1, 209
 3666         t4(nn) = 0.0d0
            call pbkt2(dg(1,iq),2,jq,t1)
            call product(t1,1,dg(1,jp),3,t4)
            call product(dg(1,ip),2,t4,4,tmp1)
  366     continue
   36   continue
c  add #1,#2,#2, with correct coefficient:
c
        call pmadd(tmp1,6,-1.0d0/2.0d0,f)
c
c  six  derivative terms ( 4f3 )
c
c
        do 1777 iq = 1,5,2
          ip = iq+1
          do 177 jq = 1,5,2
            jp = jq+1
            call pbkt2(dg(1,ip),2,jp,t1)
c  *** second term ( done with the first ).
            call pbkt2(crstrm,4,iq,temp)
            do 37 nn = 1, 83
   37         t3(nn) = 0.0d0
            call product(t1,1,dg(1,jq),2,t3)
            call product(t3,3,temp,3,dum)
c  *** end second term ******
            do 17 kq = 1,5,2
              kp = kq+1
              do 117 nn = 1,27
  117           t2(nn) = 0.0d0
              do 27 nn = 1,209
   27           t4(nn) = 0.0d0
              deriv3 = t1(kq)
              call pmadd(dg(1,iq),2,deriv3,t2)
              call product(dg(1,jq),2,t2,2,t4)
              call product(dg(1,kq),2,t4,4,t6)
   17       continue
  177     continue
 1777   continue
c *** 3rd term :
        do 88 iq = 1,5,2
          ip = iq+1
          call pbkt2(t5,5,iq,t4)
          call product(t4,4,dg(1,ip),2,t6)
   88   continue
c *** ****
c add second term with coeff 3. ):
        call pmadd(dum,6,3.0d0,t6)
c
        call pmadd(t6,6,-1.0d0/24.0d0,f)
c  End Six derivative terms
c
c  Commutator term:
        call pbkt1(g,3,g,5,temp)
        call pmadd(temp,6,-1.d0/2.0d0,f)
      endif
C  End six and 4 deriv terms **************
c
c  add in 2nd order piece
c
      f(8)=f(8)+1.d0
      f(19)=f(19)+1.d0
      f(26)=f(26)+1.d0
c
c   print out F (debug mode only)
c      write(6,*) '  The generating function F is : '
c      do 77 i=1,923
c         if (dabs(f(i)).gt.1.d-10 ) write(6,*) i,f(i)
c  77  continue
c
c=======================================================
c this is a quick hack to compare these results w/ 3.0:
c     write(6,*)'(canx) hack'
c     f(210:monoms)=0.d0
c 8/15/02 verified that, for example 2_7, the symplectic
c tracking results are the same for ML3.0 and ML5.0 with
c the above statement enabled. RDR
c
cryne 12/31/2004:
      if(norder.eq.4)f(462:monoms)=0.d0
      if(norder.eq.3)f(210:monoms)=0.d0
      if(norder.eq.2)f(84:monoms) =0.d0
c=======================================================
c
c  now compute the derivs of the generating fxn f (which define
c  the standard representation of the canonical transformation)
c
      do 5000 i=1,5,2
         do 500 k=2,imaxi
            do 5 l=1,len(imaxi)
               tmp1(l)=0.d0
               tmp2(l)=0.d0
    5       continue
            ivalue=i+1
            call pbkt2(f,k,ivalue,tmp1)
            call pbkt2(f,k,i,tmp2)
            do 50 l=1,len(imaxi)
               df(i,l)=df(i,l)+tmp1(l)
               df(ivalue,l)=df(ivalue,l)-tmp2(l)
   50       continue
  500    continue
 5000 continue
c print out the derivatives of F (debug only)
c      write(6,*) ' The derivatives of F are :'
c      do 771 i=1,6
c        write(6,772) i
c 772    format('  Derivatives with repect to z',i1, ' are:')
c        do 773 l=1,923
c          if(dabs(df(i,l)).gt.1.d-10 ) write(6,*) l,df(i,l)
c 773    continue
c 771  continue
c
c  now compute the jacobian of the momentum part of the
c  standard rep (which will be used in sub. newton to determine
c  the new momentum)
c
      ivalue=0
      do 6000 i=1,5,2
         do 6 n=1,len(imaxi)
            dum(n)=df(i,n)
    6    continue
         ivalue=ivalue+1
         jvalue=0
         do 600 j=1,5,2
            jvalue=jvalue+1
            do 60 k=2,imaxi-1
               do 58 l=1,len(imaxi)
                  temp(l)=0.d0
   58          continue
               call pbkt2(dum,k,j,temp)
               do 59 l=1,len(imaxi)
                  rjac(ivalue,jvalue,l)=rjac(ivalue,jvalue,l)+temp(l)
   59          continue
   60       continue
  600    continue
 6000 continue
      return
      end
c
***********************************************************************
c
      subroutine invrs(a,b,d)
c  subroutine to compute the inverse b and determinant d
c  of a 3x3 matrix a.  used by newton search subroutine
c  when searching for the new momentum.
c     Written by D. Douglas, ca 1983
      include 'impli.inc'
      include 'files.inc'
      dimension a(3,3),b(3,3),c(3,3)
      do 10 i=1,3
         do 1 j=1,3
            b(i,j)=0.d0
            c(i,j)=0.d0
    1    continue
   10 continue
c
c  compute determinant of a
c
      d=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)
     & -a(1,3)*a(2,2)*a(3,1)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)
c     write(jof,9100)d
c9100 format(1h ,'determinant of a is ',d21.15)
      abd=dabs(d)
      if(abd.lt.1.d-30)write(6,9000)
 9000 format(1h ,'determinant underflow in subroutine invrs')
      if(abd.lt.1.d-30)return
c
c  compute b=inverse of a
c
      b(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/d
      b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))/d
      b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))/d
      b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
      b(1,3)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
      b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
      b(2,3)=(a(2,1)*a(1,3)-a(1,1)*a(2,3))/d
      b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
      b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
c
c  check if this is a reasonable inverse
c
c     do 200 i=1,3
c        do 20 j=1,3
c           do 2 k=1,3
c              c(i,j)=c(i,j) + a(i,k)*b(k,j)
c   2       continue
c  20    continue
c 200 continue
c     write(6,9150)
c9150 format(1h ,'a * ainverse = '/)
c     do 3 i=1,3
c        write(6,9200)(c(i,j),j=1,3)
c9200    format(1h ,3(d21.15,2x))
c   3 continue
      return
      end
c
***********************************************************************
c
      subroutine newt(zi,zo,rdf,rrjac)
cryne 8/15/02      subroutine newt(zi,zo)
c  subroutine to invert map p(old)=df(q(old),p(new))/dq(old)
c  for p(new) using a newton's search procedure
c
c     Written by D. Douglas, ca 1983 and substantially modified
c  by F. Neri
c
      use lieaparam, only : monoms
      include 'impli.inc'
cryne 8/15/02      include 'deriv.inc'
cryne 8/15/02      dimension rdf(3,monom1+1),rrjac(3,3,monom2+1)
      dimension rdf(3,monoms),rrjac(3,3,monoms)
      dimension ztmp(6),zi(6),zo(6)
      dimension rlin(3,3),rmat(3,3),ri(3,3),rinv(3,3)
      dimension pi(3),crtrm(3),pimag(3),p(3),delp(3),cp(3)
      dimension pjac(3,3,84),pdf(3,84)
      dimension qvec(84),pvec(84)
      include 'files.inc'
      include 'ind.inc'
      include 'len.inc'
      include 'len3.inc'
      include 'ind3.inc'
      include 'talk.inc'
      data ri/1.,3*0.,1.,3*0.,1./
cryne 5/22/02
      save ri
c
c  initialize arrays
c
      l31=len3(imaxi-1)
      l32=len3(imaxi-2)
      ivalue=0
      do 1 i=1,3
         ivalue=ivalue+2
         pi(i)=zi(ivalue)
    1 continue
      do 2 i=1,6
         ztmp(i)=zi(i)
    2 continue
c
c generate monomials in the q's
c
      qvec(1) = 1.0d0
      ivalue=-1
      do 10 i=1,3
      ivalue=ivalue+2
      qvec(i+1)=ztmp(ivalue)
 10   continue
      do 11 i=5,len3(imaxi-1)+1
      qvec(i)=qvec(ind31(i))*qvec(ind32(i))
 11   continue
c
c  generate coefficient of polynomials in the momentum
c  variables
c
      do 1001 i1=1,3
      do 1001 i2=1,3
      pjac(i1,i2,1)=0.d0
 1001 continue
c
      do 1002 n=1,4
      do 1003 i=1,3
      pdf(i,n) = 0.d0
 1003 continue
 1002 continue
c
c generate coefficients of reduced df (pdf),function only of p:
c***********************
      do 190 i = 1,3
      do 190 ip = 5,len3(imaxi-1)+1
      pdf(i,ip) = rdf(i,ip)
  190  continue
c***** previous loops replaced by: *************
c     call dcopy(3*(len3(imaxi-1)-3),rdf(1,5),1,pdf(1,5),1)
c**********************
      n0 = len3(imaxi-1)+1
      iq2=1
      do 9 ipoq=1,imaxi-2
      iq1=iq2+1
      iq2=len3(ipoq)+1
      if(imaxi.eq.ipoq+1) l = 1
      if(imaxi.gt.ipoq+1) l = len3(imaxi-1-ipoq)+1
      do 99 iq=iq1,iq2
      qval = qvec(iq)
c********************
      do 191 i = 1,3
      do 191 ip = 1,l
      pdf(i,ip) = pdf(i,ip) + rdf(i,n0+ip)*qval
  191  continue
c**** previous loops replaced by: **********
c     call daxpy(3*l,qval,rdf(1,n0+1),1,pdf(1,1),1)
c********************
          n0 = n0 + l
 99      continue
 9      continue
c********************
      do 999 i = 1,3
      do 999 ip = 1,l31-l32
      pdf(i,1) = pdf(i,1) + rdf(i,n0+ip)*qvec(l32+1+ip)
  999  continue
c**** previous loops replaced by: **********
c     pdf(1,1)=pdf(1,1)+ddot(l31-l32,qvec(l32+2),1,rdf(1,n0+1),3)
c     pdf(2,1)=pdf(2,1)+ddot(l31-l32,qvec(l32+2),1,rdf(2,n0+1),3)
c     pdf(3,1)=pdf(3,1)+ddot(l31-l32,qvec(l32+2),1,rdf(3,n0+1),3)
c********************
c
c   genarate coefficients of reduce  jacobian (rrjac):
c********************
      do 1190 i1 = 1,3
      do 1190 i2 = 1,3
      do 1190 ip = 2,len3(imaxi-2)+1
      pjac(i1,i2,ip) = rrjac(i1,i2,ip)
 1190 continue
c**** previous loops replaced by: **********
c     call dcopy(9*len3(imaxi-2),rrjac(1,1,2),1,pjac(1,1,2),1)
c********************
      n0 =  len3(imaxi-2)+1
      iq2=1
      do 19 ipoq = 1,imaxi-2
      iq1 = iq2 + 1
      iq2 = len3(ipoq)+1
      if(imaxi.eq.ipoq+2) l = 1
      if(imaxi.gt.ipoq+2) l = len3(imaxi-2-ipoq)+1
      do 199 iq = iq1,iq2
      qval = qvec(iq)
c*********************
      do 1999 i1 = 1,3
      do 1999 i2 = 1,3
      do 1999 ip = 1,l
      pjac(i1,i2,ip) = pjac(i1,i2,ip) + rrjac(i1,i2,n0+ip)*qval
 1999 continue
c**** previous loops replaced by: **********
c     call daxpy(9*l,qval,rrjac(1,1,n0+1),1,pjac(1,1,1),1)
c********************
          n0 = n0 + l
 199  continue
 19   continue
c
c
c  loop to apply contraction mapping
c
      do 9999 index=1,12
c
         do 30 i=1,3
            crtrm(i)=0.d0
            pimag(i)=0.d0
            p(i)=0.d0
            delp(i)=0.d0
            cp(i)=0.d0
            do 3 j=1,3
               rlin(i,j)=0.d0
    3       continue
   30    continue
c   compute monomials in the momenta
c
      pvec(1) = 1.d0
c  first order monomials are equal to the p's;
      ivalue=0
      do 4 i=1,3
      ivalue=ivalue+2
      pvec(i+1)=ztmp(ivalue)
 4    continue
c
c  compute the remaining monomials as products of
c  previous ones:
c  need only term up to order imaxi-1
      do 40 i=5,len3(imaxi-1)+1
      pvec(i)=pvec(ind31(i))*pvec(ind32(i))
 40   continue
c
c     compute jacobian
c
         do 500 i=1,3
            do 50 j=1,3
c only terms of order up to imaxi-2 are present in the jacobian,
c the jacobian is produced by taking the second derivatives of a
c imaxi order polynomial.
c********************
           do 5 n =1,len3(imaxi-2)+1
                rlin(i,j) = rlin(i,j) + pjac(i,j,n)*pvec(n)
    5       continue
c**** previous loop replaced by: **********
c        rlin(i,j)=ddot(l32+1,pvec,1,pjac(i,j,1),9)
c********************
c
c
               rmat(i,j)=ri(i,j)-rlin(i,j)
   50       continue
  500    continue
         call invrs(rmat,rinv,det)
         adet=dabs(det)
         if(adet.lt.1.d-30)write(6,501)
  501    format(1h ,'determinant underflow in newt')
         if(adet.lt.1.d-30)return
         ivalue=0
         do 6 i=1,3
            ivalue=ivalue+2
            p(i)=ztmp(ivalue)
    6    continue
         do 600 i=1,3
            pimag(i)=pimag(i)+pi(i)
c only term of order 2 to imaxi-1 are present
c  so the sum only goes up to len(imaxi-1) .   not
c********************
             do 60 n=1,l31+1
                   pimag(i) = pimag(i) - pdf(i,n)*pvec(n)
   60     continue
c**** previous loop repplaced by: **********
c       pimag(i)=pimag(i)-ddot(l31+1,pvec,1,pdf(i,1),3)
c********************
c
c
            delp(i)=p(i)-pimag(i)
  600    continue
         do 70 i=1,3
            do 7 j=1,3
               crtrm(i)=crtrm(i)+rinv(i,j)*delp(j)
    7       continue
   70    continue
c         square=crtrm(1)**2 + crtrm(2)**2 + crtrm(3)**2
c         root=dsqrt(square)
         root=abs(crtrm(1)) + abs(crtrm(2)) + abs(crtrm(3))
         if(ibrief.ne.2)goto 705
         write(jof,71)index,(crtrm(i),i=1,3)
   71    format(1h ,'corrections at iteration ',i2/
     &   1h ,'crtrm(1) = ',d21.15/1h ,'crtrm(2) = ',d21.15/
     &   1h ,'crtrm(3) = ',d21.15)
c         if (1.+root.eq.1.) write(jof,72)index
         if (root .le. 1.d-12) write(jof,72)index
   72    format(1h ,'iteration has converged; i= ',i1)
  705    continue
         do 8 i=1,3
            cp(i)=p(i)-crtrm(i)
    8    continue
         ivalue=0
         do 80 i=1,3
            ivalue=ivalue+2
            ztmp(ivalue)=cp(i)
   80    continue
         do 800 i=1,6
            zo(i)=ztmp(i)
  800    continue
c         if(1.+root.eq.1.) goto 1234
         if(root .le. 1.d-12) goto 1234
 9999 continue
      write(6,906) root
  906 format(' Search did not converge; root= ',e12.5)
c      if (root.gt..1) jwarn=1
      jwarn=1
 1234 continue
      return
      end
c
***********************************************************************
c
      subroutine prod(dg,crstrm)
c  subroutine to compute crossterm in generating fxn of standard
c  rep of transfer map
c     Written by D. Douglas, ca 1982, and modified by Liam Healy
c
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension dg(6,*),crstrm(*),l(6)
      include 'expon.inc'
      include 'lims.inc'
c
      do 1 m=1,top(4)
         crstrm(m)=0.d0
    1 continue
      do 2000 i=1,5,2
         ivalue=i+1
         do 200 j=bottom(2),top(2)
            if(dg(i,j).eq.0.d0)goto 200
            do 20 k=bottom(2),top(2)
               if(dg(ivalue,k).eq.0.d0)goto 20
               do 2 m=1,6
                  l(m)=expon(m,j)+expon(m,k)
    2          continue
               n=ndex(l)
               crstrm(n)=crstrm(n) - dg(i,j)*dg(ivalue,k)/2.d0
   20       continue
  200    continue
 2000 continue
      return
      end
c
***********************************************************************
      subroutine rearr
c Written by F. Neri, ca 1984
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension j(6)
      include 'ja3.inc'
      include 'ind.inc'
      include 'deriv.inc'
      include 'len.inc'
      include 'len3.inc'
      include 'ind3.inc'
c
      do 22 i=1,6
      j(i)=0
 22   continue
      do 1 ip = 2,len3(imaxi-1)+1
      ivalue = 0
      do 11 i = 1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,ip)
 11   continue
      n = ndex(j)
      if(ip.gt.(len3(imaxi-2)+1)) goto 112
      do 111 i1 = 1,3
      do 111 i2 = 1,3
      rrjac(i1,i2,ip) = rjac(i1,i2,n)
 111  continue
 112   continue
      ivalue = -1
      do 10 i = 1,3
      ivalue = ivalue + 2
      rdf(i,ip) = df(ivalue,n)
 10   continue
 1    continue
      n1 = len3(imaxi-1)+1
      n2 = len3(imaxi-2)+1
      iq2 = 1
      do 3 ipoq = 1,imaxi-1
      iq1 = iq2+1
      iq2 = len3(ipoq)+1
      if(imaxi-2.gt.ipoq) then
      l2 = len3(imaxi-2-ipoq)+1
      l1 = len3(imaxi-1-ipoq)+1
      else if(imaxi-2.eq.ipoq) then
      l2 = 1
      l1 = len3(imaxi-1-ipoq)+1
      else
      l1 = 1
      end if
      do 33 iq = iq1,iq2
      if(imaxi-2.ge.ipoq) then
      do 333 ip = 1,l2
      ivalue = -1
      do 3333 i=1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,iq)
      j(ivalue+1) = ja3(i,ip)
 3333 continue
      n = ndex(j)
      do 4 i1=1,3
      do 4 i2 = 1,3
      rrjac(i1,i2,n2+ip) = rjac(i1,i2,n)
 4    continue
 333  continue
      n2 = n2 + l2
      endif
      do 334 ip = 1,l1
      ivalue = -1
      do 3334 i = 1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,iq)
      j(ivalue+1) = ja3(i,ip)
 3334 continue
      n = ndex(j)
      ivalue = -1
      do 44 i = 1,3
      ivalue = ivalue+2
      rdf(i,n1+ip)  = df(ivalue,n)
 44   continue
 334  continue
      n1 =  n1 + l1
 33   continue
 3    continue
      return
      end
c
***********************************************************************
c
      subroutine trace(icfile,nunit,ntaysym,norder,ntrace,nwrite,          &
     &                 jdmin,jdmax,nprecision,nraysinp,th,tmh)
c     Written by D. Douglas, ca 1982, and modified since by
c     nearly everyone at Maryland.  Could still stand improvement.
c
      use rays
      use beamdata
      use lieaparam, only : monoms
      use multitrack
      use ml_timer
      include 'impli.inc'
      include 'deriv.inc'
      include 'files.inc'
      include 'talk.inc'
      dimension th(monoms),tmh(6,6)
      character*16 ubuf
cryne 12/30/2004 moved idmin/idmax statements below past raysin
cryne because calling raysin can potentially change maxray
!     data iifile/50/
!     save iifile
c
!     write(6,*)'(trace) hello from PE ',idproc
!     if(idproc.eq.0)then
!     write(6,*)'in routine trace; icfile,nunit,norder,ntrace,nwrite='
!     write(6,*)icfile,nunit,norder,ntrace,nwrite
!     endif
c
c  read initial conditions, if requested:
      if(icfile.gt.0)then
        scaleleni=0.d0
        scalefrqi=0.d0
        scalemomi=0.d0
        call increment_timer('raysin',0)
        call raysin(icfile,nraysinp,scaleleni,scalefrqi,scalemomi,       &
     &              ubuf,jerror)
        call increment_timer('raysin',1)
        call rbcast(scaleleni)
        call rbcast(scalefrqi)
        call rbcast(scalemomi)
c  if all the scaling constants are zero, then they are not in the prologue,
c  so the sectio of code below on units conversion should be skipped:
        if(scaleleni.eq.0.d0.and.scalefrqi.eq.0.d0.and.scalemomi.eq.0.d0&
     &  )goto 555
c
c  if one of the scaling constants are zero, but not all of them,
c  then there is some sort of problem, so exit:
        scaleprod=scaleleni*scalefrqi*scalemomi
        if(scaleprod.eq.0)then
        if(idproc.eq.0)then
        write(6,*)'problem reading scaling constants in input file'
        write(6,*)'is one of them zero???'
        endif
        call myexit
        endif
c    here is where the units conversion takes place:
        iconvert=1
        if(iconvert.eq.1)then
          slratio=scaleleni/sl
          smomratio=scalemomi/p0sc
          timratio=freqscl/scalefrqi
          wldratio=(scalefrqi*scaleleni*scalemomi)/(freqscl*sl*p0sc)
          if(idproc.eq.0)then
          write(6,*)'raysin complete; scaling input data'
c         write(6,*)'scaleleni,sl=',scaleleni,sl
c         write(6,*)'scalemomi,p0sc=',scalemomi,p0sc
c         write(6,*)'scalefrqi,freqscl=',scalefrqi,freqscl
c         write(6,*)'twopi*freqscl,omegascl=',4.d0*asin(1.d0)*freqscl,       &
c    &    omegascl
c         write(6,*)'wld_i,wld=',(scalefrqi*scaleleni*scalemomi),            &
c    &    (freqscl*sl*p0sc)
c         write(6,*)' '
c         write(6,*)' '
          write(6,*)'scale length ratio=',slratio
          write(6,*)'scale momentum ratio=',smomratio
          write(6,*)'scale time ratio=',timratio
          write(6,*)'scale energy ratio=',wldratio
          endif
          do n=1,nraysp
            zblock(1,n)=zblock(1,n)*slratio
            zblock(2,n)=zblock(2,n)*smomratio
            zblock(3,n)=zblock(3,n)*slratio
            zblock(4,n)=zblock(4,n)*smomratio
            zblock(5,n)=zblock(5,n)*timratio
            zblock(6,n)=zblock(6,n)*wldratio
          enddo
        endif
  555 continue
      endif
c
!---new code to print particles in range jdmin to jdmax:
      idmin=jdmin
      idmax=jdmax
      if(idmin.eq.0)idmin=1
      if(idmax.eq.0)idmax=maxray
!---
c
c  examine various cases for norder
c  norder=-1
      if(norder.eq.-1.and.icfile.gt.0) return
      if(norder.eq.-1.and.icfile.eq.0) then
        write(6,*) 'icfile and norder have inconsistent values'
        write(12,*) 'icfile and norder have inconsistent values'
        call myexit
      endif
c
c  norder=0
c  write final conditions in full precision
!     if(norder.eq.0.and.jfcf.lt.0)then
      if(norder.eq.0.and.nunit.lt.0)then
        write(6,*)'FULL PRECISION WRITE ON FILE ',jfcf
        call pwritezfp
      endif
cryne      if(norder.eq.0.and.jfcf.lt.0)then
cryne        do 120 i=1,nraysp
cryne        do 110 j=1,6
cryne        write(-jfcf,2468)zblock(j,i)
 2468 format(1x,d32.25,1x,d32.25)
cryne  110   continue
cryne  120   continue
cryne        write(jof,130) -jfcf
cryne  130   format(1x,'final conditions written in full precision on',
cryne     #  ' file ',i4)
cryne      endif
c
c  write final conditions in standard format
!     if(norder.eq.0.and.jfcf.gt.0)call pwritez(nunit,idmin,idmax,           &
      if(norder.eq.0.and.nunit.gt.0)call pwritez(nunit,idmin,idmax,          &
     &   nprecision,0,0)
cryne      if(norder.eq.0.and.jfcf.gt.0)then
cryne        do 135,i=1,nraysp
cryne        write(jfcf,136)(zblock(j,i),j=1,6)
cryne  136   format(6(1x,1pe12.5))
cryne  135   continue
cryne      endif
      if(norder.eq.0)return
c
c  cases where norder > 0
c
c  call preliminary routines depending on type of ray trace:
cryne May 10, 2006 added "if" test on norder  ===============
      if(norder.gt.1)then
        if(ntaysym.eq.1)then
           if(multitrac.eq.0)then
             call brkts(th)
           else
             call multibrkts
           endif
        else
cryne 12/31/2004 added norder to canx argument list.
cryne This is a temporary fix to allow symplectic tracking of order < 5.
cryne It would be better if do loops did not extend to 923 in evalsr.
          if(multitrac.eq.0)then
            call canx(tmh,th,norder)
            call rearr
          else
            call multicanx
          endif
        endif
      endif
c
c  setup before entering main loop:
      ktrace=ntrace
      if(ktrace.eq.0)ktrace=1
c
c the 2 statements below seem to cause a problem and have been commented out:
c      kwrite=nwrite
c      if(kwrite.eq.0)kwrite=1
c
cryne May 22, 2006 it is very bad and confusing that the names th and tmh
cryne used in the argument list to this routine, since they might/might not be
cryne the total transfer map coefficients stored in common under these names!!!
      call increment_timer('eval',0)
      if(ntaysym.eq.1)then
        call eval(tmh,norder,ntrace,nwrite,nunit,idmin,idmax,nprecision)
        call increment_timer('eval',1)
        return
      else
        jwarn=0
        call evalsr(tmh,df,rdf,rrjac,norder,ntrace,nwrite,nunit,            &
     &              idmin,idmax,nprecision)
        call increment_timer('eval',1)
        return
      endif
cryne 8/3/2002 original code still in place for symplectic ray trace:
c main loop; do 'ktrace' ray traces (of the 'nraysp' initial rays)
      if(idproc.eq.0)write(6,*)'ERROR: CODE SHOULD NOT GET HERE!!!'
      do 300 idum=1,ktrace
      do 200 k=1,nraysp
ccc   if (nlost.ge.nrays) then
ccc     write (jof,*) 'all particles lost'
ccc     return
ccc   endif
c check to see if particle has already been lost
         if (istat(k).ne.0) goto 200
c
         do 150 l=1,6
  150    zi(l)=zblock(l,k)
c
cryne 8/3/2002         if(norder.le.4)then
cryne 8/3/2002           call eval(tmh,norder,zi,zf)
cryne 8/3/2002         endif
         if(norder.gt.4)then
c perform a symplectic ray trace
           jwarn=0
           call evalsr_old(tmh,zi,zf,df,rdf,rrjac)
cryne 8/13/02           call evalsr(tmh,zi,zf)
c          if(jfcf.lt.0)write(6,*)'back from evalsr, jfcf=',jfcf
c          write(6,*)'returned from evalsr for ray #',k
c check to see if particle was 'lost' by the symplectic ray tracer
c and mark so 'lost' particles by a a distinctive negative number
           if (jwarn.ne.0) then
c            write(6,*)'setting istat,nlost'
             istat(k)=-idum
             nlost=nlost+1
             ihist(1,nlost)=-idum
             ihist(2,nlost)=k
c            write(6,*)'done setting istat,nlost'
           endif
         endif
c  the line below has been replaced by the following two lines
c         if(mod(idum,kwrite).ne.0)goto 176
          if(nwrite.eq.0)goto 176
          if(mod(idum,nwrite).ne.0)goto 176
c          if(jfcf.lt.0)write(6,*)'did not goto 176;jfcf=',jfcf
c  end of modifications
c  write out final conditions
c  full precision case
         if(jfcf.lt.0) then
           write(6,*)'writing in full prec, idum=',idum
cryne June 4 2002           do 152 l=1,6
cryne June 4 2002           write(-jfcf,*) zf(l)
cryne June 4 2002  152      continue
        do 152 l=1,6,2
        write(-jfcf,2468)zf(l),zf(l+1)
  152 continue
           write(6,*)'donewriting in full prec, idum=',idum
         endif
c  standard format case
         if(jfcf.gt.0)then
         tblock(:,k)=zf(:)
         endif
         if(ibrief.eq.0)goto 176
c  print lengthy info:
         write(jof,155)
  155    format(/1h ,'initial conditions are: (dimensionless form)')
         write(jof,*)zi(1),zi(2)
         write(jof,*)zi(3),zi(4)
         write(jof,*)zi(5),zi(6)
         write(jof,160)
  160    format(/1h ,'final conditions are: (dimensionless form)')
         write(jof,*)zf(1),zf(2)
         write(jof,*)zf(3),zf(4)
         write(jof,*)zf(5),zf(6)
c  (end of lengthy info)
c
  176    continue
         if(ntrace.eq.0)goto 200
         zblock(:,k)=zf(:)
c
  200 continue
      if(nwrite.ne.0)then
        if(mod(idum,nwrite).eq.0 .and. jfcf.gt.0)then
        call pwritet(nunit,idmin,idmax,nprecision)
        endif
      endif
  300 continue
      call increment_timer('eval',1)
      return
      end
c---------------------------------------
      subroutine pwritez(nunit,idmin,idmax,nprecision,nphysunits,         &
     &                   iprintarc)
c routine to perform parallel write of zblock array
      use rays
      use beamdata
      use lieaparam, only : monoms  !cryne 9/23/07
      include 'impli.inc'
      include 'map.inc'   !cryne 9/23/07 need to know arclen
      integer l,nraysw
      include 'files.inc'
      character*16 :: myformat='(7(1x,1pe12.5 ))'
      character*2 :: a2
      if(idproc.eq.0)then
!---format:
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
!---
        mycount=0
        do i=1,nraysp
         mycount=mycount+1
         if(mycount.lt.idmin)cycle
         if(mycount.gt.idmax)exit
         if(nphysunits.eq.0)then
           if(iprintarc.eq.0)write(nunit,myformat)zblock(:,i)
           if(iprintarc.eq.1)write(nunit,myformat)arclen,zblock(:,i)
         else
           if(iprintarc.eq.0)                                              &
     &     write(nunit,myformat)zblock(1,i)*sl,zblock(2,i)*p0sc,           &
     &     zblock(3,i)*sl,zblock(4,i)*p0sc,zblock(5,i)/omegascl,           &
     &     zblock(6,i)*(omegascl*sl*p0sc)
           if(iprintarc.eq.1)                                              &
     &     write(nunit,myformat)arclen,zblock(1,i)*sl,zblock(2,i)*p0sc,    &
     &     zblock(3,i)*sl,zblock(4,i)*p0sc,zblock(5,i)/omegascl,           &
     &     zblock(6,i)*(omegascl*sl*p0sc)
         endif
        enddo
        do l=1,nvp-1
          call MPI_RECV(tblock,6*maxrayp,mreal,l,98,lworld,mpistat,ierr)
          call MPI_GET_COUNT(mpistat,mreal,nraysw,ierr)
          nraysw=nraysw/6
          do i=1,nraysw
           mycount=mycount+1
           if(mycount.lt.idmin)cycle
           if(mycount.gt.idmax)exit
           if(nphysunits.eq.0)then
             if(iprintarc.eq.0)write(nunit,myformat)tblock(1:6,i)
             if(iprintarc.eq.1)write(nunit,myformat)arclen,tblock(1:6,i)
           else
           if(iprintarc.eq.0)                                              &
     &     write(nunit,myformat)tblock(1,i)*sl,tblock(2,i)*p0sc,           &
     &     tblock(3,i)*sl,tblock(4,i)*p0sc,tblock(5,i)/omegascl,           &
     &     tblock(6,i)*(omegascl*sl*p0sc)
           if(iprintarc.eq.1)                                              &
     &     write(nunit,myformat)arclen,tblock(1,i)*sl,tblock(2,i)*p0sc,    &
     &     tblock(3,i)*sl,tblock(4,i)*p0sc,tblock(5,i)/omegascl,           &
     &     tblock(6,i)*(omegascl*sl*p0sc)
           endif
  100      format(6(1x,1pe12.5))
          enddo
        enddo
      else
        call MPI_SEND(zblock,6*nraysp,mreal,0,98,lworld,ierr)
      endif
      return
      end
c---------------------------------------
      subroutine pwritet(nunit,idmin,idmax,nprecision)
c routine to perform parallel write of tblock array
      use rays
!     implicit none
      include 'impli.inc'
      integer l,nraysw
      include 'files.inc'
      character*16 :: myformat='(6(1x,1pe12.5 ))'
      character*2 :: a2
      if(idproc.eq.0)then
!---format:
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
!       write(6,*)'nlength,nprecision=',nlength,nprecision
!       write(6,*)'myformat=',myformat
!---
        mycount=0
        do i=1,nraysp
         mycount=mycount+1
         if(mycount.lt.idmin)cycle
         if(mycount.gt.idmax)exit
!        write(nunit,100)tblock(:,i)
         write(nunit,myformat)tblock(:,i)
        enddo
        do l=1,nvp-1
          call MPI_RECV(tblock,6*maxrayp,mreal,l,98,lworld,mpistat,ierr)
          call MPI_GET_COUNT(mpistat,mreal,nraysw,ierr)
          nraysw=nraysw/6
!         write(6,*)'(pwrite) PE 0 received ',nraysw,' rays from PE ',l
          do i=1,nraysw
           mycount=mycount+1
           if(mycount.lt.idmin)cycle
           if(mycount.gt.idmax)exit
!          write(nunit,100)tblock(1:6,i)
           write(nunit,myformat)tblock(1:6,i)
  100      format(6(1x,1pe12.5))
          enddo
        enddo
      else
        call MPI_SEND(tblock,6*nraysp,mreal,0,98,lworld,ierr)
      endif
      return
      end
c---------------------------------------
      subroutine pwritezfp
c routine to perform parallel write of zblock array in FULL PRECISION
      use rays
!     implicit none
      include 'impli.inc'
      integer l,nraysw
      include 'files.inc'
      if(idproc.eq.0)then
        do i=1,nraysp
         write(jfcf,*)zblock(1,i),zblock(2,i)
         write(jfcf,*)zblock(3,i),zblock(4,i)
         write(jfcf,*)zblock(5,i),zblock(6,i)
        enddo
        do l=1,nvp-1
          call MPI_RECV(tblock,6*maxrayp,mreal,l,98,lworld,mpistat,ierr)
          call MPI_GET_COUNT(mpistat,mreal,nraysw,ierr)
          nraysw=nraysw/6
!         write(6,*)'(pwrite) PE 0 received ',nraysw,' rays from PE ',l
          do i=1,nraysw
           write(jfcf,*)tblock(1,i),tblock(2,i)
           write(jfcf,*)tblock(3,i),tblock(4,i)
           write(jfcf,*)tblock(5,i),tblock(6,i)
          enddo
        enddo
      else
        call MPI_SEND(zblock,6*nraysp,mreal,0,98,lworld,ierr)
      endif
      return
      end
c end of file
c---------------------------------------
      subroutine pwritetfp
c routine to perform parallel write of tblock array
      use rays
!     implicit none
      include 'impli.inc'
      integer l,nraysw
      include 'files.inc'
      if(idproc.eq.0)then
        do i=1,nraysp
         write(jfcf,*)tblock(1,i),tblock(2,i)
         write(jfcf,*)tblock(3,i),tblock(4,i)
         write(jfcf,*)tblock(5,i),tblock(6,i)
        enddo
        do l=1,nvp-1
          call MPI_RECV(tblock,6*maxrayp,mreal,l,98,lworld,mpistat,ierr)
          call MPI_GET_COUNT(mpistat,mreal,nraysw,ierr)
          nraysw=nraysw/6
!         write(6,*)'(pwrite) PE 0 received ',nraysw,' rays from PE ',l
          do i=1,nraysw
           write(jfcf,*)tblock(1,i),tblock(2,i)
           write(jfcf,*)tblock(3,i),tblock(4,i)
           write(jfcf,*)tblock(5,i),tblock(6,i)
          enddo
        enddo
      else
        call MPI_SEND(tblock,6*nraysp,mreal,0,98,lworld,ierr)
      endif
      return
      end
c---------------------------------------
c end of file
