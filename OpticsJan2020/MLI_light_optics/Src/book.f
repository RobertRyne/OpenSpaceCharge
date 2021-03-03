************************************************************************
* header:                 BOOKKEEP
*  Index manipulation, table creation, block data, startup routines
************************************************************************
c
      subroutine binom
c
c     computes the table of the binomial coefficients
c
      implicit double precision (a-h,o-z)
       integer bin(24,20)
      common /bin/ bin
      do 1 i=1,20
      bin(i,1)=i
      do 1 k=2,20
      if (i-k) 2,3,4
    2 bin(i,k)=0
      go to 1
    3 bin(i,k)=1
      go to 1
    4 ip=i-1
      kp=k-1
      bin(i,k)=bin(ip,kp)+bin(ip,k)
    1 continue
      return
      end
c
      integer function ndexn(num,j)
c  calculates the index of the long vector (dimension 209)
c  for a given set of short indices that indicate a power of
c  each variable, stored in the array j.
c  This function is a combination of DRD's subroutines
c  'indice' and 'binom' in Marylie.
c  Written by Liam Healy, June 8, 1984.
c
c  'bin(i,k)' is the binomial coefficient i select k (i>k)
cryne integer bin(24,20),j(6)
      integer bin(24,20),j(*)
      common/bin/ bin
c  calculate the index
      n=j(num)
      m=j(num)-1
      do 100 i=2,num
      ib=num+1-i
      m=m+j(ib)
      ib=m+i
      n=n+bin(ib,i)
  100 continue
      ndexn=n
      return
      end
c
      integer function ndex(j)
      integer j(6)
      ndex=ndexn(6,j)
c
      return
      end
c
***********************************************************************
c
      subroutine setup
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'ind.inc'
cryne 5/3/2006      include 'files.inc'
cryne 5/3/2006      include 'time.inc'
c
c initialize commons (other than those in block data's)
      imaxi = 6
      call initia
c
c initialize lie algebraic things
      call binom
      call tables
      call init
c
      return
      end
c
***********************************************************************
c
      subroutine cpyrt
c
      use parallel
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
cryne 5/3/2006      include 'time.inc'
      include 'ind.inc'
c
c write out copyright message at terminal
c
      if(idproc.eq.0)then
      write(jof,90)
   90 format(
     &'***************************************************************',
     &/,
     &'*                       MARYLIE/IMPACT                        *',
     &/,
     &'* Parallel Lie Algebraic Beam Dynamics Code with Space Charge *',
     &/,
     &'*              Last modified May 22, 2006  08:20PDT           *',
     &/,
     &'*     MaryLie copyright 1987, Alex J. Dragt, U. Maryland      *',
     &/,
     &'*     IMPACT copyright 2002, Robert D. Ryne, U. California    *',
     &/,
     &'***************************************************************')
      endif
c
      return
      end
c
******************************************************************************
c
      subroutine initia
c-----------------------------------------------------------------------
c  This routine initializes some constants in common blocks, which are
c  not initialized in block data
c
c  Petra Schuett  10/30/87
c  Alex Dragt 6/20/88
c-----------------------------------------------------------------------
c
cryne 7/23/2002 this is here because of the speed of light. fix later.
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pie.inc'
      include 'frnt.inc'
c
c set up constants
c
c      pi    = 3.14159265358979323846264d0
      pi = 4.d0*atan(1.d0)
      pi180 = pi/180.d0
      twopi = 2.d0*pi
      c     = 2.99792458d+08
c
c set up default values for cfbd fringe field parameters
c
      cfbgap=0.d0
      cfblk1=.5
      cfbtk1=.5
c
      return
      end
c
***********************************************************************
c
      subroutine tables
c  This subroutine creates two tables used in bookkeeping:
c    1) expon(ind,psv) is the exponent of phase space variable
c       'psv' (1 to 6)  for monomial index 'ind' (1 to top)
c    2) vblist(ind,vnum) is the variable list for each
c       monomial index number 'ind' (1 to top).  For nth
c       order terms, there will be n non-zero variable
c       numbers.
c  For example, monomial number 109 is X.PX.PX.Pt.
c  Thus, expon(109,1to6)=1,2,0,0,0,1 and vblist(109,1to4)=1,2,2,6.
c  This subroutine was written by Liam Healy on June 20,1984,
c  and is a replacement for Dave Douglas's subroutines 'addss' and 'sht'
c
      use lieaparam, only : monoms,monom1
      implicit none
c  ------Variables produced by subroutine-------
      include 'expon.inc'
      include 'vblist.inc'
      include 'prodex.inc'
      include 'order.inc'
      include 'maxcat.inc'
      include 'iprod.inc'
c  ------Variables used internally--------
c  carry = carry amount from j(6) in counting for expon
c  ind = monomial number index
c  lnzj = last non-zero j
c  psv = phase space variable 1...6 corresponds to X...Pt
c  vnum = current postion for writing variable number in 'vblist'
      integer carry,ind,lnzj,psv,vnum
c  ord, pwr = order, exponent
      integer ord,pwr
c  indpr= index for product
      integer indpr
c  declaration for the integer function ndex(j)
      integer ndex
c  other needed integers
      integer k,ind1,ind2
c
      include 'lims.inc'
c  j = array of exponents
      integer j(6)
      save j
      data j/6*0/
c
c-----------------------------------------------------
c  Sequentially create exponent table & calculate order & rearrangements
      do 150 ind=1,topcat
        carry=j(6)
        j(6)=0
        lnzj=0
        do 100 psv=1,5
          if (j(psv).gt.0) lnzj=psv
  100   continue
        if (lnzj.gt.0) j(lnzj)=j(lnzj)-1
        j(lnzj+1)=j(lnzj+1)+1+carry
        ord=0
        do 120 psv=1,6
          pwr=j(psv)
          expon(psv,ind)=pwr
          ord=ord+pwr
  120   continue
        order(ind)=ord
c
c-------------------------------------------------------
c  Create variable list table, using exponent table
        vnum=1
        do 220 psv=1,6
          do 200 k=1,j(psv)
            vblist(vnum,ind)=psv
            vnum=vnum+1
  200     continue
  220   continue
  150 continue
c
c---------------------------------------------------------
c  Create product table, based on the idea of F. Neri and the
c   method of C. Iselin.
      do 380 psv=1,6
  380   prodex(psv,0)=psv
      do 300 ord=1,ordcat-1
        do 320 psv=1,6
          indpr=bottom(ord+1)
          do 340 ind=bottom(ord),top(ord)
  360       if (expon(psv,indpr).eq.0) then
              indpr=indpr+1
              if (indpr.le.top(ord+1)) goto 360
            endif
            prodex(psv,ind)=indpr
            indpr=indpr+1
  340     continue
  320   continue
  300 continue
c
c Create table of products
      do 401 ind1 = 1,monom1
        do 402 ind2 = 1,monom1
          do 403 psv = 1,6
            j(psv) = expon(psv,ind1)+expon(psv,ind2)
  403     continue
          iprod(ind1,ind2) = ndex(j)
  402   continue
  401 continue
      iprod(0,0) = 0
      do 500 ind1 = 1,monom1
        iprod(0,ind1) = ind1
        iprod(ind1,0) = ind1
  500 continue
c
c---------------------------------------------------------
cryne 6/21/2002 (based on ryne 8/26/2001)
c  Create sum tables needed to convert maps from one set of
c  scale variables to another.
cryne 08/26/2001
      do 600 ind=1,topcat
      nxp135(ind)=expon(1,ind)+expon(3,ind)+expon(5,ind) -1
      nxp246(ind)=expon(2,ind)+expon(4,ind)+expon(6,ind) -1
      nxp13(ind)=expon(1,ind)+expon(3,ind) -1
      nxp24(ind)=expon(2,ind)+expon(4,ind) -1
!cryne===== 12/21/2004 actually what is needed is:
! n1+n3+n6-1 for scale length [use nxp13+nxp6]
! n2+n4+n6-1 for scale momentum [use nxp24+nxp6]
! n6-n5 for scale angular freq [use nxp6-nxp5]
!!!!!!nxp5(ind)=expon(5,ind) -1
!!!!!!nxp6(ind)=expon(6,ind) -1
!cryne=====
      nxp5(ind)=expon(5,ind)
      nxp6(ind)=expon(6,ind)
  600 continue
c
      return
      end
c
***********************************************************************
c
      subroutine init
c
c     produces the arrays index1 and index2
c     used to compute vector of monomials
c     as follows:
c     do 1 k = 7,length
c     vector(k) = vector(index1(k))*vector(index2(k))
c 1   continue
c     the monomial of address k is product of
c     monomial index1(k) and index2(k)
c
c     also produces array pv(k), containing
c     the sum
c                         2
c     j(1)+j(2)*16+j(3)*16 +...
c
c     of the exponents corresponding to monomial k.
c
c
c       ind31 and ind32 are the equivalent of index1() and index2()
c        for monomials in three variables.
c       the monomials in 3 varibles are numbered in the same way the
c       monomials in 6 are, except the start from 2, not from 1.
c        1 is used for the zeroth order monomial.
c       this convention is use in ind31(), ind32(), ja3(), jv3() and ip-
c       indecn(3,**,*), returns an index that starts from 1 for a first
c       monomial so a 1 has to be added to the index returned by indecn(
c       to convert the index to the convention used in the arrays built
c      jv3() is the equivalent of jv() for third order monomials.
c      ja3(i,n3) contains the ith exponent (i=1,3) of the
c      n3 monomial in 3 varibles.
c
c Written by F. Neri, Spring 1985
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'ind.inc'
      include 'pq.inc'
      include 'ind3.inc'
      include 'ja3.inc'
c
      integer j(6),j1(6),j2(6)
      integer jp(3),jq(3)
c
      call length
      do 100 i1=0,imaxi
      do 100 i2=0,imaxi-i1
      do 100 i3=0,imaxi-i1-i2
      j(1)=i1
      j(2)=i2
      j(3)=i3
      n3 = ndexn(3,j)
      n3=n3+1
      jv3(n3)=0
      do 306 k=1,3
      jv3(n3)=jv3(n3)+j(k)*(16**(k-1))
      ja3(k,n3) = j(k)
 306  continue
      ior3=i1+i2+i3
      ii2=ior3/2
      ii1=ior3-ii2
      itot=0
      do 2 k=1,3
      j1(k)=0
      j2(k)=0
 2    continue
      do 3 i=1,3
      j1(i)=j(i)
      itot=itot+j(i)
      if(itot-ii1) 3,31,32
 3    continue
 32   continue
      j1(i)=j1(i)+ii1-itot
      j2(i)=itot-ii1
 31   continue
      if(i-3) 303,304,304
 303  continue
      do 305 ii=i+1,3
      j2(ii)=j(ii)
 305  continue
 304  continue
      n1 = ndexn(3,j1)
      n2 = ndexn(3,j2)
      ind31(n3)=n1+1
      ind32(n3)=n2+1
      do 100 i4=0,imaxi-i1-i2-i3
      do 100 i5=0,imaxi-i1-i2-i3-i4
      do 100 i6=0,imaxi-i1-i2-i3-i4-i5
      iorder=i1+i2+i3+i4+i5+i6
      if(iorder-imaxi) 101,101,100
 101  continue
      if (iorder) 100,100,107
 107  continue
      j(1)=i1
      j(2)=i2
      j(3)=i3
      j(4)=i4
      j(5)=i5
      j(6)=i6
      jp(1)=j(2)
      jp(2)=j(4)
      jp(3)=j(6)
      jq(1)=j(1)
      jq(2)=j(3)
      jq(3)=j(5)
      n = ndex(j)
      n1 = ndexn(3,jp)
      n2 = ndexn(3,jq)
      ip(n)=n1+1
      iq(n)=n2+1
      jv(n) = 0
      do 33 k = 1,6
      jv(n) = jv(n)+j(k)*(16**(k-1))
 33   continue
      ii2=iorder/2
      ii1=iorder-ii2
      itot=0
      do 55 k = 1,6
      j1(k) = 0
      j2(k) = 0
 55   continue
      do 200 i=1,6
      j1(i)=j(i)
      itot=itot+j(i)
      if(itot-ii1) 200,111,112
 200  continue
 112  continue
      j1(i)=j(i)+ii1-itot
      j2(i)=itot-ii1
 111  continue
      if(i-6) 113,114,114
 113  continue
      do 115 ii=i+1,6
      j2(ii)=j(ii)
 115  continue
 114  continue
      n1   = ndex(j1)
      n2 = ndex(j2)
      index1(n)=n1
      index2(n)=n2
 100  continue
      return
      end
c
***********************************************************************
c
      subroutine length
c  Written by F. Neri, Spring 1985
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'len.inc'
      include 'len3.inc'
      include 'ind.inc'
      dimension j(6),j3(3)
c
      do 100 k =1,6
 100  j(k) = 0
      do 103 k = 1,3
 103  j3(k)=0
      do 200 k = 1, imaxi
      j(6) = k
      j3(3)=k
      n=ndex(j)
      len(k) = ndex(j)
      len3(k) = ndexn(3,j3)
 200  continue
      return
      end
c
***********************************************************************
c
      block data matrcs
c  Written by D. Douglas, ca 1982
c      implicit none
      include 'id.inc'
      include 'symp.inc'
      data ident/1.,6*0.,1.,6*0.,1.,6*0.,1.,6*0.,1.,6*0.,1./
      data jm/0.,-1.,4*0.,1.,8*0.,-1.,4*0.,1.,8*0.,-1.,4*0.,1.,0./
      end
c
***********************************************************************
c
      block data spurex
c  Written by F. Neri, Summer 1986
      use lieaparam, only : monoms
c      implicit none
      include 'sr.inc'
c
      data (srexp(i,1),i=0,2)/0,0,0/
      data (srexp(i,2),i=0,2)/1,1,0/
      data (srexp(i,3),i=0,2)/0,0,0/
      data (srexp(i,4),i=0,2)/1,0,1/
      data (srexp(i,5),i=0,2)/0,0,0/
      data (srexp(i,6),i=0,2)/0,0,0/
      data (srexp(i,7),i=0,2)/0,0,0/
      data (srexp(i,8),i=0,2)/0,0,0/
      data (srexp(i,9),i=0,2)/0,0,0/
      data (srexp(i,10),i=0,2)/1,1,0/
      data (srexp(i,11),i=0,2)/0,0,0/
      data (srexp(i,12),i=0,2)/1,0,1/
      data (srexp(i,13),i=0,2)/0,0,0/
      data (srexp(i,14),i=0,2)/1,2,0/
      data (srexp(i,15),i=0,2)/0,0,0/
      data (srexp(i,16),i=0,2)/1,0,2/
      data (srexp(i,17),i=0,2)/0,0,0/
      data (srexp(i,18),i=0,2)/1,1,1/
      data (srexp(i,19),i=0,2)/0,0,0/
      data (srexp(i,20),i=0,2)/1,1,-1/
      data (srexp(i,21),i=0,2)/0,0,0/
      data (srexp(i,22),i=0,2)/0,0,0/
      data (srexp(i,23),i=0,2)/0,0,0/
      data (srexp(i,24),i=0,2)/0,0,0/
      data (srexp(i,25),i=0,2)/0,0,0/
      data (srexp(i,26),i=0,2)/0,0,0/
      data (srexp(i,27),i=0,2)/0,0,0/
      data (srexp(i,28),i=0,2)/0,0,0/
      data (srexp(i,29),i=0,2)/0,0,0/
      data (srexp(i,30),i=0,2)/0,0,0/
      data (srexp(i,31),i=0,2)/1,1,0/
      data (srexp(i,32),i=0,2)/0,0,0/
      data (srexp(i,33),i=0,2)/1,0,1/
      data (srexp(i,34),i=0,2)/0,0,0/
      data (srexp(i,35),i=0,2)/1,2,0/
      data (srexp(i,36),i=0,2)/0,0,0/
      data (srexp(i,37),i=0,2)/1,0,2/
      data (srexp(i,38),i=0,2)/0,0,0/
      data (srexp(i,39),i=0,2)/1,1,1/
      data (srexp(i,40),i=0,2)/0,0,0/
      data (srexp(i,41),i=0,2)/1,1,-1/
      data (srexp(i,42),i=0,2)/0,0,0/
      data (srexp(i,43),i=0,2)/1,1,0/
      data (srexp(i,44),i=0,2)/0,0,0/
      data (srexp(i,45),i=0,2)/1,0,1/
      data (srexp(i,46),i=0,2)/0,0,0/
      data (srexp(i,47),i=0,2)/1,1,0/
      data (srexp(i,48),i=0,2)/0,0,0/
      data (srexp(i,49),i=0,2)/1,0,1/
      data (srexp(i,50),i=0,2)/0,0,0/
      data (srexp(i,51),i=0,2)/1,3,0/
      data (srexp(i,52),i=0,2)/0,0,0/
      data (srexp(i,53),i=0,2)/1,0,3/
      data (srexp(i,54),i=0,2)/0,0,0/
      data (srexp(i,55),i=0,2)/1,2,1/
      data (srexp(i,56),i=0,2)/0,0,0/
      data (srexp(i,57),i=0,2)/1,1,2/
      data (srexp(i,58),i=0,2)/0,0,0/
      data (srexp(i,59),i=0,2)/1,2,-1/
      data (srexp(i,60),i=0,2)/0,0,0/
      data (srexp(i,61),i=0,2)/1,-1,2/
      data (srexp(i,62),i=0,2)/0,0,0/
      data (srexp(i,63),i=0,2)/0,0,0/
      data (srexp(i,64),i=0,2)/0,0,0/
      data (srexp(i,65),i=0,2)/0,0,0/
      data (srexp(i,66),i=0,2)/0,0,0/
      data (srexp(i,67),i=0,2)/0,0,0/
      data (srexp(i,68),i=0,2)/0,0,0/
      data (srexp(i,69),i=0,2)/0,0,0/
      data (srexp(i,70),i=0,2)/0,0,0/
      data (srexp(i,71),i=0,2)/0,0,0/
      data (srexp(i,72),i=0,2)/0,0,0/
      data (srexp(i,73),i=0,2)/0,0,0/
      data (srexp(i,74),i=0,2)/0,0,0/
      data (srexp(i,75),i=0,2)/0,0,0/
      data (srexp(i,76),i=0,2)/0,0,0/
      data (srexp(i,77),i=0,2)/0,0,0/
      data (srexp(i,78),i=0,2)/0,0,0/
      data (srexp(i,79),i=0,2)/0,0,0/
      data (srexp(i,80),i=0,2)/0,0,0/
      data (srexp(i,81),i=0,2)/0,0,0/
      data (srexp(i,82),i=0,2)/0,0,0/
      data (srexp(i,83),i=0,2)/0,0,0/
      data (srexp(i,84),i=0,2)/0,0,0/
      data (srexp(i,85),i=0,2)/0,0,0/
      data (srexp(i,86),i=0,2)/0,0,0/
      data (srexp(i,87),i=0,2)/0,0,0/
      data (srexp(i,88),i=0,2)/0,0,0/
      data (srexp(i,89),i=0,2)/0,0,0/
      data (srexp(i,90),i=0,2)/1,1,0/
      data (srexp(i,91),i=0,2)/0,0,0/
      data (srexp(i,92),i=0,2)/1,0,1/
      data (srexp(i,93),i=0,2)/0,0,0/
      data (srexp(i,94),i=0,2)/1,2,0/
      data (srexp(i,95),i=0,2)/0,0,0/
      data (srexp(i,96),i=0,2)/1,0,2/
      data (srexp(i,97),i=0,2)/0,0,0/
      data (srexp(i,98),i=0,2)/1,1,1/
      data (srexp(i,99),i=0,2)/0,0,0/
      data (srexp(i,100),i=0,2)/1,1,-1/
      data (srexp(i,101),i=0,2)/0,0,0/
      data (srexp(i,102),i=0,2)/1,1,0/
      data (srexp(i,103),i=0,2)/0,0,0/
      data (srexp(i,104),i=0,2)/1,0,1/
      data (srexp(i,105),i=0,2)/0,0,0/
      data (srexp(i,106),i=0,2)/1,1,0/
      data (srexp(i,107),i=0,2)/0,0,0/
      data (srexp(i,108),i=0,2)/1,0,1/
      data (srexp(i,109),i=0,2)/0,0,0/
      data (srexp(i,110),i=0,2)/1,3,0/
      data (srexp(i,111),i=0,2)/0,0,0/
      data (srexp(i,112),i=0,2)/1,0,3/
      data (srexp(i,113),i=0,2)/0,0,0/
      data (srexp(i,114),i=0,2)/1,2,1/
      data (srexp(i,115),i=0,2)/0,0,0/
      data (srexp(i,116),i=0,2)/1,1,2/
      data (srexp(i,117),i=0,2)/0,0,0/
      data (srexp(i,118),i=0,2)/1,2,-1/
      data (srexp(i,119),i=0,2)/0,0,0/
      data (srexp(i,120),i=0,2)/1,-1,2/
      data (srexp(i,121),i=0,2)/0,0,0/
      data (srexp(i,122),i=0,2)/1,2,0/
      data (srexp(i,123),i=0,2)/0,0,0/
      data (srexp(i,124),i=0,2)/1,0,2/
      data (srexp(i,125),i=0,2)/0,0,0/
      data (srexp(i,126),i=0,2)/1,2,0/
      data (srexp(i,127),i=0,2)/0,0,0/
      data (srexp(i,128),i=0,2)/1,0,2/
      data (srexp(i,129),i=0,2)/0,0,0/
      data (srexp(i,130),i=0,2)/1,1,1/
      data (srexp(i,131),i=0,2)/0,0,0/
      data (srexp(i,132),i=0,2)/1,1,1/
      data (srexp(i,133),i=0,2)/0,0,0/
      data (srexp(i,134),i=0,2)/1,1,-1/
      data (srexp(i,135),i=0,2)/0,0,0/
      data (srexp(i,136),i=0,2)/1,-1,1/
      data (srexp(i,137),i=0,2)/0,0,0/
      data (srexp(i,138),i=0,2)/1,4,0/
      data (srexp(i,139),i=0,2)/0,0,0/
      data (srexp(i,140),i=0,2)/1,0,4/
      data (srexp(i,141),i=0,2)/0,0,0/
      data (srexp(i,142),i=0,2)/1,3,1/
      data (srexp(i,143),i=0,2)/0,0,0/
      data (srexp(i,144),i=0,2)/1,1,3/
      data (srexp(i,145),i=0,2)/0,0,0/
      data (srexp(i,146),i=0,2)/1,3,-1/
      data (srexp(i,147),i=0,2)/0,0,0/
      data (srexp(i,148),i=0,2)/1,-1,3/
      data (srexp(i,149),i=0,2)/0,0,0/
      data (srexp(i,150),i=0,2)/1,2,2/
      data (srexp(i,151),i=0,2)/0,0,0/
      data (srexp(i,152),i=0,2)/1,2,-2/
      data (srexp(i,153),i=0,2)/0,0,0/
      data (srexp(i,154),i=0,2)/0,0,0/
      data (srexp(i,155),i=0,2)/0,0,0/
      data (srexp(i,156),i=0,2)/0,0,0/
      data (srexp(i,157),i=0,2)/0,0,0/
      data (srexp(i,158),i=0,2)/0,0,0/
      data (srexp(i,159),i=0,2)/0,0,0/
      data (srexp(i,160),i=0,2)/0,0,0/
      data (srexp(i,161),i=0,2)/0,0,0/
      data (srexp(i,162),i=0,2)/0,0,0/
      data (srexp(i,163),i=0,2)/0,0,0/
      data (srexp(i,164),i=0,2)/0,0,0/
      data (srexp(i,165),i=0,2)/0,0,0/
      data (srexp(i,166),i=0,2)/0,0,0/
      data (srexp(i,167),i=0,2)/0,0,0/
      data (srexp(i,168),i=0,2)/0,0,0/
      data (srexp(i,169),i=0,2)/0,0,0/
      data (srexp(i,170),i=0,2)/0,0,0/
      data (srexp(i,171),i=0,2)/0,0,0/
      data (srexp(i,172),i=0,2)/0,0,0/
      data (srexp(i,173),i=0,2)/0,0,0/
      data (srexp(i,174),i=0,2)/0,0,0/
      data (srexp(i,175),i=0,2)/0,0,0/
      data (srexp(i,176),i=0,2)/0,0,0/
      data (srexp(i,177),i=0,2)/0,0,0/
      data (srexp(i,178),i=0,2)/0,0,0/
      data (srexp(i,179),i=0,2)/0,0,0/
      data (srexp(i,180),i=0,2)/0,0,0/
      data (srexp(i,181),i=0,2)/0,0,0/
      data (srexp(i,182),i=0,2)/0,0,0/
      data (srexp(i,183),i=0,2)/0,0,0/
      data (srexp(i,184),i=0,2)/0,0,0/
      data (srexp(i,185),i=0,2)/0,0,0/
      data (srexp(i,186),i=0,2)/0,0,0/
      data (srexp(i,187),i=0,2)/0,0,0/
      data (srexp(i,188),i=0,2)/0,0,0/
      data (srexp(i,189),i=0,2)/0,0,0/
      data (srexp(i,190),i=0,2)/0,0,0/
      data (srexp(i,191),i=0,2)/0,0,0/
      data (srexp(i,192),i=0,2)/0,0,0/
      data (srexp(i,193),i=0,2)/0,0,0/
      data (srexp(i,194),i=0,2)/0,0,0/
      data (srexp(i,195),i=0,2)/0,0,0/
      data (srexp(i,196),i=0,2)/0,0,0/
      data (srexp(i,197),i=0,2)/0,0,0/
      data (srexp(i,198),i=0,2)/0,0,0/
      data (srexp(i,199),i=0,2)/0,0,0/
      data (srexp(i,200),i=0,2)/0,0,0/
      data (srexp(i,201),i=0,2)/0,0,0/
      data (srexp(i,202),i=0,2)/0,0,0/
      data (srexp(i,203),i=0,2)/0,0,0/
      data (srexp(i,204),i=0,2)/0,0,0/
      data (srexp(i,205),i=0,2)/0,0,0/
      data (srexp(i,206),i=0,2)/0,0,0/
      data (srexp(i,207),i=0,2)/0,0,0/
      data (srexp(i,208),i=0,2)/0,0,0/
      data (srexp(i,209),i=0,2)/0,0,0/
      end
c
***************************************************************
c
      block data dpurex
c  Written by F. Neri, Summer 1986
      use lieaparam, only : monoms
c      implicit none
      include 'dr.inc'
c
      data (drexp(i,1),i=0,3)/1,1,0,0/
      data (drexp(i,2),i=0,3)/0,0,0,0/
      data (drexp(i,3),i=0,3)/1,0,1,0/
      data (drexp(i,4),i=0,3)/0,0,0,0/
      data (drexp(i,5),i=0,3)/1,0,0,1/
      data (drexp(i,6),i=0,3)/0,0,0,0/
      data (drexp(i,7),i=0,3)/0,0,0,0/
      data (drexp(i,8),i=0,3)/0,0,0,0/
      data (drexp(i,9),i=0,3)/0,0,0,0/
      data (drexp(i,10),i=0,3)/1,0,0,2/
      data (drexp(i,11),i=0,3)/0,0,0,0/
      data (drexp(i,12),i=0,3)/1,1,0,1/
      data (drexp(i,13),i=0,3)/0,0,0,0/
      data (drexp(i,14),i=0,3)/1,1,0,-1/
      data (drexp(i,15),i=0,3)/0,0,0,0/
      data (drexp(i,16),i=0,3)/1,0,1,1/
      data (drexp(i,17),i=0,3)/0,0,0,0/
      data (drexp(i,18),i=0,3)/1,0,1,-1/
      data (drexp(i,19),i=0,3)/0,0,0,0/
      data (drexp(i,20),i=0,3)/1,2,0,0/
      data (drexp(i,21),i=0,3)/0,0,0,0/
      data (drexp(i,22),i=0,3)/1,0,2,0/
      data (drexp(i,23),i=0,3)/0,0,0,0/
      data (drexp(i,24),i=0,3)/1,1,1,0/
      data (drexp(i,25),i=0,3)/0,0,0,0/
      data (drexp(i,26),i=0,3)/1,1,-1,0/
      data (drexp(i,27),i=0,3)/0,0,0,0/
      data (drexp(i,28),i=0,3)/1,0,0,1/
      data (drexp(i,29),i=0,3)/0,0,0,0/
      data (drexp(i,30),i=0,3)/1,0,0,1/
      data (drexp(i,31),i=0,3)/0,0,0,0/
      data (drexp(i,32),i=0,3)/1,0,0,3/
      data (drexp(i,33),i=0,3)/0,0,0,0/
      data (drexp(i,34),i=0,3)/1,0,0,1/
      data (drexp(i,35),i=0,3)/0,0,0,0/
      data (drexp(i,36),i=0,3)/1,1,0,0/
      data (drexp(i,37),i=0,3)/0,0,0,0/
      data (drexp(i,38),i=0,3)/1,0,1,0/
      data (drexp(i,39),i=0,3)/0,0,0,0/
      data (drexp(i,40),i=0,3)/1,1,0,2/
      data (drexp(i,41),i=0,3)/0,0,0,0/
      data (drexp(i,42),i=0,3)/1,1,0,-2/
      data (drexp(i,43),i=0,3)/0,0,0,0/
      data (drexp(i,44),i=0,3)/1,0,1,2/
      data (drexp(i,45),i=0,3)/0,0,0,0/
      data (drexp(i,46),i=0,3)/1,0,1,-2/
      data (drexp(i,47),i=0,3)/0,0,0,0/
      data (drexp(i,48),i=0,3)/1,2,0,1/
      data (drexp(i,49),i=0,3)/0,0,0,0/
      data (drexp(i,50),i=0,3)/1,2,0,-1/
      data (drexp(i,51),i=0,3)/0,0,0,0/
      data (drexp(i,52),i=0,3)/1,0,2,1/
      data (drexp(i,53),i=0,3)/0,0,0,0/
      data (drexp(i,54),i=0,3)/1,0,2,-1/
      data (drexp(i,55),i=0,3)/0,0,0,0/
      data (drexp(i,56),i=0,3)/1,1,1,1/
      data (drexp(i,57),i=0,3)/0,0,0,0/
      data (drexp(i,58),i=0,3)/1,1,1,-1/
      data (drexp(i,59),i=0,3)/0,0,0,0/
      data (drexp(i,60),i=0,3)/1,1,-1,1/
      data (drexp(i,61),i=0,3)/0,0,0,0/
      data (drexp(i,62),i=0,3)/1,1,-1,-1/
      data (drexp(i,63),i=0,3)/0,0,0,0/
      data (drexp(i,64),i=0,3)/1,1,0,0/
      data (drexp(i,65),i=0,3)/0,0,0,0/
      data (drexp(i,66),i=0,3)/1,0,1,0/
      data (drexp(i,67),i=0,3)/0,0,0,0/
      data (drexp(i,68),i=0,3)/1,1,0,0/
      data (drexp(i,69),i=0,3)/0,0,0,0/
      data (drexp(i,70),i=0,3)/1,0,1,0/
      data (drexp(i,71),i=0,3)/0,0,0,0/
      data (drexp(i,72),i=0,3)/1,3,0,0/
      data (drexp(i,73),i=0,3)/0,0,0,0/
      data (drexp(i,74),i=0,3)/1,0,3,0/
      data (drexp(i,75),i=0,3)/0,0,0,0/
      data (drexp(i,76),i=0,3)/1,2,1,0/
      data (drexp(i,77),i=0,3)/0,0,0,0/
      data (drexp(i,78),i=0,3)/1,1,2,0/
      data (drexp(i,79),i=0,3)/0,0,0,0/
      data (drexp(i,80),i=0,3)/1,2,-1,0/
      data (drexp(i,81),i=0,3)/0,0,0,0/
      data (drexp(i,82),i=0,3)/1,1,-2,0/
      data (drexp(i,83),i=0,3)/0,0,0,0/
      data (drexp(i,84),i=0,3)/0,0,0,0/
      data (drexp(i,85),i=0,3)/0,0,0,0/
      data (drexp(i,86),i=0,3)/0,0,0,0/
      data (drexp(i,87),i=0,3)/0,0,0,0/
      data (drexp(i,88),i=0,3)/0,0,0,0/
      data (drexp(i,89),i=0,3)/0,0,0,0/
      data (drexp(i,90),i=0,3)/1,0,0,2/
      data (drexp(i,91),i=0,3)/0,0,0,0/
      data (drexp(i,92),i=0,3)/1,0,0,2/
      data (drexp(i,93),i=0,3)/0,0,0,0/
      data (drexp(i,94),i=0,3)/1,0,0,4/
      data (drexp(i,95),i=0,3)/0,0,0,0/
      data (drexp(i,96),i=0,3)/1,0,0,2/
      data (drexp(i,97),i=0,3)/0,0,0,0/
      data (drexp(i,98),i=0,3)/1,1,0,3/
      data (drexp(i,99),i=0,3)/0,0,0,0/
      data (drexp(i,100),i=0,3)/1,1,0,-3/
      data (drexp(i,101),i=0,3)/0,0,0,0/
      data (drexp(i,102),i=0,3)/1,0,1,3/
      data (drexp(i,103),i=0,3)/0,0,0,0/
      data (drexp(i,104),i=0,3)/1,0,1,-3/
      data (drexp(i,105),i=0,3)/0,0,0,0/
      data (drexp(i,106),i=0,3)/1,1,0,1/
      data (drexp(i,107),i=0,3)/0,0,0,0/
      data (drexp(i,108),i=0,3)/1,1,0,-1/
      data (drexp(i,109),i=0,3)/0,0,0,0/
      data (drexp(i,110),i=0,3)/1,0,1,1/
      data (drexp(i,111),i=0,3)/0,0,0,0/
      data (drexp(i,112),i=0,3)/1,0,1,-1/
      data (drexp(i,113),i=0,3)/0,0,0,0/
      data (drexp(i,114),i=0,3)/1,2,0,0/
      data (drexp(i,115),i=0,3)/0,0,0,0/
      data (drexp(i,116),i=0,3)/1,0,2,0/
      data (drexp(i,117),i=0,3)/0,0,0,0/
      data (drexp(i,118),i=0,3)/1,2,0,2/
      data (drexp(i,119),i=0,3)/0,0,0,0/
      data (drexp(i,120),i=0,3)/1,2,0,-2/
      data (drexp(i,121),i=0,3)/0,0,0,0/
      data (drexp(i,122),i=0,3)/1,0,2,2/
      data (drexp(i,123),i=0,3)/0,0,0,0/
      data (drexp(i,124),i=0,3)/1,0,2,-2/
      data (drexp(i,125),i=0,3)/0,0,0,0/
      data (drexp(i,126),i=0,3)/1,1,1,0/
      data (drexp(i,127),i=0,3)/0,0,0,0/
      data (drexp(i,128),i=0,3)/1,1,-1,0/
      data (drexp(i,129),i=0,3)/0,0,0,0/
      data (drexp(i,130),i=0,3)/1,1,1,2/
      data (drexp(i,131),i=0,3)/0,0,0,0/
      data (drexp(i,132),i=0,3)/1,1,1,-2/
      data (drexp(i,133),i=0,3)/0,0,0,0/
      data (drexp(i,134),i=0,3)/1,1,-1,2/
      data (drexp(i,135),i=0,3)/0,0,0,0/
      data (drexp(i,136),i=0,3)/1,1,-1,-2/
      data (drexp(i,137),i=0,3)/0,0,0,0/
      data (drexp(i,138),i=0,3)/1,1,0,1/
      data (drexp(i,139),i=0,3)/0,0,0,0/
      data (drexp(i,140),i=0,3)/1,1,0,-1/
      data (drexp(i,141),i=0,3)/0,0,0,0/
      data (drexp(i,142),i=0,3)/1,0,1,1/
      data (drexp(i,143),i=0,3)/0,0,0,0/
      data (drexp(i,144),i=0,3)/1,0,1,-1/
      data (drexp(i,145),i=0,3)/0,0,0,0/
      data (drexp(i,146),i=0,3)/1,1,0,1/
      data (drexp(i,147),i=0,3)/0,0,0,0/
      data (drexp(i,148),i=0,3)/1,1,0,-1/
      data (drexp(i,149),i=0,3)/0,0,0,0/
      data (drexp(i,150),i=0,3)/1,0,1,1/
      data (drexp(i,151),i=0,3)/0,0,0,0/
      data (drexp(i,152),i=0,3)/1,0,1,-1/
      data (drexp(i,153),i=0,3)/0,0,0,0/
      data (drexp(i,154),i=0,3)/1,3,0,1/
      data (drexp(i,155),i=0,3)/0,0,0,0/
      data (drexp(i,156),i=0,3)/1,3,0,-1/
      data (drexp(i,157),i=0,3)/0,0,0,0/
      data (drexp(i,158),i=0,3)/1,0,3,1/
      data (drexp(i,159),i=0,3)/0,0,0,0/
      data (drexp(i,160),i=0,3)/1,0,3,-1/
      data (drexp(i,161),i=0,3)/0,0,0,0/
      data (drexp(i,162),i=0,3)/1,2,1,1/
      data (drexp(i,163),i=0,3)/0,0,0,0/
      data (drexp(i,164),i=0,3)/1,2,1,-1/
      data (drexp(i,165),i=0,3)/0,0,0,0/
      data (drexp(i,166),i=0,3)/1,1,2,1/
      data (drexp(i,167),i=0,3)/0,0,0,0/
      data (drexp(i,168),i=0,3)/1,1,2,-1/
      data (drexp(i,169),i=0,3)/0,0,0,0/
      data (drexp(i,170),i=0,3)/1,2,-1,1/
      data (drexp(i,171),i=0,3)/0,0,0,0/
      data (drexp(i,172),i=0,3)/1,2,-1,-1/
      data (drexp(i,173),i=0,3)/0,0,0,0/
      data (drexp(i,174),i=0,3)/1,1,-2,-1/
      data (drexp(i,175),i=0,3)/0,0,0,0/
      data (drexp(i,176),i=0,3)/1,1,-2,1/
      data (drexp(i,177),i=0,3)/0,0,0,0/
      data (drexp(i,178),i=0,3)/1,2,0,0/
      data (drexp(i,179),i=0,3)/0,0,0,0/
      data (drexp(i,180),i=0,3)/1,0,2,0/
      data (drexp(i,181),i=0,3)/0,0,0,0/
      data (drexp(i,182),i=0,3)/1,2,0,0/
      data (drexp(i,183),i=0,3)/0,0,0,0/
      data (drexp(i,184),i=0,3)/1,0,2,0/
      data (drexp(i,185),i=0,3)/0,0,0,0/
      data (drexp(i,186),i=0,3)/1,1,1,0/
      data (drexp(i,187),i=0,3)/0,0,0,0/
      data (drexp(i,188),i=0,3)/1,1,1,0/
      data (drexp(i,189),i=0,3)/0,0,0,0/
      data (drexp(i,190),i=0,3)/1,1,-1,0/
      data (drexp(i,191),i=0,3)/0,0,0,0/
      data (drexp(i,192),i=0,3)/1,1,-1,0/
      data (drexp(i,193),i=0,3)/0,0,0,0/
      data (drexp(i,194),i=0,3)/1,4,0,0/
      data (drexp(i,195),i=0,3)/0,0,0,0/
      data (drexp(i,196),i=0,3)/1,0,4,0/
      data (drexp(i,197),i=0,3)/0,0,0,0/
      data (drexp(i,198),i=0,3)/1,3,1,0/
      data (drexp(i,199),i=0,3)/0,0,0,0/
      data (drexp(i,200),i=0,3)/1,1,3,0/
      data (drexp(i,201),i=0,3)/0,0,0,0/
      data (drexp(i,202),i=0,3)/1,3,-1,0/
      data (drexp(i,203),i=0,3)/0,0,0,0/
      data (drexp(i,204),i=0,3)/1,1,-3,0/
      data (drexp(i,205),i=0,3)/0,0,0,0/
      data (drexp(i,206),i=0,3)/1,2,2,0/
      data (drexp(i,207),i=0,3)/0,0,0,0/
      data (drexp(i,208),i=0,3)/1,2,-2,0/
      data (drexp(i,209),i=0,3)/0,0,0,0/
      end
c
***********************************************************************
c
      block data srlabl
      use lieaparam, only : monoms
c      implicit none
      include 'srl.inc'
c
c Assignment of designation labels for static resonance basis:
c  Written by F. Neri, Summer 1986
c
      data(sln(i),i=1,25)/
     &'R00001','R10000','I10000','R00100','I00100',
     &'000010','R11000','R00110','R00002','R10001',
     &'I10001','R00101','I00101','R20000','I20000',
     &'R00200','I00200','R10100','I10100','R10010',
     &'I10010','100010','010010','001010','000110'/
      data(sln(i),i=26,50)/
     &'000020','000011','R11001','R00111','R00003',
     &'R10002','I10002','R00102','I00102','R20001',
     &'I20001','R00201','I00201','R10101','I10101',
     &'R10011','I10011','R21000','I21000','R00210',
     &'I00210','R10110','I10110','R11100','I11100'/
      data(sln(i),i=51,75)/
     &'R30000','I30000','R00300','I00300','R20100',
     &'I20100','R10200','I10200','R20010','I20010',
     &'R01200','I01200','200010','110010','101010',
     &'100110','100020','100011','020010','011010',
     &'010110','010020','010011','002010','001110'/
      data(sln(i),i=76,100)/
     &'001020','001011','000210','000120','000111',
     &'000030','000021','000012','R11002','R00112',
     &'R00004','R22000','R00220','R11110','R10003',
     &'I10003','R00103','I00103','R20002','I20002',
     &'R00202','I00202','R10102','I10102','R10012'/
      data(sln(i),i=101,125)/
     &'I10012','R21001','I21001','R00211','I00211',
     &'R10111','I10111','R11101','I11101','R30001',
     &'I30001','R00301','I00301','R20101','I20101',
     &'R10201','I10201','R20011','I20011','R01201',
     &'I01201','R31000','I31000','R00310','I00310'/
      data(sln(i),i=126,150)/
     &'R20110','I20110','R11200','I11200','R21100',
     &'I21100','R10210','I10210','R21010','I21010',
     &'R01210','I01210','R40000','I40000','R00400',
     &'I00400','R30100','I30100','R10300','I10300',
     &'R30010','I30010','R01300','I01300','R20200'/
      data(sln(i),i=151,175)/
     &'I20200','R20020','I20020','300010','210010',
     &'201010','200110','200020','200011','120010',
     &'111010','110110','110020','110011','102010',
     &'101110','101020','101011','100210','100120',
     &'100111','100030','100021','100012','030010'/
      data(sln(i),i=176,200)/
     &'021010','020110','020020','020011','012010',
     &'011110','011020','011011','010210','010120',
     &'010111','010030','010021','010012','003010',
     &'002110','002020','002011','001210','001120',
     &'001111','001030','001021','001012','000310'/
      data(sln(i),i=201,209)/
     &'000220','000211','000130','000121','000112',
     &'000040','000031','000022','000013'/
c
      end
c
************************************************************************
c
      block data drlabl
      use lieaparam, only : monoms
c      implicit none
      include 'drl.inc'
c
c Assignment of designation labels for dynamic resonance basis:
c  Written by F. Neri, Summer 1986
c
      data(dln(i),i=1,25)/
     &'R100000','I100000','R001000','I001000','R000010',
     &'I000010','R110000','R001100','R000011','R000020',
     &'I000020','R100010','I100010','R100001','I100001',
     &'R001010','I001010','R001001','I001001','R200000',
     &'I200000','R002000','I002000','R101000','I101000'/
      data(dln(i),i=26,50)/
     &'R100100','I100100','R110010','I110010','R001110',
     &'I001110','R000030','I000030','R000021','I000021',
     &'R100011','I100011','R001011','I001011','R100020',
     &'I100020','R100002','I100002','R001020','I001020',
     &'R001002','I001002','R200010','I200010','R200001'/
      data(dln(i),i=51,75)/
     &'I200001','R002010','I002010','R002001','I002001',
     &'R101010','I101010','R101001','I101001','R100110',
     &'I100110','R100101','I100101','R210000','I210000',
     &'R002100','I002100','R101100','I101100','R111000',
     &'I111000','R300000','I300000','R003000','I003000'/
      data(dln(i),i=76,100)/
     &'R201000','I201000','R102000','I102000','R200100',
     &'I200100','R100200','I100200','R220000','R002200',
     &'R000022','R111100','R110011','R001111','R110020',
     &'I110020','R001120','I001120','R000040','I000040',
     &'R000031','I000031','R100030','I100030','R100003'/
      data(dln(i),i=101,125)/
     &'I100003','R001030','I001030','R001003','I001003',
     &'R100021','I100021','R100012','I100012','R001021',
     &'I001021','R001012','I001012','R200011','I200011',
     &'R002011','I002011','R200020','I200020','R200002',
     &'I200002','R002020','I002020','R002002','I002002'/
      data(dln(i),i=126,150)/
     &'R101011','I101011','R100111','I100111','R101020',
     &'I101020','R101002','I101002','R100120','I100120',
     &'R100102','I100102','R210010','I210010','R210001',
     &'I210001','R002110','I002110','R002101','I002101',
     &'R101110','I101110','R101101','I101101','R111010'/
      data(dln(i),i=151,175)/
     &'I111010','R111001','I111001','R300010','I300010',
     &'R300001','I300001','R003010','I003010','R003001',
     &'I003001','R201010','I201010','R201001','I201001',
     &'R102010','I102010','R102001','I102001','R200110',
     &'I200110','R200101','I200101','R100201','I100201'/
      data(dln(i),i=176,200)/
     &'R100210','I100210','R310000','I310000','R003100',
     &'I003100','R201100','I201100','R112000','I112000',
     &'R211000','I211000','R102100','I102100','R210100',
     &'I210100','R101200','I101200','R400000','I400000',
     &'R004000','I004000','R301000','I301000','R103000'/
      data(dln(i),i=201,209)/
     &'I103000','R300100','I300100','R100300','I100300',
     &'R202000','I202000','R200200','I200200'/
c
      end
c end of file
