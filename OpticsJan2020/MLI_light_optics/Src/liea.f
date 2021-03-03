************************************************************************
* header:                 LIE ALGEBRAIC                                *
*  Lie algebraic manipulations: Poisson brackets, etc.                 *
************************************************************************
c
      subroutine brkts(h)
c
c     computes poisson brackets of total
c     lattice generator h with each
c     dynamical variable  z(i)
c     pbh(j,i=1,6) contains
c     exp(:h3:) exp(:h4:) exp(:h5:) exp(:h6:) z(i)
c     truncated to sixth order ( for marylie 5.0 )
c
c
      use rays
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'pbkh.inc'
      include 'ind.inc'
      include 'lims.inc'
cryne 7/23/2002      implicit double precision (a-h,o-z)
c     double precision tmh(6,6)
      dimension h(923),pb1(923),pb2(923),pb3(923),pb4(923),pb5(923)
cryne 7/23/2002      common /pbkh/ pbh(923,12)
cryne 7/23/2002      common/ind/imaxi,jv(923),index1(923),index2(923)
cryne 7/23/2002      integer bottom(0:12),top(0:12)
cryne 7/23/2002      common/lims/bottom,top
c
c     if(idproc.eq.0)write(6,*)'inside brkts'
c
c     initialize arrays
c
      do 10 i=1,923
      pb1(i)=0.0d0
      pb2(i)=0.0d0
      pb3(i)=0.0d0
      do 20 j=1,12
      pbh(i,j)=0.0d0
   20 continue
   10 continue
c
c     if(idproc.eq.0)write(6,*)'done with initialization'
      do i=1,27
      if(h(i).ne.h(i))write(6,*)'warning: element undefined:',i
      enddo
c
c     compute poisson brackets and store in pbh
c
      do 100 iz = 1,6
c     if(idproc.eq.0)write(6,*)'do 100;  iz=',iz
        do 110 i=1,923
          pb1(i) = 0.d0
  110   continue
        pb1(iz)=1.d0
        do 200 ideg = imaxi,3,-1
          call exphf(h,ideg,pb1,imaxi,pb2)
          do 220 i=1,top(imaxi)
            pb1(i) = pb2(i)
  220     continue
  200   continue
        do 300 i = 1,top(imaxi)
          pbh(i,iz) = pb1(i)
  300   continue
  100 continue
c      if(idproc.eq.0)write(6,*)'done with do 100'
c      debug = 0
c      if ( debug .eq. 1) then
c      do 7 i=1,6
c       do 70 j=1,6
c  70    if(tmh(i,j) .ne. 0.d0 ) write(36,*) i,j,tmh(i,j)
c       do 7 no = 2,6
c        call xform5(pbh(1,i),no,tmh,pb1)
c        do 7 j=bottom(no),top(no)
c          if(pb1(j).ne.0.d0 ) write(36,*) i,j,pb1(j)
c   7  continue
c      endif
c
cryne 8/6/2002
c     pbh6=0.d0
c     do i=7,top(imaxi)
c     do j=1,6
c       pbh6(i,j)=pbh(i,j)
c     enddo
c     enddo
c     if(idproc.eq.0)write(6,*)'here I am before alloc check'
      if(allocated(pbh6t))then
c       if(idproc.eq.0)write(6,*)'pbh6t is already allocated'
        pbh6t=0.d0
        do i=7,top(imaxi)
          do j=1,6
            pbh6t(j,i)=pbh(i,j)
          enddo
        enddo
      endif
c     if(idproc.eq.0)write(6,*)'done with final do loop'
c     if(idproc.eq.0) then
c       write(6,*) ' === liea::brkts() ==='
c       write(6,'(a)') ' pbh6t(:) ='
c       do i=1,top(imaxi)
c         write(6,123) i,(pbh6t(j,i),j=1,6)
c123      format(i4,6(1x,1pe16.9))
c       enddo
c       write(6,*) ' leaving liea::brkts()'
c       write(6,*) ' ====================='
c     end if
      return
      end
c
************************************************************************
c
      subroutine clear(h,mh)
c Clears to zero the polynomials h and the matrix mh.
c  Written by Liam Healy, June 12, 1984.
      use lieaparam, only : monoms
      double precision mh(6,6),h(monoms)
c  Clear out polynomial coefficients:
      do 120 i=1,monoms
  120 h(i)=0.
c  Clear out matrix:
      do 100 i=1,6
      do 100 j=1,6
  100 mh(i,j)=0.
      return
      end
***********************************************************************
c
      subroutine concat(fa,fm,ga,gm,ha,hm)
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  variables
c

      dimension fa(923),fm(6,6)
      dimension ga(923),gm(6,6)
      dimension ha(923),hm(6,6)
      dimension ha1(923),hm1(6,6)
c
c      write(6,*) 'first factor in concat is'
c      call pcmap(1,1,0,0,fa,fm)
c      write(6,*) 'second factor in concat is'
c      call pcmap(1,1,0,0,ga,gm)
c
      maxcat = 6
      call drcat(fa,fm,ga,gm,ha1,hm1)
c
c      write(6,*) 'result of concat is'
c      call pcmap(1,1,0,0,ha1,hm1)
c
      call mapmap(ha1,hm1,ha,hm)
c
      return
      end
c
********************************************************************
c
      subroutine cpadd(fa,ga,ha)
c this is a subroutine for computing the sum of two polynomials
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ha(monoms)
c
c perform calculation
      do 10 j=1,monoms
      ha(j)=fa(j)+ga(j)
   10 continue
c
      return
      end
c
************************************************************************
c
      subroutine cpdnf(pow,fa,fm,ga,gm)
c this subroutine computes the power of a dynamic normal form map:
c mapg=(mapf)**pow
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms)
      dimension fm(6,6),gm(6,6)
      dimension ta(monoms),t1a(monoms)
c
c clear output array
      call clear(t1a,gm)
c
c begin calculation
c
c compute phase advances
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
      cwt=fm(5,5)
      swt=fm(5,6)
      wt=atan2(swt,cwt)
c compute new quantities
      pwx=pow*wx
      cpwx=cos(pwx)
      spwx=sin(pwx)
      pwy=pow*wy
      cpwy=cos(pwy)
      spwy=sin(pwy)
      pwt=pow*wt
      cpwt=cos(pwt)
      spwt=sin(pwt)
c compute new matrix
      gm(1,1)=cpwx
      gm(1,2)=spwx
      gm(2,1)=-spwx
      gm(2,2)=cpwx
      gm(3,3)=cpwy
      gm(3,4)=spwy
      gm(4,3)=-spwy
      gm(4,4)=cpwy
      gm(5,5)=cpwt
      gm(5,6)=spwt
      gm(6,5)=-spwt
      gm(6,6)=cpwt
c compute polynomials of old normal form map in dynamic resonance basis
      call ctodr(fa,ta)
c pick out those terms which are allowed to be nonzero
      t1a(84)=ta(84)
      t1a(85)=ta(85)
      t1a(86)=ta(86)
      t1a(87)=ta(87)
      t1a(88)=ta(88)
      t1a(89)=ta(89)
c compute polynomials of new map in dynamic resonance basis
      call csmul(pow,t1a,t1a)
c transform back to cartesian basis
      call drtoc(t1a,ga)
c
      return
      end
c
********************************************************************
c
      subroutine cppb(fa,ga,ha)
c This subroutine computes the poisson brackets two polynomials
c It gives the result ha=[fa,ga]
c The aray ha is cleared upon entry
c written 5/22/02 AJD
c revised 6/14/02 AJD
c
      use lieaparam, only : monoms
      include 'impli.inc'
cryne include 'param.inc'
c
c     calling arrays
c
      dimension fa(monoms)
      dimension ga(monoms)
      dimension ha(monoms)
c
c     working arrays
      dimension ta(monoms)
      dimension tm(6,6)
c
c clear ha
      call clear(ha,tm)
c
c compute Poisson bracket
      do if=1,6 !4 changed to 6
      do ig=1,6 !4 changed to 6
      iord=if+ig-2
      if((iord .ge. 0) .and. (iord .le. 6)) then !4 changed to 6
c clear the result array ta since pbkt does not do so completely
      call clear(ta,tm)
      call pbkt(fa,if,ga,ig,ta)
      call cpadd(ha,ta,ha)
      endif
      enddo
      enddo
      return
      end
c
********************************************************************
c
      subroutine cpmul(fa,ga,ha)
c Subroutine for computing the product of polynomials
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ha(monoms)
c
c perform calculation
      do 10 nf=1,3
      do 20 ng=1,3
      if ((nf+ng).le.4) then
      call pprod(fa,nf,ga,ng,ha)
      endif
   20 continue
   10 continue
c
      return
      end
c
************************************************************************
c
      subroutine cpsnf(pow,fa,fm,ga,gm)
c this subroutine computes the power of a static normal form map:
c mapg=(mapf)**pow
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms)
      dimension fm(6,6),gm(6,6)
      dimension ta(monoms),t1a(monoms)
c
c clear output array
      call clear(t1a,gm)
c
c begin calculation
c
c compute phase advances
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
c compute momentum compaction
      wt=fm(5,6)
c compute new quantities
      pwx=pow*wx
      cpwx=cos(pwx)
      spwx=sin(pwx)
      pwy=pow*wy
      cpwy=cos(pwy)
      spwy=sin(pwy)
      pwt=pow*wt
c compute new matrix
      gm(1,1)=cpwx
      gm(1,2)=spwx
      gm(2,1)=-spwx
      gm(2,2)=cpwx
      gm(3,3)=cpwy
      gm(3,4)=spwy
      gm(4,3)=-spwy
      gm(4,4)=cpwy
      gm(5,5)=1.d0
      gm(5,6)=pwt
      gm(6,6)=1.d0
c compute polynomials of old normal form map in static resonance basis
      call ctosr(fa,ta)
c pick out those terms which are allowed to be nonzero
      t1a(28)=ta(28)
      t1a(29)=ta(29)
      t1a(30)=ta(30)
      t1a(84)=ta(84)
      t1a(85)=ta(85)
      t1a(86)=ta(86)
      t1a(87)=ta(87)
      t1a(88)=ta(88)
      t1a(89)=ta(89)
c compute polynomials of new map in static resonance basis
      call csmul(pow,t1a,t1a)
c transform back to cartesian basis
      call srtoc(t1a,ga)
c
      return
      end
c
********************************************************************
c
      subroutine csmul(scalar,fa,ga)
c this subroutine computes the scalar multiple of a polynomial
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms)
c
c perform calculation
      do 10 j=1,monoms
      ga(j)=scalar*fa(j)
   10 continue
c
      return
      end
c
***********************************************************************
c
      subroutine drcat(f,mf,g,mg,h,mh)
c     concatenates a map with linear piece
c     represented by a matrix mf and nonlinearities
c     represented by exp:f:, with a map whose
c     linear piece has matrix representation mg
c     and whose nonlinearities are represented by exp:g:
c
c     the result is a map with linearities possessing
c     a matrix representation mh=mg*mf
c     and with nonlinearities generated by exp:h:
c
c     the "f" map is assumed to occur first in the beamline
c     the "g" map is encountered second
c
      implicit double precision (a-h,o-z)
      double precision mf,mg,mh,m
      dimension f(923), g(923), h(923)
      dimension mf(6,6),mg(6,6),mh(6,6),m(6,6)
      dimension f3(923),f4(923),f5(923),f6(923)
      dimension t3(923),t4(923),t5(923),t6(923)
      include 'symp.inc'
c
c      write(6,*) 'first factor in mycat'
c      call pcmap(1,1,0,0,f,mf)
c      write(6,*) 'second factor in mycat'
c      call pcmap(1,1,0,0,g,mg)
c
      call ident(h,mh)
      do 666 i1=1,923
        t3(i1) = 0.0d0
  666 continue
c
c     compute mh=mg*mf 
c
      call mmult(mg,mf,mh)
c
c compute m = mginverse
c
      call matmat(mg,m)
      call minv(m)
c
c  compute transformed arrays
c
       call xform5(f,3,m,f3)
       call xform5(f,4,m,f4)
       call xform5(f,5,m,f5)
c
c      write(6,*) 'contents of f5 and m'
c      call pcmap(1,1,0,0,f5,m)
c
       call xform5(f,6,m,f6)
c    third order terms
       call pmadd(g,3,1.d0,h)
       call pmadd(f3,3,1.d0,h)
c   fourth order terms
       call pbkt1(f3,3,g,3,t4)
       call pmadd(g,4,1.d0,h)
       call pmadd(f4,4,1.d0,h)
       call pmadd(t4,4,.5d0,h)
c    fifth order terms
        call pmadd(g,5,1.d0,h)
        call pmadd(f5,5,1.d0,h)
        call pmadd(g,3,1.d0,t3)
        call pmadd(f3,3,1.d0/2.d0,t3)
        call pbkt1(t3,3,t4,4,t5)
        call pmadd(t5,5,-1.d0/3.d0,h)
        call pbkt1(g,3,f4,4,t5)
        call pmadd(t5,5,-1.d0,h)
c      write(6,*) 'result of drcat at end of fifth order'
c      call pcmap(1,1,0,0,h,mh)

c    sixth order terms
        call pbkt1(g,4,f4,4,t6)
        call pmadd(t6,6,-.5d0,h)
        call pmadd(g,6,1.d0,h)
        call pmadd(f6,6,1.d0,h)
        call pbkt1(g,3,f5,5,t6)
        call pmadd(t6,6,-1.d0,h)
        call pbkt1(f3,3,g,3,t4)
        call pbkt1(g,3,t4,4,t5)
        call pbkt1(f3,3,t5,5,t6)
        call pmadd(t6,6,(-1.d0/24.d0),h)
        call pbkt1(h,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(1.d0/12.d0),h)
        call pbkt1(g,3,f4,4,t5)
        call pbkt1(g,3,t5,5,t6)
        call pmadd(t6,6,(.5d0),h)
        call pbkt1(f4,4,t4,4,t6)
        call pmadd(t6,6,-.25d0,h)
        call pbkt1(g,4,t4,4,t6)
        call pmadd(t6,6,-.25d0,h)
        call pbkt1(g,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(1.d0/24.d0),h)
        call pbkt1(f3,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(-1.d0/24.d0),h)
c  This is not very efficient, but we are trying to
c  make it work before making it fast.
c
c      write(6,*) 'result of drcat'
c      call pcmap(1,1,0,0,h,mh)
c
      return
      end
c
***********************************************************************

      subroutine cpsp(job,f,g,ans1,ans2,ans3,ans4)
c  this subroutine computes the scalar product of the two polynomials
c  f and g.
c  Written by Alex Dragt 10/23/89
c
      use lieaparam, only : monoms
      include 'impli.inc'
ccccc include 'param.inc'
      include 'expon.inc'
c
c  calling arrays
      dimension f(*), g(*)
ccccc dimension f(0:monoms), g(0:monoms)
c
      integer jf(6),jg(6)
c
c  compute scalar products
c
c  procedure for USp(6) invariant scalar product
c
      if(job .eq. 1) then
c
c  first order part
c
      ans1=0.d0
      do ind=1,6
      test=f(ind)*g(ind)
      if (test .ne. 0.d0) then
      do k=1,6
      jf(k)=expon(k,ind)
      enddo
      call msp1(jf,val)
      ans1=ans1+test*val
      endif
      enddo
c
c  second order part
c
      ans2=0.d0
      do ind=7,27
      test=f(ind)*g(ind)
      if (test .ne. 0.d0) then
      do k=1,6
      jf(k)=expon(k,ind)
      enddo
      call msp1(jf,val)
      ans2=ans2+test*val
      endif
      enddo
c
c  third order part
c
      ans3=0.d0
      do ind=28,83
      test=f(ind)*g(ind)
      if (test .ne. 0.d0) then
      do k=1,6
      jf(k)=expon(k,ind)
      enddo
      call msp1(jf,val)
      ans3=ans3+test*val
      endif
      enddo
c
c  fourth order part
c
      ans4=0.d0
      do ind=84,209
      test=f(ind)*g(ind)
      if (test .ne. 0.d0) then
      do k=1,6
      jf(k)=expon(k,ind)
      enddo
      call msp1(jf,val)
      ans4=ans4+test*val
      endif
      enddo
c
      endif
c
c  procedure for integration over S^5 invariant scalar product
c
      if(job .eq. 2) then
c
c  first order part
c
      ans1=0.d0
      do 10 indf=1,6
      do 10 indg=1,6
      test=f(indf)*g(indg)
      if (test .ne. 0.d0) then
      do 12 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   12 continue
      call msp2(jf,jg,val)
      ans1=ans1+test*val
      endif
   10 continue
      ans1=ans1/3.d0
c
c  second order part
c
      ans2=0.d0
      do 20 indf=7,27
      do 20 indg=7,27
      test=f(indf)*g(indg)
      if (test .ne. 0.) then
      do 22 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   22 continue
      call msp2(jf,jg,val)
      ans2=ans2+test*val
      endif
   20 continue
      ans2=ans2/12.d0
c
c  third order part
c
      ans3=0.d0
      do 30 indf=28,83
      do 30 indg=28,83
      test=f(indf)*g(indg)
      if (test .ne. 0.d0) then
      do 32 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   32 continue
      call msp2(jf,jg,val)
      ans3=ans3+test*val
      endif
   30 continue
      ans3=ans3/60.d0
c
c  fourth order part
c
      ans4=0.d0
      do 40 indf=84,209
      do 40 indg=84,209
      test=f(indf)*g(indg)
      if (test .ne. 0.) then
      do 42 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   42 continue
      call msp2(jf,jg,val)
      ans4=ans4+test*val
      endif
   40 continue
      ans4=ans4/360.d0
c
      endif
c
      return
      end

***********************************************************************
c
      subroutine old_cpsp(f,g,ans1,ans2,ans3,ans4)
c  this subroutine computes the scalar product of the two polynomials
c  f and g.
c  Written by Alex Dragt 10/23/89
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
c
c  calling arrays
ccccc dimension f(209), g(209)
      dimension f(*), g(*)
c
      integer jf(6),jg(6)
c
c  compute inner products
c
c  first order part
c
      ans1=0.d0
      do 10 indf=1,6
      do 10 indg=1,6
      test=f(indf)*g(indg)
      if (test .ne. 0.d0) then
      do 12 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   12 continue
      call msp(jf,jg,val)
      ans1=ans1+test*val
      endif
   10 continue
      ans1=ans1/3.d0
c
c  second order part
c
      ans2=0.d0
      do 20 indf=7,27
      do 20 indg=7,27
      test=f(indf)*g(indg)
      if (test .ne. 0.) then
      do 22 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   22 continue
      call msp(jf,jg,val)
      ans2=ans2+test*val
      endif
   20 continue
      ans2=ans2/12.d0
c
c  third order part
c
      ans3=0.d0
      do 30 indf=28,83
      do 30 indg=28,83
      test=f(indf)*g(indg)
      if (test .ne. 0.d0) then
      do 32 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   32 continue
      call msp(jf,jg,val)
      ans3=ans3+test*val
      endif
   30 continue
      ans3=ans3/60.d0
c
c  fourth order part
c
      ans4=0.d0
      do 40 indf=84,209
      do 40 indg=84,209
      test=f(indf)*g(indg)
      if (test .ne. 0.) then
      do 42 k=1,6
      jf(k)=expon(k,indf)
      jg(k)=expon(k,indg)
   42 continue
      call msp(jf,jg,val)
      ans4=ans4+test*val
      endif
   40 continue
      ans4=ans4/360.d0
c
      return
      end
c
***********************************************************************
c     
      subroutine evalf(zi,h,val2,val3,val4)
c this subroutine computes the value of the function h(zi)
c Written by Alex Dragt, Fall 1986, based on work of F. Neri
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension zi(6)
      dimension h(monoms),avect(monoms)
      include 'ind.inc'
c compute vector containing values of basis monomials
cryne call evalm(zi,vect)
      call evalm(zi,avect)
c     
c the following code has been commented out since it is replaced
c by the use of evalm
c compute linear monomials
c      do 10 i=1,6
c   10 avect(i) = zi(i)
c compute higher order monomials
c      do 20 i = 7,monoms
c   20 avect(i) = avect(index1(i))*avect(index2(i))
c     
c compute value of h
      val2=0.d0
      do 30 i=1,27
   30 val2=val2+h(i)*avect(i)
      val3=val2
      do 40 i=28,83
   40 val3=val3+h(i)*avect(i)
      val4=val3
      do 50 i=84,209
   50 val4=val4+h(i)*avect(i)
      return
      end
c
************************************************************
c
      subroutine evalf_old(zi,h,val2,val3,val4)
c this subroutine computes the value of the function h(zi)
c Written by Alex Dragt, Fall 1986, based on work of F. Neri
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension zi(6)
      dimension h(monoms),avect(monoms)
      include 'ind.inc'
c compute vector containing values of basis monomials
c compute linear monomials
      do 10 i=1,6
   10 avect(i) = zi(i)
c compute higher order monomials
      do 20 i = 7,monoms
   20 avect(i) = avect(index1(i))*avect(index2(i))
c compute value of h
      val2=0.d0
      do 30 i=1,27
   30 val2=val2+h(i)*avect(i)
      val3=val2
      do 40 i=28,83
   40 val3=val3+h(i)*avect(i)
      val4=val3
      do 50 i=84,209
   50 val4=val4+h(i)*avect(i)
      return
      end
c
********************************************************************
c
      subroutine evalm(zi,vect)
c this subroutine computes the values of the basis
c monomials and stores them in the vector vect.
c Written by Alex Dragt, 1 July 1991, based on work of F. Neri
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'ind.inc'
c     
c calling arrays
      dimension zi(6)
      dimension vect(monoms)
c     
c compute linear monomials
      do 10 i=1,6
   10 vect(i) = zi(i)
c compute higher order monomials
      do i = 7,27
        vect(i) = vect(index1(i))*vect(index2(i))
      enddo
      do i = 28,209
        vect(i) = vect(index1(i))*vect(index2(i))
      enddo
      do i = 210,monoms
        vect(i) = vect(index1(i))*vect(index2(i))
      enddo
c      
      return
      end
c
      subroutine exphf(h,ideg,f,maxf,trf)
c  Applies Exp(:h:) on polynomial f.
c  h is a polynomial of degree ideg.
c  f has terms from 1 thru maxf.
c  The result is trf, which has terms 1-maxf.
c
c   Written By F. Neri 9/26/86.
c
      include 'lims.inc'
      double precision f(923),h(923),trf(923)
      double precision tmpf1(923),tmpf2(923)
c  maxf has to be .le. 6.
      integer maxf
      integer maxpow, ifact, maxord
cryne 7/23/2002      integer bottom(0:12),top(0:12)
cryne 7/23/2002      common/lims/bottom,top
      do 1 i=1,top(maxf)
        tmpf1(i) = f(i)
   1  continue
      do 3 i = 1,top(maxf)
        trf(i) = f(i)
    3 continue
      maxpow = int((maxf-1)/(ideg-2))
      ifact = 1
      do 10 n=1,maxpow
        ifact =  ifact * (-n)
        maxord = int(maxf - (ideg-2) )
        do 11 i=1,top(maxf)
          tmpf2(i) = 0.d0
   11   continue
        do 20 iord=maxord,1,-1
          call pbkt(tmpf1,iord,h,ideg,tmpf2)
          call pmadd(tmpf2,iord+ideg-2,(1.d0/ifact),trf)
  20    continue
        do 30 i=1,top(maxf)
          tmpf1(i) = tmpf2(i)
   30   continue
  10  continue
      return
      end
c
******************************************************************
c
      subroutine fxform(ga,gm,fa,ha)
c this is a subroutine for transforming a function f.
c that is, it  computes h=exp(:g:)f where f is given
c by f=f2+f3+f4 and exp(:g:) denotes a general map.
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension ga(monoms),fa(monoms),ha(monoms),t1a(monoms),
     &          t2a(monoms),t3a(monoms)
      dimension gm(6,6),tm(6,6)
c clear arrays
      call clear(t1a,tm)
      call clear(t2a,tm)
      call clear(t3a,tm)
      call clear(ha,tm)
c compute [g3,f2]
      call pbkt(ga,3,fa,2,t1a)
c compute [g3,f3]
      call pbkt(ga,3,fa,3,t2a)
c compute [g3,[g3,f2]]
      call pbkt(ga,3,t1a,3,t3a)
c compute [g4,f2]
      call pbkt(ga,4,fa,2,t1a)
c accumulate results
c set up second order part
      do 10 i=7,27
   10 t1a(i)=fa(i)
c set up third order part
      do 20 i=28,83
   20 t1a(i)=fa(i) + t1a(i)
c set up fourth order part
      do 30 i=84,209
   30 t1a(i)=fa(i) + t2a(i) + t3a(i)/2.d0 + t1a(i)
c transform results by gm
      call xform(t1a,2,gm,0,ha)
      call xform(t1a,3,gm,1,ha)
      call xform(t1a,4,gm,1,ha)
      return
      end
c
      subroutine mapmap(rh,rmh,th,tmh)
c  Written by Rob Ryne, ca 1982
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension rh(monoms),th(monoms),rmh(6,6),tmh(6,6)
      do 10 i=1,6
      do 10 j=1,6
   10 tmh(i,j)=rmh(i,j)
      do 20 i=1,monoms
   20 th(i)=rh(i)
c
      return
      end
c
***********************************************************************
c
      subroutine matify(matrix,f2)
c  Computes the matrix that corresponds to :f2:.
c  It is written in a simple-minded manner to keep execution time short.
c  Written by Liam Healy, May 29, 1985.
c
c----Variables----
c     implicit none
c matrix = matrix supplied
      double precision matrix(6,6)
c f2 = array of coefficients giving f2 values (others are ignored)
      double precision f2(*)
c
c----Routine----
      matrix(1,1)=-f2(8)
      matrix(1,2)=-2.*f2(13)
      matrix(1,3)=-f2(14)
      matrix(1,4)=-f2(15)
      matrix(1,5)=-f2(16)
      matrix(1,6)=-f2(17)
      matrix(2,1)=2.*f2(7)
      matrix(2,2)=f2(8)
      matrix(2,3)=f2(9)
      matrix(2,4)=f2(10)
      matrix(2,5)=f2(11)
      matrix(2,6)=f2(12)
      matrix(3,1)=-f2(10)
      matrix(3,2)=-f2(15)
      matrix(3,3)=-f2(19)
      matrix(3,4)=-2.*f2(22)
      matrix(3,5)=-f2(23)
      matrix(3,6)=-f2(24)
      matrix(4,1)=f2(9)
      matrix(4,2)=f2(14)
      matrix(4,3)=2*f2(18)
      matrix(4,4)=f2(19)
      matrix(4,5)=f2(20)
      matrix(4,6)=f2(21)
      matrix(5,1)=-f2(12)
      matrix(5,2)=-f2(17)
      matrix(5,3)=-f2(21)
      matrix(5,4)=-f2(24)
      matrix(5,5)=-f2(26)
      matrix(5,6)=-2.*f2(27)
      matrix(6,1)=f2(11)
      matrix(6,2)=f2(16)
      matrix(6,3)=f2(20)
      matrix(6,4)=f2(23)
      matrix(6,5)=2.*f2(25)
      matrix(6,6)=f2(26)
      return
      end
c
***********************************************************************
c
      subroutine mclear(mh)
c 
c Clears to zero the matrix mh.
c  Written by Liam Healy, June 12, 1984.
c
      double precision mh(6,6)
c
      do 100 i=1,6
      do 100 j=1,6
  100 mh(i,j)=0.d0
      return
      end
c
***********************************************************************
c
c
***********************************************************************
c
      subroutine mident(mh)
c 
c  Sets mh to the identity matrix.
c
      double precision mh(6,6)
c
      call mclear(mh)
c
      do 100 i=1,6
  100 mh(i,i)=1.d0
c
      return
      end
c
************************************************************************
c
      subroutine minv(fm)
c  Returns the inverse of the matrix fm on the
c  assumption that fm is symplectic.
c  Written by Alex Dragt, 4 October 1989.
c
      include 'impli.inc'
      include 'symp.inc'
c
c  Calling array
      dimension fm(6,6)
c
c  Temporary arrays
      dimension temp(6,6)
c
c----Routine----
c
c  Calculate (fm transpose)*jm
c
      call mclear(temp)
      do 120 j=1,6
      do 120 k=1,6
      do 120 l=1,6
  120  temp(j,k)=temp(j,k)+fm(l,j)*jm(l,k)
c
c  Calculate (jm transpose)*(fm transpose)*jm
c
      call mclear(fm)
      do 140 j=1,6
      do 140 l=1,6
      do 140 k=1,6
  140  fm(j,k)=fm(j,k)+jm(l,j)*temp(l,k)
c
      return
      end
c
***********************************************************************
c
      subroutine mycat(maxcat,f,mf,g,mg,h,mh)
c     concatenates a map with linear piece
c     represented by a matrix mf and nonlinearities
c     represented by exp:f:, with a map whose
c     linear piece has matrix representation mg
c     and whose nonlinearities are represented by exp:g:
c
c     the result is a map with linearities possessing
c     a matrix representation mh=mg*mf
c     and with nonlinearities generated by exp:h:
c
c     the "f" map is assumed to occur first in the beamline
c     the "g" map is encountered second
c
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
      include 'symp.inc'
      include 'ind.inc'
      double precision j4,mf,mg,mh,m,l3,l4
      dimension mf(6,6),mg(6,6),mh(6,6),temp(6,6),m(6,6)
      dimension f(monoms), g(monoms), h(monoms)
      dimension f3(monom2),f4(monom2),f5(monom1),f6(monoms)
      dimension t3(monom2),t4(monom2),t5(monom1),t6(monoms)
cryne 7/23/2002      implicit double precision (a-h,o-z)
cryne 7/23/2002      dimension f(923), g(923), h(923)
cryne 7/23/2002      dimension f3(209),f4(209),f5(461),f6(923)
cryne 7/23/2002      dimension t3(209),t4(209),t5(461),t6(923)
cryne 7/23/2002      common /ind/imaxi,jv(923),index1(923),index2(923)
c
      write(6,*) 'first factor in mycat'
      call pcmap(1,1,0,0,f,mf)
      write(6,*) 'second factor in mycat'
      call pcmap(1,1,0,0,g,mg)
c
      call ident(h,mh)
      do 666 i1=1,209
        t3(i1) = 0.0d0
  666 continue
      do 777 i1=1,6
      do 777 i2=1,6
      m(i1,i2)=0.d0
      mh(i1,i2) = 0.d0
      temp(i1,i2) = 0.d0
 777  continue
c
c     compute mh=mg*mf and temp=mgtransposed*jm
c
      do 40 j=1,6
      do 30 k=1,6
      do 20 l=1,6
      mh(j,k)=mh(j,k)+mg(j,l)*mf(l,k)
      temp(j,k)=temp(j,k)+mg(l,j)*jm(l,k)
   20 continue
   30 continue
   40 continue
c
c     compute m=inverse of mg=jmtransposed*temp
c
      do 41 j=1,6
      do 31 k=1,6
      do 21 l=1,6
      m(j,k)=m(j,k)+jm(l,j)*temp(l,k)
   21 continue
   31 continue
   41 continue
c
c     compute transformed arrays
c
       call xform5(f,3,m,f3)
       call xform5(f,4,m,f4)
       if(maxcat.gt.4) call xform5(f,5,m,f5)
       if(maxcat.gt.5) call xform5(f,6,m,f6)
c    third order terms
       call pmadd(g,3,1.d0,h)
       call pmadd(f3,3,1.d0,h)
c   fourth order terms
       call pbkt1(f3,3,g,3,t4)
       call pmadd(g,4,1.d0,h)
       call pmadd(f4,4,1.d0,h)
       call pmadd(t4,4,.5d0,h)
c    fifth order terms
      if(maxcat.gt.4) then
        call pmadd(g,5,1.d0,h)
        call pmadd(f5,5,1.d0,h)
        call pmadd(g,3,1.d0,t3)
        call pmadd(f3,3,1.d0/2.d0,t3)
        call pbkt1(t3,3,t4,4,t5)
        call pmadd(t5,5,-1.d0/3.d0,h)
        call pbkt1(g,3,f4,4,t5)
        call pmadd(t5,5,-1.d0,h)
      endif
c    sixth order terms
      if(maxcat.gt.5) then
        call pbkt1(g,4,f4,4,t6)
        call pmadd(t6,6,-.5d0,h)
        call pmadd(g,6,1.d0,h)
        call pmadd(f6,6,1.d0,h)
        call pbkt1(g,3,f5,5,t6)
        call pmadd(t6,6,-1.d0,h)
        call pbkt1(f3,3,g,3,t4)
        call pbkt1(g,3,t4,4,t5)
        call pbkt1(f3,3,t5,5,t6)
        call pmadd(t6,6,(-1.d0/24.d0),h)
        call pbkt1(h,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(1.d0/12.d0),h)
        call pbkt1(g,3,f4,4,t5)
        call pbkt1(g,3,t5,5,t6)
        call pmadd(t6,6,(.5d0),h)
        call pbkt1(f4,4,t4,4,t6)
        call pmadd(t6,6,-.25d0,h)
        call pbkt1(g,4,t4,4,t6)
        call pmadd(t6,6,-.25d0,h)
        call pbkt1(g,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(1.d0/24.d0),h)
        call pbkt1(f3,3,t4,4,t5)
        call pbkt1(h,3,t5,5,t6)
        call pmadd(t6,6,(-1.d0/24.d0),h)
c  This is not very efficient, but we are trying to
c  make it work before making it fast.
      endif
c
      write(6,*) 'result of mycat'
      call pcmap(1,1,0,0,h,mh)
c
      return
      end
c
***********************************************************************
c
      subroutine msp(j1,j2,val)
c  Calculated the scalar product of two monomials having exponent
c  arrays j1 and j2.
c  Written by Alex Dragt 10/22/89
c
      include 'impli.inc'
c
c  calling arrays
      integer j1(6),j2(6)
c
      integer jt(6)
c
      do 10 i=1,6
   10 jt(i)=j1(i)+j2(i)
      call s5i(jt,val)
c
      return
      end
c
***********************************************************************
c
      subroutine msp1(j1,val)
c  For a monomial having exponent array j1, calculates
c  the USp(6) invariant scalar product of the monomial
c  with itself.
c  Written by Alex Dragt 2/18/03.
c
      include 'impli.inc'
c
c  calling arrays
      integer j1(6)
c
c  procedure for computing j1!
c
      val=1.d0
      do i=1,6
      ipow=j1(i)
      call factfn(ipow,coefi)
      val=val*coefi
      enddo
c
      return
      end
c
***********************************************************************
c
      subroutine msp2(j1,j2,val)
c  Calculates the S^5 invariant scalar product of two monomials
c  having exponent arrays j1 and j2.
c  Written by Alex Dragt 10/22/89
c
      include 'impli.inc'
c
c  calling arrays
      integer j1(6),j2(6)
c
      integer jt(6)
c
      do 10 i=1,6
   10 jt(i)=j1(i)+j2(i)
      call s5i(jt,val)
c
      return
      end
***********************************************************************
c
      subroutine factfn(intg,val)
c
c This subroutine finds val=intg!
c Written by Alex Dragt 2/18/03
c
      include 'impli.inc'
c
      dimension table(0:6)
c
      data table/1.d0,1.d0,2.d0,6.d0,24.d0,120.d0,720.d0/
c
      val=table(intg)
c
      return
      end

***********************************************************************
c
      subroutine pbkt(f,ordf,g,ordg,pb)
c  Calculates the Poisson Bracket of the order ordf part of f with
c  the order ordg part of g, leaving the result in pb.
c  Written by Liam Healy, November 26, 1984.
c
      use lieaparam, only : monoms
c----Variables----
c  f,g,pb = two arrays and their Poisson bracket
      double precision f(*),g(*),pb(*)
c  ordf,ordg = order desired from f and g
      integer ordf,ordg
c  indf,indg,indpb = index in f,g and pb
      integer indf,indg,indpb
c  pr = array of exponents for the product of f(indf) and g(indg)
c  prd =copy of pr
c  prx, prp = exponents of particular x, p variables, from pr
      integer pr(6),prd(6),prx,prp
c  xvb, pvb, vbl = x variable, p variable (pvb=xvb+1), variable number
      integer xvb,pvb,vbl
      include 'expon.inc'
      include 'lims.inc'
c
c----Routine----
c  initialize array pb
      do 80 indpb=bottom(ordf+ordg-2),top(ordf+ordg-2)
   80   pb(indpb)=0.
c  pick individual indf and indg, find what element of pb it affects,
c  and calculate the new value.  Loop for all indeces in the specified
c  orders.
      do 100 indf=bottom(ordf),top(ordf)
        if (f(indf).ne.0.) then
          do 120 indg=bottom(ordg),top(ordg)
            if (g(indg).ne.0.) then
c               load sum of exponents into pr
              do 140 vbl=1,6
  140           pr(vbl)=expon(vbl,indf)+expon(vbl,indg)
              do 160 xvb=1,5,2
c                 go through three axes; for each one that pb can
c                 be taken, modify appropriate element of pb.
                pvb=xvb+1
                do 180 vbl=1,6
  180             prd(vbl)=pr(vbl)
                prx=prd(xvb)
                prp=prd(pvb)
                if (prx*prp.gt.0) then
                  prd(xvb)=prx-1
                  prd(pvb)=prp-1
                  indpb=ndex(prd)
                  pb(indpb)=pb(indpb)+f(indf)*g(indg)
     &               *(expon(xvb,indf)*expon(pvb,indg)
     &               -expon(pvb,indf)*expon(xvb,indg))
                endif
  160         continue
            endif
  120     continue
        endif
  100 continue
      return
      end
c
      subroutine pbkt1(f,ordf,g,ordg,pb)
c  Calculates the Poisson Bracket of the order ordf part of f with
c  the order ordg part of g, leaving the result in pb.
c
      use lieaparam, only : monoms,monom1
      include 'iprod.inc'
      include 'expon.inc'
      include 'lims.inc'
      include 'prodex.inc'
c----Variables----
c  f,g,pb = two arrays and their Poisson bracket
      double precision f(*),g(*),pb(*)
c  ordf,ordg = order desired from f and g
      integer ordf,ordg
c  indf,indg,indpb = index in f,g and pb
      integer indf,indg,indpb
c  pr = array of exponents for the product of f(indf) and g(indg)
c  prd =copy of pr
c  prx, prp = exponents of particular x, p variables, from pr
      integer pr(6),prd(6),prx,prp
c  xvb, pvb, vbl = x variable, p variable (pvb=xvb+1), variable number
      integer xvb,pvb,vbl
c  expon = table of exponents
cryne 7/23/2002      integer expon(6,0:923)
cryne 7/23/2002      common/expon/expon
c  bottom, top = lowest and highest monomial number for each order
cryne 7/23/2002      integer bottom(0:12),top(0:12)
cryne 7/23/2002      common/lims/bottom,top
cryne 7/23/2002      integer prodex(6,0:923)
cryne 7/23/2002      common /prodex/prodex
      integer vbf,vbg
      integer conj(6),sigj(6)
      data conj /2,1,4,3,6,5/
      data sigj /1,-1,1,-1,1,-1/
      save conj,sigj           !cryne 7/23/3003
c
c----Routine----
c  initialize array pb
      do 80 indpb=bottom(ordf+ordg-2),top(ordf+ordg-2)
   80   pb(indpb)=0.
c  pick individual indf and indg, find what element of pb it affects,
c  and calculate the new value.  Loop for all indeces in the specified orders.
      do 160 vbf = 1,6
        vbg = conj(vbf)
        sign = sigj(vbf)
        do 100 indf=bottom(ordf-1),top(ordf-1)
          if(f(prodex(vbf,indf)).eq.0.0d0 ) goto 100
          do 110 indg=bottom(ordg-1),top(ordg-1)
            indpb = iprod(indg,indf)
            pb(indpb) = pb(indpb) +
     &      sign * f(prodex(vbf,indf)) * g(prodex(vbg,indg))
     &      * (expon(vbf,indf)+1) * (expon(vbg,indg)+1)
  110     continue
  100   continue
  160 continue
      return
      end
c
      subroutine pbkt2(f,n,index,bpb)
c     computes the poisson bracket [f,z(i)] of
c     an n-th degree polynomial whose coefficients
c     are stored in an array f with the i-th component
c     of the dynamical variable array z.
c     z(1)=x
c     z(2)=px
c     .
c     .
c     .
c
c     the coefficients of the result are
c     stored in the array bpb.  note the
c     result is a polynomial of degree n-1
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      include 'vblist.inc'
      include 'len.inc'
cryne 7/23/2002      implicit double precision (a-h,o-z)
cryne 7/23/2002      integer expon(6,0:923),vblist(6,0:923)
cryne 7/23/2002      dimension f(923),bpb(923),j(6)
cryne 7/23/2002      common/expon/expon
cryne 7/23/2002      common/vblist/vblist
      dimension f(monoms),bpb(monoms),j(6)
      dimension zz(6)
cryne 7/23/2002      common /len/ len(16)
      do 10 i = 1,6
          zz(i) = 0.0d0
   10 continue
      zz(index) = 1.0d0
      call pbkt1(f,n,zz,1,bpb)
c
      return
      end
c
************************************************************************
c
      subroutine polr(fa,fm,rm,pdsm,reval,revec)
c  this subroutine finds the polar decomposition of a symplectic map
c  fa,fm is the incoming matrix, and is left unchanged
c  rm and pdsm are its orthogonal (rotation) and positive definite
c  symmetric factors
c  reval and revec are the eigenvalues and eigenvector array for pdsm
c  Written by Alex Dragt, Fall 1986
c
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms)
      dimension fm(6,6),rm(6,6),pdsm(6,6),reval(6),revec(6,6)
      dimension t1m(6,6),t2m(6,6)
      dimension ta(monoms)
c
c clear array
      call clear(ta,t1m)
c
c begin computation
c
c  computation of positive definite symmetric factor
c  along with its eigenvalues and eigenvectors
c
c  find (fm transpose)*fm
      call matmat(fm,t1m)
      call mtran(t1m)
      call mmult(t1m,fm,t2m)
c  extract square root
      call seig6(t2m,reval,revec)
      call mclear(t2m)
      do 10 i=1,6
      reval(i)=sqrt(reval(i))
      t2m(i,i)=reval(i)
   10 continue
      call matmat(revec,t1m)
      call inv(ta,t1m)
      call sndwch(ta,t1m,ta,t2m,ta,pdsm)
c
c  computation of orthogonal (rotation) factor
c
      call matmat(pdsm,t1m)
      call inv(ta,t1m)
      call concat(ta,t1m,ta,fm,ta,rm)
c
      return
      end
c
***********************************************************************
c
      subroutine pprod(a,na,b,nb,c)
c this subroutine computes the product of two polynomials: c=a*b
c Written by Alex Dragt, Fall 1986, based on work of F. Neri
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      include 'lims.inc'
      dimension a(monoms),b(monoms),c(monoms),l(6)
c
      do 200 ia=bottom(na),top(na)
           if(a(ia).eq.0.d0) goto 200
           do 20 ib = bottom(nb),top(nb)
               if(b(ib).eq.0.d0) goto 20
               do 2 m=1,6
                   l(m) = expon(m,ia) +  expon(m,ib)
   2           continue
               n = ndex(l)
               c(n) = c(n) + a(ia)*b(ib)
  20       continue
 200  continue
c
      return
      end
c
***********************************************************************
c
      subroutine s2f(revec,k1,k2,val)
c  this subroutine evaluates the symplectic 2-form assiciated with J
c  Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension revec(6,6)
c
c  begin calculation
      val=0.d0
      do 10 ip=2,6,2
      iq=ip-1
      val=val+revec(iq,k1)*revec(ip,k2)-revec(ip,k1)*revec(iq,k2)
   10 continue
c
      return
      end
c
***********************************************************************
c
      subroutine s5i(j,ans)
c Calculates the integral over S5 of the monomial having exponents
c stored in j.
c Written by Alex Dragt 10/22/89
c
      include 'impli.inc'
c
c     calling array
      integer j(6)
c
      dimension anst(6)
c
      ans=0.
      do 10 i=1,6
      tans=0.d0
      ji=j(i)
      if (ji .eq. 0) tans=1.d0
      if (ji .eq. 2) tans=1.d0/2.d0
      if (ji .eq. 4) tans=3.d0/4.d0
      if (ji .eq. 6) tans=15.d0/8.d0
      if (ji .eq. 8) tans=105.d0/16.d0
      if (tans .eq. 0.d0) return
      anst(i)=tans
   10 continue
      ans=anst(1)*anst(2)*anst(3)*anst(4)*anst(5)*anst(6)
c
      return
      end
c
***********************************************************************
c
      subroutine scncat(fa,fm,ga,gm,ha,hm)
c     this is a special routine for concatenation.
c     it calls concat and then transfers the
c     coefficients for the linear and quadratic
c     polynomials.
c     Written by Alex Dragt, Fall 1985
c
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension fa(monoms),ga(monoms),ha(monoms)
      dimension fm(6,6),gm(6,6),hm(6,6)
      call clear(ha,hm)
      call concat(fa,fm,ga,gm,ha,hm)
c  transfer coefficients for linear and quadratic polynomials
      do 10 i=1,27
      ha(i)=ga(i)
   10 continue
      return
      end
c
***********************************************************************
c
      subroutine smtof(scale,tm,ta)
c converts the symmetric matrix scale*tm to the array ta.
c written by Alex Dragt 22 May 1991
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c calling arrays
      dimension ta(monoms), tm(6,6)
c
c procedure
c clear array
      do 5 i=1,monoms
    5 ta(i)=0.d0
c put in new contents
      ta(7)=tm(1,1)
      ta(8)=tm(1,2)
      ta(9)=tm(1,3)
      ta(10)=tm(1,4)
      ta(11)=tm(1,5)
      ta(12)=tm(1,6)
      ta(13)=tm(2,2)
      ta(14)=tm(2,3)
      ta(15)=tm(2,4)
      ta(16)=tm(2,5)
      ta(17)=tm(2,6)
      ta(18)=tm(3,3)
      ta(19)=tm(3,4)
      ta(20)=tm(3,5)
      ta(21)=tm(3,6)
      ta(22)=tm(4,4)
      ta(23)=tm(4,5)
      ta(24)=tm(4,6)
      ta(25)=tm(5,5)
      ta(26)=tm(5,6)
      ta(27)=tm(6,6)
c
      do 10 i=7,27
  10  ta(i)=scale*ta(i)
c
      return
      end
c
***********************************************************************
c
      subroutine sndwch(t1a,t1m,t2a,t2m,t3a,t3m)
c this is a subroutine for sandwiching a map
c it computes t3=t1*t2*(t1 inverse)
c Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension t1a(monoms),t2a(monoms),t3a(monoms),t4a(monoms),
     &          t5a(monoms)
      dimension t1m(6,6),t2m(6,6),t3m(6,6),t4m(6,6),t5m(6,6)
      call mapmap(t1a,t1m,t4a,t4m)
      call inv(t4a,t4m)
      call concat(t1a,t1m,t2a,t2m,t5a,t5m)
      call concat(t5a,t5m,t4a,t4m,t3a,t3m)
      return
      end
c
************************************************************************
      subroutine sndwchi(t1a,t1m,t2a,t2m,t3a,t3m)
c added to liea.f by rdr and ctm on 5/22/02
cryne 08/17/2001 This version (sndwchi) reverses the order
c
c this is a subroutine for sandwiching a map
c it computes t3=t1*t2*(t1 inverse)
c Written by Alex Dragt, Fall 1986
c Modified by Alex Dragt, 18 July 1988
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c calling arrays
      dimension t1a(monoms),t2a(monoms),t3a(monoms)
      dimension t1m(6,6),t2m(6,6),t3m(6,6)
c
c working arrays
      dimension ta(monoms)
      dimension tm(6,6)
c
      call mapmap(t1a,t1m,ta,tm)
      call inv(ta,tm)
      call concat(ta,tm,t2a,t2m,ta,tm)
      call concat(ta,tm,t1a,t1m,t3a,t3m)
c
      return
      end
c
***********************************************************************
c
      subroutine svpbkt(in,ord,psv,base, pb)
c  Single Variable Poisson BracKeT
c  Takes the poisson bracket of the order 'ord' part of 'in' with the
c  phase space variable represented by 'psv'.  The result is left in pb.
c  This subroutine replaces DRD's pbkt2 and was written by Liam Healy
c  on November 29, 1984 from an idea of Christoph Iselin.
c
      use lieaparam, only : monoms
c----Variables----
c  base = lowest index of pb
      integer base
c  in, pb = incoming array, outgoing pb array
      double precision in(*),pb(base:*)
c  ord = order of 'in' to be pb-ed
      integer ord
c  psv,wrt = phase space variable in pb,
c            variable to be diff with resp to
      integer psv,wrt
c  indin,indpb = indeces for 'in' and 'pb'
      integer indin,indpb
c  sign = multiplier +1 or -1 depending on whether x or p
c         differentiating
      double precision sign
      include 'expon.inc'
      include 'lims.inc'
c
c----Routine----
      if (mod(psv,2).eq.1) then
c  psv is a coordinate, diff wrt the canonical momentum and flip sign
        wrt=psv+1
        sign=-1.
      else
c  psv is a momentum, diff wrt the canonical coordinate
        wrt=psv-1
        sign=+1.
      endif
      indin=bottom(ord)
      do 120 indpb=bottom(ord-1),top(ord-1)
  100   if (expon(wrt,indin).eq.0) then
          indin=indin+1
          if (indin.le.top(ord)) goto 100
          return
        endif
        pb(indpb)=sign*expon(wrt,indin)*in(indin)
  120   indin=indin+1
      return
      end
c
***********************************************************************
c
      subroutine xform(in,ord,matrix,matold,out)
c  Transforms the order 'ord' part of polynomial represented by array
c  'in' with the matrix 'matrix', and leaves the resulting polynomial
c  in the array 'out'.  'matold' should be set >=1 if the
c  matrix in the previous call to xform was the same as this call.
c  This subroutine replaces DRD's xform and was written by
c  Liam Healy on November 30, 1984.
c
      use lieaparam, only : monoms
c----Variables----
c in, out = original polynomial array, transformed polynomial array
      double precision in(*),out(*)
c ord = order to be transformed
      integer ord
c matrix = matrix that represents the linear transformation
c prod = product, used to accumulate matrix elements
      double precision matrix(6,6),prod
c indin, indout = indices for arrays 'in' and 'out'
      integer indin,indout
c colme(row,count) = column number of the count-th one
      integer colme(6,0:6),size
      save colme
c matold = matrix supplied was used in the last call to xform
      integer matold
c posn,psv = position in incoming variable list, phase space variable
      integer posn,psv
c code = packed information on which combination of matrix elts to
c        select
c ncombs = number of different combinations of me's that can be selected
c rem= remainder, gives the code information for the remaining positions
      integer code,ncombs,rem
c row, this = row for this position, ordinal of non-zero me at this posn
c col = column
      integer row(6),this(6),col
      include 'vblist.inc'
      include 'prodex.inc'
      include 'lims.inc'
c
c rowa, cola = row and column in matrix
c count = ordinal of the non-zero matrix elements for a given row
      integer count,rowa,cola
c
c------------------------
c Find non-zero matrix elements
c Analyzes the matrix to find its non-zero entries, which are
c recorded in the array 'colme' by row number, together with
c the total number of non-zero entries in each row in colme(rowa,0).
      if (matold.le.0) then
        do 100 rowa=1,6
          count=0
          do 110 cola=1,6
            if(matrix(rowa,cola).ne.0.) then
              count=count+1
              colme(rowa,count)=cola
            endif
  110     continue
          colme(rowa,0)=count
  100   continue
      endif
c
c---------------
c Clear the 'out' array elements of this order
      do 80 indout=bottom(ord),top(ord)
   80   out(indout)=0.
c------------------------
c Cycle through incoming array elements, for each that are non-zero
c Map them to outgoing elements
      do 120 indin=bottom(ord),top(ord)
        if (in(indin).ne.0.) then
          ncombs=1
          do 140 posn=1,ord
            row(posn)=vblist(posn,indin)
  140       ncombs=ncombs*colme(row(posn),0)
c------------------------
c Go through all possible combinations of matrix elements that
c have the correct rows, i.e. corrosponding to the variables in 'indin'
c 'Code' contains packed (by variable base) information on which me to
c select.
          do 160 code=0,ncombs-1
            rem=code
            prod=1.
            indout=0
            do 200 posn=1,ord
              size=colme(row(posn),0)
              this(posn)=mod(rem,size)+1
              rem=rem/size
              col=colme(row(posn),this(posn))
              indout=prodex(col,indout)
              prod=prod*matrix(row(posn),col)
  200       continue
            out(indout)=out(indout)+in(indin)*prod
  160     continue
c------------------------
        endif
  120 continue
      return
      end
c
******************************************************************
c
      subroutine xform5(f,no,m,l)
c
c     Transforms arguments of the polynomial f of degree no
c     by the matrix m.  The coefficients of the resultant 
c     polynomial are stored in the the array l which
c     thus contains the coefficients of f(m*z).
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      include 'vblist.inc'
      include 'prodex.inc'
      include 'len.inc'
      include 'ind.inc'
c
      double precision l,m,f
cryne 7/23/2002      dimension f(923),l(923)
cryne 7/23/2002      double precision temp(923)
cryne dimension f(monoms),l(monoms),temp(monoms)
      dimension f(*),l(*),temp(monoms)
      dimension m(6,6)
c
c      do 127 my=1,6
c      write(6,*) my, len(my),
c     & vblist(my,211), prodex(my,84)
c 127  continue
c      return
c
c      write(6,*) 'm and f coming into xform5'
c      write(6,*) 'with no =',no
c      call pcmap(1,1,0,0,f,m)
c
c     initialise arrays
c
      do 10 kp=1,len(no)
      l(kp)=0.0d0
      temp(kp) = 0.0d0
   10 continue
      do 100 n= len(no-1)+1,len(no)
        if(f(n).eq.0.0d0) goto 100
        do 101 kp=1,len(no-1)
          temp(kp) = 0.0d0
  101   continue
        do 102 k=1,6
          if(m(vblist(1,n),k).eq.0.0d0) goto 102
          temp(k) = f(n)*m(vblist(1,n),k)
  102    continue
        do 110 ior=1,no-1
          k1 = vblist(ior+1,n)
          if(ior.eq.1) then
            n1 = 1
          else
            n1=len(ior-1)+1
          endif
          do 112 k=1,6
            if(m(k1,k).eq.0.0d0) goto 112
            xm = m(k1,k)
cryne 7/21/01  cdir$ ivdep
            do 111 nn=n1,len(ior)
            temp(prodex(k,nn)) = temp(prodex(k,nn)) + temp(nn)*xm
  111       continue
  112     continue
  110   continue
  100 continue
      do 200 nn=len(no-1)+1,len(no)
        l(nn) = l(nn)+temp(nn)
  200 continue
c
c      write(6,*) 'm and l coming out of xform5'
c      write(6,*) 'with no =',no
c      call pcmap(1,1,0,0,l,m)

      return
      end
c end of file
