***********************************************************************
c
c    THIS FILE CONTAINS THE ROUTINES THAT IMPLEMENT THE SCALING 
c    SPLITTING SQUARING (SSS) ALGORITHM FOR MARYLIE5.0
c    (Written 23 April 97)
c    (Fixed  8 May 97)
c    (Modified  4 Aug 97)
***********************************************************************
c
        subroutine  sss(p,ga,gm)
c
c  This routine performs the Dragt-Finn factorization of the map
c  exp(-t:h:) using the Scaling Splitting, Squaring algorithm.
c  The arrays ga initially contains the Hamiltonian h.
c
c  Two options are possible:
c  1. ncheck = 0   ;   the calculation of the  map is done only once.
c  The number  of squarings np is fixed by the requirement that
c  |norm(JS)*.5**np| < eps, with eps specified by the user.
c  2. ncheck =/= 0 ;  the calculation of the map is repeated with np=np+1
c  and the difference between  selected generators of the map
c  calculated with np and np+1 squarings checked.
c  If the relative difference is greater than 'err' the calculation is 
c  repeated by increasing np by one.  If the relative difference increases 
c  with np (round off errors) the algorithm is stopped.
c
c  Written by M.Venturini April 23, 97.
c  Modified Aug. 4, 97  (M.V.)
c
c***************************
c  err                = error allowed for the map generators
c  integration length = t
c  np                 = # of squaring
c  eps                = max norm of the linear part of the generators
c  npse               = n. of squaring set by the user regardless of err
c                       (if npse=0 np is set based on err)
c***************************
c
       use lieaparam, only : monoms
       implicit double precision (a-h,o-z)
       include 'files.inc'
c
       parameter(err=1.d-13)
c
       character*3 kynd
       dimension p(6)
       dimension fa(monoms),fm(6,6)
       dimension ga(monoms),gm(6,6)
       dimension fa2(monoms),fm2(6,6)
       dimension fap(monoms),fmp(6,6)
       dimension faw(monoms),fmw(6,6)
c
       t=p(1)
       nmapf=nint(p(2))
       nmapg=nint(p(3))
       eps=p(4)
       ncheck=nint(p(5))
       npset=nint(p(6))
c
c Get the array fa
       if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
       if (nmapg.ge.1 .and. nmapf.le.5) then
       kynd='gtm'
       call strget(kynd,nmapf,fa,fm)
       endif
c
c Multiply fa by t
       call csmul(t,fa,fa)
c
       If(npset.gt.0) then
       np=nint(p(6))
       ncheck=0
cryne       write(jof,*) '    '
cryne       write(jof,*) 'number of squaring np=',np
       write(jodf,*) '    '
       write(jodf,*) 'number of squaring np=',np
       goto 21
       endif
c
c Compute the number of squarings based on eps
c
       call matify(fm,fa)
       call mnorm(fm,res)
       np=1
       scale=.5d0
  10   continue
       test=6*res*scale
       if (test .lt. eps) goto 20
       np=np+1
       scale=scale/2.d0
       goto 10
  20   continue
c
cryne       write(jof,*) '    '
cryne       write(jof,*) 'number of squaring np=',np
cryne       write(jof,*) 'norm of the linear part of the map generators',res
cryne       write(jof,*) 'scale',scale
cryne       write(jof,*) '    '
c
       write(jodf,*) '    '
       write(jodf,*) 'number of squaring np=',np
       write(jodf,*) 'norm of the linear part of the map generators',res
       write(jodf,*) 'scale',scale
       write(jodf,*) '    '
cryne 6/21/2002
cryne put in a diagnostic until we are sure that these routines are
cryne robust when ncheck=0:
       if(ncheck.eq.0 .and. np.lt.5)then
       write(jof,*)'WARNING FROM SUBROUTINE SSS: # OF SPLITTINGS=',np
       write(jodf,*)'WARNING FROM SUBROUTINE SSS: # OF SPLITTINGS=',np
       endif
c
       if (ncheck.ne.0) then
       call mapmap(fa,fm,fap,fmp)
       endif
c
  21   tau=1/2.d0**np
       call splitT(np,tau,fa,fm)    
c
c Concatenate
c
       do 30  i=1,np
         call concat(fa,fm,fa,fm,fa2,fm2)
         call mapmap(fa2,fm2,fa,fm)
 30    continue
c
c Want to increase np?
c
       if (ncheck.eq.0) goto 70
c
       reldifbef=1
  50   continue
       fbef463 = fa(463)
       np=np+1
       tau=1/2.d0**np
       call mapmap(fap,fmp,fa,fm)
       call splitT(np,tau,fa,fm) 
c
       do 31  i=1,np
         call concat(fa,fm,fa,fm,fa2,fm2)
         call mapmap(fa2,fm2,fa,fm)
 31    continue
c
c
       diff463 = fa(463) -fbef463
       reldif=diff463/fa(463)
c
       write(jof,*)'np=',np,'       f(463)=',fa(463)
       write(jof,*)'np=',np-1,'     f(463)=',fbef463
       write(jof,*)'difference =',diff463
       write(jof,*)'relative error =',reldif
       write(jof,*)'  '
c
       write(jodf,*)'np=',np,'       f(463)=',fa(463)
       write(jodf,*)'np=',np-1,'     f(463)=',fbef463
       write(jodf,*)'difference =',diff463
       write(jodf,*)'relative error =',reldif
       write(jodf,*)'  '

c
c Stop the calculation if the relative error increases with np
c
          if (abs(reldifbef).lt. abs(reldif)) then
c
          write(jof,*) 'Increasing np (n. of squaring) does not'
          write(jof,*) 'improve the accuracy.'
          write(jof,*) 'Stop at np=',np-1
          write(jof,*) '  '
c
          write(jodf,*) 'Increasing np (n. of squaring) does not'
          write(jodf,*) 'improve the accuracy.'
          write(jodf,*) 'Stop at np=',np-1
          write(jodf,*) '  '
c
          call mapmap(faw,fmw,fa,fm)
          goto  70
          endif
c
               call mapmap(fa,fm,faw,fmw)
               reldifbef=reldif
c
          if (abs(reldif).le. err)  then
          goto 70
          endif
c
       goto 50
c
c Set the output
c
  70   if (nmapg.ge.1 .and. nmapg.le.5) then
       kynd='stm'
       call strget(kynd,nmapg,fa,fm)
       endif
c
       if (nmapg.eq.0) call mapmap(fa,fm,ga,gm)
c
       return
       end
c
***********************************************************************
c
        subroutine splitT(np,tau,ga,gm)
c
c  The routine calculates the splitting term in the SSS algorithm
c  using the Taylor expansion of the generators through 5th order
c  in t=tau. The expansions have been evaluated using a Mathematica program. 
c  The splitting term has the structure:
c        R(t) exp(:g3(t):) exp(:g4(t):) exp(:g5(t):) exp(:g6(t):),
c  The gn(t) are  given through 5th order in t. R(t) is computed
c  using the routine 'exptay'.
c  Written by M.Venturini April 23, 97.
c  Modified Aug. 4, 97 M.V.

c
c   np = number of squarings
c
       implicit double precision (a-h,o-z)
       dimension h(923),hw(923),hw1(923),ga(923),gm(6,6),hm(6,6)
       dimension em(6,6)
c
       dimension pb23(923),pb2p23(923),pb2p33(923)
       dimension pb24(923),pb2p24(923),pb2p34(923)
       dimension pb25(923),pb2p25(923),pb2p35(923)
c
       dimension pb234(923),pb233(923),pb2p233(923),pb3b233(923)
       dimension pb34(923),pb2p2324(923),pb2p234(923),pb2324(923)
c
c  Map  the content of ga,gm  into ha,hm
         call mapmap(ga,gm,h,hm)
c
c  Reset ga,gm:  
         call clear(ga,gm)
c
c  Evaluate the linear part:
c
        do 10 i=7,27
            ga(i)= -tau*h(i)
  10    continue
c
        call matify(em,ga)
        call exptay(em,gm)
c
c
c   Evaluate some recurrent quantities
c
        tau2=tau*tau
        tau3=tau2*tau
        tau4=tau3*tau
        tau5=tau4*tau
c-----------
c       pb23= :h2: h3
        call pbkt1(h,2,h,3,pb23)
c
c       pb2p23=:h2:^2 h3
        call pbkt1(h,2,pb23,3,pb2p23)
c
c       pb2p33 = :h2:^3 h3
        call pbkt1(h,2,pb2p23,3,pb2p33)
c
c-------
c       pb24= :h2: h4
        call pbkt1(h,2,h,4,pb24)
c
c       pb2p24=:h2:^2 h4
        call pbkt1(h,2,pb24,4,pb2p24)
c
c       pb2p34 = :h2:^3 h4
        call pbkt1(h,2,pb2p24,4,pb2p34)
c--------
c
c       pb25= :h2: h5
        call pbkt1(h,2,h,5,pb25)
c
c       pb2p25=:h2:^2 h5
        call pbkt1(h,2,pb25,5,pb2p25)
c
c       pb2p35 = :h2:^3 h5
        call pbkt1(h,2,pb2p25,5,pb2p35)
c
c--------
c       pb34 = [h3,h4]
        call pbkt1(h,3,h,4,pb34)
c
c       pb234 = [[h2,h3],h4]
        call pbkt1(pb23,3,h,4,pb234)
c
c       pb233 = [[h2,h3],h3]
        call pbkt1(pb23,3,h,3,pb233)

c       pb2p233 = [:h2:^2 h3 ,h3]
        call pbkt1(pb2p23,3,h,3,pb2p233)
c
c       pb3b233 = [h3,[:h2: h3,h3]]
        call pbkt1(h,3,pb233,4,pb3b233)
c
c       pb2p234 = [:h2:^2 h3,h4]
        call pbkt1(pb2p23,3,h,4,pb2p234)
c
c       pb2p2324 = [:h2: h3, :h2: h4]
        call pbkt1(pb23,3,pb24,4,pb2324)
c
c          
c
c***************************
c Evaluate the g3 term (only the direct term is present):
c
        call dirterm(tau,h,3,ga)
c
c**************************
c Evaluate the g4 term:
c direct term
          call dirterm(tau,h,4,ga)
c
c feed-up term
           tc=tau3/12
          call pmadd(pb233,4,tc,ga)
c
           tc=tau4/24
          call pmadd(pb2p233,4,tc,ga)
c
          call pbkt1(pb2p23,3,pb23,3,hw)
           tc=tau5/120 
          call pmadd(hw,4,tc,ga)
c
          call pbkt1(pb2p33,3,h,3,hw)
           tc=tau5/80
          call pmadd(hw,4,tc,ga)
c
c**************************
c Evaluate the g5 term:
c  direct term
          call dirterm(tau,h,5,ga)
c
c  feed-up term
c1
          tc=-tau2/2
          call pmadd(pb34,5,tc,ga)
c2
          call pbkt1(h,3,pb24,4,hw)
          tc=-tau3/3
          call pmadd(hw,5,tc,ga)
c3
          tc=-tau3/6
          call pmadd(pb234,5,tc,ga)
c4
          tc=tau4/24
          call pmadd(pb3b233,5,tc,ga)
c5
          call pbkt1(h,3,pb2p24,4,hw)
          tc=-tau4/8
          call pmadd(hw,5,tc,ga)
c6
          tc=-tau4/8
          call pmadd(pb2324,5,tc,ga)
c7
          tc=-tau4/24
          call pmadd(pb2p234,5,tc,ga)
c8
          call pbkt1(h,3,pb2p233,4,hw)
          tc=tau5/45
          call pmadd(hw,5,tc,ga)
c9
          call pbkt1(h,3,pb2p34,4,hw)
          tc = -tau5/30
          call pmadd(hw,5,tc,ga)
c10
          call pbkt1(pb23,3,pb233,4,hw)
          tc = tau5/60
          call pmadd(hw,5,tc,ga)
c11
          call pbkt1(pb23,3,pb2p24,4,hw)
          tc = -tau5/20
          call pmadd(hw,5,tc,ga)
c12
          call pbkt1(pb2p23,3,pb24,4,hw)
          tc = -tau5/30
          call pmadd(hw,5,tc,ga)
c13
          call pbkt1(pb2p33,3,h,4,hw)
          tc = -tau5/120
          call pmadd(hw,5,tc,ga)
c
c
c**************************
c
c Evaluate the g6 term
c  direct term
          call dirterm(tau,h,6,ga)
c
c  feed-up term
c1
          call pbkt1(h,3,h,5,hw)
          tc = -tau2/2
          call pmadd(hw,6,tc,ga)
c2
          call pbkt1(h,3,pb34,5,hw)
          tc = -tau3/6
          call pmadd(hw,6,tc,ga)
c3
          call pbkt1(h,3,pb25,5,hw)
          tc = -tau3/3
          call pmadd(hw,6,tc,ga)
c4
          call pbkt1(pb23,3,h,5,hw)
          tc = -tau3/6
          call pmadd(hw,6,tc,ga)
c5
          call pbkt1(pb24,4,h,4,hw)
          tc = tau3/12
          call pmadd(hw,6,tc,ga)
c6
          call pbkt1(h,3,pb24,4,hw1)
          call pbkt1(h,3,hw1,5,hw)
          tc = -tau4/8
          call pmadd(hw,6,tc,ga)
c7
          call pbkt1(h,3,pb234,5,hw)
          tc = -tau4/16
          call pmadd(hw,6,tc,ga)
c8
          call pbkt1(h,3,pb2p25,5,hw)
          tc = -tau4/8
          call pmadd(hw,6,tc,ga)
c9
          call pbkt1(h,4,pb233,4,hw)
          tc = tau4/48
          call pmadd(hw,6,tc,ga)
c10
          call pbkt1(pb23,3,pb34,5,hw)
          tc = -tau4/16
          call pmadd(hw,6,tc,ga)
c11
          call pbkt1(pb23,3,pb25,5,hw)
          tc = -tau4/8
          call pmadd(hw,6,tc,ga)
c12
          call pbkt1(pb2p23,3,h,5,hw)
          tc = -tau4/24
          call pmadd(hw,6,tc,ga)
c13
          call pbkt1(pb2p24,4,h,4,hw)
          tc = tau4/24
          call pmadd(hw,6,tc,ga)
c14
          call pbkt1(h,3,pb3b233,5,hw)
          tc = tau5/80
          call pmadd(hw,6,tc,ga)
c15
          call pbkt1(h,3,pb2p24,4,hw1)
          call pbkt1(h,3,hw1,5,hw)
          tc = -tau5/20
          call pmadd(hw,6,tc,ga)
c16
          call pbkt1(h,3,pb2324,5,hw)
          tc = -tau5/20
          call pmadd(hw,6,tc,ga)
c17
          call pbkt1(h,3,pb2p234,5,hw)
          tc = -tau5/60
          call pmadd(hw,6,tc,ga)
c18
          call pbkt1(h,3,pb2p35,5,hw)
          tc = -tau5/30
          call pmadd(hw,6,tc,ga)
c19
          call pbkt1(h,4,pb2p233,4,hw)
          tc = tau5/80
          call pmadd(hw,6,tc,ga)
c20
          call pbkt1(h,3,pb24,4,hw1)
          call pbkt1(pb23,3,hw1,5,hw)
          tc = -tau5/20
          call pmadd(hw,6,tc,ga)
c21
          call pbkt1(pb23,3,pb234,5,hw)
          tc = -tau5/40
          call pmadd(hw,6,tc,ga)
c22
          call pbkt1(pb23,3,pb2p25,5,hw)
          tc = -tau5/20
          call pmadd(hw,6,tc,ga)
c23
          call pbkt1(pb24,4,pb233,4,hw)
          tc = tau5/240
          call pmadd(hw,6,tc,ga)
c24
          call pbkt1(pb2p23,3,pb34,5,hw)
          tc = -tau5/60
          call pmadd(hw,6,tc,ga)
c25
          call pbkt1(pb2p23,3,pb25,5,hw)
          tc = -tau5/30
          call pmadd(hw,6,tc,ga)
c26
          call pbkt1(pb2p24,4,pb24,4,hw)
          tc = tau5/120
          call pmadd(hw,6,tc,ga)
c27
          call pbkt1(pb2p33,3,h,5,hw)
          tc = -tau5/120
          call pmadd(hw,6,tc,ga)
c28
          call pbkt1(pb2p34,4,h,4,hw)
          tc = tau5/80
          call pmadd(hw,6,tc,ga)
c
        return
        end
c
***********************************************************************
c 
        subroutine dirterm(t,h,ideg,ga)
c
c   The routine calculates g =- sum_{i=1}^{nt} (t^i/i!) :h_2:^{i-1} h_{ideg}
c   through t^nt;
c   (direct term in the SSS algorithm).
c   Written by M.Venturini April 23, 97.
c
       implicit double precision (a-h,o-z)
       parameter(nt=5)
       dimension h(923),hw(923),ga(923),hw1(923)
       dimension isup(6),iinf(6)
       data iinf /0,0,28,84,210,462/
       data isup /0,0,83,209,461,923/
       save iinf,isup    !cryne 7/23/2002
c
c Initialize
       do 10 i=iinf(ideg),isup(ideg)
         hw(i)=h(i)
         ga(i)=0.d0
  10   continue
c
       tc=t
       call pmadd(hw,ideg,-tc,ga)
c
       do 20 i=2,nt
            call pbkt1(h,2,hw,ideg,hw1)
            tc=tc*t/float(i)
            call pmadd(hw1,ideg,-tc,ga)
c
                 do 11 j=iinf(ideg),isup(ideg)
                 hw(j)=hw1(j)
  11             continue
c           
  20   continue
c
       return
       end
c
***********************************************************************
c 
        subroutine   expJS(ga,gm)
c
c   The subroutine calculates the matrix exp. out of the 
c   Lie generator gm. For now the sub cex is called. In the 
c   future a direct algorithm should be implemented. 
c   Written by M.Venturini April 24, 97.
c
       implicit double precision (a-h,o-z)
       dimension p(6),gm(6,6),ga(923),ga2(923)
c
       do i=7,923
         ga2(i)=ga(i)
       enddo
c
       do i=28,923
         ga2(i)=0.d0
       enddo
c
       p(1)=1.d0
       p(2)=0
       p(3)=0
c
       call cex(p,ga2,gm)
c
       return
       end
c
***********************************************************************
         
