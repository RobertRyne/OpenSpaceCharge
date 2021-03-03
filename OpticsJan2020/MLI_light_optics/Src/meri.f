************************************************************************
* header                 MERIT FUNCTIONS
*  Weighted sum of sqares and user specified merit functions
***********************************************************************
c
      subroutine mrt0(pp)
c
c   This is the weighted sum of squares or the square root of
c   weighted sum of squares merit function defined by aims for
c   loop pp(2).
c
c     T. Mottershead, LANL  Feb 89; A. Dragt 7/10/98
c-----------------------------------------------------
      include 'impli.inc'
      include 'merit.inc'
c
      dimension pp(*)
c
c        gather the computed function values
c
      kind = nint(pp(1))
      loon = nint(pp(2))
c      write(6,*) 'in mrt0'
      temp = rmserr(loon)
      if (kind .eq. 1) fval = temp**2
      if (kind .eq. 2) fval = temp
      val(0) = fval
c
      return
      end
c
*****************************************************************
c
      subroutine mrt1(p)
c
      include 'impli.inc'
      include 'merit.inc'
      common/count/ifcnt,igcnt,ihcnt
c
c Calling arrays:
      dimension p(6)
c
c Local arrays:
      dimension c(10), a(10), b(10)
c
      write(6,*) 'in mrt1'
      ifcnt=ifcnt+1
c
      fcn = 0.d0
      n=6
      kc = n/2
      do 10 i= 1,n
      c(i) = i-kc
      a(i) = 2 + mod(i,3)
      b(i) = mod(i,kc)
      next = 1 + mod(i,n)
      fcn = fcn + a(i)*(p(i) + b(i)*p(next)- c(i))**2
   10 continue
      val(1) = fcn
      return
      end
c
************************************************************************
c
      subroutine mrt2(p)
c
      include 'impli.inc'
      include 'merit.inc'
      include 'parset.inc'
      include 'usrdat.inc'
c
c Calling arrays:
      dimension p(6)
c
c set up sum of squares merit function based 
c on contents of ucalc
c
c get strengths of octupole components of cfqds     
      res1=ucalc(1)
      res2=ucalc(2)
      res3=ucalc(3)
c get strengths of remaining octupoles
      res4=ucalc(4)
      res5=ucalc(5)
      res6=ucalc(6)
      res7=ucalc(7)
c form sum of squares of cfqds octupole strengths
      sscf=res1**2 + res2**2 + res3**2
c form sum of squares of remaining octupole strengths 
      sso=res4**2 + res5**2 + res6**2 + res7**2
c form full sum of squares
      sst = sscf + sso
c divide by 7 (the number of octupoles)
      wss = sst/7.d0
c square root and square result
      rwssq = sqrt(wss)
      wss = rwssq**2
c assign value
      val(2)=wss
c
      return
      end
c
************************************************************************
c
      subroutine mrt3(p)
c
      include 'impli.inc'
      include 'merit.inc'
c
c Calling arrays:
      dimension p(6)
      write(6,*) 'in mrt3'
      val(3)=p(3)
      return
      end
c
************************************************************************
c
      subroutine mrt4(p)
c
      include 'impli.inc'
      include 'merit.inc'
c
c Calling arrays:
      dimension p(6)
      write(6,*) 'in mrt4'
      val(4)=p(4)
      return
      end
c
************************************************************************
c
      subroutine mrt5(p)
c
c      provides standard set of optimization test problems
c      C. T. Mottershead  LANL AT-3  19 Mar 93
c------------------------------------------------------------
      include 'impli.inc'
      include 'merit.inc'
      include 'parset.inc'
c
c Calling arrays:
      dimension p(6), xv(6)
c      write(6,*) 'in mrt5'
      mode = nint(p(1))
      numf = nint(p(2))
      nvar = nint(p(3))
      npset = nint(p(4))
      if((npset.lt.1).or.(npset.gt.maxpst)) then
      write(6,*) ' *** mrt5 error:'
      write(6,*) npset,' is not a valid parameter set'
      return
      endif
      isend = nint(p(5))
      do 20 j = 1, nvar
      xv(j) = pst(j,npset)
  20  continue
      call optest(mode, numf, nvar, xv, fval)
      val(5) = fval
      return
      end
c
c end of file
      subroutine optest(mode,numf,nv,xv,fv)
c
c    selects and evaluates at (xv(i), i=1,nv) test function number numf
c    from the standard suite from R. Schnabel, U. Colo
c            C. T. Mottershead  AT-3   19 Mar 93
c-------------------------------------------------
      implicit double precision (a-h,o-z)
      character*24 remark(20), title
      character*6  fcname(20)
      dimension maxvar(20), xv(*)
      external fcn03,fcn04,fcn05,fcn07,fcn09,fcn12,fcn14
     *, fcn16,fcn18,fcn20,fcn21,fcn22,fcn23,fcn24,fcn25
     *, fcn26,fcn35,fcn36,fcn37,fcn38
      data maxvar /3*2,3*3,2*4,11*6,4/,  ndone /0/
      if(ndone.eq.0) call fcnlbl(num,remark,fcname)
      ndone = 1
c
c     list available test functions
c
      if(mode.eq.1) then
      write(6,*)  num,' test functions available'
      kpr = (num + 1)/2
      do 20 k=1,kpr
      k2 = k+kpr
      write(6,27) k, fcname(k), remark(k), k2, fcname(k2), remark(k2)
  27  format(2(i4,'= ',a6,': ',a24))
  20  continue
      go to 80
      endif
c
c         evaluate selected function
c
      if((numf.le.0).or.(numf.gt.20)) go to 77
      title = remark(numf)
      if(nv.gt.maxvar(numf)) nv=maxvar(numf)
      write(6,51) nv,numf,title
  51  format(i5,' variables in test function',i3,': ',a)
      go to (61,62,63,64,65,66,67,68,69,70,71,72,73,24,25,26,35,36
     * ,37,38), numf
  61  call fcn03(nv, xv, fv)
      go to 80
  62  call fcn04(nv, xv, fv)
      go to 80
  63  call fcn05(nv, xv, fv)
      go to 80
  64  call fcn07(nv, xv, fv)
      go to 80
  65  call fcn09(nv, xv, fv)
      go to 80
  66  call fcn12(nv, xv, fv)
      go to 80
  67  call fcn14(nv, xv, fv)
      go to 80
  68  call fcn16(nv, xv, fv)
      go to 80
  69  call fcn18(nv, xv, fv)
      go to 80
  70  call fcn20(nv, xv, fv)
      go to 80
  71  call fcn21(nv, xv, fv)
      go to 80
  72  call fcn22(nv, xv, fv)
      go to 80
  73  call fcn23(nv, xv, fv)
      go to 80
  24  call fcn24(nv, xv, fv)
      go to 80
  25  call fcn25(nv, xv, fv)
      go to 80
  26  call fcn26(nv, xv, fv)
      go to 80
  35  call fcn35(nv, xv, fv)
      go to 80
  36  call fcn36(nv, xv, fv)
      go to 80
  37  call fcn37(nv, xv, fv)
      go to 80
  38  call fcn38(nv, xv, fv)
  80  return
  77  write(6,*) numf, ' = invalid function number'
      return
      end
c
cccccccccccccccccccccccccc   fcnlbl  ccccccccccccccccccccccccccccccccccc
c
      subroutine fcnlbl(num,remark,fcname)
      character*24 remark(*)
      character*6  fcname(*)
      num = 0
c
      num = num + 1
      remark(num) =  'POWELLS BADLY SCALED'
      fcname(num) = 'FCN03'
c
      num = num + 1
      remark(num) =  'BROWN BADLY SCALED'
      fcname(num) = 'FCN04'
c
      num = num + 1
      remark(num) =  'BEALE'
      fcname(num) = 'FCN05'
c
      num = num + 1
      remark(num) =  'HELICAL VALLEY'
      fcname(num) = 'FCN07'
c
      num = num + 1
      remark(num) =  'GAUSSIAN'
      fcname(num) = 'FCN09'
c
      num = num + 1
      remark(num) =  'BOX 3D'
      fcname(num) = 'FCN12'
c
      num = num + 1
      remark(num) =  'WOOD'
      fcname(num) = 'FCN14'
c
      num = num + 1
      remark(num) =  'BROWN-DENNIS'
      fcname(num) = 'FCN16'
c
      num = num + 1
      remark(num) =  'BIGGS EXP'
      fcname(num) = 'FCN18'
c
      num = num + 1
      remark(num) =  'WATSON'
      fcname(num) = 'FCN20'
c
      num = num + 1
      remark(num) =  'ROSENBROCK BANANA VALLEY'
      fcname(num) = 'FCN21'
c
      num = num + 1
      remark(num) =  'POWELL SINGULAR'
      fcname(num) = 'FCN22'
c
      num = num + 1
      remark(num) =  'PENALTY I'
      fcname(num) = 'FCN23'
c
      num = num + 1
      remark(num) =  'PENALTY II'
      fcname(num) = 'FCN24'
c
      num = num + 1
      remark(num) =  'VARIABLE DIMENSION'
      fcname(num) = 'FCN25'
c
      num = num + 1
      remark(num) =  'TRIGONOMETRIC'
      fcname(num) = 'FCN26'
c
      num = num + 1
      remark(num) =  'CHEBYQUAD'
      fcname(num) = 'FCN35'
c
      num = num + 1
      remark(num) =  'ROSENBROCK BADLY SCALED'
      fcname(num) = 'FCN36'
c
      num = num + 1
      remark(num) =  'QUADRATIC FORM'
      fcname(num) = 'FCN37'
c
      num = num + 1
      remark(num) =  'BERZ'
      fcname(num) = 'FCN38'
      return
      end
c
cccccccccccccccccccccccc powell bad scale cccccccccccccccccccccccccc
c
      subroutine fcn03(n,x,f)
      implicit double precision (a-h, o-z)
c
c powell"s badly scaled function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      f=(1.d4*x(1)*x(2)-1.d0)**2 + (exp(-x(1))+exp(-x(2))-1.0001d0)**2
      return
      end
c
ccccccccccccccccccccccc  brown bad scale  ccccccccccccccccccc
c
      subroutine fcn04(n,x,f)
      implicit double precision (a-h, o-z)
c
c brown badly scaled function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      f=(x(1)-1.d6)**2 + (x(2)-2.d-6)**2 + (x(1)*x(2)-2.d0)**2
      return
      end
c
cccccccccccccccccccccccccccccccc  beale  ccccccccccccccccccccccccc
c
      subroutine fcn05(n,x,f)
      implicit double precision (a-h, o-z)
c
c beale function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      f=(1.5d0-1.0d0*x(1)*(1.d0-1.0d0*x(2)))**2
     +      + (2.25d0 -1.0d0*x(1)*(1.d0-1.0d0*x(2)*x(2)))**2
     +      + (2.625d0-1.0d0*x(1)*(1.d0-1.d0*x(2)*x(2)*x(2)))**2
      return
      end
c
ccccccccccccccccccccccccc  helical valley  cccccccccccccccccccccc
c
      subroutine fcn07(n,x,f)
      implicit double precision (a-h, o-z)
c
c helical valley function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
c      pi=3.1415926535898
      pi=3.1415926535897931d0
      if(x(1).gt. 0.d0) th= (1.d0/(2.d0*pi)) * atan(x(2)/x(1))
      if(x(1).lt. 0.d0) th= (1.d0/(2.d0*pi)) * atan(x(2)/x(1))+0.5d0
      f=100.d0*(x(3)-10.d0*th)**2
     +      +100.d0*(sqrt(x(1)**2 + x(2)**2)-1.d0)**2
     +      + x(3)*x(3)
      return
      end
c
cccccccccccccccccccccc    odd gaussian  cccccccccccccccccccccccccc
c
      subroutine fcn09(n,x,f)
      implicit double precision (a-h, o-z)
c
c gaussian function
c
      dimension x(n),y(15)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      y(1)=.0009d0
      y(2)=.0044d0
      y(3)=.0175d0
      y(4)=.0540d0
      y(5)=.1295d0
      y(6)=.2420d0
      y(7)=.3521d0
      y(8)=.3989d0
      y(9)= y(7)
      y(10)=y(6)
      y(11)=y(5)
      y(12)=y(4)
      y(13)=y(3)
      y(14)=y(2)
      y(15)=y(1)
      g=0.d0
      do 10 i=1,15
        g=g + (x(1)*exp( -x(2)*((8.d0-i)/2.d0-x(3))**2/2.d0) - y(i))**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccccc  3D box  ccccccccccccccccccccccc
c
      subroutine fcn12(n,x,f)
      implicit double precision (a-h, o-z)
c
c box 3-dimensional function(n,x,f)
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 10 i=1,n
        t=i/10.d0
        g=g + (exp(-t*x(1)) - exp(-t*x(2))
     +        - x(3)*(exp(-t)-exp(-10.d0*t)))**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccccc wood  cccccccccccccccccccccccccccc
c
      subroutine fcn14(n,x,f)
      implicit double precision (a-h, o-z)
c
c wood function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      f=100.d0*(x(2)-x(1)*x(1))**2 + (1.d0-x(1))**2
     +      + 90.d0*(x(4)-x(3)*x(3))**2 + (1.d0-x(3))**2
     +      +10.1d0*( (1.d0-x(2))**2 + (1.d0-x(4))**2)
     +      +19.8d0*(1.d0-x(2))*(1.d0-x(4))
      return
      end
c
cccccccccccccccccccccccccc   brown-dennis  cccccccccccccccccccccccccccc
c
      subroutine fcn16(n,x,f)
      implicit double precision (a-h, o-z)
c
c brown + dennis function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 10 i=1,20
        t=i/5.d0
        g=g + ((x(1)+t*x(2)-exp(t))**2
     +      +  (x(3)+x(4)*sin(t)-cos(t))**2 )**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccc   biggs exp6   cccccccccccccccccccccccc
c
      subroutine fcn18(n,x,f)
      implicit double precision (a-h, o-z)
c
c biggs exp6 function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 10 i=1,13
        t=i/10.d0
        g=g + (x(3)*exp(-t*x(1))-x(4)*exp(-t*x(2))
     +        + x(6)*exp(-t*x(5))-exp(-t)+5.d0*exp(-10.d0*t)
     +        - 3.d0*exp(-4.d0*t) )**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccc   watson   cccccccccccccccccccccccccccccc
c
      subroutine fcn20(n,x,f)
      implicit double precision (a-h, o-z)
c
c watson function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 40 i=1,29
        t=i/29.d0
        sum=0.d0
        s=1.d0
        do 20 j=2,n
          sum=sum + (j-1)*x(j)*s
          s=t*s
   20   continue
        sum2=0.d0
        s=1.d0
        do 30 j=1,n
          sum2=sum2 + x(j)*s
          s=t*s
   30   continue
        g=g + (sum - sum2*sum2 - 1.d0)**2
   40 continue
      g=g + x(1)*x(1) + (x(2)-x(1)**2-1.d0)**2
      f=g
      return
      end
c
ccccccccccccccccccccccc   extended rosenbrock  cccccccccccccccccccccc
c
      subroutine fcn21(n,x,f)
      implicit double precision (a-h, o-z)
c
c extended rosenbrock function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      j=n/2
      do 10 i=1,j
        g=g+ 100.d0*(x(2*i)- x(2*i-1)**2)**2 + (1.d0-x(2*i-1))**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccc   powell singular  cccccccccccccccccccccccc
c
      subroutine fcn22(n,x,f)
      implicit double precision (a-h, o-z)
c
c extended powell singular function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      j=n/4
      do 10 i=1,j
        g=g + (x(4*i-3) + 10.d0*x(4*i-2))**2
     +      + 5.d0*(x(4*i-1) - x(4*i))**2
     +      + ( (x(4*i-2)-2.d0*x(4*i-1))**2)**2
     +      + 10.d0*( (x(4*i-3)-x(4*i))**2)**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccccccccc   penalty I ccccccccccccccccccc
c
      subroutine fcn23(n,x,f)
      implicit double precision (a-h, o-z)
c
c penalty function i
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 20 i=1,n
        g=g + (1.d-5) * (x(i)-1.d0)**2
   20 continue
      sum=0.d0
      do 30 j=1,n
        sum=sum + x(j)*x(j)
   30 continue
      g=g + (sum-.25d0)**2
      f=g
      return
      end
c
ccccccccccccccccccccccccccccccc   penalty II  cccccccccccccccccccc
c
      subroutine fcn24(n,x,f)
      implicit double precision (a-h, o-z)
c
c penalty function ii
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=(x(1)-.2d0)**2
      do 10 i=2,n
        t=exp(i/10.d0) + exp((i-1)/10.d0)
        g=g + 1.d-5* (exp(x(i)/10.d0) + exp(x(i-1)/10.d0) - t)**2
   10 continue
      m1=n+1
      m2=2*n-1
      do 20 i=m1,m2
        g=g + 1.d-5* (exp(x(i-n+1)/10.d0) - exp(1.d0/10.d0))**2
   20 continue
      sum=0.d0
      do 30 j=1,n
        sum=sum + (n-j+1)*x(j)*x(j)
   30 continue
      g=g + (sum-1.d0)**2
      f=g
      return
      end
c
ccccccccccccccccccccccccccccc   variable dimension  ccccccccccccccccc
c
      subroutine fcn25(n,x,f)
      implicit double precision (a-h, o-z)
c
c variably dimensioned function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      do 10 i=1,n
        g=g + (x(i)-1.d0)**2
   10 continue
      s=0.d0
      do 20 i=1,n
        s=s + i*(x(i)-1.d0)
   20 continue
      g=g + s*s + (s*s)**2
      f=g
      return
      end
c
cccccccccccccccccccccc   trigonometric   cccccccccccccccccccccc
c
      subroutine fcn26(n,x,f)
      implicit double precision (a-h, o-z)
c
c trigonometric function
c
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      sum=0.d0
      do 10 j=1,n
        sum=sum + cos(x(j))
   10 continue
      do 20 i=1,n
        g=g+ (n-sum + i*(1.d0-cos(x(i))) - sin(x(i)) )**2
   20 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccc   chebyquad  cccccccccccccccccccccccccccccccccc
c
      subroutine fcn35(n,x,f)
      implicit double precision (a-h, o-z)
c
c chebyquad function
c
      dimension x(n)
c work arrays passed thru common dimensioned .ge. n
      common y1(100),y2(100),y3(100)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
c
      sum=0.d0
      do 10 j=1,n
        y1(j)=2.d0*x(j)-1.d0
        sum=sum+y1(j)
   10 continue
      g=(sum/n)**2
c
      sum=0.d0
      do 20 j=1,n
        y2(j)=2.d0*(2.d0*x(j)-1.d0)*y1(j)-1.d0
        sum=sum+y2(j)
   20 continue
      g=g+ (sum/n + 1.d0/3.d0)**2
c
      sum=0.d0
      do 30 j=1,n
        y3(j)=2.d0*(2.d0*x(j)-1.d0)*y2(j)-y1(j)
        sum=sum+y3(j)
   30 continue
      g=g+(sum/n)**2
c
      do 40 i=4,n
c i-th component
        sum=0.d0
        do 35 j=1,n
          y1(j)=y2(j)
          y2(j)=y3(j)
          y3(j)=2.d0*(2.d0*x(j)-1.d0)*y3(j)-y1(j)
          sum=sum+y3(j)
   35   continue
        k=mod(i,2)
        if(k.eq.1) g=g+(sum/n)**2
        if(k.eq.0) g=g+(sum/n + 1.d0/(i*i-1))**2
   40 continue
      f=g
      return
      end
c
ccccccccccccccccccccccc   extended rosenbrock badscale  ccccccccccccccccccc
c
      subroutine fcn36(n,x,f)
      implicit double precision (a-h,o-z)
c
c extended rosenbrock badly scaled function
c
      dimension x(n)
      parameter (alpha = 1.0d-2)
      dimension xt(10)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
      g=0.d0
      j=n/2
      do 10 i=1,j
        xt(2*i-1) = x(2*i-1)*alpha
        xt(2*i) = x(2*i)/alpha
        g=g+ 100.d0*(xt(2*i)- xt(2*i-1)**2)**2 + (1.d0-xt(2*i-1))**2
   10 continue
      f=g
      return
      end
c
ccccccccccccccccccccccccccccccc  quad  cccccccccccccccccccccccccccc
c
      subroutine fcn37(n,x,fv)
c
c    rotated quadratic form with minimum of n at x(k) = k, k=1,n
c    C. T. Mottershead  LANL AT3 22 Mar 93
c------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(n), c(10), a(10), b(10)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
c
      fv = n
      do 10 k = 1,n
      ak = k
      km = k-1
      if(km.lt.1) km = n
      akm = km
      kp = k+1
      if(kp.gt.n) kp = 1
      akp = kp
      fv = fv + ak*(x(k)-ak + x(kp)-akp - x(km)+akm)**2
   10 continue
      return
      end
c
ccccccccccccccccccccccccccccccc  berz  cccccccccccccccccccccccccccc
c
      subroutine fcn38(n,x,f)
      implicit double precision (a-h,o-z)
      dimension x(n)
      common/count/ifcnt,igcnt,ihcnt
      ifcnt=ifcnt+1
c
      f = x(1)**4 + x(2)**2 + dabs( x(2)**3 + x(4)**3) +
     $    13.0d0 * x(3)**2 + (x(2) - x(3)**2)**2
      return
      end
