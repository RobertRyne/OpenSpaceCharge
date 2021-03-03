c*******************************************************************
c  ImKmpak is a collection of portable modified Bessel function routines.
c  Written 4-96 by P. L. Walstrom Northrop Grumman
c  Modified from CLAMS and Numerical Recipes.  Should give 12-digit
c  precision for all valid arguments.
c  Should use exponentially scaled functions when x is greater than 20.
c  Function routines:
c  dbesI0s=I0
c  dbesI0es=exponentially scaled I0
c  dbesK0s=K0
c  dbesK0es=exponentially scaled K0
c  dbesI1s=I1
c  dbesI1es=exponentially scaled I1
c  dbesK1s=K1
c  dbesK1es=exponentially scaled K1
c  Subroutines:
c
c  BESSIn( M,N,x,y,KODE)
c      computes a M-member sequence I_k(x), k=N, N+1, N+2 ... N+M-1.
c   KODE=1 means unscaled, KODE=2 is exponentially scaled. Results in vector
c   y, with y(1)=I_N, y(2)=I_N+1, etc.
c
c
c  bessKn( M,N,x,y,KODE)
c      computes a M-member sequence K_k(x), k=N, N+1, N+2 ... N+M-1.
c   KODE=1 means unscaled, KODE=2 is exponentially scaled. Results in vector
c   y, with y(1)=K_N, y(2)=K_N+1, etc.
c
c
c  Exponential scaling means for the I_m
c
c  I_m (x) scaled = exp(-x) * I_m (x) unscaled
c
c  and for the K_m
c
c  K_m (x) scaled = exp(x) * K_m (x) unscaled
c
c*******************************************************************
      subroutine BESSIn(M,N,X,y,KODE)
        implicit double precision(a-h,o-z)
c Calculates an M-member sequence of modified Bessel functions of the
c  first kind.  Needs subroutines for exponentially scaled and unscaled
c  I_0(x).  Sequence contains either unscaled or scaled Bessel functions.
c  KODE=1 means unscaled.
c  KODE=2 means scaled.
c  y(k), k=1,2,..m contains I_N (x),I_N+1 (x),...I_N+M-1 (x)
c  N must be > or = 0.
      PARAMETER(IACC=40,BIGNO=1.0d10,BIGNI=1.0d-10)
        parameter(mmax=1000)
        dimension y(m),z(mmax)
C
        if(m.gt.20) go to 1005
        if(dabs(x).gt.0.1d0) go to 44
        if(dabs(x).ne.0.d0) go to 23
        do 24 l=1,m
        if(n.eq.0) y(1)=1.d0
   24 y(l)=0.d0
        return
   23 continue
c  Series expansion for I_n, I_n+1, ...I_n+m-1.
        if(m.gt.mmax) go to 1001
        qx2=0.25d0*x**2
        hfx=0.5d0*x
c  Find 1/n! and (x/2)**n.
        rnfac=1.d0
        hfxn=1.d0
        if(n.gt.0) hfxn=hfx
        if(n.lt.2) go to 1
        do 2 l=2,n
        hfxn=hfxn*hfx
    2 rnfac=rnfac/dfloat(l)
    1 y(1)=rnfac
c  find 1/ (n+l)!, l=0,1,2...m-1.
c  store in z(1),z(2)...z(m)
        z(1)=rnfac
        if(m.lt.2) go to 6
        do 3 l=2,m
        z(l)=z(l-1)/dfloat(l+n-1)
    3 y(l)=z(l)
    6 yy=1.d0
        do 5 k=1,4
        yy=yy*qx2/dfloat(k)
        do 5 l=1,m
        z(l)=z(l)/dfloat(l+k+n-1)
    5 y(l)=y(l)+yy*z(l)
        y(1)=y(1)*hfxn
        if(m.lt.2) go to 9
        yy=hfxn
        do 10 l=2,m
        yy=yy*hfx
   10 y(l)=y(l)*yy
    9 continue
        if(kode.eq.1) return
c  convert vector of I_n, I_n+1...I_n+m-1 to exponentially scaled values
        zz=dexp(-dabs(x))
        do 45 l=1,m
   45 y(l)=zz*y(l)
        return
c  Larger x- use downwards recursion
   44 n2=n+m-1
      IF (N.LT.0) PAUSE 'bad argument N < 0 in BESSIn'
        if(dabs(x).gt.dfloat(2*n2)) go to 60
      TOX=2.d0/X
      BIP=0.d0
      BI=1.d0
        do 21 k=1,m
   21 y(k)=0.d0
      MAX=2*((N2+dINT(dSQRT(dFLOAT(IACC*N2)))))
      DO 11 J=MAX,1,-1
        BIM=BIP+dFLOAT(J)*TOX*BI
        BIP=BI
        BI=BIM
c  Rescale large numbers
        IF (dABS(BI).GT.BIGNO) THEN
          BI=BI*BIGNI
          BIP=BIP*BIGNI
        jmax=N2
        jmin=j+1
        if(jmin.gt.jmax) go to 22
        jdif=n2-j
        if(jdif.gt.m) go to 17
        do 13 jj=jmin,jmax
        k=jj-jmax+m
   13 y(k)=y(k)*bigni
        go to 22
   17 do 18 k=1,m
   18 y(k)=y(k)*bigni
   22 continue
        ENDIF
c  Rescale small numbers.
        IF (dABS(BI).lT.BIGNI) THEN
          BI=BI*BIGNO
          BIP=BIP*BIGNO
        jmax=N2
        jmin=j+1
        if(jmin.gt.jmax) go to 25
        jdif=n2-j
        if(jdif.gt.m) go to 19
        do 16 jj=jmin,jmax
        k=jj-jmax+m
   16 y(k)=y(k)*bigno
        go to 25
   19 do 20 k=1,m
   20 y(k)=y(k)*bigno
   25 continue
        ENDIF
        do 14 k=1,m
        ncheck=k+n2-m
   14 if(j.eq.ncheck) y(k)=bip
11    CONTINUE
        if(KODE.eq.1) zz=dbesi0s(x)
        if(KODE.eq.2) zz=dbsi0es(x)
        do 15 k=1,m
   15 y(k)=y(k)*zz/BI
        if(n.eq.0) y(1)=zz
      RETURN
c  Very large x- use upwards recursion
   60 continue
        if(kode.eq.2) bim=dbsi0es(x)
        if(kode.eq.2) bi=dbsi1es(x)
        if(kode.eq.1) bim=dbesi0s(x)
        if(kode.eq.1) bi=dbesi1s(x)
        if(n.gt.1) go to 64
        if(n.gt.0) go to 63
c  n=0
        y(1)=bim
        if(m.lt.2) return
        y(2)=bi
        if(m.lt.3) return
        lmin=3
        go to 64
   63 y(1)=bi
        if(m.lt.2) return
        lmin=2
   64 n2=n+m-1
        if(n.gt.1) lmin=1
        tox=2.d0/x
        do 61  j=1,n2-1
        bip=bim-dfloat(j)*tox*bi
        bim=bi
        bi=bip
        do 62 l=lmin,m
        jj=n+l-2
   62 if(jj.eq.j) y(l)=bi
   61 continue
        return
 1001 write(5,200) m,mmax
  200 format(1x,'m=',i4,' > mmax=',i4,' in BESSIN-stopped')
        stop
 1005 write(6,201) m
  201 format(1x,'m in BESSIn= ',i2,' is > 20- stopped')
        stop
      END
      DOUBLE PRECISION FUNCTION DBSI0Es(X)
c  Stripped of extraneous calls, etc. to make it portable.
C***BEGIN PROLOGUE  DBSI0E
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI0E-S DBSI0E-D),BESSEL FUNCTION,
C             EXPONENTIALLY SCALED,FIRST KIND,
C             HYPERBOLIC BESSEL FUNCTION,MODIFIED BESSEL FUNCTION,
C             ORDER ZERO,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. exponentially scaled hyperbolic Bessel
C            function of the first kind of order zero.
C***DESCRIPTION
C
C DBSI0E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C zero for double precision argument X.  The result is the Bessel
C function I0(X) multiplied by EXP(-ABS(X)).
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
C Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.74E-32
C                                         log weighted error  31.56
C                               significant figures required  30.15
C                                    decimal places required  32.39
C
C Series for AI02       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.97E-32
C                                         log weighted error  31.71
C                               significant figures required  30.15
C                                    decimal places required  32.63
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCSEVLs,INITDS
C***END PROLOGUE  DBSI0E
      DOUBLE PRECISION X, BI0CS(18), AI0CS(46), AI02CS(69),
     1   Y,  DCSEVLs
      SAVE BI0 CS, AI0 CS, AI02CS, NTI0, NTAI0, NTAI02
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA AI0 CS(  1) / +.7575994494 0237959427 2987203743 8 D-1      /
      DATA AI0 CS(  2) / +.7591380810 8233455072 9297873320 4 D-2      /
      DATA AI0 CS(  3) / +.4153131338 9237505018 6319749138 2 D-3      /
      DATA AI0 CS(  4) / +.1070076463 4390730735 8242970217 0 D-4      /
      DATA AI0 CS(  5) / -.7901179979 2128946607 5031948573 0 D-5      /
      DATA AI0 CS(  6) / -.7826143501 4387522697 8898980690 9 D-6      /
      DATA AI0 CS(  7) / +.2783849942 9488708063 8118538985 7 D-6      /
      DATA AI0 CS(  8) / +.8252472600 6120271919 6682913319 8 D-8      /
      DATA AI0 CS(  9) / -.1204463945 5201991790 5496089110 3 D-7      /
      DATA AI0 CS( 10) / +.1559648598 5060764436 1228752792 8 D-8      /
      DATA AI0 CS( 11) / +.2292556367 1033165434 7725480285 7 D-9      /
      DATA AI0 CS( 12) / -.1191622884 2790646036 7777423447 8 D-9      /
      DATA AI0 CS( 13) / +.1757854916 0324098302 1833124774 3 D-10     /
      DATA AI0 CS( 14) / +.1128224463 2189005171 4441135682 4 D-11     /
      DATA AI0 CS( 15) / -.1146848625 9272988777 2963387698 2 D-11     /
      DATA AI0 CS( 16) / +.2715592054 8036628726 4365192160 6 D-12     /
      DATA AI0 CS( 17) / -.2415874666 5626878384 4247572028 1 D-13     /
      DATA AI0 CS( 18) / -.6084469888 2551250646 0609963922 4 D-14     /
      DATA AI0 CS( 19) / +.3145705077 1754772937 0836026730 3 D-14     /
      DATA AI0 CS( 20) / -.7172212924 8711877179 6217505917 6 D-15     /
      DATA AI0 CS( 21) / +.7874493403 4541033960 8390960332 7 D-16     /
      DATA AI0 CS( 22) / +.1004802753 0094624023 4524457183 9 D-16     /
      DATA AI0 CS( 23) / -.7566895365 3505348534 2843588881 0 D-17     /
      DATA AI0 CS( 24) / +.2150380106 8761198878 1205128784 5 D-17     /
      DATA AI0 CS( 25) / -.3754858341 8308744291 5158445260 8 D-18     /
      DATA AI0 CS( 26) / +.2354065842 2269925769 0075710532 2 D-19     /
      DATA AI0 CS( 27) / +.1114667612 0479285302 2637335511 0 D-19     /
      DATA AI0 CS( 28) / -.5398891884 3969903786 9677932270 9 D-20     /
      DATA AI0 CS( 29) / +.1439598792 2407526770 4285840452 2 D-20     /
      DATA AI0 CS( 30) / -.2591916360 1110934064 6081840196 2 D-21     /
      DATA AI0 CS( 31) / +.2238133183 9985839074 3409229824 0 D-22     /
      DATA AI0 CS( 32) / +.5250672575 3647711727 7221683199 9 D-23     /
      DATA AI0 CS( 33) / -.3249904138 5332307841 7343228586 6 D-23     /
      DATA AI0 CS( 34) / +.9924214103 2050379278 5728471040 0 D-24     /
      DATA AI0 CS( 35) / -.2164992254 2446695231 4655429973 3 D-24     /
      DATA AI0 CS( 36) / +.3233609471 9435940839 7333299199 9 D-25     /
      DATA AI0 CS( 37) / -.1184620207 3967424898 2473386666 6 D-26     /
      DATA AI0 CS( 38) / -.1281671853 9504986505 4833868799 9 D-26     /
      DATA AI0 CS( 39) / +.5827015182 2793905116 0556885333 3 D-27     /
      DATA AI0 CS( 40) / -.1668222326 0261097193 6450150399 9 D-27     /
      DATA AI0 CS( 41) / +.3625309510 5415699757 0068480000 0 D-28     /
      DATA AI0 CS( 42) / -.5733627999 0557135899 4595839999 9 D-29     /
      DATA AI0 CS( 43) / +.3736796722 0630982296 4258133333 3 D-30     /
      DATA AI0 CS( 44) / +.1602073983 1568519633 6551253333 3 D-30     /
      DATA AI0 CS( 45) / -.8700424864 0572298845 2249599999 9 D-31     /
      DATA AI0 CS( 46) / +.2741320937 9374811456 0341333333 3 D-31     /
      DATA AI02CS(  1) / +.5449041101 4108831607 8960962268 0 D-1      /
      DATA AI02CS(  2) / +.3369116478 2556940898 9785662979 9 D-2      /
      DATA AI02CS(  3) / +.6889758346 9168239842 6263914301 1 D-4      /
      DATA AI02CS(  4) / +.2891370520 8347564829 6692402323 2 D-5      /
      DATA AI02CS(  5) / +.2048918589 4690637418 2760534093 1 D-6      /
      DATA AI02CS(  6) / +.2266668990 4981780645 9327743136 1 D-7      /
      DATA AI02CS(  7) / +.3396232025 7083863451 5084396952 3 D-8      /
      DATA AI02CS(  8) / +.4940602388 2249695891 0482449783 5 D-9      /
      DATA AI02CS(  9) / +.1188914710 7846438342 4084525196 3 D-10     /
      DATA AI02CS( 10) / -.3149916527 9632413645 3864862961 9 D-10     /
      DATA AI02CS( 11) / -.1321581184 0447713118 7540739926 7 D-10     /
      DATA AI02CS( 12) / -.1794178531 5068061177 7943574026 9 D-11     /
      DATA AI02CS( 13) / +.7180124451 3836662336 7106429346 9 D-12     /
      DATA AI02CS( 14) / +.3852778382 7421427011 4089801777 6 D-12     /
      DATA AI02CS( 15) / +.1540086217 5214098269 1325823339 7 D-13     /
      DATA AI02CS( 16) / -.4150569347 2872220866 2689972015 6 D-13     /
      DATA AI02CS( 17) / -.9554846698 8283076487 0214494312 5 D-14     /
      DATA AI02CS( 18) / +.3811680669 3526224207 4605535511 8 D-14     /
      DATA AI02CS( 19) / +.1772560133 0565263836 0493266675 8 D-14     /
      DATA AI02CS( 20) / -.3425485619 6772191346 1924790328 2 D-15     /
      DATA AI02CS( 21) / -.2827623980 5165834849 4205593759 4 D-15     /
      DATA AI02CS( 22) / +.3461222867 6974610930 9706250813 4 D-16     /
      DATA AI02CS( 23) / +.4465621420 2967599990 1042054284 3 D-16     /
      DATA AI02CS( 24) / -.4830504485 9441820712 5525403795 4 D-17     /
      DATA AI02CS( 25) / -.7233180487 8747539545 6227240924 5 D-17     /
      DATA AI02CS( 26) / +.9921475412 1736985988 8046093981 0 D-18     /
      DATA AI02CS( 27) / +.1193650890 8459820855 0439949924 2 D-17     /
      DATA AI02CS( 28) / -.2488709837 1508072357 2054491660 2 D-18     /
      DATA AI02CS( 29) / -.1938426454 1609059289 8469781132 6 D-18     /
      DATA AI02CS( 30) / +.6444656697 3734438687 8301949394 9 D-19     /
      DATA AI02CS( 31) / +.2886051596 2892243264 8171383073 4 D-19     /
      DATA AI02CS( 32) / -.1601954907 1749718070 6167156200 7 D-19     /
      DATA AI02CS( 33) / -.3270815010 5923147208 9193567485 9 D-20     /
      DATA AI02CS( 34) / +.3686932283 8264091811 4600723939 3 D-20     /
      DATA AI02CS( 35) / +.1268297648 0309501530 1359529710 9 D-22     /
      DATA AI02CS( 36) / -.7549825019 3772739076 9636664410 1 D-21     /
      DATA AI02CS( 37) / +.1502133571 3778353496 3712789053 4 D-21     /
      DATA AI02CS( 38) / +.1265195883 5096485349 3208799248 3 D-21     /
      DATA AI02CS( 39) / -.6100998370 0836807086 2940891600 2 D-22     /
      DATA AI02CS( 40) / -.1268809629 2601282643 6872095924 2 D-22     /
      DATA AI02CS( 41) / +.1661016099 8907414578 4038487490 5 D-22     /
      DATA AI02CS( 42) / -.1585194335 7658855793 7970504881 4 D-23     /
      DATA AI02CS( 43) / -.3302645405 9682178009 5381766755 6 D-23     /
      DATA AI02CS( 44) / +.1313580902 8392397817 4039623117 4 D-23     /
      DATA AI02CS( 45) / +.3689040246 6711567933 1425637280 4 D-24     /
      DATA AI02CS( 46) / -.4210141910 4616891492 1978247249 9 D-24     /
      DATA AI02CS( 47) / +.4791954591 0828657806 3171401373 0 D-25     /
      DATA AI02CS( 48) / +.8459470390 2218217952 9971707412 4 D-25     /
      DATA AI02CS( 49) / -.4039800940 8728324931 4607937181 0 D-25     /
      DATA AI02CS( 50) / -.6434714653 6504313473 0100850469 5 D-26     /
      DATA AI02CS( 51) / +.1225743398 8756659903 4464736990 5 D-25     /
      DATA AI02CS( 52) / -.2934391316 0257089231 9879821175 4 D-26     /
      DATA AI02CS( 53) / -.1961311309 1949829262 0371205728 9 D-26     /
      DATA AI02CS( 54) / +.1503520374 8221934241 6229900309 8 D-26     /
      DATA AI02CS( 55) / -.9588720515 7448265520 3386388206 9 D-28     /
      DATA AI02CS( 56) / -.3483339380 8170454863 9441108511 4 D-27     /
      DATA AI02CS( 57) / +.1690903610 2630436730 6244960725 6 D-27     /
      DATA AI02CS( 58) / +.1982866538 7356030438 9400115718 8 D-28     /
      DATA AI02CS( 59) / -.5317498081 4918162145 7583002528 4 D-28     /
      DATA AI02CS( 60) / +.1803306629 8883929462 3501450390 1 D-28     /
      DATA AI02CS( 61) / +.6213093341 4548931758 8405311242 2 D-29     /
      DATA AI02CS( 62) / -.7692189292 7721618632 0072806673 0 D-29     /
      DATA AI02CS( 63) / +.1858252826 1117025426 2556016596 3 D-29     /
      DATA AI02CS( 64) / +.1237585142 2813957248 9927154554 1 D-29     /
      DATA AI02CS( 65) / -.1102259120 4092238032 1779478779 2 D-29     /
      DATA AI02CS( 66) / +.1886287118 0397044900 7787447943 1 D-30     /
      DATA AI02CS( 67) / +.2160196872 2436589131 4903141406 0 D-30     /
      DATA AI02CS( 68) / -.1605454124 9197432005 8446594965 5 D-30     /
      DATA AI02CS( 69) / +.1965352984 5942906039 3884807331 8 D-31     /
      DATA NTI0, NTAI0, NTAI02 / 3*0/
C***FIRST EXECUTABLE STATEMENT  DBSI0E
      IF (NTI0.NE.0) GO TO 10
c      ETA = 0.1*SNGL(D1MACH(3))
c  replace above statement
        eta=1.4d-18
      NTI0 = INITDSs (BI0CS, 18, ETA)
      NTAI0 = INITDSs (AI0CS, 46, ETA)
      NTAI02 = INITDSs (AI02CS, 69, ETA)
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBSI0Es = 1.0D0
      IF (Y.GT.1.d-15) DBSI0Es = DEXP(-Y) * (2.75D0 +
     1  DCSEVLs (Y*Y/4.5D0-1.D0, BI0CS, NTI0) )
      RETURN
C
 20   IF (Y.LE.8.D0) DBSI0Es =
     & (0.375D0 + DCSEVLs ((48.D0/Y-11.D0)/5.D0,
     1  AI0CS, NTAI0))/DSQRT(Y)
      IF (Y.GT.8.D0) DBSI0Es =
     & (0.375D0 + DCSEVLs (16.D0/Y-1.D0, AI02CS,
     1  NTAI02))/DSQRT(Y)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCSEVLs(X,A,N)
c  Stripped down version - no XERROR call- stops on wrong input.
C***BEGIN PROLOGUE  DCSEVLs
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(CSEVL-S DCSEVLs-D),CHEBYSHEV,
C             SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluate the double precision N-term Chebyshev series A
C            at X.
C***DESCRIPTION
C
C Evaluate the N-term Chebyshev series A at X.  Adapted from
C R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).
C W. Fullerton, C-3, Los Alamos Scientific Laboratory.
C
C       Input Arguments --
C X    double precision value at which the series is to be evaluated.
C A    double precision array of N terms of a Chebyshev series.  In
C      evaluating A, only half of the first coefficient is summed.
C N    number of terms in array A.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DCSEVLs
C
       DOUBLE PRECISION A(N),X,TWOX,B0,B1,B2
C***FIRST EXECUTABLE STATEMENT  DCSEVLs
       IF(N.LT.1) go to 1001
       IF(N.GT.1000) go to 1002
       IF ((X.LT.-1.D0) .OR. (X.GT.1.D0)) go to 1003
C
       TWOX = 2.0D0*X
       B1 = 0.D0
       B0=0.D0
       DO 10 I=1,N
         B2=B1
         B1=B0
         NI = N - I + 1
         B0 = TWOX*B1 - B2 + A(NI)
 10    CONTINUE
C
       DCSEVLs = 0.5D0 * (B0-B2)
C
       RETURN
 1001 write(6,201) N
  201 format(1x,'N < 1 in DCSEVLs- stopped')
        stop
 1002 write(6,202) N
  202 format(1x,'N>1000 in DCSEVLs- stopped')
        stop
 1003 write(6,203) x
  203 format(1x,'x outside interval [-1,1] in DCSEVLs- stopped')
        stop
      END
      FUNCTION INITDSs(DOS,NOS,ETA)
C***BEGIN PROLOGUE  INITDS
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(INITS-S INITDS-D),CHEBYSHEV,
C             INITIALIZE,ORTHOGONAL POLYNOMIAL,ORTHOGONAL SERIES,SERIES,
C             SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Initializes the d.p. properly normalized orthogonal
C            polynomial series to determine the number of terms needed
C            for specific accuracy.
C***DESCRIPTION
C
C Initialize the double precision orthogonal series DOS so that INITDS
C is the number of terms needed to insure the error is no larger than
C ETA.  Ordinarily ETA will be chosen to be one-tenth machine precision
C
C             Input Arguments --
C DOS    dble prec array of NOS coefficients in an orthogonal series.
C NOS    number of coefficients in DOS.
C ETA    requested accuracy of series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  INITDS
C
      DOUBLE PRECISION DOS(NOS)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS.LT.1) write(6,201)
  201 format( 'INITDS  NUMBER OF COEFFICIENTS LT 1')
C
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(DOS(I))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
C
 20   IF (I.EQ.NOS) write(6,200)
  200 format( 'INITDSs  ETA MAY BE TOO SMALL')
      INITDSs = I
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBSI1Es(X)
C***BEGIN PROLOGUE  DBSI1E
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI1E-S DBSI1E-D),BESSEL FUNCTION,
C             EXPONENTIALLY SCALED,FIRST KIND,
C             HYPERBOLIC BESSEL FUNCTION,MODIFIED BESSEL FUNCTION,
C             ORDER ONE,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. exponentially scaled modified (hyper-
C            bolic) Bessel function of the first kind of order one.
C***DESCRIPTION
C
C DBSI1E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C one for double precision argument X.  The result is I1(X)
C multiplied by EXP(-ABS(X)).
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C
C Series for AI1        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.81E-32
C                                         log weighted error  31.55
C                               significant figures required  29.93
C                                    decimal places required  32.38
C
C Series for AI12       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.83E-32
C                                         log weighted error  31.74
C                               significant figures required  29.97
C                                    decimal places required  32.66
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBSI1E
      DOUBLE PRECISION X, BI1CS(17), AI1CS(46), AI12CS(69), XMIN,
     1  XSML, Y,  DCSEVLs
      SAVE BI1 CS, AI1 CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML
      DATA BI1 CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1 CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1 CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1 CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1 CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1 CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1 CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1 CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1 CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1 CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1 CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1 CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1 CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1 CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1 CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1 CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1 CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
      DATA AI1 CS(  1) / -.2846744181 8814786741 0037246830 7 D-1      /
      DATA AI1 CS(  2) / -.1922953231 4432206510 4444877497 9 D-1      /
      DATA AI1 CS(  3) / -.6115185857 9437889822 5624991778 5 D-3      /
      DATA AI1 CS(  4) / -.2069971253 3502277088 8282377797 9 D-4      /
      DATA AI1 CS(  5) / +.8585619145 8107255655 3694467313 8 D-5      /
      DATA AI1 CS(  6) / +.1049498246 7115908625 1745399786 0 D-5      /
      DATA AI1 CS(  7) / -.2918338918 4479022020 9343232669 7 D-6      /
      DATA AI1 CS(  8) / -.1559378146 6317390001 6068096907 7 D-7      /
      DATA AI1 CS(  9) / +.1318012367 1449447055 2530287390 9 D-7      /
      DATA AI1 CS( 10) / -.1448423418 1830783176 3913446781 5 D-8      /
      DATA AI1 CS( 11) / -.2908512243 9931420948 2504099301 0 D-9      /
      DATA AI1 CS( 12) / +.1266388917 8753823873 1115969040 3 D-9      /
      DATA AI1 CS( 13) / -.1664947772 9192206706 2417839858 0 D-10     /
      DATA AI1 CS( 14) / -.1666653644 6094329760 9593715499 9 D-11     /
      DATA AI1 CS( 15) / +.1242602414 2907682652 3216847201 7 D-11     /
      DATA AI1 CS( 16) / -.2731549379 6724323972 5146142863 3 D-12     /
      DATA AI1 CS( 17) / +.2023947881 6458037807 0026268898 1 D-13     /
      DATA AI1 CS( 18) / +.7307950018 1168836361 9869812612 3 D-14     /
      DATA AI1 CS( 19) / -.3332905634 4046749438 1377861713 3 D-14     /
      DATA AI1 CS( 20) / +.7175346558 5129537435 4225466567 0 D-15     /
      DATA AI1 CS( 21) / -.6982530324 7962563558 5062922365 6 D-16     /
      DATA AI1 CS( 22) / -.1299944201 5627607600 6044608058 7 D-16     /
      DATA AI1 CS( 23) / +.8120942864 2427988920 5467834286 0 D-17     /
      DATA AI1 CS( 24) / -.2194016207 4107368981 5626664378 3 D-17     /
      DATA AI1 CS( 25) / +.3630516170 0296548482 7986093233 4 D-18     /
      DATA AI1 CS( 26) / -.1695139772 4391041663 0686679039 9 D-19     /
      DATA AI1 CS( 27) / -.1288184829 8979078071 1688253822 2 D-19     /
      DATA AI1 CS( 28) / +.5694428604 9670527801 0999107310 9 D-20     /
      DATA AI1 CS( 29) / -.1459597009 0904800565 4550990028 7 D-20     /
      DATA AI1 CS( 30) / +.2514546010 6757173140 8469133448 5 D-21     /
      DATA AI1 CS( 31) / -.1844758883 1391248181 6040002901 3 D-22     /
      DATA AI1 CS( 32) / -.6339760596 2279486419 2860979199 9 D-23     /
      DATA AI1 CS( 33) / +.3461441102 0310111111 0814662656 0 D-23     /
      DATA AI1 CS( 34) / -.1017062335 3713935475 9654102357 3 D-23     /
      DATA AI1 CS( 35) / +.2149877147 0904314459 6250077866 6 D-24     /
      DATA AI1 CS( 36) / -.3045252425 2386764017 4620617386 6 D-25     /
      DATA AI1 CS( 37) / +.5238082144 7212859821 7763498666 6 D-27     /
      DATA AI1 CS( 38) / +.1443583107 0893824464 1678950399 9 D-26     /
      DATA AI1 CS( 39) / -.6121302074 8900427332 0067071999 9 D-27     /
      DATA AI1 CS( 40) / +.1700011117 4678184183 4918980266 6 D-27     /
      DATA AI1 CS( 41) / -.3596589107 9842441585 3521578666 6 D-28     /
      DATA AI1 CS( 42) / +.5448178578 9484185766 5051306666 6 D-29     /
      DATA AI1 CS( 43) / -.2731831789 6890849891 6256426666 6 D-30     /
      DATA AI1 CS( 44) / -.1858905021 7086007157 7190399999 9 D-30     /
      DATA AI1 CS( 45) / +.9212682974 5139334411 2776533333 3 D-31     /
      DATA AI1 CS( 46) / -.2813835155 6535611063 7083306666 6 D-31     /
      DATA AI12CS(  1) / +.2857623501 8280120474 4984594846 9 D-1      /
      DATA AI12CS(  2) / -.9761097491 3614684077 6516445730 2 D-2      /
      DATA AI12CS(  3) / -.1105889387 6262371629 1256921277 5 D-3      /
      DATA AI12CS(  4) / -.3882564808 8776903934 5654477627 4 D-5      /
      DATA AI12CS(  5) / -.2512236237 8702089252 9452002212 1 D-6      /
      DATA AI12CS(  6) / -.2631468846 8895195068 3705236523 2 D-7      /
      DATA AI12CS(  7) / -.3835380385 9642370220 4500678796 8 D-8      /
      DATA AI12CS(  8) / -.5589743462 1965838068 6811252222 9 D-9      /
      DATA AI12CS(  9) / -.1897495812 3505412344 9892503323 8 D-10     /
      DATA AI12CS( 10) / +.3252603583 0154882385 5508067994 9 D-10     /
      DATA AI12CS( 11) / +.1412580743 6613781331 6336633284 6 D-10     /
      DATA AI12CS( 12) / +.2035628544 1470895072 2452613684 0 D-11     /
      DATA AI12CS( 13) / -.7198551776 2459085120 9258989044 6 D-12     /
      DATA AI12CS( 14) / -.4083551111 0921973182 2849963969 1 D-12     /
      DATA AI12CS( 15) / -.2101541842 7726643130 1984572746 2 D-13     /
      DATA AI12CS( 16) / +.4272440016 7119513542 9778833699 7 D-13     /
      DATA AI12CS( 17) / +.1042027698 4128802764 1741449994 8 D-13     /
      DATA AI12CS( 18) / -.3814403072 4370078047 6707253539 6 D-14     /
      DATA AI12CS( 19) / -.1880354775 5107824485 1273453396 3 D-14     /
      DATA AI12CS( 20) / +.3308202310 9209282827 3190335240 5 D-15     /
      DATA AI12CS( 21) / +.2962628997 6459501390 6854654205 2 D-15     /
      DATA AI12CS( 22) / -.3209525921 9934239587 7837353288 7 D-16     /
      DATA AI12CS( 23) / -.4650305368 4893583255 7128281897 9 D-16     /
      DATA AI12CS( 24) / +.4414348323 0717079499 4611375964 1 D-17     /
      DATA AI12CS( 25) / +.7517296310 8421048054 2545808029 5 D-17     /
      DATA AI12CS( 26) / -.9314178867 3268833756 8484784515 7 D-18     /
      DATA AI12CS( 27) / -.1242193275 1948909561 1678448869 7 D-17     /
      DATA AI12CS( 28) / +.2414276719 4548484690 0515390217 6 D-18     /
      DATA AI12CS( 29) / +.2026944384 0532851789 7192286069 2 D-18     /
      DATA AI12CS( 30) / -.6394267188 2690977870 4391988681 1 D-19     /
      DATA AI12CS( 31) / -.3049812452 3730958960 8488450357 1 D-19     /
      DATA AI12CS( 32) / +.1612841851 6514802251 3462230769 1 D-19     /
      DATA AI12CS( 33) / +.3560913964 3099250545 1027090462 0 D-20     /
      DATA AI12CS( 34) / -.3752017947 9364390796 6682800324 6 D-20     /
      DATA AI12CS( 35) / -.5787037427 0747993459 5198231074 1 D-22     /
      DATA AI12CS( 36) / +.7759997511 6481619619 8236963209 2 D-21     /
      DATA AI12CS( 37) / -.1452790897 2022333940 6445987408 5 D-21     /
      DATA AI12CS( 38) / -.1318225286 7390367021 2192275337 4 D-21     /
      DATA AI12CS( 39) / +.6116654862 9030707018 7999133171 7 D-22     /
      DATA AI12CS( 40) / +.1376279762 4271264277 3024338363 4 D-22     /
      DATA AI12CS( 41) / -.1690837689 9593478849 1983938230 6 D-22     /
      DATA AI12CS( 42) / +.1430596088 5954331539 8720108538 5 D-23     /
      DATA AI12CS( 43) / +.3409557828 0905940204 0536772990 2 D-23     /
      DATA AI12CS( 44) / -.1309457666 2707602278 4573872642 4 D-23     /
      DATA AI12CS( 45) / -.3940706411 2402574360 9352141755 7 D-24     /
      DATA AI12CS( 46) / +.4277137426 9808765808 0616679735 2 D-24     /
      DATA AI12CS( 47) / -.4424634830 9826068819 0028312302 9 D-25     /
      DATA AI12CS( 48) / -.8734113196 2307149721 1530978874 7 D-25     /
      DATA AI12CS( 49) / +.4045401335 6835333921 4340414242 8 D-25     /
      DATA AI12CS( 50) / +.7067100658 0946894656 5160771780 6 D-26     /
      DATA AI12CS( 51) / -.1249463344 5651052230 0286451860 5 D-25     /
      DATA AI12CS( 52) / +.2867392244 4034370329 7948339142 6 D-26     /
      DATA AI12CS( 53) / +.2044292892 5042926702 8177957421 0 D-26     /
      DATA AI12CS( 54) / -.1518636633 8204625683 7134680291 1 D-26     /
      DATA AI12CS( 55) / +.8110181098 1875758861 3227910703 7 D-28     /
      DATA AI12CS( 56) / +.3580379354 7735860911 2717370327 0 D-27     /
      DATA AI12CS( 57) / -.1692929018 9279025095 9305717544 8 D-27     /
      DATA AI12CS( 58) / -.2222902499 7024276390 6775852777 4 D-28     /
      DATA AI12CS( 59) / +.5424535127 1459696550 4860040112 8 D-28     /
      DATA AI12CS( 60) / -.1787068401 5780186887 6491299330 4 D-28     /
      DATA AI12CS( 61) / -.6565479068 7228149388 2392943788 0 D-29     /
      DATA AI12CS( 62) / +.7807013165 0611452809 2206770683 9 D-29     /
      DATA AI12CS( 63) / -.1816595260 6689797173 7933315222 1 D-29     /
      DATA AI12CS( 64) / -.1287704952 6600848203 7687559895 9 D-29     /
      DATA AI12CS( 65) / +.1114548172 9881645474 1370927369 4 D-29     /
      DATA AI12CS( 66) / -.1808343145 0393369391 5936887668 7 D-30     /
      DATA AI12CS( 67) / -.2231677718 2037719522 3244822893 9 D-30     /
      DATA AI12CS( 68) / +.1619029596 0803415106 1790980361 4 D-30     /
      DATA AI12CS( 69) / -.1834079908 8049414139 0130843921 0 D-31     /
      DATA NTI1, NTAI1, NTAI12, XMIN, XSML / 3*0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBSI1Es
      IF (NTI1.NE.0) GO TO 10
c        type *,'D1mach(3)',d1mach(3)
        eta=1.4e-18
c      ETA = 0.1*SNGL(D1MACH(3))
c        type *,'eta',eta
      NTI1 = INITDSs (BI1CS, 17, ETA)
      NTAI1 = INITDSs (AI1CS, 46, ETA)
      NTAI12 = INITDSs (AI12CS, 69, ETA)
C
      XMIN = 4.d-39
c        type *,'D1mach(1)',d1mach(1)
      XSML = 1.d-8
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBSI1Es = 0.0D0
      IF (Y.EQ.0.D0)  RETURN
C
      IF (Y.LT.XMIN) write(6,301)
  301 format(1x,'DBSI1Es  DABS(X) SO SMALL I1 UNDERFLOWS')
      IF (Y.GT.XMIN) DBSI1Es = 0.5D0*X
      IF (Y.GT.XSML) DBSI1Es =
     & X*(0.875D0 + DCSEVLs (Y*Y/4.5D0-1.D0,
     1  BI1CS, NTI1) )
      DBSI1Es = DEXP(-Y) * DBSI1Es
      RETURN
C
 20   IF (Y.LE.8.D0) DBSI1Es =
     & (0.375D0 + DCSEVLs ((48.D0/Y-11.D0)/5.D0,
     1  AI1CS, NTAI1))/DSQRT(Y)
      IF (Y.GT.8.D0) DBSI1Es =
     & (0.375D0 + DCSEVLs (16.D0/Y-1.D0, AI12CS,
     1  NTAI12))/DSQRT(Y)
      DBSI1Es = DSIGN (DBSI1Es, X)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBESI1s(X)
C***BEGIN PROLOGUE  DBESI1
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI1-S DBESI1-D),BESSEL FUNCTION,
C             FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ONE,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. modified (hyperbolic) Bessel function
C            of the first kind of order one.
C***DESCRIPTION
C
C DBESI1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order one and double precision
C argument X.
C
C Series for BI1        on the interval  0.          to  9.00000E+00
C                                        with weighted error   1.44E-32
C                                         log weighted error  31.84
C                               significant figures required  31.45
C                                    decimal places required  32.46
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBSI1E,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBESI1
      DOUBLE PRECISION X, BI1CS(17), XMAX, XMIN, XSML, Y,
     1  DCSEVLs, DBSI1Es
      SAVE BI1 CS, NTI1, XMIN, XSML, XMAX
      DATA BI1 CS(  1) / -.1971713261 0998597316 1385032181 49 D-2     /
      DATA BI1 CS(  2) / +.4073488766 7546480608 1553936520 14 D+0     /
      DATA BI1 CS(  3) / +.3483899429 9959455866 2450377837 87 D-1     /
      DATA BI1 CS(  4) / +.1545394556 3001236038 5984010584 89 D-2     /
      DATA BI1 CS(  5) / +.4188852109 8377784129 4588320041 20 D-4     /
      DATA BI1 CS(  6) / +.7649026764 8362114741 9597039660 69 D-6     /
      DATA BI1 CS(  7) / +.1004249392 4741178689 1798080372 38 D-7     /
      DATA BI1 CS(  8) / +.9932207791 9238106481 3712980548 63 D-10    /
      DATA BI1 CS(  9) / +.7663801791 8447637275 2001716813 49 D-12    /
      DATA BI1 CS( 10) / +.4741418923 8167394980 3880919481 60 D-14    /
      DATA BI1 CS( 11) / +.2404114404 0745181799 8631720320 00 D-16    /
      DATA BI1 CS( 12) / +.1017150500 7093713649 1211007999 99 D-18    /
      DATA BI1 CS( 13) / +.3645093565 7866949458 4917333333 33 D-21    /
      DATA BI1 CS( 14) / +.1120574950 2562039344 8106666666 66 D-23    /
      DATA BI1 CS( 15) / +.2987544193 4468088832 0000000000 00 D-26    /
      DATA BI1 CS( 16) / +.6973231093 9194709333 3333333333 33 D-29    /
      DATA BI1 CS( 17) / +.1436794822 0620800000 0000000000 00 D-31    /
      DATA NTI1, XMIN, XSML, XMAX / 0, 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESI1
      IF (NTI1.NE.0) GO TO 10
        eta=1.4d-18
      NTI1 = INITDSs (BI1CS, 17, eta)
      XMIN = 4.0d-39
      XSML = 1.d-8
      XMAX = 86.d0
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI1s = 0.D0
      IF (Y.EQ.0.D0)  RETURN
C
      IF (Y.LT.XMIN)
     & write(6,*) 'DBESI1s  DABS(X) SO SMALL IT UNDERFLOWS'
      IF (Y.GT.XMIN) DBESI1s = 0.5D0*X
      IF (Y.GT.XSML) DBESI1s = X*(0.875D0 +
     & DCSEVLs (Y*Y/4.5D0-1.D0,
     1  BI1CS, NTI1))
      RETURN
C
 20   IF (Y.GT.XMAX) write(6,*)'DBESI1s  DABS(X) SO BIG I1 OVERFLOWS'
C
      DBESI1s = DEXP(Y) * DBSI1Es(X)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBESI0s(X)
c  Stripped version of CLAMS DBESI0 - no machine-dependent calls.
C***BEGIN PROLOGUE  DBESI0
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI0-S DBESI0-D),BESSEL FUNCTION,
C             FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ZERO,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. hyperbolic Bessel function of the first
C            kind of order zero.
C***DESCRIPTION
C
C DBESI0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order zero and double
C precision argument X.
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBSI0E,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBESI0
      DOUBLE PRECISION X, BI0CS(18), XMAX, XSML, Y,
     1  DCSEVLs, DBSI0Es
      SAVE BI0 CS, NTI0, XSML, XMAX
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA NTI0, XSML, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESI0
      IF (NTI0.NE.0) GO TO 10
      NTI0 = INITDSs(BI0CS, 18, 1.4e-18)
c      XSML = DSQRT (8.0D0*D1MACH(3))
        xsml=7.5d-9
c      XMAX = DLOG (D1MACH(2))
        xmax=86.4d0
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI0s= 1.0D0
      IF (Y.GT.XSML) DBESI0s= 2.75D0 + DCSEVLs(Y*Y/4.5D0-1.D0, BI0CS,
     1  NTI0)
      RETURN
C
 20   IF (Y.GT.XMAX) go to 1001
C
      DBESI0s= DEXP(Y) * DBSI0Es(X)
C
      RETURN
 1001 write(6,*) 'X>XMAX in DBESI0s-stopped:XMAX=86.4,x=',x
        stop
      END
      DOUBLE PRECISION FUNCTION DBESK0s(X)
c  Stripped of machine-dependent calls .
C***BEGIN PROLOGUE  DBESK0s
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESK0-S DBESK0s-D),BESSEL FUNCTION,
C             HYPERBOLIC BESSEL FUNCTION,MODIFIED BESSEL FUNCTION,
C             ORDER ZERO,SPECIAL FUNCTIONS,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes  d.p. modified (hyperbolic) Bessel function of
C            the third kind of order zero.
C***DESCRIPTION
C
C DBESK0s(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order zero for double
C precision argument X.  The argument must be greater than zero
C but not so large that the result underflows.
C
C Series for BK0        on the interval  0.          to  4.00000E+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DBESI0,DBSK0E,DCSEVL,INITDS
C***END PROLOGUE  DBESK0s
      DOUBLE PRECISION X, BK0CS(16), XMAX,  XSML, Y,
     1  DCSEVLs, DBESI0s, DBSK0Es
      SAVE BK0 CS, NTK0, XSML, XMAX
      DATA BK0 CS(  1) / -.3532739323 3902768720 1140060063 153 D-1    /
      DATA BK0 CS(  2) / +.3442898999 2462848688 6344927529 213 D+0    /
      DATA BK0 CS(  3) / +.3597993651 5361501626 5721303687 231 D-1    /
      DATA BK0 CS(  4) / +.1264615411 4469259233 8479508673 447 D-2    /
      DATA BK0 CS(  5) / +.2286212103 1194517860 8269830297 585 D-4    /
      DATA BK0 CS(  6) / +.2534791079 0261494573 0790013428 354 D-6    /
      DATA BK0 CS(  7) / +.1904516377 2202088589 7214059381 366 D-8    /
      DATA BK0 CS(  8) / +.1034969525 7633624585 1008317853 089 D-10   /
      DATA BK0 CS(  9) / +.4259816142 7910825765 2445327170 133 D-13   /
      DATA BK0 CS( 10) / +.1374465435 8807508969 4238325440 000 D-15   /
      DATA BK0 CS( 11) / +.3570896528 5083735909 9688597333 333 D-18   /
      DATA BK0 CS( 12) / +.7631643660 1164373766 7498666666 666 D-21   /
      DATA BK0 CS( 13) / +.1365424988 4407818590 8053333333 333 D-23   /
      DATA BK0 CS( 14) / +.2075275266 9066680831 9999999999 999 D-26   /
      DATA BK0 CS( 15) / +.2712814218 0729856000 0000000000 000 D-29   /
      DATA BK0 CS( 16) / +.3082593887 9146666666 6666666666 666 D-32   /
      DATA NTK0, XSML, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESK0s
      IF (NTK0.NE.0) GO TO 10
      NTK0 = INITDSs(BK0CS, 16, 1.4e-18)
      XSML = 7.5d-9
      XMAX = 86.4
C
   10 if(x.le.0.d0) go to 1001
      IF (X.GT.2.0D0) GO TO 20
C
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBESK0s = -DLOG(0.5D0*X)*DBESI0s(X)
     & - 0.25D0 + DCSEVLs(.5D0*Y-1.D0,
     1  BK0CS, NTK0)
      RETURN
C
 20   DBESK0s = 0.D0
      IF (X.GT.XMAX) write(6,*)'DBESK0s  X SO BIG K0 UNDERFLOWS'
      IF (X.GT.XMAX) RETURN
C
      DBESK0s = DEXP(-X) * DBSK0Es(X)
C
      RETURN
 1001  write(6,*) 'DBESK0s  X IS ZERO OR NEGATIVE- stopped'
        stop
      END
      DOUBLE PRECISION FUNCTION DBSK0Es(X)
c  Stripped of machine-dependent calls.
C***BEGIN PROLOGUE  DBSK0Es
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESK0E-S DBSK0Es-D),BESSEL FUNCTION,
C             EXPONENTIALLY SCALED,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ZERO,SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. exponentially scaled modified (hyper-
C            bolic) Bessel function of the third kind of order zero.
C***DESCRIPTION
C
C DBSK0Es(X) computes the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the third kind of
C order zero for positive double precision argument X.
C
C Series for BK0        on the interval  0.          to  4.00000E+00
C                                        with weighted error   3.08E-33
C                                         log weighted error  32.51
C                               significant figures required  32.05
C                                    decimal places required  33.11
C
C Series for AK0        on the interval  1.25000E-01 to  5.00000E-01
C                                        with weighted error   2.85E-32
C                                         log weighted error  31.54
C                               significant figures required  30.19
C                                    decimal places required  32.33
C
C Series for AK02       on the interval  0.          to  1.25000E-01
C                                        with weighted error   2.30E-32
C                                         log weighted error  31.64
C                               significant figures required  29.68
C                                    decimal places required  32.40
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBESI0,DCSEVL,INITDSs,XERROR
C***END PROLOGUE  DBSK0Es
      DOUBLE PRECISION X, BK0CS(16), AK0CS(38), AK02CS(33),
     1  XSML, Y,  DCSEVLs,DBESI0s
      SAVE BK0 CS, AK0 CS, AK02CS,NTK0, NTAK0, NTAK02, XSML
      DATA BK0 CS(  1) / -.3532739323 3902768720 1140060063 153 D-1    /
      DATA BK0 CS(  2) / +.3442898999 2462848688 6344927529 213 D+0    /
      DATA BK0 CS(  3) / +.3597993651 5361501626 5721303687 231 D-1    /
      DATA BK0 CS(  4) / +.1264615411 4469259233 8479508673 447 D-2    /
      DATA BK0 CS(  5) / +.2286212103 1194517860 8269830297 585 D-4    /
      DATA BK0 CS(  6) / +.2534791079 0261494573 0790013428 354 D-6    /
      DATA BK0 CS(  7) / +.1904516377 2202088589 7214059381 366 D-8    /
      DATA BK0 CS(  8) / +.1034969525 7633624585 1008317853 089 D-10   /
      DATA BK0 CS(  9) / +.4259816142 7910825765 2445327170 133 D-13   /
      DATA BK0 CS( 10) / +.1374465435 8807508969 4238325440 000 D-15   /
      DATA BK0 CS( 11) / +.3570896528 5083735909 9688597333 333 D-18   /
      DATA BK0 CS( 12) / +.7631643660 1164373766 7498666666 666 D-21   /
      DATA BK0 CS( 13) / +.1365424988 4407818590 8053333333 333 D-23   /
      DATA BK0 CS( 14) / +.2075275266 9066680831 9999999999 999 D-26   /
      DATA BK0 CS( 15) / +.2712814218 0729856000 0000000000 000 D-29   /
      DATA BK0 CS( 16) / +.3082593887 9146666666 6666666666 666 D-32   /
      DATA AK0 CS(  1) / -.7643947903 3279414240 8297827008 8 D-1      /
      DATA AK0 CS(  2) / -.2235652605 6998190520 2309555079 1 D-1      /
      DATA AK0 CS(  3) / +.7734181154 6938582353 0061817404 7 D-3      /
      DATA AK0 CS(  4) / -.4281006688 8860994644 5214643541 6 D-4      /
      DATA AK0 CS(  5) / +.3081700173 8629747436 5001482666 0 D-5      /
      DATA AK0 CS(  6) / -.2639367222 0096649740 6744889272 3 D-6      /
      DATA AK0 CS(  7) / +.2563713036 4034692062 9408826574 2 D-7      /
      DATA AK0 CS(  8) / -.2742705549 9002012638 5721191524 4 D-8      /
      DATA AK0 CS(  9) / +.3169429658 0974995920 8083287340 3 D-9      /
      DATA AK0 CS( 10) / -.3902353286 9621841416 0106571796 2 D-10     /
      DATA AK0 CS( 11) / +.5068040698 1885754020 5009212728 6 D-11     /
      DATA AK0 CS( 12) / -.6889574741 0078706795 4171355798 4 D-12     /
      DATA AK0 CS( 13) / +.9744978497 8259176913 8820133683 1 D-13     /
      DATA AK0 CS( 14) / -.1427332841 8845485053 8985534012 2 D-13     /
      DATA AK0 CS( 15) / +.2156412571 0214630395 5806297652 7 D-14     /
      DATA AK0 CS( 16) / -.3349654255 1495627721 8878205853 0 D-15     /
      DATA AK0 CS( 17) / +.5335260216 9529116921 4528039260 1 D-16     /
      DATA AK0 CS( 18) / -.8693669980 8907538076 3962237883 7 D-17     /
      DATA AK0 CS( 19) / +.1446404347 8622122278 8776344234 6 D-17     /
      DATA AK0 CS( 20) / -.2452889825 5001296824 0467875157 3 D-18     /
      DATA AK0 CS( 21) / +.4233754526 2321715728 2170634240 0 D-19     /
      DATA AK0 CS( 22) / -.7427946526 4544641956 9534129493 3 D-20     /
      DATA AK0 CS( 23) / +.1323150529 3926668662 7796746240 0 D-20     /
      DATA AK0 CS( 24) / -.2390587164 7396494513 3598146559 9 D-21     /
      DATA AK0 CS( 25) / +.4376827585 9232261401 6571255466 6 D-22     /
      DATA AK0 CS( 26) / -.8113700607 3451180593 3901141333 3 D-23     /
      DATA AK0 CS( 27) / +.1521819913 8321729583 1037815466 6 D-23     /
      DATA AK0 CS( 28) / -.2886041941 4833977702 3595861333 3 D-24     /
      DATA AK0 CS( 29) / +.5530620667 0547179799 9261013333 3 D-25     /
      DATA AK0 CS( 30) / -.1070377329 2498987285 9163306666 6 D-25     /
      DATA AK0 CS( 31) / +.2091086893 1423843002 9632853333 3 D-26     /
      DATA AK0 CS( 32) / -.4121713723 6462038274 1026133333 3 D-27     /
      DATA AK0 CS( 33) / +.8193483971 1213076401 3568000000 0 D-28     /
      DATA AK0 CS( 34) / -.1642000275 4592977267 8075733333 3 D-28     /
      DATA AK0 CS( 35) / +.3316143281 4802271958 9034666666 6 D-29     /
      DATA AK0 CS( 36) / -.6746863644 1452959410 8586666666 6 D-30     /
      DATA AK0 CS( 37) / +.1382429146 3184246776 3541333333 3 D-30     /
      DATA AK0 CS( 38) / -.2851874167 3598325708 1173333333 3 D-31     /
      DATA AK02CS(  1) / -.1201869826 3075922398 3934621245 2 D-1      /
      DATA AK02CS(  2) / -.9174852691 0256953106 5256107571 3 D-2      /
      DATA AK02CS(  3) / +.1444550931 7750058210 4884387805 7 D-3      /
      DATA AK02CS(  4) / -.4013614175 4357097286 7102107787 9 D-5      /
      DATA AK02CS(  5) / +.1567831810 8523106725 9034899033 3 D-6      /
      DATA AK02CS(  6) / -.7770110438 5217377103 1579975446 0 D-8      /
      DATA AK02CS(  7) / +.4611182576 1797178825 3313052958 6 D-9      /
      DATA AK02CS(  8) / -.3158592997 8605657705 2666580330 9 D-10     /
      DATA AK02CS(  9) / +.2435018039 3650411278 3588781432 9 D-11     /
      DATA AK02CS( 10) / -.2074331387 3983478977 0985337350 6 D-12     /
      DATA AK02CS( 11) / +.1925787280 5899170847 4273650469 3 D-13     /
      DATA AK02CS( 12) / -.1927554805 8389561036 0034718221 8 D-14     /
      DATA AK02CS( 13) / +.2062198029 1978182782 8523786964 4 D-15     /
      DATA AK02CS( 14) / -.2341685117 5792424026 0364019507 1 D-16     /
      DATA AK02CS( 15) / +.2805902810 6430422468 1517882845 8 D-17     /
      DATA AK02CS( 16) / -.3530507631 1618079458 1548246357 3 D-18     /
      DATA AK02CS( 17) / +.4645295422 9351082674 2421633706 6 D-19     /
      DATA AK02CS( 18) / -.6368625941 3442664739 2205346133 3 D-20     /
      DATA AK02CS( 19) / +.9069521310 9865155676 2234880000 0 D-21     /
      DATA AK02CS( 20) / -.1337974785 4236907398 4500531199 9 D-21     /
      DATA AK02CS( 21) / +.2039836021 8599523155 2208896000 0 D-22     /
      DATA AK02CS( 22) / -.3207027481 3678405000 6086997333 3 D-23     /
      DATA AK02CS( 23) / +.5189744413 6623099636 2635946666 6 D-24     /
      DATA AK02CS( 24) / -.8629501497 5405721929 6460799999 9 D-25     /
      DATA AK02CS( 25) / +.1472161183 1025598552 0803840000 0 D-25     /
      DATA AK02CS( 26) / -.2573069023 8670112838 1235199999 9 D-26     /
      DATA AK02CS( 27) / +.4601774086 6435165873 7664000000 0 D-27     /
      DATA AK02CS( 28) / -.8411555324 2010937371 3066666666 6 D-28     /
      DATA AK02CS( 29) / +.1569806306 6353689393 0154666666 6 D-28     /
      DATA AK02CS( 30) / -.2988226453 0057577889 7919999999 9 D-29     /
      DATA AK02CS( 31) / +.5796831375 2168365206 1866666666 6 D-30     /
      DATA AK02CS( 32) / -.1145035994 3476813321 5573333333 3 D-30     /
      DATA AK02CS( 33) / +.2301266594 2496828020 0533333333 3 D-31     /
      DATA NTK0, NTAK0, NTAK02, XSML / 3*0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  DBSK0Es
      IF (NTK0.NE.0) GO TO 10
c        type *,d1mach(3)
      ETA = 1.d-18
      NTK0 = INITDSs (BK0CS, 16, ETA)
      NTAK0 = INITDSs (AK0CS, 38, ETA)
      NTAK02 = INITDSs (AK02CS, 33, ETA)
      XSML = 1.d-9
c        type *,xsml
C
 10   IF (X.LE.0.D0) write(6,*)'DBSK0Es  X IS ZERO OR NEGATIVE'
      IF (X.GT.2.0D0) GO TO 20
C
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBSK0Es =
     & DEXP(X)*(-DLOG(0.5D0*X)*DBESI0s(X) - 0.25D0 +
     1  DCSEVLs (.5D0*Y-1.D0, BK0CS, NTK0))
      RETURN
C
 20   IF (X.LE.8.D0) DBSK0Es =
     & (1.25D0 + DCSEVLs ((16.D0/X-5.D0)/3.D0,
     1  AK0CS, NTAK0))/DSQRT(X)
      IF (X.GT.8.D0) DBSK0Es = (1.25D0 +
     1  DCSEVLs (16.D0/X-1.D0, AK02CS, NTAK02))/DSQRT(X)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBSK1Es(X)
c  Stripped of machine-dependent calls- same as CLAMS DBSK1E
C***BEGIN PROLOGUE  DBSK1E
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESK1E-S DBSK1E-D),BESSEL FUNCTION,
C             EXPONENTIALLY SCALED,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ONE,SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the exponentially scaled,modified (hyperbolic)
C            Bessel function of the third kind of order one (double
C            precision).
C***DESCRIPTION
C
C DBSK1E(S) computes the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the third kind of order
C one for positive double precision argument X.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
C Series for AK1        on the interval  1.25000E-01 to  5.00000E-01
C                                        with weighted error   3.07E-32
C                                         log weighted error  31.51
C                               significant figures required  30.71
C                                    decimal places required  32.30
C
C Series for AK12       on the interval  0.          to  1.25000E-01
C                                        with weighted error   2.41E-32
C                                         log weighted error  31.62
C                               significant figures required  30.25
C                                    decimal places required  32.38
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DBESI1,DCSEVLs,INITDS
C***END PROLOGUE  DBSK1E
      DOUBLE PRECISION X, BK1CS(16), AK1CS(38), AK12CS(33), XMIN,
     1  XSML, Y,  DCSEVLs, DBESI1s
      SAVE BK1 CS, AK1 CS, AK12CS, NTK1, NTAK1, NTAK12, XMIN, XSML
      DATA BK1 CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1 CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1 CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1 CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1 CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1 CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1 CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1 CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1 CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1 CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1 CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1 CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1 CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1 CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1 CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1 CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
      DATA AK1 CS(  1) / +.2744313406 9738829695 2576662272 66 D+0     /
      DATA AK1 CS(  2) / +.7571989953 1993678170 8923781492 90 D-1     /
      DATA AK1 CS(  3) / -.1441051556 4754061229 8531161756 25 D-2     /
      DATA AK1 CS(  4) / +.6650116955 1257479394 2513854770 36 D-4     /
      DATA AK1 CS(  5) / -.4369984709 5201407660 5808450891 67 D-5     /
      DATA AK1 CS(  6) / +.3540277499 7630526799 4171390085 34 D-6     /
      DATA AK1 CS(  7) / -.3311163779 2932920208 9826882457 04 D-7     /
      DATA AK1 CS(  8) / +.3445977581 9010534532 3114997709 92 D-8     /
      DATA AK1 CS(  9) / -.3898932347 4754271048 9819374927 58 D-9     /
      DATA AK1 CS( 10) / +.4720819750 4658356400 9474493390 05 D-10    /
      DATA AK1 CS( 11) / -.6047835662 8753562345 3735915628 90 D-11    /
      DATA AK1 CS( 12) / +.8128494874 8658747888 1938379856 63 D-12    /
      DATA AK1 CS( 13) / -.1138694574 7147891428 9239159510 42 D-12    /
      DATA AK1 CS( 14) / +.1654035840 8462282325 9729482050 90 D-13    /
      DATA AK1 CS( 15) / -.2480902567 7068848221 5160104405 33 D-14    /
      DATA AK1 CS( 16) / +.3829237890 7024096948 4292272991 57 D-15    /
      DATA AK1 CS( 17) / -.6064734104 0012418187 7682103773 86 D-16    /
      DATA AK1 CS( 18) / +.9832425623 2648616038 1940046506 66 D-17    /
      DATA AK1 CS( 19) / -.1628416873 8284380035 6666201156 26 D-17    /
      DATA AK1 CS( 20) / +.2750153649 6752623718 2841203370 66 D-18    /
      DATA AK1 CS( 21) / -.4728966646 3953250924 2810695680 00 D-19    /
      DATA AK1 CS( 22) / +.8268150002 8109932722 3920503466 66 D-20    /
      DATA AK1 CS( 23) / -.1468140513 6624956337 1939648853 33 D-20    /
      DATA AK1 CS( 24) / +.2644763926 9208245978 0858948266 66 D-21    /
      DATA AK1 CS( 25) / -.4829015756 4856387897 9698688000 00 D-22    /
      DATA AK1 CS( 26) / +.8929302074 3610130180 6563327999 99 D-23    /
      DATA AK1 CS( 27) / -.1670839716 8972517176 9977514666 66 D-23    /
      DATA AK1 CS( 28) / +.3161645603 4040694931 3686186666 66 D-24    /
      DATA AK1 CS( 29) / -.6046205531 2274989106 5064106666 66 D-25    /
      DATA AK1 CS( 30) / +.1167879894 2042732700 7184213333 33 D-25    /
      DATA AK1 CS( 31) / -.2277374158 2653996232 8678400000 00 D-26    /
      DATA AK1 CS( 32) / +.4481109730 0773675795 3058133333 33 D-27    /
      DATA AK1 CS( 33) / -.8893288476 9020194062 3360000000 00 D-28    /
      DATA AK1 CS( 34) / +.1779468001 8850275131 3920000000 00 D-28    /
      DATA AK1 CS( 35) / -.3588455596 7329095821 9946666666 66 D-29    /
      DATA AK1 CS( 36) / +.7290629049 2694257991 6799999999 99 D-30    /
      DATA AK1 CS( 37) / -.1491844984 5546227073 0240000000 00 D-30    /
      DATA AK1 CS( 38) / +.3073657387 2934276300 7999999999 99 D-31    /
      DATA AK12CS(  1) / +.6379308343 7390010366 0048853410 2 D-1      /
      DATA AK12CS(  2) / +.2832887813 0497209358 3503028470 8 D-1      /
      DATA AK12CS(  3) / -.2475370673 9052503454 1454556673 2 D-3      /
      DATA AK12CS(  4) / +.5771972451 6072488204 7097662576 3 D-5      /
      DATA AK12CS(  5) / -.2068939219 5365483027 4553319655 2 D-6      /
      DATA AK12CS(  6) / +.9739983441 3818041803 0921309788 7 D-8      /
      DATA AK12CS(  7) / -.5585336140 3806249846 8889551112 9 D-9      /
      DATA AK12CS(  8) / +.3732996634 0461852402 2121285473 1 D-10     /
      DATA AK12CS(  9) / -.2825051961 0232254451 3506575492 8 D-11     /
      DATA AK12CS( 10) / +.2372019002 4841441736 4349695548 6 D-12     /
      DATA AK12CS( 11) / -.2176677387 9917539792 6830166793 8 D-13     /
      DATA AK12CS( 12) / +.2157914161 6160324539 3956268970 6 D-14     /
      DATA AK12CS( 13) / -.2290196930 7182692759 9155133815 4 D-15     /
      DATA AK12CS( 14) / +.2582885729 8232749619 1993956522 6 D-16     /
      DATA AK12CS( 15) / -.3076752641 2684631876 2109817344 0 D-17     /
      DATA AK12CS( 16) / +.3851487721 2804915970 9489684479 9 D-18     /
      DATA AK12CS( 17) / -.5044794897 6415289771 1728250880 0 D-19     /
      DATA AK12CS( 18) / +.6888673850 4185442370 1829222399 9 D-20     /
      DATA AK12CS( 19) / -.9775041541 9501183030 0213248000 0 D-21     /
      DATA AK12CS( 20) / +.1437416218 5238364610 0165973333 3 D-21     /
      DATA AK12CS( 21) / -.2185059497 3443473734 9973333333 3 D-22     /
      DATA AK12CS( 22) / +.3426245621 8092206316 4538880000 0 D-23     /
      DATA AK12CS( 23) / -.5531064394 2464082325 0124800000 0 D-24     /
      DATA AK12CS( 24) / +.9176601505 6859954037 8282666666 6 D-25     /
      DATA AK12CS( 25) / -.1562287203 6180249114 4874666666 6 D-25     /
      DATA AK12CS( 26) / +.2725419375 4843331323 4943999999 9 D-26     /
      DATA AK12CS( 27) / -.4865674910 0748279923 7802666666 6 D-27     /
      DATA AK12CS( 28) / +.8879388552 7235025873 5786666666 6 D-28     /
      DATA AK12CS( 29) / -.1654585918 0392575489 3653333333 3 D-28     /
      DATA AK12CS( 30) / +.3145111321 3578486743 0399999999 9 D-29     /
      DATA AK12CS( 31) / -.6092998312 1931276124 1600000000 0 D-30     /
      DATA AK12CS( 32) / +.1202021939 3698158346 2399999999 9 D-30     /
      DATA AK12CS( 33) / -.2412930801 4594088413 8666666666 6 D-31     /
      DATA NTK1, NTAK1, NTAK12, XMIN, XSML / 3*0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBSK1Es
      IF (NTK1.NE.0) GO TO 10
c      ETA = 0.1*SNGL(D1MACH(3))
        eta=1.4e-18
      NTK1 = INITDSs (BK1CS, 16, ETA)
      NTAK1 = INITDSs (AK1CS, 38, ETA)
      NTAK12 = INITDSs (AK12CS, 33, ETA)
C
c      XMIN = DEXP (DMAX1(DLOG(D1MACH(1)), -DLOG(D1MACH(2))) + 0.01D0)
        xmin=5.94d-39
c      XSML = DSQRT (4.0D0*D1MACH(3))
        xsml=7.5d-9
C
 10   IF (X.LE.0.D0) write(6,*)'DBSK1Es  X IS ZERO OR NEGATIVE'
      IF (X.GT.2.0D0) GO TO 20
C
      IF (X.LT.XMIN) write(6,*)'DBSK1Es  X SO SMALL K1 OVERFLOWS'
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBSK1Es = DEXP(X)*(DLOG(0.5D0*X)*DBESI1s(X) + (0.75D0 +
     1  DCSEVLs (0.5D0*Y-1.D0, BK1CS, NTK1))/X )
      RETURN
C
 20   IF (X.LE.8.D0) DBSK1Es =
     & (1.25D0 + DCSEVLs ((16.D0/X-5.D0)/3.D0,
     1  AK1CS, NTAK1))/DSQRT(X)
      IF (X.GT.8.D0) DBSK1Es = (1.25D0 +
     1  DCSEVLs (16.D0/X-1.D0, AK12CS, NTAK12))/DSQRT(X)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBESK1s(X)
C***BEGIN PROLOGUE  DBESK1
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESK1-S DBESK1-D),BESSEL FUNCTION,
C             HYPERBOLIC BESSEL FUNCTION,MODIFIED BESSEL FUNCTION,
C             ORDER ONE,SPECIAL FUNCTIONS,THIRD KIND
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the dp modified Bessel function of the third kind
C            of order one.
C***DESCRIPTION
C
C DBESK1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order one for double precision
C argument X.  The argument must be large enough that the result does
C not overflow and small enough that the result does not underflow.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBESI1,DBSK1E,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBESK1
      DOUBLE PRECISION X, BK1CS(16), XMAX, XMIN, XSML, Y,
     1  DCSEVLs, DBESI1s, DBSK1Es
      SAVE BK1 CS, NTK1, XMIN, XSML, XMAX
      DATA BK1 CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1 CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1 CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1 CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1 CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1 CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1 CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1 CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1 CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1 CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1 CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1 CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1 CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1 CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1 CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1 CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
      DATA NTK1, XMIN, XSML, XMAX / 0, 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESK1s
      IF (NTK1.NE.0) GO TO 10
      NTK1 = INITDSs (BK1CS, 16, 1.4e-18)
c      XMIN = DEXP (DMAX1(DLOG(D1MACH(1)), -DLOG(D1MACH(2))) + 0.01D0)
        xmin=5.94d-39
        xsml=7.5e-9
c      XSML = DSQRT (4.0D0*D1MACH(3))
c      XMAX = -DLOG(D1MACH(1))
c      XMAX = XMAX - 0.5D0*XMAX*DLOG(XMAX)/(XMAX+0.5D0)
        xmax=86.4
C
 10   IF (X.LE.0.D0) go to 1001
      IF (X.GT.2.0D0) GO TO 20
C
      IF (X.LT.XMIN) write(6,*)'DBESK1s  X SO SMALL K1 OVERFLOWS'
        Y=0.d0
      IF (X.GT.XSML) Y = X*X
      DBESK1s = DLOG(0.5D0*X)*DBESI1s(X) +
     & (0.75D0 + DCSEVLs(.5D0*Y-1.D0,
     1  BK1CS, NTK1))/X
      RETURN
C
 20   DBESK1s = 0.D0
      IF (X.GT.XMAX) write(6,*)'DBESK1s  X SO BIG K1 UNDERFLOWS'
      IF (X.GT.XMAX) RETURN
C
      DBESK1s = DEXP(-X) * DBSK1Es(X)
C
      RETURN
 1001 write(6,*)'DBESK1s  X IS ZERO OR NEGATIVE- stopped'
        stop
      END
        subroutine BESSKn(M,N,X,y,KODE)
        implicit double precision(a-h,o-z)
c Calculates an M-member sequence of modified Bessel functions of the
c  2nd kind.  Needs subroutines for exponentially scaled and unscaled
c  K_0(x) and K1(x).
c Sequence contains either unscaled or scaled Bessel functions.
c  KODE=1 means unscaled.
c  KODE=2 means scaled.
c  y(k), k=1,2,..m contains K_N (x),K_N+1 (x),...K_N+M-1 (x)
c  N must be > or = 0.
        dimension y(m)
      IF (N.LT.0) go to 1001
       if(kode.eq.2) bkm=dbsk0es(x)
       if(kode.eq.2) bk=dbsk1es(x)
       if(kode.eq.1) bkm=dbesk0s(x)
       if(kode.eq.1) bk=dbesk1s(x)
       if(n.gt.1) go to 64
       if(n.gt.0) go to 63
c  n=0
       y(1)=bkm
       if(m.lt.2) return
       y(2)=bk
       if(m.lt.3) return
       lmin=3
       go to 64
c  n=1
   63 y(1)=bk
       if(m.lt.2) return
       lmin=2
   64 n2=n+m-1
        if(n.gt.1) lmin=1
       tox=2.d0/x
       do 61  j=1,n2-1
       bkp=bkm+dfloat(j)*tox*bk
       bkm=bk
       bk=bkp
       do 62 l=lmin,m
       jj=n+l-2
   62 if(jj.eq.j) y(l)=bk
   61 continue
        return
 1001 write(6,*) 'Stopped in BESSKn- N<0:  N=',N
        stop
        end
