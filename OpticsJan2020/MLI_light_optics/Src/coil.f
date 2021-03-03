c
c
c  A few changes made (like changing .0 to 0.d0) by P. Walstrom 6/04
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine coil(pa)
c
c  Parse multipole input:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use acceldata
       include 'impli.inc'
       include 'multipole.inc'
       include 'parset.inc'
cryneneriwalstrom       include 'elmnts.inc'
        parameter(maxgauss=100)
        common/gauss/xgn(maxgauss,maxgauss),wgn(maxgauss,maxgauss)
c  Gaussian nodes and weights- used with some thick coils
        dimension xg(maxgauss),wg(maxgauss)
       double precision pa(*)
c    F. Neri  3/23/90  
c    Revised  7/12/90
c    Version  6/26/91
c    Including spacers...
c    F. Neri
       integer np
ctm
ctm    dump input parameters
ctm
cryne 8/5/2004       write(6,16) (pa(i),i=1,5)
cryne 8/5/2004  16   format('  coil: ',5f12.5)
cryne 8/5/2004       write(6,17) (pa(i),i=6,11)
cryne 8/5/2004  17   format(' shape: ',6f10.4)
       it = nint(pa(5))
c
c If itype == 0 initilize commons bincof,coeffs, and gauss
ctm /gauss/ initialization added for types 13 through 17
c
       if ( it .eq. 0 ) then
          ncoil = 0
          ztotal = 0.0d0
          call bincof
          call coeffs
c
c  Precompute gaussian weights and nodes on [-1,1]-
c  these are scaled linearily for other intervals..
c
          x1=-1.d0
          x2=1.d0
          do 13 ng=1,maxgauss
            call gauleg(x1,x2,xg,wg,ng)
          do 13 n=1,ng
            xgn(n,ng)=xg(n)
            wgn(n,ng)=wg(n)
  13      continue
          return
       endif
c
c If itype == -1 spacer.
c
       if ( it .eq. -1 ) then
          if ( ncoil .eq. 0 ) then
            xldr = pa(1)
          else
            ztotal = ztotal + pa(1)
          endif
          return
       endif
c  Increase coil number
       ncoil = ncoil + 1
       if ( ncoil .gt. maxcoils ) then
          write (6,*) ' Exceeded maximum number of Coils!'
          call myexit
       endif
c  Stack name:
       lblcoil(ncoil) = lmnlbl(inmenu)
c  Set Length, Strength and Multipole Number:       
       alcoil(ncoil) = pa(1)
       glprod(ncoil) = -pa(2)*pa(1)
       mcoil(ncoil) = nint(pa(3))
       acoil(ncoil) = pa(4)
c   Set shape(*,1)
c   For type 1: Quartic = slope
c       type 2: Halback = a2
c       type 3: Lambertson = xlmin
c       type 4  Lamb1rsth  = xlmin
c       type 5  Square = N.A.
c       type 6  User = ifile
c       type 7  FlatTop = wflat
c       type 8  Thin Halbach = N.A.
c
ctm :  load all 6 shape components from pa(6) through pa(11) (4/99)
       do 20 j=1,6
         shape(ncoil,j) = pa(5+j)
  20   continue
c
c  Increase length of String:
       ztotal = ztotal + pa(1)
c  Store position of center of coil:   
       zcoil(ncoil) = ztotal - pa(1)/2.d0
c  Get Type Number:       
       itype(ncoil) = nint(pa(5))
c   Move to position of next coil:
c       ztotal = ztotal + pa(6)  NOW done by spacers
c       
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Plot multipole configuration::::
c   Option -1 in integ.
c   Tlie output in file 19(?)
c   F. Neri 7/28/91
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine mulplt(zint)
       use acceldata
       include 'impli.inc'
       include 'files.inc'
       include 'multipole.inc'
cryneneriwalstrom       include 'elmnts.inc'
       character*11 names(10)
c
       character*11 types(10)
c
c  F. Neri 3/28/90
c  Names output 7/26/91
c  F. Neri
c  Tlie Output 7/27/91
c
       names(1) =  ' dipole'
       names(2) =  ' quadrupole'
       names(3) =  ' sextupole'
       names(4) =  ' octupole'
c
c       type 1: Quartic
c       type 2: Halback
c       type 3: Lambertson
c       type 4  Lamb1rsth
c       type 5  Square
c       type 6  User
c       type 7  FlatTop
c       type 8  Thin Halbach
c
       types(1) = ' Quartic'
       types(2) = ' RECM'
       types(3) = ' Lambertson'
       types(4) = ' Lamb1sth'
       types(5) = ' Square'
       types(6) = ' User'
       types(7) = ' Flattop'
       types(8) = ' Halbach'
c
       zi = -xldr
       zf = zint-xldr
       write(jof , 101) zi, zf 
       write(jodf, 101) zi, zf
  101  format(1x,' Integration from ',f16.8, ' to ', f16.8)
       do 2 nc = 1, ncoil
          mm = mcoil(nc)
          z0 = zcoil(nc)
          ar = acoil(nc)
          al = alcoil(nc)
          z1 = z0 - (al/2.d0)
          z2 = z0 + (al/2.d0)
          if (mm .lt. 0 ) then
            mm = -mm
            write(jof ,*) lblcoil(nc)
            write(jodf,*) lblcoil(nc)
c            
            write(jof , 102) names(mm), ar
            write(jodf, 102) names(mm), ar
  102       format(' Skew',a11,' of radius ', f16.8, ' from')
            write(jof , 103) z1, z2
            write(jodf, 103) z1, z2
  103       format(' z = ',f16.8, ' to  z = ', f16.8)
c Write Tlie file
            write(19,*) lblcoil(nc),': multipole, ',types(itype(nc)),','
            write(19,*) ' L = ',al,', Skew, K = ', -glprod(nc)/al,','
            write(19,*) ' position = ',z1,','
            write(19,*) ' radius = ',ar,', m = ',mm
          else
            write(jof ,*) lblcoil(nc)
            write(jodf,*) lblcoil(nc)
c
            write(jof , 104) names(mm), ar
            write(jodf, 104) names(mm), ar
  104       format(' ',a11,' of radius ', f16.8, ' from')
c Write Tlie file
            write(jof , 103) z1, z2
            write(jodf, 103) z1, z2
            write(19,*) lblcoil(nc),': multipole, ',types(itype(nc)),','
            write(19,*) ' L = ',al,', K = ', -glprod(nc)/al,','
            write(19,*) ' position = ',z1,','
            write(19,*) ' radius = ',ar,', m = ',mm
c
          endif
  2      continue
         write(19,*) ' : multipolelist =(',(lblcoil(nc),',',nc=1,ncoil)
  4      continue
         write(19,*) ' zmin = ',zi,', zmax = ',zf,')'
         return
         end        
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine multcfq(zlength,fa,fm)
c
c  option 0 in integ: translate coils to cfqd.
c  this version by F. Neri, 6/27/91
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
cryneneriwalstrom      include 'parm.inc'
cryneneriwalstrom      include 'param.inc'
      include 'dip.inc'
      include 'pie.inc'
      include 'gronax.inc'
      include 'multipole.inc'
c  calling arrays
      dimension fa(monoms), fm(6,6)
c  internal arrays
      parameter (MAXX = 1000)
      parameter (TINY=1.d-10)
      dimension xv(MAXX)
      dimension a(100), b(100)
c
      dimension qa(monoms), qm(6,6)
c
      xv(1) = 0.d0
      nx = 1
c
      nx = nx+1
      if( nx.gt. MAXX) goto 111
      xv(nx) =  zlength
c
      do 1 j = 1, ncoil
        do 12 i = 1,nx
          if(abs(xldr+zcoil(j)-alcoil(j)/2.d0-xv(i)).lt.TINY) goto 100
 12     continue
c
        nx = nx+1
        if( nx.gt. MAXX) goto 111
        xv(nx) = xldr+zcoil(j)-alcoil(j)/2.d0
c
 100    continue
c
        do 22 i = 1,nx
          if(abs(xldr+zcoil(j)+alcoil(j)/2.d0-xv(i)).lt.TINY) goto 200
 22     continue
c
        nx = nx+1
        if( nx.gt. MAXX) goto 111
        xv(nx) = xldr+zcoil(j)+alcoil(j)/2.d0
c
 200    continue
 1    continue
c
c  Sort xv
c
      call piksrt(nx,xv)
c
      call clear(fa,fm)
c
      do 444 j = 1,6
        fm(j,j) = 1.d0
 444  continue
c

      do 3 i = 1, nx-1
        alen = xv(i+1)-xv(i)
        xpos = (xv(i+1)+xv(i))/2.d0
        do 32 m = 1, 4
          a(m) = 0.0d0
          b(m) = 0.0d0
 32     continue
c
        do 31 j = 1, ncoil
          if(xpos.gt.zcoil(j)-alcoil(j)/2.d0+xldr
     #  .and.xpos.lt.zcoil(j)+alcoil(j)/2.d0+xldr) then
            if (mcoil(j) .gt. 0 ) then
              b(mcoil(j)) = b(mcoil(j))-glprod(j)/alcoil(j)
            else
              a(-mcoil(j)) = a(-mcoil(j))-glprod(j)/alcoil(j)
            endif
          endif
 31     continue
         write(6,113) alen,b(2),b(3),b(4)
 113     format(1x,4(1x,1pg16.8))
        call gcfqd(alen, b, a, qa, qm, 1, 1)
c        call pcmap(3, 0, 0 ,0, qa, qm) 
        call concat(fa,fm,qa,qm,fa,fm)
 3    continue
      return
 111  write(6,*) ' MAXX exceeded in integ CFQD interpreter!'
      call myexit
      return
      end
c
c
c
      SUBROUTINE PIKSRT(N,ARR)
      include 'impli.inc'
      DIMENSION ARR(N)
      DO 12 J=2,N
        A=ARR(J)
        DO 11 I=J-1,1,-1
          IF(ARR(I).LE.A)GO TO 10
          ARR(I+1)=ARR(I)
11      CONTINUE
        I=0
10      ARR(I+1)=A
12    CONTINUE
      RETURN
      END
c
c
      subroutine a3(q,x0,y0)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
cryneneriwalstrom      include 'parm.inc'
cryneneriwalstrom      include 'param.inc'
      include 'vecpot.inc'
      include 'gronax.inc'
      include 'prodex.inc'
c
      double precision AAx(0:10,0:10),AAy(0:10,0:10),AAz(0:10,0:10)
      double precision a(0:10,0:10)
c
      do 101 i = 0,10
        do 101 j = 0,10
          AAx(i,j) = 0.0d0
          AAy(i,j) = 0.0d0
          AAz(i,j) = 0.0d0
  101  continue
c
      AAz(1,0)=-(q*gn(1,0))
      AAz(0,1)=q*gs(1,0)
      AAx(2,0)=q*gn(1,1)
      AAz(2,0)=-(q*gn(2,0))
      AAx(1,1)=-(q*gs(1,1))
      AAy(1,1)=q*gn(1,1)
      AAz(1,1)=2*q*gs(2,0)
      AAy(0,2)=-(q*gs(1,1))
      AAz(0,2)=q*gn(2,0)
      AAx(3,0)=q*gn(2,1)/2
      AAz(3,0)=3*q*gn(1,2)/8 - q*gn(3,0)
      AAx(2,1)=-(q*gs(2,1))
      AAy(2,1)=q*gn(2,1)/2
      AAz(2,1)=-3*q*gs(1,2)/8 + 3*q*gs(3,0)
      AAx(1,2)=-(q*gn(2,1))/2
      AAy(1,2)=-(q*gs(2,1))
      AAz(1,2)=3*q*gn(1,2)/8 + 3*q*gn(3,0)
      AAy(0,3)=-(q*gn(2,1))/2
      AAz(0,3)=-3*q*gs(1,2)/8 - q*gs(3,0)
      AAx(4,0)=-(q*gn(1,3))/8 + q*gn(3,1)/3
      AAz(4,0)=q*gn(2,2)/6 - q*gn(4,0)
      AAx(3,1)=q*gs(1,3)/8 - q*gs(3,1)
      AAy(3,1)=-(q*gn(1,3))/8 + q*gn(3,1)/3
      AAz(3,1)=-(q*gs(2,2))/3 + 4*q*gs(4,0)
      AAx(2,2)=-(q*gn(1,3))/8 - q*gn(3,1)
      AAy(2,2)=q*gs(1,3)/8 - q*gs(3,1)
      AAz(2,2)=6*q*gn(4,0)
      AAx(1,3)=q*gs(1,3)/8 + q*gs(3,1)/3
      AAy(1,3)=-(q*gn(1,3))/8 - q*gn(3,1)
      AAz(1,3)=-(q*gs(2,2))/3 - 4*q*gs(4,0)
      AAy(0,4)=q*gs(1,3)/8 + q*gs(3,1)/3
      AAz(0,4)=-(q*gn(2,2))/6 - q*gn(4,0)
      AAx(5,0)=-(q*gn(2,3))/24 + q*gn(4,1)/4
      AAz(5,0)=-5*q*gn(1,4)/192 + 5*q*gn(3,2)/48 - q*gn(5,0)
      AAx(4,1)=q*gs(2,3)/12 - q*gs(4,1)
      AAy(4,1)=-(q*gn(2,3))/24 + q*gn(4,1)/4
      AAz(4,1)=5*q*gs(1,4)/192 - 5*q*gs(3,2)/16 + 5*q*gs(5,0)
      AAx(3,2)=-3*q*gn(4,1)/2
      AAy(3,2)=q*gs(2,3)/12 - q*gs(4,1)
      AAz(3,2)=-5*q*gn(1,4)/96 - 5*q*gn(3,2)/24 + 10*q*gn(5,0)
      AAx(2,3)=q*gs(2,3)/12 + q*gs(4,1)
      AAy(2,3)=-3*q*gn(4,1)/2
      AAz(2,3)=5*q*gs(1,4)/96 - 5*q*gs(3,2)/24 - 10*q*gs(5,0)
      AAx(1,4)=q*gn(2,3)/24 + q*gn(4,1)/4
      AAy(1,4)=q*gs(2,3)/12 + q*gs(4,1)
      AAz(1,4)=-5*q*gn(1,4)/192 - 5*q*gn(3,2)/16 - 5*q*gn(5,0)
      AAy(0,5)=q*gn(2,3)/24 + q*gn(4,1)/4
      AAz(0,5)=5*q*gs(1,4)/192 + 5*q*gs(3,2)/48 + q*gs(5,0)
      AAx(6,0)=q*gn(1,5)/192 - q*gn(3,3)/48 + q*gn(5,1)/5
      AAz(6,0)=-(q*gn(2,4))/128 + 3*q*gn(4,2)/40 - q*gn(6,0)
      AAx(5,1)=-(q*gs(1,5))/192 + q*gs(3,3)/16 - q*gs(5,1)
      AAy(5,1)=q*gn(1,5)/192 - q*gn(3,3)/48 + q*gn(5,1)/5
      AAz(5,1)=q*gs(2,4)/64 - 3*q*gs(4,2)/10 + 6*q*gs(6,0)
      AAx(4,2)=q*gn(1,5)/96 + q*gn(3,3)/24 - 2*q*gn(5,1)
      AAy(4,2)=-(q*gs(1,5))/192 + q*gs(3,3)/16 - q*gs(5,1)
      AAz(4,2)=-(q*gn(2,4))/128 - 3*q*gn(4,2)/8 + 15*q*gn(6,0)
      AAx(3,3)=-(q*gs(1,5))/96 + q*gs(3,3)/24 + 2*q*gs(5,1)
      AAy(3,3)=q*gn(1,5)/96 + q*gn(3,3)/24 - 2*q*gn(5,1)
      AAz(3,3)=q*gs(2,4)/32 - 20*q*gs(6,0)
      AAx(2,4)=q*gn(1,5)/192 + q*gn(3,3)/16 + q*gn(5,1)
      AAy(2,4)=-(q*gs(1,5))/96 + q*gs(3,3)/24 + 2*q*gs(5,1)
      AAz(2,4)=q*gn(2,4)/128 - 3*q*gn(4,2)/8 - 15*q*gn(6,0)
      AAx(1,5)=-(q*gs(1,5))/192 - q*gs(3,3)/48 - q*gs(5,1)/5
      AAy(1,5)=q*gn(1,5)/192 + q*gn(3,3)/16 + q*gn(5,1)
      AAz(1,5)=q*gs(2,4)/64 + 3*q*gs(4,2)/10 + 6*q*gs(6,0)
      AAy(0,6)=-(q*gs(1,5))/192 - q*gs(3,3)/48 - q*gs(5,1)/5
      AAz(0,6)=q*gn(2,4)/128 + 3*q*gn(4,2)/40 + q*gn(6,0)
      AAx(7,0)=q*gn(2,5)/768 - q*gn(4,3)/80 + q*gn(6,1)/6
      AAz(7,0)=7*q*gn(1,6)/9216 - 7*q*gn(3,4)/1920 + 
     -  7*q*gn(5,2)/120
      AAx(6,1)=-(q*gs(2,5))/384 + q*gs(4,3)/20 - q*gs(6,1)
      AAy(6,1)=q*gn(2,5)/768 - q*gn(4,3)/80 + q*gn(6,1)/6
      AAz(6,1)=-7*q*gs(1,6)/9216 + 7*q*gs(3,4)/640 - 
     -  7*q*gs(5,2)/24
      AAx(5,2)=q*gn(2,5)/768 + q*gn(4,3)/16 - 5*q*gn(6,1)/2
      AAy(5,2)=-(q*gs(2,5))/384 + q*gs(4,3)/20 - q*gs(6,1)
      AAz(5,2)=7*q*gn(1,6)/3072 + 7*q*gn(3,4)/1920 - 
     -  21*q*gn(5,2)/40
      AAx(4,3)=-(q*gs(2,5))/192 + 10*q*gs(6,1)/3
      AAy(4,3)=q*gn(2,5)/768 + q*gn(4,3)/16 - 5*q*gn(6,1)/2
      AAz(4,3)=-7*q*gs(1,6)/3072 + 7*q*gs(3,4)/384 + 
     -  7*q*gs(5,2)/24
      AAx(3,4)=-(q*gn(2,5))/768 + q*gn(4,3)/16 + 5*q*gn(6,1)/2
      AAy(3,4)=-(q*gs(2,5))/192 + 10*q*gs(6,1)/3
      AAz(3,4)=7*q*gn(1,6)/3072 + 7*q*gn(3,4)/384 - 
     -  7*q*gn(5,2)/24
      AAx(2,5)=-(q*gs(2,5))/384 - q*gs(4,3)/20 - q*gs(6,1)
      AAy(2,5)=-(q*gn(2,5))/768 + q*gn(4,3)/16 + 5*q*gn(6,1)/2
      AAz(2,5)=-7*q*gs(1,6)/3072 + 7*q*gs(3,4)/1920 + 
     -  21*q*gs(5,2)/40
      AAx(1,6)=-(q*gn(2,5))/768 - q*gn(4,3)/80 - q*gn(6,1)/6
      AAy(1,6)=-(q*gs(2,5))/384 - q*gs(4,3)/20 - q*gs(6,1)
      AAz(1,6)=7*q*gn(1,6)/9216 + 7*q*gn(3,4)/640 + 
     -  7*q*gn(5,2)/24
      AAy(0,7)=-(q*gn(2,5))/768 - q*gn(4,3)/80 - q*gn(6,1)/6
      AAz(0,7)=-7*q*gs(1,6)/9216 - 7*q*gs(3,4)/1920 - 
     -  7*q*gs(5,2)/120
      AAx(8,0)=-(q*gn(1,7))/9216 + q*gn(3,5)/1920 - q*gn(5,3)/120
      AAz(8,0)=q*gn(2,6)/5760 - q*gn(4,4)/480 + q*gn(6,2)/21
      AAx(7,1)=q*gs(1,7)/9216 - q*gs(3,5)/640 + q*gs(5,3)/24
      AAy(7,1)=-(q*gn(1,7))/9216 + q*gn(3,5)/1920 - q*gn(5,3)/120
      AAz(7,1)=-(q*gs(2,6))/2880 + q*gs(4,4)/120 - 2*q*gs(6,2)/7
      AAx(6,2)=-(q*gn(1,7))/3072 - q*gn(3,5)/1920 + 
     -  3*q*gn(5,3)/40
      AAy(6,2)=q*gs(1,7)/9216 - q*gs(3,5)/640 + q*gs(5,3)/24
      AAz(6,2)=q*gn(2,6)/2880 + q*gn(4,4)/120 - 2*q*gn(6,2)/3
      AAx(5,3)=q*gs(1,7)/3072 - q*gs(3,5)/384 - q*gs(5,3)/24
      AAy(5,3)=-(q*gn(1,7))/3072 - q*gn(3,5)/1920 + 
     -  3*q*gn(5,3)/40
      AAz(5,3)=-(q*gs(2,6))/960 + q*gs(4,4)/120 + 2*q*gs(6,2)/3
      AAx(4,4)=-(q*gn(1,7))/3072 - q*gn(3,5)/384 + q*gn(5,3)/24
      AAy(4,4)=q*gs(1,7)/3072 - q*gs(3,5)/384 - q*gs(5,3)/24
      AAz(4,4)=q*gn(4,4)/48
      AAx(3,5)=q*gs(1,7)/3072 - q*gs(3,5)/1920 - 3*q*gs(5,3)/40
      AAy(3,5)=-(q*gn(1,7))/3072 - q*gn(3,5)/384 + q*gn(5,3)/24
      AAz(3,5)=-(q*gs(2,6))/960 - q*gs(4,4)/120 + 2*q*gs(6,2)/3
      AAx(2,6)=-(q*gn(1,7))/9216 - q*gn(3,5)/640 - q*gn(5,3)/24
      AAy(2,6)=q*gs(1,7)/3072 - q*gs(3,5)/1920 - 3*q*gs(5,3)/40
      AAz(2,6)=-(q*gn(2,6))/2880 + q*gn(4,4)/120 + 2*q*gn(6,2)/3
      AAx(1,7)=q*gs(1,7)/9216 + q*gs(3,5)/1920 + q*gs(5,3)/120
      AAy(1,7)=-(q*gn(1,7))/9216 - q*gn(3,5)/640 - q*gn(5,3)/24
      AAz(1,7)=-(q*gs(2,6))/2880 - q*gs(4,4)/120 - 2*q*gs(6,2)/7
      AAy(0,8)=q*gs(1,7)/9216 + q*gs(3,5)/1920 + q*gs(5,3)/120
      AAz(0,8)=-(q*gn(2,6))/5760 - q*gn(4,4)/480 - q*gn(6,2)/21
      AAx(9,0)=-(q*gn(2,7))/46080 + q*gn(4,5)/3840 - 
     -  q*gn(6,3)/168
      AAz(9,0)=-(q*gn(1,8))/81920 + q*gn(3,6)/15360 - 
     -  3*q*gn(5,4)/2240
      AAx(8,1)=q*gs(2,7)/23040 - q*gs(4,5)/960 + q*gs(6,3)/28
      AAy(8,1)=-(q*gn(2,7))/46080 + q*gn(4,5)/3840 - 
     -  q*gn(6,3)/168
      AAz(8,1)=q*gs(1,8)/81920 - q*gs(3,6)/5120 + 3*q*gs(5,4)/448
      AAx(7,2)=-(q*gn(2,7))/23040 - q*gn(4,5)/960 + q*gn(6,3)/12
      AAy(7,2)=q*gs(2,7)/23040 - q*gs(4,5)/960 + q*gs(6,3)/28
      AAz(7,2)=-(q*gn(1,8))/20480 + 3*q*gn(5,4)/280
      AAx(6,3)=q*gs(2,7)/7680 - q*gs(4,5)/960 - q*gs(6,3)/12
      AAy(6,3)=-(q*gn(2,7))/23040 - q*gn(4,5)/960 + q*gn(6,3)/12
      AAz(6,3)=q*gs(1,8)/20480 - q*gs(3,6)/1920
      AAx(5,4)=-(q*gn(4,5))/384
      AAy(5,4)=q*gs(2,7)/7680 - q*gs(4,5)/960 - q*gs(6,3)/12
      AAz(5,4)=-3*q*gn(1,8)/40960 - q*gn(3,6)/2560 + 
     -  3*q*gn(5,4)/160
      AAx(4,5)=q*gs(2,7)/7680 + q*gs(4,5)/960 - q*gs(6,3)/12
      AAy(4,5)=-(q*gn(4,5))/384
      AAz(4,5)=3*q*gs(1,8)/40960 - q*gs(3,6)/2560 - 
     -  3*q*gs(5,4)/160
      AAx(3,6)=q*gn(2,7)/23040 - q*gn(4,5)/960 - q*gn(6,3)/12
      AAy(3,6)=q*gs(2,7)/7680 + q*gs(4,5)/960 - q*gs(6,3)/12
      AAz(3,6)=-(q*gn(1,8))/20480 - q*gn(3,6)/1920
      AAx(2,7)=q*gs(2,7)/23040 + q*gs(4,5)/960 + q*gs(6,3)/28
      AAy(2,7)=q*gn(2,7)/23040 - q*gn(4,5)/960 - q*gn(6,3)/12
      AAz(2,7)=q*gs(1,8)/20480 - 3*q*gs(5,4)/280
      AAx(1,8)=q*gn(2,7)/46080 + q*gn(4,5)/3840 + q*gn(6,3)/168
      AAy(1,8)=q*gs(2,7)/23040 + q*gs(4,5)/960 + q*gs(6,3)/28
      AAz(1,8)=-(q*gn(1,8))/81920 - q*gn(3,6)/5120 - 
     -  3*q*gn(5,4)/448
      AAy(0,9)=q*gn(2,7)/46080 + q*gn(4,5)/3840 + q*gn(6,3)/168
      AAz(0,9)=q*gs(1,8)/81920 + q*gs(3,6)/15360 + 
     -  3*q*gs(5,4)/2240
      AAx(10,0)=q*gn(1,9)/737280 - q*gn(3,7)/138240 + 
     -  q*gn(5,5)/6720
      AAz(10,0)=-(q*gn(2,8))/442368 + q*gn(4,6)/32256 - 
     -  5*q*gn(6,4)/5376
      AAx(9,1)=-(q*gs(1,9))/737280 + q*gs(3,7)/46080 - 
     -  q*gs(5,5)/1344
      AAy(9,1)=q*gn(1,9)/737280 - q*gn(3,7)/138240 + 
     -  q*gn(5,5)/6720
      AAz(9,1)=q*gs(2,8)/221184 - q*gs(4,6)/8064 + 
     -  5*q*gs(6,4)/896
      AAx(8,2)=q*gn(1,9)/184320 - q*gn(5,5)/840
      AAy(8,2)=-(q*gs(1,9))/737280 + q*gs(3,7)/46080 - 
     -  q*gs(5,5)/1344
      AAz(8,2)=-(q*gn(2,8))/147456 - q*gn(4,6)/10752 + 
     -  65*q*gn(6,4)/5376
      AAx(7,3)=-(q*gs(1,9))/184320 + q*gs(3,7)/17280
      AAy(7,3)=q*gn(1,9)/184320 - q*gn(5,5)/840
      AAz(7,3)=q*gs(2,8)/55296 - q*gs(4,6)/4032 - 5*q*gs(6,4)/672
      AAx(6,4)=q*gn(1,9)/122880 + q*gn(3,7)/23040 - q*gn(5,5)/480
      AAy(6,4)=-(q*gs(1,9))/184320 + q*gs(3,7)/17280
      AAz(6,4)=-(q*gn(2,8))/221184 - q*gn(4,6)/2304 + 
     -  5*q*gn(6,4)/384
      AAx(5,5)=-(q*gs(1,9))/122880 + q*gs(3,7)/23040 + 
     -  q*gs(5,5)/480
      AAy(5,5)=q*gn(1,9)/122880 + q*gn(3,7)/23040 - q*gn(5,5)/480
      AAz(5,5)=q*gs(2,8)/36864 - 5*q*gs(6,4)/192
      AAx(4,6)=q*gn(1,9)/184320 + q*gn(3,7)/17280
      AAy(4,6)=-(q*gs(1,9))/122880 + q*gs(3,7)/23040 + 
     -  q*gs(5,5)/480
      AAz(4,6)=q*gn(2,8)/221184 - q*gn(4,6)/2304 - 
     -  5*q*gn(6,4)/384
      AAx(3,7)=-(q*gs(1,9))/184320 + q*gs(5,5)/840
      AAy(3,7)=q*gn(1,9)/184320 + q*gn(3,7)/17280
      AAz(3,7)=q*gs(2,8)/55296 + q*gs(4,6)/4032 - 5*q*gs(6,4)/672
      AAx(2,8)=q*gn(1,9)/737280 + q*gn(3,7)/46080 + 
     -  q*gn(5,5)/1344
      AAy(2,8)=-(q*gs(1,9))/184320 + q*gs(5,5)/840
      AAz(2,8)=q*gn(2,8)/147456 - q*gn(4,6)/10752 - 
     -  65*q*gn(6,4)/5376
      AAx(1,9)=-(q*gs(1,9))/737280 - q*gs(3,7)/138240 - 
     -  q*gs(5,5)/6720
      AAy(1,9)=q*gn(1,9)/737280 + q*gn(3,7)/46080 + 
     -  q*gn(5,5)/1344
      AAz(1,9)=q*gs(2,8)/221184 + q*gs(4,6)/8064 + 
     -  5*q*gs(6,4)/896
      AAy(0,10)=-(q*gs(1,9))/737280 - q*gs(3,7)/138240 - 
     -  q*gs(5,5)/6720
      AAz(0,10)=q*gn(2,8)/442368 + q*gn(4,6)/32256 + 
     -  5*q*gn(6,4)/5376
c
cryneneriwalstrom 16 July, 2004 changed "4" to maxorder
cryneneriwalstrom and set maxorder to 6
      maxorder=6
c
      call exp2f1(x0,y0,AAx,a,10)
      ind1 = 0
      do 1 i = 0, maxorder
        ind2 = ind1
        do 2 j = 0, maxorder-i
          Ax(ind2) = a(i,j)
          if (j .lt. maxorder-i) then
            ind2 = prodex(3,ind2)
          endif
 2      continue
        if ( i .lt. maxorder ) ind1 = prodex(1,ind1)
 1    continue
c
      call exp2f1(x0,y0,AAy,a,10)
      ind1 = 0
      do 10 i = 0, maxorder
        ind2 = ind1
        do 20 j = 0, maxorder-i
          Ay(ind2) = a(i,j)
          if (j .lt. maxorder-i) then
            ind2 = prodex(3,ind2)
          endif
 20     continue
        if ( i .lt. maxorder ) ind1 = prodex(1,ind1)
 10   continue
c
      call exp2f1(x0,y0,AAz,a,10)
      ind1 = 0
      do 100 i = 0, maxorder
        ind2 = ind1
        do 200 j = 0, maxorder-i
          Az(ind2) = a(i,j)
          if (j .lt. maxorder-i) then
            ind2 = prodex(3,ind2)
          endif
 200    continue
        if ( i .lt. maxorder ) ind1 = prodex(1,ind1)
 100  continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Reexpands a Polynomial in x and y where the coefficient
c  of x^i y^j is a(i,j), and a() is dimensioned as a(0:md,0:md),
c  around the point x0, y0.
c  The coefficients of the reexpanded polynomial are in the array
c  b(0:md,0:md).
c
c  F. Neri 3/29/91.
c
      subroutine exp2f1(x0,y0,a,b,md)
      implicit double precision (a-h,o-z)
      parameter (maxd = 100)
      dimension a(0:md,0:md), b(0:md,0:md)
c
      dimension c(0:maxd,0:maxd), d(0:maxd), e(0:maxd)
c
      do 1 i = 0, md
        call expoly(a(0,i), md, x0, c(0,i), md)
 1    continue
c
      do 2 i = 0, md
        do 20 j = 0, md
          d(j) = c(i,j)
 20     continue
        call expoly(d, md, y0 ,e, md)
        do 200 j = 0, md
          b(i,j) = e(j)
 200    continue
 2      continue
        return
        end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE EXPOLY(C,NC,X,PD,ND)
c
c Given the NC+1 coefficients of a polynomial of degree NC as an array C
c with C(0) being the constant term, and a given a value X, and given a
c value ND>0, this routine returns the polynomial evaluated at X as PD(0)
c and ND derivatives times ND!, as PD(1)  . . . PD(ND).
c From Numerical Recipes, pag. 138 (1986 FORTRAN/Pascal version).
c
      implicit double precision (a-h,o-z)  
      DIMENSION C(0:NC),PD(0:ND)
      if ( x .ne. 0.d0 ) then
      PD(0)=C(NC)
      DO 11 J=1,ND
        PD(J)=0.d0
 11    CONTINUE
      DO 13 I=NC-1,0,-1
        NND=MIN(ND,NC-I)
        DO 12 J=NND,1,-1
          PD(J)=PD(J)*X+PD(J-1)
 12      CONTINUE
        PD(0)=PD(0)*X+C(I)
 13    CONTINUE
      else
       do 14 i = 0, nd
         pd(i) = c(i)
 14    continue
      endif
      RETURN
      END
c
c**************************************************
c
	subroutine gauleg(x1,x2,x,w,n)
	implicit double precision(a-h,o-z)
	dimension x(n),w(n)
c  Routine from Numerical Recipes to calculate Gaussian weights and node
c  for Gaussian quadrature on the interval [x1,x2].
c  x1=lower end of integration interval
c  x2=upper  "   "     "          "
c  x=array of gauss nodes
c  w=array of gauss weights
c  n=gaussian order
	parameter (eps=3.d-14)
	parameter( half=0.5d0)
	pi=2.d0*dasin(1.d0)
	m=(n+1)/2
	xm=half*(x2+x1)
	xl=half*(x2-x1)
	do 12 i=1,m
	z=dcos(pi*(i-0.25d0)/(n+half))
    1 continue
	p1=1.d0
	p2=0.d0
	do 11 j=1,n
	p3=p2
	p2=p1
	p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
   11 continue
	pp=n*(z*p1-p2)/(z*z-1.d0)
	z1=z
	z=z1-p1/pp
	if(dabs(z-z1).gt.eps) go to 1
	x(i)=xm-xl*z
	x(n+1-i)=xm+xl*z
	w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
	w(n+1-i)=w(i)
   12 continue
	return
	end
