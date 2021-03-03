c IMPACT RF Gap routines
c Copyright 2001 University of California
c
      subroutine rfgap(zedge,zmap,nstep,rffreq,escale,theta0,           &
     &t00,gam00,itype,h,mh)
cryne use Dataclass
c h(28-209) = 2nd order and 3rd order map, mh(6,6)=transfer matrix
      use parallel, only : idproc
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
cryne include 'para.inc'
      include 'map.inc'
ctm   include 'arcblk.inc'
      include 'usrdat.inc'
      double precision zedge,zmap,rffreq,escale,theta0
      integer nstep,itype,npusr
      double precision h(monoms),mh(6,6)
      common/rfcdata/zdat(20500),edat(20500),epdat(20500),Nsize
!     write(6,*)'(rfgap) ref(6),brho=',reftraj(6),brho
! nonlinear map for rf gap not implemented yet
      do 10 i=1,monoms
      h(i)=0.
   10 continue
!!!!!!!!!!!!! fix these hardwired values later:
!     qmcc=1./938.28e6
cryne 7/21/02      qmcc=1./939.294308e6
      qmcc=1.d0/pmass
c
      xl=sl
      xk=1./xl
      ww=rffreq/(c/xl)
c     write(6,*) "==> ww = ",ww
c g0 is the quadrupole gradient, assumed to be zero
c (combined rfgap+quad will be implemented later if requested)
      g0=0.
c
      if(zmap.lt.0.)then
        write(6,*)'error(rfgap): zmap cannot be < 0 in this version'
      endif
      tau=zmap
c
      engin=gam00*pmass
c     ucalc(11)=engin
      if(idproc.eq.0)write(36,*)arclen,engin
      call myflush(36)
      call mapgap(mh,arclen,tau,nstep,qmcc,xk,xl,ww,zedge,              &
     &                    escale,theta0,g0,t00,gam00,itype)
! update brho, beta, and gamma:
cryne 5/1/2006 No. don't update here. Moved to routine lmnt in afro.f
!      engfin=gam00*pmass
cryne-abell Mon Nov 24 18:39:30 PST 2003
!      ucalc(12)=engfin
      return
      end
 
      subroutine mapgap(xm,t,tau,mapstp,qmcc,xk,xl,ww,zedge,escale,     &
     &                    theta0,g0,t00,gam00,itype)
! Compute the linear map for an rf gap + quadrupole
! 2 design orbit eqns + 12 matrix coeffs = 14 eqns
! all the following arrays are serial arrays:
cryne== implicit none
      include 'impli.inc'
cryne==
        dimension xm(6,6), y(14)
ctm        integer, intent(in) :: mapstp,itype
ctm        double precision, intent(in) :: t,tau,qmcc,xk,xl,ww,zedge,         &
ctm     &                                escale,theta0,g0
ctm        double precision, intent(inout) :: t00,gam00
ctm        double precision, dimension(6,6), intent(out) :: xm
ctm        double precision, dimension(14) :: y
ctm        double precision :: h,t0,gi,betai,ui,squi,sqi3,ein,e1in,e2in,      &
ctm     &                    sinthi
ctm        double precision :: costhi,uprimi,qpwi,dlti,thli,gf,betaf,uf,      &
ctm     &                    squf,sqf3
ctm        double precision :: tfin,ef,e1f,e2f,sinthf,costhf,uprimf,qpwf,     &
ctm     &                    dltf,thlf
!       xm = 0.0
        do 11 i=1,6
        do 10 j=1,6
        xm(i,j)=0.d0
   10   continue
   11   continue
c       y(1)=mpasyn(5)
c       y(2)=mpasyn(6)
        y(1)=t00
        y(2)=-gam00
        y(3)=1.d0
        y(4)=0.d0
        y(5)=0.d0
        y(6)=1.d0
        y(7)=1.d0
        y(8)=0.d0
        y(9)=0.d0
        y(10)=1.d0
        y(11)=1.d0
        y(12)=0.d0
        y(13)=0.d0
        y(14)=1.d0
        h=tau/mapstp
        t0=t
c       gi=-mpasyn(6)
        gi=gam00
        betai=sqrt((gi-1.d0)*(gi+1.d0))/gi
        ui=gi*betai
        squi=sqrt(ui)
        sqi3=sqrt(ui)**3
        call eintrp(t,ein,e1in,e2in,zedge,escale,itype,ww,theta0,y)
        sinthi=sin(ww*t00+theta0)
        costhi=cos(ww*t00+theta0)
        uprimi=qmcc*ein/betai*costhi
        qpwi=0.5d0*qmcc*xl/(ui*ww)
        dlti=xl*(0.5d0*uprimi/ui-qpwi*e1in*sinthi)
        thli=1.5d0*xl*(uprimi/ui)
!       call rk6i(h,mapstp,t0,y,qmcc,xk,xl,ww,zedge,escale,g0,           &
!    &            theta0,itype)
        call adam11rf(h,mapstp,t0,y,qmcc,xk,xl,ww,zedge,escale,g0,       &
     &            theta0,itype)
c       mpasyn(5)=y(1)
c       mpasyn(6)=y(2)
c       gf=-mpasyn(6)
        t00=y(1)
        gam00=-y(2)
        gf=gam00
        betaf=sqrt((gf-1.d0)*(gf+1.d0))/gf
        uf=gf*betaf
        squf=sqrt(uf)
        sqf3=sqrt(uf)**3
        tfin=t+tau
        call eintrp(tfin,ef,e1f,e2f,zedge,escale,itype,ww,theta0,y)
        sinthf=sin(ww*t00+theta0)
        costhf=cos(ww*t00+theta0)
        uprimf=qmcc*ef/betaf*costhf
        qpwf=0.5d0*qmcc*xl/(uf*ww)
        dltf=xl*(0.5d0*uprimf/uf-qpwf*e1f*sinthf)
        thlf=1.5d0*xl*(uprimf/uf)
        xm(1,1)= y(3)*squi/squf+y(5)*squi/squf*dlti
        xm(2,1)=(y(4)-y(3)*dltf)*squi*squf+(y(6)-y(5)*dltf)*squi*squf*   &
     &           dlti
        xm(1,2)= y(5)/(squi*squf)
        xm(2,2)=(y(6)-y(5)*dltf)*squf/squi
        xm(3,3)= y(7)*squi/squf+y(9)*squi/squf*dlti
        xm(4,3)=                                                         &
     &  (y(8)-y(7)*dltf)*squi*squf+(y(10)-y(9)*dltf)*squi*squf*dlti
        xm(3,4)= y(9)/(squi*squf)
        xm(4,4)=(y(10)-y(9)*dltf)*squf/squi
        xm(5,5)= y(11)*sqi3/sqf3+y(13)*sqi3/sqf3*thli
        xm(6,5)=                                                         &
     &  (y(12)-y(11)*thlf)*sqi3*sqf3+(y(14)-y(13)*thlf)*sqi3*sqf3*thli
        xm(5,6)= y(13)/(sqi3*sqf3)
        xm(6,6)=(y(14)-y(13)*thlf)*sqf3/sqi3
 
        end    !subroutine mapgap
! 
      subroutine eintrp(z,ez1,ezp1,ezpp1,zedge,escale,itype,ww,theta0,y)
      use parallel
      use beamdata, only : pmass
      include 'impli.inc'
      include 'fitbuf.inc'
      dimension y(14)
      common/rfcdata/zdat(20500),edat(20500),epdat(20500),Nsize
!     if(itype.eq.0 .or. itype.eq.1)then
!       ez1=0.
!       ezp1=0.
!       ezpp1=0.
!       goto 1000
!     endif
      zz=z-zedge
      checklo=zz-zdat(1)
      checkhi=zdat(Nsize)-zz
      eps=1.d-10
      if(checklo.lt.0.d0)then
        if(abs(checklo).lt.eps)then
          zz=zdat(1)
        else
          if(idproc.eq.0)then 
          write(6,*)'interpolation error in routine eintrp'
          write(6,*)'zz,zdat(1)=',zz,zdat(1)
          endif
          call myexit
        endif
      endif
      if(checkhi.lt.0.d0)then
        if(abs(checkhi).lt.eps)then
          zz=zdat(Nsize)
        else
          if(idproc.eq.0)then 
          write(6,*)'interpolation error in routine eintrp'
          write(6,*)'zz,zdat(Nsize)=',zz,zdat(Nsize)
          endif
          call myexit
        endif
      endif
!============================
      klo=1
      khi=Nsize
    1 if (khi-klo.gt.1)then
         k=(khi+klo)/2
         if(zdat(k).gt.zz)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=zdat(khi)-zdat(klo)
      if(h.eq.0.)then
        write(6,*)'bad data in routine eintrp'
        call myexit
      endif
!============================
!start linear interpolation.
        slope1=(edat(khi)-edat(klo))/h
        ez1 =edat(klo)+slope1*(zz-zdat(klo))
        slope2=(epdat(khi)-epdat(klo))/h
        ezp1=epdat(klo)+slope2*(zz-zdat(klo))
        ez1=ez1*escale
        ezp1=ezp1*escale
        ezpp1=0.
cryne
cryne   if(idproc.eq.0) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c modification to only write when fitting has converged:
 1000 continue
cryne 5/9/2006      write(42,*) z,ez1,ezp1,klo,khi,min(zz-zdat(klo),zdat(khi)-zz)
c       if(kfit.eq.1 .and. idproc.eq.0)then
        if(idproc.eq.0)then
cccc      write(49,101)z,zz,y(1),(-y(2)-1.)*pmass,ez1,ezp1,1.0*m
cccc      call myflush(49)
          m=klo
c         if(m.eq.435 .or. m.eq.920)then
          if( m.gt.1 .and. m.lt.Nsize)then
c          if( (edat(m-1).lt.edat(m)).and.(edat(m+1).lt.edat(m)) )then
           if( ((edat(m-1).lt.edat(m)).and.(edat(m+1).lt.edat(m))).or.  &
     &         ((edat(m-1).gt.edat(m)).and.(edat(m+1).gt.edat(m))) )then
           twopi=4.*asin(1.0d0)
           phasprnt=360./twopi*mod(ww*y(1)+theta0,twopi)
           if(phasprnt.lt.10.)phasprnt=phasprnt+360.
           prxrad=360./twopi*(ww*y(1)+theta0)
           preng=(-y(2)-1.)*939.294308
           write(59,101)z,1.0*m,phasprnt,prxrad,theta0*360./twopi,preng
           call myflush(59)
           endif
          endif
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cryne   endif
  101 format(7(1x,1pe12.5))
        return
        end    !subroutine eintrp
 
        double precision function func1(x,pi52)
cryne== implicit none
      include 'impli.inc'
cryne==
ctm        double precision, intent(in) :: x, pi52
 
        func1 = exp(-(4.*x)**4)*cos(pi52*tanh(5.*x))
 
        end    !function func1
c************************************************************
        double precision function func2(x,pi52)
cryne== implicit none
      include 'impli.inc'
cryne==
ctm        double precision, intent(in) :: x, pi52
 
        func2=-5.*pi52*exp(-(4.*x)**4)*sin(pi52*tanh(5.*x))/              &
     &  cosh(5.*x)**2-16.*(4.*x)**3*exp(-(4.*x)**4)*cos(pi52*tanh(5.*x))
 
        end    !function func2
 
c***********************************************
        subroutine rk6i(h,ns,t,y,qmcc,xk,xl,ww,zedge,escale,g0,           &
     &                  theta0,itype)
cryne== implicit none
      include 'impli.inc'
cryne==
ctm        integer, intent(in) :: ns,itype
ctm        double precision, intent(in) :: qmcc,xk,xl,ww,zedge,escale,       &
ctm     &                                g0,theta0
ctm        double precision, intent(inout) :: t
ctm        double precision, intent(in) :: h
ctm        double precision, dimension(14), intent(inout) :: y
ctm        double precision, dimension(14) :: yt,a,b,c,d,e,f,g,o,p
ctm        double precision:: tint,tt,blg
ctm        integer :: i
        dimension y(14),yt(14),a(14),b(14),c(14),d(14),e(14),f(14),          &
     &  g(14),o(14),p(14)
        blg=0.
        tint=t
        do i=1,ns
          call evalrf(t,y,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,       &
     &              itype)
      do 10 j=1,14
          a(j)=h*f(j)
          yt(j)=y(j)+a(j)/9.d0
   10 continue
          tt=t+h/9.d0
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 20 j=1,14
          b(j)=h*f(j)
          yt(j)=y(j) + (a(j) + 3.d0*b(j))/24.d0
   20 continue
          tt=t+h/6.d0
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 30 j=1,14
          c(j)=h*f(j)
          yt(j)=y(j)+(a(j)-3.d0*b(j)+4.d0*c(j))/6.d0
   30 continue
          tt=t+h/3.d0
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 40 j=1,14
          d(j)=h*f(j)
          yt(j)=y(j) +                                                  &
     &(278.d0*a(j) - 945.d0*b(j) + 840.d0*c(j) + 99.d0*d(j))/544.d0
   40 continue
          tt=t+.5d0*h
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 50 j=1,14
          e(j)=h*f(j)
          yt(j)=y(j) +                                                  &
     &(-106.d0*a(j)+273.d0*b(j)-104.d0*c(j)-107.d0*d(j)+48.d0*e(j))/6.d0
   50 continue
          tt = t+2.d0*h/3.d0
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 60 j=1,14
          g(j)=h*f(j)
          yt(j) = y(j)+(110974.d0*a(j)-236799.d0*b(j)+68376.d0*c(j)+    &
     &            103803.d0*d(j)-10240.d0*e(j) + 1926.d0*g(j))/45648.d0
   60 continue
          tt = t + 5.d0*h/6.d0
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 70 j=1,14
          o(j)=h*f(j)
          yt(j) = y(j)+(-101195.d0*a(j)+222534.d0*b(j)-71988.d0*c(j)-   &
     &    26109.d0*d(j)-2.d4*e(j)-72.d0*g(j)+22824.d0*o(j))/25994.d0
   70 continue
          tt = t + h
          call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     &              itype)
      do 80 j=1,14
          p(j)=h*f(j)
          y(j) = y(j)+(41.d0*a(j)+216.d0*c(j)+27.d0*d(j)+               &
     &           272.d0*e(j)+27.d0*g(j)+216.d0*o(j)+41.d0*p(j))/840.d0
   80 continue
          t=tint+i*h
crynedabell
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2)
cryne 5/9/2006      call myflush(41)
crynedabell
!         write(20,101)t,y(1),(-y(2)-1.0)*938.28
cryne     write(20,101)t,y(1),(-y(2)-1.0)*939.294308
  101     format(3(1pe12.5,1x))
cryne     write(21,102)t,y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),   &
cryne&                 y(12)
  102     format(11(1pe12.5,1x))
        enddo
        end    !subroutine rk6i
c********************************************************************
        subroutine evalrf(t,y,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,          &
     &                  theta0,itype)
        use beamdata, only : pmass
cryne== implicit none
      include 'impli.inc'
cryne==
ctm        integer, intent(in) :: itype
ctm        double precision, intent(in) :: t,qmcc,xk,xl,ww,zedge,blg,        &
ctm     &                                escale,g0,theta0
ctm        double precision, dimension(14), intent(in) :: y
ctm        double precision, dimension(14), intent(out) :: f
ctm        double precision :: clite,ez1,ezp1,ezpp1,gamma0,beta0,            &
ctm     &                   gbet,sinphi,cosphi
ctm        double precision :: rfdsgn,pmass,brho,s11tmp,s11,s33,s55
cryne 7/21/02
c
       dimension f(14), y(14)
        clite=299792458.d0
        call eintrp(t,ez1,ezp1,ezpp1,zedge,escale,itype,ww,theta0,y)
! synchronous particle:
        gamma0=-y(2)
        beta0=sqrt((gamma0-1.d0)*(gamma0+1.d0))/gamma0
        gbet=beta0*gamma0
        sinphi=sin(ww*y(1)+theta0)
        cosphi=cos(ww*y(1)+theta0)
        rfdsgn=ez1*cosphi
        f(1)=xk/beta0
        f(2)=-qmcc*rfdsgn
cryneabell:09.Nov.05
cryne 5/9/2006        write(44,*) t,ww*y(1)+theta0,ez1,ezp1
cryneabell----------
! matrix elements
!       pmass=938.28d6
! mistaken comment was here previously; now set to H- mass   Sept 5, 2000
!       pmass=939.294308d6
!!!************************************************************************
!!!cryne 11/06/2003 eventually this should be fixed to use pmass in common**
!       pmass=938.27231d6
!!!************************************************************************
!!!************************************************************************
! put a negative sign here for H-.
!       brho=-gbet/clite*pmass
        brho=gbet/clite*pmass
        s11tmp=0.5d0*(1.d0+0.5d0*gamma0**2)*                            &
     &   (qmcc*rfdsgn/beta0**2/gamma0**2)**2 +                          &
     &    qmcc*xk*ww*0.5d0/beta0**3/gamma0**3*ez1*sinphi
!          qmcc*xk*0.5d0/beta0**3/gamma0**3*ez1*sinphi
        s11=(s11tmp+g0/brho)*xl
        s33=(s11tmp-g0/brho)*xl
        s55=                                                            &
     &   -1.5d0*qmcc/beta0**2/gamma0*ezp1*cosphi                        &
     &   +(beta0**2+0.5d0)/beta0**3/gamma0*qmcc*xk*ww*ez1*sinphi        &
     &   +1.5d0*(1.d0-0.5d0*gamma0**2)*                                 &
     &   (qmcc*rfdsgn/beta0**2/gamma0**2)**2
!         +(beta0**2+0.5d0)/beta0**3/gamma0*qmcc*xk*ez1*sinphi
        s55=xl*s55
        f(3)=y(4)/xl
        f(4)=-s11*y(3)
        f(5)=y(6)/xl
        f(6)=-s11*y(5)
        f(7)=y(8)/xl
        f(8)=-s33*y(7)
        f(9)=y(10)/xl
        f(10)=-s33*y(9)
        f(11)=y(12)/xl
        f(12)=-s55*y(11)
        f(13)=y(14)/xl
        f(14)=-s55*y(13)
cryneabell---9.Nov.05
cryne 5/9/2006      write(43,*) t,f(1)*ww,f(2)/ww,f(1:2),f(7:14)
cryneabell
        end    !subroutine evalrf
 
c==========================================================
      subroutine read_RFdata(nunit,numdata,zlen)
      use parallel
      include 'impli.inc'
      common/rfcdata/zdat(20500),edat(20500),epdat(20500),Nsize
      if(idproc.eq.0)then
        numdata=0
        do i = 1, 999999
        read(nunit,*,end=999)zdat(i),edat(i),epdat(i)
        numdata=numdata+1
        enddo
  999   continue
cryne note: this sets the left hand edge to zero:
        zdat1=zdat(1)
        do i=1,numdata
        zdat(i)=zdat(i)-zdat1
        enddo
c       if(idproc.eq.0)then
c         write(6,*)'Done reading RF data. Closing file unit ',nunit
c       endif
        close(nunit)
      endif
      call MPI_BCAST(numdata,1,mntgr,0,lworld,ierr)
      call MPI_BCAST(zdat,numdata,mreal,0,lworld,ierr)
      call MPI_BCAST(edat,numdata,mreal,0,lworld,ierr)
      call MPI_BCAST(epdat,numdata,mreal,0,lworld,ierr)
      Nsize = numdata
      zlen=zdat(numdata)-zdat(1)
      return
      end
************************************************************************
c
!     subroutine adam11rf(h,ns,nf,t,y,ne)
      subroutine adam11rf(h,ns,t,y,qmcc,xk,xl,ww,zedge,escale,g0,       &
     &            theta0,itype)
c  Written by Rob Ryne, Spring 1986, based on a routine of Alex Dragt
c  This integration routine makes local truncation errors at each
c  step of order h**11.  That is, it is locally correct through
c  order h**10.  Due to round off errors, its true precision is
c  realized only when more than 64 bits are used.
      use lieaparam, only : monoms
      use beamdata
      use phys_consts
      include 'impli.inc'
c     character*6 nf
cryneneriwalstrom fix later to use monoms instead of hardwire:
      dimension y( 14),yp( 14),yc( 14),f1( 14),f2( 14),f3( 14),f4( 14),
     & f5( 14),f6( 14),f7( 14),f8( 14),f9( 14),f10( 14),f11( 14)
      dimension a(10),am(10),b(10),bm(10)
c
      data (a(i),i=1,10)/57281.d0,-583435.d0,2687864.d0,
     & -7394032.d0,13510082.d0,-17283646.d0,16002320.d0,
     & -11271304.d0,9449717.d0,2082753.d0/
      data (b(i),i=1,10)/-2082753.d0,20884811.d0,-94307320.d0,
     & 252618224.d0,-444772162.d0,538363838.d0,-454661776.d0,
     & 265932680.d0,-104995189.d0,30277247.d0/
cryne 7/23/2002
      save a,b
c
cryneabell 11/8/2005
      ne=14
c     nf='start'
cryneabell 11/8/2005
c
cryne 1 August 2004      ne=monoms+15
c
      nsa=ns
c     if (nf.eq.'cont') go to 20
c  rk start
      iqt=5
      qt=float(iqt)
      hqt=h/qt
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f1,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f2,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f3,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f4,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f5,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f6,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f7,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f8,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f9,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,itype)
      call rk78iirf(hqt,iqt,t,y,qmcc,xk,xl,ww,zedge,blg,escale,g0,      &
     & theta0,itype)
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
      call evalrf(t,y,f10,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      nsa=ns-9
      hdiv=h/7257600.0d+00
      do 10 i=1,10
      am(i)=hdiv*a(i)
  10  bm(i)=hdiv*b(i)
  20  tint=t
      do 100 i=1,nsa
      do 30 j=1,ne
      yp(j)=y(j)+bm(1)*f1(j)+bm(2)*f2(j)+bm(3)*f3(j)                    &
     &+bm(4)*f4(j)+bm(5)*f5(j)+bm(6)*f6(j)+bm(7)*f7(j)                  &
     & +bm(8)*f8(j)+bm(9)*f9(j)+bm(10)*f10(j)
   30 continue
      call evalrf(t+h,yp,f11,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,  &
     & itype)
      do 40 j=1,ne
      yp(j)=y(j)+am(1)*f2(j)+am(2)*f3(j)+am(3)*f4(j)+am(4)*f5(j)        &
     & +am(5)*f6(j)+am(6)*f7(j)+am(7)*f8(j)+am(8)*f9(j)+am(9)*f10(j)
  40  yc(j)=yp(j)+am(10)*f11(j)
  41  call evalrf(t+h,yc,f11,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,  &
     & itype)
      do 50 j=1,ne
   50 y(j)=yp(j)+am(10)*f11(j)
      do 60 j=1,ne
      f1(j)=f2(j)
      f2(j)=f3(j)
      f3(j)=f4(j)
      f4(j)=f5(j)
      f5(j)=f6(j)
      f6(j)=f7(j)
      f7(j)=f8(j)
      f8(j)=f9(j)
      f9(j)=f10(j)
   60 f10(j)=f11(j)
      t=tint+i*h
cryne 5/9/2006      write(41,*) t,y(1:2),-y(2),-y(2)*omegascl*sl/c_light
cryne 5/9/2006      call myflush(41)
  100 continue
      return
      end
c
***********************************************************************
c
c     subroutine rk78iirf(h,ns,t,y,ne)
      subroutine rk78iirf(h,ns,t,y,qmcc,xk,xl,ww,zedge,blg,escale,      &
     & g0,theta0,itype)
c  Written by Rob Ryne, Spring 1986, based on a routine of
c  J. Milutinovic.
c  For a reference, see page 76 of F. Ceschino and J Kuntzmann,
c  Numerical Solution of Initial Value Problems, Prentice Hall 1966.
c  This integration routine makes local truncation errors at each
c  step of order h**7.
c  That is, it is locally correct through terms of order h**6.
c  Each step requires 8 function evaluations.

      use lieaparam, only : monoms
      include 'impli.inc'
cryneneriwalstrom fix later to use monoms instead of hardwire:
      dimension y( 14),yt( 14),f( 14),a( 14),b( 14),c( 14),d( 14),      &
     &e( 14),g( 14),o( 14),p( 14)
cryne 1 August 2004      ne=monoms+15
c
      ne=14
c
      tint=t
      do 200 i=1,ns
      call evalrf(t,y,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,       &
     & itype)
      do 10 j=1,ne
   10 a(j)=h*f(j)
      do 20 j=1,ne
   20 yt(j)=y(j)+a(j)/9.d+0
      tt=t+h/9.d+0
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 30 j=1,ne
   30 b(j)=h*f(j)
      do 40 j=1,ne
   40 yt(j)=y(j) + (a(j) + 3.d+0*b(j))/24.d+0
      tt=t+h/6.d+0
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 50 j=1,ne
   50 c(j)=h*f(j)
      do 60 j=1,ne
  60  yt(j)=y(j)+(a(j)-3.d+0*b(j)+4.d+0*c(j))/6.d+0
      tt=t+h/3.d+0
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 70 j=1,ne
  70  d(j)=h*f(j)
      do 80 j=1,ne
   80 yt(j)=y(j) + (-5.d+0*a(j) + 27.d+0*b(j) -
     &  24.d+0*c(j) + 6.d+0*d(j))/8.d+0
      tt=t+.5d+0*h
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 90 j=1,ne
   90 e(j)=h*f(j)
      do 100 j=1,ne
  100 yt(j)=y(j) + (221.d+0*a(j) - 981.d+0*b(j) +
     &  867.d+0*c(j)- 102.d+0*d(j) + e(j))/9.d+0
      tt = t+2.d+0*h/3.d+0
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 110 j=1,ne
  110 g(j)=h*f(j)
      do 120 j=1,ne
  120 yt(j) = y(j)+(-183.d+0*a(j)+678.d+0*b(j)-472.d+0*c(j)-
     &  66.d+0*d(j)+80.d+0*e(j) + 3.d+0*g(j))/48.d+0
      tt = t + 5.d+0*h/6.d+0
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 130 j=1,ne
  130 o(j)=h*f(j)
      do 140 j=1,ne
  140 yt(j) = y(j)+(716.d+0*a(j)-2079.d+0*b(j)+1002.d+0*c(j)+
     & 834.d+0*d(j)-454.d+0*e(j)-9.d+0*g(j)+72.d+0*o(j))/82.d+0
      tt = t + h
      call evalrf(tt,yt,f,qmcc,xk,xl,ww,zedge,blg,escale,g0,theta0,     &
     & itype)
      do 150 j=1,ne
  150 p(j)=h*f(j)
      do 160 j=1,ne
  160 y(j) = y(j)+(41.d+0*a(j)+216.d+0*c(j)+27.d+0*d(j)+
     &  272.d+0*e(j)+27.d+0*g(j)+216.d+0*o(j)+41.d+0*p(j))/840.d+0
      t=tint+i*h
  200 continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine multirfsetup(imap5,jmap6)
      use multitrack
      include 'impli.inc'
      include 'map.inc'   ! this contains the reftraj array (need 5 and 6)
      integer imap5,jmap6
      real*8 hx,hy,hxi,hyi
      integer nx,ny,i,j,n,ierror
      integer imin,imax,jmin,jmax
      real*8, dimension(4) :: diag,rdiag
!
! if arrays exist, make sure they have the correct size
      if( allocated(tlistin) .and.                                         &
     &    ((imap5.ne.imaps) .or. (jmap6.ne.jmaps)) )then
        if(idproc.eq.0)write(6,*)'(rfsetup) deallocating multiarrays...'
        call del_multiarrays
        imaps=imap5
        jmaps=jmap6
        if(idproc.eq.0)write(6,*)'(rfsetup) ...allocating multiarrays'
        call new_multiarrays
      endif
! if arrays do not exist, create them:
      if( .not.allocated(tlistin) )then
        imaps=imap5
        jmaps=jmap6
        if(idproc.eq.0)write(6,*)'(rfsetup) allocating multiarrays'
        call new_multiarrays
      endif
!
! save the reference t,pt (these are changed by rfgap, so need to store them)
      refentry5=reftraj(5)
      refentry6=reftraj(6)
!
! Set up grid of initial t-pt values:
! Note: In theory the 5th variable should be periodic. But in most cases it
! will not fill the full 2pi. And it is not worth the hassle to deal with the
! special case that it does fill the full 2pi. So don't bother with periodic:
!
! pack the values into an array of length 4 to minimize MPI calls:
      diag(1)= maxval(zblock(5,1:nraysp))          !tmax
      diag(2)= maxval(zblock(6,1:nraysp))          !ptmax
      diag(3)=-minval(zblock(5,1:nraysp))          !tmin
      diag(4)=-minval(zblock(6,1:nraysp))          !ptmin
      call MPI_ALLREDUCE(diag,rdiag,4,mreal,mpimax,lworld,ierror)
      xmax= rdiag(1)
      ymax =rdiag(2)
      xmin=-rdiag(3)
      ymin=-rdiag(4)
!
! nx, ny are just local variables (could use imaps and jmaps instead)
      nx=imaps
      ny=jmaps
! Note that, in the following, tlist and ptlist are *not* deviation variables,
! as can be seen from the addition of refentry5 and refentry6 in the formulas.
! The reason is that the non-deviation variables must be passed to rfgap.
!
! Compute the tlistin grid points and array inear:
      if( (xmin.eq.xmax).and.(nx.gt.1) )then
        if(idproc.eq.0)                                                   &
     &  write(6,*)'warning: phase spread is zero and imaps.ne.1'
        tlistin(1:imaps)=xmin+refentry5
        inear(1:nraysp)=1
      elseif( (xmin.ne.xmax).and.(nx.gt.1) )then
        hx=(xmax-xmin)/(nx-1)
        hxi=1./hx
        do i=1,nx
          tlistin(i)=xmin+(i-1)*hx+refentry5
        enddo
        do n=1,nraysp
          inear(n)=nint((zblock(5,n)-xmin)*hxi)+1 !nearest grid point
        enddo
      else
!       tlistin(1)=0.5d0*(xmin+xmax)+refentry5
        tlistin(1)=refentry5
        inear(1:nraysp)=1
      endif
! Compute the ptlistin grid points and array jnear:
      if( (ymin.eq.ymax).and.(ny.gt.1) )then
        if(idproc.eq.0)                                                  &
     &  write(6,*)'warning: energy spread is zero and jmaps.ne.1'
        ptlistin(1:jmaps)=ymin+refentry6
        jnear(1:nraysp)=1
      elseif( (ymin.ne.ymax).and.(ny.gt.1) )then
        hy=(ymax-ymin)/(ny-1)
        hyi=1./hy
        do j=1,ny
          ptlistin(j)=ymin+(j-1)*hy+refentry6
        enddo
        do n=1,nraysp
          jnear(n)=nint((zblock(6,n)-ymin)*hyi)+1 !nearest grid point
        enddo
      else
!       ptlistin(1)=0.5d0*(ymin+ymax)+refentry6
        ptlistin(1)=refentry6
        jnear(1:nraysp)=1
      endif
!
! compute tdelt, ptdelt arrays:
      do n=1,nraysp
        tdelt(n)= zblock(5,n)+refentry5- tlistin(inear(n))
        ptdelt(n)=zblock(6,n)+refentry6-ptlistin(jnear(n))
      enddo
!
! In the future, use rho56 to determine whether or not to compute a map:
      rho56=0.d0
      do n=1,nraysp
        rho56(inear(n),jnear(n))=rho56(inear(n),jnear(n))+1.d0
      enddo
!-----
! just for checking/debugging
!     imin=minval(inear)
!     imax=maxval(inear)
!     jmin=minval(jnear)
!     jmax=maxval(jnear)
!     if(idproc.eq.0)write(6,*)'imin,imax=',imin,imax
!     if(idproc.eq.0)write(6,*)'jmin,jmax=',jmin,jmax
!-----
      return
      end
