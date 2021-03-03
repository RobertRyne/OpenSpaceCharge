c******************* ML/I DIAGNOSTIC ROUTINES*******************
      subroutine writereftraj(sarclen,nfile,nprecision,nunits)
      use parallel, only: idproc
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'map.inc'
      include 'usrdat.inc'
      integer :: nfile,nprecision,nunits
      real*8 sarclen
      character*16 :: myformat='(9(1x,1pe12.5 ))'
      character*2 :: a2
      common/envdata/env(6),envold(6),emap(6,6) 
!---format:
      nlength=nprecision+7
      call num2string(nlength,a2,2)
      myformat(10:11)=a2
      call num2string(nprecision,a2,2)
      myformat(13:14)=a2
!
      xl=sl
      wld=omegascl*xl*p0sc
      z=sarclen
! reference trajecory (t-pt portion):
      energy=(gamma-1.)*pmass
      if(idproc.eq.0)then
        if(nunits.eq.0)then
          write(nfile,myformat)z,reftraj(1:6),energy
          write(17,myformat)z,env(1:6)
          call myflush(17)
        elseif(nunits.eq.1)then
          write(nfile,myformat)z,reftraj(1)*xl,reftraj(2)*p0sc,             &
     &    reftraj(3)*xl,reftraj(4)*p0sc,                                    &
     &    reftraj(5)/omegascl,reftraj(6)*wld,energy
          write(17,myformat)z,env(1:6)
          call myflush(17)
        endif
      endif
      ucalc(100)=energy
      return
      end
c
      subroutine writemom2d(sarclen,nfile1,nfile2,nfile3,                 &
     &                      nprecision,nunits,ncorr,includepi,ncent)
cryne 8/11/2001 routine to print some 2nd moments
c arguments:
c   sarclen = current arclength
c   nfile1 = unit # for output xfile
c   nfile2 = unit # for output yfile
c   nfile3 = unit # for output tfile
c   nprecision = precision at which to write data
c   nunits = 0/1 for dimensionless/physical variables in output
c   ncorr = 0/1 to report cross-terms (<xpx>,...) as is/as ratios
c   includepi = 0/1 to no/yes divide emittances by pi
c   ncent = 0/1 to remove/keep beam centroid in emittance computation
      use beamdata
      use rays
      use ml_timer
      include 'impli.inc'
      include 'usrdat.inc'
!
      integer :: nfile1,nfile2,nfile3,nprecision,nunits,ncorr,ncent
      real*8 sarclen
!
      logical msk
      dimension msk(nraysp)
      dimension avg(6),ravg(6)
      dimension diag(9),rdiag(9)
      character*16 :: myformat='(6(1x,1pe12.5 ))'
      character*2 :: a2
      real*8 :: zero=0.
!
      call increment_timer('moments',0)
!
c set up particle mask
      msk(1:nraysp)=.true.
      do j=1,nraysp
        if(istat(j).ne.0)msk(j)=.false.
      enddo
!
c compute moments
c   scale factors
      den1=1./nrays
      den2=den1*den1
      slp0=sl*p0sc
      slp1=slp0
      wld=omegascl*slp0
      clyte=299792458.d0
      pi=2.d0*asin(1.d0)
      z=sarclen
c   compute beam centroid
      ravg(1:6)=0.d0
      if(ncent.eq.0)then
        avg(1)=sum(zblock(1,1:nraysp),msk)*den1
        avg(2)=sum(zblock(2,1:nraysp),msk)*den1
        avg(3)=sum(zblock(3,1:nraysp),msk)*den1
        avg(4)=sum(zblock(4,1:nraysp),msk)*den1
        avg(5)=sum(zblock(5,1:nraysp),msk)*den1
        avg(6)=sum(zblock(6,1:nraysp),msk)*den1
        call MPI_ALLREDUCE(avg,ravg,6,mreal,mpisum,lworld,ierror)
      endif
c   give readable names to coordinate entries
      xbar=ravg(1)
      ybar=ravg(3)
      tbar=ravg(5)
c   compute second-order moments
      diag(1)=sum((zblock(1,1:nraysp)-ravg(1))                          &
     &           *(zblock(1,1:nraysp)-ravg(1)),msk)
      diag(2)=sum((zblock(2,1:nraysp)-ravg(2))                          &
     &           *(zblock(2,1:nraysp)-ravg(2)),msk)
      diag(3)=sum((zblock(3,1:nraysp)-ravg(3))                          &
     &           *(zblock(3,1:nraysp)-ravg(3)),msk)
      diag(4)=sum((zblock(4,1:nraysp)-ravg(4))                          &
     &           *(zblock(4,1:nraysp)-ravg(4)),msk)
      diag(5)=sum((zblock(5,1:nraysp)-ravg(1))                          &
     &           *(zblock(5,1:nraysp)-ravg(2)),msk)
      diag(6)=sum((zblock(6,1:nraysp)-ravg(3))                          &
     &           *(zblock(6,1:nraysp)-ravg(4)),msk)
      diag(7)=sum((zblock(1,1:nraysp)-ravg(1))                          &
     &           *(zblock(2,1:nraysp)-ravg(2)),msk)
      diag(8)=sum((zblock(3,1:nraysp)-ravg(3))                          &
     &           *(zblock(4,1:nraysp)-ravg(4)),msk)
      diag(9)=sum((zblock(5,1:nraysp)-ravg(5))                          &
     &           *(zblock(6,1:nraysp)-ravg(6)),msk)
      call MPI_ALLREDUCE(diag,rdiag,9,mreal,mpisum,lworld,ierror)
c   give readable names to entries, and divide by number of particles
      sq1=rdiag(1)*den1
      sq2=rdiag(2)*den1
      sq3=rdiag(3)*den1
      sq4=rdiag(4)*den1
      sq5=rdiag(5)*den1
      sq6=rdiag(6)*den1
      xpx=rdiag(7)*den1
      ypy=rdiag(8)*den1
      tpt=rdiag(9)*den1
c   compute derived quantities
      epsx2=(sq1*sq2-xpx*xpx)
      epsy2=(sq3*sq4-ypy*ypy)
      epst2=(sq5*sq6-tpt*tpt)
      xrms=sqrt(sq1)
      yrms=sqrt(sq3)
      trms=sqrt(sq5)
      pxrms=sqrt(sq2)
      pyrms=sqrt(sq4)
      ptrms=sqrt(sq6)
      epsx=sqrt(max(epsx2,zero))
      epsy=sqrt(max(epsy2,zero))
      epst=sqrt(max(epst2,zero))
      if(includepi.eq.1)then
        epsx=epsx/pi
        epsy=epsy/pi
        epst=epst/pi
      endif
      if(ncorr.eq.1)then
        xpxfac=0.d0
        ypyfac=0.d0
        tptfac=0.d0
        if(xrms.ne.0. .and. pxrms.ne.0.)xpxfac=1./(xrms*pxrms)
        if(yrms.ne.0. .and. pyrms.ne.0.)ypyfac=1./(yrms*pyrms)
        if(trms.ne.0. .and. ptrms.ne.0.)tptfac=1./(trms*ptrms)
        xpx=xpx*xpxfac
        ypy=ypy*ypyfac
        tpt=tpt*tptfac
        slp1=1.d0
      endif
!
c output results
      if(idproc.eq.0)then
c   set format for amount of data and desired precision
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
c   write 2D (transverse 2nd moments)
        if(nunits.eq.0)then
          write(nfile1,myformat)z,xrms,pxrms,xpx,epsx,xbar
          write(nfile2,myformat)z,yrms,pyrms,ypy,epsy,ybar
        else
          write(nfile1,myformat)z,xrms*sl,pxrms*p0sc,xpx*slp1,epsx*slp0,&
     &                            xbar*sl
          write(nfile2,myformat)z,yrms*sl,pyrms*p0sc,ypy*slp1,epsy*slp0,&
     &                            ybar*sl
        endif
      endif
      ucalc(1)=xrms
      ucalc(2)=pxrms
      ucalc(3)=xpx
      ucalc(4)=yrms
      ucalc(5)=pyrms
      ucalc(6)=ypy
c   write 3D (t-pt 2nd moments)
      if(idproc.eq.0)then
        if(nunits.eq.0)then
          write(nfile3,myformat)z,trms,ptrms,tpt,epst,tbar
        else
          write(nfile3,myformat)z,trms*beta*clyte/omegascl,ptrms*wld,   &
     &                          tpt*slp1,epst*slp0,tbar/omegascl
        endif
      endif
      ucalc(7)=trms
      ucalc(8)=ptrms
      ucalc(9)=tpt
!
      call increment_timer('moments',1)
      return
      end
!
      subroutine writeemit(sarclen,nfile,nprecision,nunits,
     &                     ne2,ne4,ne6,ncent)
cabell 8/29/2004 routine to write 2- 4- and/or 6-d emittances
c arguments:
c   sarclen = current arclength
c   nfile = unit # for output file
c   nprecision = precision at which to write data
c   nunits = 0/1 for dimensionless/physical variables in output
c   ne2 = flag to compute 2-d emittances
c   ne4 = flag to compute 4-d transverse emittance
c   ne6 = flag to compute 6-d emittance
c   ncent = 0/1 to remove/keep beam centroid in emittance computation
      use beamdata
      use rays
      use ml_timer
      include 'impli.inc'
!
      integer :: nfile,nprecision,nunits,ne2,ne4,ne6,ncent
      real*8 sarclen
!
      logical msk
      dimension msk(nraysp)
      dimension avg(6),ravg(6)
      dimension diag(21),rdiag(21)
      dimension sigma6(6,6)
      character*16 :: myformat='(6(1x,1pe12.5 ))'
      character*1 :: a1
      character*2 :: a2
      real*8 :: zero=0.
!
      call increment_timer('moments',0)
!
!      if(idproc.eq.0) then
!        write(6,*) "writeemit()::"
!        if(nunits.eq.0) write(6,*) "  using scaled variables"
!        if(nunits.eq.1) write(6,*) "  using physical variables"
!        if(ne2.eq.1) write(6,*)    "    --2-D"
!        if(ne4.eq.1) write(6,*)    "    --4-D"
!        if(ne6.eq.1) write(6,*)    "    --6-D"
!      endif
!
c set up particle mask
      msk(1:nraysp)=.true.
      do j=1,nraysp
        if(istat(j).ne.0)msk(j)=.false.
      enddo
!
c compute emittances
c   scale factors
      den1=1./nrays
      slp0=sl*p0sc
      clyte=299792458.d0
      wld=omegascl*slp0
      z=sarclen
c   compute beam centroid
      ravg(1:6)=0.d0
      if(ncent.eq.0)then
        avg(1)=sum(zblock(1,1:nraysp),msk)*den1
        avg(2)=sum(zblock(2,1:nraysp),msk)*den1
        avg(3)=sum(zblock(3,1:nraysp),msk)*den1
        avg(4)=sum(zblock(4,1:nraysp),msk)*den1
        avg(5)=sum(zblock(5,1:nraysp),msk)*den1
        avg(6)=sum(zblock(6,1:nraysp),msk)*den1
        call MPI_ALLREDUCE(avg,ravg,6,mreal,mpisum,lworld,ierror)
      endif
c   compute entries of beam sigma matrix
c   entries of diag() follow diagonals of the beam matrix
c   starting from the main diagonal and going out to UR corner
c   compute only required entries
      diag(1:21)=0
      diag( 1)=sum((zblock(1,1:nraysp)-ravg(1))                         &
     &            *(zblock(1,1:nraysp)-ravg(1)),msk)
      diag( 2)=sum((zblock(2,1:nraysp)-ravg(2))                         &
     &            *(zblock(2,1:nraysp)-ravg(2)),msk)
      diag( 3)=sum((zblock(3,1:nraysp)-ravg(3))                         &
     &            *(zblock(3,1:nraysp)-ravg(3)),msk)
      diag( 4)=sum((zblock(4,1:nraysp)-ravg(4))                         &
     &            *(zblock(4,1:nraysp)-ravg(4)),msk)
      diag( 7)=sum((zblock(1,1:nraysp)-ravg(1))                         &
     &            *(zblock(2,1:nraysp)-ravg(2)),msk)
      diag( 9)=sum((zblock(3,1:nraysp)-ravg(3))                         &
     &            *(zblock(4,1:nraysp)-ravg(4)),msk)
      if(ne2.eq.1.or.ne6.eq.1)then
        diag( 5)=sum((zblock(5,1:nraysp)-ravg(5))                       &
     &              *(zblock(5,1:nraysp)-ravg(5)),msk)
        diag( 6)=sum((zblock(6,1:nraysp)-ravg(6))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
        diag(11)=sum((zblock(5,1:nraysp)-ravg(5))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
      endif
      if(ne4.eq.1.or.ne6.eq.1)then
        diag( 8)=sum((zblock(2,1:nraysp)-ravg(2))                       &
     &              *(zblock(3,1:nraysp)-ravg(3)),msk)
        diag(12)=sum((zblock(1,1:nraysp)-ravg(1))                       &
     &              *(zblock(3,1:nraysp)-ravg(3)),msk)
        diag(13)=sum((zblock(2,1:nraysp)-ravg(2))                       &
     &              *(zblock(4,1:nraysp)-ravg(4)),msk)
        diag(16)=sum((zblock(1,1:nraysp)-ravg(1))                       &
     &              *(zblock(4,1:nraysp)-ravg(4)),msk)
      endif
      if(ne6.eq.1)then
        diag(10)=sum((zblock(4,1:nraysp)-ravg(4))                       &
     &              *(zblock(5,1:nraysp)-ravg(5)),msk)
        diag(14)=sum((zblock(3,1:nraysp)-ravg(3))                       &
     &              *(zblock(5,1:nraysp)-ravg(5)),msk)
        diag(15)=sum((zblock(4,1:nraysp)-ravg(4))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
        diag(17)=sum((zblock(2,1:nraysp)-ravg(2))                       &
     &              *(zblock(5,1:nraysp)-ravg(5)),msk)
        diag(18)=sum((zblock(3,1:nraysp)-ravg(3))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
        diag(19)=sum((zblock(1,1:nraysp)-ravg(1))                       &
     &              *(zblock(5,1:nraysp)-ravg(5)),msk)
        diag(20)=sum((zblock(2,1:nraysp)-ravg(2))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
        diag(21)=sum((zblock(1,1:nraysp)-ravg(1))                       &
     &              *(zblock(6,1:nraysp)-ravg(6)),msk)
      endif
      call MPI_ALLREDUCE(diag,rdiag,21,mreal,mpisum,lworld,ierror)
c  give readable names to entries, and divide by number of particles
        xx=rdiag( 1)*den1
      pxpx=rdiag( 2)*den1
        yy=rdiag( 3)*den1
      pypy=rdiag( 4)*den1
        tt=rdiag( 5)*den1
      ptpt=rdiag( 6)*den1
       xpx=rdiag( 7)*den1
       pxy=rdiag( 8)*den1
       ypy=rdiag( 9)*den1
       pyt=rdiag(10)*den1
       tpt=rdiag(11)*den1
        xy=rdiag(12)*den1
      pxpy=rdiag(13)*den1
        yt=rdiag(14)*den1
      pypt=rdiag(15)*den1
       xpy=rdiag(16)*den1
       pxt=rdiag(17)*den1
       ypt=rdiag(18)*den1
        xt=rdiag(19)*den1
      pxpt=rdiag(20)*den1
       xpt=rdiag(21)*den1
c epsX2=det(2x2 beam sigma matrix for horizontal plane)
c epsY2=det(2x2 beam sigma matrix for vertical plane)
c epsT2=det(2x2 beam sigma matrix for longitudinal plane)
      if(ne2.eq.1)then
        epsX2=(xx*pxpx-xpx*xpx)
        epsY2=(yy*pypy-ypy*ypy)
        epsT2=(tt*ptpt-tpt*tpt)
        epsX=sqrt(max(epsX2,zero))
        epsY=sqrt(max(epsY2,zero))
        epsT=sqrt(max(epsT2,zero))
        if(nunits.eq.1)then
          epsX=epsX*slp0
          epsY=epsY*slp0
          epsT=epsT*slp0
        endif
      endif
c epsXY2=det(4x4 beam sigma matrix for transverse planes)
      if(ne4.eq.1)then
        epsXY2=xx*yy*pxpx*pypy                                          &
     &   + xpx**2*ypy**2 + xy**2*pxpy**2 + xpy**2*pxy**2                &
     &   - xx*(pxpx*ypy**2 + yy*pxpy**2 + pypy*pxy**2)                  &
     &   - pxpx*(yy*xpy**2 + pypy*xy**2) - yy*pypy*xpx**2               &
     &   + 2*(xx*pxy*pxpy + pxpx*xy*xpy)*ypy                            &
     &   + 2*(yy*xpy*pxpy + pypy*xy*pxy)*xpx                            &
     &   - 2*xpx*(xpy*pxy + xy*pxpy)*ypy - 2*xy*xpy*pxy*pxpy
        epsXYth=max(epsXY2,zero)**0.25
        if(nunits.eq.1)then
          epsXYth=epsXYth*slp0
        endif
      endif
c epsXYT2=det(6x6 beam sigma matrix)
c  this COULD be coded by hand, but the det contains 388 terms!
      if(ne6.eq.1)then
        sigma6(1,1)=xx
        sigma6(2,1)=xpx
        sigma6(3,1)=xy
        sigma6(4,1)=xpy
        sigma6(5,1)=xt
        sigma6(6,1)=xpt
        sigma6(2,2)=pxpx
        sigma6(3,2)=pxy
        sigma6(4,2)=pxpy
        sigma6(5,2)=pxt
        sigma6(6,2)=pxpt
        sigma6(3,3)=yy
        sigma6(4,3)=ypy
        sigma6(5,3)=yt
        sigma6(6,3)=ypt
        sigma6(4,4)=pypy
        sigma6(5,4)=pyt
        sigma6(6,4)=pypt
        sigma6(5,5)=tt
        sigma6(6,5)=tpt
        sigma6(6,6)=ptpt
        do i=1,6
          do j=i+1,6
            sigma6(i,j)=sigma6(j,i)
          enddo
        enddo
cabell:need to compute epsXYT2=det(sigma6)
        epsXYT2=0.0
        epsXYTth=(max(epsXYT2,zero))**(1.d0/6.d0)
        if(nunits.eq.1)then
          epsXYTth=epsXYTth*slp0
        endif
      endif
!
c output results
      if(idproc.eq.0)then
c   set format for amount of data and desired precision
        nout=1
        if(ne2.eq.1)nout=nout+3
        if(ne4.eq.1)nout=nout+1
        if(ne6.eq.1)nout=nout+1
        call num2string(nout,a1,1)
        myformat(2:2)=a1
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
c   write results
        if(ne2.eq.1.and.ne4.eq.0.and.ne6.eq.0)then
c   2d only
          write(nfile,myformat)z,epsX,epsY,epsT
        else if(ne2.eq.0.and.ne4.eq.1.and.ne6.eq.0)then
c   4d only
          write(nfile,myformat)z,epsXYth
        else if(ne2.eq.0.and.ne4.eq.0.and.ne6.eq.1)then
c   6d only
          write(nfile,myformat)z,epsXYTth
        else if(ne2.eq.1.and.ne4.eq.1.and.ne6.eq.0)then
c   2d and 4d
          write(nfile,myformat)z,epsX,epsY,epsT,epsXYth
        else if(ne2.eq.1.and.ne4.eq.0.and.ne6.eq.1)then
c   2d and 6d
          write(nfile,myformat)z,epsX,epsY,epsT,epsXYTth
        else if(ne2.eq.0.and.ne4.eq.1.and.ne6.eq.1)then
c   4d and 6d
          write(nfile,myformat)z,epsXYth,epsXYTth
        else
c   2d and 4d and 6d
          write(nfile,myformat)z,epsX,epsY,epsT,epsXYth,epsXYTth
        endif
      endif
!
      call increment_timer('moments',1)
      return
      end
!
!
      subroutine writemaxsize(sarclen,nfile1,nfile2,nfile3,nfile4,        &
     &                      nprecision)
cryne 8/11/2001 routine to print min/max beam size
cryne various mods on 11/26/03
      use beamdata
      use rays
      include 'impli.inc'
      integer :: nfile1,nfile2,nfile3,nfile4,nprecision
      real*8 sarclen
      logical msk
      dimension msk(nraysp)
      dimension mloca(1),diag(16),rdiag(16)
!
      character*16 :: myformat='(6(1x,1pe12.5 ))'
      character*2 :: a2
      if(idproc.eq.0)then
!---format:
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
      endif

!
!     if(idproc.eq.0)write(6,*)'inside routine printmaxsize'
      msk(1:nraysp)=.true.
      do j=1,nraysp
        if(istat(j).ne.0)msk(j)=.false.
      enddo
!
!------------ min/max quantities1:nrays----------------
!     pxmax=maxval(abs(zblock(2,1:nrays)))/gambet*econ
!     pymax=maxval(abs(zblock(4,1:nrays)))/gambet*econ
! minvals and maxvals:
      diag(1)=maxval(zblock(1,1:nraysp))
      diag(2)=maxval(zblock(2,1:nraysp))
      diag(3)=maxval(zblock(3,1:nraysp))
      diag(4)=maxval(zblock(4,1:nraysp))
      diag(5)=maxval(zblock(5,1:nraysp))
      diag(6)=maxval(zblock(6,1:nraysp))
      diag(7)=-minval(zblock(1,1:nraysp))
      diag(8)=-minval(zblock(2,1:nraysp))
      diag(9)=-minval(zblock(3,1:nraysp))
      diag(10)=-minval(zblock(4,1:nraysp))
      diag(11)=-minval(zblock(5,1:nraysp))
      diag(12)=-minval(zblock(6,1:nraysp))
      call MPI_ALLREDUCE(diag,rdiag,12,mreal,mpimax,lworld,ierror)
      clyte=299792458.d0
      xl=sl
      xk=1./xl
      z=sarclen
      econ=1.0
      wld=omegascl*xl*p0sc
      xmax=rdiag(1)*xl
      pxmax=rdiag(2)*p0sc
      ymax=rdiag(3)*xl
      pymax=rdiag(4)*p0sc
      tmax=rdiag(5)/omegascl
      ptmax=rdiag(6)*wld
!     zmax=rdiag(5)*beta/xk
      zmax=rdiag(5)*beta*clyte/omegascl
      emax=rdiag(6)/(gamma**3*beta**2)*econ
      xmin=-rdiag(7)*xl
      pxmin=-rdiag(8)*p0sc
      ymin=-rdiag(9)*xl
      pymin=-rdiag(10)*p0sc
      tmin=-rdiag(11)/omegascl
      ptmin=-rdiag(12)*wld
!     zmin=-rdiag(11)*beta/xk
      zmin=-rdiag(11)*beta*clyte/omegascl
      emin=-rdiag(12)/(gamma**3*beta**2)*econ
! location of min/max
!     diag(15)=diag(8)
!     diag(13)=diag(7)
!     diag(11)=diag(6)
!     diag(9)=diag(5)
!     diag(7)=diag(4)
!     diag(5)=diag(3)
!     diag(3)=diag(2)
! x min/max1:nrays
!     mloca=maxloc(zblock(1,1:nraysp))
!     diag(2)=mloca(1)+maxrayp*idproc
!     mloca=maxloc(zblock(3,1:nraysp))
!     diag(4)=mloca(1)+maxrayp*idproc
!     mloca=maxloc(zblock(5,1:nraysp))
!     diag(6)=mloca(1)+maxrayp*idproc
!     mloca=maxloc(zblock(6,1:nraysp))
!     diag(8)=mloca(1)+maxrayp*idproc
!     mloca=minloc(zblock(1,1:nraysp))
!     diag(10)=mloca(1)+maxrayp*idproc
!     mloca=minloc(zblock(3,1:nraysp))
!     diag(12)=mloca(1)+maxrayp*idproc
!     mloca=minloc(zblock(5,1:nraysp))
!     diag(14)=mloca(1)+maxrayp*idproc
!     mloca=minloc(zblock(6,1:nraysp))
!     diag(16)=mloca(1)+maxrayp*idproc
!     call MPI_ALLREDUCE(diag,rdiag,8,m2real,mpimaxloc,lworld,ierror)
!     maxpcl1=rdiag(2)
!     maxpcl3=rdiag(4)
!     maxpcl5=rdiag(6)
!     maxpcl6=rdiag(8)
!     minpcl1=rdiag(10)
!     minpcl3=rdiag(12)
!     minpcl5=rdiag(14)
!     minpcl6=rdiag(16)
      if(idproc.eq.0)then
        write(nfile1,myformat)z,xmin,xmax,pxmin,pxmax
        write(nfile2,myformat)z,ymin,ymax,pymin,pymax
        write(nfile3,myformat)z,tmin,tmax,ptmin,ptmax
        write(nfile4,myformat)z,zmin,zmax,emin,emax
!!!     write(nfile4,myformat)z,gamma*zmin,gamma*zmax,minpcl5,maxpcl5
!       call myflush(nfile1)
!       call myflush(nfile2)
!       call myflush(nfile3)
!       call myflush(nfile4)
      endif
      return
      end
!
      subroutine profile1d(ncol,nbins,rwall,sarclen,nfile,fname,          &
     &                     nprecision)
      use beamdata
      use rays
      use lieaparam, only : monoms
      use ml_timer
      include 'impli.inc'
      include 'map.inc'
      real*8, dimension(nbins) :: rho1d,rhotmp
      real*8, dimension(2) :: bndy,rbndy
      character*16 fname
      call increment_timer('profile1d',0)
!     write(6,*)'PE#',idproc,' is in profile1d w/ nraysp=',nraysp
!
      if(idproc.eq.0)then
!     write(6,*)'s= ',sarclen,';writing 1D profile on file ',fname
      endif
      if(rwall.gt.0.d0)then
        xmin=-rwall
        xmax=+rwall
      else
!     determine the range from the particle coords:
        bndy(1)= maxval(zblock(ncol,:))
        bndy(2)=-minval(zblock(ncol,:))
        call MPI_ALLREDUCE(bndy,rbndy,2,mreal,mpimax,lworld,ierr)
        xmax=rbndy(1)
        xmin=-rbndy(2)
!DO  express in physical units (i.e. meters):
        xmin=xmin*sl
        xmax=xmax*sl
!       write(6,1122)sl
!1122   format(1x,'sl=',d21.14)
! a little bigger to make sure nothing falls off the end:
        xmin=xmin*(1.d0+1.d-8)
        xmax=xmax*(1.d0+1.d-8)
      endif
!     write(6,*)'inside profile1d w/ min,max=',xmin,xmax
! note that nbins is REALLY the number of grid points, not bins:
      hx=(xmax-xmin)/(nbins-1)
      hxi=1.d0/hx
      rho1d=0.d0
      do n=1,nraysp
        indx=(zblock(ncol,n)*sl-xmin)*hxi + 1
        if(indx.lt.1 .or. indx.gt.nbins-1)then
!         write(6,*)'particle #',n,' is out of range; indx=',indx
!         write(6,*)'zblock(' ,ncol, ',' ,n, ')=' ,zblock(ncol,n)
          cycle
        endif
        ab=((xmin-zblock(ncol,n)*sl)+indx*hx)*hxi
        rho1d(indx)=rho1d(indx)+ab
        rho1d(indx+1)=rho1d(indx+1)+(1.d0-ab)
      enddo
!
!-------------------------------------------
! now that every processor has its own 1d profile, sum them together
! on processor 0:
      if(idproc.eq.0)then
        do l=1,nvp-1
          call MPI_RECV(rhotmp,nbins,mreal,l,95,lworld,mpistat,ierr)
!         write(6,*)'PE#0 has received rhotmp data from proc# ',l
          rho1d(:)=rho1d(:)+rhotmp(:)
        enddo
      else
        call MPI_SEND(rho1d,nbins,mreal,0,95,lworld,ierr)
!       write(6,*)'PE#',idproc,' has sent is rho1d data.'
      endif
!-------------------------------------------
!
      if(idproc.eq.0)then
        do n=1,nbins
          xval=xmin+(n-1)*hx
!         write(nfile,999,err=1001)xval,rho1d(n),arclen
          write(nfile,999)xval,rho1d(n)
        enddo
  999   format(3(1pe14.7,1x))
        total=sum(rho1d)
        write(6,*)'s= ',sarclen,' ;total deposited=',total,              &
     &            '; 1D profile to file ',fname
        call increment_timer('profile1d',1)
        return
 1001   write(6,*)'PE 0: TROUBLE WITH WRITE STMT IN PROFILE1D'
        return
      endif
      call increment_timer('profile1d',1)
!     write(6,*)'PE# ',idproc,' is returning from profile1d'
      return
      end
!
      subroutine profilerad(ncol,nbins,rwall,sarclen,nfile,fname,          &
     &                     nprecision)
      use beamdata
      use rays
      use lieaparam, only : monoms
      use ml_timer
      include 'impli.inc'
      include 'map.inc'
      real*8, dimension(nbins) :: rho1d,rhotmp,rhotest
      real*8, dimension(2) :: bndy,rbndy
      real*8, dimension(nraysp) :: radii
      character*16 fname
      call increment_timer('profile1d',0)
!     write(6,*)'PE#',idproc,' is in profile1d w/ nraysp=',nraysp
!
      if(idproc.eq.0)then
      write(6,*)'s= ',sarclen,';writing radial profile on file ',fname
      endif
      clyte=299792458.d0
      t2z=beta*clyte/omegascl
      do n=1,nraysp
        radii(n)=sqrt(                                                    &
     &  (sl*zblock(1,n))**2+(sl*zblock(3,n))**2+(t2z*zblock(5,n))**2 )
      enddo
      if(rwall.gt.0.d0)then
        xmin=0.d0
        xmax=+rwall
      else
!     determine the range from the particle coords:
        bndy(1)=maxval(radii(1:nraysp))
        bndy(2)=0.d0
        if(idproc.eq.0)then
        write(6,*)'ready to call MPI_ALLREDUCE'
        endif
        call MPI_ALLREDUCE(bndy,rbndy,2,mreal,mpimax,lworld,ierr)
        if(idproc.eq.0)then
        write(6,*)'returned from MPI_ALLREDUCE'
        endif
        xmax=rbndy(1)
        xmin=0.d0
!DO  express in physical units (i.e. meters):
!       xmin=xmin*sl
!       xmax=xmax*sl
!       write(6,1122)sl
!1122   format(1x,'sl=',d21.14)
! a little bigger to make sure nothing falls off the end:
!!!!!!!!xmin=xmin*(1.d0+1.d-8)
        xmax=xmax*(1.d0+1.d-8)
      endif
!     write(6,*)'inside profilerad w/ min,max=',xmin,xmax
! (old comment) note that nbins is REALLY the number of grid points, not bins.
! ABOVE COMMENT IS WRONG: in this routine nbins really IS the # of bins
!!!!!!hx=(xmax-xmin)/(nbins-1)
      hx=xmax/nbins
      hxi=1.d0/hx
      bcon=hx*8.*asin(1.d0)*nbins
      rho1d=0.d0
!cryne This routine is a hacked version of the 1d histogram routine.
!cryne As such, its treatment of the value at r=0 is fishy, since I
!cryne will end up dividing by the radial values.
      do n=1,nraysp
        indx=(radii(n)-xmin)*hxi + 1
        if(indx.lt.1 .or. indx.gt.nbins)then
!       if(indx.lt.1 .or. indx.gt.nbins-1)then  !applies for linear weighting???
!         write(6,*)'particle #',n,' is out of range; indx=',indx
!         write(6,*)'zblock(' ,ncol, ',' ,n, ')=' ,zblock(ncol,n)
          if(rwall.gt.0.d0)then
           cycle
          else
            write(6,*)'TROUBLE IN PROFILERAD, PE=',idproc
            write(6,*)'xmin,xmax=',xmin,xmax
            write(6,*)'n,radii(n)=',n,radii(n)
            write(6,*)'indx=',indx
            call myexit
          endif
        endif
!nearest grid point:
!       rho1d(indx)=rho1d(indx)+(radii(n))**2
!!      rho1d(indx)=rho1d(indx)+(hx*(n-1))**2
!!!     rho1d(indx)=rho1d(indx)+(hx*(n-1)+0.5*hx)**2
        rho1d(indx)=rho1d(indx)+1.d0
!linear weighting:
!       ab=((xmin-radii(n))+indx*hx)*hxi
!       rho1d(indx)=rho1d(indx)+ab
!       rho1d(indx+1)=rho1d(indx+1)+(1.d0-ab)
      enddo
!
!-------------------------------------------
! now that every processor has its own 1d profile, sum them together
! on processor 0:
      if(idproc.eq.0)then
        do l=1,nvp-1
          call MPI_RECV(rhotmp,nbins,mreal,l,95,lworld,mpistat,ierr)
!         write(6,*)'PE#0 has received rhotmp data from proc# ',l
          rho1d(:)=rho1d(:)+rhotmp(:)
        enddo
      else
        call MPI_SEND(rho1d,nbins,mreal,0,95,lworld,ierr)
!       write(6,*)'PE#',idproc,' has sent is rho1d data.'
      endif
!-------------------------------------------
      if(idproc.eq.0)then
      write(6,*)'done with SEND and RECV'
      endif
!
!cryne Here is the code for dealing the the radial values;
!cryne For now, I simply use the radial value at the bin center:
!cryne (Also, use rhotmp to hold the unscaled values)
!cryne (Also, use rhotest to scale by the radial value at the bin edge)
      do n=1,nbins
        rhotmp(n)=rho1d(n)
        redge=hx*(n-1)
        rcent=hx*(n-1)+0.5*hx
        rho1d(n)=rho1d(n)/(bcon*rcent**2)
        if(n.eq.1)rhotest(n)=0.d0
        if(n.ne.1)rhotest(n)=rhotmp(n)/(bcon*redge**2)
      enddo
      if(idproc.eq.0)then
      write(6,*)'done computing rval and recomputing rho1d'
      endif
!
      if(idproc.eq.0)then
        do n=1,nbins
          xedge=xmin+(n-1)*hx
          xcent=xmin+(n-1)*hx+0.5*hx
!         write(nfile,999,err=1001)xcent,rho1d(n),arclen
          write(nfile,999)xcent,rho1d(n),xedge,rhotest(n),rhotmp(n)
        enddo
  999   format(5(1pe14.7,1x))
        total=sum(rho1d)
        write(6,*)'Comment from Ryne: the following is not correct'
        write(6,*)'for the radial histogram calculation. Fix later.'
        write(6,*)'s= ',sarclen,' ;total deposited=',total,              &
     &            '; 1D profile to file ',fname
        call increment_timer('profile1d',1)
        return
 1001   write(6,*)'PE 0: TROUBLE WITH WRITE STMT IN PROFILE1D'
        return
      endif
      call increment_timer('profile1d',1)
!     write(6,*)'PE# ',idproc,' is returning from profile1d'
      return
      end
!
      subroutine writeenv2d(sarclen,nfile1,nfile2,nfile3,                 &
     &                      nprecision,nunits,ncorr)
cryne routine to print 2nd moments obtained from env array
      use parallel, only: idproc
      use beamdata
      include 'impli.inc'
      include 'usrdat.inc'
      integer :: nfile1,nfile2,nfile3,nprecision,nunits
      real*8 sarclen
      common/envdata/env(6),envold(6),emap(6,6) 
      common/emitdata/emxn2,emyn2,emtn2
!
      character*16 :: myformat='(6(1x,1pe12.5 ))'
      character*2 :: a2

      if(idproc.eq.0)then
!---format:
        nlength=nprecision+7
        call num2string(nlength,a2,2)
        myformat(10:11)=a2
        call num2string(nprecision,a2,2)
        myformat(13:14)=a2
      endif
!
      xl=sl
      xk=1./xl
      xlp0=xl*p0sc
      xlp1=xlp0
      gambet=gamma*beta
      clyte=299792458.d0
      wld=omegascl*xl*p0sc
      z=sarclen
! 
      xrms=env(1)
      xpx=env(2)*xrms
      pxrms=0.d0
      if(xrms.ne.0.d0)pxrms=sqrt(emxn2+xpx**2)/xrms
      epsx=sqrt(emxn2)
!
      yrms=env(3)
      ypy=env(4)*yrms
      pyrms=0.d0
      if(yrms.ne.0.d0)pyrms=sqrt(emyn2+ypy**2)/yrms
      epsy=sqrt(emyn2)
!
      trms=env(5)
      tpt=env(6)*trms
      ptrms=0.d0
      if(trms.ne.0.d0)ptrms=sqrt(emtn2+tpt**2)/trms
      epst=sqrt(emtn2)
      if(ncorr.eq.1)then
        xpxfac=0.
        ypyfac=0.
        tptfac=0.
        if(xrms.ne.0. .and. pxrms.ne.0.)xpxfac=1./(xrms*pxrms)
        if(yrms.ne.0. .and. pyrms.ne.0.)ypyfac=1./(yrms*pyrms)
        if(trms.ne.0. .and. ptrms.ne.0.)tptfac=1./(trms*ptrms)
        xpx=xpx*xpxfac
        ypy=ypy*ypyfac
        tpt=tpt*tptfac
        xlp1=1.d0
      endif
      if(idproc.eq.0)then
        if(nunits.eq.0)then
          write(nfile1,myformat)z,xrms,pxrms,xpx,epsx
          write(nfile2,myformat)z,yrms,pyrms,ypy,epsy
        else
          write(nfile1,myformat)z,xrms*xl,pxrms*p0sc,xpx*xlp1,epsx*xlp0
          write(nfile2,myformat)z,yrms*xl,pyrms*p0sc,ypy*xlp1,epsy*xlp0
        endif
      endif
      ucalc(31)=xrms
      ucalc(32)=pxrms
      ucalc(33)=xpx
      ucalc(34)=yrms
      ucalc(35)=pyrms
      ucalc(36)=ypy
!
! 3D (t-pt 2nd moments)
      if(idproc.eq.0)then
        if(nunits.eq.0)then
        write(nfile3,myformat)z,trms,ptrms,tpt,epst
        else
        write(nfile3,myformat)z,                                        &
     &  trms*beta*clyte/omegascl,ptrms*wld,tpt*xlp1,epst*xlp0
        endif
      endif
      ucalc(37)=trms
      ucalc(38)=ptrms
      ucalc(39)=tpt
      return
      end
!
