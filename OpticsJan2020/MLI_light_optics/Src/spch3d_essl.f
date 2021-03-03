c IMPACT 3D space charge routines
c Copyright 2001 University of California
c
c at this point np is the # of particles on a processor:
      subroutine spch3d(c,ex,ey,ez,msk,np,ntot,                         &
     &                  nx,ny,nz,n1,n2,n3,n3a,nadj)
      use parallel
      use ml_timer
      implicit double precision(a-h,o-z)
      real*8, dimension(6,np) :: c
      logical, dimension(np) :: msk
      real*8, dimension(np) :: ex,ey,ez
      real*8, dimension(nx,ny,nz) :: rho,exg,eyg,ezg,rhosum
!     real*8,dimension(n1,n2,n3a) :: rho2
!     complex*16,dimension(n3a,n2,n1) :: grnxtr,rho2xtr
!
      real*8,dimension(n1+2,n2,n3a) :: rho2
      complex*16,dimension(n1/2+1,n2,n3a) :: grnxtr,rho2xtr
!
      integer, parameter :: naux=120000
c rho=charge density on the grid
c rho2=...on doubled grid
c rho2xtr=...xformed (by fft) and transposed
c grnxtr=green function, xformed, transposed
c weights, indices associated with area weighting:
      dimension ab(np),de(np),gh(np),indx(np),jndx(np),kndx(np),
     &indxp1(np),jndxp1(np),kndxp1(np)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/accum/at10,at21,at32,at43,at54,at65,at76,at70
      common/showme/iverbose
      real*8 aux(naux)
!     if(idproc.eq.0)write(6,*)'inside spchg; nadj=',nadj
c
!     call system_clock(count_rate=ihertz)
!     hertz=ihertz
!     call system_clock(count=iticks0)
c
! compute the Green function on the grid (use rho2 for storage):
      call increment_timer('greenf3d',0)
      call greenf3d(rho2,nx,ny,nz,n1,n2,n3,n3a,nadj)
      call increment_timer('greenf3d',1)
!     call system_clock(count=iticks3)
!     if(idproc.eq.0)write(6,*)'returned from greenf'
! FFT the Green function:
      scale=1.0
!     write(6,*)'(greenf) calling fft3dhpf; grn='
      iunity=1
      izero=0
      scale=1.d0
      ijsign=1
!     if(idproc.eq.0)then
!     write(6,*)'calling dcft w/ n1,n2*n1,n3a,n2*n3a,n1,n2,n3a='
!     write(6,*)n1,n2*n1
!     write(6,*)n3a,n2*n3a
!     write(6,*)n1,n2,n3a
!     endif
      call drcft3(rho2,n1+2,n2*(n1+2),                                     &
     & grnxtr,n1/2+1,n2*(n1/2+1),n1,n2,n3a,ijsign,scale,aux,naux)
!     call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,grnxtr)
!     if(idproc.eq.0)write(6,*)'rtrnd frm dcft3 grnf'
!
! deposit charge on the grid:
      call increment_timer('rhoslo3d',0)
!     if(idproc.eq.0)write(6,*)'calling rhoslo3d'
      call rhoslo3d(c,rho,msk,np,ab,de,gh,indx,jndx,kndx,
     &indxp1,jndxp1,kndxp1,nx,ny,nz,nadj)
!     if(idproc.eq.0)write(6,*)'returned from rhoslo3d'
      call MPI_ALLREDUCE(rho,rhosum,nx*ny*nz,mreal,mpisum,lworld,ierror)
!     if(idproc.eq.0)write(6,*)'returned from MPI_ALLREDUCE'
      rho=rhosum
      call increment_timer('rhoslo3d',1)
!     glrhochk=sum(rho)
!     if(idproc.eq.0)write(6,*)'global sum of rho = ',glrhochk
cryne august 1, 2002
cryne moved normalization out of rhoslo to here:
      rho=rho/ntot
c     gnrhochk=sum(rho)
c     if(iverbose.ge.0)then
c       write(6,*)'normalized global rhosum=',gnrhochk
c     endif
c---------------
      idebug=0
      hxyzi=hxi*hyi*hzi
      if(idebug.eq.1)then
      write(6,*)'writing charge density (note mult by hxi*hyi*hzi)'
      do i=1,nx
      write(65,*)xmin+(i-1)*hx,rho(i,ny/2,nz/2)*hxi*hyi*hzi
      enddo
      endif
c---------------
      call system_clock(count=iticks1)
! store rho in lower left quadrant of doubled grid:
      do k=1,n3a
      do j=1,n2
      do i=1,n1
      rho2(i,j,k)=0.d0
      enddo
      enddo
      enddo
cryne forall(i=1:nx,j=1:ny,k=1:nz)rho2(i,j,k)=cmplx(rho(i,j,k),0.)
      do k=1,nz
      do j=1,ny
      do i=1,nx
      rho2(i,j,k)=rho(i,j,k)
      enddo
      enddo
      enddo
!     if(idproc.eq.0)then
!     write(6,*)'finished storing in lower left quadrant; start conv'
!     endif
!     call system_clock(count=iticks2)
!---------convolution------------
!
! fft the charge density:
      ijsign=1
!     if(idproc.eq.0)then
!     write(6,*)'(rho)calling dcft w/ n1,n2*n1,n3a,n2*n3a,n1,n2,n3a='
!     write(6,*)n1,n2*n1
!     write(6,*)n3a,n2*n3a
!     write(6,*)n1,n2,n3a
!     endif
      call drcft3(rho2,n1+2,n2*(n1+2),                                   &
     & rho2xtr,n1/2+1,n2*(n1/2+1),n1,n2,n3a,ijsign,scale,aux,naux)
!     call fft3dhpf(n1,n2,n3a,1,scale,0,nadj,rho2,rho2xtr)
!     call system_clock(count=iticks4)
!     if(idproc.eq.0)write(6,*)'returned from dcft3 of rho'
! multiply transformed charge density and transformed Green function:
!     rho2xtr=rho2xtr*grnxtr/(n1*n2*n3a)
      rho2xtr=rho2xtr*grnxtr
      scale=1.d0/(n1*n2*n3a)
! inverse fft:
!     if(idproc.eq.0)then
!     write(6,*)'calling inv. dcft w/ n1/2+1,n2*(n1/2+1),n1,n2,n3a='
!     write(6,*)n1/2+1,n2*(n1/2+1)
!     write(6,*)n1,n2,n3a
!     endif
      ijsign=-1
      call dcrft3(rho2xtr,n1/2+1,n2*(n1/2+1),                            &
     & rho2,n1+2,n2*(n1+2),n1,n2,n3a,ijsign,scale,aux,naux)
!!!! & rho2,n1,n2*n1,n1/2+1,n2,n3a,ijsign,scale,aux,naux)
!     if(idproc.eq.0)write(6,*)'returned from inverse dcft3'
!     call fft3dhpf(n3a,n2,n1,-1,scale,0,nadj,rho2xtr,rho2)
!     call system_clock(count=iticks5)
!!!!  rho2=rho2*hx*hy*hz
!----done with convolution-------
! store physical data back on grid of correct (not doubled) size:
cryne forall(i=1:nx,j=1:ny,k=1:nz)rho(i,j,k)=real(rho2(i,j,k))
      do k=1,nz
      do j=1,ny
      do i=1,nx
      rho(i,j,k)=rho2(i,j,k)
      enddo
      enddo
      enddo
c---------------
      idebug=0
      if(idebug.eq.1)then
      write(6,*)'writing out potential'
      do i=1,nx
      write(66,*)xmin+(i-1)*hx,rho(i,ny/2,nz/2)
      enddo
      write(6,*)'stopping due to debug statement in spch3d'
      call myexit
      endif
c---------------
!     if(idproc.eq.0)write(6,*)'obtaining electric fields'
! obtain the electric fields:
      exg=0.5*hxi*(cshift(rho,-1,1)-cshift(rho,1,1))
      eyg=0.5*hyi*(cshift(rho,-1,2)-cshift(rho,1,2))
      ezg=0.5*hzi*(cshift(rho,-1,3)-cshift(rho,1,3))
      call system_clock(count=iticks6)
!     if(idproc.eq.0)write(6,*)'interpolating electric fields'
! interpolate electric field at particle postions:
      call increment_timer('ntrslo3d',0)
      call ntrslo3d(exg,eyg,ezg,ex,ey,ez,msk,ab,de,gh,indx,jndx,kndx,
     &indxp1,jndxp1,kndxp1,nx,ny,nz,np)
      call increment_timer('ntrslo3d',1)
!     if(idproc.eq.0)write(6,*)'done interpolating electric fields'
!     call system_clock(count=iticks7)
!     t10=(iticks1-iticks0)/hertz
!     t21=(iticks2-iticks1)/hertz
!     t32=(iticks3-iticks2)/hertz
!     t43=(iticks4-iticks3)/hertz
!     t54=(iticks5-iticks4)/hertz
!     t65=(iticks6-iticks5)/hertz
!     t76=(iticks7-iticks6)/hertz
!     t70=(iticks7-iticks0)/hertz
!     at10=at10+t10
!     at21=at21+t21
!     at32=at32+t32
!     at43=at43+t43
!     at54=at54+t54
!     at65=at65+t65
!     at76=at76+t76
!     at70=at70+t70
!     write(3,2468)t10,t21,t32,t43,t54,t65,t76,t70
!     write(4,2468)at10,at21,at32,at43,at54,at65,at76,at70
!2468 format(8(1pe9.3,1x))
!     call flush_(3)
!     call flush_(4)
      return
      end

      subroutine greenf3d(g,nx,ny,nz,n1,n2,n3,n3a,nadj)
      implicit double precision(a-h,o-z)
      real*8 g
!     dimension g(n1,n2,n3a)
      dimension g(n1+2,n2,n3a)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!     write(6,*)'inside greenf3d'
!     write(6,*)'nx,ny,nz=',nx,ny,nz
!     write(6,*)'n1,n2,n3,n3a=',n1,n2,n3,n3a
!     write(6,*)'nadj=',nadj
      if(nadj.eq.1)goto 50
!     write(6,*)'(greenf): isolated boundary conditions'
! zero adjacent bunches (isolated boundary conditions):
      do k=1,nz+1
      do j=1,ny+1
      do i=1,nx+1
      g(i,j,k)   =   (hx*(i-1))**2   +(hy*(j-1))**2   +(hz*(k-1))**2
      enddo
      enddo
      enddo
      g(1,1,1)=1.
      g(1:nx+1,1:ny+1,1:nz+1)=1.d0/sqrt(g(1:nx+1,1:ny+1,1:nz+1))
      g(1,1,1)=g(1,1,2)
      do k=1,nz
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j,k)=g(n1-i+2,j,k)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j,k)=g(i,n2-j+2,k)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1,ny
      do i=1,nx
      g(i,j,k)=g(i,j,n3-k+2)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j,k)=g(i,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j,k)=g(n1-i+2,j,n3-k+2)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j,k)=g(n1-i+2,n2-j+2,k)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j,k)=g(n1-i+2,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      goto 200
! adjacent bunches (periodic boundary conditions):
   50 continue
      d=hz*nz
      do 100 k=-nz/2,nz/2-1
      do 99 j=-ny,ny-1
      do 98 i=-nx,nx-1
      if(i.eq.0 .and. j.eq.0 .and. k.eq.0)goto 98
      g(1+mod(i+n1,n1),1+mod(j+n2,n2),1+mod(k+nz,nz))=                  &
     & 1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k)**2)                         &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+1.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+2.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+3.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+4.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+5.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+6.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+7.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k+8.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-1.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-2.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-3.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-4.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-5.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-6.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-7.*d)**2)                    &
     &+1.0/sqrt( (hx*i)**2+(hy*j)**2+(hz*k-8.*d)**2)
   98 continue
   99 continue
  100 continue
      g(1,1,1)=g(1,1,2)
  200 continue
!     write(6,*)'leaving greenf3d'
      return
      end

      subroutine rhoslo3d(coord,rho,msk,np,ab,de,gh,indx,jndx,kndx,
     &indxp1,jndxp1,kndxp1,nx,ny,nz,nadj)
cryne 08/24/2001      use hpf_library
      implicit double precision(a-h,o-z)
      logical msk
      dimension coord(6,np),msk(np),vol(np)
!hpf$ distribute coord(*,block)
!hpf$ align (:) with coord(*,:) :: msk,vol
      dimension ab(np),de(np),gh(np),indx(np),jndx(np),kndx(np),
     &indxp1(np),jndxp1(np),kndxp1(np)
!hpf$ align (:) with coord(*,:) :: ab,de,gh,indx,jndx,kndx
!hpf$ align (:) with coord(*,:) :: indxp1,jndxp1,kndxp1
      dimension rho(nx,ny,nz),tmp(nx,ny,nz)
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: tmp
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/showme/iverbose
      if(idproc.eq.0) write(6,*)'inside rhoslo3d'
!       write(6,*)'xmin,xmax=',xmin,xmax
!       write(6,*)'nx,ny,nz=',nx,ny,nz
!       write(6,*)'nadj=',nadj
!       xminnn=minval(coord(1,:))
!       xmaxxx=maxval(coord(1,:))
!       write(6,*)'xminnn,xmaxxx=',xminnn,xmaxxx
cryne 6/25/2002      do i=1,np
      do j=1,np
      indx(j)=(coord(1,j)-xmin)*hxi + 1
      jndx(j)=(coord(3,j)-ymin)*hyi + 1
      kndx(j)=(coord(5,j)-zmin)*hzi + 1
      enddo
      do j=1,np
      indxp1(j)=indx(j)+1
      jndxp1(j)=jndx(j)+1
      kndxp1(j)=kndx(j)+1
      enddo
!-------
      imin=minval(indx,1,msk)
      imax=maxval(indx,1,msk)
      jmin=minval(jndx,1,msk)
      jmax=maxval(jndx,1,msk)
      kmin=minval(kndx,1,msk)
      kmax=maxval(kndx,1,msk)
      if((imin.lt.1).or.(imax.gt.nx-1))then
        write(6,*)'error in rhoslo3d: imin,imax=',imin,imax
        call myexit
      endif
      if((jmin.lt.1).or.(jmax.gt.ny-1))then
        write(6,*)'error in rhoslo3d: jmin,jmax=',jmin,jmax
        call myexit
      endif
      if(nadj.eq.0)then
       if((kmin.lt.1).or.(kmax.gt.nz-1))then
        write(6,*)'error in rhoslo3d (nadj=0): kmin,kmax=',kmin,kmax
        call myexit
       endif
      endif
      if(nadj.eq.1)then
       if((kmin.lt.1).or.(kmax.gt.nz))then
        write(6,*)'error in rhoslo3d (nadj=1): kmin,kmax=',kmin,kmax
        call myexit
       endif
      endif
!-------
      if(nadj.eq.1)then
        do n=1,np
        if(kndxp1(n).eq.nz+1)kndxp1(n)=1
        enddo
      endif
      ab=((xmin-coord(1,:))+indx*hx)*hxi
      de=((ymin-coord(3,:))+jndx*hy)*hyi
      gh=((zmin-coord(5,:))+kndx*hz)*hzi
      rho=0.
!1 (i,j,k):
      vol=ab*de*gh
      do 100 n=1,np
      rho(indx(n),jndx(n),kndx(n))=                                     &
     &rho(indx(n),jndx(n),kndx(n))+vol(n)
  100 continue
!2 (i,j+1,k):
      vol=ab*(1.-de)*gh
      do 200 n=1,np
      rho(indx(n),jndxp1(n),kndx(n))=                                   &
     &rho(indx(n),jndxp1(n),kndx(n))+vol(n)
  200 continue
!3 (i,j+1,k+1):
      vol=ab*(1.-de)*(1.-gh)
      do 300 n=1,np
      rho(indx(n),jndxp1(n),kndxp1(n))=                                 &
     &rho(indx(n),jndxp1(n),kndxp1(n))+vol(n)
  300 continue
!4 (i,j,k+1):
      vol=ab*de*(1.-gh)
      do 400 n=1,np
      rho(indx(n),jndx(n),kndxp1(n))=                                   &
     &rho(indx(n),jndx(n),kndxp1(n))+vol(n)
  400 continue
!5 (i+1,j,k+1):
      vol=(1.-ab)*de*(1.-gh)
      do 500 n=1,np
      rho(indxp1(n),jndx(n),kndxp1(n))=                                 &
     &rho(indxp1(n),jndx(n),kndxp1(n))+vol(n)
  500 continue
!6 (i+1,j+1,k+1):
      vol=(1.-ab)*(1.-de)*(1.-gh)
      do 600 n=1,np
      rho(indxp1(n),jndxp1(n),kndxp1(n))=                               &
     &rho(indxp1(n),jndxp1(n),kndxp1(n))+vol(n)
  600 continue
!7 (i+1,j+1,k):
      vol=(1.-ab)*(1.-de)*gh
      do 700 n=1,np
      rho(indxp1(n),jndxp1(n),kndx(n))=                                 &
     &rho(indxp1(n),jndxp1(n),kndx(n))+vol(n)
  700 continue
!8 (i+1,j,k):
      vol=(1.-ab)*de*gh
      do 800 n=1,np
      rho(indxp1(n),jndx(n),kndx(n))=                                   &
     &rho(indxp1(n),jndx(n),kndx(n))+vol(n)
  800 continue
!
cryne august 1, 2002:
cccc      ngood=count(msk)
cccc      write(6,*)'ngood=',ngood
cccc      rho=rho/ngood
!wrong     rho=rho/ngood*hxi*hyi*hzi
!     if(idproc.eq.0)write(6,*)'leaving rhoslo3d'
!     rhochk=sum(rho)
!     write(6,*)'[rhoslo3d]sum(rho)=',rhochk
      return
      end

      subroutine ntrslo3d(exg,eyg,ezg,ex,ey,ez,msk,ab,de,gh,
     &indx,jndx,kndx,indxp1,jndxp1,kndxp1,nx,ny,nz,np)
      use parallel
      implicit double precision(a-h,o-z)
      logical msk
      dimension exg(nx,ny,nz),eyg(nx,ny,nz),ezg(nx,ny,nz)
!hpf$ distribute exg(*,*,block)
!hpf$ align (*,*,:) with exg(*,*,:) :: eyg,ezg
      dimension ex(np),ey(np),ez(np),msk(np)
      dimension ab(np),de(np),gh(np),indx(np),jndx(np),kndx(np),
     &indxp1(np),jndxp1(np),kndxp1(np)
!hpf$ template t1(np)
!hpf$ distribute t1(block)
!hpf$ align (:) with t1(:) :: ex,ey,ez,msk,ab,de,gh
!hpf$ align (:) with t1(:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!     if(idproc.eq.0)write(6,*)'inside ntrslo3d'
!     if(idproc.eq.0)then
!       do n=1,np
!         if(indx(n).lt.1 .or. indx(n).gt.nx)write(6,*)'error: indx(n)'
!       enddo
!       do n=1,np
!         if(jndx(n).lt.1 .or. jndx(n).gt.ny)write(6,*)'error: jndx(n)'
!       enddo
!       do n=1,np
!         if(kndx(n).lt.1 .or. kndx(n).gt.nz)write(6,*)'error: kndx(n)'
!       enddo
!       do n=1,np
!         if(indxp1(n).lt.1.or.indxp1(n).gt.nx)write(6,*)'err:indxp1(n)'
!       enddo
!       do n=1,np
!         if(jndxp1(n).lt.1.or.jndxp1(n).gt.ny)write(6,*)'err:jndxp1(n)'
!       enddo
!       do n=1,np
!         if(kndxp1(n).lt.1.or.kndxp1(n).gt.nz)write(6,*)'err:kndxp1(n)'
!       enddo
!     endif
cryne forall(n=1:np)ex(n)=
      do 100 n=1,np
      ex(n)=                                                            &
     & exg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     &+exg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     &+exg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     &+exg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     &+exg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     &+exg(indxp1(n),
     &jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     &+exg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     &+exg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)
  100 continue

cryne forall(n=1:np)ey(n)=
      do 200 n=1,np
      ey(n)=                                                            &
     & eyg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     &+eyg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     &+eyg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     &+eyg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     &+eyg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     &+eyg(indxp1(n),
     &jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     &+eyg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     &+eyg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)
  200 continue

cryne forall(n=1:np)ez(n)=
      do 300 n=1,np
      ez(n)=                                                            &
     & ezg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     &+ezg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     &+ezg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     &+ezg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     &+ezg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     &+ezg(indxp1(n),
     &jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     &+ezg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     &+ezg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)
  300 continue
!     if(idproc.eq.0)write(6,*)'leaving ntrslo3d'
      return
      end
c
c     block data etimes
c     common/accum/at10,at21,at32,at43,at54,at65,at76,at70
c     data at10,at21,at32,at43,at54,at65,at76,at70/0.,0.,0.,0.,0.,0.,0.,&
c    &     0./
c     end
