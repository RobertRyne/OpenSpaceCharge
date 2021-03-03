      subroutine setbound(c6,msk,np,gblam,gamma,nx0,ny0,nz0,noresize,   &
     &                    nadj0)
      use rays
      use parallel
      include 'impli.inc'
      real*8 hmax
      logical msk
      dimension c6(6,np),msk(np)
!hpf$ distribute (*,block) :: c6
!hpf$ align (:) with c6(*,:) :: msk
      real*8          Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Hx,Hy,Hz,Hxi,Hyi,Hzi
      common/GRIDSZ3D/Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Hx,Hy,Hz,Hxi,Hyi,Hzi
      common/gxyzsave/xmin0,xmax0,ymin0,ymax0,zmin0,zmax0,kfixbdy,madegr
      integer           IDirectFieldCalc,IDensityFunction,ISolve
     &                 ,AnagPatchSize,AnagRefRatio
      common/NEWPOISSON/IDirectFieldCalc,IDensityFunction,ISolve
     &                 ,AnagPatchSize,AnagRefRatio
      integer       IVerbose
      common/SHOWME/IVerbose
      n0min=16
      n0max=1024
   37 continue
c     if(idproc.eq.0)then
c       write(6,*)'inside setbound; gamma=',gamma
c     endif
      call boundp3d(c6,msk,np,gblam,nx0,ny0,nz0,nadj0)
c     if(idproc.eq.0)then
c     write(6,*)'returned from boundp3d'
c     write(6,*)'xmin,xmax=',xmin,xmax
c     write(6,*)'ymin,ymax=',ymin,ymax
c     write(6,*)'zmin,zmax=',zmin,zmax
c     endif

! special case for ANAG Infinite Domain Poisson solver which requires
! uniform grid spacing, only applicable in variable grid case
      if( kfixbdy .NE. 1 .AND. ISolve/10 .EQ. 1 )then
        ! set all h's to max and adjust the physical domain boundaries
        ! without moving the center of the domain
        hmax = MAX( hx,hy,hz )
        if( IVerbose .GT. 4 .AND. IdProc .EQ. 0 .AND.
     &      (hmax.NE.hx .OR. hmax.NE.hy .OR. hmax.NE.hz) )then
          write(6,*) 'info: SETBOUND: enforcing isotropic grid spacing '
     &              ,'for ANAG Poisson solver; old h[xyz] = ',hx,hy,hz
     &              ,', new h = ',hmax
        endif
        if( hx .LT. hmax )then
          ! new_xmin is center minus 1/2 new X domain length
          hx = hmax
          xmin0 = 0.5 * ( xmax0 + xmin0 ) - 0.5 * ( hmax * (nx0-1) )
          xmax0 = xmin0 + ( hmax * (nx0-1) )
        endif
        if( hy .LT. hmax )then
          hy = hmax
          ymin0 = 0.5 * ( ymax0 + ymin0 ) - 0.5 * ( hmax * (ny0-1) )
          ymax0 = ymin0 + ( hmax * (ny0-1) )
        endif
        if( hz .LT. hmax )then
          hz = hmax
          zmin0 = 0.5 * ( zmax0 + zmin0 ) - 0.5 * ( hmax * (nz0-1) )
          zmax0 = zmin0 + ( hmax * (nz0-1) )
        endif
        hxi = 1.0d0 / hx
        hyi = 1.0d0 / hy
        hzi = 1.0d0 / hz
      endif

c===================================================================
c 12/4/2002 new code to allow for fixed grid size:
      if(kfixbdy.eq.1)then
cryne 11/27/03        if(idproc.eq.0)write(6,*)'kfixbdy.eq.1'
        if(xmin0.lt.xmin)then
          xmin=xmin0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside xmin; halting.'
          write(6,*)'xmin0,xmin=',xmin0,xmin
          endif
          call myexit
        endif
        if(xmax0.gt.xmax)then
          xmax=xmax0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside xmax; halting.'
          write(6,*)'xmax0,xmax=',xmax0,xmax
          endif
          call myexit
        endif
        if(ymin0.lt.ymin)then
          ymin=ymin0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside ymin; halting.'
          write(6,*)'ymin0,ymin=',ymin0,ymin
          endif
          call myexit
        endif
        if(ymax0.gt.ymax)then
          ymax=ymax0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside ymax; halting.'
          write(6,*)'ymax0,ymax=',ymax0,ymax
          endif
          call myexit
        endif
        hx=(xmax-xmin)/(nx0-1)
        hy=(ymax-ymin)/(ny0-1)
        hxi=1./hx
        hyi=1./hy
      endif
      if(kfixbdy.eq.1 .and. nadj0.eq.0)then
!       if(idproc.eq.0)write(6,*)'kfixbdy.eq.1 .and. nadj0.eq.0'
cryne 11/27/03 from now on, zmin is the grid size IN THE BUNCH FRAME,
cryne          WHERE THE POISSON EQUATION IS ACTUALLY SOLVED
cryne        if(gamma*zmin0.lt.zmin)then
cryne          zmin=gamma*zmin0
        if(zmin0.lt.zmin)then
          zmin=zmin0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside zmin; halting.'
cryne          write(6,*)'zmin0,zmin=',gamma*zmin0,zmin
          write(6,*)'zmin0,zmin=',zmin0,zmin
          endif
          call myexit
        endif
cryne        if(gamma*zmax0.gt.zmax)then
cryne          zmax=gamma*zmax0
        if(zmax0.gt.zmax)then
          zmax=zmax0
        else
          if(idproc.eq.0)then
          write(6,*)'ERROR: particle range falls outside zmax; halting.'
cryne          write(6,*)'zmax0,zmax=',gamma*zmax0,zmax
          write(6,*)'zmax0,zmax=',zmax0,zmax
          endif
          call myexit
        endif
        hz=(zmax-zmin)/(nz0-1)
        hzi=1./hz
      endif
!===================================================================
! check grid sizes:
! aspect is the max allowed ratio between any *two* grid sizes
      aspect=1.9
      aspinv=1./aspect
      ierrx=0
      ierry=0
      ierrz=0
      if(hx.lt.hy .and. hx.lt.hz)then
        if(hx.lt.aspinv*min(hy,hz))ierrx=1
      else if(hy.lt.hx .and. hy.lt.hz)then
        if(hy.lt.aspinv*min(hx,hz))ierry=1
      else if(hz.lt.hx .and. hz.lt.hy)then
        if(hz.lt.aspinv*min(hx,hy))ierrz=1
      endif
      if(ierrx.eq.1 .and. idproc.eq.0)write(12,1111)hx,hy,hz
      if(ierry.eq.1 .and. idproc.eq.0)write(12,1112)hx,hy,hz
      if(ierrz.eq.1 .and. idproc.eq.0)write(12,1113)hx,hy,hz
 1111 format('(aspect ratio)should increase hx; hx,hy,hz=',3(1pe9.3,1x))
 1112 format('(aspect ratio)should increase hy; hx,hy,hz=',3(1pe9.3,1x))
 1113 format('(aspect ratio)should increase hz; hx,hy,hz=',3(1pe9.3,1x))
      if(hx.gt.hy .and. hx.gt.hz)then
        if(hx.gt.aspect*max(hy,hz))ierrx=-1
      else if(hy.gt.hx .and. hy.gt.hz)then
        if(hy.gt.aspect*max(hx,hz))ierry=-1
      else if(hz.gt.hx .and. hz.gt.hy)then
        if(hz.gt.aspect*max(hx,hy))ierrz=-1
      endif
      if(ierrx.eq.-1 .and. idproc.eq.0)write(12,1114)hx,hy,hz
      if(ierry.eq.-1 .and. idproc.eq.0)write(12,1115)hx,hy,hz
      if(ierrz.eq.-1 .and. idproc.eq.0)write(12,1116)hx,hy,hz
 1114 format('(aspect ratio)should decrease hx; hx,hy,hz=',3(1pe9.3,1x))
 1115 format('(aspect ratio)should decrease hy; hx,hy,hz=',3(1pe9.3,1x))
 1116 format('(aspect ratio)should decrease hz; hx,hy,hz=',3(1pe9.3,1x))
c flush file 12 just in case the code crashes or runs out of time:
      call myflush(12)
c===================================================================
c     if(ierrx.ne.0 .or. ierry.ne.0 .or. ierrz.ne.0)
c    #write(6,*)'warning from setbound: ierr (x,y,z)=',ierrx,ierry,ierrz
        if(noresize.eq.1 .or. kfixbdy.eq.1)goto 38
c resize the grid if needed:
        if((ierrx.eq.1 .or. ierry.eq.1 .or. ierrz.eq.1) .and.
     &     (nx0.le.n0min .or. ny0.le.n0min .or. nz0.le.n0min))then
           write(6,*)'cannot decrease grid size; doubling all sizes'
           nx0=2*nx0
           ny0=2*ny0
           nz0=2*nz0
           goto 37
        endif
        if(ierrx.eq.1)then
          nx0=nx0/2
          if(idproc.eq.0)write(6,*)'#DECREASING NX0 to ',nx0
        else if(ierry.eq.1)then
          ny0=ny0/2
          if(idproc.eq.0)write(6,*)'#DECREASING NY0 to ',ny0
        else if(ierrz.eq.1)then
          nz0=nz0/2
          if(idproc.eq.0)write(6,*)'#DECREASING NZ0 to ',nz0
        endif
        if((ierrx.eq.-1 .or. ierry.eq.-1 .or. ierrz.eq.-1) .and.
     &     (nx0.ge.n0max .or. ny0.ge.n0max .or. nz0.ge.n0max))then
           write(6,*)'cannot increase grid size; stopping.'
           stop
        endif
        if(ierrx.eq.-1)then
          nx0=2*nx0
          if(idproc.eq.0)write(6,*)'#INCREASING NX0 to ',nx0
        else if(ierry.eq.-1)then
          ny0=2*ny0
          if(idproc.eq.0)write(6,*)'#INCREASING NY0 to ',ny0
        else if(ierrz.eq.-1)then
          nz0=2*nz0
          if(idproc.eq.0)write(6,*)'#INCREASING NZ0 to ',nz0
        endif
        if(abs(ierrx)+abs(ierry)+abs(ierrz).ne.0)then
          write(6,'(3i5,9x,3i5)')nx0,ny0,nz0
          goto 37
        endif
   38 continue
c==========
      if(idproc.eq.0.and.iverbose.ge.5)then
        write(6,*)'setbound: [xyz]min=',xmin,ymin,zmin
     &           ,' [xyz]max=',xmax,ymax,zmax, ' h[xyz]='
     &           ,hx,hy,hz
      endif
      return
      end
