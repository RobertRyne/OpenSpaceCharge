!***********************************************************************
!
! IMPACT 3D space charge routines
! Copyright 2001 University of California
!
!     subroutine SPCH3D(c,ex,ey,ez,msk,Np,Ntot,
!    &                  Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,rparams,cparams)
! Arguments in SPCH3D(...)
!   c        in:out (1:5:2,*) are x,y,z coordinates of particles
!     NB: in the rest of MLI, col.5 denotes time-of-flight (i.e., phase),
!         but prior to calling the space-charge routines, it is
!         converted to longitudinal position z
!                   (2:6:4,*) are momenta and not used here
!   ex,ey,ez :out electric field at particles
!   msk      :in  flag indicating valid particles
!   Np       :in  number of particles (good and bad) on this processor
!   Ntot     :in  number of (good) particles on all processors
!   Nx,Ny,Nz :in  dimensions of grid on the regular sized grid
!   N1,N2,N3 :in  dimensions of the bigger (usually doubled) grid
!   N3a      :in  =Nz for longitudinal periodic BCs, =2*Nz for open
!   Nadj     :in  =0 for long. open or Dirichlet BCs, >=1 for periodic
!
! Globals:
!   ISolve   :in  flag denoting Poisson solver to use
!                   1 for default Infinite Domain
!                   1x for ANAG ID
!                   2x for ANAG Homogenous Dirichlet
!                   30 for Chombo AMR Homogenous Dirichlet
!                   4x for Chombo AMR Infinite Domain
!                  here "x" denotes the discretization to use
!                    0 for Spectral (MLI-style)
!                    1 for Laplace (Chombo-style 7-point)
!                    2 for Mehrstellen O(h^6)
!                    3 for Mehrstellen O(h^4)
!
!***********************************************************************
      subroutine spch3d(c,ex,ey,ez,msk,np,ntot,                         &
     &                  nx,ny,nz,n1,n2,n3,n3a,nadj,rparams,cparams)
      use parallel
      implicit none
!Arguments
      integer          np,ntot,nx,ny,nz,n1,n2,n3,n3a,nadj
      double precision c(6,np)
      double precision ex(np),ey(np),ez(np)
      logical          msk(np)
      double precision rparams
      character*16 cparams
      dimension rparams(60),cparams(60)
!Local variables
      double precision phi(nx,ny,nz)
      integer i,j,k ,ip
!Globals
      integer       iverbose
      common/showme/iverbose
      integer           idirectfieldcalc,idensityfunction,isolve,idum
      common/newpoisson/idirectfieldcalc,idensityfunction,isolve,idum(2)   
      integer          idirich
      common/dirichlet/idirich
      doubleprecision xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!Externals

!Execute

!-----------------------------------------------------------------------

! save particle data in real coordinates for the Chombo programs
      if (iverbose.GT.9.AND.idproc.EQ.0) then
        write(6,*) 'Writing particle coords to file MLI.Pxyz.dat ...'
        open(90,file='MLI.Pxyz.dat')
        write(90,9001) ((c(i,ip),i=1,5,2),ip=1,Np)
        close(90)
 9001   format(3(1x,E16.9))
      endif

!-----------------------------------------------------------------------

!dbs -- added other solvers and formatted output of phi
      if( iverbose .ge. 5 .AND. idproc .eq. 0 )
     &    write(6,*) 'in SPCH3D with solver,dirich = ',isolve,idirich

      if( isolve .ge. 1 .AND. isolve .lt. 10 )then
        ! default open BC (infinite domain) solver
        ![NOTE: the value of isolve isnt important because SPCH3D1() 
        !       ignores the fft_stencil input keyword because it only
        !       implements the spectral stencil.]
        if( idirich .ne. 0 )then
          write(6,*) 'error: SPCH3D1: Dirichlet BC not supported'
          stop 'SPCH3D1a'
        endif
        if( iverbose .ge. 6 .AND. idproc .eq. 0 )
     &      write(6,*) 'SPCH3D: calling SPCH3D1'
        call SPCH3D1( c ,ex,ey,ez ,msk                                  &
     &               ,np,ntot,nx,ny,nz,n1,n2,n3,n3a,nadj,phi )

      elseif( isolve .ge. 10 .and. isolve .lt. 20 )then
        ! ANAG infinite domain solver
        if( idirectfieldcalc .ne. 0 )then
          write(6,*) 'error: SPCH3D2: solving_for=E not supported'
          stop 'SPCH3D2a'
        endif
        if( idirich .ne. 0 )then
          write(6,*) 'error: SPCH3D2: Dirichlet BC not supported'
          stop 'SPCH3D2b'
        endif
        if( iverbose .ge. 6 .AND. idproc .eq. 0 )
     &      write(6,*) 'SPCH3D: calling SPCH3D2'
        call SPCH3D2( c ,ex,ey,ez ,msk                                  &
     &               ,np,ntot,nx,ny,nz,n1,n2,n3,n3a,nadj,phi )

      elseif( isolve .ge. 20 .and. isolve .lt. 30 )then
        ! ANAG solver for homogenous Dirichlet BCs
        if( idirectfieldcalc .ne. 0 )then
          write(6,*) 'error: SPCH3DBC0: solving_for=E not supported'
          stop 'SPCH3D0a'
        endif
        if( idirich .EQ. 0 )then
          write(6,*) 'warning: SPCH3D: non-Dirichlet BC is inconsistent'
     &              ,' with SPCH3DBC0 solver.'
          stop 'SPCH3D0b'
        endif
        if( iverbose .ge. 6 .AND. idproc .eq. 0 )
     &      write(6,*) 'SPCH3D: calling SPCH3DBC0'
        call SPCH3DBC0( c ,ex,ey,ez ,msk                                &
     &                 ,np,ntot,nx,ny,nz,n1,n2,n3,n3a,nadj,phi )

      elseif( isolve .ge. 30 .and. isolve .lt. 50 )then
        ! ANAG Chombo AMR solver for infinite domain (open) or homogenous Dirichlet BCs
        if( idirectfieldcalc .ne. 0 )then
          write(6,*) 'error: SPCH3DAMR: solving_for=E not supported'
          stop 'SPCH3DEAb'
        endif
        if( iverbose .ge. 6 .AND. idproc .eq. 0 )
     &      write(6,*) 'SPCH3D: calling SPCH3DAMR'
        call SPCH3DAMR( c ,ex,ey,ez ,msk                                &
     &                 ,np,ntot,nx,ny,nz,n1,n2,n3,n3a,nadj,phi )

      else
        write(6,*) 'error: unknown solver = ',isolve
        stop 'SPCH3D'
      endif

!-----------------------------------------------------------------------

      if( iverbose .ge. 10 .and. idproc .eq. 0 .AND.
     &    idirectfieldcalc .EQ. 0 )then  !no phi is solving_for=E
          write(6,*) 'Writing grid phi data to unit 86 ...'
          write(86,8601)' phi on grid (',nx,ny,nz,'), h=' ,hx,hy,hz
          write(86,*)'  I    J    K     X         Y        Z        phi'
          do k = 1 ,nz
          do j = 1 ,ny
          do i = 1 ,nx
              write(86,8602) i,j,k
     &            ,xmin+(i-1)*hx,ymin+(j-1)*hy,zmin+(k-1)*hz ,phi(i,j,k)
          enddo
          enddo
          enddo
      elseif( iverbose .ge. 8 .and. idproc .eq. 0 )then
          write(6,*) 'Writing grid centerline phi data to unit 86 ...'
          write(86,8601)' phi on grid (',nx,ny,nz,'), h=' ,hx,hy,hz
          write(86,*)'  I    J    K     X         Y        Z        phi'
          i=nx/2+1
          j=ny/2+1
          do k = 1 ,nz
              write(86,8602) i,j,k
     &            ,xmin+(i-1)*hx,ymin+(j-1)*hy,zmin+(k-1)*hz ,phi(i,j,k)
          enddo
          i=nx/2+1
          k=nz/2+1
          do j = 1 ,ny
              write(86,8602) i,j,k
     &            ,xmin+(i-1)*hx,ymin+(j-1)*hy,zmin+(k-1)*hz ,phi(i,j,k)
          enddo
          j=ny/2+1
          k=nz/2+1
          do i = 1 ,nx
              write(86,8602) i,j,k
     &            ,xmin+(i-1)*hx,ymin+(j-1)*hy,zmin+(k-1)*hz ,phi(i,j,k)
          enddo
          write(6,*) 'Writing E on particles 1-10 to unit 87 ...'
          write(87,8701)'p     X      Y      Z      Ex       Ey      Ez'
          do ip = 1,10
            write(87,8702) ip,(c(i,ip),i=1,5,2),ex(ip),ey(ip),ez(ip)
          enddo
      endif
      if( iverbose .ge. 9 .and. idproc .eq. 0 )then
          write(6,*) 'Writing particle XYZ,E data to MLI.Pxyze.dat ...'
          open(90,file='MLI.Pxyze.dat')
          write(90,9002)((c(i,ip),i=1,5,2),ex(ip),ey(ip),ez(ip),ip=1,Np)
          close(90)
      endif
 8601 format(A,1x,3(I3,1x),A,1x,1P,3(E11.5,1x))
 8602 format(1x,3(1X,I3),1P,3(1X,E15.8),4(1X,E15.8))
 8701 format(1x,A)
 8702 format(1x,I2,1P,6(1x,E16.9))
 9002 format(1P,6(1x,E16.9))
!dbs
!Done
      if( IVerbose .GT. 19 ) stop 'SPCHDONE'
      if( iverbose .ge. 5  .AND. idproc .eq. 0 )
     &    write(6,*) 'leaving SPCH3D'

      return
      end
!
!***********************************************************************
      subroutine spch3d1(c,ex,ey,ez,msk,np,ntot,                        &
     &                   nx,ny,nz,n1,n2,n3,n3a,nadj,rho)
      use parallel
      use spchdata 
      use intgreenfn 
      use ml_timer
      implicit double precision(a-h,o-z)
      real*8, dimension(6,np) :: c
      logical, dimension(np) :: msk
      real*8, dimension(np) :: ex,ey,ez
      real*8, dimension(nx,ny,nz) :: rho,exg,eyg,ezg,rhosum
      type (griddata) :: griddims
      integer :: idebug=0
cryne 8/4/2004   integer :: idebug=1
!in module spchdata   complex*16,dimension(n1,n2,n3a) :: rho2
!in module spchdata   complex*16,dimension(n3a,n2,n1) :: grnxtr,rho2xtr

! rho=charge density on the grid
! rho2=...on doubled grid
! rho2xtr=...xformed (by fft) and transposed
! grnxtr=green function, xformed, transposed
! weights, indices associated with area weighting:
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!XXX     common/accum/at10,at21,at32,at43,at54,at65,at76,at70
      common/showme/iverbose
      common/gxyzsave/xmin0,xmax0,ymin0,ymax0,zmin0,zmax0,kfixbdy,madegr
      common/gridxtra/xsml,xbig,ysml,ybig,zsml,zbig
      common/newpoisson/idirectfieldcalc,idensityfunction,idum(3)
      save idebug

!     if(idproc.eq.0)write(6,*)'inside spch3d'
      iexactrho=0

      iunity=1
      munity=-1
      izero=0
      scale=1.d0
      hxyz=hx*hy*hz
      hxyzi=hxi*hyi*hzi
      n123a=n1*n2*n3a

      if (idproc.eq.0.and.iverbose.ge.3) then
        write(6,*) '  In spch3d1() ...'
        write(6,*) '     {np ntot} = {',np,ntot,'}'
        write(6,*) '     {nx ny nz} = {',nx,ny,nz,'}'
        write(6,*) '     {hx hy hz} = {',hx,hy,hz,'}'
        write(6,*) '     {n1 n2 n3} = {',n1,n2,n3,'}'
        write(6,*) '     {hxi hyi hzi} = {',hxi,hyi,hzi,'}'
        write(6,*) '     {n3a nadj} = {',n3a,nadj,'}'
        write(6,*) '     hxyz = ',hxyz
        write(6,*) '     hxyzi = ',hxyzi
        write(6,*) '     nxyz = ',nx*ny*nz
        write(6,*) '     n123 = ',n1*n2*n3
        write(6,*) '     n123a = ',n1*n2*n3a
      end if
!======================================================================
! deposit charge on the grid:
      call increment_timer('rhoslo3d',0)
      call rhoslo3d(c,rho,msk,np,nx,ny,nz,nadj)
      call MPI_ALLREDUCE(rho,rhosum,nx*ny*nz,mreal,mpisum,lworld,ierror)
      rho=rhosum
      call increment_timer('rhoslo3d',1)
!--------
      if(iexactrho.eq.1)then
        glrhochk=sum(rho)
!       if(idproc.eq.0)write(6,*)'(exact rho): global sum = ',glrhochk
        call getrms(c,xrms,yrms,zrms,np)
!       if(idproc.eq.0)write(6,*)'xrms,yrms,zrms=',xrms,yrms,zrms
        xmac=sqrt(5.d0)*xrms
        ymac=sqrt(5.d0)*yrms
        zmac=sqrt(5.d0)*zrms
!       if(idproc.eq.0)write(6,*)'xmac,ymac,zmac=',xmac,ymac,zmac
        rho=0.d0
        do k=1,nz
          z=zmin+(k-1)*hz
          do j=1,ny
            y=ymin+(j-1)*hy
            do i=1,nx
              x=xmin+(i-1)*hx
!       if((x/xbig)**2+(y/ybig)**2+(z/zbig)**2 .le.1.d0)rho(i,j,k)=1.d0
        if((x/xmac)**2+(y/ymac)**2+(z/zmac)**2 .le.1.d0)rho(i,j,k)=1.d0
            enddo
          enddo
        enddo
        glrhonew=sum(rho)
!       if(idproc.eq.0)write(6,*)'glrhonew=',glrhonew
        rho=rho*glrhochk/glrhonew
        glrhochk=sum(rho)
!       if(idproc.eq.0)write(6,*)'new global sum of rho = ',glrhochk
      endif
!ryne august 1, 2002
!ryne moved normalization out of rhoslo to here:
      rho=rho/(hx*hy*hz*ntot)
      call system_clock(count=iticks1)
!--------
      ! print rho on the grid
      if(iverbose.GT.9.AND.idproc.EQ.0)then
          write(6,*) 'Writing grid rho data to unit 85 ...'
          write(85,*)' SPCH3D rho on grid (',nx,ny,nz,'), h=',hx,hy,hz
          write(85,*)'  I    J    K     X         Y        Z        rho'
          do k = 1 ,nz
          do j = 1 ,ny
          do i = 1 ,nx
              write(85,8501) i,j,k
     $            ,xmin+(i-1)*hx,ymin+(j-1)*hy,zmin+(k-1)*hz
     $            ,rho(i,j,k)
          enddo
          enddo
          enddo
 8501     format(1x,3(1X,I3),1P,3(1X,E15.8),1X,E15.8)
      endif
!--------
! keep a copy of rho for diagnostics:
      if(iexactrho.eq.1)rhosum=rho
      checknorm=hx*hy*hz*sum(rhosum)
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
! fft the charge density:
      call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,rho2xtr)
!======================================================================
!HERE IS THE OLD IMPLEMENTATION THAT USES PHI AND CSHIFTS:
      if(idirectfieldcalc.ne.1)then
!       if(idproc.eq.0)write(6,*)'******obtaining scalar potential******'
        call increment_timer('greenf3d',0)
        if(madegr.eq.0 .or. kfixbdy.eq.0)then
!         if (idproc.eq.0) then
!           write(12,*) ' spch3d::spch3d1:'
!           write(12,*) '   computing green function for phi'
!         end if
          griddims%NxIntrvl=nx-1
          griddims%NyIntrvl=ny-1
          griddims%NzIntrvl=nz-1
          griddims%Nx=nx; griddims%Nx2=n1; griddims%hx=hx
          griddims%Ny=ny; griddims%Ny2=n2; griddims%hy=hy
          griddims%Nz=nz; griddims%Nz2=n3; griddims%hz=hz
          griddims%Vh=hxyz
          griddims%Vhinv=hxyzi
          if (idensityfunction.eq.0) then
!           write(12,*) '   using function greenphi'
            call greenphi(rho2,nx,ny,nz,n1,n2,n3,n3a,nadj)
          else if (idensityfunction.eq.1) then
!           write(12,*) '   using function intgreenphi'
            call intgreenphi(griddims,rho2)
          end if
c ===DBG===
!         write(6,*) ' === G_phi array ==='
!         write(6,*) ' {hx,hy,hz} = {',hx,',',hy,',',hz,'}'
!         gmn=rho2(1,1,1)
!         gmx=gmn
!         agmn=abs(gmn)
!         agmx=agmn
!         do iz=1,n1
!           do iy=1,n2
!             do ix=1,n1
!               g=rho2(ix,iy,iz)
!               ag=abs(g)
!               if(g.lt.gmn) gmn=g
!               if(g.gt.gmx) gmx=g
!               if(ag.lt.agmn) agmn=ag
!               if(ag.gt.agmx) agmx=ag
!             end do
!           end do
!         end do
!         write(6,*) ' min(G_phi) =',gmn
!         write(6,*) ' max(G_phi) =',gmx
!         write(6,*) ' min(abs(G_phi)) =',agmn
!         write(6,*) ' max(abs(G_phi)) =',agmx
!         do iz=1,4
!           do iy=1,4
!             do ix=1,4
!               write(6,*) ix,iy,iz,rho2(ix,iy,iz)
!             end do
!           end do
!         end do
!         do iz=1,nz+1
!           do iy=1,ny+1
!             do ix=1,nx+1
!               g=rho2(ix,iy,iz)
!               if (g.lt.0) write(6,'(3(1x,i3),2x,1pe12.5))')           &
!    &            ix-1,iy-1,iz-1,g
!             end do
!           end do
!         end do
c =========
          call increment_timer('greenf3d',1)
          call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,xgrnxtr)
          madegr=1
        else
!         if(idproc.eq.0)write(12,*)'using saved green function'
          call increment_timer('greenf3d',1)
        endif
        grnxtr=rho2xtr*xgrnxtr/n123a
        call fft3dhpf(n3a,n2,n1,munity,scale,izero,nadj,grnxtr,rho2)
        rho(1:nx,1:ny,1:nz)=hxyz*real(rho2(1:nx,1:ny,1:nz))
! obtain the electric fields:
!XXX        exg=0.d0
!XXX        eyg=0.d0
!XXX        ezg=0.d0
        exg=cshift(rho,-1,1)
        exg=exg-cshift(rho,1,1)
        exg=exg*(0.5*hxi)
        eyg=cshift(rho,-1,2)
        eyg=eyg-cshift(rho,1,2)
        eyg=eyg*(0.5*hyi)
!XXX        exg=0.5*hxi*(cshift(rho,-1,1)-cshift(rho,1,1))
!XXX        eyg=0.5*hyi*(cshift(rho,-1,2)-cshift(rho,1,2))
!       ezg=0.5*hzi*(cshift(rho,-1,3)-cshift(rho,1,3))
! the previous statement (ezg=...) causes a crash on seaborg. Replace with: 
        ezg=cshift(rho,-1,3)
        ezg=ezg-cshift(rho,1,3)
        ezg=ezg*(0.5*hzi)
      else
        !=======================================================
        ! New implementation that solves for E directly:
        ! compute the Green function on the grid (use rho2 for storage):
        call increment_timer('greenf3d',0)
        if(madegr.eq.0 .or. kfixbdy.eq.0)then
!         if(idproc.eq.0)write(12,*)'computing green function'
          call greenfx(rho2,nx,ny,nz,n1,n2,n3,n3a,nadj)
          call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,xgrnxtr)
          call greenfy(rho2,nx,ny,nz,n1,n2,n3,n3a,nadj)
          call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,ygrnxtr)
          call greenfz(rho2,nx,ny,nz,n1,n2,n3,n3a,nadj)
          call fft3dhpf(n1,n2,n3a,iunity,scale,izero,nadj,rho2,zgrnxtr)
          madegr=1
        else
!       if(idproc.eq.0)write(12,*)'using saved green function'
        endif
        call increment_timer('greenf3d',1)
        !========================================================
!       if(idproc.eq.0)write(6,*)'starting convolution'
        !---------convolution------------
        grnxtr=rho2xtr*xgrnxtr/n123a
        call fft3dhpf(n3a,n2,n1,munity,scale,izero,nadj,grnxtr,rho2)
        exg(1:nx,1:ny,1:nz)=hxyz*real(rho2(1:nx,1:ny,1:nz))

        grnxtr=rho2xtr*ygrnxtr/n123a
        call fft3dhpf(n3a,n2,n1,munity,scale,izero,nadj,grnxtr,rho2)
        eyg(1:nx,1:ny,1:nz)=hxyz*real(rho2(1:nx,1:ny,1:nz))

        grnxtr=rho2xtr*zgrnxtr/n123a
        call fft3dhpf(n3a,n2,n1,munity,scale,izero,nadj,grnxtr,rho2)
        ezg(1:nx,1:ny,1:nz)=hxyz*real(rho2(1:nx,1:ny,1:nz))
        !----done with convolution-------
      endif

!     if(idproc.eq.0)write(6,*)'done with convolution'
!---------------
      if((idebug.eq.1.or.iverbose.ge.8) .and. idproc.eq.0)then
        write(6,*)'writing out potential and fields'
        do i=1,nx
          j=ny/2+1
          k=nz/2+1
          xval=xmin+(i-1)*hx
          yval=ymin+(j-1)*hy
          zval=zmin+(k-1)*hz
          del2=0.d0
          if(i.gt.1 .and. i.lt.nx)then
            del2=(exg(i+1,j,k)-exg(i-1,j,k))*0.5*hxi+                   &
     &           (eyg(i,j+1,k)-eyg(i,j-1,k))*0.5*hyi+                   &
     &           (ezg(i,j,k+1)-ezg(i,j,k-1))*0.5*hzi
          endif
          write(61,1001)xval,yval,zval,                                 &
     &                 rho(i,j,k),exg(i,j,k),eyg(i,j,k),ezg(i,j,k)      &
     &                ,rhosum(i,j,k),del2
        enddo
        write(61,*)' '
        call myflush(61)
        do j=1,ny
          i=nx/2+1
          k=nz/2+1
          xval=xmin+(i-1)*hx
          yval=ymin+(j-1)*hy
          zval=zmin+(k-1)*hz
          del2=0.d0
          if(j.gt.1 .and. j.lt.ny)then
            del2=(exg(i+1,j,k)-exg(i-1,j,k))*0.5*hxi+                   &
     &           (eyg(i,j+1,k)-eyg(i,j-1,k))*0.5*hyi+                   &
     &           (ezg(i,j,k+1)-ezg(i,j,k-1))*0.5*hzi
          endif
          write(62,1001)xval,yval,zval,                                 &
     &                 rho(i,j,k),exg(i,j,k),eyg(i,j,k),ezg(i,j,k)      &
     &                ,rhosum(i,j,k),del2
        enddo
        write(62,*)' '
        call myflush(62)
        do k=1,nz
          i=nx/2+1
          j=ny/2+1
          xval=xmin+(i-1)*hx
          yval=ymin+(j-1)*hy
          zval=zmin+(k-1)*hz
          del2=0.d0
          if(k.gt.1 .and. k.lt.nz)then
            del2=(exg(i+1,j,k)-exg(i-1,j,k))*0.5*hxi+                   &
     &           (eyg(i,j+1,k)-eyg(i,j-1,k))*0.5*hyi+                   &
     &           (ezg(i,j,k+1)-ezg(i,j,k-1))*0.5*hzi
          endif
          write(63,1001)xval,yval,zval,                                 &
     &                 rho(i,j,k),exg(i,j,k),eyg(i,j,k),ezg(i,j,k)      &
     &                ,rhosum(i,j,k),del2
        enddo
        write(63,*)' '
        call myflush(63)
 1001   format(9(1pe13.6,1x))

      endif
      if( IVerbose .GT. 11 .AND. idproc .EQ. 0 )then
        write(6,*) 'Writing grid E data to unit 88 ...'
        write(88,8601)' phi on grid (',nx,ny,nz,'), h=' ,hx,hy,hz
        write(88,8701)'  I    J   K       X               Y            '
     &            ,'   Z              Ex               '
     &            ,'Ey                 Ez'
        do k = 1 ,nz
        do j = 1 ,ny
        do i = 1 ,nx
          write(88,8602) i,j,k,xmin+(i-1)*hx,ymin+(j-1)*hy
     &                  ,zmin+(k-1)*hz
     &                  ,exg(i,j,k),eyg(i,j,k),ezg(i,j,k)
        enddo
        enddo
        enddo
      endif
 8601 format(A,1x,3(I3,1x),A,1x,1P,3(E11.5,1x))
 8602 format(1x,3(1X,I3),1P,3(1X,E15.8),4(1X,E15.8))
 8701 format(1x,A,A,A)


      if(idebug.eq.1)then
       write(6,*)'PE',idproc,':hxyz,n123a=',hxyz,n123a
       write(6,*)'PE',idproc,':rho2xtr(4,6,8),=',rho2xtr(4,6,8)
       write(6,*)'PE',idproc,':xgrnxtr(4,6,8),=',xgrnxtr(4,6,8)
      endif
      if(idebug.eq.1)idebug=0
!---------------
!
!
!     if(idproc.eq.0)write(6,*)'interpolating electric fields'
! interpolate electric field at particle postions:
      call increment_timer('ntrslo3d',0)
      call ntrslo3d(c,exg,eyg,ezg,ex,ey,ez,msk,nx,ny,nz,np,nadj)
      call increment_timer('ntrslo3d',1)
!     if(idproc.eq.0)write(6,*)'done interpolating; leaving spch3d'
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine rhoslo3d(coord,rho,msk,np,nx,ny,nz,nadj)
cryne 08/24/2001      use hpf_library
      implicit double precision(a-h,o-z)
      logical msk
      dimension coord(6,np),msk(np)
      dimension rho(nx,ny,nz)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/showme/iverbose
!     if(idproc.eq.0) write(6,*)'inside rhoslo3d'
!       write(6,*)'xmin,xmax=',xmin,xmax
!       write(6,*)'nx,ny,nz=',nx,ny,nz
!       write(6,*)'nadj=',nadj
!       xminnn=minval(coord(1,:))
!       xmaxxx=maxval(coord(1,:))
!       write(6,*)'xminnn,xmaxxx=',xminnn,xmaxxx
cryne 6/25/2002      do i=1,np
      rho=0.
      do n=1,np
c INSERT STATEMENT HERE TO SKIP IF MSK(N)=.FALSE.
      indx=(coord(1,n)-xmin)*hxi + 1
      jndx=(coord(3,n)-ymin)*hyi + 1
      kndx=(coord(5,n)-zmin)*hzi + 1
      indxp1=indx+1
      jndxp1=jndx+1
      kndxp1=kndx+1
      ab=((xmin-coord(1,n))+indx*hx)*hxi
      de=((ymin-coord(3,n))+jndx*hy)*hyi
      gh=((zmin-coord(5,n))+kndx*hz)*hzi
!-------
      imin=indx
      imax=indx
      jmin=jndx
      jmax=jndx
      kmin=kndx
      kmax=kndx
      if((imin.lt.1).or.(imax.gt.nx-1))then
        write(6,*)'error in rhoslo3d: imin,imax=',imin,imax
        write(6,*)'nx,xmin,xmax,hx=',nx,xmin,xmax,hx
        write(6,*)'nadj=',nadj
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
        if(kndxp1.eq.nz+1)kndxp1=1
      endif
!1 (i,j,k):
      rho(indx,jndx,kndx)=rho(indx,jndx,kndx)+ab*de*gh
!2 (i,j+1,k):
      rho(indx,jndxp1,kndx)=rho(indx,jndxp1,kndx)+ab*(1.-de)*gh
!3 (i,j+1,k+1):
      rho(indx,jndxp1,kndxp1)=rho(indx,jndxp1,kndxp1)+ab*(1.-de)*(1.-gh)
!4 (i,j,k+1):
      rho(indx,jndx,kndxp1)=rho(indx,jndx,kndxp1)+ab*de*(1.-gh)
!5 (i+1,j,k+1):
      rho(indxp1,jndx,kndxp1)=rho(indxp1,jndx,kndxp1)+(1.-ab)*de*(1.-gh)
!6 (i+1,j+1,k+1):
      rho(indxp1,jndxp1,kndxp1)=                                           &
     &rho(indxp1,jndxp1,kndxp1)+(1.-ab)*(1.-de)*(1.-gh)
!7 (i+1,j+1,k):
      rho(indxp1,jndxp1,kndx)=rho(indxp1,jndxp1,kndx)+(1.-ab)*(1.-de)*gh
!8 (i+1,j,k):
      rho(indxp1,jndx,kndx)=rho(indxp1,jndx,kndx)+(1.-ab)*de*gh
      enddo
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

      subroutine ntrslo3d(coord,exg,eyg,ezg,ex,ey,ez,msk,nx,ny,nz,np,   &
     &                    nadj)
      use parallel
      implicit double precision(a-h,o-z)
      logical msk
      dimension coord(6,np)
      dimension exg(nx,ny,nz),eyg(nx,ny,nz),ezg(nx,ny,nz)
!hpf$ distribute exg(*,*,block)
!hpf$ align (*,*,:) with exg(*,*,:) :: eyg,ezg
      dimension ex(np),ey(np),ez(np),msk(np)
      dimension abz(np),dez(np),ghz(np),indz(np),jndz(np),kndz(np),
     &indzp1(np),jndzp1(np),kndzp1(np)
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
      indx=(coord(1,n)-xmin)*hxi + 1
      jndx=(coord(3,n)-ymin)*hyi + 1
      kndx=(coord(5,n)-zmin)*hzi + 1
      indxp1=indx+1
      jndxp1=jndx+1
      kndxp1=kndx+1
!cryne August 4, 2004 -------
      if(nadj.eq.1)then
      if(kndxp1.eq.nz+1)kndxp1=1
      endif
!cryne -------
      ab=((xmin-coord(1,n))+indx*hx)*hxi
      de=((ymin-coord(3,n))+jndx*hy)*hyi
      gh=((zmin-coord(5,n))+kndx*hz)*hzi
      ex(n)=                                                            &
     & exg(indx,jndx,kndx)*ab*de*gh
     &+exg(indx,jndxp1,kndx)*ab*(1.-de)*gh
     &+exg(indx,jndxp1,kndxp1)*ab*(1.-de)*(1.-gh)
     &+exg(indx,jndx,kndxp1)*ab*de*(1.-gh)
     &+exg(indxp1,jndx,kndxp1)*(1.-ab)*de*(1.-gh)
     &+exg(indxp1,jndxp1,kndxp1)*(1.-ab)*(1.-de)*(1.-gh)
     &+exg(indxp1,jndxp1,kndx)*(1.-ab)*(1.-de)*gh
     &+exg(indxp1,jndx,kndx)*(1.-ab)*de*gh

      ey(n)=                                                            &
     & eyg(indx,jndx,kndx)*ab*de*gh
     &+eyg(indx,jndxp1,kndx)*ab*(1.-de)*gh
     &+eyg(indx,jndxp1,kndxp1)*ab*(1.-de)*(1.-gh)
     &+eyg(indx,jndx,kndxp1)*ab*de*(1.-gh)
     &+eyg(indxp1,jndx,kndxp1)*(1.-ab)*de*(1.-gh)
     &+eyg(indxp1,jndxp1,kndxp1)*(1.-ab)*(1.-de)*(1.-gh)
     &+eyg(indxp1,jndxp1,kndx)*(1.-ab)*(1.-de)*gh
     &+eyg(indxp1,jndx,kndx)*(1.-ab)*de*gh

      ez(n)=                                                            &
     & ezg(indx,jndx,kndx)*ab*de*gh
     &+ezg(indx,jndxp1,kndx)*ab*(1.-de)*gh
     &+ezg(indx,jndxp1,kndxp1)*ab*(1.-de)*(1.-gh)
     &+ezg(indx,jndx,kndxp1)*ab*de*(1.-gh)
     &+ezg(indxp1,jndx,kndxp1)*(1.-ab)*de*(1.-gh)
     &+ezg(indxp1,jndxp1,kndxp1)*(1.-ab)*(1.-de)*(1.-gh)
     &+ezg(indxp1,jndxp1,kndx)*(1.-ab)*(1.-de)*gh
     &+ezg(indxp1,jndx,kndx)*(1.-ab)*de*gh
  100 continue
!     if(idproc.eq.0)write(6,*)'leaving ntrslo3d'
      return
      end
c
c     block data etimes
c     common/accum/at10,at21,at32,at43,at54,at65,at76,at70
c     data at10,at21,at32,at43,at54,at65,at76,at70/0.,0.,0.,0.,0.,0.,0.,&
c    &     0./
c     end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine getrms(cblock,xrms,yrms,zrms,np)
      use beamdata
      use rays
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'map.inc'
      real*8, dimension(6,np) :: cblock
!
      dimension diag(16),rdiag(16)
!     if(1.gt.0)return
! rms quantities:
      den1=1./nrays
      den2=den1*den1
      econ=1.0
      xl=sl
      xk=1./xl
      gambet=gamma*beta
      z=arclen
!     xbar=sum(cblock(1,1:nraysp))*xl*den1
!     ybar=sum(cblock(3,1:nraysp))*xl*den1
!     zbar=sum(cblock(5,1:nraysp))*(beta/xk)*den1
!     sq1=sum(cblock(1,1:nraysp)*cblock(1,1:nraysp))*xl**2
!     sq2=sum(cblock(2,1:nraysp)*cblock(2,1:nraysp))/gambet**2*econ**2
!     sq3=sum(cblock(3,1:nraysp)*cblock(3,1:nraysp))*xl**2
!     sq4=sum(cblock(4,1:nraysp)*cblock(4,1:nraysp))/gambet**2*econ**2
!     xpx=sum(cblock(1,1:nraysp)*cblock(2,1:nraysp))*xl/gambet*econ
!     ypy=sum(cblock(3,1:nraysp)*cblock(4,1:nraysp))*xl/gambet*econ
      diag(1)=sum(cblock(1,1:nraysp))*den1
      diag(2)=sum(cblock(3,1:nraysp))*den1
      diag(3)=sum(cblock(5,1:nraysp))*den1
      diag(4)=sum(cblock(1,1:nraysp)*cblock(1,1:nraysp))
      diag(5)=sum(cblock(2,1:nraysp)*cblock(2,1:nraysp))
      diag(6)=sum(cblock(3,1:nraysp)*cblock(3,1:nraysp))
      diag(7)=sum(cblock(4,1:nraysp)*cblock(4,1:nraysp))
      diag(8)=sum(cblock(5,1:nraysp)*cblock(5,1:nraysp))
      diag(9)=sum(cblock(6,1:nraysp)*cblock(6,1:nraysp))
      diag(10)=sum(cblock(1,1:nraysp)*cblock(2,1:nraysp))
      diag(11)=sum(cblock(3,1:nraysp)*cblock(4,1:nraysp))
      diag(12)=sum(cblock(5,1:nraysp)*cblock(6,1:nraysp))
      call MPI_ALLREDUCE(diag,rdiag,12,mreal,mpisum,lworld,ierror)
c-------
      xbar=rdiag(1)
      ybar=rdiag(2)
      zbar=rdiag(3)
      sq1=rdiag(4)
      sq2=rdiag(5)
      sq3=rdiag(6)
      sq4=rdiag(7)
      sq5=rdiag(8)
      sq6=rdiag(9)
      xpx=rdiag(10)
      ypy=rdiag(11)
      zpz=rdiag(12)
c-------
      epsx2=(sq1*sq2-xpx*xpx)*den2
      epsy2=(sq3*sq4-ypy*ypy)*den2
      epsz2=(sq5*sq6-zpz*zpz)*den2
      xrms=sqrt( sq1*den1 )
      yrms=sqrt( sq3*den1 )
      zrms=sqrt( sq5*den1 )
      pxrms=sqrt( sq2*den1 )
      pyrms=sqrt( sq4*den1 )
      pzrms=sqrt( sq6*den1 )
      zero=0.
      epsx=sqrt(max(epsx2,zero))
      epsy=sqrt(max(epsy2,zero))
      epsz=sqrt(max(epsz2,zero))
      xpx=xpx*den1
      ypy=ypy*den1
      zpz=zpz*den1
      xpxfac=0.
      ypyfac=0.
      zpzfac=0.
      if(xrms.ne.0. .and. pxrms.ne.0.)xpxfac=1./(xrms*pxrms)
      if(yrms.ne.0. .and. pyrms.ne.0.)ypyfac=1./(yrms*pyrms)
      if(zrms.ne.0. .and. pzrms.ne.0.)zpzfac=1./(zrms*pzrms)
      return
      end
c
c--------------------------------------------------------------
c--------------------------------------------------------------

      subroutine greenphi(g,nx,ny,nz,n1,n2,n3,n3a,nadj)
! green function routine.
      implicit double precision(a-h,o-z)
      complex*16 g
      dimension g(n1,n2,n3a)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      gfun(x,y,z)=1.d0/sqrt(x**2+y**2+z**2+1.d-20)
!
!*******
! 3/21/2003 include correct factor of 1/(4pi) in the Green function:
      fourpi=8.d0*(asin(1.d0))
      fourpinv=1.d0/fourpi
!*******
!     write(6,*)'inside greenphi'
      if(nadj.eq.1)goto 50
!     write(6,*)'(greenf): isolated boundary conditions'
! zero adjacent bunches (isolated boundary conditions):
      do k=1,nz+1
      z=(k-1)*hz
      do j=1,ny+1
      y=(j-1)*hy
      do i=1,nx+1
      x=(i-1)*hx
      g(i,j,k)=fourpinv*gfun(x,y,z)
      enddo
      enddo
      enddo
!OTHER OPTIONS POSSIBLE HERE:
!!!   g(1,1,1)=0.d0
      g(1,1,1)=g(1,1,2)
!
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
!     nbunch=200
      nbunch=20
      d=hz*nz
      do 100 k=-nz/2,nz/2-1
      do 99 j=-ny,ny-1
      do 98 i=-nx,nx-1
      if(i.eq.0 .and. j.eq.0 .and. k.eq.0)goto 98
      tmp=0.d0
      do n=-nbunch,nbunch
        tmp=tmp+1.0/(fourpi*sqrt( (hx*i)**2+(hy*j)**2+(hz*k+n*d)**2))
      enddo
      g(1+mod(i+n1,n1),1+mod(j+n2,n2),1+mod(k+nz,nz))=tmp
   98 continue
   99 continue
  100 continue
      g(1,1,1)=g(1,1,2)
  200 continue
!     write(6,*)'leaving greenphi'
      return
      end
c--------------------------------------------------------------
c--------------------------------------------------------------

      subroutine greenfx(g,nx,ny,nz,n1,n2,n3,n3a,nadj)
! green function routine.
      implicit double precision(a-h,o-z)
      complex*16 g
      dimension g(n1,n2,n3a)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      gfun(x,y,z)=x/sqrt(x**2+y**2+z**2+1.d-20)**3
!
!*******
! 3/21/2003 include correct factor of 1/(4pi) in the Green function:
      fourpi=8.d0*(asin(1.d0))
      fourpinv=1.d0/fourpi
!*******
!     write(6,*)'inside greenfx'
      if(nadj.eq.1)goto 50
!     write(6,*)'(greenf): isolated boundary conditions'
! zero adjacent bunches (isolated boundary conditions):
      do k=1,nz+1
      z=(k-1)*hz
      do j=1,ny+1
      y=(j-1)*hy
      do i=1,nx+1
      x=(i-1)*hx
      g(i,j,k)=fourpinv*gfun(x,y,z)
      enddo
      enddo
      enddo
      g(1,1,1)=0.d0
!
      do k=1,nz
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j,k)=-g(n1-i+2,j,k)
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
      g(i,j,k)=-g(n1-i+2,j,n3-k+2)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j,k)=-g(n1-i+2,n2-j+2,k)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j,k)=-g(n1-i+2,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      goto 200
! adjacent bunches (periodic boundary conditions):
   50 continue
      if(idproc.eq.0)write(6,*)'(greenfx)THIS PORTION OF CODE IS BROKEN'
!      nbunch=200
      nbunch=20
      d=hz*nz
      do 100 k=-nz/2,nz/2-1
      do 99 j=-ny,ny-1
      do 98 i=-nx,nx-1
      if(i.eq.0 .and. j.eq.0 .and. k.eq.0)goto 98
      tmp=0.d0
      do n=-nbunch,nbunch
      tmp=tmp+(hx*i)/(fourpi*sqrt((hx*i)**2+(hy*j)**2+(hz*k+n*d)**2))**3
      enddo
      g(1+mod(i+n1,n1),1+mod(j+n2,n2),1+mod(k+nz,nz))=tmp
   98 continue
   99 continue
  100 continue
!cryne July 7, 2004      g(1,1,1)=g(1,1,2)
      g(1,1,1)=0.d0
  200 continue
!     write(6,*)'leaving greenfx'
      return
      end

c--------------------------------------------------------------
c--------------------------------------------------------------

      subroutine greenfy(g,nx,ny,nz,n1,n2,n3,n3a,nadj)
! green function routine.
      implicit double precision(a-h,o-z)
      complex*16 g
      dimension g(n1,n2,n3a)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      gfun(x,y,z)=y/sqrt(x**2+y**2+z**2+1.d-20)**3
!
!*******
! 3/21/2003 include correct factor of 1/(4pi) in the Green function:
      fourpi=8.d0*(asin(1.d0))
      fourpinv=1.d0/fourpi
!*******
!     write(6,*)'inside greenfy'
      if(nadj.eq.1)goto 50
!     write(6,*)'(greenf): isolated boundary conditions'
! zero adjacent bunches (isolated boundary conditions):
      do k=1,nz+1
      z=(k-1)*hz
      do j=1,ny+1
      y=(j-1)*hy
      do i=1,nx+1
      x=(i-1)*hx
      g(i,j,k)=fourpinv*gfun(x,y,z)
      enddo
      enddo
      enddo
      g(1,1,1)=0.d0
!case2         g(1,1,1)=2.d0*g(1,1,2)
!case3         g(1,1,1)=g(2,1,1)
!
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
      g(i,j,k)=-g(i,n2-j+2,k)
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
      g(i,j,k)=-g(i,n2-j+2,n3-k+2)
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
      g(i,j,k)=-g(n1-i+2,n2-j+2,k)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j,k)=-g(n1-i+2,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      goto 200
! adjacent bunches (periodic boundary conditions):
   50 continue
      if(idproc.eq.0)write(6,*)'(greenfy)THIS PORTION OF CODE IS BROKEN'
!      nbunch=200
      nbunch=20
      d=hz*nz
      do 100 k=-nz/2,nz/2-1
      do 99 j=-ny,ny-1
      do 98 i=-nx,nx-1
      if(i.eq.0 .and. j.eq.0 .and. k.eq.0)goto 98
      tmp=0.d0
      do n=-nbunch,nbunch
      tmp=tmp+(hy*j)/(fourpi*sqrt((hx*i)**2+(hy*j)**2+(hz*k+n*d)**2))**3
      enddo
      g(1+mod(i+n1,n1),1+mod(j+n2,n2),1+mod(k+nz,nz))=tmp
   98 continue
   99 continue
  100 continue
!cryne July 7, 2004      g(1,1,1)=g(1,1,2)
      g(1,1,1)=0.d0
  200 continue
!     write(6,*)'leaving greenfy'
      return
      end

c--------------------------------------------------------------
c--------------------------------------------------------------

      subroutine greenfz(g,nx,ny,nz,n1,n2,n3,n3a,nadj)
      use parallel, only : idproc
! green function routine.
      implicit double precision(a-h,o-z)
      complex*16 g
      dimension g(n1,n2,n3a)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      gfun(x,y,z)=z/sqrt(x**2+y**2+z**2+1.d-20)**3
!
!*******
! 3/21/2003 include correct factor of 1/(4pi) in the Green function:
      fourpi=8.d0*(asin(1.d0))
      fourpinv=1.d0/fourpi
!*******
      if(idproc.eq.0)write(6,*)'inside greenfz'
      if(nadj.eq.1)goto 50
!     write(6,*)'(greenf): isolated boundary conditions'
! zero adjacent bunches (isolated boundary conditions):
      do k=1,nz+1
      z=(k-1)*hz
      do j=1,ny+1
      y=(j-1)*hy
      do i=1,nx+1
      x=(i-1)*hx
      g(i,j,k)=fourpinv*gfun(x,y,z)
      enddo
      enddo
      enddo
      g(1,1,1)=0.d0
!case2         g(1,1,1)=2.d0*g(1,1,2)
!case3         g(1,1,1)=g(2,1,1)
!
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
      g(i,j,k)=-g(i,j,n3-k+2)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j,k)=-g(i,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      do k=1+nz,nz+nz
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j,k)=-g(n1-i+2,j,n3-k+2)
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
      g(i,j,k)=-g(n1-i+2,n2-j+2,n3-k+2)
      enddo
      enddo
      enddo
      goto 200
! adjacent bunches (periodic boundary conditions):
   50 continue
      if(idproc.eq.0)write(6,*)'adjacent bunches'
!      nbunch=200
      if(idproc.eq.0)write(6,*)'(greenfz)THIS PORTION OF CODE IS BROKEN'
      nbunch=20
      d=hz*nz
      do 100 k=-nz/2,nz/2-1
      do 99 j=-ny,ny-1
      do 98 i=-nx,nx-1
      if(i.eq.0 .and. j.eq.0 .and. k.eq.0)goto 98
      tmp=0.d0
      do n=-nbunch,nbunch
      tmp=tmp+(hz*k+n*d)/                                                &
     &    (fourpi*sqrt( (hx*i)**2+(hy*j)**2+(hz*k+n*d)**2))**3
      enddo
      g(1+mod(i+n1,n1),1+mod(j+n2,n2),1+mod(k+nz,nz))=tmp
   98 continue
   99 continue
  100 continue
!cryne July 7, 2004      g(1,1,1)=g(1,1,2)
      g(1,1,1)=0.d0
  200 continue
      if(idproc.eq.0)write(6,*)'leaving greenfz'
      return
      end

