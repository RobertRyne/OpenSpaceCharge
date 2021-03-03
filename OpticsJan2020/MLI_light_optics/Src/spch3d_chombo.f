! IMPACT 3D space charge routines
! Copyright 2001 University of California
!
! Using Chombo AMR Poisson solver for Infinite Domain or Homogenous Dirichlet BCs
!
!     subroutine SPCH3DAMR( c,ex,ey,ez,msk,Np,Ntot
!    &                     ,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,phi )

! Arguments:
!   c      in:out (1:5:2,*) are x,y,z coordinates of particles
!          note: in the rest of MLI column 5 is time-of-flight (i.e. phase),
!          but prior to calling the space charge routines it is
!          converted to longitudinal position z.
!                 (2:6:2,*) are momenta and not used here.
!   ex,ey,ez :out electric field at particles
!   Msk      :in  flag indicating valid particles
!   Np       :in  number of particles (both good and bad) on this processor
!   Ntot     :in  number of (good) particles on all processors
!   Nx,Ny,Nz :in  dimensions of grid on the regular sized grid
!   N1,N2,N3 :in  dimensions of the bigger (usually doubled) grid
!   N3a      :in  =Nz for periodic bc longitudinally, =2*Nz for open bc longit.
!   Nadj     :in  =0 for open or Dirichlet BCs, >=1 for longitudinally periodic BCs
!   phi      :out solution to Poisson equation
!
! Globals:
!  Common /NEWPOISSON/
!    ISolve :in choose which Poisson solver to use -- 1 for default Infinite Domain,
!               3x for Chombo AMR Infinite Domain solver
!               4x for Chombo Hom.Dirichlet
!               x0 = spectral discretization
!               x1 = Mehrstellen "
!               x2 = Laplacian   " (7-point)
!               x3 = Mehrstellen (4th order)
!  Common /SHOWME/
!     IVerbose :in
!   Common /GRIDSZ3D/
!     Xmin,Xmax   :in  location of corners of physical grid
!     Ymin,Ymax   :in
!     Zmin,Zmax   :in
!     Hx,Hy,Hz    :in  grid spacings
!     Hxi,Hyi,Hzi :in  1.0 / Hx,Hy,Hz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SPCH3DAMR( c,ex,ey,ez,Msk,Np,Ntot                      &
     &                     ,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj,phi )
      use parallel ,ONLY : NVP,IDPROC
      use ml_timer
      implicit none
!Constants
      integer ,parameter :: InterpType = 1
      integer            :: Nxyz                  !set below
      real*8             :: FourPi, Hxyz, Neg4Pi  ! "   "
!Arguments
      integer Np,Ntot,Nx,Ny,Nz,N1,N2,N3,N3a,Nadj
      real*8,  dimension(6,Np)     :: c
      logical, dimension(Np)       :: Msk
      real*8,  dimension(Np)       :: ex,ey,ez
      real*8,  dimension(Nx,Ny,Nz) :: phi
!Locals
      character filename*32
      integer i,j,k,ierror,status ,stride
      integer CHPHandle ,maxlvls ,maxp ,tagsbuffer ,solvertype
      integer domain(3,2) ,subdomains(3,2,1) ,bcflags(3,2)
      integer refratios(10) 
      logical ltmp
!!!      real*8, dimension(Nx,Ny,Nz) :: exg,eyg,ezg,rhosum
      real*8       x0(3) ,fillratio,tolerance ,charge ,bcvals(3,2)
      real*8       checknorm ,scale ,xval,yval,zval,del2
      integer      iteration ,plotiters
      save         iteration
      data         iteration/0/
      character*12 stencil(0:3) !input names of stencils, indexed by ISolve
      data stencil/'spectral' ,'mehrstellen' ,'laplace' ,'4mehrstellen'/
      save stencil

!Globals
      real*8          Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Hx,Hy,Hz,Hxi,Hyi,Hzi
      common/GRIDSZ3D/Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,Hx,Hy,Hz,Hxi,Hyi,Hzi
      integer       IVerbose
      common/SHOWME/IVerbose
      integer           IDirectFieldCalc,IDensityFunction,ISolve
      common/NEWPOISSON/IDirectFieldCalc,IDensityFunction,ISolve
      character*256     ChomboFilename
      common/chombochar/ChomboFilename
!Externals
      integer  LENC
      logical        HANDLE_ERROR
      external LENC ,HANDLE_ERROR

!======================================================================

!Execution

      if(idproc.eq.0 .AND. iverbose.ge.5)
     &   write(6,*)'info: SPCH3DAMR: solving for phi with '
     &            ,TRIM( stencil( MOD( ISolve,10 ) ) ) ,' stencil'

      !! Set constants
      FourPi = 4.0 * ACOS( -1.0d0 )
      Neg4Pi = -FourPi
      Nxyz = Nx*Ny*Nz
      Hxyz = Hx*Hy*Hz

!-----------------------------------------------------------------------

! check for valid conditions for this solver

      if( ISolve.LT.30 .OR. ISolve.GT.42 )then
        write(6,*) 'error: SPCH3DAMR: invalid Chombo solver type. '
     &            ,' Should be 30--42, not ',ISolve
        stop 'SPCH3DAMR'
      endif
      if( ISolve - 30 .LT. 10 .AND. ( Nx.NE.Ny .OR. Ny.NE.Nz ) )then
        write(6,*) 'error: SPCH3DAMR: non-cubic grids are not supported'
     &            ,'.  Nx,Ny,Nz = ' ,Nx,Ny,Nz
        stop 'SPCH3DAMR'
      endif
      if( Hx.NE.Hy .OR. Hy.NE.Hz )then
        write(6,*) 'error: SPCH3DAMR: anisotropic grids are not '
     &            ,'supported.  Hx,Hy,Hz = ' ,Hx,Hy,Hz
        stop 'SPCH3DAMR'
      endif
      if( IDirectFieldCalc .EQ. 1 )then
        write(6,*) 'error: SPCH3DAMR: direct Efield calculation is not '
     &            ,'supported (i.e. IDirectFieldCalc should be 0).'
        stop 'SPCH3DAMR'
      endif
      if( Nadj .NE. 0 )then
        write(6,*) 'error: SPCH3DAMR: periodic longitudinal BC is not '
     &            ,'supported (i.e. Nadj should be 0).'
        stop 'SPCH3DAMR'
      endif
      if( NVP .NE. 1 )then
        write(6,*) 'error: SPCH3DAMR: parallel not implemented.'
        stop 'NVP>1'
      endif
      if( Ntot .NE. NP )then
        write(6,*) 'error: SPCH3DAMR: all particles not on process 1'
        stop 'NP!=NPtot'
      endif
      if( COUNT( Msk ) .NE. NP )then
        write(6,*) 'error: SPCH3DAMR: masked particles not supported.'
        stop 'Msk!=NP'
      endif
      if( ChomboFilename .EQ. ' ' )then
        write(6,*) 'error: SPCH3DAMR: chombo input filename is blank.'
        stop 'SPCH3DAMR'
      endif

!-----------------------------------------------------------------------

      !! keep track of how many calls
      iteration = iteration + 1

      !! read the Chombo input file from
      call CH_READINFILE( ChomboFilename ,maxlvls ,refratios ,maxp
     &                   ,tagsbuffer ,fillratio ,tolerance ,solvertype
     &                   ,bcvals ,plotiters )

      !! instantiate a ChomboPIC object and return a handle to it
      call CHP_CREATE( CHPHandle ,status )
      ltmp = HANDLE_ERROR( status, 'CHP_CREATE' ) !ignore warnings

      call CHP_SETDEBUG( CHPHandle, IVerbose )

      ! set the computational domain (in grid points)
      ! and the per-processor subdomains
!XXX for now, assume one cpu
      ! domain is in nodes
      domain(:,1) = 0
      domain(1,2) = Nx-1 ; domain(2,2) = Ny-1 ; domain(3,2) = Nz-1
      subdomains(:,:,1) = domain
      x0(1) = Xmin ; x0(2) = Ymin ; x0(3) = Zmin

      call CHP_SETGRIDPARAMS( CHPHandle, x0,Hx,domain,subdomains,maxlvls
     &                       ,refratios ,maxp ,tagsbuffer ,fillratio
     &                       ,status )
      ltmp = HANDLE_ERROR( status, 'CHP_SETGRID' )
      if( ltmp ) stop 'SPCH3A02'

      if( IVerbose .GE. 5 .AND. iteration .EQ. 1 )then
        write(6,*) 'SPCH3D_CHOMBO parameters:'
        write(6,*) ' MaxLevels',maxlvls
        write(6,*) ' RefRatios',refratios
        write(6,*) ' MaxParticlesPerCell',MaxP
        write(6,*) ' TagsBufferCells',tagsbuffer
        write(6,*) ' FillRatio',fillratio
        write(6,*) ' Tolerance',Tolerance
        write(6,*) ' SolverType',solvertype
        write(6,*) ' BoundaryValues' ,bcvals
        write(6,*) ' PlotIters' ,plotiters
      endif

      !! set solver parameters
      bcflags = 0  !default dirichlet
      if( ISolve .LT. 40 )then
        !! infinite domain boundary condition on all faces
        bcflags = 5  !5==Chombo::PICOpen
!XXX -- discretization type is hardcoded        
!XXX        !! the infinite domain solver has 3 flavors for the type of
!XXX        !! discretization used on the level_0 FFT solver, so extract
!XXX        !! this from the ISolve global paramter and pack it into solvertype
!XXX        !! (default: spectral, +10 for Mehrstellan, +20 for laplace
!XXX        solvertype = solvertype + 10*(ISolve-30)
        if( ISolve-30 .NE. 1 )then
          write(6,*) 'warning: SPCH3DAMR: ignoring FFT discretization '
     &              ,'type ' ,ISolve-30
          write(6,*) 'info: using Mehrstellen discretization instead.'
        endif
!XXX
      endif
      
      call CHP_SETSOLVEPARAMS( CHPHandle ,tolerance ,bcflags,bcvals
     &                      ,InterpType ,IVerbose,solvertype,status )
      ltmp = HANDLE_ERROR( status, 'CHP_SETSOLVE' )
      if( ltmp ) stop 'SPCH3A03'

!-----------------------------------------------------------------------

! pass particles to ChomboPIC

      !! divide charge by number of particles so rho will be same as MLI
      !![NOTE: Chombo::PIC divides by cell volume, so dont need to do that here.]
      !![NOTE: poisson3d solver assumes rho is negative.]
      call increment_timer('rhoslo3d',0)
      charge = -1.0d0 / Ntot 
      stride = 6
      call CHP_PUTPARTICLES( CHPHandle,1,charge,NP,stride
     &                      ,c(1,1),c(3,1),c(5,1)       !x,y,z
     &                      ,status )
      ltmp = HANDLE_ERROR( status, 'CHP_PUTPARTICLES' )
      call increment_timer('rhoslo3d',1)
      if( ltmp ) stop 'SPCH3A04'

!-----------------------------------------------------------------------

! solve the Poisson equation using the ChomboPIC solver

      call increment_timer('fft',0)
      call CHP_SOLVE( CHPHandle ,status )
      ltmp = HANDLE_ERROR( status, 'CHP_SOLVE' )
      call increment_timer('fft',1)
      if( ltmp ) stop 'SPCH3A05'

!-----------------------------------------------------------------------

! extract results from Chombo::PIC

      call increment_timer('timer1',0)

      !! get the solution on the base grid
      call CHP_GETPHIGRID0( CHPHandle ,phi ,status )
      ltmp = HANDLE_ERROR( status, 'CHP_GETPHIGRID' )
      call increment_timer('timer1',1)
      if( ltmp ) stop 'SPCH3A08'

      !! get the electric field at the particles
      call CHP_GETEFIELD( CHPHandle, 1,NP,stride ,c(1,1),c(3,1),c(5,1)
     &                   ,1 ,ex ,ey ,ez ,status )  !stride for exyz is 1, not 6
      ltmp = HANDLE_ERROR( status, 'CHP_GETEFIELD' )
      if( ltmp ) stop 'SPCH3A06'

      ![NOTE: this may need a factor of 'h' as well]
!!      ex = ex * FourPi
!!      ey = ey * FourPi
!!      ez = ez * FourPi
      ex = -ex
      ey = -ey
      ez = -ez

      call increment_timer('timer1',1)

!-----------------------------------------------------------------------

! write data files

      call increment_timer('timer2',0)

      if( iteration .EQ. 1 .OR. MOD( iteration,plotiters ) .EQ. 0 )then

        !! poisson solution on all grids
        filename = 'phi_0000.hdf5 '
        write(filename(5:8),'(I4.4)') iteration
        if( iverbose .ge. 4 )
     &      write(6,*) 'SPCH3DAMR: writing phi to ',TRIM(filename)
        call CHP_WRITEPHI( CHPHandle,filename,status )
        ltmp = HANDLE_ERROR( status, 'CHP_WRITEPHI' )
        ! ignore errors or warnings

        !! all data including particles on all grids
        filename = 'soln_0000.hdf5 '
        write(filename(6:9),'(I4.4)') iteration
        if( iverbose .ge. 4 )
     &      write(6,*) 'SPCH3DAMR: writing solution to ',TRIM(filename)
        call CHP_WRITESOLN( CHPHandle ,filename ,status )
        ltmp = HANDLE_ERROR( status, 'CHP_WRITESOLN' )
        ! ignore errors or warnings

      endif
      
      call increment_timer('timer2',1)

!======================================================================

! Done
      call CHP_DESTROY( CHPHandle, status )
      ltmp = HANDLE_ERROR( status, 'CHP_DESTROY' ) !ignore warnings
      if(idproc.eq.0 .AND. iverbose.gt.5) write(6,*)'leaving SPCH3DAMR'
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function HANDLE_ERROR( Status ,IdString )
      integer Status
      character IdString*(*)
      if( Status .EQ. 0 )then
          HANDLE_ERROR = .FALSE.
      elseif( Status .LT. 0 )then
          write(6,*) 'FATALERROR(' ,Status ,') from '
     $        ,IdString(1:LEN(IdString))
          stop
      else
          write(6,*) 'WARNING(' ,Status ,') from '
     $        ,IdString(1:LEN(IdString))
          HANDLE_ERROR = .TRUE. 
      endif
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CH_READINFILE( FileName ,maxlvls ,refratios ,maxp
     &                        ,tagsbuffer ,fillratio ,tolerance
     &                        ,solvertype ,bcvals ,plotiters )

!! Input file format:  (filename: chombo.input)
!!   Version <int>  -- (V2 or later)
!!   MaxLevels <int> -- number of AMR levels to use
!!   RefRatios <int[]> -- refinement ratios for levels > 0 (none if MaxLevels==1)
!!   MaxParticlesPerCell <int> -- threshold for making tags for MeshRefine
!!   TagsBufferCells <int>  -- number of cells around tags to include in grids (V3 or later)
!!   FillRatio <Real> -- portion of cells in a refined box that must be tagged (V4 or later)
!!   Tolerance <Real> -- AMRNodeElliptic solver tolerance (V4 or later)
!!   SolverType <int> -- AMRNodeElliptic solver type (0=default, 1=alternate)
!!   BoundaryValues <RealVect> -- homogenous BC values for each direction (both faces) (V2 or later)
!!   PlotIters <int> -- number of iterations between plot files (V5 or later) (def: 1)          
!! Notes:
!!
!! Defaults:
!!  MaxLevels: 1
!!  RefRatios: none
!!  MaxParticlesPerCell: 10          
!!  TagsBufferCells: 1
!!  FillRatio: .75
!!  Tolerance: 1e-12
!!  SolverType: 0          
!!  BoundaryValues: 0 0 0
!!  PlotIters: 1
!!
!!=========================================================================
      implicit none
!Arguments
      character FileName*(*)
      integer maxlvls ,refratios(*) ,maxp ,tagsbuffer ,solvertype
      integer plotiters
      real*8 x0(3) ,fillratio ,tolerance ,bcvals(3,2)
!Locals
      character word*32
      integer i ,version 
!Globals
      integer       IVerbose
      common/SHOWME/IVerbose
!Execute
      version = 4
      maxlvls = 1
      maxp = 10
      tagsbuffer = 1
      fillratio = 0.75
      tolerance = 1.0d-12
      solvertype = 0
      bcvals = 0
      plotiters = 1

      open(2,file=FileName,err=99)
      read(2,*,err=99,end=99) word ,version
      if( word .NE. 'Version' )then
          write(6,*) 'error: Chombo input: expecting Version, got '
     &              ,version
          stop 'CHOMBOIN'
      endif
      if( version .LT. 4 )then
        write(6,*) 'error: Chombo input: input file versions <4 '
     &            ,' are not supported.'
        stop 'CHOMBOIN'
      endif
      read(2,*) word, maxlvls
      write(6,*) 'Expecting MaxLevels, got [',TRIM(word),']=',maxlvls
      read(2,*) word, (refratios(i),i=1,maxlvls-1)
      write(6,*) 'Expecting RefRatios, got [',TRIM(word),']='
     &          ,(refratios(i),i=1,maxlvls-1)
      read(2,*) word, maxp
      write(6,*) 'Expecting MaxParticlesPerCell, got [',TRIM(word),']='
     &          ,maxp
      read(2,*) word,tagsbuffer
      write(6,*) 'Expecting TagsBufferCells, got [',TRIM(word),']='
     &          ,tagsbuffer
      read(2,*) word,fillratio
      write(6,*) 'Expecting FillRatio, got [',TRIM(word),']=',fillratio
      read(2,*) word,tolerance
      write(6,*) 'Expecting Tolerance, got [',TRIM(word),']=',tolerance
      read(2,*) word,solvertype
      write(6,*) 'Expecting SolverType, got [',TRIM(word),']='
     &          ,solvertype
      read(2,*) word,(bcvals(i,1),i=1,3)
      do i = 1,3
        bcvals(i,2) = bcvals(i,1)
      enddo
      if( version .GE. 5 )then
        read(2,*) word,plotiters
        write(6,*) 'Expecting PlotIters, got [',TRIM(word),']='
     &            ,plotiters
      endif
!Done
      close(2)
      return

!Errors
      ! handle missing file
   99 write(6,*) 'info: Chombo: input file [' ,TRIM( FileName )
     &          ,'] is missing or empty so using default input values.'
      return
      end
