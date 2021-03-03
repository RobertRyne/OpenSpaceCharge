!    wakefld.f
!
!    Contains subroutines for computing wakefield forces
!    written by Roman Samulyak (BNL), July 2002 
!               rosamu@bnl.gov, (631) 344-3304
!
!***********************************************************************
      module wakefld_params
      double precision :: z0 = 377.0  ! Free space impedance

      character*8 ::      wktype     !Wakefield type

! Global defaults
      integer ::          gnmodes   !Number of modes
      integer ::          nturns    !Number of turns for collecting
                                    !wakefield forces (long range wakes)
      double precision :: gcond     !Conductivity
      double precision :: gradius   !Average radius of chamber elements
      

! elements creating wake fields

      type WKFLD_WALL
        double precision :: length
        double precision :: conduct
        double precision :: radius
        integer ::          nmodes
      end type WKFLD_WALL

      type WKFLD_RLC
        double precision, dimension(6) :: res_fr
        double precision, dimension(6) :: r_shunt
        double precision, dimension(6) :: quality
        integer ::                        nmodes
      end type WKFLD_RLC

      type WKFLD_PILLBOX
        double precision, dimension(3) :: length
        double precision, dimension(3) :: depth
        double precision, dimension(3) :: radius
        integer ::                        nmodes
      end type WKFLD_PILLBOX

      type WKFLD_DATAFILE
        character*31 ::                   fname
        integer ::                        nmodes
      end type WKFLD_DATAFILE

      type(WKFLD_WALL), dimension(100) :: rwall_param 
      type(WKFLD_RLC), dimension(100) :: rlc_param   
      type(WKFLD_PILLBOX), dimension(100) :: pbx_param  
      type(WKFLD_DATAFILE), dimension(100) :: file_param  
 
      integer :: nrwall = 0 ! Order number of rwall element model
      integer :: nrlc = 0   ! Order number of rlc element model
      integer :: npbx = 0   ! Order number of pbx element model
      integer :: nfile = 0  ! Order number of datafile element

      integer :: ncurel = 0 ! Order number of the current wakefield   
!                             element in the ring
      integer :: nwkel = 0  ! total number of wakefield el. in the list
      character*8, dimension(500) :: wk_el_list ! contains element names 
      integer, dimension(500,3) :: wkfld_work 
! wkfld_work(i,j) :
! i = element # in the list of elements creating wakes
! wkfld_work(i,1) : model # for ith elemnt: 
!                                 0 = read wake function data from file
!                                 1 = res. wall
!                                 2 = rlc resonator
!                                 3 = pbx model
!                                 4 = 
! wkfld_work(i,2) : element # in this group of elements
! wkfld_work(i,3) : # of modes for this el.


! variables used for the beam moments calculation 
      double precision, dimension(:,:), allocatable :: moments
      double precision, dimension(:,:,:), allocatable :: density
      type LOCAL_COORD
        double precision :: x
        integer          :: k
        integer          :: nn ! entry # in the moments array
      end type LOCAL_COORD
      type(LOCAL_COORD) :: lcoord

      type WKGRID
        integer          :: nx,ny,nz  ! the number of nodes.
! users should enter number of mesh cells (2^n numbers)
! number of nodes = number of mesh cells + 1 
        double precision :: hx,hy,hz
        double precision :: xmin,xmax,ymin,ymax,zmin,zmax
        integer :: ndx,ndy ! number of subdomains in x and y directions
      end type WKGRID
      type(WKGRID) :: grid 

! variables used for wake field force calculation 
      double precision, dimension(:,:,:), allocatable :: wkfx,wkfy,wkfz 
      double precision,dimension(:,:,:),allocatable::lwkfx,lwkfy,lwkfz 

! an array for tabulated wake function values 
      double precision, dimension(:,:), allocatable :: wkfunc_table

! temp arrays needed to deposit particles on the grid and interp. fields
      double precision, dimension(:,:), allocatable :: p_data
      integer, dimension(:), allocatable :: indx,jndx,kndx
      integer, dimension(:), allocatable :: indxp1,jndxp1,kndxp1
      double precision, dimension(:), allocatable :: ab,de,gh
      logical, dimension(:), allocatable :: msk
      integer :: flg

      save

      end module wakefld_params

!***********************************************************************

      subroutine wake_defaults(line)

      use wakefld_params, only: wktype,gnmodes,nturns,gradius,gcond,grid
      use parallel, only : nvp,idproc   !cryne

      character*80 :: line
      character*16 :: cbuf
      logical keypres, numpres
      double precision, dimension(1) :: bufr
      integer :: mmax
      mmax=80

      call getparm(line,mmax,'type=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
         wktype=cbuf
      endif

      call getparm(line,mmax,'nmodes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         gnmodes=int(bufr(1))
      endif

      call getparm(line,mmax,'nturns=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         nturns=int(bufr(1))
      endif

      call getparm(line,mmax,'r=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         gradius=bufr(1)
      endif

      call getparm(line,mmax,'conduct=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         gcond=bufr(1)
      endif
! another name
      call getparm(line,mmax,'cond=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         gcond=bufr(1)
      endif

! new line containinf wakefield grid parameters
      call readin(line,leof)

      call getparm(line,mmax,'nx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         grid%nx=int(bufr(1))+1
      endif

      call getparm(line,mmax,'ny=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         grid%ny=int(bufr(1))+1
      endif

      call getparm(line,mmax,'nz=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         grid%nz=int(bufr(1))+1
      endif

      call getparm(line,mmax,'ndx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         grid%ndx=int(bufr(1))
      endif

      call getparm(line,mmax,'ndy=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         grid%ndy=int(bufr(1))
      endif
!cryne---
      if(grid%ndx*grid%ndy.ne.nvp)then
        if(idproc.eq.0)then
          write(6,*)'ERROR (in wake_defaults) : ndx*ndy .ne. nvp'
          write(6,*)'ndx, ndy, nvp = ',grid%ndx,grid%ndy,nvp
          write(6,*)'Halting.'
        endif
        call myexit
      endif
!cryne---

      return
      end subroutine wake_defaults

!***********************************************************************

      subroutine wake_init(cbufin)

      use wakefld_params
      use acceldata
      include 'files.inc'

      character*16 :: cbufin, cbuf
      integer :: mmax, i
      double precision, dimension(1) :: bufr
      logical keypres, numpres, leof
      character*80 :: line
      double precision :: l,r,sigma
      mmax=80
      nwkel=nwkel+1

      if ((cbufin.eq.'rwall').or.(cbufin.eq.'res_wall')) then 
         go to 110
      elseif (cbufin.eq.'rlc') then 
         go to 120
      elseif ((cbufin.eq.'pillbox').or.(cbufin.eq.'pbx')) then
         go to 130
      elseif ((cbufin.eq.'datafile').or.(cbufin.eq.'file')) then 
         go to 140
      else
         write(6,*)'Wakefield model',cbufin,'is not supported'
      call myexit
      endif              

!***  res_wall
! Input format:
! wakedata: nmodes=  l=  r=  conduct=

 110  nrwall = nrwall+1       
      backspace lf
! Defaults:
      call readin(line,leof)
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         rwall_param(nrwall)%length=bufr(1)
      endif
      rwall_param(nrwall)%nmodes = gnmodes
      rwall_param(nrwall)%radius = gradius         
      rwall_param(nrwall)%conduct = gcond

      wk_el_list(nwkel) = lmnlbl(na)
      wkfld_work(nwkel,1) = 1
      wkfld_work(nwkel,2) = nrwall
      wkfld_work(nwkel,3) = gnmodes

! End defaults

      call readin(line,leof)
      if(line(1:9) .ne. 'wakedata:')then
         write(6,*)'No user specified wakefield data. Using defaults'  
         backspace lf
         return
      endif

      call getparm(line,mmax,'nmodes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        rwall_param(nrwall)%nmodes=int(bufr(1))
      endif

      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        rwall_param(nrwall)%length=bufr(1)
      endif

      call getparm(line,mmax,'r=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         rwall_param(nrwall)%radius=bufr(1)
      endif

      call getparm(line,mmax,'conduct=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         rwall_param(nrwall)%conduct=bufr(1)
      endif

      wkfld_work(nwkel,3) = rwall_param(nrwall)%nmodes

      if (ibrief.eq.2) then
         write(6,*)'\n Initialized parameters for rwall:'
         write(6,*)'nrwall = ',nrwall
         write(6,*)'nmodes = ',rwall_param(nrwall)%nmodes
         write(6,*)'l = ',rwall_param(nrwall)%length
         write(6,*)'r = ',rwall_param(nrwall)%radius
         write(6,*)'conduct = ',rwall_param(nrwall)%conduct,'\n'
      endif
      return

!***  rlc
! Input format:
! wakedata: nmodes= 
! mode0: fr=  q=  r=
! mode1: fr=  q=  r=
! ............
! the word 'mode' is optional

 120  nrlc = nrlc + 1
      call readin(line,leof)
      if(line(1:9) .ne. 'wakedata:')then
         write(6,*)'\n Error:'
         write(6,*)'There are no defaults for rlc wakefield model.'
         write(6,*)'Please enter wakefield model data: '
         write(6,*)'the number of modes and the resonant frequency, '
         write(6,*)'shunt impedance, and quality factor for each mode' 
         write(6,*)'in the following format' 
         write(6,*)'wakedata: nmodes= ' 
         write(6,*)'mode0: fr=  q=  r='
         write(6,*)'mode1: fr=  q=  r='
         write(6,*)' ............\n'
         call myexit
         return
      endif

! If there is wake data, get the number of modes
      call getparm(line,mmax,'nmodes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         rlc_param(nrlc)%nmodes=int(bufr(1))
      endif

! Get parameters for each mode
      do i=1,rlc_param(nrlc)%nmodes
         call readin(line,leof)

         call getparm(line,mmax,'fr=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            rlc_param(nrlc)%res_fr(i)=bufr(1)
         endif
         
         call getparm(line,mmax,'r=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            rlc_param(nrlc)%r_shunt(i)=bufr(1)
         endif
         
         call getparm(line,mmax,'q=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            rlc_param(nrlc)%quality(i)=bufr(1)
         endif
      enddo

      wk_el_list(nwkel) = lmnlbl(na)
      wkfld_work(nwkel,1) = 2
      wkfld_work(nwkel,2) = nrlc
      wkfld_work(nwkel,3) = rlc_param(nrlc)%nmodes

      if (ibrief.eq.2) then
      write(6,*)'\n Initialized parameters for rlc:'
      write(6,*)'nrlc = ',nrlc, 'nmodes = ',rlc_param(nrlc)%nmodes
      write(6,*)'fr = ',rlc_param(nrlc)%res_fr(1:rlc_param(nrlc)%nmodes)
      write(6,*)'R = ',rlc_param(nrlc)%r_shunt(1:rlc_param(nrlc)%nmodes)
      write(6,*)'Q = ',rlc_param(nrlc)%quality(1:rlc_param(nrlc)%nmodes)
      write(6,*)'\n'
      endif
      return

!*** pillbox
! Input format:
! wakedata: nmodes= 
! mode0: l=  h=  r=
! mode1: l=  h=  r=
! ............
! the word 'mode' is optional

 130  npbx = npbx + 1 
      call readin(line,leof)
      if(line(1:9) .ne. 'wakedata:')then
         write(6,*)'\n Error:'
         write(6,*)'There are no defaults for pillbox wakefield model.'
         write(6,*)'Please enter wakefield model data: '
         write(6,*)'the number of modes and the length, depth,'
         write(6,*)'and radius for each mode in the following format' 
         write(6,*)'wakedata: nmodes= ' 
         write(6,*)'mode0: l=  h=  r= '
         write(6,*)'mode1: l=  h=  r= '
         write(6,*)' ............\n'
         call myexit
         return
      endif

! If there is wake data, get the number of modes
      call getparm(line,mmax,'nmodes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         pbx_param(npbx)%nmodes=int(bufr(1))
      endif

! Get parameters for each mode
      do i=1,pbx_param(npbx)%nmodes
         call readin(line,leof)
         call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            pbx_param(npbx)%length(i)=bufr(1)
         endif
         
         call getparm(line,mmax,'r=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            pbx_param(npbx)%radius(i)=bufr(1)
         endif
         
         call getparm(line,mmax,'h=',bufr,keypres,numpres,0,cbuf)
         if(keypres.and.numpres)then
            pbx_param(npbx)%depth(i)=bufr(1)
         endif
      enddo

      wk_el_list(nwkel) = lmnlbl(na)
      wkfld_work(nwkel,1) = 3
      wkfld_work(nwkel,2) = npbx
      wkfld_work(nwkel,3) = pbx_param(npbx)%nmodes

      if (ibrief.eq.2) then
      write(6,*)'\n Initialized parameters for pillbox:'
      write(6,*)'npbx = ',npbx, 'nmodes = ',pbx_param(npbx)%nmodes
      write(6,*)'l = ',pbx_param(npbx)%length(1:pbx_param(npbx)%nmodes)
      write(6,*)'r = ',pbx_param(npbx)%radius(1:pbx_param(npbx)%nmodes)
      write(6,*)'h = ',pbx_param(npbx)%depth(1:pbx_param(npbx)%nmodes)
      write(6,*)'\n'
      endif

      return

!*** datafile
! Input format:
! wakedata: file=file_name nmodes=number_of_modes

 140  nfile = nfile + 1
      call readin(line,leof)
      if(line(1:9) .ne. 'wakedata:')then
         write(6,*)'\n Error:'
         write(6,*)'There are no defaults for tabulated wakefield data.'
         write(6,*)'Please enter the file name and the number of modes'
         write(6,*)'in the following format' 
         write(6,*)'wakedata: file=file_name nmodes=number_of_modes\n' 
         call myexit
         return
      endif
! If there is wake data, get parameters
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
         file_param(nfile)%fname=cbuf
      endif

      call getparm(line,mmax,'nmodes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
         file_param(nfile)%nmodes=int(bufr(1))
      endif

      wk_el_list(nwkel) = lmnlbl(na)
      wkfld_work(nwkel,1) = 0
      wkfld_work(nwkel,2) = nfile
      wkfld_work(nwkel,3) = file_param(nfile)%nmodes

      if (ibrief.eq.2) then
      write(6,*)'\nInitialized parameters for tabulated wakefield data:'
      write(6,*)'nfile = ',nfile
      write(6,*)'file name = ',file_param(nfile)%fname
      write(6,*)'nmodes = ',file_param(nfile)%nmodes,'\n'
      endif

      return

      end subroutine wake_init


!***********************************************************************

      subroutine wkfld_srange(nslices,elname,tau)

      use wakefld_params, only: wktype,ncurel,nwkel,wk_el_list

      integer :: nslices, i
      character*8 :: elname
      double precision :: tau

      if (wktype.ne.'short') return

! Identify the current element
      do ncurel=1,nwkel
         if (elname .eq. wk_el_list(ncurel)) exit
      enddo
      
      call initialize_wk_structures
      call compute_beam_moments
      call compute_wk_forces
      call apply_wk_forces(tau)
      call free_wk_structures

      end subroutine wkfld_srange

!***********************************************************************

      subroutine initialize_wk_structures

      use wakefld_params
      use rays, only : nraysp,maxrayp

      implicit double precision(a-h,o-z)
      integer :: nx,ny,nz,nmodes,j
!      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0,nspchset

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      allocate(density(nx,ny,nz))
      allocate(p_data(6,nraysp))

      nmodes = wkfld_work(ncurel,3)
      allocate(moments(nz,2*nmodes+1))

      allocate(wkfx(nx,ny,nz))
      allocate(wkfy(nx,ny,nz))
      allocate(wkfz(nx,ny,nz))

      allocate(indx(nraysp)) 
      allocate(jndx(nraysp)) 
      allocate(kndx(nraysp)) 
      allocate(indxp1(nraysp)) 
      allocate(jndxp1(nraysp)) 
      allocate(kndxp1(nraysp)) 
      allocate(ab(nraysp)) 
      allocate(de(nraysp)) 
      allocate(gh(nraysp)) 
      allocate(msk(maxrayp)) 
 
      end subroutine initialize_wk_structures

!***********************************************************************

      subroutine free_wk_structures

      use wakefld_params

      deallocate(density)
      deallocate(moments)
      deallocate(p_data)
      deallocate(wkfx)
      deallocate(wkfy)
      deallocate(wkfz)

      deallocate(indx) 
      deallocate(jndx) 
      deallocate(kndx) 
      deallocate(indxp1) 
      deallocate(jndxp1) 
      deallocate(kndxp1) 
      deallocate(ab) 
      deallocate(de) 
      deallocate(gh) 
      deallocate(msk) 
      
      end subroutine free_wk_structures

!***********************************************************************

      subroutine wkfld_lrange

      use wakefld_params

      if (wktype.ne.'long') return

      end subroutine wkfld_lrange

!***********************************************************************
! INTEGRATOR
     
      subroutine compute_beam_moments

      use wakefld_params
      use rays, only : zblock,nraysp,maxrayp
      use beamdata
      use parallel

      double precision,dimension(:,:,:), allocatable :: dens_local    
      double precision, dimension(:,:), allocatable :: mom_local
      integer :: i,j,nx,ny,nz,nmodes,nz_start,nz_end,n_zone,noresize
      double precision :: ss,bbyk,gblam
      double precision :: xmin,xmax,ymin,ymax,zmin,zmax
      double precision :: hx,hy,hz,hxi,hyi,hzi
      double precision :: toopi,c5j,ac5j

      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0,nspchset

      external ysimpsonint

      nmodes = wkfld_work(ncurel,3)
      nx = grid%nx
      ny = grid%ny
      nz = grid%nz
 
      allocate(dens_local(nx,ny,nz))
      allocate(mom_local(nz,2*nmodes+1))

      bbyk = gamma*beta*c/omegascl
      gblam = gamma*beta*c/bfreq   

      noresize = 1

      do j = 1,nraysp
         p_data(1,j) = zblock(1,j)
         p_data(2,j) = zblock(2,j)
         p_data(3,j) = zblock(3,j)
         p_data(4,j) = zblock(4,j)
         p_data(5,j) = zblock(5,j)
         p_data(6,j) = zblock(6,j)
      enddo

!      if(nadj0.ne.0)then
!         toopi=4.d0*asin(1.d0)
!         do j=1,nraysp
!            c5j=p_data(5,j)
!            ac5j=abs(c5j)
!            if(c5j.gt.0.)p_data(5,j)=mod(c5j,toopi)
!            if(c5j.lt.0.)p_data(5,j)=toopi-mod(ac5j,toopi)
!         enddo!
!
!c  place the value which is in(-pi,pi) back in the zblock array:
!         do j=1,nraysp
!            zblock(5,j)=p_data(5,j)-0.5d0*toopi
!         enddo
!      endif

c convert phase to z (bbyk) and transform to the beam frame (gamma)
c also multiply by scale length (sl) to convert x and y to units of meters
      do j=1,nraysp
         p_data(5,j)=p_data(5,j)*gamma*bbyk
         p_data(1,j)=p_data(1,j)*sl
         p_data(3,j)=p_data(3,j)*sl
      enddo
      
c mask off lost particles
      do j=1,nraysp
         msk(j)=.true.
      enddo
      if(nraysp.lt.maxrayp)then
         do j=nraysp+1,maxrayp
            msk(j)=.false.
         enddo
      endif

      call setbound(p_data,msk,nraysp,gblam,gamma,nx,ny,nz,noresize,0)

      grid%xmin = xmin
      grid%xmax = xmax
      grid%ymin = ymin
      grid%ymax = ymax
      grid%zmin = zmin
      grid%zmax = zmax
      grid%hx = (xmax - xmin)/(nx - 1)
      grid%hy = (ymax - ymin)/(ny - 1)
      grid%hz = (zmax - zmin)/(nz - 1)
 
! Deposit local chunk of particles on every processor on the entire grid
      call rhoslo3d_loc(p_data,dens_local,msk,nraysp,ab,de,gh,indx,jndx,&
     &kndx,indxp1,jndxp1,kndxp1,nx,ny,nz,0) 
      call MPI_ALLREDUCE(dens_local,density,nx*ny*nz,mreal,mpisum,      &
     &lworld,ierror)

      n_zone = int(nz/nvp)
      nz_start = 1 + n_zone*idproc
      nz_end = n_zone*(idproc + 1)
      if (idproc .eq. nvp - 1) nz_end = nz 

      do i = nz_start,nz_end
        lcoord%k = i
          do j = 1, 2*nmodes+1
            lcoord%nn = j
            call qsimpsonx(ysimpsonint,xmin,xmax,ss)
            mom_local(i,j) = ss
          enddo
      enddo

      call MPI_ALLREDUCE(mom_local,moments,nz*(2*nmodes+1),mreal,mpisum,&
     &lworld,ierror)

      deallocate(dens_local)
      deallocate(mom_local)

      return
      end subroutine compute_beam_moments

!***********************************************************************
      
      function ysimpsonint(xx)

      use wakefld_params
      implicit double precision(a-h,o-z)

      external func_rho
      ymin = grid%ymin
      ymax = grid%ymax

      lcoord%x = xx
      call qsimpsony(func_rho,ymin,ymax,sst)
      ysimpsonint = sst

      return
      end function ysimpsonint

!***********************************************************************

      function func_rho(yy)

      use wakefld_params, only: grid, lcoord, density
      implicit double precision(a-h,o-z)
      integer i,j,k

      xmin = grid%xmin
      xmax = grid%xmax
      ymin = grid%ymin
      ymax = grid%ymax
      hx = grid%hx
      hy = grid%hy

      xx = lcoord%x
      i = int((xx - xmin)/hx) + 1
      j = int((yy - ymin)/hy) + 1
      k = lcoord%k

! mode # 0
      if (lcoord%nn .eq. 1) then
        func_rho = density(i,j,k)
      endif
! mode # 1
      if (lcoord%nn .eq. 2) then
        func_rho = density(i,j,k)*xx
      endif

      if (lcoord%nn .eq. 3) then
        func_rho = density(i,j,k)*yy
      endif
! mode # 2
      if (lcoord%nn .eq. 4) then
        func_rho = density(i,j,k)*(xx**2 - yy**2)
      endif

      if (lcoord%nn .eq. 5) then
        func_rho = density(i,j,k)*2*xx*yy
      endif
! mode # 3
      if (lcoord%nn .eq. 6) then
        func_rho = density(i,j,k)*(xx**3 - 3*xx*yy**2)
      endif

      if (lcoord%nn .eq. 7) then
        func_rho = density(i,j,k)*(3*xx**2*yy - yy**3)
      endif
! mode # 4
      if (lcoord%nn .eq. 8) then
        func_rho = density(i,j,k)
      endif

      if (lcoord%nn .eq. 9) then
        func_rho = density(i,j,k)
      endif
! mode # 5
      if (lcoord%nn .eq. 10) then
        func_rho = density(i,j,k)
      endif

      if (lcoord%nn .eq. 11) then
        func_rho = density(i,j,k)
      endif

      
      if (lcoord%nn .gt. 11) then
        write(6,*)'Error in the beam moment computation:'
        write(6,*)'only 5 modes are supported'
        call myexit
      endif

      return
      end function func_rho

!***********************************************************************

      subroutine qsimpsonx(function,start,end,dint)

      implicit double precision(a-h,o-z)
      integer :: j,jmax
      external function

      eps=1.e-9
      jmax = 10
      ost=-1.e30
      os =-1.e30

      do j=1,jmax
       call trapzdx(function,start,end,st,j)
       dint=(4.*st-ost)/3.
       if (abs(dint-os) .lt. eps*abs(os)) return
       if (dint .eq. 0. .and. os .eq. 0. .and. j .gt. 6) return
       os = dint
       ost = st
      enddo
      end subroutine qsimpsonx

!***********************************************************************
 
      subroutine qsimpsony(function,start,end,dint)

      implicit double precision(a-h,o-z)
      integer :: j,jmax
      external function

      eps=1.e-9
      jmax = 10
      ost=-1.e30
      os =-1.e30

      do j=1,jmax
       call trapzdy(function,start,end,st,j)
       dint=(4.*st-ost)/3.
       if (abs(dint-os) .lt. eps*abs(os)) return
       if (dint .eq. 0. .and. os .eq. 0. .and. j .gt. 6) return
       os = dint
       ost = st
      enddo
      end subroutine qsimpsony

!***********************************************************************

      subroutine trapzdx(func,a,b,s,n)

      implicit double precision(a-h,o-z)
      integer :: n,it,j
      external func

      if (n .eq. 1) then
        s = 0.5*(b-a)*(func(a)+func(b))
      else
        it = 2**(n-2)
        tnm = it
        del = (b-a)/tnm
        x = a+0.5*del
        sum = 0.
        do j = 1,it
          sum = sum + func(x)
          x = x +del
        enddo
        s = 0.5*(s+(b-a)*sum/tnm)
      endif
      return
      end subroutine trapzdx

!***********************************************************************

      subroutine trapzdy(func,a,b,s,n)

      implicit double precision(a-h,o-z)
      integer :: n,it,j
      external func

      if (n .eq. 1) then
        s = 0.5*(b-a)*(func(a)+func(b))
      else
        it = 2**(n-2)
        tnm = it
        del = (b-a)/tnm
        x = a+0.5*del
        sum = 0.
        do j = 1,it
          sum = sum + func(x)
          x = x +del
        enddo
        s = 0.5*(s+(b-a)*sum/tnm)
      endif
      return
      end subroutine trapzdy

!***********************************************************************
      
      subroutine compute_wk_forces

      use wakefld_params
      use parallel

! 2d domain decomposition in x-y plane. ndx, ndy = number of subdomains
! in x and y directions (declared in wakefld_params)
! nxst,nxend,nyst,nyend = start and end grid nodes in x and y directions
! nnx,nny = number of grids in x and y directions
! ii,jj = indices of subdomains:
!                      
! ---------------------|--------------------
!   idproc = 4         | idproc = 5     
!                      |
!   ii = 1, jj = 3     | ii = 2, jj = 3
! ---------------------|--------------------
!   idproc = 2         | idproc = 3     
!                      |
!   ii = 1, jj = 2     | ii = 2,  jj = 2   
! ---------------------|--------------------
!   idproc = 0         | idproc = 1 
!                      |
!   ii = 1, jj = 1     | ii = 2, jj = 1
! ---------------------|--------------------
!                      

      double precision :: rho, z
      integer :: nxst,nxend,nyst,nyend,ii,jj,nnx,nny
      integer :: i,j,k1,k2,nx,ny,nz,nm

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz
      hz = grid%hz
      ndx = grid%ndx
      ndy = grid%ndy
 
      allocate(lwkfx(nx,ny,nz))
      allocate(lwkfy(nx,ny,nz))
      allocate(lwkfz(nx,ny,nz))

      jj = int(idproc/ndx) + 1
      ii = idproc - (jj - 1)*ndx + 1
      nnx = int(nx/ndx)
      nny = int(ny/ndy)

      nxst = 1 + (ii - 1)*nnx
      nxend = ii*nnx
      if(ii .eq. ndx) nxend = nx 

      nyst = 1 + (jj - 1)*nny
      nyend = jj*nny
      if (jj .eq. ndy) nyend = ny 

      do k1 = 1,nz
        do j = nyst,nyend
          do i = nxst,nxend

! For every grid point include wake fields created by
! particles in the head of the bunch           
            do k2 = k1,nz
              z = (k1 - k2)*hz
              rho = density(i,j,k1)
              nm = wkfld_work(ncurel,3)

              if     (wkfld_work(ncurel,1) .eq. 0) then
                 call compute_wkfrc_from_table
              elseif (wkfld_work(ncurel,1) .eq. 1) then
                 call compute_wkfrc_reswall(rho,z,i,j,k1,k2,nm)
              elseif (wkfld_work(ncurel,1) .eq. 2) then
                 call compute_wkfrc_rlc(rho,z,i,j,k1,k2,nm)
              elseif (wkfld_work(ncurel,1) .eq. 3) then
                 call compute_wkfrc_pbx(rho,z,i,j,k1,k2,nm)
              else
                 write(6,*)'Error in compute_wk_forces:'
                 write(6,*)'unimplemented wakefield model'
                 call myexit
              endif              
            enddo
          
          enddo
        enddo
      enddo
      

      call MPI_ALLREDUCE(lwkfx,wkfx,nx*ny*nz,mreal,mpisum,lworld,ierror)
      call MPI_ALLREDUCE(lwkfy,wkfy,nx*ny*nz,mreal,mpisum,lworld,ierror)
      call MPI_ALLREDUCE(lwkfz,wkfz,nx*ny*nz,mreal,mpisum,lworld,ierror)

      deallocate(lwkfx)
      deallocate(lwkfy)
      deallocate(lwkfz)

      ncurel  = ncurel + 1
      if (ncurel .gt. nwkel) then
        ncurel = 1
      endif
      return
      end subroutine compute_wk_forces

!***********************************************************************


      subroutine apply_wk_forces(tau)

      use wakefld_params, only: wkfx,wkfy,wkfz,ab,de,gh,indx,jndx,kndx, &
     &msk,indxp1,jndxp1,kndxp1,grid
      use rays, only: zblock, nraysp
      use beamdata

      double precision, dimension (:), allocatable :: fx,fy,fz
      double precision :: tau,fourpi,perv,xycon,tcon,ratio3,gbi,fpei
      integer :: j

      allocate(fx(nraysp))
      allocate(fy(nraysp))
      allocate(fz(nraysp))

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

! Interpolate forces to positions of particles      

      call ntrslo3d_loc(wkfx,wkfy,wkfz,fx,fy,fz,msk,ab,de,gh,indx,jndx, &
     &kndx,indxp1,jndxp1,kndxp1,nx,ny,nz,nraysp)

      fpei=c*c*1.d-7
      perv=2.d0*bcurr*fpei/(brho*beta*beta*c*c*gamma*gamma)
      fourpi=8.d0*(asin(1.d0))
      xycon=fourpi*0.5*perv*gamma*beta*(beta*c/bfreq)
      tcon=beta*xycon*gamma**2
      ratio3=omegascl*sl/299792458.d0
      tcon=tcon/ratio3

      gbi=1./(gamma*beta)

      if(lflagmagu)then
         xycon=xycon*gbi
         tcon=tcon*gbi
      endif

      open(77,file='wake_data')

      do j=1,nraysp

         write(77,88)zblock(2,j),tau*xycon*fx(j),zblock(4,j),           &
     &               xycon*fy(j),zblock(6,j),tau*tcon*fz(j)
 88      format(6d12.3)

         zblock(2,j)=zblock(2,j) + tau*xycon*fx(j)
         zblock(4,j)=zblock(4,j) + tau*xycon*fy(j)
         zblock(6,j)=zblock(6,j) + tau*tcon*fz(j)
      enddo
      
      deallocate(fx)
      deallocate(fy)
      deallocate(fz)
      
      end subroutine apply_wk_forces

!***********************************************************************

      subroutine compute_wkfrc_reswall(rho,z,ii,jj,kk1,kk2,nmod)

! Computes wake forces acting on the test particle due to wake wields
! created by the source particle using the resistive wake field model
! Writes forces to arrays lwkfx,lwkfy,lwkfz in wakefld_params.mod
!
! z is the distance between the source and test particles
! ii,jj,kk1 are discrete coordinates of the test particle
! kk2 is the discrete z coord. of the source particle
! nmod is the number of modes

      use wakefld_params
      implicit double precision(a-h,o-z)
      integer ii,jj,kk1,kk2,nmod,nrw 

      x = grid%xmin + grid%hx*(ii-1)
      y = grid%ymin + grid%hy*(jj-1)

      nrw = wkfld_work(ncurel,2)
      call wkfunction_rwall(z,1,1,nrw,wkfl)
      lwkfx(ii,jj,kk1) = 0.
      lwkfy(ii,jj,kk1) = 0.
      lwkfz(ii,jj,kk1) = -moments(kk2,1)*wkfl

      if (nmod .gt. 1) then
         call wkfunction_rwall(z,2,0,nrw,wkft)
         call wkfunction_rwall(z,2,1,nrw,wkfl)

         lwkfx(ii,jj,kk1) = -moments(kk2,2)*wkft

         lwkfy(ii,jj,kk1) = -moments(kk2,3)*wkft

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &                   wkfl*(x*moments(kk2,2) + y*moments(kk2,3)) 
      endif
      
      if (nmod .gt. 2) then
         call wkfunction_rwall(z,3,0,nrw,wkft)
         call wkfunction_rwall(z,3,1,nrw,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   2*wkft*(x*moments(kk2,4) + y*moments(kk2,5))        

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   2*wkft*(-y*moments(kk2,4) + x*moments(kk2,5))       

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**2-y**2)*moments(kk2,4)+2*x*y*moments(kk2,5))               
      endif
      
      if (nmod .gt. 3) then
         call wkfunction_rwall(z,4,0,nrw,wkft)
         call wkfunction_rwall(z,4,1,nrw,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   3*wkft*((x**2-y**2)*moments(kk2,6)+2*x*y*moments(kk2,7))

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   3*wkft*(-2*x*y*moments(kk2,6)+(x**2-y**2)*moments(kk2,7))

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**3-3*x*y**2)*moments(kk2,6) +                           &
     &   (3*x**2*y-y**3)*moments(kk2,7))
      endif

!      if (nmod .gt. 4) then
!         call wkfunction_rwall(z,5,0,nrw,wkft)
!         call wkfunction_rwall(z,5,1,nrw,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1)
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                  
!      endif

!      if (nmod .gt. 5) then
!         call wkfunction_rwall(z,6,0,nrw,wkft)
!         call wkfunction_rwall(z,6,1,nrw,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)      
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) 
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                     
!      endif

      return
      end subroutine compute_wkfrc_reswall

!***********************************************************************

      subroutine compute_wkfrc_rlc(rho,z,ii,jj,kk1,kk2,nmod)

! Computes wake forces acting on the test particle due to wake wields
! created by the source particle using the rlc resonator model
! Writes forces to arrays lwkfx,lwkfy,lwkfz in wakefld_params.mod
!
! z is the distance between the source and test particles
! ii,jj,kk1 are discrete coordinates of the test particle
! kk2 is the discrete z coord. of the source particle
! nmod is the number of modes

      use wakefld_params
      implicit double precision(a-h,o-z)
      integer ii,jj,kk1,kk2,nmod,nrf

      x = grid%xmin + grid%hx*(ii-1)
      y = grid%ymin + grid%hy*(jj-1)

      nrf = wkfld_work(ncurel,2)
      call wkfunction_rlc(z,1,1,nrf,wkfl)
      lwkfx(ii,jj,kk1) = 0.
      lwkfy(ii,jj,kk1) = 0.
      lwkfz(ii,jj,kk1) = -moments(kk2,1)*wkfl

      if (nmod .gt. 1) then
         call wkfunction_rlc(z,2,0,nrf,wkft)
         call wkfunction_rlc(z,2,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = -moments(kk2,2)*wkft

         lwkfy(ii,jj,kk1) = -moments(kk2,3)*wkft

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &                   wkfl*(x*moments(kk2,2) + y*moments(kk2,3)) 
      endif
      
      if (nmod .gt. 2) then
         call wkfunction_rlc(z,3,0,nrf,wkft)
         call wkfunction_rlc(z,3,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   2*wkft*(x*moments(kk2,4) + y*moments(kk2,5))        

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   2*wkft*(-y*moments(kk2,4) + x*moments(kk2,5))       

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**2-y**2)*moments(kk2,4)+2*x*y*moments(kk2,5))               
      endif
      
      if (nmod .gt. 3) then
         call wkfunction_rlc(z,4,0,nrf,wkft)
         call wkfunction_rlc(z,4,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   3*wkft*((x**2-y**2)*moments(kk2,6)+2*x*y*moments(kk2,7))

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   3*wkft*(-2*x*y*moments(kk2,6)+(x**2-y**2)*moments(kk2,7))

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**3-3*x*y**2)*moments(kk2,6) +                           &
     &   (3*x**2*y-y**3)*moments(kk2,7))
      endif

!      if (nmod .gt. 4) then
!         call wkfunction_rlc(z,5,0,nrf,wkft)
!         call wkfunction_rlc(z,5,1,nrf,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1)
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                  
!      endif

!      if (nmod .gt. 5) then
!         call wkfunction_rlc(z,6,0,nrf,wkft)
!         call wkfunction_rlc(z,6,1,nrf,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)      
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) 
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                     
!      endif

      return
      end subroutine compute_wkfrc_rlc

!***********************************************************************

      subroutine compute_wkfrc_pbx(rho,z,ii,jj,kk1,kk2,nmod)

! Computes wake forces acting on the test particle due to wake wields
! created by the source particle using the pbx resonator model
! Writes forces to arrays lwkfx,lwkfy,lwkfz in wakefld_params.mod
!
! z is the distance between the source and test particles
! ii,jj,kk1 are discrete coordinates of the test particle
! kk2 is the discrete z coord. of the source particle
! nmod is the number of modes

      use wakefld_params
      implicit double precision(a-h,o-z)
      integer ii,jj,kk1,kk2,nmod,nrf

      x = grid%xmin + grid%hx*(ii-1)
      y = grid%ymin + grid%hy*(jj-1)

      nrf = wkfld_work(ncurel,2)
      call wkfunction_pbx(z,1,1,nrf,wkfl)
      lwkfx(ii,jj,kk1) = 0.
      lwkfy(ii,jj,kk1) = 0.
      lwkfz(ii,jj,kk1) = -moments(kk2,1)*wkfl

      if (nmod .gt. 1) then
         call wkfunction_pbx(z,2,0,nrf,wkft)
         call wkfunction_pbx(z,2,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = -moments(kk2,2)*wkft

         lwkfy(ii,jj,kk1) = -moments(kk2,3)*wkft

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &                   wkfl*(x*moments(kk2,2) + y*moments(kk2,3)) 
      endif
      
      if (nmod .gt. 2) then
         call wkfunction_pbx(z,3,0,nrf,wkft)
         call wkfunction_pbx(z,3,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   2*wkft*(x*moments(kk2,4) + y*moments(kk2,5))        

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   2*wkft*(-y*moments(kk2,4) + x*moments(kk2,5))       

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**2-y**2)*moments(kk2,4)+2*x*y*moments(kk2,5))               
      endif
      
      if (nmod .gt. 3) then
         call wkfunction_pbx(z,4,0,nrf,wkft)
         call wkfunction_pbx(z,4,1,nrf,wkfl)

         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1) -                            &
     &   3*wkft*((x**2-y**2)*moments(kk2,6)+2*x*y*moments(kk2,7))

         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) -                            &
     &   3*wkft*(-2*x*y*moments(kk2,6)+(x**2-y**2)*moments(kk2,7))

         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1) -                            & 
     &   wkfl*((x**3-3*x*y**2)*moments(kk2,6) +                           &
     &   (3*x**2*y-y**3)*moments(kk2,7))
      endif

!      if (nmod .gt. 4) then
!         call wkfunction_pbx(z,5,0,nrf,wkft)
!         call wkfunction_pbx(z,5,1,nrf,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1)
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                  
!      endif

!      if (nmod .gt. 5) then
!         call wkfunction_pbx(z,6,0,nrf,wkft)
!         call wkfunction_pbx(z,6,1,nrf,wkfl)
!         lwkfx(ii,jj,kk1) = lwkfx(ii,jj,kk1)      
!         lwkfy(ii,jj,kk1) = lwkfy(ii,jj,kk1) 
!         lwkfz(ii,jj,kk1) = lwkfz(ii,jj,kk1)                     
!      endif

      return
      end subroutine compute_wkfrc_pbx

!***********************************************************************

      subroutine compute_wkfrc_from_table

      end subroutine compute_wkfrc_from_table

!***********************************************************************

      subroutine wkfunction_rlc(z,m,ntype,nel,wkfunction)

! z is the longitudinal coordinate
! m is the mode number
! ntype = 0 for the transverse wake function
! ntype = 1 for the longitudinal wake function
! nel is the order number of this rlc element model
! wkfunction is the wake function value 

      use wakefld_params
      use beamdata, only: c

      implicit double precision(a-h,o-z)
      integer m,ntype,nel
      common/pie/pi,pi180,twopi

      Q = rlc_param(nel)%quality(m)
      omega = rlc_param(nel)%res_fr(m)
      Rs = rlc_param(nel)%r_shunt(m)
      alpha = 0.5*omega/Q
      bomega = sqrt(omega**2 - alpha**2)

      if (ntype .eq. 0) then ! transverse wake function
        if (z .gt. 0.) then
           write(6,*)'Error: positive z in subroutine wkfunction_rlc'
           call myexit
        else
         wkfunction=c*Rs*omega*exp(alpha*z/c)*sin(bomega*z/c)/(Q*bomega)
        endif
      endif

      if (ntype .eq. 1) then ! longitudinal wake function
        if (z .gt. 0.) then
           write(6,*)'Error: positive z in subroutine wkfunction_rlc'
           call myexit
        elseif (z .eq. 0.) then
           wkfunction=alpha*Rs
        else
           wkfunction=2*alpha*Rs*exp(alpha*z/c)*(cos(bomega*z/c)+       &
     &                alpha*sin(bomega*z/c)/bomega)
        endif
      endif

      if (ntype .gt. 1) then ! invalid ntype
        write(6,*)'Invalid ntype of the rlc wake function'
        call myexit
      endif
      return

      end subroutine wkfunction_rlc

!***********************************************************************

      subroutine wkfunction_rwall(z,m,ntype,nel,wkfunction)

! z is the distance between the source and test particles (with - sign) 
! m is the mode number
! ntype = 0 for the transverse wake function
! ntype = 1 for the longitudinal wake function
! nel is the order number of this res. wall element model
! wkfunction is the wake function value 

      use wakefld_params
      use beamdata, only : c

      implicit double precision(a-h,o-z)
      integer m,ntype,nel
      common/pie/pi,pi180,twopi

      if (z .eq. 0) then
         wkfunction = 0
         return
      endif

      if (m .eq. 0) then
         delta = 1.
      else
         delta = 0.
      endif

      if (ntype .eq. 0) then ! transverse wake function
        wkfunction = c*sqrt(z0)/(pi*(rwall_param(nel)%radius)**(m+1)*   &
     &               (1+delta)*sqrt(pi*rwall_param(nel)%conduct*abs(z)))
      endif

      if (ntype .eq. 1) then ! longitudinal wake function
        wkfunction = c*sqrt(z0)/(twopi*(rwall_param(nel)%radius)**(m+1)*&
     &             (1+delta)*z*sqrt(pi*rwall_param(nel)%conduct*abs(z)))
      endif

      if (ntype .gt. 1) then ! invalid ntype
        write(6,*)'Invalid ntype of the resistive wall wake function'
        call myexit
      endif
      return
      end subroutine wkfunction_rwall

!***********************************************************************

      subroutine wkfunction_pbx(z,m,ntype,nel,wkfunction)

! z is the longitudinal coordinate
! m is the mode number
! ntype = 0 for the transverse wake function
! ntype = 1 for the longitudinal wake function
! nel is the order number of this rlc element model
! wkfunction is the wake function value 

      use wakefld_params
      use beamdata, only: c

      implicit double precision(a-h,o-z)
      integer m,ntype,nel
      common/pie/pi,pi180,twopi

      l = pbx_param(nel)%length(m)
      d = pbx_param(nel)%depth(m)
      r = pbx_param(nel)%radius(m)

      if (m .eq. 1) then
         dlt = 2
      else
         dlt = 1
      endif

      if (ntype .eq. 0) then ! transverse wake function
        if (z .gt. 0.) then
           write(6,*)'Error: positive z in subroutine wkfunction_pbx'
           call myexit
        else
         wkfunction=2*z0*c*sqrt(-2.d0*l*z)/(dlt*pi**2*r**(2*m-1))
        endif
      endif

      if (ntype .eq. 1) then ! longitudinal wake function
        if (z .gt. 0.) then
           write(6,*)'Error: positive z in subroutine wkfunction_pbx'
           call myexit
        elseif (z .eq. 0.) then
           wkfunction=0
        else
           wkfunction=z0*c*sqrt(2.d0*l)/(dlt*pi**2*r**(2*m-1)*sqrt(-z))
        endif
      endif

      if (ntype .gt. 1) then ! invalid ntype
        write(6,*)'Invalid ntype of the pbx wake function'
        call myexit
      endif
      return

      end subroutine wkfunction_pbx

!***********************************************************************


      subroutine locate_index(x,n,j)
      
!  Given an array wkfunc_table(:,:), and given a value x, returns a value 
!  j  such that  x is between wkfunc_table(j,1) and wkfunc_table(j+1,1). 
!  wkfunc_table(1:n,1) must be monotonic. j = 0 or j = n is rturned to 
!  indicate that x is out of range

      use wakefld_params, only: wkfunc_table

      implicit double precision(a-h,o-z)
      integer j,n
      integer jl,jm,ju
!temp
      allocate(wkfunc_table(100,2))

      jl = 0
      ju = n+1
 10   if (ju-jl .gt. 1) then
         jm = (ju+jl)/2
         if ((wkfunc_table(n,1) .ge. wkfunc_table(1,1)) .eqv.           &
     &       (x .ge. wkfunc_table(jm,1))) then
!XXX <dbs> gfortran compiler complains about uninitialized variable
!XXX            jl = lm
            write(6,*) 'error: locate_index(): uninitialized var: lm'
            stop 'LMERROR'
!XXX
         else
            ju = jm
         endif
      goto 10
      endif
      if (x .eq. wkfunc_table(1,1)) then
         j = 1
      else if (x .eq. wkfunc_table(n,1)) then
         j = n-1
      else
         j = jl
      endif
!temp
      deallocate(wkfunc_table)


      return
      end subroutine locate_index

!***********************************************************************

      subroutine rhoslo3d_loc(coord,rho,msk,np,ab,de,gh,indx,jndx,kndx, &
     &indxp1,jndxp1,kndxp1,nx,ny,nz,nadj)
cryne 08/24/2001      use hpf_library
      implicit double precision(a-h,o-z)
      logical msk
      dimension coord(6,np),msk(np),vol(np)
!hpf$ distribute coord(*,block)
!hpf$ align (:) with coord(*,:) :: msk,vol
      dimension ab(np),de(np),gh(np),indx(np),jndx(np),kndx(np),        &
     &indxp1(np),jndxp1(np),kndxp1(np)
!hpf$ align (:) with coord(*,:) :: ab,de,gh,indx,jndx,kndx
!hpf$ align (:) with coord(*,:) :: indxp1,jndxp1,kndxp1
      dimension rho(nx,ny,nz),tmp(nx,ny,nz)
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: tmp
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

!*************************************************************


      subroutine ntrslo3d_loc(exg,eyg,ezg,ex,ey,ez,msk,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1,nx,ny,nz,np)
      use parallel
      implicit double precision(a-h,o-z)
      logical msk
      dimension exg(nx,ny,nz),eyg(nx,ny,nz),ezg(nx,ny,nz)
!hpf$ distribute exg(*,*,block)
!hpf$ align (*,*,:) with exg(*,*,:) :: eyg,ezg
      dimension ex(np),ey(np),ez(np),msk(np)
      dimension ab(np),de(np),gh(np),indx(np),jndx(np),kndx(np),
     #indxp1(np),jndxp1(np),kndxp1(np)
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
!     if(idproc.eq.0)write(6,*)'leaving ntrslo3d_loc'
      return
      end










