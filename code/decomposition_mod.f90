module decomposition_mod

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128


type index_struct
 integer :: n(3)  = 0   ! number of points in each direction x, y, z
 integer :: lo(3) = 0   ! Low index
 integer :: hi(3) = 0   ! high index = n + low - 1
end type


type domain_decomposition_struct
  integer :: rank                       ! Rank of this subdomain. 
  integer :: n_process                  ! Maximum rank  
  integer :: communicator               ! MPI communicator
  integer :: idecomp                    ! decomposition style 
  type(index_struct) :: process         ! process indices
  type(index_struct) :: global          ! global grid indices
  type(index_struct) :: local           ! local grid indices
end type


contains

subroutine allocate_real_3d(data, index)
type(index_struct) :: index
real(dp), allocatable :: data(:,:,:)

allocate(data(index%lo(1):index%hi(1), &
              index%lo(2):index%hi(2), &
              index%lo(3):index%hi(3)))

end subroutine

subroutine allocate_complex_3d(data, index)
type(index_struct) :: index
complex(dp), allocatable :: data(:,:,:)

allocate(data(index%lo(1):index%hi(1), &
              index%lo(2):index%hi(2), &
              index%lo(3):index%hi(3)))

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+   
!
!
! Input:
!   domain    -- domain_decomposition_struct 
!            %rank
!            %n_process
!            %global
!			 %idecomp, optional
!          
! Output:
!   domain    -- domain_decomposition_struct 
!                %process (if not previously set)
!                %local

subroutine init_domain_decomposition(domain,  communicator)
type (domain_decomposition_struct) :: domain
integer, optional :: communicator 

if (present(communicator)) domain%communicator = communicator
  
! Check if processes have been decomposed
if (domain%process%n(1) == 0) then
  call procgriddecomp(domain%n_process, &
               domain%global%lo(1), domain%global%hi(1),  &
               domain%global%lo(2), domain%global%hi(2),  &
               domain%global%lo(3), domain%global%hi(3),  & 
               domain%idecomp, &
               domain%process%n(1), domain%process%n(2), domain%process%n(3))
endif 

! Set lo, hi indices
call decompose(domain%rank, domain%n_process, &
               domain%global%lo(1), domain%global%hi(1),  &
               domain%global%lo(2), domain%global%hi(2),  &
               domain%global%lo(3), domain%global%hi(3),  & 
               domain%idecomp, &
               domain%process%n(1), domain%process%n(2), domain%process%n(3), &
               domain%local%lo(1), domain%local%hi(1), &
               domain%local%lo(2), domain%local%hi(2), &
               domain%local%lo(3), domain%local%hi(3) )

! Set for consistency:
domain%process%lo = 1
domain%process%hi = domain%process%n 
domain%local%n = domain%local%hi - domain%local%lo +1
domain%global%n = domain%global%hi - domain%global%lo +1

end subroutine


subroutine print_domain_decomposition(d)
type (domain_decomposition_struct) :: d

print *, 'rank: ', d%rank
print *, 'n_process: ', d%n_process
print *, 'communicator: ', d%communicator
print *, 'idecomp: ', d%idecomp
print *, 'process indices: '
call print_index(d%process)
print *, 'global indices: '
call print_index(d%global)
print *, 'local indices: '
call print_index(d%local)

end subroutine

subroutine print_index(ix)
type(index_struct) :: ix
print *, '  n: ', ix%n
print *, '  lo: ', ix%lo
print *, '  hi: ', ix%hi
end subroutine

! -----------------------------------------------------------------------
! partition array for all possibilities
! nx = global # in fast dimension
! ny = global # in mid dimension
! nz = global # in slow dimension

! idecomp = 0 = xyz (blocks)
!           1 = xy, 2 = yz, 3 = xz
!           4 = x, 5 = y, 6 = z


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+      
! modified Plimpton routine to allow arbitrary lower bounds (nxlo,nylo,nzlo) and upper bounds (nxhi,nyhi,nzhi)
! also, this version allows the user to specify npx,npy,npz on input (in which case idecomp must equal -1)
subroutine decompose(node,mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,   &
                     idecomp,npxarg,npyarg,npzarg,ilo,ihi,jlo,jhi,klo,khi)
implicit none

integer :: node,mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,idecomp,npxarg,npyarg,npzarg,ilo,ihi,jlo,jhi,klo,khi
integer :: nx,ny,nz,npx,npy,npz,node_x,node_y,node_z
!consistency checks:
if(.not.(idecomp.ge.-1.or.idecomp.le.6))then
  write(6,*)'(decompose) idecomp error: idecomp=',idecomp
!       call myexit
  stop
endif
if((idecomp.eq.-1).and.(npxarg*npyarg*npzarg.ne.mprocs))then !this could be commented out if it is checked in the main program
  write(6,*)'error: idecomp.eq.-1 but npx*npy*npz .ne. mprocs'
!       call myexit
  stop
endif
!end of consistency checks


nx=nxhi-nxlo+1
ny=nyhi-nylo+1
nz=nzhi-nzlo+1
if(idecomp.eq.-1)then
  npx=npxarg
  npy=npyarg
  npz=npzarg
else
  npx=1
  npy=1
  npz=1
endif

if(idecomp.eq.0)call proc2grid3d(mprocs,nx,ny,nz,npx,npy,npz) !xyz decomp
if(idecomp.eq.1)call proc2grid2d(mprocs,nx,ny,npx,npy)        !xy
if(idecomp.eq.2)call proc2grid2d(mprocs,ny,nz,npy,npz)        !yz
if(idecomp.eq.3)call proc2grid2d(mprocs,nx,nz,npx,npz)        !xz

if(npx.ne.1.and.npy.ne.1.and.npz.ne.1)then !xyz decomp
  node_x = mod(node,npx)
  node_y = mod(node/npx,npy)
  node_z = node/(npx*npy)
  ilo = node_x*nx/npx + 1
  ihi = (node_x+1)*nx/npx
  jlo = node_y*ny/npy + 1
  jhi = (node_y+1)*ny/npy
  klo = node_z*nz/npz + 1
  khi = (node_z+1)*nz/npz
elseif(npx.ne.1.and.npy.ne.1.and.npz.eq.1)then !xy decomp
  node_x = mod(node,npx)
  node_y = node/npx
  ilo = node_x*nx/npx + 1
  ihi = (node_x+1)*nx/npx
  jlo = node_y*ny/npy + 1
  jhi = (node_y+1)*ny/npy
  klo = 1
  khi = nz
elseif(npx.eq.1.and.npy.ne.1.and.npz.ne.1)then !yz decomp
  node_y = mod(node,npy)
  node_z = node/npy
  ilo = 1
  ihi = nx
  jlo = node_y*ny/npy + 1
  jhi = (node_y+1)*ny/npy
  klo = node_z*nz/npz + 1
  khi = (node_z+1)*nz/npz
elseif(npx.ne.1.and.npy.eq.1.and.npz.ne.1)then !xz decomp
  node_x = mod(node,npx)
  node_z = node/npx
  ilo = node_x*nx/npx + 1
  ihi = (node_x+1)*nx/npx
  jlo = 1
  jhi = ny
  klo = node_z*nz/npz + 1
  khi = (node_z+1)*nz/npz
elseif(npx.eq.mprocs)then !x decomp
  node_x=nx/mprocs
  ilo = node*node_x + 1
  ihi = (node+1)*node_x
  jlo = 1
  jhi = ny
  klo = 1
  khi = nz
elseif(npy.eq.mprocs)then !y decomp
  node_y=ny/mprocs
  ilo = 1
  ihi = nx
  jlo = node*node_y + 1
  jhi = (node+1)*node_y
  klo = 1
  khi = nz
elseif(npz.eq.mprocs)then !z decomp
  node_z=nz/mprocs
  ilo = 1
  ihi = nx
  jlo = 1
  jhi = ny
  klo = node*node_z + 1
  khi = (node+1)*node_z
else
  write(6,*)'(decompose) npx,npy,npz error'
  !call myexit
  stop
endif

ilo=ilo+nxlo-1
ihi=ihi+nxlo-1
jlo=jlo+nylo-1
jhi=jhi+nylo-1
klo=klo+nzlo-1
khi=khi+nzlo-1


end subroutine 
      
      
      
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+

subroutine procgriddecomp(mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,idecomp,npx,npy,npz)
implicit none
integer :: mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,idecomp,npx,npy,npz
integer :: nx,ny,nz

nx=nxhi-nxlo+1
ny=nyhi-nylo+1
nz=nzhi-nzlo+1
npx=1
npy=1
npz=1
if (idecomp.eq.0) then
  call proc2grid3d(mprocs,nx,ny,nz,npx,npy,npz)
elseif (idecomp.eq.1) then
  call proc2grid2d(mprocs,nx,ny,npx,npy)
elseif (idecomp.eq.2) then
  call proc2grid2d(mprocs,ny,nz,npy,npz)
elseif (idecomp.eq.3) then
  call proc2grid2d(mprocs,nx,nz,npx,npz)
elseif (idecomp.eq.4) then
  npx=nx/mprocs
elseif (idecomp.eq.5) then
  npy=ny/mprocs
elseif (idecomp.eq.6) then
  npz=nz/mprocs
else
  write(6,*)'(procgriddecomp) error, this routine is called only when 0 <= idecomp <= 6, but idecomp=',idecomp
endif

end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! assign mprocs to 3-d nx,ny,nz grid so as to minimize surface area
! return 3-d px,py,pz grid of procs
subroutine proc2grid3d(mprocs,nx,ny,nz,px,py,pz)
implicit none

integer mprocs
integer nx,ny,nz
integer px,py,pz
integer surf,boxx,boxy,boxz
integer bestsurf,bestboxx,bestboxy,bestboxz
integer ipx,ipy,ipz,nremain
bestsurf = 2 * (nx*ny + ny*nz + nz*nx)
bestboxx = 0
bestboxy = 0
bestboxz = 0
! loop thru all possible factorizations of mprocs
! surf = surface area of largest proc sub-domain
! innermost if test minimizes surface area and surface/volume ratio
ipx = 1
do while (ipx.le.mprocs)
  if (mod(mprocs,ipx).eq.0) then
    nremain = mprocs/ipx
    ipy = 1
    do while (ipy.le.nremain)
      if (mod(nremain,ipy).eq.0) then
        ipz = nremain/ipy
        boxx = nx/ipx
        if (mod(nx,ipx).ne.0) boxx = boxx + 1
        boxy = ny/ipy
        if (mod(ny,ipy).ne.0) boxy = boxy + 1
        boxz = nz/ipz
        if (mod(nz,ipz).ne.0) boxz = boxz + 1
        surf = boxx*boxy + boxy*boxz + boxz*boxx
        if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and.                    &
             boxx*boxy*boxz.gt.bestboxx*bestboxy*bestboxz)) then
          bestsurf = surf
          bestboxx = boxx
          bestboxy = boxy
          bestboxz = boxz
          px = ipx
          py = ipy
          pz = ipz
        endif
      endif
      ipy = ipy + 1
    enddo
  endif
  ipx = ipx + 1
enddo

if (px*py*pz.ne.mprocs) then
  write (6,*) 'Bad result in proc2grid3d'
  !call exit(0)
  stop
endif

end subroutine 

      
      
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+      
! assign mprocs to 2d nx,ny grid so as to minimize surface area
! return 2d px,py grid of procs

subroutine proc2grid2d(mprocs,nx,ny,px,py)
implicit none

integer mprocs
integer nx,ny
integer px,py

integer surf,boxx,boxy
integer bestsurf,bestboxx,bestboxy
integer ipx,ipy

bestsurf = 2 * (nx + ny)
bestboxx = 0
bestboxy = 0

! loop thru all possible factorizations of mprocs
! surf = surface area of largest proc sub-domain
! innermost if test minimizes surface area and surface/volume ratio

ipx = 1
do while (ipx.le.mprocs)
  if (mod(mprocs,ipx).eq.0) then
    ipy = mprocs/ipx
    boxx = nx/ipx
    if (mod(nx,ipx).ne.0) boxx = boxx + 1
    boxy = ny/ipy
    if (mod(ny,ipy).ne.0) boxy = boxy + 1
    surf = boxx + boxy
    if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and. &
         boxx*boxy.gt.bestboxx*bestboxy)) then
      bestsurf = surf
      bestboxx = boxx
      bestboxy = boxy
      px = ipx
      py = ipy
    endif
  endif
  ipx = ipx + 1
enddo

if (px*py.ne.mprocs) then
  write (6,*) 'Bad result in proc2grid2d'
  !call exit(0)
  stop
endif


end subroutine

end module

