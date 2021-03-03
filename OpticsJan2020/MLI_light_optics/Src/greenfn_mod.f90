!***********************************************************************
!
! intgreenfn: module for computing integrated Green functions
!
! Description: This module implements the derived types and subroutines
! for computing 'effective integrated Green functions'.
!
! Version: 0.1
! Author: D.T.Abell, Tech-X Corp., Jan.2007
!
! Comments
!   24.Jan.07 DTA: This module currently implements only the 3-D
!     integrated Green functions for open boundary conditions.
!
!***********************************************************************
      module intgreenfn
        !use parallel, only : idproc
        implicit none
!
! data
!
      double precision, parameter, private :: dtiny=1.d-13
!
! derived types
!
      type griddata
        integer :: NxIntrvl,NyIntrvl,NzIntrvl
        integer :: Nx,Ny,Nz,Nx2,Ny2,Nz2
        double precision  :: xmin,xmax,ymin,ymax,zmin,zmax
        double precision  :: hx,hy,hz
        double precision  :: Vh,Vhinv
      end type griddata
!
! functions and subroutines
!
      contains
!
!***********************************************************************
      subroutine intgreenphi(grid,g)
      use math_consts
      implicit none
      type(griddata), intent(in) :: grid
!     complex*16, dimension(:,:,:), intent(out) :: g
      complex*16,dimension(grid%Nx2,grid%Ny2,grid%Nz2),intent(out) :: g

      integer :: ix,nx1,nx12
      integer :: iy,ny1,ny12
      integer :: iz,nz1,nz12
      double precision :: invfourpi,gg,vh2i

      invfourpi=1/fourpi
      vh2i = grid%vhinv ** 2
      nx1 = grid%Nx + 1; nx12 = 2 * nx1
      ny1 = grid%Ny + 1; ny12 = 2 * ny1
      nz1 = grid%Nz + 1; nz12 = 2 * nz1
      do iz = 1,nz1
        do iy = 1,ny1
          do ix = 1,nx1
            call geff3phi(grid%hx,grid%hy,grid%hz,ix-1,iy-1,iz-1,gg)
            !g(ix,iy,iz) = vh2i * gg
            g(ix,iy,iz) = invfourpi * gg
          end do
          do ix = nx1+1,grid%Nx2
            g(ix,iy,iz) = g(nx12-ix,iy,iz)
          end do
        end do
        do iy = ny1+1,grid%Ny2
          do ix = 1,grid%Nx2
            g(ix,iy,iz) = g(ix,ny12-iy,iz)
          end do
        end do
      end do
      do iz = nz1+1,grid%Nz2
        do iy = 1,grid%Ny2
          do ix = 1,grid%Nx2
            g(ix,iy,iz) = g(ix,iy,nz12-iz)
          end do
        end do
      end do

      return
      end subroutine intgreenphi
!
!***********************************************************************
      subroutine geff3phi(hx,hy,hz,i,j,k,g)
      implicit none
      double precision, intent(in) :: hx,hy,hz
      integer, intent(in) :: i,j,k
      double precision, intent(out) :: g

      double precision :: g000,g001,g010,g011,g100,g101,g110,g111

      call intg3(hx,hy,hz,i,j,k,0,0,0,g000)
      call intg3(hx,hy,hz,i,j,k,0,0,1,g001)
      call intg3(hx,hy,hz,i,j,k,0,1,0,g010)
      call intg3(hx,hy,hz,i,j,k,0,1,1,g011)
      call intg3(hx,hy,hz,i,j,k,1,0,0,g100)
      call intg3(hx,hy,hz,i,j,k,1,0,1,g101)
      call intg3(hx,hy,hz,i,j,k,1,1,0,g110)
      call intg3(hx,hy,hz,i,j,k,1,1,1,g111)

      g = g000 - (g001 + g010 + g100) + (g011 + g101 + g110) - g111
      g = g/(hx*hy*hz)**2

      return
      end subroutine geff3phi
!
!***********************************************************************
      subroutine intg3(hx,hy,hz,i,j,k,r,s,t,g)
      implicit none
      double precision, intent(in) :: hx,hy,hz
      integer, intent(in) :: i,j,k,r,s,t
      double precision, intent(out) :: g

      double precision :: a,b,c,xl,xu,yl,yu,zl,zu

      a = hx*(i - 1 + 2*r); xl = hx*(i - 1 + r); xu = hx*(i + r)
      b = hy*(j - 1 + 2*s); yl = hy*(j - 1 + s); yu = hy*(j + s)
      c = hz*(k - 1 + 2*t); zl = hz*(k - 1 + t); zu = hz*(k + t)

      call defintg3(a,b,c,xl,xu,yl,yu,zl,zu,g)

      return
      end subroutine intg3
!
!***********************************************************************
      subroutine defintg3(a,b,c,xl,xu,yl,yu,zl,zu,g)
      implicit none
      double precision, intent(in) :: a,b,c
      double precision, intent(in) :: xl,xu,yl,yu,zl,zu
      double precision, intent(out) :: g

      double precision :: aa,ba,ca
      double precision :: lll,llu,lul,ull,uul,ulu,luu,uuu

      aa=abs(a); ba=abs(b); ca=abs(c)
      if (aa.lt.dtiny.and.ba.lt.dtiny.and.ca.lt.dtiny) then
        call indef000(xl,yl,zl,lll)
        call indef000(xl,yl,zu,llu)
        call indef000(xl,yu,zl,lul)
        call indef000(xl,yu,zu,luu)
        call indef000(xu,yl,zl,ull)
        call indef000(xu,yl,zu,ulu)
        call indef000(xu,yu,zl,uul)
        call indef000(xu,yu,zu,uuu)
      else if (aa.lt.dtiny.and.ba.lt.dtiny) then
        call indef00c(c,xl,yl,zl,lll)
        call indef00c(c,xl,yl,zu,llu)
        call indef00c(c,xl,yu,zl,lul)
        call indef00c(c,xl,yu,zu,luu)
        call indef00c(c,xu,yl,zl,ull)
        call indef00c(c,xu,yl,zu,ulu)
        call indef00c(c,xu,yu,zl,uul)
        call indef00c(c,xu,yu,zu,uuu)
      else if (aa.lt.dtiny.and.ca.lt.dtiny) then
        call indef0b0(b,xl,yl,zl,lll)
        call indef0b0(b,xl,yl,zu,llu)
        call indef0b0(b,xl,yu,zl,lul)
        call indef0b0(b,xl,yu,zu,luu)
        call indef0b0(b,xu,yl,zl,ull)
        call indef0b0(b,xu,yl,zu,ulu)
        call indef0b0(b,xu,yu,zl,uul)
        call indef0b0(b,xu,yu,zu,uuu)
      else if (aa.lt.dtiny) then
        call indef0bc(b,c,xl,yl,zl,lll)
        call indef0bc(b,c,xl,yl,zu,llu)
        call indef0bc(b,c,xl,yu,zl,lul)
        call indef0bc(b,c,xl,yu,zu,luu)
        call indef0bc(b,c,xu,yl,zl,ull)
        call indef0bc(b,c,xu,yl,zu,ulu)
        call indef0bc(b,c,xu,yu,zl,uul)
        call indef0bc(b,c,xu,yu,zu,uuu)
      else if (ba.lt.dtiny.and.ca.lt.dtiny) then
        call indefa00(a,xl,yl,zl,lll)
        call indefa00(a,xl,yl,zu,llu)
        call indefa00(a,xl,yu,zl,lul)
        call indefa00(a,xl,yu,zu,luu)
        call indefa00(a,xu,yl,zl,ull)
        call indefa00(a,xu,yl,zu,ulu)
        call indefa00(a,xu,yu,zl,uul)
        call indefa00(a,xu,yu,zu,uuu)
      else if (ba.lt.dtiny) then
        call indefa0c(a,c,xl,yl,zl,lll)
        call indefa0c(a,c,xl,yl,zu,llu)
        call indefa0c(a,c,xl,yu,zl,lul)
        call indefa0c(a,c,xl,yu,zu,luu)
        call indefa0c(a,c,xu,yl,zl,ull)
        call indefa0c(a,c,xu,yl,zu,ulu)
        call indefa0c(a,c,xu,yu,zl,uul)
        call indefa0c(a,c,xu,yu,zu,uuu)
      else if (ca.lt.dtiny) then
        call indefab0(a,b,xl,yl,zl,lll)
        call indefab0(a,b,xl,yl,zu,llu)
        call indefab0(a,b,xl,yu,zl,lul)
        call indefab0(a,b,xl,yu,zu,luu)
        call indefab0(a,b,xu,yl,zl,ull)
        call indefab0(a,b,xu,yl,zu,ulu)
        call indefab0(a,b,xu,yu,zl,uul)
        call indefab0(a,b,xu,yu,zu,uuu)
      else
        call indefabc(a,b,c,xl,yl,zl,lll)
        call indefabc(a,b,c,xl,yl,zu,llu)
        call indefabc(a,b,c,xl,yu,zl,lul)
        call indefabc(a,b,c,xl,yu,zu,luu)
        call indefabc(a,b,c,xu,yl,zl,ull)
        call indefabc(a,b,c,xu,yl,zu,ulu)
        call indefabc(a,b,c,xu,yu,zl,uul)
        call indefabc(a,b,c,xu,yu,zu,uuu)
      end if
      g = uuu - (uul + ulu + luu) + (llu + lul + ull) - lll

      return
      end subroutine defintg3
!
!***********************************************************************
      subroutine aux23(a,b,v,r)
      implicit none
      double precision, intent(in) :: a,b,v
      double precision, intent(out) :: r

      double precision :: d

      d = 2*v - 3*b
      r = 0
      if (abs(d).gt.dtiny) r = d*atan(v*(3*v - 4*b)/(4*a*d))
      
      return
      end subroutine aux23
!
!***********************************************************************
      subroutine aux34(b,c,v,w,v2,w2,vpr,r)
      implicit none
      double precision, intent(in) :: b,c,v,w,v2,w2,vpr
      double precision, intent(out) :: r

      double precision :: d

      d = 3*w - 4*c
      r = 0
      if (abs(d).gt.0.d0) r = d*log((v2+w2)*(4*vpr/(b*d*w*v2*w2))**2)
      
      return
      end subroutine aux34
!
!***********************************************************************
      subroutine indef000(u,v,w,g)
      implicit none
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      g = (u**2 + v**2 + w**2)**2.5d0/15.d0

      return
      end subroutine indef000
!
!***********************************************************************
      subroutine indef00c(c,u,v,w,g)
      implicit none
      double precision, intent(in) :: c
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2
      double precision :: u2v2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny) then
        g = (4*w2 - 5*c*w)*w2*wa/60.d0
      else if (va.lt.dtiny) then
        g = (32*r2**2.5d0 + 5*c*(3*u2*(3*r2 + w2) - 4*r*w*(5*u2 + 2*    &
     &       w2)) - 60*c*u2**2*log(wpr))/480.d0
      else if (wa.lt.dtiny) then
        g = (32*r2**2.5d0 + 45*c*u2*(r2 + v2) - 30*c*r2**2*log(r2))/    &
     &      480.d0
      else
        u2v2 = u2 + v2
        g = (15*c*u2*(3*(r2 + v2) + w2) + 4*r*(8*r2**2 - 5*c*w*(2*r2 +  &
     &       3*u2v2)) - 60*c*u2*(u2v2 + v2)*log(wpr) - 30*c*v2**2*      &
     &       log((v2 + w2) * ((4*wpr)/(3*c*v2**2*w2))**2))/480.d0
      end if

      return
      end subroutine indef00c
!
!***********************************************************************
      subroutine indef0b0(b,u,v,w,g)
      implicit none
      double precision, intent(in) :: b
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.wa.lt.dtiny) then
        g = v2*(4*v2 - 5*b*v)*va/60.d0
      else if (va.lt.dtiny) then
        g = (32*r2**2.5d0 + 15*b*u2*(r2 + w2) - 30*b*r2**2*log(r2))/    &
     &      480.d0
      else if (wa.lt.dtiny) then
        g = (32*r2**2.5d0 - 20*b*r*v*(5*u2 + 2*v2) + 15*b*u2**2*(1 - 4* &
     &       log(vpr)))/480.d0
      else
        g = (4*r*(8*r2**2 - 5*b*v*(5*(u2 + w2) + 2*v2)) + 15*b*u2*(u2 + &
     &       2*w2)*(1 - 4*log(vpr)) - 30*b*w2**2*log((16*(v2 +          &
     &       w2)*vpr**2)/(9*b**2*v2**2*w2**4)))/480.d0
      end if

      return
      end subroutine indef0b0
!
!***********************************************************************
      subroutine indef0bc(b,c,u,v,w,g)
      implicit none
      double precision, intent(in) :: b,c
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2
      double precision :: aux1,aux2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny.and.wa.lt.dtiny) then
        g = 0.d0
      else if (ua.lt.dtiny.and.va.lt.dtiny) then
        g = w2*w*((16*w - 20*c)*wa - 15*b*(w - 4*c)*log(w2))/240.d0
      else if (ua.lt.dtiny.and.wa.lt.dtiny) then
        g = v2*v*((16*v - 20*b)*va - 5*c*(3*v - 4*b)*log(v2))/240.d0
      else if (va.lt.dtiny) then
        g = (15*u2*(c*(3*r2 + w2) + b*(r2 + w2 - 8*c*w)) + 4*r*(8*      &
     &       r2**2 - 5*c*w*(5*u2 + 2*w2)) - 30*b*r2*(r2 - 4*c*w)*       &
     &       log(r2) - 60*c*u2**2*log(wpr))/480.d0
      else if (wa.lt.dtiny) then
        g = (5*u2*(b*(3*u2 - 56*c*v) + 9*c*(r2 + v2)) + 4*r*(8*r2**2 -  &
     &       5*b*v*(5*u2 + 2*v2)) + 160*b*c*u2*u*atan(v/u) - 10*c*(3*   &
     &       r2**2 - 4*b*v*(r2 + 2*u2))*log(r2) - 60*b*u2**2*log(vpr))/ &
     &      480.d0
      else if (ua.lt.dtiny) then
        call aux34(b,c,v,w,v2,w2,vpr,aux1)
        call aux34(c,b,w,v,w2,v2,wpr,aux2)
        g = (2*r*(8*r2**2 + 5*(8*b*c*v*w - b*v*(5*w2 + 2*v2) - c*w*(5*  &
     &       v2 + 2*w2))) + 40*b*c*w*w2*log(w2) - 5*b*w*w2*aux1 - 5*c*  &
     &       v*v2*aux2)/240.d0
      else
        call aux34(b,c,v,w,v2,w2,vpr,aux1)
        call aux34(c,b,w,v,w2,v2,wpr,aux2)
        g = (u2*(3*c*(3*(r2 + v2) + w2) + b*(3*u2 + 6*w2 - 8*c*(7*v +   &
     &       3*w))) + 0.8d0*r*(8*r2**2 - 5*(b*v*(2*r2 + 3*(u2 + w2)) +  &
     &       c*w*(2*r2 + 3*(u2 + v2)) - 8*b*c*v*w)) + 32*b*c*u*u2*      &
     &       (atan(v/u) - atan((v*w)/(u*r))) - 12*u2*(b*(u2 + 2*w*(w -  &
     &       2*c))*log(vpr) + c*(u2 + 2*v*(v - 2*b))*log(wpr)) + 16*b*  &
     &       c*w*w2*log(u2 + w2) - 2*b*w*w2*aux1 - 2*c*v*v2*aux2)/96.d0
      end if

      return
      end subroutine indef0bc
!
!***********************************************************************
      subroutine indefa00(a,u,v,w,g)
      implicit none
      double precision, intent(in) :: a
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2
      double precision :: v2w2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny.and.wa.lt.dtiny) then
        g = 0.d0
      else if (va.lt.dtiny.and.wa.lt.dtiny) then
        g = r*u2*(4*u2 - 5*a*u)/60.d0
      else
        v2w2 = v2 + w2
        g = (r*(8*r2**2 - 5*a*u*(2*r2 + 3*v2w2)) - 15*a*v2w2**2*        &
     &       log(upr))/120.d0
      end if

      return
      end subroutine indefa00
!
!***********************************************************************
      subroutine indefa0c(a,c,u,v,w,g)
      implicit none
      double precision, intent(in) :: a,c
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2
      double precision :: u2v2,v2w2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny.and.wa.lt.dtiny) then
        g = 0.d0
      else if (ua.lt.dtiny.and.va.lt.dtiny) then
        g = w2*w*((16*w - 20*c)*wa - 5*a*(3*w - 4*c)*log(w2))/240.d0
      else if (va.lt.dtiny.and.wa.lt.dtiny) then
        g = u2*u*((96*u - 120*a)*ua + 5*c*((27*u - 40*a) - 6*(3*u - 4*  &
     &            a)*log(u2)))/1440.d0
      else if (va.lt.dtiny) then
        g = (5*c*u*(9*u*(3*r2 + w2) - 8*a*(5*u2 + 9*w2)) +              &
     &       12*r*(8*r2**2 + 5*(8*a*c*u*w - a*u*(2*r2 + 3*w2) - c*w*(2* &
     &       r2 + 3*u2))) + 60*a*w*(4*c - 3*w)*w2*log(upr) + 60*c*u*(4* &
     &       a - 3*u)*u2*log(wpr))/1440.d0
      else if (wa.lt.dtiny) then
        g = (5*c*u*(27*u*(r2 + v2) - 8*a*(5*r2 + 16*v2)) + 12*r*(8*     &
     &       r2**2 - 5*a*u*(2*r2 + 3*v2)) + 480*a*c*v*v2*(atan(u/v) +   &
     &       atan((3*v)/(8*a))) - 30*c*(3*r2**2 - 4*a*u*(u2 + 3*v2))*   &
     &       log(r2) - 180*a*v2**2*log(upr))/1440.d0
      else
        u2v2 = u2 + v2
        v2w2 = v2 + w2
        g = (c*u*(9*u*(3*(r2 + v2) + w2) - 8*a*(5*r2 + 16*v2 + 4*w2)) + &
     &       2.4d0*r*(8*r2**2 + 40*a*c*u*w - 5*a*u*(2*r2 + 3*v2w2) - 5* &
     &       c*w*(2*r2 + 3*u2v2)) + 96*a*c*v*v2*(atan(u/v) - atan((8*a* &
     &       u*w - 3*r*v2)/(v*(8*a*r + 3*u*w)))) - 12*a*(3*v2w2**2 -    &
     &       4*c*w*(3*v2 + w2))*log(upr) - 12*c*u*(3*u*(u2v2 + v2) -    &
     &       4*a*(u2v2 + 2*v2))*log(wpr) - 18*c*v2**2*log((v2w2/(64*    &
     &       a**2 + 9*v2))*((4*wpr)/(c*v*v2*w2))**2))/288.d0
      end if

      return
      end subroutine indefa0c
!
!***********************************************************************
      subroutine indefab0(a,b,u,v,w,g)
      implicit none
      double precision, intent(in) :: a,b
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny.and.wa.lt.dtiny) then
        g = 0.d0
      else if (ua.lt.dtiny.and.wa.lt.dtiny) then
        g = v2*v*((16*v - 20*b)*va - 5*a*(3*v - 4*b)*log(v2))/240.d0
      else if (va.lt.dtiny.and.wa.lt.dtiny) then
        g = u2*u*((96*u - 120*a)*ua + 5*b*((9*u - 16*a) - 6*(3*u - 4*  &
     &            a)*log(u2)))/1440.d0
      else if (wa.lt.dtiny) then
        g = (5*b*u2*(9*u2 - 16*a*u) + 12*r*(8*r2**2 + 40*a*b*u*v - 5*   &
     &       (a*u*(2*u2 + 5*v2) + b*v*(5*u2 + 2*v2))) + 60*(a*(4*b*v -  &
     &       3*v2)*v2*log(upr) + b*(4*a*u - 3*u2)*u2*log(vpr)))/1440.d0
      else if (va.lt.dtiny) then
        g = (5*b*u*(9*u*(r2 + w2) - 16*a*(r2 + 5*w2)) + 12*r*(8*r2**2 - &
     &       5*a*u*(2*u2 + 5*w2)) + 480*a*b*w*w2*(atan(u/w) + atan((3*  &
     &       w)/(8*a))) - 30*b*(3*r2**2 - 4*a*u*(r2 + 2*w2))*           &
     &       log(r2) - 180*a*w2**2*log(upr))/1440.d0
      else
        g = (5*b*u*(9*u*(u2 + 2*w2) - 16*a*(u2 + 6*w2)) + 12*r*(8*      &
     &       r2**2 + 40*a*b*u*v - 5*(a*u*(2*r2 + 3*(v2 + w2)) + b*v*(2* &
     &       r2 + 3*(u2 + w2)))) + 480*a*b*w*w2*(atan(u/w) - atan((8*   &
     &       a*u*v - 3*r*w2)/((8*a*r + 3*u*v)*w))) + 60*a*(4*b*v*(v2 +  &
     &       3*w2) - 3*(v2 + w2)**2)*log(upr) - 60*b*u*(3*u*(u2 + 2*    &
     &       w2) - 4*a*(u2 + 3*w2))*log(vpr) - 90*b*w2**2*log(((v2 +    &
     &       w2)/(64*a**2 + 9*w2))*((4*vpr)/(b*w*v2*w2))**2))/1440.d0
      end if

      return
      end subroutine indefab0
!
!***********************************************************************
      subroutine indefabc(a,b,c,u,v,w,g)
      implicit none
      double precision, intent(in) :: a,b,c
      double precision, intent(in) :: u,v,w
      double precision, intent(out) :: g

      double precision :: r,r2,ua,upr,u2,va,vpr,v2,wa,wpr,w2
      double precision :: aux1,aux2

      ua=abs(u); va=abs(v); wa=abs(w);
      u2=u**2; v2=v**2; w2=w**2;
      r2=u2+v2+w2; r=sqrt(r2);
      upr=u+r; vpr=v+r; wpr=w+r;

      if (ua.lt.dtiny.and.va.lt.dtiny.and.wa.lt.dtiny) then
        g = 0.d0
      else if (ua.lt.dtiny.and.va.lt.dtiny) then
        call aux23(a,c,w,aux2)
        g = w2*(4*wa*(4*w2 - 5*c*w) + 5*(8*a*b*aux2 + w*(4*c*(a + 3*    &
     &          b) - 3*w*(a + b))*log(w2)))/240.d0
      else if (ua.lt.dtiny.and.wa.lt.dtiny) then
        call aux23(a,b,v,aux1)
        g = v2*(4*va*(4*v2 - 5*b*v) + 5*(8*a*c*aux1 + (a + c)*(4*b*v -  &
     &          3*v2)*log(v2)))/240.d0
      else if (va.lt.dtiny.and.wa.lt.dtiny) then
        g = u*u2*(24*(4*u - 5*a)*ua + 5*(9*(b + 3*c)*u - 8*a*(2*b +     &
     &            5*c) + 6*(b + c)*(4*a - 3*u)*log(u2)))/1440.d0
      else if (va.lt.dtiny) then
        call aux23(a,c,w,aux2)
        g = (12*r*(8*r2**2 + 40*a*c*u*w - 5*(a*u*(2*u2 + 5*w2) + c*w*   &
     &       (5*u2 + 2*w2))) + 5*u*(9*u*((b + 3*c)*u2 + 2*(b + 2*c)*    &
     &       w2) - 8*a*(2*b + 5*c)*u2 + 72*b*c*(4*a - u)*w - 24*a*(4*   &
     &       b + 3*c)*w2) + 240*a*b*w2*((2*w - 6*c)*atan(u/w) + aux2) + &
     &       30*b*(4*a*u*(u2 - 6*c*w + 3*w2) + 3*r2*(4*c*w - r2))*      &
     &       log(r2) + 60*(a*w*(4*c - 3*w)*w2*log(upr) + c*u*(4*a - 3*  &
     &       u)*u2*log(wpr)))/ 1440.d0
      else if (wa.lt.dtiny) then
        call aux23(a,b,v,aux1)
        g = (5*u*(3*u*(3*(b + 3*c)*u2 + c*(18*v2 - 56*b*v)) - 8*a*((2*  &
     &       b + 5*c)*u2 + c*(21*v2 - 54*b*v))) + 12*r*(8*r2**2 + 40*a* &
     &       b*u*v - 5*(a*u*(2*u2 + 5*v2) + b*v*(5*u2 + 2*v2))) - 240*  &
     &       c*(b*(3*a - 2*u)*u2*atan(v/u) + a*v2*((3*b - 2*v)*         &
     &       atan(u/v) - aux1)) + 30*c*(4*(a*u*(u2 - 3*(b - v)*v) + b*  &
     &       v*(v2 - 3*(a - u)*u)) - 3*r2**2)*log(r2) + 60*(a*v*(4*b -  &
     &       3*v)*v2*log(upr) + b*u*(4*a - 3*u)*u2*log(vpr)))/1440.d0
      else if (ua.lt.dtiny) then
        call aux23(a,b,v,aux1)
        call aux23(a,c,w,aux2)
        g = (12*r*(8*r2**2 + 40*b*c*v*w - 5*(b*v*(2*v2 + 5*w2) + c*w*   &
     &       (5*v2 + 2*w2))) + 30*a*(4*b*v*(v2 + 3*w*(w - c)) + 4*c*w*  &
     &       (w2 + 3*v*(v - b)) - 3*r2**2)*log(r2) + 240*b*c*w*w2*      &
     &       log(w2) + 240*a*(b*w2*aux2 + c*v2*aux1) + 30*b*w*(4*c - 3* &
     &       w)*w2*log((1/(16*a**2*(3*c - 2*w)**2 + (4*c - 3*w)**2*     &
     &       w2))*((4*r*vpr)/(b*v2*w2))**2) + 30*c*v*(4*b - 3*v)*v2*    &
     &       log((1/(16*a**2*(3*b - 2*v)**2 + (4*b - 3*v)**2*v2))*((4*  &
     &       r*wpr)/(c*v2*w2))**2))/1440.d0
      else
        g = (2.4d0*r*(8*r2**2 + 40*(a*b*u*v + a*c*u*w + b*c*v*w) - 5*   &
     &       (a*u*(2*u2 + 5*(v2 + w2)) + b*v*(5*(u2 + w2) + 2*v2) + c*w &
     &       *(5*(u2 + v2) + 2*w2))) + u*(3*u*(3*(b + 3*c)*u2 + 2*c*(9* &
     &       v - 28*b)*v - 24*b*c*w + 6*(b + 2*c)*w2) - 8*a*(2*b*(u2 -  &
     &       27*c*v - 18*c*w + 6*w2) + c*(5*u2 + 21*v2 + 9*w2))) - 144* &
     &       a*b*c*w2*atan(u/w) - 48*b*c*(3*a - 2*u)*u2*(atan(v/u) -    &
     &       atan((v*w)/(u*r))) - 48*a*c*(3*b - 2*v)*v2*(atan(u/v) -    &
     &       atan((4*a*u*(3*b - 2*v)*w - (4*b - 3*v)*v2*r)/(v*(u*(4*    &
     &       b - 3*v)*w + 4*a*(3*b - 2*v)*r)))) - 48*a*b*(3*c - 2*w)*   &
     &       w2*(atan(u/w) - atan((4*a*u*v*(3*c - 2*w) - (4*c - 3*w)*   &
     &       w2*r)/(w*(u*v*(4*c - 3*w) + 4*a*(3*c - 2*w)*r)))) + 12*a*  &
     &       (4*b*v*(v2 + 3*w*(w - c)) + 4*c*w*(w2 + 3*v*(v - b)) - 3*  &
     &       (v2 + w2)**2)*log(upr) - 12*b*u*(3*u*(u2 + 2*w*(w - 2*     &
     &       c)) - 4*a*(u2 + 3*w*(w - 2*c)))*log(vpr) - 12*c*u*(3*u*    &
     &       (u2 + 2*v*(v - 2*b)) - 4*a*(u2 + 3*v*(v - 2*b)))*          &
     &       log(wpr) + 48*b*c*w*w2*log(u2 + w2) + 6*b*w*(4*c - 3*w)*   &
     &       w2*log(((v2 + w2)/(16*a**2*(3*c - 2*w)**2 + (4*c - 3*w)**  &
     &       2*w2))*((4*vpr)/(b*v2*w2))**2) + 6*c*v*(4*b - 3*v)*v2*     &
     &       log(((v2 + w2)/(16*a**2*(3*b - 2*v)**2 + (4*b - 3*v)**2*   &
     &       v2))*((4*wpr)/(c*v2*w2))**2))/288.d0
      end if

      return
      end subroutine indefabc
!
!***********************************************************************
      end module intgreenfn

