c IMPACT 2D space charge routines
c Copyright 2001 University of California 
c
      subroutine spch2dphi(c,ex,ey,msk,np,ntot,nx,ny,n1,n2)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      complex*16 grnxtr,erhox,erhoxtr,rho2,rho2xtr
! rho=charge density on the grid
! rho2=...on doubled grid
! rho2xtr=...xformed (by fft) and transposed
! grnxtr=green function, xformed, transposed
      dimension c(4,maxrayp),msk(maxrayp),ex(maxrayp),ey(maxrayp)
      dimension rho(nx,ny),exg(nx,ny),eyg(nx,ny),rhosum(nx,ny)
      dimension rho2(n1,n2),rho2xtr(n2,n1),grnxtr(n2,n1)
! weights, indices associated with area weighting:
      dimension ab(maxrayp),cd(maxrayp)
      dimension indx(maxrayp),jndx(maxrayp)
      dimension indxp1(maxrayp),jndxp1(maxrayp)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/newpoisson/idirectfieldcalc,idensityfunction,idummy(4)
! compute grid quantities hx,hy,...
      call boundp2d(c,msk,np,nx,ny)
! deposit charge on the grid:
      call rhoslo2d(c,rho,msk,np,ab,cd,indx,jndx,indxp1,jndxp1,nx,ny)
      call MPI_ALLREDUCE(rho,rhosum,nx*ny,mreal,mpisum,lworld,ierror)
      rho=rhosum
! store rho in lower left quadrant of doubled grid:
      rho2=(0.,0.)
      do j=1,ny
      do i=1,nx
      rho2(i,j)=cmplx(rho(i,j),0.)
      enddo
      enddo
!---------convolution------------
      twopi=4.d0*asin(1.d0)
! compute FFT of the Green function on the grid:
      if(idensityfunction.eq.0)then
!       if(idproc.eq.0)write(6,*)'in spch2d (via phi); NOT using slic'
        call greenf2d(grnxtr,nx,ny,n1,n2)
        grnxtr=0.5d0*grnxtr/(twopi*nrays)
      else
!       if(idproc.eq.0)write(6,*)'in spch2d, using slic (for phi)'
        call greenphi2d(grnxtr,rho,nx,ny,n1,n2)
        grnxtr=-1.d0*grnxtr/(hx*hy*twopi*nrays)
      endif
! fft the charge density:
      scale=1.0
      call fft2dhpf(n1,n2,1,0,scale,0,rho2,rho2xtr)
! multiply transformed charge density and transformed Green function:
      rho2xtr=rho2xtr*grnxtr/(n1*n2)
! inverse fft:
      call fft2dhpf(n2,n1,-1,0,scale,0,rho2xtr,rho2)
!----done with convolution-------
! store physical data back on grid of correct (not doubled) size:
      do j=1,ny
      do i=1,nx
      rho(i,j)=real(rho2(i,j))
      enddo
      enddo
! obtain the electric fields:
      exg=0.5*hxi*(cshift(rho,-1,1)-cshift(rho,1,1))
      eyg=0.5*hyi*(cshift(rho,-1,2)-cshift(rho,1,2))
! interpolate electric field at particle postions:
      call ntrslo2d(exg,eyg,ex,ey,msk,ab,cd,indx,jndx,indxp1,jndxp1,       &
     &              nx,ny,np)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     do i=1,nx
!      j=ny/2
!      x=xmin+(i-1)*hx
!      y=ymin+(j-1)*hy
!      write(58,*)x,y,exg(i,j),eyg(i,j)
!     enddo
!     do j=1,ny
!      i=nx/2
!      x=xmin+(i-1)*hx
!      y=ymin+(j-1)*hy
!      write(59,*)x,y,exg(i,j),eyg(i,j)
!     enddo
!     if(nx.ne.12345)call myexit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end
!
      subroutine boundp2d(c,msk,np,nx,ny)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      dimension c(4,maxrayp),msk(maxrayp)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
! transverse:
      xbig=maxval(c(1,:),1,msk)
      xsml=minval(c(1,:),1,msk)
      ybig=maxval(c(3,:),1,msk)
      ysml=minval(c(3,:),1,msk)
! note: this would be a good place to check which particles
! are inside the beampipe, then mask off those that are not.
      xmin=xsml-0.05*(xbig-xsml)
      xmax=xbig+0.05*(xbig-xsml)
      ymin=ysml-0.05*(ybig-ysml)
      ymax=ybig+0.05*(ybig-ysml)
      hx=(xmax-xmin)/(nx-1)
      hy=(ymax-ymin)/(ny-1)
      hxi=1./hx
      hyi=1./hy
      return
      end
!
      subroutine greenf2d(grnxtr,nx,ny,n1,n2)
      implicit double precision(a-h,o-z)
      complex*16 grn,grnxtr
      dimension grn(n1,n2)
      dimension grnxtr(n2,n1)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      eps=1.d-25
      do j=-ny,ny-1
      do i=-nx,nx-1
      grn(1+mod(i+n1,n1),1+mod(j+n2,n2))= -log((hx*i)**2+(hy*j)**2+eps)
      enddo
      enddo
      grn(1,1)=grn(1,2)
! compute FFT of the Green function
      scale=1.d0
      call fft2dhpf(n1,n2,1,0,scale,0,grn,grnxtr)
      return
      end
!
!
      subroutine rhoslo2d(c,rho,msk,np,ab,cd,indx,jndx,indxp1,           &
     &                    jndxp1,nx,ny)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      dimension c(4,maxrayp),msk(maxrayp),vol(maxrayp)
      dimension ab(maxrayp),cd(maxrayp)
      dimension indx(maxrayp),jndx(maxrayp)
      dimension indxp1(maxrayp),jndxp1(maxrayp)
      dimension rho(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      do j=1,np
      indx(j)=(c(1,j)-xmin)*hxi + 1
      jndx(j)=(c(3,j)-ymin)*hyi + 1
      enddo
      do j=1,np
      indxp1(j)=indx(j)+1
      jndxp1(j)=jndx(j)+1
      enddo
!-------
      imin=minval(indx,1,msk)
      imax=maxval(indx,1,msk)
      jmin=minval(jndx,1,msk)
      jmax=maxval(jndx,1,msk)
      if((imin.lt.1).or.(imax.gt.nx-1))then
        write(6,*)'error in rhoslo2d: imin,imax=',imin,imax
        call myexit
      endif
      if((jmin.lt.1).or.(jmax.gt.ny-1))then
        write(6,*)'error in rhoslo3d: jmin,jmax=',jmin,jmax
        call myexit
      endif
!-------
      do j=1,np
      ab(j)=((xmin-c(1,j))+indx(j)*hx)*hxi
      cd(j)=((ymin-c(3,j))+jndx(j)*hy)*hyi
      enddo
      rho=0.
!1 (i,j):
      vol=ab*cd
      do n=1,np
      rho(indx(n),jndx(n))=rho(indx(n),jndx(n))+vol(n)
      enddo
!2 (i,j+1):
      vol=ab*(1.-cd)
      do n=1,np
      rho(indx(n),jndxp1(n))=rho(indx(n),jndxp1(n))+vol(n)
      enddo
!3 (i+1,j):
      vol=(1.-ab)*cd
      do n=1,np
      rho(indxp1(n),jndx(n))=rho(indxp1(n),jndx(n))+vol(n)
      enddo
!4 (i+1,j+1):
      vol=(1.-ab)*(1.-cd)
      do n=1,np
      rho(indxp1(n),jndxp1(n))=rho(indxp1(n),jndxp1(n))+vol(n)
      enddo
!     ngood=count(msk)
!     rho=rho/ngood
      return
      end
!
      subroutine ntrslo2d(exg,eyg,ex,ey,msk,ab,cd,indx,jndx,indxp1,
     &jndxp1,nx,ny,np)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      dimension exg(nx,ny),eyg(nx,ny)
      dimension ex(maxrayp),ey(maxrayp),msk(maxrayp)
      dimension ab(maxrayp),cd(maxrayp)
      dimension indx(maxrayp),jndx(maxrayp)
      dimension indxp1(maxrayp),jndxp1(maxrayp)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      do n=1,np
      ex(n)=
     &       exg(indx(n),jndx(n))*ab(n)*cd(n)
     &      +exg(indx(n),jndxp1(n))*ab(n)*(1.-cd(n))
     &      +exg(indxp1(n),jndx(n))*(1.-ab(n))*cd(n)
     &      +exg(indxp1(n),jndxp1(n))*(1.-ab(n))*(1.-cd(n))
      enddo
      do n=1,np
      ey(n)=
     &       eyg(indx(n),jndx(n))*ab(n)*cd(n)
     &      +eyg(indx(n),jndxp1(n))*ab(n)*(1.-cd(n))
     &      +eyg(indxp1(n),jndx(n))*(1.-ab(n))*cd(n)
     &      +eyg(indxp1(n),jndxp1(n))*(1.-ab(n))*(1.-cd(n))
      enddo
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
!                2D SPACE CHARGE ROUTINES BEGIN HERE
!=======================================================================
      subroutine spch2ddirect(c,ex,ey,msk,np,ntot,nx,ny,n1,n2)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      double complex grnxtr
      double complex rho2,rho2xtr,rho2tmp
! rho=charge density on the grid
! rho2=...on doubled grid
! rho2xtr=...xformed (by fft) and transposed
! grnxtr=green function, xformed, transposed
      dimension c(4,maxrayp),msk(maxrayp),ex(maxrayp),ey(maxrayp)
      dimension rho(nx,ny),exg(nx,ny),eyg(nx,ny),rhosum(nx,ny)
      dimension rho2(n1,n2),rho2xtr(n2,n1),rho2tmp(n2,n1)
      dimension grnxtr(n2,n1)
! weights, indices associated with area weighting:
      dimension ab(maxrayp),cd(maxrayp)
      dimension indx(maxrayp),jndx(maxrayp)
      dimension indxp1(maxrayp),jndxp1(maxrayp)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/newpoisson/idirectfieldcalc,idensityfunction,idummy(4)
!     if(idproc.eq.0)write(6,*)'inside spch2ddirect; np,ntot,nx,ny='
!     if(idproc.eq.0)write(6,*)np,ntot,nx,ny

! compute grid quantities hx,hy,...
      call boundp2d(c,msk,np,nx,ny)
! deposit charge on the grid:
      call rhoslo2d(c,rho,msk,np,ab,cd,indx,jndx,indxp1,jndxp1,nx,ny)
      call MPI_ALLREDUCE(rho,rhosum,nx*ny,mreal,mpisum,lworld,ierror)
      rho=rhosum
!
! store rho in lower left quadrant of doubled grid:
      rho2=(0.,0.)
      do j=1,ny
        do i=1,nx
          rho2(i,j)=cmplx(rho(i,j),0.)
        enddo
      enddo
!---------convolution------------
! fft the charge density:
      scale=1.0
      call fft2dhpf(n1,n2,1,0,scale,0,rho2,rho2xtr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute FFT of the x-Green function on the grid:
      twopi=4.d0*asin(1.d0)
      if(idensityfunction.eq.0)then
!       if(idproc.eq.0)write(6,*)'in spch2ddirect; NOT using slic'
        call greenfxold(grnxtr,rhosum,nx,ny,n1,n2)
        grnxtr=grnxtr/(twopi*nrays)
      else
!       if(idproc.eq.0)write(6,*)'in spch2ddirect; using slic'
        call greenfx2d(grnxtr,rhosum,nx,ny,n1,n2)
        grnxtr=grnxtr/(twopi*nrays)
      endif
! multiply transformed charge density and transformed Green function:
      rho2tmp=rho2xtr*grnxtr/(n1*n2)
! inverse fft:
      call fft2dhpf(n2,n1,-1,0,scale,0,rho2tmp,rho2)
! store physical data back on grid of correct (not doubled) size:
      do j=1,ny
        do i=1,nx
          rho(i,j)=real(rho2(i,j))
        enddo
      enddo
! obtain the x-electric field:
!     exg=rho/(4.d0*asin(1.d0))
      exg=rho/(hx*hy)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute FFT of the y-Green function on the grid:
      if(idensityfunction.eq.0)then
        call greenfy2d(grnxtr,rhosum,nx,ny,n1,n2)
        grnxtr=grnxtr/(twopi*nrays)
      else
        call greenfyold(grnxtr,rhosum,nx,ny,n1,n2)
        grnxtr=grnxtr/(twopi*nrays)
      endif
! multiply transformed charge density and transformed Green function:
      rho2tmp=rho2xtr*grnxtr/(n1*n2)
! inverse fft:
      call fft2dhpf(n2,n1,-1,0,scale,0,rho2tmp,rho2)
! store physical data back on grid of correct (not doubled) size:
      do j=1,ny
        do i=1,nx
          rho(i,j)=real(rho2(i,j))
        enddo
      enddo
! obtain the y-electric field:
!     eyg=rho/(4.d0*asin(1.d0))
      eyg=rho/(hx*hy)
!
! interpolate electric field at particle postions:
      call ntrslo2d(exg,eyg,ex,ey,msk,ab,cd,indx,jndx,indxp1,jndxp1,       &
     &              nx,ny,np)
!     if(idproc.eq.0)write(6,*)'done with interpolation'
!     do j=1,np
!       ex(j)=0.5*ex(j)
!     enddo
!     do j=1,np
!       ey(j)=0.5*ey(j)
!     enddo
      return
      end
!
!
      subroutine spch2dold(c,ex,ey,msk,np,exg,eyg,nx,ny,n1,n2)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      double complex grnxtr
      double complex rho2,rho2xtr,rho2tmp
! rho=charge density on the grid
! rho2=...on doubled grid
! rho2xtr=...xformed (by fft) and transposed
! grnxtr=green function, xformed, transposed
      dimension c(4,maxrayp),msk(maxrayp),ex(maxrayp),ey(maxrayp)
      dimension rho(nx,ny),exg(nx,ny),eyg(nx,ny),rhosave(nx,ny)
      dimension rho2(n1,n2),rho2xtr(n2,n1),rho2tmp(n2,n1)
      dimension grnxtr(n2,n1)
! weights, indices associated with area weighting:
      dimension ab(maxrayp),cd(maxrayp)
      dimension indx(maxrayp),jndx(maxrayp)
      dimension indxp1(maxrayp),jndxp1(maxrayp)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      if(idproc.eq.0)write(6,*)'inside spch2dold'
! deposit charge on the grid:
!     call rho2d(c,rho,msk,ab,de,indx,jndx,indxp1,jndxp1,nx,ny)
      call rhoslo2d(c,rho,msk,np,ab,cd,indx,jndx,indxp1,jndxp1,nx,ny)
!
! store rho in lower left quadrant of doubled grid:
      rhosave(:,:)=rho(:,:)
      rho2=(0.,0.)
      do j=1,ny
        do i=1,nx
          rho2(i,j)=cmplx(rho(i,j),0.)
        enddo
      enddo
!---------convolution------------
! fft the charge density:
      scale=1.0
      call fft2dhpf(n1,n2,1,0,scale,0,rho2,rho2xtr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute FFT of the x-Green function on the grid:
      if(idproc.eq.0)write(6,*)'calling greenfxold'
      call greenfxold(grnxtr,rhosave,nx,ny,n1,n2)
      if(idproc.eq.0)write(6,*)'done with greenfxold'
! multiply transformed charge density and transformed Green function:
      rho2tmp=rho2xtr*grnxtr/(n1*n2)
! inverse fft:
      call fft2dhpf(n2,n1,-1,0,scale,0,rho2tmp,rho2)
      if(idproc.eq.0)write(6,*)'done with inverse fft'
! store physical data back on grid of correct (not doubled) size:
      do j=1,ny
        do i=1,nx
          rho(i,j)=real(rho2(i,j))
        enddo
      enddo
! obtain the x-electric field:
      exg=rho/(4.d0*asin(1.d0))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute FFT of the y-Green function on the grid:
      if(idproc.eq.0)write(6,*)'calling greenfyold'
      call greenfyold(grnxtr,rhosave,nx,ny,n1,n2)
      if(idproc.eq.0)write(6,*)'done with greenfyold'
! multiply transformed charge density and transformed Green function:
      rho2tmp=rho2xtr*grnxtr/(n1*n2)
! inverse fft:
      call fft2dhpf(n2,n1,-1,0,scale,0,rho2tmp,rho2)
      if(idproc.eq.0)write(6,*)'done with inverse fft'
! store physical data back on grid of correct (not doubled) size:
      do j=1,ny
        do i=1,nx
          rho(i,j)=real(rho2(i,j))
        enddo
      enddo
! obtain the y-electric field:
      eyg=rho/(4.d0*asin(1.d0))
!
! interpolate electric field at particle postions:
      call ntrslo2d(exg,eyg,ex,ey,msk,ab,cd,indx,jndx,indxp1,jndxp1,       &
     &              nx,ny,np)
!     if(idproc.eq.0)write(6,*)'done with interpolation'
      return
      end
!
!=============================================================
!
      subroutine greenphi2d(gxtr,rho,nx,ny,n1,n2)
! green function routine.
      implicit double precision(a-h,o-z)
      double complex g,gxtr
      dimension g(n1,n2),gxtr(n2,n1)
      dimension rho(nx,ny),exg(nx,ny),fxg(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi

      xint1(x,y,a,b,hx,hy,e2)=
     &  (((a - hx)*x**3)/9. - x**4/16. + 
     &    (x**2*(28*b - 28*hy - 9*y)*y)/24. - 
     &    ((a - hx)*x*(18*b - 18*hy - 7*y)*y)/6. + 
     &    ((3*a*b*y**2 - 3*b*hx*y**2 - 3*a*hy*y**2 + 3*hx*hy*y**2 - 
     &         2*a*y**3 + 2*hx*y**3)*atan(x/y))/3. - 
     &    ((b - hy)*x**2*(-3*a + 3*hx + 2*x)*atan(y/x))/3. + 
     &    (x*(-4*a*x**2 + 4*hx*x**2 + 3*x**3 + 24*a*b*y - 
     &         24*b*hx*y - 24*a*hy*y + 24*hx*hy*y - 12*b*x*y + 
     &         12*hy*x*y - 12*a*y**2 + 12*hx*y**2 + 6*x*y**2)*
     &       Log(x**2 + y**2))/24. + 
     &    ((-4*b*y**3 + 4*hy*y**3 + 3*y**4)*Log(x**2 + y**2))/24.)/
     &  2.

      xint2(x,y,a,b,hx,hy,e2)=
     &  (((-a - hx)*x**3)/9. + x**4/16. + 
     &    ((a + hx)*x*(18*b - 18*hy - 7*y)*y)/6. + 
     &    (x**2*y*(-28*b + 28*hy + 9*y))/24. + 
     &    ((-3*a*b*y**2 - 3*b*hx*y**2 + 3*a*hy*y**2 + 
     &         3*hx*hy*y**2 + 2*a*y**3 + 2*hx*y**3)*atan(x/y))/3.
     &     + ((b - hy)*x**2*(-3*a - 3*hx + 2*x)*atan(y/x))/3. - 
     &    (x*(-4*a*x**2 - 4*hx*x**2 + 3*x**3 + 24*a*b*y + 
     &         24*b*hx*y - 24*a*hy*y - 24*hx*hy*y - 12*b*x*y + 
     &         12*hy*x*y - 12*a*y**2 - 12*hx*y**2 + 6*x*y**2)*
     &       Log(x**2 + y**2))/24. + 
     &    ((4*b*y**3 - 4*hy*y**3 - 3*y**4)*Log(x**2 + y**2))/24.)/2.

      xint3(x,y,a,b,hx,hy,e2)=
     &  (((-a + hx)*x**3)/9. + x**4/16. + 
     &    ((a - hx)*x*(18*b + 18*hy - 7*y)*y)/6. + 
     &    (x**2*y*(-28*b - 28*hy + 9*y))/24. + 
     &    ((-3*a*b*y**2 + 3*b*hx*y**2 - 3*a*hy*y**2 + 
     &         3*hx*hy*y**2 + 2*a*y**3 - 2*hx*y**3)*atan(x/y))/3.
     &     + ((b + hy)*x**2*(-3*a + 3*hx + 2*x)*atan(y/x))/3. - 
     &    (x*(-4*a*x**2 + 4*hx*x**2 + 3*x**3 + 24*a*b*y - 
     &         24*b*hx*y + 24*a*hy*y - 24*hx*hy*y - 12*b*x*y - 
     &         12*hy*x*y - 12*a*y**2 + 12*hx*y**2 + 6*x*y**2)*
     &       Log(x**2 + y**2))/24. + 
     &    ((4*b*y**3 + 4*hy*y**3 - 3*y**4)*Log(x**2 + y**2))/24.)/2.

      xint4(x,y,a,b,hx,hy,e2)=
     &  (((a + hx)*x**3)/9. - x**4/16. + 
     &    (x**2*(28*b + 28*hy - 9*y)*y)/24. - 
     &    ((a + hx)*x*(18*b + 18*hy - 7*y)*y)/6. + 
     &    ((3*a*b*y**2 + 3*b*hx*y**2 + 3*a*hy*y**2 + 3*hx*hy*y**2 - 
     &         2*a*y**3 - 2*hx*y**3)*atan(x/y))/3. - 
     &    ((b + hy)*x**2*(-3*a - 3*hx + 2*x)*atan(y/x))/3. + 
     &    (x*(-4*a*x**2 - 4*hx*x**2 + 3*x**3 + 24*a*b*y + 
     &         24*b*hx*y + 24*a*hy*y + 24*hx*hy*y - 12*b*x*y - 
     &         12*hy*x*y - 12*a*y**2 - 12*hx*y**2 + 6*x*y**2)*
     &       Log(x**2 + y**2))/24. + 
     &    ((-4*b*y**3 - 4*hy*y**3 + 3*y**4)*Log(x**2 + y**2))/24.)/
     &  2.
!
      xintlimit11=
     &  (-25*hx**2*hy**2 + 8*hx*hy**3*atan(hx/hy) + 
     &    8*hx**3*hy*atan(hy/hx) + hx**4*Log(hx**2) + 
     &    hy**4*Log(hy**2) - hx**4*Log(hx**2 + hy**2) + 
     &    6*hx**2*hy**2*Log(hx**2 + hy**2) - 
     &    hy**4*Log(hx**2 + hy**2))/12.

      xintlimit21=
     &  (-50*hx**2*hy**2 - 16*hx*hy**3*atan(hx/hy) + 
     &    16*hx*hy**3*atan((2*hx)/hy) + 
     &    64*hx**3*hy*atan(hy/(2.*hx)) - 
     &    16*hx**3*hy*atan(hy/hx) - 2*hx**4*Log(hx**2) + 
     &    16*hx**4*Log(4*hx**2) - hy**4*Log(hy**2) + 
     &    2*hx**4*Log(hx**2 + hy**2) - 
     &    12*hx**2*hy**2*Log(hx**2 + hy**2) + 
     &    2*hy**4*Log(hx**2 + hy**2) - 
     &    16*hx**4*Log(4*hx**2 + hy**2) + 
     &    24*hx**2*hy**2*Log(4*hx**2 + hy**2) - 
     &    hy**4*Log(4*hx**2 + hy**2))/24.
!
      xintlimit12=
     &  (-50*hx**2*hy**2 + 64*hx*hy**3*atan(hx/(2.*hy)) - 
     &    16*hx*hy**3*atan(hx/hy) - 16*hx**3*hy*atan(hy/hx) + 
     &    16*hx**3*hy*atan((2*hy)/hx) - hx**4*Log(hx**2) - 
     &    2*hy**4*Log(hy**2) + 16*hy**4*Log(4*hy**2) + 
     &    2*hx**4*Log(hx**2 + hy**2) - 
     &    12*hx**2*hy**2*Log(hx**2 + hy**2) + 
     &    2*hy**4*Log(hx**2 + hy**2) - hx**4*Log(hx**2 + 4*hy**2) + 
     &    24*hx**2*hy**2*Log(hx**2 + 4*hy**2) - 
     &    16*hy**4*Log(hx**2 + 4*hy**2))/24.
!
      xintlimit22=
     &  (-100*hx**2*hy**2 - 128*hx*hy**3*atan(hx/(2.*hy)) + 
     &    160*hx*hy**3*atan(hx/hy) - 
     &    32*hx*hy**3*atan((2*hx)/hy) - 
     &    128*hx**3*hy*atan(hy/(2.*hx)) + 
     &    160*hx**3*hy*atan(hy/hx) - 
     &    32*hx**3*hy*atan((2*hy)/hx) + 2*hx**4*Log(hx**2) - 
     &    16*hx**4*Log(4*hx**2) + 2*hy**4*Log(hy**2) - 
     &    16*hy**4*Log(4*hy**2) - 4*hx**4*Log(hx**2 + hy**2) + 
     &    24*hx**2*hy**2*Log(hx**2 + hy**2) - 
     &    4*hy**4*Log(hx**2 + hy**2) + 
     &    32*hx**4*Log(4*hx**2 + hy**2) - 
     &    48*hx**2*hy**2*Log(4*hx**2 + hy**2) + 
     &    2*hy**4*Log(4*hx**2 + hy**2) + 
     &    2*hx**4*Log(hx**2 + 4*hy**2) - 
     &    48*hx**2*hy**2*Log(hx**2 + 4*hy**2) + 
     &    32*hy**4*Log(hx**2 + 4*hy**2) - 
     &    16*hx**4*Log(4*hx**2 + 4*hy**2) + 
     &    96*hx**2*hy**2*Log(4*hx**2 + 4*hy**2) - 
     &    16*hy**4*Log(4*hx**2 + 4*hy**2))/48.
!
!     eps2=1.d-20
      eps2=0.d0
! compute integrals involving the Green function for Ex:
!       do j=1,ny+1
        do j=1,ny
          y0=hy*(j-1)
          y1=hy*j
!         do i=1,nx+1
          do i=1,nx
            x0=hx*(i-1)
            x1=hx*i
!-----------
      fmp=
     & xint1(x0-hx,y0-hy,x0,y0,hx,hy,eps2)+xint1(x0,y0,x0,y0,hx,hy,eps2)
     &-xint1(x0-hx,y0,x0,y0,hx,hy,eps2)-xint1(x0,y0-hy,x0,y0,hx,hy,eps2)
     &+xint2(x0,y0-hy,x0,y0,hx,hy,eps2)+xint2(x0+hx,y0,x0,y0,hx,hy,eps2)
     &-xint2(x0,y0,x0,y0,hx,hy,eps2)-xint2(x0+hx,y0-hy,x0,y0,hx,hy,eps2)
     &+xint3(x0-hx,y0,x0,y0,hx,hy,eps2)+xint3(x0,y0+hy,x0,y0,hx,hy,eps2)
     &-xint3(x0-hx,y0+hy,x0,y0,hx,hy,eps2)-xint3(x0,y0,x0,y0,hx,hy,eps2)
     &+xint4(x0,y0,x0,y0,hx,hy,eps2)+xint4(x0+hx,y0+hy,x0,y0,hx,hy,eps2)
     &-xint4(x0,y0+hy,x0,y0,hx,hy,eps2)-xint4(x0+hx,y0,x0,y0,hx,hy,eps2)
      if(fmp.ne.fmp)then
      write(12,123)i,j,x0,y0
  123 format(1x,'(greenfx)fmp=NaN at i,j,x0,y0=',2(i5,1x),2(1pe12.5,1x))
      call myflush(12)
      if(i.eq.1.and.j.eq.1)fmp=xintlimit11
      if(i.eq.1.and.j.eq.2)fmp=xintlimit12
      if(i.eq.2.and.j.eq.1)fmp=xintlimit21
      if(i.eq.2.and.j.eq.2)fmp=xintlimit22
      write(12,*)'new limiting value of fmp =',fmp
      call myflush(12)
      endif
      g(i,j)=fmp/(hx*hy)
!-----------
          enddo
        enddo
!reflections:
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j)=g(n1-i+2,j)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j)=g(i,n2-j+2)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j)=g(n1-i+2,n2-j+2)
      enddo
      enddo
! compute FFT of the Green function
      scale=1.0
!     write(6,*)'calling fft2dhpf in greenfx'
      call fft2dhpf(n1,n2,1,0,scale,0,g,gxtr)
!     write(6,*)'returned from fft2dhpf in greenfx'
      return
      end
!
!=============================================================
!
!
      subroutine greenfx2d(gxtr,rho,nx,ny,n1,n2)
! green function routine.
      implicit double precision(a-h,o-z)
      double complex g,gxtr
      dimension g(n1,n2),gxtr(n2,n1)
      dimension rho(nx,ny),exg(nx,ny),fxg(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi

      xint1(x,y,a,b,hx,hy,e2)=
     &  -((-a + hx)*x**2)/4. - x**3/9. + 
     &  (x*y*(-3*b + 3*hy + 2*y))/6. - 
     &  ((-3*b*y**2 + 3*hy*y**2 + 2*y**3)*atan(x/(y+e2)))/6. - 
     &  ((b - hy)*x*(-2*a + 2*hx + x)*atan(y/(x+e2)))/2. + 
     &  (x**2*(-3*a + 3*hx + 2*x)*log(x**2 + y**2 +e2))/12. - 
     &  ((-2*a*b*y + 2*b*hx*y + 2*a*hy*y - 2*hx*hy*y + 
     &       a*y**2 - hx*y**2)*log(x**2 + y**2 +e2))/4.

      xint2(x,y,a,b,hx,hy,e2)=
     &  -((a + hx)*x**2)/4. + x**3/9. - 
     &  (x*y*(-3*b + 3*hy + 2*y))/6. - 
     &  ((3*b*y**2 - 3*hy*y**2 - 2*y**3)*atan(x/(y+e2)))/6. + 
     &  ((b - hy)*x*(-2*a - 2*hx + x)*atan(y/(x+e2)))/2. - 
     &  (x**2*(-3*a - 3*hx + 2*x)*log(x**2 + y**2 +e2))/12. - 
     &  ((2*a*b*y + 2*b*hx*y - 2*a*hy*y - 2*hx*hy*y - 
     &       a*y**2 - hx*y**2)*log(x**2 + y**2 +e2))/4.

      xint3(x,y,a,b,hx,hy,e2)=
     &  ((-a + hx)*x**2)/4. + x**3/9. - 
     &  (x*y*(-3*b - 3*hy + 2*y))/6. + 
     &  ((-3*b*y**2 - 3*hy*y**2 + 2*y**3)*atan(x/(y+e2)))/6. + 
     &  ((b + hy)*x*(-2*a + 2*hx + x)*atan(y/(x+e2)))/2. - 
     &  (x**2*(-3*a + 3*hx + 2*x)*log(x**2 + y**2 +e2))/12. + 
     &  ((-2*a*b*y + 2*b*hx*y - 2*a*hy*y + 2*hx*hy*y + 
     &       a*y**2 - hx*y**2)*log(x**2 + y**2 +e2))/4.

      xint4(x,y,a,b,hx,hy,e2)=
     &  ((a + hx)*x**2)/4. - x**3/9. + 
     &  (x*y*(-3*b - 3*hy + 2*y))/6. + 
     &  ((3*b*y**2 + 3*hy*y**2 - 2*y**3)*atan(x/(y+e2)))/6. - 
     &  ((b + hy)*x*(-2*a - 2*hx + x)*atan(y/(x+e2)))/2. + 
     &  (x**2*(-3*a - 3*hx + 2*x)*log(x**2 + y**2 +e2))/12. + 
     &  ((2*a*b*y + 2*b*hx*y + 2*a*hy*y + 2*hx*hy*y - 
     &       a*y**2 - hx*y**2)*log(x**2 + y**2 +e2))/4.
!
!     write(6,*)'inside greenfx2d w/ hx,hy=',hx,hy
      xintlimit21=
     &  (-2*hy**3*atan(hx/hy) + hy**3*atan((2*hx)/hy) + 
     &    12*hx**2*hy*atan(hy/(2.*hx)) - 
     &    6*hx**2*hy*atan(hy/hx) - hx**3*Log(hx**2) + 
     &    4*hx**3*Log(4*hx**2) + hx**3*Log(hx**2 + hy**2) - 
     &    3*hx*hy**2*Log(hx**2 + hy**2) - 
     &    4*hx**3*Log(4*hx**2 + hy**2) + 
     &    3*hx*hy**2*Log(4*hx**2 + hy**2))/3.
!
      xintlimit22=
     &  (-16*hy**3*atan(hx/(2.*hy)) + 12*hy**3*atan(hx/hy) - 
     &    2*hy**3*atan((2*hx)/hy) - 
     &    24*hx**2*hy*atan(hy/(2.*hx)) + 
     &    36*hx**2*hy*atan(hy/hx) - 
     &    12*hx**2*hy*atan((2*hy)/hx) + hx**3*Log(hx**2) - 
     &    4*hx**3*Log(4*hx**2) - 2*hx**3*Log(hx**2 + hy**2) + 
     &    6*hx*hy**2*Log(hx**2 + hy**2) + 
     &    8*hx**3*Log(4*hx**2 + hy**2) - 
     &    6*hx*hy**2*Log(4*hx**2 + hy**2) + 
     &    hx**3*Log(hx**2 + 4*hy**2) - 
     &    12*hx*hy**2*Log(hx**2 + 4*hy**2) - 
     &    4*hx**3*Log(4*hx**2 + 4*hy**2) + 
     &    12*hx*hy**2*Log(4*hx**2 + 4*hy**2))/6.
!
!     eps2=1.d-20
      eps2=0.d0
! compute integrals involving the Green function for Ex:
!       do j=1,ny+1
        do j=1,ny
          y0=hy*(j-1)
          y1=hy*j
!         do i=1,nx+1
          do i=1,nx
            x0=hx*(i-1)
            x1=hx*i
!-----------
      fmp=
     & xint1(x0-hx,y0-hy,x0,y0,hx,hy,eps2)+xint1(x0,y0,x0,y0,hx,hy,eps2)
     &-xint1(x0-hx,y0,x0,y0,hx,hy,eps2)-xint1(x0,y0-hy,x0,y0,hx,hy,eps2)
     &+xint2(x0,y0-hy,x0,y0,hx,hy,eps2)+xint2(x0+hx,y0,x0,y0,hx,hy,eps2)
     &-xint2(x0,y0,x0,y0,hx,hy,eps2)-xint2(x0+hx,y0-hy,x0,y0,hx,hy,eps2)
     &+xint3(x0-hx,y0,x0,y0,hx,hy,eps2)+xint3(x0,y0+hy,x0,y0,hx,hy,eps2)
     &-xint3(x0-hx,y0+hy,x0,y0,hx,hy,eps2)-xint3(x0,y0,x0,y0,hx,hy,eps2)
     &+xint4(x0,y0,x0,y0,hx,hy,eps2)+xint4(x0+hx,y0+hy,x0,y0,hx,hy,eps2)
     &-xint4(x0,y0+hy,x0,y0,hx,hy,eps2)-xint4(x0+hx,y0,x0,y0,hx,hy,eps2)
      if(fmp.ne.fmp)then
      write(12,123)i,j,x0,y0
  123 format(1x,'(greenfx)fmp=NaN at i,j,x0,y0=',2(i5,1x),2(1pe12.5,1x))
      call myflush(12)
      if(i.eq.1.and.j.eq.1)fmp=0.d0
      if(i.eq.1.and.j.eq.2)fmp=0.d0
      if(i.eq.2.and.j.eq.1)fmp=xintlimit21
      if(i.eq.2.and.j.eq.2)fmp=xintlimit22
      write(12,*)'new limiting value of fmp =',fmp
      call myflush(12)
      endif
      g(i,j)=fmp/(hx*hy)
!-----------
          enddo
        enddo
!reflections:
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,j)
!     g(i,j)=-g(n1-i+1,j)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j)=g(i,n2-j+2)
!     g(i,j)=g(i,n2-j+1)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,n2-j+2)
!     g(i,j)=-g(n1-i+1,n2-j+1)
      enddo
      enddo
! compute FFT of the Green function
      scale=1.0
!     write(6,*)'calling fft2dhpf in greenfx'
      call fft2dhpf(n1,n2,1,0,scale,0,g,gxtr)
!     write(6,*)'returned from fft2dhpf in greenfx'
      return
      end
!
!=============================================================
!
      subroutine greenfy2d(gxtr,rho,nx,ny,n1,n2)
! green function routine.
      implicit double precision(a-h,o-z)
      double complex g,gxtr
      dimension g(n1,n2),gxtr(n2,n1)
      dimension rho(nx,ny),fyg(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      yint1(x,y,a,b,hx,hy,e2)=
     &  -((a - hx)*x*(2*b - 2*hy + y))/2. + 
     &  (x**2*(3*b - 3*hy + 4*y))/12. + 
     &  ((2*a*b*y - 2*b*hx*y - 2*a*hy*y + 2*hx*hy*y - 
     &       a*y**2 + hx*y**2)*atan(x/(y+e2)))/2. - 
     &  (x**2*(-3*a + 3*hx + 2*x)*atan(y/(x+e2)))/6. - 
     &  ((b - hy)*x*(-2*a + 2*hx + x)*log(x**2 + y**2 +e2))/4. + 
     &  ((-3*b*y**2 + 3*hy*y**2 + 2*y**3)*log(x**2 + y**2 +e2))/
     &   12.

      yint2(x,y,a,b,hx,hy,e2)=
     &  (x**2*(-3*b + 3*hy - 4*y))/12. + 
     &  ((a + hx)*x*(2*b - 2*hy + y))/2. + 
     &  ((-2*a*b*y - 2*b*hx*y + 2*a*hy*y + 2*hx*hy*y + 
     &       a*y**2 + hx*y**2)*atan(x/(y+e2)))/2. + 
     &  (x**2*(-3*a - 3*hx + 2*x)*atan(y/(x+e2)))/6. + 
     &  ((b - hy)*x*(-2*a - 2*hx + x)*log(x**2 + y**2 +e2))/4. + 
     &  ((3*b*y**2 - 3*hy*y**2 - 2*y**3)*log(x**2 + y**2 +e2))/
     &   12.

      yint3(x,y,a,b,hx,hy,e2)=
     &  (x**2*(-3*b - 3*hy - 4*y))/12. + 
     &  ((a - hx)*x*(2*b + 2*hy + y))/2. + 
     &  ((-2*a*b*y + 2*b*hx*y - 2*a*hy*y + 2*hx*hy*y + 
     &       a*y**2 - hx*y**2)*atan(x/(y+e2)))/2. + 
     &  (x**2*(-3*a + 3*hx + 2*x)*atan(y/(x+e2)))/6. + 
     &  ((b + hy)*x*(-2*a + 2*hx + x)*log(x**2 + y**2 +e2))/4. + 
     &  ((3*b*y**2 + 3*hy*y**2 - 2*y**3)*log(x**2 + y**2 +e2))/
     &   12.

      yint4(x,y,a,b,hx,hy,e2)=
     &  -((a + hx)*x*(2*b + 2*hy + y))/2. + 
     &  (x**2*(3*b + 3*hy + 4*y))/12. + 
     &  ((2*a*b*y + 2*b*hx*y + 2*a*hy*y + 2*hx*hy*y - 
     &       a*y**2 - hx*y**2)*atan(x/(y+e2)))/2. - 
     &  (x**2*(-3*a - 3*hx + 2*x)*atan(y/(x+e2)))/6. - 
     &  ((b + hy)*x*(-2*a - 2*hx + x)*log(x**2 + y**2 +e2))/4. + 
     &  ((-3*b*y**2 - 3*hy*y**2 + 2*y**3)*log(x**2 + y**2 +e2))/
     &   12.
!
!     write(6,*)'inside greenfy2d'

      yintlimit12=
     &  (12*hx*hy**2*atan(hx/(2.*hy)) - 
     &    6*hx*hy**2*atan(hx/hy) - 2*hx**3*atan(hy/hx) + 
     &    hx**3*atan((2*hy)/hx) - hy**3*Log(hy**2) + 
     &    4*hy**3*Log(4*hy**2) - 3*hx**2*hy*Log(hx**2 + hy**2) + 
     &    hy**3*Log(hx**2 + hy**2) + 
     &    3*hx**2*hy*Log(hx**2 + 4*hy**2) - 
     &    4*hy**3*Log(hx**2 + 4*hy**2))/3.

      yintlimit22=
     &  (-24*hx*hy**2*atan(hx/(2.*hy)) + 
     &    36*hx*hy**2*atan(hx/hy) - 
     &    12*hx*hy**2*atan((2*hx)/hy) - 
     &    16*hx**3*atan(hy/(2.*hx)) + 12*hx**3*atan(hy/hx) - 
     &    2*hx**3*atan((2*hy)/hx) + hy**3*Log(hy**2) - 
     &    4*hy**3*Log(4*hy**2) + 6*hx**2*hy*Log(hx**2 + hy**2) - 
     &    2*hy**3*Log(hx**2 + hy**2) - 
     &    12*hx**2*hy*Log(4*hx**2 + hy**2) + 
     &    hy**3*Log(4*hx**2 + hy**2) - 
     &    6*hx**2*hy*Log(hx**2 + 4*hy**2) + 
     &    8*hy**3*Log(hx**2 + 4*hy**2) + 
     &    12*hx**2*hy*Log(4*hx**2 + 4*hy**2) - 
     &    4*hy**3*Log(4*hx**2 + 4*hy**2))/6.

!     eps2=1.d-20
      eps2=0.d0
! compute integrals involving the Green function for Ey:
        do j=1,ny+1
          y0=hy*(j-1)
          y1=hy*j
          do i=1,nx+1
            x0=hx*(i-1)
            x1=hx*i
!-----------
      fmp=
     & yint1(x0-hx,y0-hy,x0,y0,hx,hy,eps2)+yint1(x0,y0,x0,y0,hx,hy,eps2)
     &-yint1(x0-hx,y0,x0,y0,hx,hy,eps2)-yint1(x0,y0-hy,x0,y0,hx,hy,eps2)
     &+yint2(x0,y0-hy,x0,y0,hx,hy,eps2)+yint2(x0+hx,y0,x0,y0,hx,hy,eps2)
     &-yint2(x0,y0,x0,y0,hx,hy,eps2)-yint2(x0+hx,y0-hy,x0,y0,hx,hy,eps2)
     &+yint3(x0-hx,y0,x0,y0,hx,hy,eps2)+yint3(x0,y0+hy,x0,y0,hx,hy,eps2)
     &-yint3(x0-hx,y0+hy,x0,y0,hx,hy,eps2)-yint3(x0,y0,x0,y0,hx,hy,eps2)
     &+yint4(x0,y0,x0,y0,hx,hy,eps2)+yint4(x0+hx,y0+hy,x0,y0,hx,hy,eps2)
     &-yint4(x0,y0+hy,x0,y0,hx,hy,eps2)-yint4(x0+hx,y0,x0,y0,hx,hy,eps2)
      if(fmp.ne.fmp)then
      write(12,123)i,j,x0,y0
  123 format(1x,'(greenfy)fmp=NaN at i,j,x0,y0=',2(i5,1x),2(1pe12.5,1x))
      call myflush(12)
      if(i.eq.1.and.j.eq.1)fmp=0.d0
      if(i.eq.1.and.j.eq.2)fmp=yintlimit12
      if(i.eq.2.and.j.eq.1)fmp=0.d0
      if(i.eq.2.and.j.eq.2)fmp=yintlimit22
      write(12,*)'new limiting value of fmp =',fmp
      call myflush(12)
      endif
      g(i,j)=fmp/(hx*hy)
!-----------
          enddo
        enddo
!reflections:
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j)=g(n1-i+2,j)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j)=-g(i,n2-j+2)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,n2-j+2)
      enddo
      enddo
! compute FFT of the Green function
      scale=1.0
!     write(6,*)'calling fft2dhpf in greenfy'
      call fft2dhpf(n1,n2,1,0,scale,0,g,gxtr)
!     write(6,*)'returned from fft2dhpf in greenfy'
      return
      end
!
      subroutine greenfxold(gxtr,rho,nx,ny,n1,n2)
! green function routine.
      implicit double precision(a-h,o-z)
      double complex g,gxtr
      dimension g(n1,n2),gxtr(n2,n1)
      dimension rho(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!
      gfunx(x,y,e2)=x/(x**2+y**2+e2)
!
!     write(6,*)'inside greenfxold'
      eps2=1.d-20
! compute integrals involving the Green function for Ex:
        do j=1,ny+1
          y0=hy*(j-1)
          y1=hy*j
          do i=1,nx+1
            x0=hx*(i-1)
            x1=hx*i
!-----------
      fmp=gfunx(x0,y0,eps2)
      g(i,j)=fmp*hx*hy
!-----------
          enddo
        enddo
!reflections:
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,j)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j)=g(i,n2-j+2)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,n2-j+2)
      enddo
      enddo
!
! compute FFT of the Green function
      scale=1.0
!     write(6,*)'calling fft2dhpf in greenfx'
      call fft2dhpf(n1,n2,1,0,scale,0,g,gxtr)
!     write(6,*)'returned from fft2dhpf in greenfx'
      return
      end
!
!=============================================================
!
      subroutine greenfyold(gxtr,rho,nx,ny,n1,n2)
! green function routine.
      implicit double precision(a-h,o-z)
      double complex g,gxtr
      dimension g(n1,n2),gxtr(n2,n1)
      dimension rho(nx,ny),fyg(nx,ny)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
!
      gfuny(x,y,e2)=y/(x**2+y**2+e2)
!
!     write(6,*)'inside greenfyold'
      eps2=1.d-20
! compute integrals involving the Green function for Ey:
        do j=1,ny+1
          y0=hy*(j-1)
          y1=hy*j
          do i=1,nx+1
            x0=hx*(i-1)
            x1=hx*i
!-----------
      fmp=gfuny(x0,y0,eps2)
      g(i,j)=fmp*hx*hy
!-----------
          enddo
        enddo
!reflections:
      do j=1,ny
      do i=1+nx,nx+nx
      g(i,j)=g(n1-i+2,j)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1,nx
      g(i,j)=-g(i,n2-j+2)
      enddo
      enddo
      do j=1+ny,ny+ny
      do i=1+nx,nx+nx
      g(i,j)=-g(n1-i+2,n2-j+2)
      enddo
      enddo
! compute FFT of the Green function
      scale=1.0
!     write(6,*)'calling fft2dhpf in greenfy'
      call fft2dhpf(n1,n2,1,0,scale,0,g,gxtr)
!     write(6,*)'returned from fft2dhpf in greenfy'
      return
      end
