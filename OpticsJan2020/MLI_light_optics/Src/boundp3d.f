      subroutine boundp3d(coord,msk,np,gblam,nx,ny,nz,nadj)
      use rays
      implicit double precision(a-h,o-z)
      logical msk
      dimension coord(6,np),msk(np)
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/gridxtra/xsml,xbig,ysml,ybig,zsml,zbig
      dimension bndy(6),rbndy(6)
!hpf$ distribute coord(*,block)
!hpf$ align (:) with coord(*,:) :: msk
!     if(idproc.eq.0)write(6,*)'inside boundp'
! transverse:
      bndy(1)=maxval(coord(1,:),1,msk)
      bndy(2)=maxval(coord(3,:),1,msk)
      bndy(3)=maxval(coord(5,:),1,msk)
      bndy(4)=-minval(coord(1,:),1,msk)
      bndy(5)=-minval(coord(3,:),1,msk)
      bndy(6)=-minval(coord(5,:),1,msk)
      call MPI_ALLREDUCE(bndy,rbndy,6,mreal,mpimax,lworld,ierror)
      xbig=rbndy(1)
      ybig=rbndy(2)
      zbig=rbndy(3)
      xsml=-rbndy(4)
      ysml=-rbndy(5)
      zsml=-rbndy(6)
!     write(6,*)'xsml,xbig=',xsml,xbig
!     write(6,*)'ysml,ybig=',ysml,ybig
!     write(6,*)'zsml,zbig=',zsml,zbig
! note: this would be a good place to check which particles
! are inside the beampipe, then mask off those that are not.
      eps=0.05
      xmin=xsml-eps*(xbig-xsml)
      xmax=xbig+eps*(xbig-xsml)
      ymin=ysml-eps*(ybig-ysml)
      ymax=ybig+eps*(ybig-ysml)
!     if(idproc.eq.0)then
!       write(6,*)'xmin,xmax=',xmin,xmax
!       write(6,*)'ymin,ymax=',ymin,ymax
!     endif
      hx=(xmax-xmin)/(nx-1)
      hy=(ymax-ymin)/(ny-1)
      hxi=1.d0/hx
      hyi=1.d0/hy
! longitudinal:
      delz=zbig-zsml
!     if(idproc.eq.0)then
!       write(6,*)'delz,gblam=',delz,gblam
!       write(6,*)'nadj=',nadj
!     endif
      if(delz.gt.gblam .and. nadj.eq.1)then
        if(idproc.eq.0)then
        write(6,*)'error: delz.gt.gam*beta*lambda in routine boundp'
        write(6,*)'zbig,zsml=',zbig,zsml
        write(6,*)'delz,gblam=',zbig-zsml,gblam
        endif
        call myexit
      endif
      if(nadj.eq.0)then
        zmin=zsml-eps*(zbig-zsml)
        zmax=zbig+eps*(zbig-zsml)
!       if(idproc.eq.0)write(6,*)'(nadj.eq.0) zmin,zmax=',zmin,zmax
      else
        zmin=zsml-0.5d0*(gblam-delz)
        zmax=zbig+0.5d0*(gblam-delz)
!       if(idproc.eq.0)write(6,*)'(nadj.ne.0) zmin,zmax=',zmin,zmax
      endif
      if(nadj.eq.0)hz=(zmax-zmin)/(nz-1)
      if(nadj.eq.1)hz=gblam/nz
      hzi=1.d0/hz
!     if(idproc.eq.0)write(6,*)'leaving boundp3d with hz=',hz
      return
      end
