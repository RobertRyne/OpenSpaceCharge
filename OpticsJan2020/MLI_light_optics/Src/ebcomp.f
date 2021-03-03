************************************************************************
      function angle(x,y)
c
c  compute the angle defined by the points (x,y), (0,0), and (1,0)
c  returns a value in [0,twopi)
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      use math_consts
      implicit none
      double precision angle,x,y
      double precision t
c
      t=datan2(y,x)
      if(t.lt.0.d0) t=t+twopi
      angle=t
c
      return
      end
c
************************************************************************
      subroutine cegengrad(fnin,ftype,infiles,nfile,kfile,jfile,        &
     &                     zi,zf,nz,f,r,kmx,nk,nprec)
c
c Compute generalized gradients for an rf cavity.
      use parallel, only : idproc
      use math_consts
      use phys_consts
      use e_gengrad
      implicit none
      character(80), intent(in) :: fnin  ! filename or base filename
      character(8), intent(in) :: ftype  ! file type
      integer, intent(in) :: infiles  ! number of input files
      integer, intent(in) :: nfile  ! file unit number for output
      integer, intent(in) :: kfile  ! file unit number for output of e0
      integer, intent(in) :: jfile  ! file unit number for output of efld
      double precision, intent(inout) :: zi, zf  ! initial and final z
      integer, intent(in) :: nz  ! number of z intervals
      double precision, intent(inout) :: f  ! cavity frequency
      double precision, intent(in) :: r  ! extract data at this radius
      double precision, intent(in) :: kmx  ! maximum wave number k
      integer, intent(in) :: nk  ! number of k intervals
      integer, intent(in) :: nprec  ! precision at which to write data
c-----!----------------------------------------------------------------!
      type (Efield_data) efield
      type (charfn) e0
c
      nullify(efield%zvals)
      nullify(efield%Ezdata)
      !nullify(efield%Erdata)
      !nullify(efield%Btdata)
      nullify(e0%e0r)
      nullify(e0%e0i)
c
      if (trim(ftype).eq."t7") then
        call read_efield_t7(fnin,infiles,zi,zf,f,r,efield)
      else if (trim(ftype).eq."ez") then
        call read_efield_ez(fnin,infiles,jfile,zi,zf,f,r,efield)
      else
        if (idproc.eq.0) then
          write(6,*) "<*** ERROR ***> in ebcomp::cegengrad:"
          write(6,*) "  File type ",trim(ftype)," not recognized!"
        end if
        call myexit()
      end if
      call comp_charfn(efield,kmx,nk,kfile,e0)
      call comp_egengrads(nfile,e0,zi,zf,nz)
      close(nfile)
c
      if(associated(efield%zvals)) deallocate(efield%zvals)
      if(associated(efield%Ezdata)) deallocate(efield%Ezdata)
      !if(associated(efield%Erdata)) deallocate(efield%Erdata)
      !if(associated(efield%Btdata)) deallocate(efield%Btdata)
      if(associated(e0%e0r)) deallocate(e0%e0r)
      if(associated(e0%e0i)) deallocate(e0%e0i)
c
      return
      end
c
************************************************************************
      subroutine nlrfcav(zed,zlc,zh,nst,tm,gm,h,mh)
c
c Compute nonlinear map for slice of an rf cavity.
      use parallel, only : idproc
      use math_consts
      use phys_consts
      use beamdata
      use lieaparam, only : monoms
      use e_gengrad
      implicit none
      double precision, intent(in) :: zed  ! z at cavity entrance (edge)
      double precision, intent(in) :: zlc  ! z at slice entrance
      double precision, intent(in) :: zh   ! length of slice
      integer, intent(in) :: nst  ! number of z-steps for this slice
      double precision, intent(inout) :: tm  ! w_scl*t_g at slice
      double precision, intent(inout) :: gm  ! rel. gamma at slice
      ! Lie generators, h, and linear map, mh, for slice
      double precision, dimension(monoms), intent(out) :: h
      double precision, dimension(6,6), intent(out) :: mh
c-----!----------------------------------------------------------------!
      include 'map.inc'
      integer :: i,j
      double precision :: betgam,engin,engfin
c
      if (idproc.eq.0) then
        write (6,*) "ebcomp::nlrfcav"
        write (6,*) "  ...under construction..."
      end if
c
c sanity check
      if (zh.le.0.) then
        if (idproc.eq.0) then
          write(6,*) "<*** ERROR ***> in ebcomp::nlrfcav:"
          write(6,*) "  The slice length zh (=",zh,") must exceed zero!"
        endif
        call myexit()
      endif
c
c initialize map to identity
      do i=1,6
        do j=1,6
          mh(i,j)=0.d0
        end do
        mh(i,i)=1.d0
      end do
      do i=1,monoms
        h(i)=0.d0
      end do
c
c print initial map
      !if (idproc.eq.0) then
      !  write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%'
      !  write(6,*) '==initial transfer map=='
      !  write(6,*) '  matrix part:'
      !  do i=1,6
      !    write(6,*) (mh(i,j),j=1,6)
      !  end do
      !  write(6,*) '  nonlinear generators:'
      !  do i=1,monoms
      !    if (h(i).ne.0.d0) then
      !      write(6,*) 'h(',i,')=',h(i)
      !    end if
      !  end do
      !  write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%'
      !end if
c
      engin=gm*pmass
      if (idproc.eq.0) then
        write(36,*) arclen,engin
        call myflush(36)
      endif
      if (idproc.eq.0) then
        write (6,*) 'calling subroutine mapnlrfcav() ...'
        write (6,*) 'zed,zlc,zh,nst,tm,gm='
        write (6,*)  zed,zlc,zh,nst,tm,gm
      end if
      call mapnlrfcav(zed,zlc,zh,nst,tm,gm,h,mh)
      if (idproc.eq.0) then
        write (6,*) 'returned from subroutine mapnlrfcav() ...'
        write (6,*) 'zed,zlc,zh,nst,tm,gm='
        write (6,*)  zed,zlc,zh,nst,tm,gm
      end if
      ! update beamdata: gamma, gamm1, beta, and brho
      gamma=gm
      gamm1=gamma-1.d0
      betgam=sqrt(gamm1*(gamma+1.d0))
      beta=betgam/gamma
      brho=pmass*betgam/c_light
      engfin=pmass*gamma
      if (idproc.eq.0) then
        write(36,*) arclen,engfin
        write(36,*) ' '
        call myflush(36)
      endif
c
      return
      end
c
************************************************************************
      subroutine mapnlrfcav(zed,zlc,zh,nst,tm,gm,h,mh)
c
c Compute nonlinear map for slice of an rf cavity.
      use parallel, only : idproc
      use math_consts
      use phys_consts
      use beamdata
      use lieaparam, only : monoms
      use e_gengrad
      implicit none
      double precision, intent(in) :: zed  ! z at cavity entrance (edge)
      double precision, intent(in) :: zlc  ! z at slice entrance
      double precision, intent(in) :: zh   ! length of slice
      integer, intent(in) :: nst  ! number of z-steps for this slice
      double precision, intent(inout) :: tm  ! w_scl*t_g at slice
      double precision, intent(inout) :: gm  ! rel. gamma at slice
      ! Lie generators, h, and linear map, mh, for slice
      double precision, dimension(monoms), intent(inout) :: h
      double precision, dimension(6,6), intent(inout) :: mh
c-----!----------------------------------------------------------------!
      !!!! module anyone? !!!!
      integer iflag
      include 'hmflag.inc'
      !!!! ^^^^^^^^^^^^^^ !!!!
      integer, parameter :: neqn=monoms+15
      integer :: i,j
      character(len=6) :: stctd
      double precision :: dz
      double precision, dimension(neqn) :: y
c
      if (idproc.eq.0) then
        write (6,*) "ebcomp::mapnlrfcav"
        write (6,*) "  ...under construction..."
      end if
c
c array y comprises the quantities to integrate
c   y(1-6) = given (design) trajectory
c   y(7-42) = matrix (6x6)
c   y(43-98) = f3 coefficients
c   y(99-224) = f4 coefficients
c   y(225-476) = f5 coefficients
c   y(477-938) = f6 coefficients
      y=0.d0
c initialize design trajectory
c   x = p_x = y = p_y = 0
c   tau = omega_ref t
c   p_tau = -m gamma c^2/(omega_ref sl p_ref)
      y(1:4)=0.d0
      y(5)=tm
      y(6)=-gm*pmass/(omegascl*sl*p0sc)
c initialize rest to identity
      do i=1,6
        y(7*i)=1.d0
      end do
c
c set flag that says "rf cavity" (really!, this approach seems antique)
c the subroutine feval uses this flag to switch to the correct RHSs
c as the differential equations for the map
      iflag=5
c now integrate map equations of motion
      !if (zlc.eq.zed) then
      !  stctd='start '
      !else
      !  stctd='cont  '
      !end if
      stctd='start '
      dz=zh/nst
      if (idproc.eq.0) then
        write (6,*) "calling adam11() with arguments"
        write (6,*) " eggrdata%dz =",eggrdata%dz
        write (6,*) " eggrdata%nz_intrvl =",eggrdata%nz_intrvl
        write (6,*) " dz      =",dz
        write (6,*) " nst     =",nst
        write (6,*) " stctd   =",stctd
        write (6,*) " zlc     =",zlc
        write (6,*) " size(y) =",size(y)
        write (6,*) " neqn    =",neqn
      endif
      !call adam11(eggrdata%dz,eggrdata%nz_intrvl,stctd,zlc,y,neqn)
      call adam11(dz,nst,stctd,zlc,y,neqn)
      write (6,*) "...returned from adam11()"
      call errchk(y(7))
      call putmap(y,h,mh)
      tm=y(5)
      gm=-y(6)*(omegascl*sl*p0sc)/pmass
c
      return
      end
c
************************************************************************
      subroutine nlrfvecpot(z,y)
c
c Use generalized gradients in eggr to compute vector potential at z.
c The vector potential is multiplied by rf_escale*(e/p_ref).
      use parallel, only : idproc
      use beamdata, only : sl,p0sc,omegascl,pmass
      use math_consts
      use phys_consts
      use e_gengrad
      use curve_fit
      implicit none
      double precision, intent(in) :: z  ! determine vec. pot. here
      double precision, intent(in) :: y(:)  ! y(1:6)=ref. ptcl. (scaled)
c-----!----------------------------------------------------------------!
      include 'map.inc'
      integer, save :: iz=0
      integer :: i,j
      double precision, parameter :: tiny=2.0d-16
      double precision :: wl,wr,mwr2,phiz,cphiz,sphiz,cfwrl,sfwl
      double precision :: sclfac
      double precision :: sl2,sl3,sl4,sl5,sl6
      double precision, dimension(:,:), allocatable :: Crc,Czc
c
      !if (idproc.eq.0) then
      !  write (6,*) "(ebcomp)::nlrfvecpot(",z,")"
      !  write (6,*) "  ...under construction..."
      !end if

      ! allocate temporary arrays
      allocate(Crc(0:eggrdata%maxj,0:eggrdata%maxm))
      allocate(Czc(0:eggrdata%maxj,0:eggrdata%maxm))

c find index iz ofzvals nearest z (zvals has index origin 0)
      call cf_searchNrst(eggrdata%zvals,z,iz,iorigin=0)
      if (dabs(eggrdata%zvals(iz)-z).le.tiny) then
        Czc(0,0)=eggrdata%G3c(iz,0,0)
        do j=1,eggrdata%maxj
          Crc(j,0)=eggrdata%G1c(iz,j,0)
          Czc(j,0)=eggrdata%G3c(iz,j,0)
        end do
      else
        call cf_polyInterp3(eggrdata%zvals,eggrdata%G3c(:,0,0),
     &                      z,Czc(0,0),iz,iorigin=0)
        do j=1,eggrdata%maxj
          call cf_polyInterp3(eggrdata%zvals,eggrdata%G1c(:,j,0),
     &                        z,Crc(j,0),iz,iorigin=0)
          call cf_polyInterp3(eggrdata%zvals,eggrdata%G3c(:,j,0),
     &                        z,Czc(j,0),iz,iorigin=0)
        end do
      end if
c compute phase phi(z) and other factors
      wl=eggrdata%angfrq
      wr=wl/omegascl
      mwr2=-wr**2
      phiz=wr*y(5)+rf_phase
      cphiz=dcos(phiz); cfwrl=cphiz*wr/wl
      sphiz=dsin(phiz); sfwl=sphiz/wl
      sl2=sl**2; sl3=sl2*sl; sl4=sl3*sl; sl5=sl4*sl; sl6=sl5*sl
cryneabell:09.Nov.05
      write(44,*) z,phiz,y(5:6),y(7:12)
      write(45,*) z,rf_escale*Czc(0,0),-rf_escale*Crc(1,0)
cryneabell----------
c A_x
      vp%Ax(1)   = -Crc(1,0)*sfwl*sl/2.d0
      vp%Ax(11)  = -Crc(1,0)*cfwrl*sl/2.d0
      vp%Ax(28)  = -Crc(2,0)*sfwl*sl3/32.d0
      vp%Ax(39)  = vp%Ax(28)
      vp%Ax(46)  = mwr2*vp%Ax(1)/2.d0
      vp%Ax(88)  = -Crc(2,0)*cfwrl*sl3/32.d0
      vp%Ax(122) = vp%Ax(88)
      vp%Ax(136) = mwr2*vp%Ax(11)/6.d0
      vp%Ax(210) = -Crc(3,0)*sfwl*sl5/1152.d0
      vp%Ax(221) = 2.d0*vp%Ax(210)
      vp%Ax(228) = mwr2*vp%Ax(39)/2.d0
      vp%Ax(301) = vp%Ax(210)
      vp%Ax(308) = vp%Ax(228)
      vp%Ax(331) = mwr2*vp%Ax(46)/12.d0
      vp%Ax(466) = -Crc(3,0)*cfwrl*sl5/1152.d0
      vp%Ax(500) = 2.d0*vp%Ax(466)
      vp%Ax(514) = mwr2*vp%Ax(122)/6.d0
      vp%Ax(660) = vp%Ax(466)
      vp%Ax(674) = vp%Ax(514)
      vp%Ax(708) = mwr2*vp%Ax(136)/20.d0
c A_y
      vp%Ay(3)   = vp%Ax(1)
      vp%Ay(20)  = vp%Ax(11)
      vp%Ay(30)  = vp%Ax(28)
      vp%Ay(64)  = vp%Ax(39)
      vp%Ay(71)  = vp%Ax(46)
      vp%Ay(97)  = vp%Ax(88)
      vp%Ay(177) = vp%Ax(122)
      vp%Ay(191) = vp%Ax(136)
      vp%Ay(212) = vp%Ax(210)
      vp%Ay(246) = vp%Ax(221)
      vp%Ay(253) = vp%Ax(228)
      vp%Ay(406) = vp%Ax(301)
      vp%Ay(413) = vp%Ax(308)
      vp%Ay(436) = vp%Ax(331)
      vp%Ay(475) = vp%Ax(466)
      vp%Ay(555) = vp%Ax(500)
      vp%Ay(569) = vp%Ax(514)
      vp%Ay(842) = vp%Ax(660)
      vp%Ay(856) = vp%Ax(674)
      vp%Ay(890) = vp%Ax(708)
c A_z
      vp%Az(0)   = -Czc(0,0)*sfwl
      vp%Az(5)   = -Czc(0,0)*cfwrl
      vp%Az(7)   = -Czc(1,0)*sfwl*sl2/4.d0
      vp%Az(18)  = vp%Az(7)
      vp%Az(25)  = mwr2*vp%Az(0)/2.d0
      vp%Az(32)  = -Czc(1,0)*cfwrl*sl2/4.d0
      vp%Az(66)  = vp%Az(32)
      vp%Az(80)  = mwr2*vp%Az(5)/6.d0
      vp%Az(84)  = -Czc(2,0)*sfwl*sl4/64.d0
      vp%Az(95)  = 2.d0*vp%Az(84)
      vp%Az(102) = mwr2*vp%Az(18)/2.d0
      vp%Az(175) = vp%Az(84)
      vp%Az(182) = vp%Az(102)
      vp%Az(205) = mwr2*vp%Az(25)/12.d0
      vp%Az(214) = -Czc(2,0)*cfwrl*sl4/64.d0
      vp%Az(248) = 2.d0*vp%Az(214)
      vp%Az(262) = mwr2*vp%Az(66)/6.d0
      vp%Az(408) = vp%Az(214)
      vp%Az(422) = vp%Az(262)
      vp%Az(456) = mwr2*vp%Az(80)/20.d0
      vp%Az(462) = -Czc(3,0)*sfwl*sl6/2304.d0
      vp%Az(473) = 3.d0*vp%Az(462)
      vp%Az(480) = mwr2*vp%Az(175)/2.d0
      vp%Az(553) = vp%Az(473)
      vp%Az(560) = 2.d0*vp%Az(480)
      vp%Az(583) = mwr2*vp%Az(182)/12.d0
      vp%Az(840) = vp%Az(462)
      vp%Az(847) = vp%Az(480)
      vp%Az(870) = vp%Az(583)
      vp%Az(917) = mwr2*vp%Az(205)/30.d0
c scale vector potential by rf_escale/(p_ref/e) = rf_escale*c/(mc^2/e)
      sclfac=rf_escale*c_light/pmass
      vp%Ax=sclfac*vp%Ax
      vp%Ay=sclfac*vp%Ay
      vp%Az=sclfac*vp%Az
c
      !if (idproc.eq.0) then
      !  write(6,*) "nlrf vector potential components"
      !  write(6,*) "  =scale factors="
      !  write(6,*) " sclfac =",sclfac
      !  write(6,*) " wl, wr =",wl,wr
      !  write(6,*) " phiz =",phiz
      !  write(6,*) " Crc0j(z) =",(Crc(j,0),j=1,3)
      !  write(6,*) " Czc0j(z) =",(Czc(j,0),j=0,3)
      !  write(6,*) "  =x="
      !  do i=0,monoms
      !    if (vp%Ax(i).ne.0.d0) then
      !      write(6,*) 'Ax(',i,')=',vp%Ax(i)
      !    end if
      !  end do
      !  write(6,*) "  =y="
      !  do i=0,monoms
      !    if (vp%Ay(i).ne.0.d0) then
      !      write(6,*) 'Ay(',i,')=',vp%Ay(i)
      !    end if
      !  end do
      !  write(6,*) "  =z="
      !  do i=0,monoms
      !    if (vp%Az(i).ne.0.d0) then
      !      write(6,*) 'Az(',i,')=',vp%Az(i)
      !    end if
      !  end do
      !end if
c
      ! deallocate temporary arrays
      deallocate(Crc)
      deallocate(Czc)
c
      return
      end
c
!***********************************************************************
      subroutine hmltnRF(s,y,h)
! Given longitudinal position s and phase-space coordinates y(1:6) of
! the reference particle, compute an expansion of the Hamiltonian h for
! an RF cavity.
! Based on the approach used in hmltn3 (by F. Neri): compute the
! hamiltonian using a polynomial expansion of the square root.
! (Later, we can use this to test a hard-wired version.)
      use math_consts
      use phys_consts
      use lieaparam, only : monoms
      use beamdata
      use e_gengrad
      use parallel, only : idproc
      implicit none
      double precision, intent(in) :: s  ! longitudinal coordinate
      double precision, intent(inout) :: y(:)  ! y(1:6)=ref. ptcl.
      double precision, intent(out) :: h(:)  ! rf cavity Hamiltonian
!-----!----------------------------------------------------------------!
************************************************************************
      interface
       subroutine nlrfvecpot(s,y)
       implicit none
       double precision, intent(in) :: s  ! longitudinal coordinate
       double precision, intent(in) :: y(:)  ! y(1:6)=ref. ptcl. data
       end subroutine nlrfvecpot
      end interface
************************************************************************
      include 'expon.inc'
      integer, parameter :: maxord=6
      integer :: i,j
      double precision :: bg,gammag,wlbyc
      double precision, dimension(0:maxord) :: a
      double precision, dimension(0:monoms) :: pkx1,pky1,pkx2,pky2
      double precision, dimension(0:monoms) :: x,yy
c
      !if (idproc.eq.0) then
      !  print '(a,1pe16.9,a)',"(ebcomp)::hmltnRF[ s=",s,"]="
      !  do j=1,monoms
      !    if(h(j).ne.0.d0) then
      !      print '(2x,a,i3,a,1pe16.9,2x,3(2i1,1x))',                   &
    ! &!            'h(',j,')=',h(j),(expon(i,j),i=1,6)
      !    end if
      !  end do
      !  do j=1,monoms
      !    if(y(j).ne.0.d0) then
      !      print '(2x,a,i3,a,1pe16.9,2x,3(2i1,1x))',                   &
    ! &!            'y(',j,')=',y(j),(expon(i,j),i=1,6)
      !    end if
      !  end do
      !end if
c
c coefficients in expansion of 1-sqrt(1+x)
      a(0) = 0.d0
      a(1) = -1.d0/2.d0
      do j=2,maxord
        a(j)=a(j-1)*(1.5d0/j-1.d0)
      end do
c
c factors
      wlbyc=omegascl*sl/c_light
      gammag=-wlbyc*y(6)
      bg=sqrt((gammag+1.d0)*(gammag-1.d0))
c
c compute rf cavity vector potential
      !write(6,*) "calling nlrfvecpot ..."
      call nlrfvecpot(s,y)
      !write(6,*) "returned from nlrfvecpot ..."
c
c construct Hamiltonian
      pkx1=0.d0; pkx1(2)=1.d0; pkx1=pkx1-vp%Ax
      pky1=0.d0; pky1(4)=1.d0; pky1=pky1-vp%Ay
      call pmult(pkx1,pkx1,pkx2,maxord)
      call pmult(pky1,pky1,pky2,maxord)
      x=0.d0; x(6)=-2.d0*gammag*wlbyc; x(27)=wlbyc**2
      x=(x-pkx2-pky2)/(bg**2)
c yy = 1-sqrt(1+x)
      call poly1(maxord,a,x,yy,maxord)
c h = {bg*[1-sqrt(1+x)]-Az}/sl, but don't include constant
c or linear terms because we've already removed them
      h=0.d0
      do j=7,monoms
        h(j)=(bg*yy(j)-vp%Az(j))/sl
      end do
c
      !if (idproc.eq.0) then
      !  print '(a,1pe16.9,a)',"(ebcomp)::hmltnRF[ s=",s,"]="
      !  do j=1,monoms
      !    if(h(j).ne.0.d0) then
      !      print '(2x,a,i3,a,1pe16.9,2x,3(2i1,1x))',                   &
    ! &!            'h(',j,')=',h(j),(expon(i,j),i=1,6)
      !    end if
      !  end do
      !end if
c
      return
      end subroutine hmltnRF
c
************************************************************************
      function finterp(x,xx,a,b,c,n)
c
c  This subroutine interpolates at x the value of a function described
c  by the n-element arrays a, b, and c.  These arrays hold the quadratic
c  fit parameter determined at the n locations xx by subroutine parfit.
c  NB: The values in xx MUST increase (or decrease) monotonically.
c-----!----------------------------------------------------------------!
      implicit none
      double precision finterp
      integer n,k
      double precision a(n),b(n),c(n),xx(n),x
c
      call locate(xx,n,x,k)
      if(k.lt.1) then
        k=1
      else if(k.gt.(n-1))  then
        k=n-1
      endif
      finterp=a(k)+x*(b(k)+x*c(k))
c
      return
      end
c
************************************************************************
      function finterA(ang)
c
c  This subroutine interpolates in the azimuthal direction a field value
c  for some component of E or B.  It first uses 'locate' to find where
c  ang lies in the array phin; it then uses coefficients determined by
c  'parfit' to interpolate a value for that field component.
c  NB: The angles in phin must increase (or decrease) monotonically.
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      include 'ebdata.inc'
c
      call locate(phin,maxna,ang,n)
      if(n.lt.1) then
        n=1
      else if(n.gt.(maxna-1))  then
        n=maxna-1
      endif
      finterA=aia(n)+ang*(bia(n)+ang*cia(n))
c
      return
      end
c
************************************************************************
        function finterZ(z)
c
c  This subroutine interpolates in the longitudinal direction a field
c  value for some component of E or B.  It first uses 'locate' to find
c  where z lies in the array zn; it then uses coefficients determined by
c 'parfit' to interpolate a value for that field component.
c  NB: The z's in zn must increase (or decrease) monotonically.
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      include 'ebdata.inc'
c
      call locate(zn,maxnz,z,n)
      if(n.lt.1) then
        n=1
      else if(n.gt.(maxnz-1))  then
        n=maxnz-1
      endif
      finterZ=aiz(n)+z*(biz(n)+z*ciz(n))
c
      return
      end
c
************************************************************************
      subroutine filonint(xmin,xmax,y,nsteps,fname,sinint,cosint)
c
c  Generic Filon integrator:
c  This subroutine uses Filon's method to evaluate the integrals from
c  xmin to xmax of
c
c            fname(x) sin y*x dx
c       and  fname(x) cos y*x dx
c
c  Filon's method allows one to evaluate the integral of an oscillating
c  function without having to follow every wiggle with many evaluation
c  points.  Here the number of function evaluations equals 2*nsteps+3.
c  This version has only a single frequency, y.
c  (See F.B.Hildebrand, Introduction to Numerical Analysis, 2nd ed.,
c  McGraw-Hill, 1974 (Dover reprint, 1987), Sect. 3.10.  See, also,
c  Abramowitz & Stegun, p.890.)  For y=0, the cosine integral reduces
c  to the extended form of Simpson's rule.
c
c  Written by Peter Walstrom.
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      external fname
c
      parameter(half=0.5d0,one=1.d0,zero=0.d0)
      dimension dsinint(3),dcosint(3)
c
c  compute step-size
      h=(xmax-xmin)*half/dfloat(nsteps)
c  get Filon weights
      theta=h*y
      call filon_wts(theta,alf,bet,gam)
c      write(20,200) theta,alf,bet,gam
c  200 format(1x,4(1pd11.4,1x))
c
c  initialize intermediate arrays to zero
      do 111 ievod=1,3
        dsinint(ievod)=zero
  111   dcosint(ievod)=zero
c
c  perform Filon integration:
c
c  step through odd, then even, x values
      do 2 ievod=1,2
c  ievod=1 -->  odd points
c  ievod=2 --> even points
c  ievod=3 -->  end points
        nstp=nsteps
c  extra x value in set of even points
        if(ievod.eq.2) nstp=nsteps+1
        do 2 n=1,nstp
          if(ievod.eq.2) go to 3
c  odd points
          x=h*dfloat(2*n-1)+xmin
          wt=one
          go to 4
c  even points
    3     wt=one
          if(n.eq.1) wt=half
          if(n.eq.nstp) wt=half
          x=h*dfloat(2*n-2)+xmin
    4     continue
          f=wt*fname(x)
          xy=x*y
          cosxy=dcos(xy)
          sinxy=dsin(xy)
          dsinint(ievod)=dsinint(ievod)+f*sinxy
    2     dcosint(ievod)=dcosint(ievod)+f*cosxy
c  end points (upper end first)
      x=xmax
      wt=one
      do 5 iend=1,2
        xy=y*x
        cosxy=dcos(xy)
        sinxy=dsin(xy)
        f=wt*fname(x)
        dsinint(3)=dsinint(3)-f*cosxy
        dcosint(3)=dcosint(3)+f*sinxy
        x=xmin
    5   wt=-one
c  now sum the Filon-weighted contributions from ievod=1,2,3
      sinint=h*(alf*dsinint(3)+bet*dsinint(2)+gam*dsinint(1))
      cosint=h*(alf*dcosint(3)+bet*dcosint(2)+gam*dcosint(1))
c
      return
      end
c
************************************************************************
      subroutine filonarr(xn,fn,y,npts,sinint,cosint)
c
c  Generic Filon integrator---array version:
c  This subroutine uses Filon's method to evaluate the integrals from
c  xn(1) to xn(npts) of
c
c            fn(x) sin y*x dx
c       and  fn(x) cos y*x dx
c
c  where fn is also an array of length npts.  Filon's method allows one
c  to evaluate the integral of an oscillating function without having to
c  follow every wiggle with many evaluation points.  This version has
c  only a single frequency, y.  (See F.B.Hildebrand, Introduction to
c  Numerical Analysis, 2nd ed., McGraw-Hill, 1974 (Dover reprint, 1987),
c  Sect. 3.10.  See, also, Abramowitz & Stegun, p.890.)  For y=0, the
c  cosine integral reduces to the extended form of Simpson's rule.
c
c  NB: The input array xn must have equally spaced points.  Moreover, it
c  should contain an _odd_ number of points.  (If xn contains an even
c  number of points, the rightmost point is ignored.)  Also note that
c  the usual definitions of Filon's formulas use index origin zero, so
c  that the evaluation points are {x_0, x_1, ..., x_2n}.  But the
c  implementation below uses index origin one, so that the "even" points
c  are those indexed 1, 3, 5, ..., while the "odd" points are those
c  indexed 2, 4, 6, ....  Sigh, ...!
c
c  Written by Peter Walstrom.
c-----!----------------------------------------------------------------!
      include 'impli.inc'
c
c  calling arrays
      parameter(maxpts=4001)
      dimension xn(maxpts),fn(maxpts)
c
      parameter(half=0.5d0,one=1.d0,zero=0.d0)
      dimension dsinint(3),dcosint(3)
c
c
c  check parity of npts; if even, reduce it by 1
      npoints=npts
      if(2*(npts/2).eq.npts) then
        write(6,*) 'WARNING from subroutine filonarr():'
        write(6,*) '  must have an odd number of data points;'
        write(6,*) '  ignoring last point!'
        npoints=npts-1
      endif
      nhalf=npoints/2
c  note step-size
      h=xn(2)-xn(1)
c  get Filon weights
      theta=h*y
      call filon_wts(theta,alf,bet,gam)
c      write(20,200) theta,alf,bet,gam
c  200 format(1x,4(1pd11.4,1x))
c
c  initialize intermediate arrays to zero
      do 111 ievod=1,3
        dsinint(ievod)=zero
  111   dcosint(ievod)=zero
c
c  perform Filon integration:
c
c  step through odd, then even, x values
      do 2 ievod=1,2
c  ievod=1 -->  odd points
c  ievod=2 --> even points
c  ievod=3 -->  end points
        nstp=nhalf
c  extra x value in set of even points
        if(ievod.eq.2) nstp=nhalf+1
        do 2 n=1,nstp
          if(ievod.eq.2) go to 3
c  "odd" points (n=2,4,6,...,npts-1)
          x=xn(2*n)
          f=fn(2*n)
          wt=one
          go to 4
c  "even" points (n=1,3,5, ... npts)
    3     wt=one
          if(n.eq.1) wt=half
          if(n.eq.nstp) wt=half
          x=xn(2*n-1)
          f=fn(2*n-1)
    4     continue
          f=wt*f
          xy=x*y
          cosxy=dcos(xy)
          sinxy=dsin(xy)
          dsinint(ievod)=dsinint(ievod)+f*sinxy
    2     dcosint(ievod)=dcosint(ievod)+f*cosxy
c  end points (upper end first)
      n=npoints
      wt=one
      do 5 iend=1,2
        x=xn(n)
        xy=y*x
        cosxy=dcos(xy)
        sinxy=dsin(xy)
        f=wt*fn(n)
        dsinint(3)=dsinint(3)-f*cosxy
        dcosint(3)=dcosint(3)+f*sinxy
        n=1
    5   wt=-one
c  now sum the Filon-weighted contributions from ievod=1,2,3
      sinint=h*(alf*dsinint(3)+bet*dsinint(2)+gam*dsinint(1))
      cosint=h*(alf*dcosint(3)+bet*dcosint(2)+gam*dcosint(1))
c
      return
      end
c
************************************************************************
      subroutine filonarr_io(xn,fn,k,io,npts,sinint,cosint)
c
c  Generic Filon integrator---array version with variable index origin:
c  This subroutine uses Filon's method to evaluate the integrals from
c  xn(io) to xn(io+npts-1) of
c
c            fn(x) sin k*x dx
c       and  fn(x) cos k*x dx
c
c  where fn is also an array of length npts with the same index origin
C  io.  Filon's method allows one to evaluate the integral of an
c  oscillating function without having to follow every wiggle with many
c  evaluation points.  This version has only a single frequency, k.
c  (See F.B.Hildebrand, Introduction to Numerical Analysis, 2nd ed.,
c  McGraw-Hill, 1974 (Dover reprint, 1987), Sect. 3.10.  See, also,
c  Abramowitz & Stegun, p.890.)  For k=0, the cosine integral reduces to
c  the extended form of Simpson's rule.
c
c  NB: The input array xn must have equally spaced points.  Moreover, it
c  should contain an _odd_ number of points---or an even number of
c  intervals.  (If xn contains an even number of points, the rightmost
c  point is ignored.)  Also note that the usual definitions of Filon's
c  formulas use index origin zero, so that the evaluation points are
c  {x_0, x_1, ..., x_2n}.  But many Fortran arrays have index origin
c  one, which means the "even" points are those indexed 1, 3, 5, ...,
c  while the "odd" points are those indexed 2, 4, 6, ....  Sigh, ...!
c  The implementation given here allows for an arbitrary index origin.
c
c  Written by Peter Walstrom.
c  Jan.2005 (DTA): Modified to allow for arbitrary index origin, and
c                  changed to implicit none.
c-----!----------------------------------------------------------------!
      implicit none
c
c  arguments
      double precision, dimension(:), intent(in) :: xn,fn
      double precision, intent(in) :: k
      integer, intent(in) :: npts,io
      double precision, intent(out) :: sinint,cosint
c
c  local variables
      double precision, parameter :: half=0.5d0,one=1.d0,zero=0.d0
      integer :: ii,ix,istart,iend,ievod,npoints
      double precision :: h,theta,alf,bet,gam
      double precision :: kx,coskx,sinkx,f,wt
      double precision, dimension(3) :: dsinint,dcosint
c
c  check parity of npts; if even, emit warning and reduce it by 1
      npoints=npts
      if(2*(npts/2).eq.npts) then
        write(6,*) '<*** WARNING ***> from subroutine filonarr():'
        write(6,*) '  must have an odd number of data points;'
        write(6,*) '  ignoring last point!'
        npoints=npts-1
      endif
c
c  set array endpoints
      istart=io
      iend=io+npoints-1
c  note step-size and get Filon weights
      h=xn(istart+1)-xn(istart)
      theta=h*k
      call filon_wts(theta,alf,bet,gam)
c      write(*,200) "filon weights:",theta,alf,bet,gam
c  200 format(1x,4(1pd11.4,1x))
c
c  initialize intermediate arrays to zero
      do ievod=1,3
        dsinint(ievod)=zero
        dcosint(ievod)=zero
      end do
c
c  perform Filon integration:
c
c  ievod=1 -->  odd points
c  ievod=2 --> even points
c  ievod=3 -->  end points
      ievod=2
      do ix=istart,iend
        wt=one
        if(ix.eq.istart.or.ix.eq.iend) wt=half
        kx=k*xn(ix)
        f=wt*fn(ix)
        sinkx=dsin(kx)
        coskx=dcos(kx)
        dsinint(ievod)=dsinint(ievod)+f*sinkx
        dcosint(ievod)=dcosint(ievod)+f*coskx
        if(ievod.eq.2) then
          ievod=1
        else
          ievod=2
        endif
      end do
c  end points (upper end first)
      do ii=1,2
        if(ii.eq.1) then
          ix=iend
          wt=one
        else
          ix=istart
          wt=-one
        endif
        kx=k*xn(ix)
        f=wt*fn(ix)
        dsinint(3)=dsinint(3)-f*dcos(kx)
        dcosint(3)=dcosint(3)+f*dsin(kx)
      end do
c  now sum the Filon-weighted contributions from ievod=1,2,3
      sinint=h*(alf*dsinint(3)+bet*dsinint(2)+gam*dsinint(1))
      cosint=h*(alf*dcosint(3)+bet*dcosint(2)+gam*dcosint(1))
c
      return
      end
c
************************************************************************
      subroutine filon_wts(theta,alpha,beta,gamma)
c
c  This subroutine computes the weights for Filon's integration method.
c-----!----------------------------------------------------------------!
      implicit none
c
c  arguments
      double precision :: theta
      double precision :: alpha,beta,gamma
c
c  local variables
      double precision, parameter :: a1=2.d0/4.5d1,a2=2.d0/3.15d2,      &
     &                               a3=2.d0/4.725d3,a4=8.d0/4.67775d5
      double precision, parameter :: b0=2.d0/3.d0,b1=2.d0/1.5d1,        &
     &                               b2=4.d0/1.05d2,b3=2.d0/5.67d2,     &
     &                               b4=4.d0/2.2275d4,b5=4.d0/6.75675d5
      double precision, parameter :: c0=4.d0/3.0d0,c1=2.d0/1.5d1,       &
     &                               c2=1.d0/2.1d2,c3=1.d0/1.134d4,     &
     &                               c4=1.d0/9.9792d5,c5=1.d0/1.297296d8
      double precision, parameter :: one=1.d0,two=2.d0,smal1=1.d-1
      double precision :: t2,t3
      double precision :: onoth,tuoth2,costh,sinc,cs
c
c  small-angle approximation
      if(dabs(theta).gt.smal1) go to 11
      t2=theta**2
      t3=theta*t2
      alpha=t3*(a1-t2*(a2-t2*(a3-t2*a4)));
      beta=b0+t2*(b1-t2*(b2-t2*(b3-t2*(b4-t2*b5))));
      gamma=c0-t2*(c1-t2*(c2-t2*(c3-t2*(c4-t2*c5))));
      return
c
c  full computation
 11   onoth=one/theta
      tuoth2=two*onoth*onoth
      costh=dcos(theta)
      sinc=dsin(theta)*onoth
      cs=costh-two*sinc
      alpha=onoth*(one+sinc*cs)
      beta=tuoth2*(one+costh*cs)
      gamma=two*tuoth2*(sinc-costh)
      return
c
      end
c
************************************************************************
      subroutine intsimpodd(n,f,out)
c
c Use Simpson's three-point rule to integrate a function whose values at
c n equally-spaced locations are given in array f; put result in "out".
c [Numerical Recipes (1992), p.128]
c Notes: 1) The array f must have an odd number of entries.
c        2) To obtain the desired integral, one MUST, after calling
c           intsimpodd, multiply the result "out" by the stepsize.  MV
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      include 'impli.inc'
      integer n
      double precision f(n),out
c
      sum1=0
      sum2=0
      nhalf=n/2
c
      sum2=sum2+f(2)
      do j=2,nhalf
        sum2=sum2 + f(2*j)
        sum1=sum1 + f(2*j-1)
      enddo
      out = (f(1) + 2.d0*sum1 + 4.d0*sum2 + f(n))/3.d0
c
      return
      end
c
************************************************************************
      subroutine intsimp(n,f,out)
c
c Use Simpson's rule to integrate a function whose values at n equally-
c spaced locations are given in the array f; put result in "out".
c [Numerical Recipes (1992), p.128]
c
c For odd n, use the extended three-point rule (subroutine intsimpodd).
c For even n, apply the three-point rule to the interval i=1,2,...,n-3;
c apply the four-point rule to the remaining interval, i=n-3,...,n; and
c then add the results.  Overall error ~ h^{-4}.  MV
c
c NB: To obtain the desired integral, one MUST, after calling intsimp,
c     multiply the result "out" by the stepsize.
c---x-!--1----x----2----x----3----x----4----x----5----x----6----x----7-!
      include 'impli.inc'
      integer n
      double precision f(n),out
c
      if(n.eq.(2*(n/2)+1)) then
        call intsimpodd(n,f,out)
      else
        call intsimpodd(n-3,f,out)
        xtra=3.d0*(f(n-3) + 3.d0*(f(n-2) + f(n-1)) + f(n))/8.d0
        out=out+xtra
       endif
c
       return
       end
c
************************************************************************
      subroutine locate(xx,n,x,j)
c
c  This subroutine takes an array xx(1..n) and a value x, and it returns
c  the index j such that x lies between xx(j) and x(j+1).  If the value
c  returned is either j=0 or j=n, then x lies outside the range.
c  NB: The array xx must be monotonic, either increasing or decreasing.
c  From Press, et al., Numerical Recipes, 2nd ed. (1992), p.111.
c-----!----------------------------------------------------------------!
      integer j,n
      double precision x,xx(n)
c
c local variables
      integer jl,jm,ju
c
      jl=0
      ju=n+1
c
 10   if(ju-jl.gt.1) then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        goto 10
      endif
      j=jl
c
      return
      end
c
************************************************************************
cryne subroutine name changed to parfit_eb to avoid conflict
cryne with parfit in magnet.f ; fix later. R.D. Ryne June 16, 2004
      subroutine parfit_eb(npoint,xi,yi,ai,bi,ci)
c
c  This subroutine takes data pairs (x(i),y(i)), i=1..npoint and fits a
c  parabola of the form a+b*x+c*x**2 to each triple of adjacent points.
c  It returns the lists of parabola coefficients a(i), b(i), and c(i).
c  Except at the ends, the parabola coefficients a(i), b(i), c(i) are
c  the average of those calculated using the points (i-1,i,i+1) and
c  those calculated using the points (i,i+1,i+2).
c-----!----------------------------------------------------------------!
      include 'impli.inc'
c
c  calling arrays
      dimension xi(npoint),yi(npoint),ai(npoint),bi(npoint),ci(npoint)
c
c  local variables
      npm1=npoint-1
      npm2=npoint-2
c
c  initialize the fitting
      x1=xi(1)
      x2=xi(2)
      x3=xi(3)
      y1=yi(1)
      y2=yi(2)
      y3=yi(3)
      h1=x2-x1
      h2=x3-x2
      h3=x3-x1
      d1= y1/(h1*h3)
      d2=-y2/(h1*h2)
      d3= y3/(h2*h3)
c  and fit a parabola to the first three points
      a1=  d1*x2*x3   + d2*x3*x1   + d3*x1*x2
      b1= -d1*(x2+x3) - d2*(x3+x1) - d3*(x1+x2)
      c1=  d1         + d2         + d3
      ai(1)=a1
      bi(1)=b1
      ci(1)=c1
c
c  compute the next set of parabola coefficients
c  and average with the previous set
      do 1 i=2,npm2
        ip2=i+2
        x1=x2
        x2=x3
        x3=xi(ip2)
        y1=y2
        y2=y3
        y3=yi(ip2)
        h1=h2
        h2=x3-x2
        h3=x3-x1
        d1= y1/(h1*h3)
        d2=-y2/(h1*h2)
        d3= y3/(h2*h3)
        a2=  d1*x2*x3   + d2*x3*x1   + d3*x1*x2
        b2= -d1*(x2+x3) - d2*(x3+x1) - d3*(x1+x2)
        c2=  d1         + d2         + d3
        ai(i)=0.5d0*(a1+a2)
        bi(i)=0.5d0*(b1+b2)
        ci(i)=0.5d0*(c1+c2)
        a1=a2
        b1=b2
        c1=c2
  1   continue
c
c  rightmost interval
      ai(npm1)=a2
      bi(npm1)=b2
      ci(npm1)=c2
c
      return
      end
c
************************************************************************
