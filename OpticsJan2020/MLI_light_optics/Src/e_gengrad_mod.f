!***********************************************************************
!
! e_gengrad: module for E-field generalized gradients
!
! Description: This module implements the derived types and subroutines
! for computing generalized gradients from E-field data given on the
! surface of a cylinder.
!
! Version: 0.1
! Author: D.T.Abell, Tech-X Corp., Jan.2005
!
! Comments
!   12.Jan.05 DTA: Only azimuthally symmetric case implemented.
!
!***********************************************************************
!
      module e_gengrad
        use parallel, only : idproc
        use lieaparam, only : monoms
        use gengrad_data
        implicit none
!
! data
!
      integer, parameter, private :: max_j = 3  ! (N+1)/2 for MaryLie-N
      ! extra z values required by adam11
      integer, parameter, private :: rkxtra = 306
      !double precision :: radius !don't think this is used
      double precision :: rf_phase
      double precision :: rf_escale
!
! derived types
!
      type Efield_data
        integer :: nz_intrvl
        double precision :: zmin, zmax
        double precision :: radius
        double precision :: freq
        ! must allocate the following arrays dim(0:nz_intrvl)
        double precision, dimension(:), pointer :: zvals
        double precision, dimension(:), pointer :: Ezdata
        !double precision, dimension(:), pointer :: Erdata
        !double precision, dimension(:), pointer :: Btdata
      end type Efield_data
!
      type charfn
        integer :: nk_intrvl
        double precision :: w_ovr_c
        double precision :: kmax, dk
        ! must allocate e0r(0:nk_intrvl), e0i(0:nk_intrvl)
        double precision, dimension(:), pointer :: e0r
        double precision, dimension(:), pointer :: e0i
      end type charfn
!
      type vecpot
        double precision, dimension(0:monoms) :: Ax,Ay,Az
      end type vecpot
!
! generalized gradient data to share
!
      type (gengrad), save :: eggrdata
      ! set G1[c,s](iz,j,m) = Crmj^{c,s}(z)
      !     G2[c,s](iz,j,m) = Cfmj^{c,s}(z)
      !     G3[c,s](iz,j,m) = Czmj^{c,s}(z)
!
! vector potential data to share
!
      type (vecpot), save :: vp
!
! functions and subroutines
!
      contains
!
!***********************************************************************
      subroutine read_efield_t7(fn,nf,zi,zf,f,r,efield)
! Read one or more 't7' files and extract the field data at the surface
! of a cylinder of specified radius r.  A t7 file contains rf data for
! Ez, Er, Etot, and Htheta on a uniform grid in r and z.  For a more on
! the t7 format, see, for example, p.35 of the Parmela manual.  If the
! input 'f', for rf cavity frequency, has a non-zero value, it will
! override the value in the field-data files; it is otherwise replaced
! by the value from the field-data files.
      use math_consts
      use parallel
      implicit none
      character(16), intent(in) :: fn  ! file name or base file name
      integer, intent(in) :: nf  ! file number range; 0 => use 'fname'
      double precision, intent(in) :: zi,zf  ! initial and final z
      double precision, intent(inout) :: f  ! cavity frequency
      double precision, intent(in) :: r  ! extract data at this radius
      type (Efield_data), intent(out) :: efield
!-----!----------------------------------------------------------------!
      character(29), parameter :: estrng="e_gengrad_mod::read_efield_t7"
      double precision, parameter :: tiny=1.d-6
      character(16) :: fname
      logical :: singleF
      integer :: ierr,lng0,lng,nfc,nfile,offset
      integer :: i,iu,ii,jj,kk,ll,irad,ir,iz,izt
      integer :: nr,nz,nvals,nztot
      double precision :: frq,r1,r2,ri,z1,z2,dr,dz
      double precision :: ez,er,et,hth
cccccc
c     if (idproc.eq.0) then
c       print *," "
c       print *,"(read_efield_t7): ..."
c       print *,"(",fn,nf,zi,zf,r,")"
c     end if
cccccc
      ! are we reading just a single file
      ! with a given filename?
      if (nf.eq.0) then
        nfc=1
        singleF=.true.
      else
        nfc=nf
        singleF=.false.
        offset=int(log10(real(nfc)))
      end if

      ! loop over the files to read
      nztot=0
      do nfile=1,nfc

        ! first open the E-field data file
        fname=fn
        if (singleF) then  ! use 'fname' as is
          call fnamechk(fname,iu,ierr,estrng)
        else  ! build name of current 't7' file
          lng0=len_trim(fname)  ! DTA: should test length
          if (nfile.ge.1.and.nfile.le.9) then
            lng=lng0+offset
            do i=1,offset
              fname(lng0+i:lng0+i)="0"
            end do
            fname(lng+1:lng+1)=char(nfile+48)
            fname(lng+2:lng+4)=".t7"
            call fnamechk(fname,iu,ierr,estrng)
            if(ierr.ne.0) then
              fname(lng+2:lng+4)=".T7"
              call fnamechk(fname,iu,ierr,estrng)
            end if
          else if (nfile.ge.10.and.nfile.le.99) then
            lng=lng0+offset-1
            do i=1,offset-1
              fname(lng0+i:lng0+i)="0"
            end do
            ii = nfile/10
            jj = nfile - ii*10
            fname(lng+1:lng+1) = char(ii+48)
            fname(lng+2:lng+2) = char(jj+48)
            fname(lng+3:lng+5)=".t7"
            call fnamechk(fname,iu,ierr,estrng)
            if(ierr.ne.0) then
              fname(lng+3:lng+5)=".T7"
              call fnamechk(fname,iu,ierr,estrng)
            end if
          else if (nfile.ge.100.and.nfile.le.999) then
            lng=lng0+offset-2
            do i=1,offset-2
              fname(lng0+i:lng0+i)="0"
            end do
            ii = nfile/100
            jj = nfile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            fname(lng+1:lng+1) = char(ii+48)
            fname(lng+2:lng+2) = char(kk+48)
            fname(lng+3:lng+3) = char(ll+48)
            fname(lng+4:lng+6)=".t7"
            call fnamechk(fname,iu,ierr,estrng)
            if(ierr.ne.0) then
              fname(lng+4:lng+6)=".T7"
              call fnamechk(fname,iu,ierr,estrng)
            end if
          else if (nfile.gt.999) then
            ierr=1
            if (idproc.eq.0) then
              write(6,*) "<*** ERROR ***> in ",estrng,":"
              write(6,*) "   range of file numbers too large!"
            end if
          end if
        end if
        if(ierr.ne.0) then
          if (idproc.eq.0) then
            write(6,*) "<*** ERROR ***> in ",estrng,":"
            write(6,*) "   can't read data file",fname
            write(6,*) "   exiting now..."
          end if
          call myexit
        end if
        if (idproc.eq.0) print *,fname," opened"
        ! file unit iu should now be open for reading
        ! ---but only by processor zero

        ! now read data file
        if (idproc.eq.0) then
          ! read data description: first three lines of t7 file describe
          ! longitudinal range, frequency, and radial range
          ! longitudinal range given in cm; convert to m
          read(iu,*) z1,z2,nz
          print *,"z-in: ",z1,z2,nz
          z1=z1*1.0d-2
          z2=z2*1.0d-2
          ! if argument f is non-zero, use that value;
          ! else read frequency in MHz, and convert to Hz
          read(iu,*) frq
          print *,"frq: ",f
          if (f.eq.0.d0) then
            f=frq*1.0d6
          endif
          ! radial range given in cm; convert to m
          read(iu,*) r1,r2,nr
          print *,"r-in: ",r1,r2,nr
          r1=r1*1.0d-2
          r2=r2*1.0d-2
          print *,"frequency: ",f
          print *,"z-range: ",z1,z2,nz
          print *,"r-range: ",r1,r2,nr
        end if

        call MPI_BCAST(nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(z1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(z2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 
        if (nfile.eq.1) then
          efield%zmin = z1
          efield%zmax = z2
          efield%radius = r
          efield%freq = f
          ! deallocate old data, then allocate space we need
          ! === not sure this is the right approach (DTA) ===
          ! DTA: if (nz(file2) > nz(file1)) problem exists!
          if(associated(efield%zvals)) deallocate(efield%zvals)
          if(associated(efield%Ezdata)) deallocate(efield%Ezdata)
          !if(associated(efield%Erdata)) deallocate(efield%Erdata)
          !if(associated(efield%Btdata)) deallocate(efield%Btdata)
          nvals=nfc*nz  ! best guess for total number of intervals
          allocate(efield%zvals(0:nvals))
          allocate(efield%Ezdata(0:nvals))
          !allocate(efield%Erdata(0:nvals))
          !allocate(efield%Btdata(0:nvals))
        else
          ! DTA: should do some sanity checks here
          efield%zmax=z2
        end if

        ! compute and test radial index
        ri=nr*(r-r1)/(r2-r1)
        irad=int(ri+0.5)
        if (abs(ri-irad).gt.tiny.or.                                    &
     &      irad.lt.0.or.irad.gt.nr) then
          if (idproc.eq.0) then
            print *,"<*** ERROR ***> in ",estrng,":"
            print *,"   desired radius ",r," does not appear in ",fname
            print *,"   r1, r2, nr, dr =",r1,r2,nr,(r2-r1)/nr
            print *,"   ri, irad =",ri,irad
            print *,"   exiting now..."
          end if
          call myexit
        end if

        ! read field data at desired radius
        if (idproc.eq.0) then
          do ir=1,irad  ! skip inner radii
            do iz=0,nz
              read(iu,*) ez,er,et
              read(iu,*) hth
            end do
          end do
          dz=(z2-z1)/nz
          do iz=0,nz
            read(iu,*) ez,er,et
            read(iu,*) hth
            izt=nztot+iz
            efield%zvals(izt)=z1+iz*dz
            ! E-field given in MV/m, convert to V/m
            efield%Ezdata(izt)=ez*1.0d6
            !efield%Erdata(izt)=er*1.0d6
            ! H-field given in A/m, convert to B-field in T
            !efield%Btdata(izt)=hth*mu_o
          end do
          close(iu)
        end if
        nztot=nztot+nz

      end do ! loop over nfile
      efield%nz_intrvl = nztot

      ! check (zmin,zmax) v. (zi,zf)
      if (idproc.eq.0) then
        if (.not.(zi.eq.0.d0.and.zf.eq.0.d0)) then  ! shift z values
          if (abs((efield%zmax-efield%zmin)-(zf-zi)).gt.tiny) then
            if (idproc.eq.0) then
              print *," <*** WARNING ***> from ",estrng
              print *,"   desired range of z values is incompatible"
              print *,"   with the data in file(s) ",fn
            end if
          else
            dz=efield%zmin-zi
            if (dz.ne.0.d0) then
              do iz=0,nztot
                efield%zvals(iz)=efield%zvals(iz)-dz
              end do
            end if
          end if
        end if
      end if
cccccc
c     if (idproc.eq.0) then
c       print *," "
c       print *," Efield data at radius ",r
c       do iz=0,nztot
c         print *,efield%zvals(iz),efield%Ezdata(iz)
c       end do
c     end if
cccccc

      ! broadcast E-field information to other processors
      ! zmin,zmax,nz_intrvl,freq,w_ovr_c already done;
      ! need to broadcast just zvals and Ezdata
      call MPI_BCAST(efield%zvals(0),efield%nz_intrvl+1,                &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(efield%Ezdata(0),efield%nz_intrvl+1,               &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      end subroutine read_efield_t7
!
!***********************************************************************
      subroutine read_efield_ez(fn,nf,jf,zi,zf,f,r,efield)
!
! Read one or more data files and extract the field data at the surface
! of a cylinder of specified radius r.  This subroutine reads field data
! for Ez on a uniform grid in r and z.  The first three lines contain
! meta-data describing the file:
!
!   r_min/m  r_max/m  N_r  (N_r := number of r intervals)
!   z_min/m  z_max/m  N_z  (N_z := number of z intervals)
!   rf-freq/Hz
!
! The remaining lines contain the data E_z(r_i,z_j), where i denotes
! the faster changing index.
!
! If the input 'f', for rf cavity frequency, has a non-zero value, it
! will override the value in the field-data files; it is otherwise
! replaced by the value from the field-data files.
!
      use math_consts
      use parallel
      implicit none
      character(16), intent(in) :: fn  ! file name or base file name
      integer, intent(in) :: nf  ! file number range; 0 => use 'fn'
      integer, intent(in) :: jf  ! unit number for E-field diagnostic
      double precision, intent(inout) :: zi,zf  ! initial and final z
      double precision, intent(inout) :: f  ! cavity frequency
      double precision, intent(in) :: r  ! extract data at this radius
      type (Efield_data), intent(out) :: efield
!-----!----------------------------------------------------------------!
      character(29), parameter :: estrng="e_gengrad_mod::read_efield_ez"
      double precision, parameter :: tiny=1.d-6
      double precision, parameter :: eps=1.d-13
      character(16) :: fname,numstr
      logical :: singleF
      integer :: ierr,lng,ndig,nfc,nfile
      integer :: i,iu,ii,jj,kk,ll,irad,ir,iz,izt
      integer :: nr,nz,nvals,nztot
      double precision :: dr,dz,frq,r1,r2,ri,z1,z2
      double precision :: ez,er,et,hth
cccccc
c     if (idproc.eq.0) then
c       print *," "
c       print *,estrng,": ..."
c       print *,"  (",fn,",",nf,",",jf,",",zi,",",zf,",",r,",efield)"
c     end if
cccccc
      ! are we reading just a single file with a given filename?
      if (nf.eq.0) then
        nfc=1
        singleF=.true.
      else
        nfc=nf
        singleF=.false.
      end if

      ! loop over the files to read
      ndig=1+int(log10(real(nfc)))
      nztot=0
      do nfile=1,nfc

        ! first open the E-field data file
        fname=fn
        if (singleF) then  ! use 'fname' as is
          call fnamechk(fname,iu,ierr,estrng)
        else  ! build name of current data file
          lng=len_trim(fname)  ! DTA: should test length
          call num2string(nfile,numstr,ndig)
          fname=fname(1:lng)//numstr(1:ndig)//".ez"
          call fnamechk(fname,iu,ierr,estrng)
        end if
        if(ierr.ne.0) then
          if (idproc.eq.0) then
            write(6,*) "<*** ERROR ***> in ",estrng,":"
            write(6,*) "   can't read data file ",fname
            write(6,*) "   exiting now..."
          end if
          call myexit
        end if
        if (idproc.eq.0) print *,fname," opened"
        ! file unit iu should now be open for reading
        ! ---but only by processor zero

        ! now read data file
        if (idproc.eq.0) then
          ! read data description: first three lines of ez file describe
          !    r_min/m  r_max/m  N_r_intrvls
          !    z_min/m  z_max/m  N_z_intrvls
          !    rf-frequency/Hz
          read(iu,*) r1,r2,nr
          read(iu,*) z1,z2,nz
          read(iu,*) frq
          ! if argument f is non-zero, use that value instead
          if (f.eq.0.d0) then
            f=frq
          else
            print *," NB: not using rf-frequency ",frq," from data file"
          endif
          print *,"r-range: ",r1,r2,nr
          print *,"z-range: ",z1,z2,nz
          print *,"frequency: ",f
        end if

        call MPI_BCAST(nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(z1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(z2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(f,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 
        if (nfile.eq.1) then
          efield%zmin = z1
          efield%zmax = z2
          efield%radius = r ! fix (below) if requested and actual differ
          efield%freq = f
          ! deallocate old data, then allocate space we need
          ! === not sure this is the right approach (DTA) ===
          ! DTA: if (nz(file_i) > nz(file_1)) problem exists!
          if(associated(efield%zvals)) deallocate(efield%zvals)
          if(associated(efield%Ezdata)) deallocate(efield%Ezdata)
          !if(associated(efield%Erdata)) deallocate(efield%Erdata)
          !if(associated(efield%Btdata)) deallocate(efield%Btdata)
          nvals=nfc*nz  ! best guess for total number of intervals
          allocate(efield%zvals(0:nvals))
          allocate(efield%Ezdata(0:nvals))
          !allocate(efield%Erdata(0:nvals))
          !allocate(efield%Btdata(0:nvals))
        else
          ! DTA: should do some sanity checks here
          efield%zmax=z2
        end if

        ! compute and test radial index
        if (nr.ne.0) then 
          dr=(r2-r1)/nr
          ri=(r-r1)/dr
        else ! we have data at a single radius
          dr=0
          ri=0
        end if
        irad=nint(ri)
        efield%radius=abs(r1+irad*dr) ! here store actual radius
        if (abs(ri-irad).gt.tiny.or.                                    &
     &      irad.lt.0.or.irad.gt.nr) then
          if (idproc.eq.0) then
            print *,"<*** ERROR ***> in ",estrng,":"
            print *,"   desired radius ",r," does not appear in ",fname
            print *,"   r1, r2, nr, dr = ",r1,r2,nr,dr
            print *,"   ri, irad = ",ri,irad
            print *,"   exiting now..."
          end if
          call myexit
        else
          print *,"extracting E-field data a radius ",efield%radius
        end if

        ! read field data at desired radius
        if (idproc.eq.0) then
          dz=(z2-z1)/nz
          do iz=0,nz
            ! skip inner radii
            do ir=0,irad-1
              read(iu,*) ez
            end do
            ! read and store next value (but don't duplicate boundary)
            read(iu,*) ez
            izt=nztot+iz ! overwrites last value at file boundary (iz=0)
            efield%zvals(izt)=z1+iz*dz
            efield%Ezdata(izt)=ez
            ! skip outer radii
            do ir=irad+1,nr
              read(iu,*) ez
            end do
          end do
          close(iu)
        end if
        nztot=nztot+nz

      end do ! loop over nfile
      efield%nz_intrvl = nztot

      ! check (zmin,zmax) v. (zi,zf)
      if (zi.eq.0.d0.and.zf.eq.0.d0) then
        ! assign zi and zf from file
        zi=efield%zmin
        zf=efield%zmax
      else if (abs((zf-zi)-(efield%zmax-efield%zmin)).lt.eps) then
        ! shift data
        dz=zi-efield%zmin
        if (dz.ne.0.d0) then
          if (idproc.eq.0) then
            do iz=0,nztot
              efield%zvals(iz)=efield%zvals(iz)+dz
            end do
          end if
          efield%zmin=efield%zmin+dz
          efield%zmax=efield%zmax+dz
        end if
        ! at this point, zi.eq.efield%zmin
        ! and zf lies within eps of efield%zmax
        if (zf.gt.efield%zmax) zf=efield%zmax
      else if ((zf-zi).gt.(efield%zmax-efield%zmin)) then
        ! complain, and assign zi and zf from file
        if (idproc.eq.0) then
          print *," <*** WARNING ***> from ",estrng
          print *,"   desired z-range [",zi,",",zf,"]"
          print *,"   is incompatible with the data in file(s) ",fn
          print *,"   using z-range [",efield%zmin,",",efield%zmax,"]"
        end if
        zi=efield%zmin
        zf=efield%zmax
      else
        ! deal with cases (zf-zi).lt.(efield%zmax-efield%zmin)
        if (zi.lt.efield%zmin) then
          if (idproc.eq.0) then
            print *," <*** WARNING ***> from ",estrng
            print *,"   changing desired z-range from [",zi,",",zf,"]"
            print *,"   to [",efield%zmin,",",zf,"]"
          end if
          zi=efield%zmin
        else if (zf.lt.efield%zmax) then
          if (idproc.eq.0) then
            print *," <*** WARNING ***> from ",estrng
            print *,"   changing desired z-range from [",zi,",",zf,"]"
            print *,"   to [",zi,",",efield%zmin,"]"
          end if
          zf=efield%zmax
        end if
      end if

      ! write E-field(s) at radius r
      if (idproc.eq.0) then
        if (jf.ne.0) then
          do iz=0,nztot
            write(jf,*) efield%zvals(iz),efield%Ezdata(iz)
          end do
        end if
      end if

      ! broadcast E-field information to other processors
      ! zmin,zmax,nz_intrvl,freq,w_ovr_c already done;
      ! need to broadcast just zvals and Ezdata
      call MPI_BCAST(efield%zvals(0),efield%nz_intrvl+1,                &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(efield%Ezdata(0),efield%nz_intrvl+1,               &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      end subroutine read_efield_ez
!
!***********************************************************************
      subroutine comp_charfn(efield,kmx,nk_intrvl,kfile,e0)
! Compute the characteristic function \tilde{e}_0(k) from a given set
! of electric field data at radius r.
      use math_consts
      use phys_consts
      use parallel, only : idproc
      implicit none
      type (Efield_data), intent(in) :: efield
      double precision, intent(in) :: kmx ! maximum value of k
      integer, intent(in) :: nk_intrvl ! number of intervals in k
      integer, intent(in) :: kfile ! file unit number for output of e0
      type (charfn), intent(out) :: e0
!-----!----------------------------------------------------------------!
      double precision, parameter :: tiny=1.d-6
      integer :: ik,iz,izb,izt
      double precision :: bessR,ei,er,invrt2pi
      double precision :: k,kapsq,klsq,kR,dz,dz1
      double precision, dimension(1) :: bv
cccccc
c     if (idproc.eq.0) then
c       print *,"(comp_charfn): ..."
c     end if
cccccc
      e0%nk_intrvl=nk_intrvl
      e0%w_ovr_c=twopi*efield%freq/c_light
      e0%kmax=kmx
      e0%dk=kmx/nk_intrvl
      klsq=e0%w_ovr_c**2
      invrt2pi=1.d0/sqrt(twopi)

      ! allocate memory for e0%e0r and e0%e0i
      ! deallocate old data, then allocate space we need
      ! === not sure this is the best approach (DTA) ===
      if(associated(e0%e0r)) deallocate(e0%e0r)
      if(associated(e0%e0i)) deallocate(e0%e0i)
      allocate(e0%e0r(0:nk_intrvl))
      allocate(e0%e0i(0:nk_intrvl))
      e0%e0r=0.d0
      e0%e0i=0.d0

      ! dz may not be uniform across the E-field data files, but the
      ! filon integrator requires dz = constant: perform z integration
      ! in discrete sections
      izt=0  ! top index of current section

      do while (izt < efield%nz_intrvl)
        ! set izb, and initialize search for next izt
        izb=izt; izt=izb+1
        dz=efield%zvals(izt)-efield%zvals(izb)
        do while (izt.lt.efield%nz_intrvl.and.                          &
     &           dabs(efield%zvals(izt+1)-efield%zvals(izt)-dz).lt.tiny)
          izt=izt+1
        end do
        ! compute real and imaginary parts of e0(k), k in [0,kmax]
cccccc
c       if (idproc.eq.0) then
c         print *,"iz-range:",izb,izt,"(",izt-izb+1,")"
c       end if
cccccc
        do ik=0,nk_intrvl
          k=ik*e0%dk
          call besselR0(k,e0%w_ovr_c,efield%radius,bessR)
          call filon_io(efield%zvals,efield%Ezdata,k,                   &
     &                  izb,izt-izb+1,ei,er)
          e0%e0r(ik)=e0%e0r(ik)+invrt2pi*er/bessR
          e0%e0i(ik)=e0%e0i(ik)-invrt2pi*ei/bessR
        end do  ! loop over k
      end do  ! loop over sections w/ dz=constant

      ! if asked, write characteristic function to file
      if (idproc.eq.0.and.kfile.ne.0) then
        write(kfile,*) "# ",e0%nk_intrvl
        write(kfile,*) "# ",e0%w_ovr_c
        write(kfile,*) "# ",e0%kmax,e0%dk
        do ik=0,nk_intrvl
          k=ik*e0%dk
          write(kfile,*) k,e0%e0r(ik),e0%e0i(ik)
        end do
      end if
!
      end subroutine comp_charfn
!
!***********************************************************************
      subroutine read_charfn(iu,e0)
! Read characteristic function from file.
      use parallel
      implicit none
      integer, intent(in) :: iu ! file unit number for data file
      type (charfn), intent(out) :: e0 ! characteristic function
!-----!----------------------------------------------------------------!
      character(26), parameter :: estrng="e_gengrad_mod::read_charfn"
      character(90) :: string
      integer :: ierr,ik,ist,len,myerr
      real*8 :: k

!      if (idproc.eq.0) then
!        print *," reading characteristic functions ..."
!      end if

      ! check file unit number
      ! NB: only PE0 knows the argument iu
      myerr=0
      if (idproc.eq.0) then
        if (iu.eq.0) then
          myerr=1
          write(6,*) "<*** ERROR ***> in ",estrng,": iu=0"
          write(6,*) '<**ERROR**> e_gengrad_mod::read_charfn: iu=0'
        endif
      endif
      call ibcast(myerr)
      if (myerr.ne.0) call myexit()

      ! read characteristic function from file
      !   first read and broadcast parameters from top three lines
      !   (skipping possible comment character '#' at start of line)
      if (idproc.eq.0) then
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) e0%nk_intrvl
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) e0%w_ovr_c
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) e0%kmax,e0%dk
      end if
      call MPI_BCAST(e0%nk_intrvl,1,                                    &
     &               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(e0%w_ovr_c,1,                                      &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(e0%kmax,1,                                         &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(e0%dk,1,                                           &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      !   allocate arrays
      if(associated(e0%e0r)) deallocate(e0%e0r)
      if(associated(e0%e0i)) deallocate(e0%e0i)
      allocate(e0%e0r(0:e0%nk_intrvl))
      allocate(e0%e0i(0:e0%nk_intrvl))
      !   then read data and close file
      if (idproc.eq.0) then
        do ik=0,e0%nk_intrvl
          read(iu,*) k,e0%e0r(ik),e0%e0i(ik)
        end do
        close(iu)
      end if
      ! broadcast char. fn. information to other processors
      call MPI_BCAST(e0%e0r(0:e0%nk_intrvl),(e0%nk_intrvl+1),           &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(e0%e0i(0:e0%nk_intrvl),(e0%nk_intrvl+1),           &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
      end subroutine read_charfn
!
!***********************************************************************
      subroutine comp_egengrads(iu,e0,zi,zf,nz_intrvl)
! Compute electric field generalized gradients and write to unit iu.
! [NB: currently handles only azimuthally symmetric case.]
      use math_consts
      use phys_consts
      use parallel, only : idproc
      implicit none
      integer, intent(in) :: iu ! file unit number for output
      type (charfn), intent(in) :: e0
      double precision, intent(in) :: zi,zf ! initial and final z
      integer, intent(in) :: nz_intrvl ! number of z intervals
!-----!----------------------------------------------------------------!
      double precision, parameter :: tiny=1.d-6
      integer :: ik,iz,j,nkpts
      double precision :: rt2ovrpi,z
      double precision :: cosint,sinint
      double precision, allocatable, dimension(:) :: kvals,sk,skj,skjr, &
     &                                               integrand
cccccc
c     if (idproc.eq.0) then
c       print *,"(comp_egengrads): ..."
c     end if
      print *,idproc,": (comp_egengrads): ..."
cccccc
      ! initialize constants in eggrdata
      eggrdata%maxj=max_j
      eggrdata%maxm=0
      eggrdata%zmin=zi
      eggrdata%zmax=zf
      eggrdata%nz_intrvl=nz_intrvl
      eggrdata%dz=(zf-zi)/nz_intrvl
      eggrdata%angfrq=e0%w_ovr_c*c_light
cccccc
c     if (idproc.eq.0) then
c       print *,"   allocating memory for zvals and gen.grad. arrays..."
c     end if
      print *," proc.",idproc,                                          &
     &        ": allocating memory for zvals and gen.grad. arrays..."
cccccc
      ! allocate memory for eggrdata%{zvals,G1c,G3c}
      ! deallocate old data, then allocate space we need
      ! === not sure this is the best approach (DTA) ===
      if(associated(eggrdata%zvals)) deallocate(eggrdata%zvals)
      if(associated(eggrdata%G1c)) deallocate(eggrdata%G1c)
      !if(associated(eggrdata%G2c)) deallocate(eggrdata%G2c)
      if(associated(eggrdata%G3c)) deallocate(eggrdata%G3c)
      allocate(eggrdata%zvals(0:nz_intrvl+rkxtra))
      allocate(eggrdata%G1c(0:nz_intrvl+rkxtra,1:eggrdata%maxj,         &
     &                                         0:eggrdata%maxm))
      !allocate(eggrdata%G2c(0:nz_intrvl+rkxtra,1:eggrdata%maxj,0:0))
      allocate(eggrdata%G3c(0:nz_intrvl+rkxtra,0:eggrdata%maxj,         &
     &                                         0:eggrdata%maxm))
cccccc
c     if (idproc.eq.0) then
c       print *,": done allocating zvals and gen.grad. arrays..."
c     end if
      print *," proc.",idproc,                                          &
     &        ": done allocating zvals and gen.grad. arrays..."
cccccc

      ! define scalar variables we need
      nkpts=e0%nk_intrvl+1
      rt2ovrpi=sqrt(2.d0/pi)

      ! allocate vectors we need
      allocate(kvals(0:e0%nk_intrvl))
      allocate(sk(0:e0%nk_intrvl))
      allocate(skj(0:e0%nk_intrvl))
      allocate(skjr(0:e0%nk_intrvl))
      allocate(integrand(0:e0%nk_intrvl))

      ! initialize vectors
      do ik=0,e0%nk_intrvl
        kvals(ik)=ik*e0%dk
      end do
      sk=kvals**2-e0%w_ovr_c**2
      skj=1.d0

      ! compute generalized gradients
      call zs_adam11(zi,zf,nz_intrvl,eggrdata%zvals)
      if (idproc.eq.0) then
        print *," Computing E-field generalized gradients ..."
      end if
      do j=0,eggrdata%maxj
        if (j.ne.0) then
          ! compute Crc(0,j,z) = eggrdata%G1c(z,j,0)
          if (idproc.eq.0) then
            print *,"      ... Crc(0,",j,",z)"
          end if
          skjr=j*kvals*skj
          integrand=skjr*e0%e0r
          do iz=0,nz_intrvl+rkxtra
            z=eggrdata%zvals(iz)
            call filon_io(kvals,integrand,z,0,nkpts,sinint,cosint)
            eggrdata%G1c(iz,j,0)=sinint
          end do
          integrand=skjr*e0%e0i
          do iz=0,nz_intrvl+rkxtra
            z=eggrdata%zvals(iz)
            call filon_io(kvals,integrand,z,0,nkpts,sinint,cosint)
            eggrdata%G1c(iz,j,0)=rt2ovrpi*(eggrdata%G1c(iz,j,0)+cosint)
          end do
        end if
        ! compute Czc(0,j,z) = eggrdata%G3c(z,j,0)
        if (idproc.eq.0) then
          print *,"      ... Czc(0,",j,",z)"
        end if
        if (j.ne.0) then
          skj=skj*sk
        end if
        integrand=skj*e0%e0r
        do iz=0,nz_intrvl+rkxtra
          z=eggrdata%zvals(iz)
          call filon_io(kvals,integrand,z,0,nkpts,sinint,cosint)
          eggrdata%G3c(iz,j,0)=cosint
        end do
        integrand=skj*e0%e0i
        do iz=0,nz_intrvl+rkxtra
          z=eggrdata%zvals(iz)
          call filon_io(kvals,integrand,z,0,nkpts,sinint,cosint)
          eggrdata%G3c(iz,j,0)=rt2ovrpi*(eggrdata%G3c(iz,j,0)-sinint)
        end do
      end do  ! loop over j

      ! write generalized gradients to file
      if(idproc.eq.0.and.iu.ne.0) then
        write(iu,*) "# ",eggrdata%nz_intrvl,eggrdata%maxj,eggrdata%maxm
        write(iu,*) "# ",eggrdata%zmin,eggrdata%zmax,eggrdata%dz
        write(iu,*) "# ",eggrdata%angfrq
        do iz=0,nz_intrvl+rkxtra
          z=eggrdata%zvals(iz)
          write(iu,*) z,eggrdata%G3c(iz,0,0),                           &
     &                 (eggrdata%G1c(iz,j,0),eggrdata%G3c(iz,j,0),      &
     &                  j=1,eggrdata%maxj)
        end do
      end if

      ! write z, Ez, dEz/dz
      !if(idproc.eq.0.and.iude.ne.0) then
      if(idproc.eq.0) then
        do iz=0,nz_intrvl+rkxtra
          write(81,*) z,eggrdata%G3c(iz,0,0),-eggrdata%G1c(iz,1,0)
    !     write(81,'(3(1x,1pe22.15))') z,eggrdata%G3c(iz,0,0),          &
    !&                                  -eggrdata%G1c(iz,1,0)
        end do
      end if

      ! clean up
      deallocate(kvals)
      deallocate(sk)
      deallocate(skj)
      deallocate(skjr)
      deallocate(integrand)
!
      end subroutine comp_egengrads
!
!***********************************************************************
      subroutine read_egengrads(iu)
! Read generalized gradients from file.
      use parallel
      implicit none
      integer, intent(in) :: iu ! file unit number for data file
!-----!----------------------------------------------------------------!
      character(29), parameter :: estrng="e_gengrad_mod::read_egengrads"
      character(90) :: string
      integer :: ierr,ist,iz,j,len,myerr
      double precision :: z

      if (idproc.eq.0) then
        print *," reading generalized gradients ..."
      end if

      ! check file unit number
      ! NB: only PE0 knows the argument iu
      myerr=0
      if (idproc.eq.0) then
        if (iu.eq.0) then
          myerr=1
          write(6,*) "<*** ERROR ***> in ",estrng,": iu=0"
        endif
      endif
      call ibcast(myerr)
      if (myerr.ne.0) call myexit()

      ! read generalized gradients from file
      !   first read and broadcast parameters from top three lines
      !   (skipping possible comment character '#' at start of line)
      if(idproc.eq.0.and.iu.ne.0) then
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) eggrdata%nz_intrvl,eggrdata%maxj,     &
     &                            eggrdata%maxm
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) eggrdata%zmin,eggrdata%zmax,          &
     &                            eggrdata%dz
        read(iu,'(a)') string
        len=len_trim(string)
        ist=index(trim(string),'#')
        read(string(ist+1:len),*) eggrdata%angfrq
      end if
      call MPI_BCAST(eggrdata%nz_intrvl,1,                              &
     &               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%maxj,1,                                   &
     &               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%maxm,1,                                   &
     &               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%zmin,1,                                   &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%zmax,1,                                   &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%dz,1,                                     &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%angfrq,1,                                 &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      !   allocate arrays
      if(associated(eggrdata%zvals)) deallocate(eggrdata%zvals)
      if(associated(eggrdata%G3c)) deallocate(eggrdata%G3c)
      if(associated(eggrdata%G1c)) deallocate(eggrdata%G1c)
      allocate(eggrdata%zvals(0:eggrdata%nz_intrvl+rkxtra))
      allocate(eggrdata%G1c(0:eggrdata%nz_intrvl+rkxtra,                &
     &                      1:eggrdata%maxj,0:eggrdata%maxm))
      allocate(eggrdata%G3c(0:eggrdata%nz_intrvl+rkxtra,                &
     &                      0:eggrdata%maxj,0:eggrdata%maxm))
      !   then read data and close file
      if(idproc.eq.0.and.iu.ne.0) then
        do iz=0,eggrdata%nz_intrvl+rkxtra
          read(iu,*) eggrdata%zvals(iz),eggrdata%G3c(iz,0,0),           &
     &               (eggrdata%G1c(iz,j,0),eggrdata%G3c(iz,j,0),        &
     &                j=1,eggrdata%maxj)
        end do
        close(iu)
      end if
      ! broadcast gen. grad. information to other processors
      call MPI_BCAST(eggrdata%zvals(0:eggrdata%nz_intrvl+rkxtra),       &
     &               (eggrdata%nz_intrvl+rkxtra+1),                     &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%G1c(0:eggrdata%nz_intrvl+rkxtra,          &
     &                            1:eggrdata%maxj,0:0),                 &
     &               (eggrdata%nz_intrvl+rkxtra+1)*(eggrdata%maxj),     &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eggrdata%G3c(0:eggrdata%nz_intrvl+rkxtra,          &
     &                            0:eggrdata%maxj,0:0),                 &
     &               (eggrdata%nz_intrvl+rkxtra+1)*(eggrdata%maxj+1),   &
     &               MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
      end subroutine read_egengrads
!
!***********************************************************************
      subroutine besselR0(k,kl,r,b)
! Compute, the values for the order zero "Bessel" function
! R_0(k,kl,r) := { J_0(sqrt(|k^2-kl^2|)*r), sgn(k^2-kl^2) < 0;
!                { I_0(sqrt(|k^2-kl^2|)*r), otherwise.
! In words, R_0 switches from the regular to the modified Bessel
! function of the first kind as |k| crosses |kl|.
      implicit none
      double precision, intent(in) :: k,kl,r
      double precision, intent(out) :: b
!-----!----------------------------------------------------------------!
      double precision :: kapsq,kR
      double precision, dimension(1) :: bv
!
      kapsq=k**2-kl**2
      if (kapsq.lt.0.d0) then  ! regular Bessel function
        kR=sqrt(-kapsq)*r
        call dbessj0(kR,b)
      else  ! modified Bessel function
        kR=sqrt(kapsq)*r
        call BESSIn(1,0,kR,bv,1)
        b=bv(1)
      endif
!
      end subroutine besselR0
!
!***********************************************************************
      subroutine besselR(m,k,kl,r,bv)
! Compute, for integer orders 0...m, values for the "Bessel"
! function R_m(k,kl,r) := { J_m(sqrt(|k^2-kl^2|)*r), sgn(k^2-kl^2) < 0;
!                         { I_m(sqrt(|k^2-kl^2|)*r), otherwise.
! In words, R_m switches from the regular to the modified Bessel
! function of the first kind as |k| crosses |kl|.
      implicit none
      integer, intent(in) :: m  ! desired orders 0...m-1
      double precision, intent(in) :: k,kl,r 
      double precision, dimension(:), intent(out) :: bv ! bv(0:m)
!-----!----------------------------------------------------------------!
      double precision :: kapsq,kR
      double precision, dimension(:), allocatable :: bessR
!
      allocate(bessR(m+1))
!
      kapsq=k**2-kl**2
      if (kapsq.lt.0.d0) then  ! regular Bessel function
        kR=sqrt(-kapsq)*r
        call dbessjm(m+1,kR,bessR)
      else  ! modified Bessel function
        kR=sqrt(kapsq)*r
        call BESSIn(m+1,0,kR,bv,1)  ! unscaled I_0 ... I_m
        bessR=bv(1)
      endif
      bv=bessR
!
      deallocate(bessR)
!
      end subroutine besselR
!
!***********************************************************************
      subroutine filon_io(xn,fn,k,io,npts,sinint,cosint)
!
!  Generic Filon integrator---array version with variable index origin:
!  This subroutine uses Filon's method to evaluate the integrals from
!  xn(io) to xn(io+npts-1) of
!
!            fn(x) sin k*x dx
!       and  fn(x) cos k*x dx
!
!  where fn is also an array of length npts with the same index origin
!  io.  Filon's method allows one to evaluate the integral of an
!  oscillating function without having to follow every wiggle with many
!  evaluation points.  This version has only a single frequency, k.
!  (See F.B.Hildebrand, Introduction to Numerical Analysis, 2nd ed.,
!  McGraw-Hill, 1974 (Dover reprint, 1987), Sect. 3.10.  See, also,
!  Abramowitz & Stegun, p.890.)  For k=0, the cosine integral reduces to
!  the extended form of Simpson's rule.
!
!  NB: The input array xn must have equally spaced points.  Moreover, it
!  should contain an _odd_ number of points---or an even number of
!  intervals.  (If xn contains an even number of points, the rightmost
!  point is ignored.)  Also note that the usual definitions of Filon's
!  formulas use index origin zero, so that the evaluation points are
!  {x_0, x_1, ..., x_2n}.  But many Fortran arrays have index origin
!  one, which means the "even" points are those indexed 1, 3, 5, ...,
!  while the "odd" points are those indexed 2, 4, 6, ....  Sigh, ...!
!  The implementation given here allows for an arbitrary index origin.
!
!  Written by Peter Walstrom.
!  Jan.2005 (DTA): Rewritten so that one may integrate over a portion
!                  of the data range.  Also changed to implicit none.
!-----!----------------------------------------------------------------!
      use parallel, only : idproc
      implicit none
!
!  arguments
      double precision, dimension(0:), intent(in) :: xn,fn
      double precision, intent(in) :: k
      integer, intent(in) :: npts,io
      double precision, intent(out) :: sinint,cosint
!
!  local variables
      double precision, parameter :: half=0.5d0,one=1.d0,zero=0.d0
      integer :: ii,ix,istart,iend,ievod,npoints
      double precision :: h,theta,alf,bet,gam
      double precision :: kx,coskx,sinkx,f,wt
      double precision, dimension(3) :: dsinint,dcosint
!
!  check parity of npts; if even, emit warning and reduce it by 1
      npoints=npts
      if (2*(npts/2).eq.npts) then
        if (idproc.eq.0) then
          write(6,*) '<*** WARNING ***> from subroutine filon_io():'
          write(6,*) '  must have an odd number of data points;'
          write(6,*) '  ignoring last point!'
        end if
        npoints=npts-1
      end if
!
!  set array endpoints
      istart=io
      iend=io+npoints-1
!
!  note step-size and get Filon weights
      h=xn(istart+1)-xn(istart)
      theta=h*k
      call filon_wts(theta,alf,bet,gam)
c      if (idproc.eq.0) then
c        write(*,200) "filon weights:",theta,alf,bet,gam
c      end if
c  200 format(1x,4(1pd11.4,1x))
!
!  initialize intermediate arrays to zero
      do ievod=1,3
        dsinint(ievod)=zero
        dcosint(ievod)=zero
      end do
!
!  perform Filon integration:
!
!  ievod=1 -->  odd points
!  ievod=2 --> even points
!  ievod=3 -->  end points
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
!  end points (upper end first)
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
!  now sum the Filon-weighted contributions from ievod=1,2,3
      sinint=h*(alf*dsinint(3)+bet*dsinint(2)+gam*dsinint(1))
      cosint=h*(alf*dcosint(3)+bet*dcosint(2)+gam*dcosint(1))
!
      return
      end subroutine filon_io
!
!***********************************************************************
      subroutine zs_adam11(zmin,zmax,nzsteps,zvals)
! Compute for given (zmin,zmax,nzsteps) the z-values required by adam11
! Note: This routine requires size(zvals) >= nzsteps+rkxtra.
      implicit none
      double precision, intent(in) :: zmin,zmax ! end-points in z
      integer, intent(in) :: nzsteps ! number of z integration steps
      double precision, dimension(0:), intent(out) :: zvals
!-----!----------------------------------------------------------------!
      double precision, dimension(7) :: rkvec
      integer :: ih,ir,iz,nzvals
      double precision :: dz,hq,z
c
      ! define the step-size fractions used in the initial Runge-Kutta
      ! steps required by the 11th-order Adams integration routine
      rkvec(1)=0.d0
      rkvec(2)=1.d0/9.d0
      rkvec(3)=1.d0/6.d0
      rkvec(4)=1.d0/3.d0
      rkvec(5)=1.d0/2.d0
      rkvec(6)=2.d0/3.d0
      rkvec(7)=5.d0/6.d0
c
      ! compute z values
      dz=(zmax-zmin)/nzsteps
      hq=dz/5.d0
      nzvals=-1
      do iz=0,nzsteps
        if (iz.lt.9) then
          do ih=0,4
            do ir=1,7
              nzvals=nzvals+1
              zvals(nzvals)=zmin+iz*dz+(ih+rkvec(ir))*hq
            end do
          end do
        else
          nzvals=nzvals+1
          zvals(nzvals)=zmin+iz*dz
        endif
      end do
cccccc
c     write(6,*) "size(zvals)=",size(zvals)
c     write(6,*) "zmin, zmax, nzsteps =",zmin,zmax,nzsteps
c     write(6,*) "dz, hq =",dz,hq
c     do iz=0,nzvals
c       write(6,*) iz,zvals(iz)
c     end do
cccccc
      return
      end subroutine zs_adam11

      end module e_gengrad

