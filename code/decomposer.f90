! -----------------------------------------------------------------------
! partition array for all possibilities
! nx = global # in fast dimension
! ny = global # in mid dimension
! nz = global # in slow dimension

! idecomp = 0 = xyz (blocks)
!           1 = xy, 2 = yz, 3 = xz
!           4 = x, 5 = y, 6 = z

      subroutine decompose(node,mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,   &
     &                     idecomp,npxarg,npyarg,npzarg,ilo,ihi,jlo,jhi,klo,khi)
! modified Plimpton routine to allow arbitrary lower bounds (nxlo,nylo,nzlo) and upper bounds (nxhi,nyhi,nzhi)
! also, this version allows the user to specify npx,npy,npz on input (in which case idecomp must equal -1)
      implicit none
      integer :: node,mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,idecomp,npxarg,npyarg,npzarg,ilo,ihi,jlo,jhi,klo,khi
      integer :: nx,ny,nz,npx,npy,npz,node_x,node_y,node_z
!consistency checks:
      if(.not.(idecomp.ge.-1.or.idecomp.le.6))then
        write(6,*)'(decompose) idecomp error: idecomp=',idecomp
        stop
      endif
      if((idecomp.eq.-1).and.(npxarg*npyarg*npzarg.ne.mprocs))then !this could be commented out if it is checked in the main program
        write(6,*)'error: idecomp.eq.-1 but npx*npy*npz .ne. mprocs'
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
        stop
      endif
      ilo=ilo+nxlo-1
      ihi=ihi+nxlo-1
      jlo=jlo+nylo-1
      jhi=jhi+nylo-1
      klo=klo+nzlo-1
      khi=khi+nzlo-1
      return
      end
! -----------------------------------------------------------------------
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
      return
      end
!
! ----------------------------------------------------------------------
! assign mprocs to 3-d nx,ny,nz grid so as to minimize surface area
! return 3-d px,py,pz grid of procs
      subroutine proc2grid3d(mprocs,nx,ny,nz,px,py,pz)
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
     &             boxx*boxy*boxz.gt.bestboxx*bestboxy*bestboxz)) then
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
        stop
      endif
      return
      end
! ----------------------------------------------------------------------
! assign mprocs to 2d nx,ny grid so as to minimize surface area
! return 2d px,py grid of procs

      subroutine proc2grid2d(mprocs,nx,ny,px,py)
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
          if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and.                         &
     &         boxx*boxy.gt.bestboxx*bestboxy)) then
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
        stop
      endif

      return
      end

! -----------------------------------------------------------------------
      function iowner(iloa,ihia,jloa,jhia,kloa,khia,mprocs,             &
     &                nxlo,nxhi,nylo,nyhi,nzlo,nzhi,i,j,k)
! determine the identity (iowner) of the PE that
! "owns" the grid point with global indices (i,j,k)
      implicit none
      integer :: iowner
      integer :: mprocs,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,i,j,k
      integer :: im,jm,km,n
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer, dimension(0:mprocs-1) :: iloa,ihia,jloa,jhia,kloa,khia
      if(i.lt.nxlo.or.i.gt.nxhi)write(6,'("iowner_input_error:",3i7)') &
     &   i,nxlo,nxhi
      if(j.lt.nylo.or.j.gt.nyhi)write(6,'("iowner_jnput_error:",3i7)') &
     &   j,nylo,nyhi
      if(k.lt.nzlo.or.k.gt.nzhi)write(6,'("iowner_knput_error:",3i7)') &
     &   k,nzlo,nzhi
      im=i
      jm=j
      km=k
!periodic boundary conditions:
!     if(i.eq.0)im=nx
!     if(i.eq.nx+1)im=1
!     if(j.eq.0)jm=ny
!     if(j.eq.ny+1)jm=1
!     if(k.eq.0)km=nz
!     if(k.eq.nz+1)km=1
      do n=0,mprocs-1
!         write(6,'(i6,3x,2i5,7x,2i5,7x,2i5)')n,iloa(n),ihia(n),jloa(n),&
!    &                                        jhia(n),kloa(n),khia(n)
!      call decompose(n,mprocs,nx,ny,nz,idecomp,ilo,ihi,jlo,jhi,klo,khi)
       if(im.lt.iloa(n) .or. im.gt.ihia(n))cycle
       if(jm.lt.jloa(n) .or. jm.gt.jhia(n))cycle
       if(km.lt.kloa(n) .or. km.gt.khia(n))cycle
       iowner=n
       return
      enddo
      write(6,'("iowner_error",3i7)')i,j,k
      iowner=-1
      if(iowner.eq.-1)stop
      return
      end function iowner
!
!============================================================
! I'm including subroutine pmove in this file. It doesn't really belong here
! but I'm putting it here for convenience.  R.D. Ryne
!
      subroutine pmove(y,iownerptcl,n1,n2,nraysp,nbufsize,ifail)
! Given a 2D array y(n1,n2), and a 1D array iownerptcl(n2) that describes who owns each row,
! move the rows to the proc that owns them.
! This routine could be greatly improved, but it appears to work.  R.D. Ryne
      USE mpi
      implicit none
      integer :: n1,n2,nraysp,nbufsize,ifail,n,nsbuf,nrecv,nshift,ntarget,nsender,nraysgbl
      real(8), dimension(n1,n2) :: y
      integer, dimension(1:n2) :: iownerptcl
      real(8), dimension(n1,nbufsize) :: sbuf,rbuf
      integer, dimension(2) :: msid
      real(8), dimension(n1,n2) :: ynew
      integer :: nrayspnew
      integer, dimension(MPI_STATUS_SIZE) :: mpistat
      integer :: mprocs,myrank,ierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!
      call MPI_REDUCE(nraysp,nraysgbl,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ifail=0
      nrayspnew=0
!#    do nshift=1,mprocs-1 !the nshift'th nearest neighbor in order from 1 to mprocs-1
      do nshift=0,mprocs-1 !if storing in a "new" array, then have to start with nshift=0
        ntarget=mod(myrank+nshift,mprocs)       !send to processor ntarget
        nsender=mod(myrank-nshift+mprocs,mprocs)   !receive from processor nsender
! buffer the particles to be sent:
        nsbuf=0;n=1 !# in the send buffer; counter that goes through particles on this proc
     4  continue
        if(n.gt.nraysp)goto 5
        if(iownerptcl(n).eq.ntarget)then
          nsbuf=nsbuf+1
          if(nsbuf.gt.nbufsize)then
            write(6,*)'send buffer overflow';ifail=-1;return
          endif
          sbuf(1:n1,nsbuf)=y(1:n1,n)
          y(1:n1,n)=y(1:n1,nraysp) !swap last good particle with location n
          iownerptcl(n)=iownerptcl(nraysp)
          nraysp=nraysp-1
          goto 4
        endif
        n=n+1
        goto 4
    5   continue
! send and receive data:
        call MPI_IRECV(rbuf,n1*nbufsize,MPI_DOUBLE_PRECISION,nsender,nshift,MPI_COMM_WORLD,msid(1),ierr) !tag value is nshift
        call MPI_ISEND(sbuf,n1*nsbuf,MPI_DOUBLE_PRECISION,ntarget,nshift,MPI_COMM_WORLD,msid(2),ierr)
        call MPI_WAIT(msid(1),mpistat,ierr)
        call MPI_GET_COUNT(mpistat,MPI_DOUBLE_PRECISION,nrecv,ierr)
        call MPI_WAIT(msid(2),mpistat,ierr)
!! for some reason, under some conditions the following causes the code to freeze:
!!      call MPI_SEND(sbuf,n1*nsbuf,MPI_DOUBLE_PRECISION,ntarget,nshift,MPI_COMM_WORLD,ierr)
!!      call MPI_RECV(rbuf,n1*nbufsize,MPI_DOUBLE_PRECISION,nsender,nshift,MPI_COMM_WORLD,mpistat,ierr) !tag value is nshift
!!      call MPI_GET_COUNT(mpistat,MPI_DOUBLE_PRECISION,nrecv,ierr)
        nrecv=nrecv/n1
        if(nrecv.eq.0)cycle
! store the received data:
        if(nrayspnew+nrecv.gt.n2)then
          write(6,'("(not enough storage)",1x,5(i9,1x))')myrank,nrayspnew,nrecv,nrayspnew+nrecv,n2
        endif
        ynew(1:n1,nrayspnew+1:nrayspnew+nrecv)=rbuf(1:n1,1:nrecv)
        nrayspnew=nrayspnew+nrecv
      enddo
      nraysp=nrayspnew
      y(1:n1,1:nraysp)=ynew(1:n1,1:nraysp)
      call MPI_REDUCE(nraysp,nraysgbl,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      return
      end
!
Subroutine dmrgrnk(XDONT, IRNGT, sizexdont,sizeirngt)
      implicit none
      integer :: sizexdont,sizeirngt
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      real(8), dimension (*), intent (in) :: XDONT     !changed (:) to (*)
      integer, dimension (*), intent (out) :: IRNGT   !changed (:) to (*)
! __________________________________________________________
      real(8) :: XVALA, XVALB
!
!!!!!!Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      integer, dimension (sizeirngt) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
!     write(6,*)'size(xdont)=',size(xdont)
!     write(6,*)'size(irngt)=',size(irngt)
!!!!!!NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      NVAL = Min (sizexdont, sizeirngt)
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine dmrgrnk
