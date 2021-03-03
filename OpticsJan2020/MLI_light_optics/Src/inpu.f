************************************************************************
* header:                 INPUT AND OUTPUT                             *
* Input and output for maps, matrices, files, arrays, and parameter    *
* sets                                                                 *
************************************************************************
c
      subroutine cf(p)
c  subroutine to close files
c  Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      dimension p(6),ip(6)
      include 'files.inc'
c
c set up control indices
      do 10 j=1,6
      ip(j)=nint(p(j))
   10 continue
c
c close indicated files
      do 20 j=1,6
      n=ip(j)
      if( n.gt.0) close(unit=n, err=30)
      go to 20
   30 write(jof,*) 'error in closing file unit ',n
   20 continue
c
      return
      end
c
***********************************************************************
c
      subroutine fnamechk(fname,iu,ierr,str)
cryne
c R. Ryne 7/16/2002
c check file names and associated unit number
c this is needed now that it is possible to open files w/ given names
c NB: At the end of this routine, only PE0 knows the file unit number iu,
c and only PE0 has an connection to file fname via file unit number iu.
c All PEs know the error flag ierr.
      use parallel
      logical ex,op,exu,opu,lname
      character*16 fname,uname
      character (len=*) str
      integer n
      common/bfileinfo/ifileinfo
!cryne===== 28 July, 2004      if(idproc.ne.0)return
      if(idproc.eq.0)then
!cryne=====
!
!
!     write(6,*)'here I ama in fnamemchk, iu=',iu
      nmax=len_trim(str)
      ierr=0
cryne----- Jan 3, 2003:
c if the requested unit number =0, then the code should obtain
c a valid, unused unit number.
c For historical purposes (original Marylie used multiple files <~20),
c the unit numbers will be searched between 31 on upward.
c This should be irrelevant to the user, as these numbers will not be seen.
      if(iu.eq.0)then
        iumin=31
        iumax=1000
        do iuu=iumin,iumax
          INQUIRE(UNIT=iuu,EXIST=exu,OPENED=opu,NAMED=lname,NAME=uname)
          if(exu .and. .not.opu)then
            iu=iuu
            exit
          endif
        enddo
        if(iu.eq.0)then
          if(idproc.eq.0)write(6,*)'unable to find an available unit #'
          ierr=1
          goto 9999
        endif
      endif
cryne-----
      INQUIRE(UNIT=iu,EXIST=exu,OPENED=opu,NAMED=lname,NAME=uname)
      INQUIRE(FILE=fname,EXIST=ex,OPENED=op,NUMBER=numb)
c if the file is already open and connected to the unit specified
c by the user (or the default value in use), then everything is OK:
      if(op.and.(numb.eq.iu))then
        if(idproc.eq.0.and.fname.ne.' ')then
          if(ifileinfo.ne.0)write(6,*)'(',str(1:nmax),') ',                 &
     &    'file ',fname(1:16), ' will stay connected to unit ',iu
        endif
        goto 9999
      endif
c the file needs to be created and/or opened.
c first check to see if iu is a valid unit number:
      if(.not.exu)then
       write(6,*)'Error: file unit ',iu,' is not valid for this system.'
       write(6,*)'It cannot be assigned to file ',fname
       write(6,*)'command ',str(1:nmax),' will be ignored'
       ierr=1
       goto 9999
      endif
c does the file exist?
      if(.not.ex)then
c if not, create/open it:
        if(fname.ne.' ')then
          if(idproc.eq.0)then
           if(ifileinfo.ne.0)write(6,*)'(',str(1:nmax),') ','creating ',    &
     &     'file ',fname(1:16), ' and connecting to unit ',iu
          endif
          open(unit=iu,file=fname,status='new')
        else
cg95      open(unit=iu,status='new')
          if(idproc.eq.0)write(6,*)'code should not get here'
          call myexit
        endif
        goto 9999
      endif
c file exists. Is it open?
c if not, open it, but first make sure unit iu is not in use.
      if(.not.op)then
        if(.not.opu)then
          if(idproc.eq.0.and.ifileinfo.ne.0)then
            write(6,*)'(',str(1:nmax),') ','opening ',                     &
     &     'existing file ',fname(1:16), ' and connecting to unit ',iu
            write(6,*)'If written to, the file will be overwritten'
          endif
          open(unit=iu,file=fname,status='old')
          goto 9999
        else
          write(6,*)'file unit ',iu,' is already in use.'
          write(6,*)'It cannot be assigned to file ',fname
c         write(6,*)'command ',str(1:nmax),' will be ignored'
c         ierr=1
          write(6,*)'instead the following file will be used: ',uname
          goto 9999
        endif
      else
c the file exists and is open. There must be a problem with
c conflicting unit numbers [numb.ne.iu] Fix it:
c       write(6,*)'(',str(1:nmax),') Problem w/ file unit specification'
c       write(6,*)'for file with name ',fname
c       write(6,*)'File is already connected to unit ',numb
c       write(6,*)'but user has specified (or default is) unit ',iu
c       write(6,*)'Existing unit number will be used'
c
        if(ifileinfo.ne.0)then
        write(6,*)'(',str(1:nmax),') File named ',fname
        write(6,*)'will continue to stay connected to unit ',numb
        endif
        iu=numb
        goto 9999
      endif
c
!cryne===== 28 July 2004 only PE 0 executed the above code
 9999 continue
      endif
!cryne===== 28 July 2004 now all PEs execute the following:
      call ibcast(ierr)
cryne-abell: delete next four lines when read_egengrads works in ||
c     call ibcast(iu)
c     if (idproc.ne.0.and.iu.ne.0) then
c       open(unit=iu,file=fname,status='old')
c     endif
      return
      end
c
***********************************************************************
c
      subroutine mapin(nopt,nskp,h,mh)
c read nonzero matrix elements and monomials from file unit mpi
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision h,mh
      include 'files.inc'
      dimension mh(6,6),h(monoms)
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
c
      if(nopt.eq.0)goto 5
      rewind mpi
      write(jof,210)mpi
    5 continue
c     initialize arrays:
      do 10 j=1,monoms
   10 h(j)=0.d0
      do 20 k=1,6
      do 20 l=1,6
   20 mh(k,l)=0.
c
c  skip nskp maps:
      ns=nskp
   55 if(ns.eq.0)goto 100
   60 read(mpi,*)i,j,temp
      if(i.eq.6.and.j.eq.6)goto 80
      goto 60
   80 read(mpi,*)k,temp
      if(k.eq.monoms)goto 90
      goto 80
   90 ns=ns-1
      goto 55
c
c  now read in the map:
  100 continue
  160 read(mpi,*)i,j,temp
      mh(i,j)=temp
      if(i.eq.6.and.j.eq.6)goto 180
      goto 160
  180 read(mpi,*)k,temp
      h(k)=temp
      if(k.eq.monoms)goto 200
      goto 180
  200 continue
      write(6,205) mpi,nskp
  205 format(1x,'map read in from file ',i3,'; ',i3,
     &' record(s) skipped')
  210 format(1x,'file unit ',i2,' rewound')
      return
      end
c
***********************************************************************
c
      subroutine mapout(nopt,h,mh)
c  output present matrix and polynomials (nonzero values)
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1984
c   modified Oct 89 by Tom Mottershead to work without requiring
c   a prior map file.
c The parameter nopt is not currently used
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
c
c calling arrays
c
      double precision h,mh
      dimension mh(6,6),h(monoms)
c
c check to see if file assignment makes sense
      if(mpo .le. 0) return
c
      nout=mpo
c
      read(nout,*,end=5,err=19)dummy
    5 rewind(nout)
c
c count the # of lines in the file:
      ncount=0
   10 read(nout,*,end=15)dummy
      ncount=ncount+1
      goto 10
   15 continue
      write(6,*)'found ',ncount,' lines in file ',nout
      rewind(nout)
      do n=1,ncount
        read(nout,*)dummy
      enddo
      write(jof,*)'adding current map to file on unit ',nout
      goto 25
c
c     read error means map file does not exist, so start a new one
c
  19  write(jof,21) nout
  21  format(' starting new map file on unit',i3)
c
c  write out map
c
  25  continue
      do 30 i=1,6
      do 30 j=1,6
   30 if(mh(i,j).ne.0.)write(nout,*)i,j,mh(i,j)
      if(mh(6,6).eq.0.)write(nout,*)6,6,mh(6,6)
      do 40 k=1,monoms
   40 if(h(k).ne.0.)write(nout,*)k,h(k)
      if(h(monoms).eq.0.)write(nout,*)monoms,h(monoms)
c      write(jof,100) nout
c  100 format(1x,'map written on file ',i2)
      return
      end
c
***********************************************************************
c
      subroutine mapsnd(iopt,nmap,ta,tm,ha,hm)
c this is a subroutine for sending a map to some buffer or file
c  Written by Alex Dragt, Spring 1987
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'buffer.inc'
      character*3 kynd
c
      dimension ta(monoms),ha(monoms)
      dimension tm(6,6),hm(6,6)
c
c procedure when iopt = 0
      if (iopt.eq.0) then
      if (nmap.ge.1 .and. nmap.le.5) then
      kynd='stm'
      call strget(kynd,nmap,ta,tm)
      endif
c
      if (nmap.eq.0) call mapmap(ta,tm,ha,hm)
c
      if (nmap.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmap.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmap.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmap.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmap.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
      endif
c
c procedure when iopt .gt. 0
      if (iopt.gt.0) then
      mpot=mpo
      mpo=iopt
      call mapout(0,ta,tm)
      mpo=mpot
      endif
c
      return
      end
c
************************************************************************
c
      subroutine numfile(p)
c This is a subroutine for numbering lines in a file
c Written by Alex Dragt, Spring 1987
      use rays
      include 'impli.inc'
      dimension p(6)
      write(6,*)'THIS ROUTINE (NUMFILE) NEEDS TO BE MODIFIED TO'
      write(6,*)'EXECUTE PROPERLY IN PARALLEL'
c
c set up control indices
      iopt=nint(p(1))
      nfile=nint(p(2))
      ifirst=nint(p(3))
      istep=nint(p(4))
c
c procedure for writing out a file with numbered lines
      line=ifirst
      do 10 i=1,nraysp
      write(nfile,100) line,
     & zblock(1,i), zblock(2,i), zblock(3,i),
     & zblock(4,i), zblock(5,i), zblock(6,i)
  100 format(1x,i5,6(1x,1pe11.4))
      line=line+istep
   10 continue
c
      return
      end
c
***********************************************************************
c
      subroutine of(p)
c  subroutine to open files
c  Written by Alex Dragt, Spring 1987
      include 'impli.inc'
      character*6 unit
      dimension p(6),ip(6)
      dimension unit(50)
      include 'files.inc'
c
c set up file names
c
      data (unit(i),i=1,50)/
     &'unit01','unit02','unit03','unit04','unit05',
     &'unit06','unit07','unit08','unit09','unit10',
     &'unit11','unit12','unit13','unit14','unit15',
     &'unit16','unit17','unit18','unit19','unit20',
     &'unit21','unit22','unit23','unit24','unit25',
     &'unit26','unit27','unit28','unit29','unit30',
     &'unit31','unit32','unit33','unit34','unit35',
     &'unit36','unit37','unit38','unit39','unit40',
     &'unit41','unit42','unit43','unit44','unit45',
     &'unit46','unit47','unit48','unit49','unit50'/
cryne 7/23/2002
      save unit
c
c set up control indices
      do 10 j=1,6
      ip(j)=nint(p(j))
   10 continue
c
c open indicated files
      do 20 j=1,6
      n=ip(j)
      if( n.gt.0 .and. n.le.50 .and. n.ne.lf .and. n.ne.jof
     & .and. n.ne.jodf) then
      open(unit=n, file=unit(n), status='unknown', err=30)
      endif
      go to 20
   30 write(jof,*) 'error in opening file unit ',n
   20 continue
c
      return
      end
c
***********************************************************************
c
      subroutine pcmap(n1,n2,n3,n4,fa,fm)
c  routine to print m,f3,f4 and t,u.
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
      use parallel
      use lieaparam, only : monoms
      include 'impli.inc'
      integer colme(6,0:6)
      include 'expon.inc'
      include 'pbkh.inc'
      include 'files.inc'
      dimension fa(monoms),fm(6,6)
      dimension t(monoms),u(monoms),u2(monoms)
      if(idproc.ne.0)return
c
c  test for matrix write
      if(n1.eq.0) goto 20
c
c  procedure for writing out matrix
c  write matrix at terminal
      if(n1.eq.1.or.n1.eq.3)then
        write(jof,13)
   13   format(/1h ,'matrix for map is :'/)
        write(jof,15)((fm(k,i),i=1,6),k=1,6)
   15   format(6(1x,1pe12.5))
        write(jof,*)
        write(jof,*)'nonzero matrix elements in full precision:'
        do k=1,6
        do i=1,6
        if(fm(k,i).ne.0.d0)write(jof,155)k,i,fm(k,i)
  155   format(i1,2x,i1,2x,1pg21.14)
        enddo
        enddo
        write(jof,*)
      endif
c  write matrix on file 12
      if(n1.eq.2.or.n1.eq.3)then
        write(jodf,13)
        write(jodf,15)((fm(k,i),i=1,6),k=1,6)
        write(jodf,*)
        write(jodf,*)'nonzero matrix elements in full precision:'
        do k=1,6
        do i=1,6
        if(fm(k,i).ne.0.d0)write(jodf,155)k,i,fm(k,i)
        enddo
        enddo
        write(jodf,*)
      endif
c
c  test for polynomial write
   20 continue
      if(n2.eq.0)goto 30
c
c  procedure for writing out polynomial
c  write polynomial at terminal
      if(n2.eq.1.or.n2.eq.3)then
        write(jof,22)
   22   format(/1h ,'nonzero elements in generating polynomial are :'/)
        do 25 i=1,monoms
ccccc   if(fa(i).eq.0.0d0) goto 25
        if(fa(i).eq.0.0d0) goto 25
        write(jof,27)i,(expon(j,i),j=1,6),fa(i)
cryne 6/21/2002   27   format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',d21.14)
   27 format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',1pg21.14)
   25   continue
      endif
c  write polynomial on file 12
      if(n2.eq.2.or.n2.eq.3)then
        write(jodf,22)
        do 26 i=1,monoms
c06jan2019  if(fa(i).eq.0.0d0) goto 26
        write(jodf,27)i,(expon(j,i),j=1,6),fa(i)
   26   continue
      endif
c
c  prepare for higher order matrix write if required
   30 continue
cryne write(6,*)'inside pcmap'
      if(n3.gt.0.or.n4.gt.0)then
cryne   write(6,*)'calling brkts'
        call brkts(fa)
cryne   write(6,*)'returned from brkts'
      endif
c
c  test for t matrix write
      if(n3.eq.0) goto 40
c
c  procedure for writing t matrix
c  write out heading
      if(n3.eq.1.or.n3.eq.3)write(jof,32)
      if(n3.eq.2.or.n3.eq.3)write(jodf,32)
   32   format(/1h ,'nonzero elements in second order matrix are :'/)
c  write out contents
        do 35 i=1,6
        call xform(pbh(1,i),2,fm,i-1,t)
        do 36 n=7,27
        if(t(n).eq.0.0d0) goto 36
        if(n3.eq.1.or.n3.eq.3)
     &  write(jof,38) i,n,i,(expon(j,n),j=1,6),t(n)
        if(n3.eq.2.or.n3.eq.3)
     &  write(jodf,38) i,n,i,(expon(j,n),j=1,6),t(n)
cryne6/21/02 38format(2x,'t',i1,'(',i3,')','=t',i1,'( ',3(2i1,1x),')=',d21.14)
   38 format(2x,'t',i1,'(',i3,')','=t',i1,'( ',3(2i1,1x),')=',1pg21.14)
   36   continue
   35   continue
c

c  test for u matrix write
   40 continue
      if(n4.eq.0) goto 50
c
c  procedure for writing u matrix
c  write out heading
      if(n4.eq.1.or.n4.eq.3)write(jof,42)
      if(n4.eq.2.or.n4.eq.3)write(jodf,42)
   42 format(/1h ,'nonzero elements in third order matrix are :'/)
c  write out contents
         do  44 i=1,6
         call xform(pbh(1,i),3,fm,i-1,u)
         call xform(pbh(1,i+6),3,fm,1,u2)
         do  45 n=28,83
         u(n)=u(n)+u2(n)/2.d0
         if(u(n).eq.0.0d0) goto 45
         if(n4.eq.1.or.n4.eq.3)
     &   write(jof,46) i,n,i,(expon(j,n),j=1,6),u(n)
         if(n4.eq.2.or.n4.eq.3)
     &   write(jodf,46) i,n,i,(expon(j,n),j=1,6),u(n)
cryne 6/21/02 46format(2x,'u',i1,'(',i3,')','=u',i1,'( ',3(2i1,1x),')=',d21.14)
   46 format(2x,'u',i1,'(',i3,')','=u',i1,'( ',3(2i1,1x),')=',1pg21.14)
   45    continue
   44    continue
c
c  procedure if all n's are zero or are faulty
   50 continue
c
      return
      end
c
***********************************************************************
c
      subroutine pset(p,k)
c  subroutine to read in parameter set values
c  the integer k labels the parameter set
c  Written by Alex Dragt, Fall 1986
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'parset.inc'
      dimension p(6)
c
c test on k
      if (k.gt.maxpst .or. k.le.0) then
      write(6,*) 'improper attempt to store parameters in
     & a parameter set with k=',k
      call myexit
      endif
c
c store parameter values
      do 10 i=1,6
   10 pst(i,k)=p(i)
c
      return
      end
c
***********************************************************************
c
      subroutine psrmap(n1,n2,fa,fm)
c  routine to print mij in the cartesian basis
c  and the f's in the static resonance basis.
c  Written by Alex Dragt, Fall 1986
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'srl.inc'
      include 'files.inc'
      dimension fa(monoms),fm(6,6)
c
c  beginning of routine
c
c  writing out matrix
      if(n1.eq.0)goto 20
      if(n1.eq.1.or.n1.eq.3)write(jof,33)
      if(n1.eq.2.or.n1.eq.3)write(jodf,33)
   33 format(/1h ,'matrix for map is :'/)
      if(n1.eq.1.or.n1.eq.3)write(jof,35)((fm(k,i),i=1,6),k=1,6)
      if(n1.eq.2.or.n1.eq.3)write(jodf,35)((fm(k,i),i=1,6),k=1,6)
c   35 format(6(1x,e12.5))
   35 format(6(1x,1pe12.5))
c
c writing out f's
   20 if(n2.eq.0)goto 30
      if(n2.eq.1.or.n2.eq.3)write(jof,55)
      if(n2.eq.2.or.n2.eq.3)write(jodf,55)
   55 format(/1h ,'nonzero elements in generating polynomial in',/,
     &1x,'the static resonance basis are :'/)
      do 22 i=1,monoms
      if (fa(i).eq.0.0d0) go to 22
      if(n2.eq.1.or.n2.eq.3)write(jof,60)i,sln(i),fa(i)
      if(n2.eq.2.or.n2.eq.3)write(jodf,60)i,sln(i),fa(i)
   60 format(1x,'f(',i3,')=f( ',a6,' )=',d21.14)
   22 continue
c
   30 continue
      return
      end
c
***********************************************************************
c
      subroutine pdrmap(n1,n2,fa,fm)
c  routine to print mij in the cartesian basis
c  and the f's in the dynamic resonance basis.
c    F. Neri  6/3/1986
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'drl.inc'
      dimension fa(monoms),fm(6,6)
c
c  beginning of routine
c
c  writing out matrix
      if(n1.eq.0)goto 20
      if(n1.eq.1.or.n1.eq.3)write(jof,33)
      if(n1.eq.2.or.n1.eq.3)write(jodf,33)
   33 format(/1h ,'matrix for map is :'/)
      if(n1.eq.1.or.n1.eq.3)write(jof,35)((fm(k,i),i=1,6),k=1,6)
      if(n1.eq.2.or.n1.eq.3)write(jodf,35)((fm(k,i),i=1,6),k=1,6)
c   35 format(6(1x,e12.5))
   35 format(6(1x,1pe12.5))
c
c  writing out f's
   20 if(n2.eq.0)goto 30
      if(n2.eq.1.or.n2.eq.3)write(jof,55)
      if(n2.eq.2.or.n2.eq.3)write(jodf,55)
   55 format(/1h ,'nonzero elements in generating polynomial in',/,
     &1x,'the dynamic resonance basis are :'/)
      do 22 i=1,monoms
      if (fa(i).eq.0.0d0) go to 22
      if(n2.eq.1.or.n2.eq.3)write(jof,60)i,dln(i),fa(i)
      if(n2.eq.2.or.n2.eq.3)write(jodf,60)i,dln(i),fa(i)
   60 format(1x,'f(',i3,')=f( ',a7,' )=',d21.14)
   22 continue
c
   30 continue
      return
      end
c
*****************************************************************
c
      subroutine pmif(iu,itype,fname)
c  subroutine to write out the master input file
c  written by Rob Ryne ca 1984
ctm  modified 9/01 to write unexpanded #include files
c
      use parallel
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'codes.inc'
      include 'incmif.inc'
      include 'files.inc'
c
      character*79 line
      logical noitem, nolab
      character*16 str1,fname
c
cryne 7/23/2002 don't know why this save statement was put here,
cryne but I will leave it for now.
       save
       if(idproc.ne.0)return
c------------------
c start routine
c------------------
ctm    itype = 0 to echo input file as is
ctm    itype = 1 to write input file from full internal commons.
ctm    itype = 2 to write only #labor part
ctm    itype = 3 to write from commons without expanding the #includes
c-------------------------------------------------------
c-------------------------------------------------------
c
c      initialize loop limits for full common dump
c
      noitem = .true.
      if(nb.gt.0) noitem = .false.
      nolab = .true.
      if(noble.gt.0) nolab = .false.
      kbgn = 1
      kend = na
      ibgn = 1
      iend = nb
      jbgn = 1
      jend = noble
      goto(4,90,30)itype
c
c  echo exact contents of file lf: (for out of bounds itype)
c
      rewind lf
    1 read(lf,1002,end=3)line
 1002 format(a)
      write(iu,1003)line
 1003 format(1x,a)
      goto 1
    3 continue
      write(iu,*)
      return
c
ctm     reset loop limits to exclude #includes
c
  30  continue
      if(ninc.eq.0) go to 4
      kend = na1
      noitem = .true.
      if((nb1.gt.0).and.(nb2.eq.nb)) then
         ibgn = 1
         iend = nb1
         noitem = .false.
      endif
      if((nb2.gt.0).and.(nb.gt.nb2)) then
         ibgn = nb2 + 1
         iend = nb
         noitem = .false.
      endif
      jbgn = noble2 + 1
c
c
c    write from internal commons
c
    4 continue
c
c  comments
      if(np.eq.0)goto 5
      write(iu,530)ling(1)
      write(iu,6)(mline(i),i=1,np)
c mline is character*80, but is written out with an a79 format
c so that the last character (which should be a blank) does not
c cause carriage returns (and hence blank lines) on laser printers
    6 format(a79)
c    6 format(a80)
c  beam
    5 write(iu,530)ling(2)
      write(iu,*)brho
      write(iu,*)gamm1
      write(iu,*)achg
      write(iu,*)sl
c  menu
      write(iu,530)ling(3)
cryne quick hack to check for long names in the menu:
      longnm=0
      do 10 k=kbgn,kend
         str1=lmnlbl(k)
         if(str1(9:9).ne.' ')longnm=1
         if(longnm.eq.0)write(iu,600)lmnlbl(k),ltc(nt1(k),nt2(k))
         if(longnm.eq.1)write(iu,6000)trim(lmnlbl(k)),ltc(nt1(k),nt2(k))
         imax=nrp(nt1(k),nt2(k))
         if(imax.eq.0)goto 10
         write(iu,603)(pmenu(i+mpp(k)),i=1,imax)
   10 continue
c
ctm  explicit #include
c
      if(itype.eq.3) then
        if(ninc.eq.0) go to 50
        write(iu,530)ling(8)
        do 35 jj = 1, ninc
        write(iu,33) incfil(jj)
  33    format(2x,a)
  35    continue
      endif
c  lines,lumps,loops
  50  if(noitem)goto 90
      do 40 ii=2,4
        write(iu,530)ling(ii+2)
        do 20 k=ibgn,iend
          if(ityp(k).ne.ii)goto 20
          write(iu,790)ilbl(k)
          if(longnm.eq.0)write(iu,791)(irep(l,k),icon(l,k),l=1,ilen(k))
          if(longnm.eq.1)                                               &
     &    write(iu,7910)(irep(l,k),trim(icon(l,k)),l=1,ilen(k))
   20   continue
   40 continue
c  labor
   90 if(nolab)goto 999
      write(iu,530)ling(7)
      do 100 j=jbgn,jend
  100 write(iu,800)num(j),latt(j)
c  529 format(1h ,'#comments')
  530 format(1h ,a8)
  600 format(1h ,1x,a8,1x,a8)
 6000 format(1h ,1x,a,1x,a8)
c  603 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  603 format((1h ,3(1x,1pg22.15)))
cryne  790 format(1h ,1x,a8)
  790 format(1h ,1x,a)
  791 format((1h ,1x,5(i5,'*',a8),1x,:'&'))
cryne 7910 format((1h ,1x,3(i5,'*',a28),1x,:'&'))
 7910 format((1h ,1x,5(i5,'*',a),1x,:'&'))
cryne  800 format(1h ,1x,i4,'*',a8)
  800 format(1h ,1x,i4,'*',a)
  999 continue
c      write(iu,*)
      return
      end
c
******************************************************************
c
      subroutine randin(nfile,kt1,kt2,p)
c  This is a subroutine for reading in and checking the parameters
c  for random elements.
c  Written by Alex Dragt, ca 1985, and modified by F. Neri
c  and Alex Dragt on 29 January 1988
      use lieaparam, only : monoms
      include 'impli.inc'
      double precision p(6)
      include 'files.inc'
      include 'codes.inc'
      include 'parset.inc'
c
      if(nfile.gt.0) then
c  Read in parameters from file:
   2    read(nfile,*,end=4,err=6)(p(i),i=1,nrp(kt1,kt2))
        goto 8
   4    rewind nfile
        goto 2
   6    write(6,7) nfile,kt1,kt2
   7    format(1x,'nfile=',i4,'kt1,kt2=',2i4,'error in randin')
        call myexit
   8    continue
      else if(nfile.lt.0) then
c  Read in parameters form parameter sets:
        npar = -nfile
        if(npar.gt.maxpst) then
          write(jof,*) ' num of pset too large in rnd lmnt.'
          call myexit
        endif
        do 7008 i=1,6
          p(i) = pst(i,npar)
 7008   continue
      else
        write(jof,*) ' zero pset number in rnd lmnt.'
        call myexit
      endif
c
c  Check on parameter values corresponding to control indices:
      goto(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     &170),kt2
      return
c
c  drift:
  10  continue
      return
c
c  normal entry bend:
  20  call rcheck(0,1,p(3),'lfrn','nbnd',nfile)
      call rcheck(0,1,p(4),'tfrn','nbnd',nfile)
      return
c
c  parallel faced bend:
  30  continue
      return
c
c  general bending magnet:
  40  call rcheck(0,1,p(5),'lfrn','gbnd',nfile)
      call rcheck(0,1,p(6),'tfrn','gbnd',nfile)
      return
c
c  leading or trailing rotation for a parallel faced bend:
  50  call rcheck(0,1,p(2),'kind','prot',nfile)
      return
c
c  body of a general bending magnet:
  60  continue
      return
c
c  hard edge fringe fields of a normal entry bend:
  70  call rcheck(0,1,p(2),'iedg','frng',nfile)
      return
c
c  combined function bend:
  80  continue
      return
c
c  quadrupole:
  90  call rcheck(0,1,p(3),'lfrn','quad',nfile)
      call rcheck(0,1,p(4),'tfrn','quad',nfile)
      return
c
c  sextupole:
 100  continue
      return
c
c  mag. octupole:
 110  continue
      return
c
c  elec. octupole:
 120  continue
      return
c
c  short rf cavity
 130  continue
      return
c
c  axial rotation:
 140  continue
      return
c
c  linear matrix via twiss parameters:
 150  call rcheck(1,3,p(1),'ipla','twsm',nfile)
      return
c
c  thin lens low order multipole:
 160  continue
      return
c
c  "compressed" low order multipole:
 170  continue
c
      return
      end
c
c
***************************************************************
c
      subroutine raysin(icfile,nraysinp,scaleleni,scalefrqi,scalemomi,   &
     &ubuf,jerr)
c  routine to read in initial conditions of rays to be traced
c  and to initialize the arrays istat and ihist
c  Written by Robert Ryne ca 1984, and modified by Alex Dragt
cryne 2/16/2002 mods to handle missing file and file with a bad record

      use rays
!     use parallel
      include 'impli.inc'
      include 'files.inc'
      include 'parset.inc'
      character*80 fname,line
      character*80 string
      character*16 symb(50),istrtsymb(50)
      character*16 ubuf
      logical keypres,numpres,leof
      dimension bufr(1)
      dimension tmpvec(6)
c
c procedure for reading a single set of initial conditions from
c a parameter set
c
      if(icfile.lt.0) then
cryne 3/10/2004 this option does not work sensibly on multiple processors.
cryne Whatever is in ipset goes in zblock, so all procs get the same partcles.
      ipset=-icfile
      if(ipset.gt.maxpst) then
      if(idproc.eq.0)write(jof,*)'(raysin)parameter icfile out of range'
      call myexit
      endif
      nrays=1
      do 10 j=1,6
   10 zblock(j,1)=pst(j,ipset)
      goto 130
      endif
c
      if(idproc.eq.0)then
c procedure for reading initial conditions from a file
c first check to see if there is a header (scaling info) in the input file:
        write(6,*)'Rewinding/reading particle data in file connected',    &
     &  ' to unit ',icfile
   20   continue
        nrays=0
        rewind icfile
        ncomm=0
   22   continue
        read(icfile,'(a)',end=110,err=120)line
cryne quick hack: % should be in column 1, or user might put an extra
cryne space or two:
        if(line(1:1).eq.'%'   .or.                                        &
     &     line(1:2).eq.' %'  .or.                                        &
     &     line(1:3).eq.'  %')then
          call low(line)
          write(6,'(a)')line
          ncomm=ncomm+1
          ubuf=' '
         call getscaleinfo(line,ubuf,scaleleni,scalefrqi,scalemomi,jerr)
         if(jerr.ne.0)write(6,*)'(RAYSIN)ERROR RETURN FROM GETSCALEINFO'
        else
          write(6,*)'found ',ncomm, 'lines of header info'
          goto 123
        endif
        goto 22
c EOF:
  110   continue
        write(6,*)'trouble:particle data file may be missing or empty.'
        goto 122
c ERR:
  120   write(jof,121)
  121   format(1x,'trouble in raysin while trying to read header info')
  122   continue
c       try opening another file:
        write(6,*)'type a file name or <cr> to stop'
        read(5,210)fname
        if(fname.eq.' ')call myexit
  210   format(a16)
        open(icfile,file=fname,status='old',err=300)
        goto 20
  300   continue
        write(6,*)'file still does not exist. Halting.'
        call myexit

  123   continue
c count particle data
        if(nraysinp.gt.0)then
        nrays=nraysinp
        goto 125
        endif
c
        nrays=0
        rewind icfile
        if(ncomm.ne.0)then
          do i=1,ncomm
cryne could put an end= and err= here, but not really needed.
cryne     read(icfile,'(a)',end=987,err=987)line
          read(icfile,'(a)')line
          enddo
        endif
cryne  222   read(icfile,'(a)',end=910,err=920)line
cryne modified March 10, 2004 to perform unformmated read of 6 numbers
cryne instead of a formatted read of a character string.
cryne this was done because there may be blank lines (e.g. at end), which
cryne lead to an incorrect count unless reading 6 numbers (not a char string)
  222   read(icfile,*,end=910,err=920)tmpvec(1:6)
        nrays=nrays+1
        goto 222
  910   continue
        if(nrays.ne.0)then
          write(6,*)'found ',nrays,' particles in data file'
          goto 125
        endif
        write(6,*)'trouble: no particles found in particle data file'
        call myexit
  920 continue
        write(6,*)'error trying to read particle data'
        call myexit
c Processor 0 knows value of nrays. Now broadcast it to all processors
  125   continue
        do l=1,nvp-1
        call MPI_SEND(nrays,1,mntgr,l,96,lworld,ierr)
        enddo
      else
        call MPI_RECV(nrays,1,mntgr,0,96,lworld,mpistat,ierr)
      endif
c
cryne 1 August, 2004
      if(nrays.lt.nvp)then
        if(idproc.eq.0)then
          write(6,*)'ERROR: # of particles is < # of processors'
          write(6,*)'Execution will be terminated'
        endif
        call myexit
      endif
c
      if( (nrays.le.maxray) .and. allocated(zblock) )goto 127
      if( (nrays.le.maxray) .and. .not.allocated(zblock) )then
        call new_particledata
        goto 127
      endif
c
c trouble: nrays.gt.maxray
c has zblock already been allocated? yes...
      if(allocated(zblock))then
        if(idproc.eq.0)then
          write(6,*)'error: zblock has been allocated but is not large'
          write(6,*)'enough to contain the particle array to be read in'
         write(6,*)'Rerun with maxray .ge. ',nrays,' in beam definition'
        endif
        call myexit
      endif
c ...no, zblock has not been allocated
      if(idproc.eq.0)then
        write(6,*)'WARNING: # of rays in input file = ',nrays
        write(6,*)'         current value of maxray = ',maxray
        write(6,*)'setting maxray to ',nrays,' and allocating zblock'
      endif
      maxray=nrays
      call new_particledata
c=========== data file counted; now proc 0 reads in the rays===========
  127 continue
      if(idproc.eq.0)then
        rewind icfile
        write(6,*)'PE0 has rewound file with unit number ',icfile
        if(ncomm.gt.0)then
          do n=1,ncomm
          read(icfile,*)line(1:80)
          enddo
        endif
cryne May 19 2006        nraysp=(nrays-1)/nvp + 1
        nraysp0=nrays/nvp
        nremain=nrays-nraysp0*nvp
        nraysp=nraysp0
        if(nremain.gt.0.and.idproc.le.nremain-1)nraysp=nraysp0+1
c1221   continue
c       if(nraysp*(nvp-1).ge.nrays)then
c         nraysp=nraysp-1
c         if(nraysp.eq.0)then
c         write(6,*)'trouble: nraysp=0. something wrong. halting.'
c         endif
c         goto 1221
c       endif
c       write(6,*)'nraysp=',nraysp
        do k=1,nraysp
c       write(6,*)k
        read(icfile,*,end=9991,err=9992)zblock(:,k)
c       write(6,1841)zblock(1:6,k)
c1841   format(6(1pe12.5,1x))
        enddo
        write(6,*)'PE0 has read its own data'
cryne May 19, 2006        nraysw=nraysp
        do l=1,nvp-1
cryne May 19, 2006        if (l.eq.nvp-1) nraysw = nrays - nraysp*(nvp-1)
        nraysw=nraysp0
        if(nremain.gt.0.and.l.le.nremain-1)nraysw=nraysp0+1
c       write(6,*)'PE0 is reading nraysw=',nraysw,' for PE#',l,'...'
        do k=1,nraysw
        read(icfile,*,end=9991,err=9992)tblock(:,k)
        enddo
c       write(6,*)'...done'
        if(l.lt.nvp)then
          call MPI_SEND(tblock,6*nraysw,mreal,l,99,lworld,ierr)
        endif
        enddo
c
        izero=0
        write(6,*)'processor ',izero,' retained ',nraysp,' rays'
c
      else
        call MPI_RECV(zblock,6*maxrayp,mreal,0,99,lworld,mpistat,ierr)
        call MPI_GET_COUNT(mpistat,mreal,nraysp,ierr)
        nraysp=nraysp/6
        write(6,*)'processor ',idproc,' received ',nraysp,' rays'
      endif

  130 continue
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nraysp
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
      return
cryne This does not work properly, since only PE0 would get here and exit.
cryne To do this right all the procs should find out that there is a problem,
cryne and all should exit. E.G. PE0 could send a special number in tblock(1,1),
cryne and all the procs could check to see if they received it, and if so, exit
cryne Too much trouble. MPI is a pain in the neck.
 9991 continue
      if(idproc.eq.0)write(6,*)'error reading particle data file'
      call myexit
 9992 continue
      if(idproc.eq.0)write(6,*)'EOF encountered reading ptcl data file'
      call myexit
      end
***********************************************************************
c
      subroutine rcheck(min,max,arg,ivar,itype,nfile)
c  This is a subroutine for checking that control
c  parameters for random elements lie within the allowed range.
c  Written by Alex Dragt ca 1985
      double precision arg
      character*4 ivar,itype
      iarg=nint(arg)
      if (iarg.lt.min.or.iarg.gt.max) goto 10
      return
   10 write (6,100) ivar,iarg,itype,nfile
  100 format(1x,a4,'=',i4,1x,'trouble with random element',1x,a4,1x,
     &'read from file',1x,i4)
      call myexit
      return
      end
c
************************************************************************
c
      subroutine wcl(pp)
c  subroutine to write out contents of a loop
c  Written by Alex Dragt, 23 August 1988
c  Modified 19 June 1998 AJD
c  Based on the subroutines cqlate and pmif
c
      use beamdata
      use acceldata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c
c common blocks
c
      include 'codes.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
c
      dimension pp(6)
c
c local variables
c
      character*8 string(5),str
      logical     ljof,ljodf
c
c  set up control indices
c
      iopt=nint(pp(1))
      ifile=nint(pp(2))
      isend=nint(pp(3))
c
c  start routine
c
c  see if a loop exists
      if(nloop.le.0) then
      write(jof ,*) ' error from wcl: no loop has been specified'
      write(jodf,*) ' error from wcl: no loop has been specified'
      return
      endif
c
c  procedure when iopt=1 (write only names of loop contents
c  at the terminal and/or file 12)
      if (iopt.eq.1) then
      ljof  = isend.eq.1 .or. isend.eq.3
      ljodf = isend.eq.2 .or. isend.eq.3
      if (ljof .or. ljodf) then
c  write loop name
      if(ljof ) write(jof ,510) ilbl(nloop), joy
      if(ljodf) write(jodf,510) ilbl(nloop), joy
  510 format(/,1h ,'contents of loop ',a8,';',i6,' items:')
c  write loop contents
      do 1 jtot=0,joy,5
       kmax=min(5,joy-jtot)
       do 2 k=1,kmax
        jk1=jtot+k
c element
        if(mim(jk1).lt.0) then
          string(k)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(k)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(k)=ilbl(inuse(mim(jk1)))
        endif
   2  continue
      if(ljof)  write(jof ,511)(string(k),k=1,kmax)
      if(ljodf) write(jodf,511)(string(k),k=1,kmax)
  511 format(' ',5(1x,a8))
   1  continue
      endif
      endif
c
c  Procedure when iopt=2 (write out names of loop contents with & signs
c  on file ifile.  Each line is preceded by a blank space.)
      if (iopt .eq. 2 .and. ifile .gt. 0) then
c  write loop contents
      do  jtot=0,joy,5
       kmax=min(5,joy-jtot)
       jcheck=joy-jtot
       do  k=1,kmax
        jk1=jtot+k
c element
        if(mim(jk1).lt.0) then
          string(k)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(k)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(k)=ilbl(inuse(mim(jk1)))
        endif
      end do
      if (jcheck .gt. 5) write(ifile,612) (string(k),k=1,kmax)
      if (jcheck .le. 5) write(ifile,613) (string(k),k=1,kmax)
  612 format(' ',' ',5(1x,a8),' &')
  613 format(' ',' ',5(1x,a8))
      end do
      endif
c
c  procedure when iopt=3 (write out names of loop contents
c  with % and & signs on file ifile)
      if (iopt .eq. 3 .and. ifile .gt. 0) then
c  write loop contents
      do  jtot=0,joy,5
       kmax=min(5,joy-jtot)
       jcheck=joy-jtot
       do  k=1,kmax
        jk1=jtot+k
c element
        if(mim(jk1).lt.0) then
          string(k)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(k)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(k)=ilbl(inuse(mim(jk1)))
        endif
      end do
      if (jcheck .gt. 5) write(ifile,712) (string(k),k=1,kmax)
      if (jcheck .eq. 5) write(ifile,713) (string(k),k=1,kmax)
      if (jcheck .lt. 5 .and. kmax .gt. 0) then
      write(ifile,714) (string(k),k=1,kmax)
      endif
  712 format(' ',' ',5(1x,'% ',a8),' &')
  713 format(' ',' ',5(1x,'% ',a8),' %')
  714 format(' ',' ',5(1x,'% ',a8))
      end do
      endif
c
c  procedure when iopt=4 or iopt=5 (write out full loop contents)
      if (ifile .gt. 0 .and. (iopt.eq.4 .or. iopt.eq.5)) then
c  comments
      write(ifile,530) ling(1)
c  write loop name
      write(ifile,512) ilbl(nloop)
  512 format(1h ,' contents of loop ',a8)
c  beam
      write(ifile,530) ling(2)
      write(ifile,*) brho
      write(ifile,*) gamm1
      write(ifile,*) achg
      write(ifile,*) sl
c  biglist heading
      write(ifile,127)
  127 format(1h ,'#biglist')
c
c  contents of biglist
c  write loop contents
      do 137 jk1=1,joy
c element
        if(mim(jk1).lt.0) then
          string(1)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(1)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(1)=ilbl(inuse(mim(jk1)))
        endif
      call lookup(string(1),itype,item)
c      write(6,513) string(1)
c  513 format(1x,a8)
c      write(6,*) 'itype and item are ',itype, item
c procedure for a menu item
      if(itype.eq.1) then
      k=item
c case where iopt=4
      if (iopt.eq.4) then
      write(ifile,600)lmnlbl(k),ltc(nt1(k),nt2(k))
  600 format(1h ,1x,a8,1x,a8)
         imax=nrp(nt1(k),nt2(k))
         if(imax.eq.0)goto 137
      write(ifile,603)(pmenu(i+mpp(k)),i=1,imax)
c  603 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  603 format((1h ,3(1x,1pg22.15)))
      endif
c case where iopt=5
      if (iopt.eq.5) then
      imax=nrp(nt1(k),nt2(k))
      write(ifile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
  605 format(1h ,1x,a8,1x,a8,1x,i5,1x,i5,1x,i5)
         if(imax.eq.0)goto 137
      write(ifile,607)(pmenu(i+mpp(k)),i=1,imax)
c  607 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  607 format((1h ,3(1x,1pg22.15)))
      endif
      endif
c procedure for a lump
      if(itype.eq.3) then
c case where iopt=4
      if (iopt .eq. 4) then
      write(ifile,514) string(1)
  514 format(1x,1x,a8,1x,'lump')
      endif
c case where iopt=5
      if (iopt .eq. 5) then
      write(ifile,515) string(1),mim(jk1)
  515 format(1x,1x,a8,1x,'lump',9x,'0',4x,'-1',2x,i4)
      endif
      endif
  137 continue
      endif
c
  530 format(1h ,a8)
c
      return
      end
c
************************************************************************
c
      subroutine whst(p)
c  subroutine to write out history of beam loss
c  Written by Alex Dragt, Fall 1986
      use rays
      include 'impli.inc'
      dimension p(6)
      write(6,*)'THIS ROUTINE (WHST) NEEDS TO BE MODIFIED TO'
      write(6,*)'EXECUTE PROPERLY IN PARALLEL'
c  begin routine
      ifile=nint(p(1))
      job=nint(p(2))
c  determine what job is to be done
      if (job.eq.2) goto 200
c  procedure for writing out istat
c
      do 5 k=1,nraysp
      write (ifile,50) k, istat(k)
   50 format (1h ,i10,i10)
    5 continue
      if(nraysp.gt.0)then
         write(6,*)'istat written on file ',ifile,' for PE# ',idproc
      endif
      return
c
c  procedure for writing out ihist
c
  200 continue
      do 10 k=1,nlost
      write (ifile,100) k, ihist(1,k), ihist(2,k)
  100 format (1h ,i10,i10,i10)
   10 continue
      write(6,*) 'nlost =',nlost
      if(nlost.eq.0) write(6,*) 'ihist not written out'
      if(nlost.gt.0) write(6,*) 'ihist written on file ',ifile
      return
      end
c
***********************************************************************
c
      subroutine wps(ipset,isend)
c subroutine for writing out values in a parameter set
c Written by Alex Dragt, 30 January 1988
c
      include 'impli.inc'
      include 'files.inc'
      include 'parset.inc'
c
c Check to see that ipset is within range
      if ((ipset.lt.1) .or. (ipset.gt.maxpst)) then
        write (jof,*) 'WARNING: ipset out of range in command',
     &  ' with typecode wps'
        return
      endif
c
c Write out values of parameters
c
      if ((isend.eq.1) .or. (isend.eq.3)) then
      write (jof,*) 'values of parameters in the parameter set',ipset
      write (jof,100)( pst(j,ipset), j=1,6)
      endif
      if ((isend.eq.2) .or. (isend.eq.3)) then
      write (jodf,*) 'values of parameters in the parameter set',ipset
      write (jodf,100)( pst(j,ipset), j=1,6)
      endif
  100 format((1h ,3(1x,1pg22.15)))
c
      return
      end
c
c--------------
c
      subroutine getfirststring(line,m1,mlast,istart,iend,ierr)
      implicit none
      character*80 line
      integer m1,mlast,istart,iend,ierr
      integer m
      ierr=0
      do m=m1,mlast
      if(line(m:m).ne.' ')then
        istart=m
        goto 100
      endif
      enddo
      ierr=1
      return
  100 continue
      do m=istart+1,mlast
      if(line(m:m).eq.' ')then
        iend=m-1
        return
      endif
      enddo
      iend=80
      return
      end
c
c--------------
c
      subroutine getscaleinfo(line,ubuf,scaleleni,scalefrqi,scalemomi,   &
     &jerr)
      implicit none
c arguments:
      character*80 line
      character*16 ubuf
      real*8 scaleleni,scalefrqi,scalemomi
      integer jerr,iostat,numerr
c other variables:
      real*8 bufr(1),value
      character*80 string
      integer m0,istart,iend,ierr
cryne!!!!!!!!!!!!!!!!!!!!!!!!!!!!!July 5, 2004 needs cleaning up!!!
      ierr=0
cryne!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c units:

          m0=index(line,'units')
          if(m0.ne.0)then
            write(6,*)'found units'
            m0=index(line,'=')
            call getfirststring(line,m0+1,80,istart,iend,ierr)
            if(ierr.eq.0)then
c             write(6,*)'found istart,iend=',istart,iend
c             write(6,*)'char value is ',line(istart:iend)
              ubuf=line(istart:iend)
            else
              write(6,*)'error return from getfirststring'
              jerr=1
            endif
c           write(6,*)'character string following UNITS is: ',ubuf
            goto 23
          endif
c scale length:
          m0=index(line,'length')
          if(m0.ne.0)then
            write(6,*)'found scale_length'
            m0=index(line,'=')
            call getfirststring(line,m0+1,80,istart,iend,ierr)
            if(ierr.eq.0)then
c             write(6,*)'found istart,iend=',istart,iend
c             write(6,*)'char value is ',line(istart:iend)
              string=line(istart:iend)
            else
              write(6,*)'error return from getfirststring'
              jerr=1
            endif
            read(string,*,iostat=numerr)value
            if(numerr.eq.0)then
              bufr(1)=value
            else
              write(6,*)'ERROR reading SCALE LENGTH info'
            endif
            scaleleni=bufr(1)
            write(6,*)'scale length in particle data file is ',scaleleni
          goto 23
          endif
c scale frequency:
          m0=index(line,'freq')
          if(m0.ne.0)then
            write(6,*)'found scale_frequency'
            m0=index(line,'=')
            call getfirststring(line,m0+1,80,istart,iend,ierr)
            if(ierr.eq.0)then
c             write(6,*)'found istart,iend=',istart,iend
c             write(6,*)'char value is ',line(istart:iend)
              string=line(istart:iend)
            else
              write(6,*)'error return from getfirststring'
              jerr=1
            endif
            read(string,*,iostat=numerr)value
            if(numerr.eq.0)then
              bufr(1)=value
            else
              write(6,*)'ERROR reading SCALE FREQENCY info'
            endif
            scalefrqi=bufr(1)
            write(6,*)'scale freq in particle data file is ',scalefrqi
            goto 23
          endif
c scale momentum:
          m0=index(line,'momentum')
          if(m0.ne.0)then
            write(6,*)'found scale_momentum'
            m0=index(line,'=')
            call getfirststring(line,m0+1,80,istart,iend,ierr)
            if(ierr.eq.0)then
c             write(6,*)'found istart,iend=',istart,iend
c             write(6,*)'char value is ',line(istart:iend)
              string=line(istart:iend)
            else
              write(6,*)'error return from getfirststring'
              jerr=1
            endif
            read(string,*,iostat=numerr)value
            if(numerr.eq.0)then
              bufr(1)=value
            else
              write(6,*)'ERROR reading SCALE MOMENTUM info'
            endif
            scalemomi=bufr(1)
            write(6,*)'scale momentum in input file is ',scalemomi
            goto 23
          endif
   23 continue
      jerr=ierr
      return
      end

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c==============================================================================
c//////////////////////////////////////////////////////////////////////////////

      subroutine wrtmap(iu,refinit,reffin,arclen,h,mh)
c------------------------------------------------------------------------------
c
c  wrtmap:  writes a map to file, including the scale factors used to generate
c           the map, the initial and final reference particle trajectories, the
c           arclength of the map, and the matrix and non-linear parts of the
c           map.  Map is written in internal units, and the scale factors
c           allow one to convert from these units to any other units by
c           rescaling.  Also, only the non-zero components of the map are
c           written to file.
c
c  input:
c     iu :        file unit number
c     refinit :   initial reference trajectory
c     reffin :    final reference trajectory
c     arclen :    length of map element
c     h :         non-linear part of map
c     mh :        linear (matrix) part of map
c
c  misc:
c     scale factors for length, momentum, and frequency are stored in the
c     beamdata module found in 'afro_mod.f90'
c
c  KMP - 9 Nov 2006
c------------------------------------------------------------------------------
      use parallel, only : idproc
      use lieaparam, only : monoms
      use beamdata, only : sl,p0sc,freqscl
      implicit none
      double precision arclen,refinit(6),reffin(6),mh(6,6),h(monoms)
      integer i,j,iu
c
c
c Begin Map File delimiter:
      write(iu,'(a)') '#BeginMap'
c
c Scale Factors:
      write(iu,100) '#ScaleLen=',sl
      write(iu,100) '#ScaleMom=',p0sc
      write(iu,100) '#ScaleFrq=',freqscl
 100  format(A10,1X,ES21.13E3)
c
c Arc Length of Map Element:
      write(iu,150) '#ArcLen=',arclen
 150  format(A8,1X,ES21.13E3)
      write(iu,160) '#Monoms=',monoms
 160  format(A8,1X,I8)
c
c Initial/Final Reference Trajectory:
      write(iu,200) '#InitialRefTraj=',refinit
      write(iu,200) '#FinalRefTraj=  ',reffin
 200  format(A16,6(1X,ES21.13E3))
c
c Matrix/Linear part of Transfer Map:
      do i=1,6
         do j=1,6
            if(abs(mh(i,j)).gt.0.)then
               write(iu,300) i,j,mh(i,j)
 300           format(I1,1X,I1,1X,ES21.13E3)
            endif
         enddo
      enddo
c
c Non-linear part of Transfer Map (output 'monoms' in case it changes at some point)
      do i=28,monoms
         if(abs(h(i)).gt.0.)then
            write(iu,400) i,h(i)
 400        format(I8,1X,ES21.13E3)
         endif
      enddo
c
c End Map File delimiter:
      write(iu,'(a)') '#EndMap'
c
c Automatically flush write to file so map can be read and used immediately
      call myflush(iu)
      end

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c==============================================================================
c//////////////////////////////////////////////////////////////////////////////

      subroutine rdmap(iu,darc,dreftraj,h,mh)
c------------------------------------------------------------------------------
c
c  readmap: reads a map, or series of maps, from file, including the scale
c           factors used to generate the map, the initial and final reference
c           particle trajectories, the arclength of the map, and the matrix
c           and non-linear parts of the map.  The map is read in assuming
c           internal units, and the scale factors allow one to convert from
c           these units to any other units by rescaling.  The map, once read,
c           is immediately applied.
c
c  input:
c     iu :        file unit number
c
c
c  misc:
c     scale factors for length, momentum, and frequency are stored in the
c     beamdata module found in 'afro_mod.f90'
c
c  KMP - 15 Dec 2006
c------------------------------------------------------------------------------
      use parallel, only : idproc
      use lieaparam, only : monoms
      use beamdata, only : sl,p0sc,freqscl
      implicit none
      include 'map.inc'
      double precision refinit(6),reffin(6),dreftraj(6)
      double precision mh(6,6),h(monoms),mh2(6,6),h2(monoms)
      double precision lsctmp,psctmp,fsctmp,darc
      double precision d1,d2
      integer i,j,iu,monomstmp
      character*160 cline
      character*20 cword
      character ch
c
c Read a line from the file
 100  read(iu,'(A160)',end=300) cline
c      write(*,*) cline
c
c Look for #BeginMap delimiter
      if(cline(1:9).eq.'#BeginMap')then
         read(iu,*) cword, lsctmp
         if(lsctmp.eq.0.0)then
            write(*,*) 'ERROR (rdmap): length scale is zero!'
            write(*,*) '      keeping current length scale'
            lsctmp = sl
         endif
c         write(*,*) cword, lsctmp
         read(iu,*) cword, psctmp
         if(psctmp.eq.0.0)then
            write(*,*) 'ERROR (rdmap): momentum scale is m*c'
            write(*,*) '      keeping current momentum scale'
            psctmp = p0sc
         endif
c         write(*,*) cword, psctmp
         read(iu,*) cword, fsctmp
         if(fsctmp.eq.0.0)then
            write(*,*) 'ERROR (rdmap): frequency scale is zero!'
            write(*,*) '      keeping current frequency scale'
            fsctmp = freqscl
         endif
c         write(*,*) cword, fsctmp
         read(iu,*) cword, darc
c         write(*,*) cword, darc
         read(iu,*) cword, monomstmp
c         write(*,*) cword, monomstmp
         if(monoms.ne.monomstmp)then
            if(idproc.eq.0)then
               write(*,*) 'WARNING (rdmap): Lie Algebraic order of ',
     &                    'map read from file is ',monomstmp,' but '
               write(*,*) '    the current order in the simulation ',
     &                    'is ',monoms,'.'
            endif
         endif
         read(iu,*) cword, refinit
c         write(*,*) cword, refinit
c
c Rescale initial reference trajectory
         refinit(1) = refinit(1) * lsctmp / sl
         refinit(2) = refinit(2) * psctmp / p0sc
         refinit(3) = refinit(3) * lsctmp / sl
         refinit(4) = refinit(4) * psctmp / p0sc
         refinit(5) = refinit(5) * freqscl / fsctmp
         refinit(6) = refinit(6) * (fsctmp*lsctmp*psctmp)
     &                           / (freqscl*sl*p0sc)
c
c Compare initial and current reference trajectories
         do i=1,6
            if(i.ne.5)then
            d1 = reftraj(i) - refinit(i)
            d2 = abs(reftraj(i)) + abs(refinit(i))
            if((d2.gt.0.0).and.(abs(d1).gt.0.05*d2))then
               if(idproc.eq.0)then
               write(*,*) 'WARNING (rdmap): Map reference trajectory ',
     &                    'coordinate',i,'differs from simulation '
               write(*,*) 'reference trajectory by greater than 10%!',
     &                    ' Results could be inaccurate!'
               endif
            endif
            endif
         enddo
         read(iu,*) cword, reffin
c         write(*,*) cword, reffin
c Rescale final reference trajectory
         reffin(1) = reffin(1) * lsctmp / sl
         reffin(2) = reffin(2) * psctmp / p0sc
         reffin(3) = reffin(3) * lsctmp / sl
         reffin(4) = reffin(4) * psctmp / p0sc
         reffin(5) = reffin(5) * freqscl / fsctmp
         reffin(6) = reffin(6) * (fsctmp*lsctmp*psctmp)
     &                           / (freqscl*sl*p0sc)
c
c Find change in reference trajectory
         dreftraj = reffin - refinit
c
c Initialize map and read in non-zero coefficients
         h = 0.0
         mh = 0.0
 200     read(iu,'(A160)') cline
c         write(*,*) '1: ',cline
         if(cline(1:1).eq.'#')then
            if(cline(1:7).eq.'#EndMap')goto 400
         else
            read(cline,*) i
c            write(*,*) '2: ',i
            if(i.gt.6)then
               read(cline,*) i,d1
c               write(*,*) '3: ',i,d1
               if(i.le.monoms)h(i) = d1
            else
               read(cline,*) i,j,d1
c               write(*,*) '4: ',i,j,d1
               mh(i,j) = d1
            endif
         endif
         goto 200
      else
         goto 100
      endif
c
c Error trapping: Only reaches here if no map found in file.
 300  continue
      write(*,*) 'ERROR (rdmap): No map found in file.  ',
     &           'Returning the identity map.'
      call ident(h,mh)
      darc = 0.0
      dreftraj = 0.0
      return
c
 400  continue
c Scale map to with current scale factors
      call rescale_map(lsctmp,psctmp,fsctmp,h,mh)
c
c Check return values:
!       if(idproc.eq.0)then
!          write(*,*) 'ERROR CHECKING:'
!          write(*,*) 'darc = ',darc
!          write(*,*) 'mh:'
!          do i=1,6
!             write(*,*) '   ',(mh(i,j),j=1,6)
!          enddo
!          write(*,*) 'h:'
!          do i=1,monoms
!             if(h(i).ne.0.0)write(*,*) '   ',i,': ',h(i)
!          enddo
!       endif
      return
      end

