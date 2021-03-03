      subroutine dumpin
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use parallel, only : idproc
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
c
cryne 5/4/2006 added this to allowing printing of messages from #labor:
      integer, parameter :: lmaxmsg=1000
      character*256 lattmsg(lmaxmsg)
      common/lattmsgb/lattmsg,lmsgpoi
c
cryne 5/4/2006 Added MADX capability. Initialized to MAD8 below.
c   MADX is enabled by putting ";" on the #comment line after #comment
      logical MAD8,MADX
      common/mad8orx/MAD8,MADX
c
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*16 strarr(40),string
c     character*80 line
      character*256 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
c-----------------------------------------------------------------------
      lmsgpoi=0
      MAD8=.true.
      MADX=.false.
c Here are some things that need to be initialized before the data are read in:
      isymbdef=-1
c default units are "static" (="magnetic") with scale length=1 and omega*l/c=1
      sl=1.d0
      clite=299792458.d0
      ts=sl/clite
      omegascl=1.d0/ts
      freqscl=omegascl/( 4.d0*asin(1.d0) )
      p0sc=0.d0
      magunits=1
      iverbose=0
c
      slicetype='none'
cryne if(idproc.eq.0)write(6,*)'AUG8, set default slicetype to none'

c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
      call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
      ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
      read(lf,2000,end=2001,err=2001) line
 2000 format(a)
      if(idproc.eq.0)                                                     &
     &write(6,*)'reading from file associated with unit 11'
      goto 4320
 2001 continue
      open(lf,file='mli.in',status='old',err=357)
      read(lf,2000,end=357,err=357) line
      if(idproc.eq.0)                                                     &
     &write(6,*)'reading from file mli.in'
      goto 4320
  357 continue
      write(6,*)'master input file does not exist or is empty'
      write(6,*)'type filename or <cr> to halt'
      read(5,2002)fname
 2002 format(a16)
      if(fname.eq.' ')call myexit
      open(lf,file=fname,status='old',err=3000)
      read(lf,2000,end=3000,err=3000) line
      goto 4320
 3000 continue
      write(6,*)'file still does not exist or is empty. Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
      rewind lf
      msegm = -1
      curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      if (leof) goto 1000
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin; line=')
      write(jof,*)TRIM(line)
      write(jodf,*)TRIM(line)
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin): ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne--- 1 August, 2004:
      omegascl=1.d0/ts
      freqscl=omegascl/( 4.d0*asin(1.d0) )
      pmass=brho/(gamma*beta/clite)
      p0sc=gamma*beta*pmass/clite
      lflagmagu=.true.
      lflagdynu=.false.
cryne---
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
c     call txtnum(line,80,4,nget,bufr)
      call txtnum(line,LEN(line),4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin near #beam'/        &
     &'Note: if using MAD-style input for beam info, it should be'/       &
     &'preceded by #menu, not #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(iverbose.eq.2)write(6,*)'#menu;',line(1:70)
      if(iverbose.eq.2)write(6,*)'#menu;',itot,strarr(1),strarr(2),msegm
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
      if(iverbose.eq.2)write(6,*)'got past 301'
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin'
c       call getsymb60(line,80,symb,istrtsym,nsymb)
        call getsymb60(line,LEN(line),symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        lftemp2=lf
        lf=42
        write(6,*)'potential bug in dumpin! Need to replace hardwired'
        write(6,*)'file unit # (42) with code-selected #. fix later'
        call dumpin2(symb(3))
        close(lf)
        lf=lftemp2
        write(6,*)'returned from dumpin2'
        goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      if(iverbose.eq.2)write(6,*)'checking for doubly defined names'
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      if(iverbose.eq.2)write(6,*)'starting do 325'
      do 325 m = 1,9
      do 325 n = 1,nrpmax
cryne Dec 1, 2002      if(string(1:8).eq.ltc(m,n)) then
c eventually the max string length should not be hardwired like this
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        if(iverbose.eq.2)then
        write(6,*)'(3344) ',na,nt1(na),nt2(na),string(1:16)
        endif
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      if(idproc.eq.0)then
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      endif
      na=na-1
      goto 300
 3344 continue
c     write(6,*)'menu element;string,m,n=',string(1:16),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
cryne 12/31/2005 imaxold is now incremented for thick elements if autoslicing.
cryne But note that mpprpoi is incremented by imax, not imaxold.
cryne This allows use of MaryLie names, but to have additional
cryne parameters if read in using the SIF (MAD) format.
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
cryne 5/4/2006
      if(line(1:1).eq.'>')then
        lmsgpoi=lmsgpoi+1
        if(lmsgpoi.gt.lmaxmsg)then
          if(idproc.eq.0)                                               &
     &    write(6,*)'error: too many output messages (>) ini #labor'
          call myexit
        endif
c       latt(noble)='>'
        latt(noble)=line(1:3)
        num(noble)=lmsgpoi
        lattmsg(lmsgpoi)=line
        goto 700
      endif
c
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      return
      end
c
************************************************************************
      subroutine dumpin2(uname)
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*60 uname
      character*16 strarr(40),string
      character*80 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
cryne 9/26/02
c this could go somewhere else too, but here is OK:
c     isymbdef=-1
c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
c     call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
c     ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
c     read(lf,2000,end=2001,err=2001) line
 2000 format(a)
c     goto 4320
c2001 continue
c     write(6,*)'master input file does not exist or is empty'
c     write(6,*)'type filename or <cr> to halt'
c     read(5,2002)fname
c2002 format(a16)
c     if(fname.eq.' ')call myexit
      open(lf,file=uname,status='old',err=3000)
      write(6,*)'successfully opened file ',uname
      goto 4320
 3000 continue
      write(6,*)'(dumpin2) file does not exist. File name=.'
      write(6,*)uname
      write(6,*)'Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
c     rewind lf
cryne 11/04/02      msegm = -1
      msegm=3
c     curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      write(6,*)'dumpin2: done reading line='
      write(6,*)line
      write(6,*)'msegm,itot,leof=',msegm,itot,leof
      if(strarr(1).eq.'call')then
        write(6,*)'found call statement after 10 in dumpin2'
        goto 301
      endif
      if (leof) goto 1000
      if(msegm.eq.3)goto 301
      if(msegm.eq.4 .or. msegm.eq.5 .or. msegm.eq.6)goto 410
      if(msegm.eq.9)goto 901
c check for the statement '#labor'
      if(msegm.eq.7)goto 700
      write(6,*)'error: read the first line of input file but'
      write(6,*)'cannot determine where to go in routine dumpin2'
      call myexit
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin2; line=')
      write(jof,*)line
      write(jodf,*)line
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin2): ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
      call txtnum(line,80,4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin2)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin2 near #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!     write(6,*)line(1:80)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin2'
        call getsymb60(line,80,symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        lftemp3=lf
        lf=43
        call dumpin3(symb(3))
        close(lf)
        lf=lftemp3
        goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin2 in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      do 325 m = 1,9
      do 325 n = 1,nrpmax
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin2 error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
 3344 continue
!     write(6,*)'found a menu element;string,m,n=',string(1:8),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin2) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
c     write(6,*)'USING STANDARD INPUT FORMAT. LINE='
c     write(6,*)line
c     write(6,*)'strarr(1),strarr(2)='
c     write(6,*)strarr(1)
c     write(6,*)strarr(2)
ccc   write(6,*)'will read params using SIF; imax,icmax=',imax,icmax
ccc   write(6,*)'current values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
c     write(6,*)'TRIMMED LINE='
c     write(6,*)line
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin2:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin2 error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin2:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin2:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      close(44)
      write(6,*)'returning from routine dumpin2'
      return
      end
c
c***********************************************************************
************************************************************************
      subroutine dumpin3(uname)
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*60 uname
      character*16 strarr(40),string
      character*80 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
cryne 9/26/02
c this could go somewhere else too, but here is OK:
c     isymbdef=-1
c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
c     call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
c     ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
c     read(lf,2000,end=2001,err=2001) line
 2000 format(a)
c     goto 4320
c2001 continue
c     write(6,*)'master input file does not exist or is empty'
c     write(6,*)'type filename or <cr> to halt'
c     read(5,2002)fname
c2002 format(a16)
c     if(fname.eq.' ')call myexit
      open(lf,file=uname,status='old',err=3000)
      write(6,*)'successfully opened file ',uname
      goto 4320
 3000 continue
      write(6,*)'(dumpin3) file does not exist. File name=.'
      write(6,*)uname
      write(6,*)'Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
c     rewind lf
cryne 11/04/02      msegm = -1
      msegm=3
c     curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      write(6,*)'dumpin3: done reading line='
      write(6,*)line
      write(6,*)'msegm,itot,leof=',msegm,itot,leof
      if(strarr(1).eq.'call')then
        write(6,*)'found call statement after 10 in dumpin3'
        goto 301
      endif
      if (leof) goto 1000
      if(msegm.eq.3)goto 301
      if(msegm.eq.4 .or. msegm.eq.5 .or. msegm.eq.6)goto 410
      if(msegm.eq.9)goto 901
c check for the statement '#labor'
      if(msegm.eq.7)goto 700
      write(6,*)'error: read the first line of input file but'
      write(6,*)'cannot determine where to go in routine dumpin3'
      call myexit
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin3; line=')
      write(jof,*)line
      write(jodf,*)line
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin3: ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
      call txtnum(line,80,4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin3)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin3 near #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!     write(6,*)line(1:80)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin3'
        call getsymb60(line,80,symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        lftemp4=lf
        lf=44
        call dumpin4(symb(3))
        close(lf)
        lf=lftemp4
        goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      do 325 m = 1,9
      do 325 n = 1,nrpmax
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
 3344 continue
!     write(6,*)'found a menu element;string,m,n=',string(1:8),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
c     write(6,*)'USING STANDARD INPUT FORMAT. LINE='
c     write(6,*)line
c     write(6,*)'strarr(1),strarr(2)='
c     write(6,*)strarr(1)
c     write(6,*)strarr(2)
ccc   write(6,*)'will read params using SIF; imax,icmax=',imax,icmax
ccc   write(6,*)'current values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
c     write(6,*)'TRIMMED LINE='
c     write(6,*)line
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      close(44)
      write(6,*)'returning from routine dumpin3'
      return
      end
c
c***********************************************************************
************************************************************************
      subroutine dumpin4(uname)
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*60 uname
      character*16 strarr(40),string
      character*80 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
cryne 9/26/02
c this could go somewhere else too, but here is OK:
c     isymbdef=-1
c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
c     call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
c     ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
c     read(lf,2000,end=2001,err=2001) line
 2000 format(a)
c     goto 4320
c2001 continue
c     write(6,*)'master input file does not exist or is empty'
c     write(6,*)'type filename or <cr> to halt'
c     read(5,2002)fname
c2002 format(a16)
c     if(fname.eq.' ')call myexit
      open(lf,file=uname,status='old',err=3000)
      write(6,*)'successfully opened file ',uname
      goto 4320
 3000 continue
      write(6,*)'(dumpin4) file does not exist. File name=.'
      write(6,*)uname
      write(6,*)'Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
c     rewind lf
cryne 11/04/02      msegm = -1
      msegm=3
c     curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      write(6,*)'dumpin4: done reading line='
      write(6,*)line
      write(6,*)'msegm,itot,leof=',msegm,itot,leof
      if(strarr(1).eq.'call')then
        write(6,*)'found call statement after 10 in dumpin4'
        goto 301
      endif
      if (leof) goto 1000
      if(msegm.eq.3)goto 301
      if(msegm.eq.4 .or. msegm.eq.5 .or. msegm.eq.6)goto 410
      if(msegm.eq.9)goto 901
c check for the statement '#labor'
      if(msegm.eq.7)goto 700
      write(6,*)'error: read the first line of input file but'
      write(6,*)'cannot determine where to go in routine dumpin4'
      call myexit
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin4; line=')
      write(jof,*)line
      write(jodf,*)line
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin4): ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
      call txtnum(line,80,4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin4)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin4 near #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!     write(6,*)line(1:80)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin4'
        call getsymb60(line,80,symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        lftemp5=lf
        lf=45
        call dumpin5(symb(3))
        close(lf)
        lf=lftemp5
        goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin4:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin4 in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      do 325 m = 1,9
      do 325 n = 1,nrpmax
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin4 error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
 3344 continue
!     write(6,*)'found a menu element;string,m,n=',string(1:8),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin4) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
c     write(6,*)'USING STANDARD INPUT FORMAT. LINE='
c     write(6,*)line
c     write(6,*)'strarr(1),strarr(2)='
c     write(6,*)strarr(1)
c     write(6,*)strarr(2)
ccc   write(6,*)'will read params using SIF; imax,icmax=',imax,icmax
ccc   write(6,*)'current values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
c     write(6,*)'TRIMMED LINE='
c     write(6,*)line
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin4:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin4 error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin4:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin4:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      close(44)
      write(6,*)'returning from routine dumpin4'
      return
      end
c
c***********************************************************************
************************************************************************
      subroutine dumpin5(uname)
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*60 uname
      character*16 strarr(40),string
      character*80 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
cryne 9/26/02
c this could go somewhere else too, but here is OK:
c     isymbdef=-1
c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
c     call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
c     ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
c     read(lf,2000,end=2001,err=2001) line
 2000 format(a)
c     goto 4320
c2001 continue
c     write(6,*)'master input file does not exist or is empty'
c     write(6,*)'type filename or <cr> to halt'
c     read(5,2002)fname
c2002 format(a16)
c     if(fname.eq.' ')call myexit
      open(lf,file=uname,status='old',err=3000)
      write(6,*)'successfully opened file ',uname
      goto 4320
 3000 continue
      write(6,*)'(dumpin5) file does not exist. File name=.'
      write(6,*)uname
      write(6,*)'Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
c     rewind lf
cryne 11/04/02      msegm = -1
      msegm=3
c     curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      write(6,*)'dumpin5: done reading line='
      write(6,*)line
      write(6,*)'msegm,itot,leof=',msegm,itot,leof
      if(strarr(1).eq.'call')then
        write(6,*)'found call statement after 10 in dumpin5'
        goto 301
      endif
      if (leof) goto 1000
      if(msegm.eq.3)goto 301
      if(msegm.eq.4 .or. msegm.eq.5 .or. msegm.eq.6)goto 410
      if(msegm.eq.9)goto 901
c check for the statement '#labor'
      if(msegm.eq.7)goto 700
      write(6,*)'error: read the first line of input file but'
      write(6,*)'cannot determine where to go in routine dumpin5'
      call myexit
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin5; line=')
      write(jof,*)line
      write(jodf,*)line
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin5): ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
      call txtnum(line,80,4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin5)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin5 near #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!     write(6,*)line(1:80)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin5'
        call getsymb60(line,80,symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        lftemp6=lf
        lf=46
        call dumpin6(symb(3))
        close(lf)
        lf=lftemp6
        goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin5:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin5 in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      do 325 m = 1,9
      do 325 n = 1,nrpmax
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin5 error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
 3344 continue
!     write(6,*)'found a menu element;string,m,n=',string(1:8),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin5) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
c     write(6,*)'USING STANDARD INPUT FORMAT. LINE='
c     write(6,*)line
c     write(6,*)'strarr(1),strarr(2)='
c     write(6,*)strarr(1)
c     write(6,*)strarr(2)
ccc   write(6,*)'will read params using SIF; imax,icmax=',imax,icmax
ccc   write(6,*)'current values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
c     write(6,*)'TRIMMED LINE='
c     write(6,*)line
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin5:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin5 error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin5:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin5:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      close(44)
      write(6,*)'returning from routine dumpin5'
      return
      end
c
c***********************************************************************
************************************************************************
      subroutine dumpin6(uname)
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c               October 21, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'incmif.inc'
      include 'codes.inc'
      include 'files.inc'
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
      common/symbdef/isymbdef
      dimension bufr(4)
c-----------------------------------------------------------------------
c local variables:
      character*60 uname
      character*16 strarr(40),string
      character*80 line
      character*16 fname
      integer narr(40)
      character*60 symb(40)
      integer istrtsym(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
cryne 9/26/02
c this could go somewhere else too, but here is OK:
c     isymbdef=-1
c
cryne 7/20/2002 initialize character array so that, even if input is
c in the original MaryLie format, the mppc array will be filled with ' '
c so that files are opened properly. (code checks for name=' ')
cryne 7/23/2002 done in new_acceldata      cmenu(1:mnumax)=' '
cryne 7/7/2002
cryne initialize #const
cryne this should go somewhere else (near start of afro.f), but I am
cryne trying not to change MaryLie too much.
c     call initcons
cryne 7/7/2002 additional mods to deal w/ huge number of comments
cryne eventually it would be a good idea to let the user specify whether
cryne or not comments should be printed in the pmif command
cryne
c start
c     ignorcom=0
cryne----- 15 Sept 2000 modified to check for input file:
ctm   open(lf,file='fort.11',status='old',err=357)
c     read(lf,2000,end=2001,err=2001) line
 2000 format(a)
c     goto 4320
c2001 continue
c     write(6,*)'master input file does not exist or is empty'
c     write(6,*)'type filename or <cr> to halt'
c     read(5,2002)fname
c2002 format(a16)
c     if(fname.eq.' ')call myexit
      open(lf,file=uname,status='old',err=3000)
      write(6,*)'successfully opened file ',uname
      goto 4320
 3000 continue
      write(6,*)'(dumpin6) file does not exist. File name=.'
      write(6,*)uname
      write(6,*)'Halting.'
      call myexit
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
c     rewind lf
cryne 11/04/02      msegm = -1
      msegm=3
c     curr0=-99999.
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,*) itot,' = itot after first REAREC'
      write(6,*)'dumpin6: done reading line='
      write(6,*)line
      write(6,*)'msegm,itot,leof=',msegm,itot,leof
      if(strarr(1).eq.'call')then
        write(6,*)'found call statement after 10 in dumpin6'
        goto 301
      endif
      if (leof) goto 1000
      if(msegm.eq.3)goto 301
      if(msegm.eq.4 .or. msegm.eq.5 .or. msegm.eq.6)goto 410
      if(msegm.eq.9)goto 901
c check for the statement '#labor'
      if(msegm.eq.7)goto 700
      write(6,*)'error: read the first line of input file but'
      write(6,*)'cannot determine where to go in routine dumpin6'
      call myexit
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800,900),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin6; line=')
      write(jof,*)line
      write(jodf,*)line
      call myexit
c--------------------
c  #comment
c
  100 continue
!     write(6,*)'here I am at #comment'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
      if(ignorcom.eq.1)goto 100
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'warning (dumpin6): ',                                 &
     & 'entries in #comment beyond line ',i6,'  will be ignored')
        ignorcom=1
        npcom=maxcmt
        goto 100
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
!     write(6,*)'here I am at #beam'
      read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
cryne 08/14/2001 the beam component may contain other info:
c first set defaults
      magunits=1
      iverbose=0
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!?    if(msegm.ne.2)goto 1
      if (npound) goto 1
      if(msegm.eq.3)goto 301
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if(msegm.ne.2) goto 1
!?
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c need to delete the following code, which is no longer useful or works.
      if(msegm.ne.12345)then
       write(6,*)'reached bad section of code in DUMPIN. Stopping.'
       stop
      endif
cryne there is other info in the #beam component:
      call txtnum(line,80,4,nget,bufr)
      if(nget.ne.4)then
        write(6,*)'(dumpin6)error: could not read 4 more #beam numbers'
        stop
      endif
c     write(6,*)'bufr(1),bufr(2),bufr(3)=',bufr(1),bufr(2),bufr(3)
cryne current:
      curr0=bufr(1)
      if(curr0.ne.-99999.)write(6,*)'beam current = ',curr0
cryne units:
      if(nint(bufr(3)).ne.0)then
        magunits=0
        sclfreq=4.*asin(1.d0)*bufr(3)
        write(6,*)'input frequency=',bufr(3),' Hz'
        write(6,*)'scale frequency=',sclfreq,' rad/sec'
        write(6,*)'dynamic units (not magnetostatic units) will be used'
        slnew=c/sclfreq
        write(6,*)'scale length specified in the input file is',sl,'m'
        write(6,*)'resetting the scale length to c/scalefreq=',slnew,'m'
        if(abs(sl-slnew)/sl .gt. 1.e-3)then
          write(6,*)'WARNING: YOU SPECIFIED A SCALE LENGTH THAT IS'
          write(6,*)'SIGNIFICANTLY DIFFERENT FROM THE NEW VALUE.'
          write(6,*)'MAKE SURE THAT YOUR INITIAL CONDITIONS ARE'
          write(6,*)'SPECIFIED IN DYNAMIC UNITS PRIOR TO TRACING RAYS'
        endif
        sl=slnew
      endif
cryne verbose to show progress:
      iverbose=bufr(4)
cryne autoslicing (thick elements need an extra parameter):
      iautosl=nint(bufr(2))
      if(iautosl.ne.0)then
        write(6,*)'autoslicing of thick elements will be enabled.'
      endif
      if(iautosl.gt.0)then
        write(6,*)'fixed # of slices/element will be = ',iautosl
      endif
      if(iautosl.lt.0)then
        write(6,*)'variable # of slices/element will be used'
        if(na.ne.0)then
         write(6,*)'error: input file specifies variable # of'
         write(6,*)'slices/element, but some elements have already'
         write(6,*)'been read in. stopping.'
         stop
        endif
        nrp(1,1)=nrp(1,1)+1
        nrp(1,2)=nrp(1,2)+1
        nrp(1,3)=nrp(1,3)+1
        nrp(1,4)=nrp(1,4)+1
        nrp(1,6)=nrp(1,6)+1
        nrp(1,8)=nrp(1,8)+1
        nrp(1,9)=nrp(1,9)+1
        nrp(1,10)=nrp(1,10)+1
        nrp(1,11)=nrp(1,11)+1
        nrp(1,12)=nrp(1,12)+1
        nrp(1,18)=nrp(1,18)+1
        nrp(1,20)=nrp(1,20)+1
        nrp(1,24)=nrp(1,24)+1
        nrp(1,30)=nrp(1,30)+1
      endif
      goto 10
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c error exit:
  290 continue
      write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin6 near #beam')
      call myexit
c--------------------
c  #menu
c
  300 continue
!     write(6,*)'here I am at #menu'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
!     write(6,*)line(1:80)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.4)goto 410
      if (msegm .ne. 3) goto 1
  301 continue
c
c 11/4/02 new code to deal with 'call' statements:
      if(trim(strarr(1)).eq.'call')then
        write(6,*)'found a call statement in dumpin6'
        call getsymb60(line,80,symb,istrtsym,nsymb)
        write(6,*)'symb(1), symb(2), symb(3)='
        write(6,*)symb(1)
        write(6,*)symb(2)
        write(6,*)symb(3)
        write(6,*)'this exceeds the maximum amount of nested calls'
        write(6,*)'in the dumpin subroutines; halting.'
        call myexit
c for more nesting, uncomment the following and add routine dumpin7
cccc    lftemp7=lf
cccc    lf=47
cccc    call dumpin7(symb(3))
cccc    close(lf)
cccc    lf=lftemp7
cccc    goto 300
      endif
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin6:',                                  &
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin6 in #menu: name ',             &
     & a16,' is doubly defined'/                                         &
!    &         ' the second definition will be ignored')
!       na = na-1
!       goto 300
     &         ' the first definition will be ignored')
        lmnlbl(indx)(1:1)='*'
      endif
cryne
      if((trim(strarr(1)).ne.'beam').and.                               &
     &(trim(strarr(1)).ne.'units'))then
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a16,' has no type code!')
        call myexit
      endif
      endif
cryne July 4 2002      if ( itot .gt. 2 ) then
cryne July 4 2002        write(jof ,1397) strarr(1)
cryne July 4 2002        write(jodf ,1397) strarr(1)
c1397   format(' error in #menu: ',a16,' has more than one type code!')
cryne July 4 2002        call myexit
cryne July 4 2002      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
      initparms=1  !cryne 12/17/2004 change to 0 in the case of element reuse
c string is name of element/command type, look up element indices
cryne      string=strarr(2)
      if((trim(strarr(1)).eq.'beam').or.                                &
     &(trim(strarr(1)).eq.'units'))then
        string=strarr(1)
      else
        string=strarr(2)
      endif
      do 325 m = 1,9
      do 325 n = 1,nrpmax
      if(string(1:16).eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
        goto 3344
      endif
  325 continue
c=============
cryne Sept 17, 2003
c before declaring this an unknown element/command name,
c check to see if it is the name of an existing menu item
c since this is a functionality that MAD allows
      call lookup(strarr(2),itype,indx)
      if(itype.ne.1)goto 3226
c this *is* the name of an existing menu item.
      if(idproc.eq.0)then
        write(6,*)'Element ',strarr(1),' is derived from ',strarr(2)
      endif
      m = nt1(indx)
      n = nt2(indx)
      nt1(na) = m
      nt2(na) = n
      imax=nrp(m,n)
      if(imax.ne.0)then
        do i=1,imax
         pmenu(i+mpprpoi)=pmenu(i+mpp(indx))
        enddo
      endif
      icmax=ncp(m,n)
      if(icmax.ne.0)then
        do i=1,icmax
         cmenu(i+mppcpoi)=cmenu(i+mppc(indx))
        enddo
      endif
      initparms=0 !cryne 12/17/2004 skip parameter initialization in stdinf.f
      goto 3344
 3226 continue
c=============
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin6 error at ',a16,': type code ',a16,             &
     &' not found.'/                                                    &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
 3344 continue
!     write(6,*)'found a menu element;string,m,n=',string(1:8),m,n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        icmax=ncp(m,n)
        imaxold=nrpold(m,n)
c       write(6,*)'ready to read parameters; imax,icmax=',imax,icmax
cryne 9/26/02        if(imax.eq.0.and.icmax.eq.0)goto 300
        if(imax.eq.0.and.icmax.eq.0.and.itot.eq.2.and.
     &  index(line,':').eq.0)goto 300
cryne July 4, 2002
c if using the Standard Input Format, the remaining parameters
c are stored in "line"; otherwise use the MaryLie input format:
c
c mpp(na) points to where the real    parameters will be stored
c mppc(na) points to where the char*16 parameters will be stored
      mpp(na) = mpprpoi
      mppc(na) = mppcpoi
cryne 7/9/2002      if(itot.eq.2)then
c original MaryLie input format (icmax not relevant here):
      if( (itot.eq.2) .and. (index(line,':').eq.0) )then
c       write(6,*)'reading ',imax,' params for ',strarr(1),' - ',string
        if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imaxold)
c       if(imax.ne.0)read(lf,*,err=390) (pmenu(i+mpp(na)),i=1,imax)
c       if(icmax.ne.0)read(lf,*,err=390)(cmenu(i+mppc(na)),ic=1,imax)
        mpprpoi = mpprpoi + imax
        mppcpoi = mppcpoi + icmax
        goto 300
c error in parameter input
  390   write(jof ,397)lmnlbl(na)
        write(jodf,397)lmnlbl(na)
  397   format(1h ,'(dumpin6) parameter input error at element ',a16)
        call myexit
      endif
cryne July 4, 2002
c if the code gets here, the Standard Input Format is being used
c to specify the elements in the menu:
c first delete the element name and type from the character string:
c     write(6,*)'USING STANDARD INPUT FORMAT. LINE='
c     write(6,*)line
c     write(6,*)'strarr(1),strarr(2)='
c     write(6,*)strarr(1)
c     write(6,*)strarr(2)
ccc   write(6,*)'will read params using SIF; imax,icmax=',imax,icmax
ccc   write(6,*)'current values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
cryne this could be more easily done with the f90 intrinsic len_trim()
c     kkk1=len_trim(strarr(1))-1
      string=strarr(1)
      kkk1=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk1=kkk1-1
      enddo
      do i=1,80
        if(line(i:i+kkk1).eq.strarr(1))then
          line(i:i+kkk1)=' '
          exit
        endif
      enddo
cryne kkk2=len_trim(strarr(2))-1
      string=strarr(2)
      kkk2=len(string)-1
      do i=1,len(string)
      if(string(i:i).eq.' ')kkk2=kkk2-1
      enddo
c
      if((trim(strarr(1)).ne.'beam').and.                                &
     &(trim(strarr(1)).ne.'units'))then
        do i=1,80
          if(line(i:i+kkk2).eq.strarr(2))then
            line(i:i+kkk2)=' '
            exit
          endif
        enddo
      endif
c get rid of other unneccesary characters:
      do i=1,80
c       if(line(i:i).eq.',')line(i:i)=' '
        if(line(i:i).eq.':')line(i:i)=' '
      enddo
c     write(6,*)'TRIMMED LINE='
c     write(6,*)line
      call stdinf(line,na,m,n,initparms,strarr(1))
      mpprpoi = mpprpoi + imax
      mppcpoi = mppcpoi + icmax
c       write(6,*)'read parameters using SIF; imax,icmax=',imax,icmax
c       write(6,*)'new values of mpprpoi,mppcpoi=',mpprpoi,mppcpoi
      goto 300
c
c--------------------
c  #lines,#lumps,#loops
c
  400 continue
!     write(6,*)'here I am at #lines,lumps,loops'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if(msegm.eq.9)goto 901
      if(msegm.eq.3)goto 301
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin6:',                                  &
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin6 error in ',                                   &
     &  a8,': name ',a16,' is doubly defined'/                          &
!    &         ' the second definition will be ignored')
!       na = na-1
!       if(msegm .eq.4) then
!         goto 400
!       else if(msegm.eq.5) then
!         goto 500
!       else if(msegm.eq.6) then
!         goto 600
!       endif
     &         ' the first definition will be ignored')
        ilbl(indx)(1:1)='*'
      endif
c new item
      ilbl(nb)=strarr(1)
      imin=0
cryne July 2002 Standard Input Format option: components might be on
cryne this record; check now
cryne Note: this assume that the 2nd string after the name is 'line',
cryne and so the 2nd string is ignored
      if( (itot.eq.2) .and. (.not.(lcont)) )then
       write(6,*)'error parsing line:'
       write(6,*)'this should be a name alone, or'
       write(6,*)'a name followed by LINE= followed by the components'
       write(6,*)'input line ='
       write(6,*)line
       stop
      endif
      if(itot.gt.2)then
c       write(6,*)'ITOT=',itot
c       do i=1,itot
c         write(6,*)i,narr(i),strarr(i)
c       enddo
        do 424 i=3,itot
          icon(imin+i-2,nb)=strarr(i)
          irep(imin+i-2,nb)=narr(i)
  424   continue
        imin=imin+itot-2
        if(lcont)then
          goto 420
        else
          goto 426
        endif
      endif
cryne remaining code is from original version:
c read components of item
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin6:',                                  &
     & ' too many entries (> itmln = ',i6,') in ',a16)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
  426 continue
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 continue
!     write(6,*)'here I am at #labor'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin6:',                                  &
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
c--------------------
c  #call (formerly #include)
  800 continue
      write(6,*)' Start #call    after user name: ',strarr(1)
c
ctm 9/01 modified to open and save up to 32 long #include file names
c
      ninc = 0
  810 read(lf,2000,end=1000) line
c
c    check for next segment
c
      if(index(line,'#').ne.0) then
        itot = -1
        na2 = na
        nb2 = nb
        noble2 = noble
        write(6,813) na2,nb2,noble2
  813   format(' After #include:',i5,' menu',i5,' items',i5,' tasks')
        write(6,817) line
  817   format(' End #include with: ',a)
         go to 10
      endif
c
ctm       save file name and pointers before opening 1st include file
c
      ninc = ninc + 1
      if(ninc.eq.1) then
         na1 = na
         nb1 = nb
         noble1 = noble
         write(6,877) na1,nb1,noble1
 877    format(' Before #include:',i5,' menu',i5,' items',i5,' tasks')
      endif
c
ctm    skip leading blanks
c
      lc = 1
 880  if(line(lc:lc).eq.' ') then
        lc = lc + 1
        go to 880
      endif
      incfil(ninc) = line(lc:)
      write(6,888) ninc,incfil(ninc)
 888  format(' include',i3,' : ',a)
      call mlfinc(incfil(ninc))
      go to 810
ctm      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      if (leof) goto 1000
ctm      goto 1
c--------------------
c  #const
  900 continue
!     write(6,*)'here I am at #const'
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (npound) goto 1
      if (leof) goto 1000
      if((msegm.eq.9).and.(itot.ne.1))then
        write(6,*)'input error (#const)'
        write(6,*)'trying to read a definition, but did not find a'
        write(6,*)'character string at the beginning of this record:'
        write(6,*)line(1:80)
        stop
      endif
      if(msegm.eq.3)goto 301
      if(msegm.eq.4)goto 410
      if (msegm .ne. 9) goto 1
  901 continue
      call getconst(line,strval,nreturn)
      if(nreturn.ne.1)then
        write(6,*)'trouble parsing the following line:'
        write(6,*)line
        write(6,*)'continuing...'
        goto 900
      endif
      nconst=nconst+1
      if(nconst.gt.nconmax)then
        write(6,*)'too many constants defined in the input file'
        stop
      endif
      constr(nconst)=strarr(1)
      conval(nconst)=strval
c     write(6,*)'nconst,constr(nconst),conval(nconst)='
c     write(6,*)nconst,constr(nconst),conval(nconst)
      goto 900
c normal return at end of file
 1000 continue
      close(44)
      write(6,*)'returning from routine dumpin6'
      return
      end
c
