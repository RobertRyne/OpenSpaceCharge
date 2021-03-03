cryne The bulk of this file (sif.f) is devoted to reading in
cryne parameters in the Standard Input Format
c
      subroutine getconst(line,strval,nreturn)
cryne read & parse a symbolic constant defined in the input file
cryne Note: this assumes only one definition per input record,
cryne i.e. can't have multiple definitions separated by ; if using MADX format.
cryne If using MADX format, the record terminates with the first ;
      use parallel
      include 'impli.inc'
      character line*(*)
cryne 5/4/2006      character linetmp*256
      character linetmp*800
      character cdummy*1
      integer   lentmp
      logical keypres,numpres
      dimension bufr(1)
c
cryne 5/4/2006      linetmp=line
cryne 5/4/2006      lentmp = LEN( linetmp )
c I don't remember exactly why I added readin800.
c I think it was to process symbolic constants defined on continued lines.
c
      call readin800(line,linetmp,mmax)
      lentmp=mmax
      nreturn=0
c

cryne---5/4/2006 I don't know exactly why I commented this out
cryne It must be that this (and other) data processing is handled in readin800
c     n1=index(line,'!')
c     if(n1.ne.0)then
c       do n=n1,lentmp
c         linetmp(n:n)=' '
c       enddo
c     endif
c for now assume that there is only one definition per record:
c     n1=index(line,';')
c     if(n1.ne.0)then
c       do n=n1,lentmp
c         linetmp(n:n)=' '
c       enddo
c     endif
cryne---
      call getparm(linetmp,lentmp,'=',bufr,keypres,numpres,0,cdummy)
      if(numpres)then
        nreturn=1
        strval=bufr(1)
cc      write(6,*)'returning from getparm w/ strval=',strval
      endif
      return
      end
c
      subroutine initcons
      use acceldata
      include 'impli.inc'
c initialize the first few values of #const
      pi=2.d0*asin(1.d0)
      twopi=4.d0*asin(1.d0)
      clite=299792458.d0
      constr(1)='pi'
      conval(1)=pi
      constr(2)='twopi'
      conval(2)=twopi
      constr(3)='clite'
      conval(3)=clite
      nconst=3
c     write(6,*)'nconst=',nconst
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stdinf(line80,na,kt1,kt2,initparms,lmntname)
c read a record in the Standard Input Format
c input:
c   line80=input line
c   na=na'th entry in the menu
c   kt1, kt2 = indices describing entry category and type
cryne 5/4/2006: removed a lot of code from stdinf and put it in readin800
cryne 5/4/2006: na appears unused. formerly for debugging?
cryne 12/17/2004:
cryne added initparms (controls whether or not to initialize parameters)
cryne added lmntname (element label; used for diagnostic info & debugging)
      include 'impli.inc'
      character line80*(*)
      character line*800
      character*16 lmntname
      call readin800(line80,line,mmax)
      call lmntparm(line,mmax,kt1,kt2,initparms,lmntname)
      return
      end
c
      subroutine readin800(line80,line,mmax)
cryne 5/4/2006 store multiline input into a buffer of length 800
cryne and do some minimal processing of the line
c   line80=first 80 character input line
c   line=full 800 characture buffer
c   mmax=length of resulting line (i.e. nonblank character length)
      include 'impli.inc'
      character line80*(*)
      character line*800
      logical leof
      dimension bufr(1)
      integer ilen
      common/showme/iverbose
c
cryne 5/4/2006
      logical MAD8,MADX
      common/mad8orx/MAD8,MADX
c
c     write(6,*)'inside readin800'
c
c
      leof = .FALSE.

!     if(iverbose.ge.6 .AND. idproc.eq.0)then
!       write(6,*)'entering stdinf; kt1,kt2=',kt1,kt2
!       write(6,*)'line80(1:78) = ',line80(1:78)
!     endif

!XXX -- I rewrote this <dbs> Aug04
      ! read input lines until all the data for this entry is read;
      ! and pack it into 'line'
      i1 = 0
      i2 = 0
      do !until eof or entry is terminated
        ! first iteration uses line80 from caller; rest get input from file
        if( i1 .NE. 0 )then
!         write(6,*)'calling readin'
          call readin(line80,leof)
!         write(6,*)'line80(1:78) = ',line80(1:78)
          if( leof )then
            if( iverbose.ge.6 .AND. idproc.eq.0 )then
              write(6,*) 'info: STDINF(): EOF happend with line80='
     &                  ,line80
            endif
            exit
          endif
        endif
        line80 = ADJUSTL( line80 ) !remove leading blanks
        ilen = LEN_TRIM( line80 )  !dont look at trailing blanks
        ! skip additional blank lines;
        ! line80 initially blank is ok but doesnt need processing
        if( ilen .EQ. 0 )then
          if(i1 .GT. 0 )cycle
        else
          ! debug
          if( iverbose.ge.6 .AND. idproc.eq.0 )then
            write(6,*)'info: STDINF(): line80='
            write(6,*) line80(:ilen)
          endif
          i1 = i2 + 1  !start inserting after end of previous line80
          ! find comments and ignore them
          j1 = INDEX( line80(:ilen) ,'!' )
          if( j1 .GT. 0 )then
            ! skip if nothing on the line except the comment
            if( j1 .EQ. 1 ) cycle
            ! strip trailing blanks before the '!'
            ilen = LEN_TRIM(line80(:j1-1))
            ! skip if non-comment part of the line is blank
            if( ilen .EQ. 0 ) cycle
          endif
          ! check for not enough buffer space
          i2 = i1 + ilen - 1
          if( i2 .GT. LEN( line ) )then
            write(6,*) 'error: STDINF(): input too long for line buffer'
            write(6,*) 'info: input is [',line80(:ilen),']'
            stop
          endif
          ! convert to lowercase
          ![NOTE: readin() does this, so this is probably redundant. <dbs>]
          call LOW( line80(:ilen) )
          if( iverbose.ge.6 .AND. idproc.eq.0 )then
            write(6,*) 'info: STDINF(): after low(), input line is '
            write(6,*) line80(:ilen)
          endif
cryne 5/4/2006
cryne in MADX style, there could be space between ; and !, so deal with it:
          do iq=ilen,1,-1
          ilast=iq
          if(line80(ilast:ilast).ne.' ')exit
          enddo
cryne
          ! save this line
          line(i1:i2) = line80(:ilen)
          ! check if there are more lines to read
          if( (MAD8)  .and.  line80(ilen:ilen) .EQ. '&' )then
            ! add a space to the buffer to make sure words are separated
            i2 = i2 + 1
            line(i2:i2) = ' '
            ! continue to next line
            cycle
          elseif( (MADX) .and. line80(ilast:ilast) .NE. ';' )then
          write(6,*)'continuing...'
!         write(6,*)line80(1:ilen)
            i2 = i2 + 1
            line(i2:i2) = ' '
            cycle
          endif
        endif
        ! line is blank or, for ML and MAD8, this is the end,
        ! but for MAD9 we should keep going until a ";" is found:
!cryne  if( MAD9 ) cycle
        exit
      enddo

      ! finished with packing all the input lines into the buffer
      mmax = i2
      if( iverbose.ge.6 .AND. idproc.eq.0 )then
        write(6,*) 'info: STDINF(): completed input buffer ='
        write(6,*) line(:mmax)
      endif
cryne 5/4/2006
cryne skip the rest of there is no equal sign on this line
c     write(6,*)'nearly finished in readin800 having found:'
c     write(6,*)line(1:mmax)
      if(index(line,'=').eq.0)then
        write(6,*)'leaving readin800 without cleanup up near = sign'
        goto 999
      endif

      ! remove blanks around "=";
      ! move each character to the left
      ! by the number of blanks found up to that character
      move = 0
      i1 = 2                    !can skip the 1st char
      do while( i1 .LE. mmax )  !loop up to the last non-blank
        ![NOTE: use a 'while' loop because i1 gets changed inside the loop]
        if( line(i1:i1) .EQ. '=' )then
          ! count spaces to the left (dont go past the start of the string)
          do i2 = 1 ,i1-1
            if( line(i1-i2:i1-i2) .NE. ' ' ) exit
          enddo
          ! increase the move distance by the number of spaces found
          move = move + i2 - 1
          ! do the move, if necessary
          if( move .GT. 0 ) line(i1-move:i1-move) = line(i1:i1)
          ! count spaces to the right (dont go past the end)
          do i2 = 1 ,mmax-i1
            if( line(i1+i2:i1+i2) .NE. ' ' ) exit
          enddo
          ! increase the move distance by the number of spaces found
          move = move + i2 - 1
          ! skip over the spaces on the right so we dont move them
          i1 = i1 + i2 - 1
        else
          ! not '=', so just move it, if necessary
          if( move .GT. 0 ) line(i1-move:i1-move) = line(i1:i1)
        endif
        i1 = i1 + 1
      enddo
      ! clean up chars at end of the string that were moved but not overwritten
      if( move .GT. 0 ) line(mmax-move+1:mmax) = ' '

      if( iverbose.ge.6 .AND. idproc.eq.0 )then
        write(6,*) 'info: STDINF(): completed parsing, line ='
        write(6,*) line(:mmax)
      endif
cryne 5/4/2006
      mmax=mmax-move
cryne now lmntparm is called separately from the above parsing
cryne call lmntparm(line,mmax,kt1,kt2,initparms,lmntname)
cryne call lmntparm(line,mmax-move,kt1,kt2,initparms,lmntname)
cryne debugging:
  999 continue
c     write(6,*)'leaving readin800 at end having found:'
c     write(6,*)line(1:mmax)
      return
      end
c
c
      subroutine lmntparm(line,mmax,kt1,kt2,initparms,lmntname)
cryne 12/17/2004 added initparms
cryne If initparms=1 the parameters get initialized, otherwise they do not.
cryne This is needed when a menu element is definied in terms of a previous
cryne element, in which case the initial values are inherited, not set here.
cryne 12/17/2004 also added lmntname: used for diagnostic info & debugging
      use rays
      use acceldata
      use beamdata
      include 'impli.inc'
      include 'pie.inc'
      include 'codes.inc'   ! added Dec 1, 2002 to check nrp
      include 'setref.inc'  ! added by RDR, April 18, 2004 (used in BEAM)
      include 'files.inc'!RDR 7/29/04 (to backspace(lf) [wake code after BEAM])
c     character*800 line
      character line*(*)
      character*16 cbuf
      character*16 lmntname
cryne 5/4/2006      logical keypres,keypres1,keypres2,numpres
      logical keypres,keypres1,keypres2,numpres,leof
      dimension bufr(1)
      common/mlunits/sclfreq,magunits
      common/showme/iverbose
!     common/autslice/iautosl
!     common/autoslice/islicetype,isliceprecedence,islicevalue
      common/symbdef/isymbdef
      data multmsg/1/
      save multmsg
c
      if(idproc.eq.0.and.iverbose.ge.6)then
        write(6,*)'(lmntparm) mmax,kt1,kt2=',mmax,kt1,kt2
        write(6,*)'line=',line(1:mmax)
      endif
c
c files.inc contains the unit # (lf) connected to the master input file.
c this is needed so that that wake_init reads from the correct file
      iunitnum=lf
cryne=== 9/16/2004 new code to deal with "at='
cryne=== (needed to emulate MAD "sequence")
c note:  eventually, it might be better not to have a separate array
c        called atarray, but instead to store the "at" info in pmenu
      call getparm(line,mmax,'at=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)atarray(1+mpp(na))=bufr(1)
cryne===
cryne===
c
c eventually "c" should be replaced by "clite" in common/parm/
      clite=c
c Select appropriate action depending on value of kt1,kt2
c
      go to (11,12,13,14,15,16,17,18,19), kt1
c
c     1: simple elements *************************************
c
11    continue
c     write(6,*)'line(1:80)='
c     write(6,*)line(1:80)
c     write(6,*)'(lmntparm) kt1,kt2=',kt1,kt2
c          'drft    ','nbnd    ','pbnd    ','gbnd    ','prot    ',
      go to(101,       102,       103,       104,       105,
c          'gbdy    ','frng    ','cfbd    ','quad    ','sext    ',
     &      106,       107,       108,       109,       110,
c          'octm    ','octe    ','srfc    ','arot    ','twsm    ',
     &      111,       112,       113,       114,       115,
c          'thlm    ','cplm    ','cfqd    ','dism    ','sol     ',
     &      116,       117,       118,       119,       120,
c          'mark    ','jmap    ','dp      ','recm    ','spce    ',
     &      121,       122,       123,       124,       125,
c          'cfrn    ','coil    ','intg    ','rmap    ','arc     ',
     &      126,       127,       128,       129,       130,
c          'rfgap   ','confoc  ','spare1  ','spare2  ','spare3  ',
     &      1135,      1136,      133,       134,       135,
c          'spare4  ','spare5  ','spare6  ','spare7  ','spare8  ',
     &      136,       137,       138,       139,       140,
c          'marker  ','drift   ','rbend   ','sbend   ','gbend   ',
     &      141,       142,       143,       1045,      145,
c          'quadrupo','sextupol','octupole','multipol','solenoid',
     &      146,       147,       148,       149,       150,
c          'hkicker ','vkicker ','kicker  ','rfcavity','elsepara',
     &      151,       152,       153,       154,       155,
c          'hmonitor','vmonitor','monitor ','instrume','sparem1 ',
     &      156,       157,       158,       159,       160,
c          'rcollima','ecollima','yrot    ','srot    ','prot3   ',
     &      161,       162,       163,       164,       165,
c          'beambeam','matrix  ','profile1d','yprofile','tprofile',
     &      166,       167,       168,       169,       170,
c          'hkick   ','vkick   ','kick    ','sparem6 ','nlrf    '/
     &      171,       172,       173,       174,       175),kt2
c
c=====================================================================
c                                DRFT
c=====================================================================
c 'drft    ': drift
101   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,1).eq.1)then
         write(6,*)'DRFT input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         pmenu(2+mpp(na))=1.
         if(numpres)pmenu(2+mpp(na))=bufr(1)
        endif
      endif
      return
c
c
c=====================================================================
c                                NBND
c=====================================================================
c 'nbnd':
102   continue
c     if(idproc.eq.0)write(6,*)'nbnd element is: ',lmntname
c defaults:
      if(initparms.eq.1)then
      angle=0.
      pmenu(1+mpp(na))=angle             !angle or angdeg
      pmenu(2+mpp(na))=0.                !gap or hgap
      pmenu(3+mpp(na))=0.5               !fint
      pmenu(4+mpp(na))=0.                !b
      pmenu(5+mpp(na))=0.                !lfrn
      pmenu(6+mpp(na))=0.                !tfrn
cryne 12/21/2004 need to add slices to MaryLie element parsed in MAD format
      if(nrp(1,2).eq.7)then
c       if(idproc.eq.0)write(6,*)'initializing nbnd slices to 1'
        pmenu(7+mpp(na))=1.
      endif
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=angdeg
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning: nbnd bend angle not found'
        endif
      endif
c
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(4+mpp(na))=brho*angle/bufr(1)
        else
          write(6,*)'Error: nbnd L or B not found'
          stop
        endif
      endif
c
c gap, field integral:
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(2+mpp(na))=2.d0*bufr(1)
      endif
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c leading fringe, trailing fringe:
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,2).eq.6)then
         write(6,*)' NBND input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         if(numpres)then
c          if(idproc.eq.0)write(6,*)'resetting nbnd slices to ',bufr(1)
           pmenu(7+mpp(na))=bufr(1)
         endif
        endif
      endif
      return
c
c=====================================================================
c                                PBND
c=====================================================================
c  'pbnd':
103   continue
c defaults:
      if(initparms.eq.1)then
      angle=0.
      pmenu(1+mpp(na))=angle
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.5
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=1.   !slices
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=bufr(1)
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning: pbnd bend angle not found'
        endif
      endif
c
c gap, field integral:
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(2+mpp(na))=2.d0*bufr(1)
      endif
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(4+mpp(na))=2.d0*brho*sin(0.5d0*angle)/bufr(1)
        else
          write(6,*)'Error: pbnd L or B not found'
          stop
        endif
      endif
c slices:
cryne 3/17/2004 this is not quite right.
cryne need to clarify how to handle mixed MaryLie and SIF input
cryne and implement it correctly throughout. This will do for now:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,3).eq.4)then
         write(6,*)' PBND input warning: SLICES will be ignored.'
         write(6,*)'To autoslice elements specified in the original'
         write(6,*)'MaryLie input style, place an autoslice element'
         write(6,*)'in the input file BEFORE any MaryLie elements'
        else
         if(numpres)pmenu(5+mpp(na))=bufr(1)
        endif
      endif
      return
c
c=====================================================================
c                                GBND
c=====================================================================
c 'gbnd    ': general bending magnet
104   continue
c defaults:
c 6 Marylie parameters = angdeg, e1deg, e2deg, gap, fint, B field
      if(initparms.eq.1)then
      angle=0.d0
      pmenu(1+mpp(na))=angle
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.5
      pmenu(6+mpp(na))=0.
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=angdeg
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning: gbnd bend angle not found'
        endif
      endif
c
c entry angle:
      call getparm(line,mmax,'e1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e1deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(2+mpp(na))=bufr(1)
        else
          write(6,*)'warning: gbnd entry angle not found'
        endif
      endif
c
c exit angle:
      call getparm(line,mmax,'e2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e2deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(3+mpp(na))=bufr(1)
        else
          write(6,*)'warning: gbnd exit angle not found'
        endif
      endif
c gap and field integral:
c
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=2.d0*bufr(1)
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(6+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          write(6,*)'Reading gbnd parameters. Performing conversion'
          write(6,*)'from length to B field assuming sbend geometry'
          pmenu(6+mpp(na))=brho*angle/bufr(1)
        else
          write(6,*)'Error: gbnd L or B not found'
          stop
        endif
      endif
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,4).eq.6)then
         write(6,*)' GBND input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         pmenu(7+mpp(na))=1.
         if(numpres)pmenu(7+mpp(na))=bufr(1)
        endif
      endif
      return
c=====================================================================
c                                PROT
c=====================================================================
c
c 'prot    ': rotation of reference plane
105   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)/pi180
      call getparm(line,mmax,'kind=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                GBDY
c=====================================================================
c 'gbdy    ': body of a general bending magnet
106   continue
c defaults:
      if(initparms.eq.1)then
      angle=0.
      pmenu(1+mpp(na))=angle
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=angdeg
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning: gbnd bend angle not found'
        endif
      endif
c
c entry angle:
      call getparm(line,mmax,'e1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e1deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(2+mpp(na))=bufr(1)
        else
          write(6,*)'warning: gbnd entry angle not found'
        endif
      endif
c
c exit angle:
      call getparm(line,mmax,'e2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e2deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(3+mpp(na))=bufr(1)
        else
          write(6,*)'warning: gbnd exit angle not found'
        endif
      endif
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          write(6,*)'Reading gbdy parameters. Performing conversion'
          write(6,*)'from length to B field assuming sbend geometry'
          pmenu(4+mpp(na))=brho*angle/bufr(1)
        else
          write(6,*)'Error: gbnd L or B not found'
          stop
        endif
      endif
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,6).eq.4)then
         write(6,*)' GBDY input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         pmenu(5+mpp(na))=1.
         if(numpres)pmenu(5+mpp(na))=bufr(1)
        endif
      endif
      return
c
c=====================================================================
c                                FRNG
c=====================================================================
c 'frng    ': hard edge dipole fringe fields
107   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.5
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=2.d0*bufr(1)
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'iedge=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      return
c
c=====================================================================
c           CFBD  :  PROBABLY NOT FUNCTIONAL; USE 'SBEND' INSTEAD
c=====================================================================
c 'cfbd    ': combined function bend (normal entry and exit)
108   continue
      write(6,*)'*****CFBD input is not functional*************'
      write(6,*)'USE SBEND INSTEAD'
cryne 12/17/2004 what was I thinking here??? this element is all screwed up!
c combined function bend with arbitrary entrance/exit angles?
      call getparm(line,mmax,'e1=',bufr,keypres1,numpres,0,cbuf)
      call getparm(line,mmax,'e2=',bufr,keypres1,numpres,0,cbuf)
      if(keypres1.or.keypres2)then
!       (make sure nt1,nt2,na are present in commons)
!       nt1(na)=???
!       nt2(na)=???
        goto 1045
      endif
c defaults:
      if(initparms.eq.1)then
      angle=0.
      pmenu(1+mpp(na))=angle
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.
      pmenu(6+mpp(na))=0.
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(numpres)pmenu(1+mpp(na))=bufr(1)
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
          angle=bufr(1)
          if(numpres)pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'WARNING: no bend angle specified for cfbd'
        endif
      endif
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(numpres)pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
          if(numpres)pmenu(2+mpp(na))=brho*bufr(1)/angle
        else
          write(6,*)'WARNING: no bfield or length specified for cfbd'
        endif
      endif
c leading fringe, trailing fringe, iopt, ipset:
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'ipset=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                SBEND
c=====================================================================
c 'sbend    ': combined function bend (arbitrary entry and exit)
1045  continue
c defaults (in MaryLie units, i.e. degrees for angles) :
c (1)bend angle in degrees, (2)B field,
c (3)entry angle in degrees, (4)exit angle in degrees,
c (5)entry fringe, (6)exit fringe, [DEFAULT FOR BOTH = 3, i.e. turned on]
c (7)gap size for leading (entry) fringe field calculation,
c (8)gap size for trailing (exit) fringe field calculation,
c (9)normalized field integral for leading (entry) edge fringe field,
c (10)normalized field integral for trailing (exit) edge fringe field,
c (11)iopt[interpretation of multipole coeffs;default=3=same meaning as MAD],
c (12)ipset [if multipole coeffs specified in pset],
c (13-18)=P1,P2,P3,P4,P5,P6
c        =BQD, AQD,BSEX,ASEX,BOCT,AOCT if iopt=1
c        =Tay1,AQD,Tay2,ASEX,Tay3,AOCT if iopt=2
c        =Tay1/brho,AQD/brho,Tay2/brho,ASEX/brho,Tay3/brho,AOCT/brho if iopt=3
c (19)axial rotation angle ["TILT" in MAD]
c (20)order
c (21)number of slices
      if(initparms.eq.1)then
c The phrases "entry angle" and "exit angle" are WRONG!!!!!!!!!!!
c They should be "entry pole face rotation angle" and "exit pole face
c rotation angle" (This is what E1 and E2 mean in MAD notation.)
c These are only comments in the code, but they will cause confusion.
c FIX LATER!!!!!!!!!
      angle=0.
      pmenu(1+mpp(na))=angle
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=3.  !fixed on 11/3/02 (formerly default was 1)
      pmenu(6+mpp(na))=3.  !ditto
      pmenu(7+mpp(na))=0.d0
      pmenu(8+mpp(na))=0.d0
      pmenu(9+mpp(na))=0.d0
      pmenu(10+mpp(na))=0.d0
      pmenu(11+mpp(na))=3.d0
      pmenu(12+mpp(na))=0.
      pmenu(13+mpp(na))=0.
      pmenu(14+mpp(na))=0.
      pmenu(15+mpp(na))=0.
      pmenu(16+mpp(na))=0.
      pmenu(17+mpp(na))=0.
      pmenu(18+mpp(na))=0.
      pmenu(19+mpp(na))=0.
      pmenu(20+mpp(na))=5. ! order
      pmenu(21+mpp(na))=1. ! slices
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=angdeg
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning: sbend bend angle not found'
        endif
      endif
c
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(2+mpp(na))=brho*angle/bufr(1)
        else
          write(6,*)'Error: sbend L or B not found'
          stop
        endif
      endif
c
c entry angle:
      call getparm(line,mmax,'e1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e1deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(3+mpp(na))=bufr(1)
        else
c         write(6,*)'warning: sbend entry angle not found'
        endif
      endif
c
c exit angle:
      call getparm(line,mmax,'e2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e2deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(4+mpp(na))=bufr(1)
        else
c         write(6,*)'warning: sbend exit angle not found'
        endif
      endif
c
c leading fringe, trailing fringe:
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c gap sizes:
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(7+mpp(na))=bufr(1)
        pmenu(8+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=2.d0*bufr(1)
        if(keypres.and.numpres)pmenu(8+mpp(na))=2.d0*bufr(1)
c
        call getparm(line,mmax,'gap1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
        call getparm(line,mmax,'gap2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
c
        call getparm(line,mmax,'hgap1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=2.d0*bufr(1)
        call getparm(line,mmax,'hgap2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(8+mpp(na))=2.d0*bufr(1)
      endif
c     write(6,*)'sbend gap/hgap: found the following:'
c     write(6,*)pmenu(7+mpp(na)),pmenu(8+mpp(na))
c normalized field integrals:
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(9+mpp(na))=bufr(1)
        pmenu(10+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'fint1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
        call getparm(line,mmax,'fint2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
      endif
c check consistency of input parameters regarding fringe fields:
      if( (pmenu(7+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
       write(6,*)'warning(sbend): input specifies a nonzero gap size'
       write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(8+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
       write(6,*)'warning(sbend): input specifies a nonzero gap size'
       write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
      if( (pmenu(9+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
      write(6,*)'warning(sbend):input specifies nonzero field integral'
      write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(10+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
      write(6,*)'warning(sbend):input specifies nonzero field integral'
      write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
      if( (pmenu(3+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
      write(6,*)'warning(sbend):input specifies nonzero entrance angle'
      write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(4+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
      write(6,*)'warning(sbend):input specifies nonzero exit angle'
      write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
c
c
c iopt, ipset:
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
      call getparm(line,mmax,'ipset=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
c multipole coefficients:
      call getparm(line,mmax,'k1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(13+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'s1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(14+mpp(na))=bufr(1)
c for compatability w/ MAD:
      call getparm(line,mmax,'ks=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(14+mpp(na))=bufr(1)
      call getparm(line,mmax,'k2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(15+mpp(na))=bufr(1)
      call getparm(line,mmax,'s2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(16+mpp(na))=bufr(1)
      call getparm(line,mmax,'k3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(17+mpp(na))=bufr(1)
      call getparm(line,mmax,'s3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(18+mpp(na))=bufr(1)
c tilt:
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(19+mpp(na))=bufr(1)
c order:
      call getparm(line,mmax,'order=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(20+mpp(na))=bufr(1)
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(21+mpp(na))=bufr(1)
c for debugging:
      if(idproc.eq.0 .and. iverbose.ge.1)then
      write(6,*)'done reading SBEND parameters,'
      write(6,*)'angle,b=',pmenu(1+mpp(na)),pmenu(2+mpp(na))
      write(6,*)'e1,e2=',pmenu(3+mpp(na)),pmenu(4+mpp(na))
      endif
      return
c
c
c=====================================================================
c                                QUAD
c=====================================================================
c 'quad    ': quadrupole
109   continue
      if(idproc.eq.0)then
      write(6,*)'reading parameters for a QUAD in MAD format.'
      write(6,*)'**NOTE WELL** QUAD is used for backward compatibility'
      write(6,*)'with MaryLie and corresponds to the original MaryLie'
      write(6,*)'code. Suggest you use QUADRUPOLE instead.'
      endif
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
cryne 12/17/2004 added initialization of  # of slices too,
cryne            commented out below at end of QUAD section
      pmenu(5+mpp(na))=1.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g1=',bufr,keypres1,numpres,0,cbuf)
      if(keypres1.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k1=',bufr,keypres2,numpres,0,cbuf)
        if(keypres2.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho
        endif
      endif
      if((.not.keypres1).and.(.not.keypres2))then
      write(6,*)lmntname,'WARNING: neither g1 nor k1 has been specified'
      endif
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,9).eq.4)then
         write(6,*)' QUAD input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
cryne 12/17/2004         pmenu(5+mpp(na))=1.
         if(numpres)pmenu(5+mpp(na))=bufr(1)
        endif
      endif
      write(6,*)'done reading QUAD parameters:'
      write(6,*)pmenu(1+mpp(na)),pmenu(2+mpp(na))
      write(6,*)pmenu(3+mpp(na)),pmenu(4+mpp(na))
      return
c
c=====================================================================
c                                SEXT
c=====================================================================
c 'sext    ': sextupole
110   continue
c     write(6,*)'(LMNTPARM) SEXT:'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho/2.d0
        endif
      endif
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,10).eq.2)then
         write(6,*)' SEXT input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         pmenu(3+mpp(na))=1.
         if(numpres)pmenu(3+mpp(na))=bufr(1)
        endif
      endif
      return
c
c=====================================================================
c                                OCTM
c=====================================================================
c 'octm    ': mag. octupole
111   continue
c     write(6,*)'(LMNTPARM) OCTM:'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k3=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho/6.d0
        endif
      endif
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(nrp(1,11).eq.2)then
         write(6,*)' SEXT input warning: SLICES will be ignored.'
         write(6,*)'To autoslice individual elements specified in the'
         write(6,*)'original MaryLie input style, use the SLICEMARYLIE'
         write(6,*)'argument of AUTOSLICE and place first in the menu'
        else
         pmenu(3+mpp(na))=1.
         if(numpres)pmenu(3+mpp(na))=bufr(1)
        endif
      endif
      return
c
c=====================================================================
c                                OCTE
c=====================================================================
c 'octe    ': elec. octupole
112   continue
      if(initparms.eq.1)then
c     ...
      endif
      write(6,*)'MAD-style input for MaryLie octe not implemented'
      stop
c
c=====================================================================
c                                SRFC
c=====================================================================
c 'srfc    ': short rf cavity
113   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.    !volts (in Volts)
      pmenu(2+mpp(na))=0.    !freq  (in Hz)
      endif
      call getparm(line,mmax,'volts=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'freq=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                RFGAP
c=====================================================================
c 'rfgap  ': rf gap
1135  continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=-1.0  !length
      pmenu(2+mpp(na))=0.    !frequency
      pmenu(3+mpp(na))=0.    !escale
      pmenu(4+mpp(na))=0.    !phasedeg
      pmenu(5+mpp(na))=0.    !unit (integer N for data file called rfdataN)
      pmenu(6+mpp(na))=0.    !steps
      pmenu(7+mpp(na))=0.    !e1deg
      pmenu(8+mpp(na))=0.    !notused
      pmenu(9+mpp(na))=1.    !tdim
      pmenu(10+mpp(na))=1.   !ptdim
      pmenu(11+mpp(na))=1.   !slices
      cmenu(1+mppc(na))=' '  !file = name of data file
      cmenu(2+mppc(na))=' '  !wake
      endif
c     write(6,*)'READING RFGAP PARAMETERS'
c length (or flag), freq,volts, phase, #, int_steps, slices
c
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'freq=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'volts=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(idproc.eq.0)then
        write(6,*)'ERROR (RFGAP): the keywork VOLTS has been replaced'
        write(6,*)'with ESCALE. Modify your input file and re-run'
        endif
        call myexit
      endif
      call getparm(line,mmax,'escale=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'phasedeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'unit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'steps=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'e1deg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'notused=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'tdim=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'ptdim=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
c file:
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
c wake:
      call getparm(line,mmax,'wake=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))='wake'
        write(6,*)'found rfgap wake; name of wake =',cbuf
        call wake_init(cbuf)
      endif
! check for over-specification:
      if(pmenu(5+mpp(na)).ne.0. .and. cmenu(1+mppc(na)).ne.' ')then
      if(idproc.eq.0)then
        write(6,*)'rfgap input error:'
        write(6,*)'specified file rfgapN, where N =',pmenu(5+mpp(na))
        write(6,*)'specified file name, where name=',cmenu(1+mppc(na))
        write(6,*)'cannot specify file by both number and name'
      endif
      call myexit
      endif
      return
c
c=====================================================================
c                                CONFOC
c=====================================================================
c 'confoc  ': "constant focusing" element   rdr 08/29/2001
1136  continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !length
      pmenu(2+mpp(na))=0.d0  !k1
      pmenu(3+mpp(na))=0.d0  !k2
      pmenu(4+mpp(na))=0.d0  !k3
      pmenu(5+mpp(na))=1.d0  !slices
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'k1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'k2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'k3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
!
      if(idproc.eq.0)then
      write(6,*)'results of confoc input:'
      write(6,*)'l=',pmenu(1+mpp(na))
      write(6,*)'k1=',pmenu(2+mpp(na))
      write(6,*)'k2=',pmenu(3+mpp(na))
      write(6,*)'k3=',pmenu(4+mpp(na))
      write(6,*)'slices=',pmenu(5+mpp(na))
      endif
      return
c
c=====================================================================
c                                AROT
c=====================================================================
c 'arot    ': axial rotation
114   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)/pi180
      return
      stop
c
c=====================================================================
c                                TWSM
c=====================================================================
c 'twsm    ': linear matrix via twiss parameters
115   continue
      if(initparms.eq.1)then
c     ...
      endif
      write(6,*)'MAD-style input for ML twsm not implemented'
      stop
c
c=====================================================================
c                                THLM
c=====================================================================
c 'thlm    ': thin lens low order multipole
116   continue
c     if(idproc.eq.0)write(6,*)'reading parameters for thlm'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.
      pmenu(6+mpp(na))=0.
      endif
c values provided by user:
c normal quadrupole (n=0):
      call getparm(line,mmax,'bqd=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
c       if(idproc.eq.0)write(6,*)'found bqd'
         if(numpres)then
          pmenu(1+mpp(na))=bufr(1)
c         if(idproc.eq.0)write(6,*)'param1=',pmenu(1+mpp(na))
         endif
      else
        call getparm(line,mmax,'k1l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
c       if(idproc.eq.0)write(6,*)'found k1l'
         if(numpres)then
           pmenu(1+mpp(na))=bufr(1)*brho
c          if(idproc.eq.0)write(6,*)'param1=',pmenu(1+mpp(na))
         endif
        endif
      endif
c skew quadrupole:
c fix later.
c
c normal sextupole:
      call getparm(line,mmax,'bsex=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
c       if(idproc.eq.0)write(6,*)'found bsex'
        if(numpres)then
          pmenu(3+mpp(na))=bufr(1)
c         if(idproc.eq.0)write(6,*)'param3=',pmenu(3+mpp(na))
        endif
      else
        call getparm(line,mmax,'k2l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
c       if(idproc.eq.0)write(6,*)'found k2l'
         if(numpres)then
           pmenu(3+mpp(na))=bufr(1)*brho
c          if(idproc.eq.0)write(6,*)'param3=',pmenu(3+mpp(na))
         endif
        endif
      endif
c skew sextupole:
c fix later.
c
c normal octupole:
      call getparm(line,mmax,'boct=',bufr,keypres,numpres,0,cbuf)
      if(keypres)then
        if(numpres)pmenu(5+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k3l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
          if(numpres)pmenu(5+mpp(na))=bufr(1)*brho
        endif
      endif
c skew sextupole:
c fix later.
c higher order stuff from MAD:
c fix later.
c printout for debugging:
c     if(idproc.eq.0)then
c       write(6,*)'multipole parameters have been set to:'
c       do ii=1,6
c         write(6,*)pmenu(ii+mpp(na))
c       enddo
c     endif
      return
c
c=====================================================================
c                                CPLM
c=====================================================================
c 'cplm    ': "compressed" low order multipole
117   continue
      write(6,*)'MAD-style input for MaryLie cplm not implemented'
      stop
c
c=====================================================================
c                                CFQD
c=====================================================================
c 'cfqd    ': combined function quadrupole
118   continue
      write(6,*)'MAD-style input for MaryLie cfqd not implemented'
      stop
c
c dispersion matrix (dism)
119   continue
      write(6,*)'MAD-style input for MaryLie dism not implemented'
      stop
c
c 'sol     ': solenoid
120   continue
      write(6,*)'MAD-style input for MaryLie sol not implemented'
      write(6,*)'Use solenoid instead of sol'
      stop
c
c 'mark    ': marker
121   continue
      write(6,*)'MAD-style input for MaryLie mark not implemented'
      stop
c
c 'jmap    ': j mapping
122   continue
      write(6,*)'MAD-style input for MaryLie jmap not implemented'
      stop
c
c 'dp      ': data point
123   continue
      write(6,*)'MAD-style input for MaryLie dp not implemented'
      stop
c
c 'recm    ': REC multiplet
124   continue
      write(6,*)'MAD-style input for MaryLie recm not implemented'
      stop
c
c 'spce    ': space
125   continue
      write(6,*)'MAD-style input for MaryLie spce not implemented'
      stop
c
c 'cfrn    ': change/write fringe field params for comb function dipole
126   continue
      write(6,*)'MAD-style input for MaryLie cfrn not implemented'
      stop
c
c 'coil'
127   continue
      write(6,*)'MAD-style input for MaryLie coil not implemented'
      stop
c
c 'intg'
128   continue
      write(6,*)'MAD-style input for MaryLie intg not implemented'
      stop
c
129   continue
c 'rmap'
      write(6,*)'MAD-style input for MaryLie rmap not implemented'
      stop
c
c 'arc'
130   continue
      write(6,*)'MAD-style input for MaryLie arc not implemented'
      stop
c
c 'transit'
133   continue
      write(6,*)'(sif) transit'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.      !l (length)
      pmenu(2+mpp(na))=1.      !n
      pmenu(3+mpp(na))=1.      !slices
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'n=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      return
c 'interface'
134   write(6,*)'(sif) interface'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.      !n1
      pmenu(2+mpp(na))=1.      !n2
      pmenu(3+mpp(na))=0.      !b2
      pmenu(4+mpp(na))=0.      !b4
      pmenu(5+mpp(na))=0.      !b6
      endif
c values provided by user:
      call getparm(line,mmax,'n1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'n2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'b2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'b4=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'b6=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      return
c 'rootmap'
135   write(6,*)'(sif) rootmap'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.      !n
      pmenu(2+mpp(na))=0.      !b2
      pmenu(3+mpp(na))=0.      !b4
      pmenu(4+mpp(na))=0.      !b6
      cmenu(1+mppc(na))='false'    !invert= : flag to compute inverted map
      endif
c values provided by user:
      call getparm(line,mmax,'n=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'b2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'b4=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'b6=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'inverse=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      return
c 'optirot'
136   write(6,*)'(sif) optirot'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)/pi180
      call getparm(line,mmax,'kind=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'spare5'
137   write(6,*)'spare5 just a placeholder (not a valid type code)'
      return
c 'spare6'
138   write(6,*)'spare6 just a placeholder (not a valid type code)'
      return
c 'spare7'
139   write(6,*)'spare7 just a placeholder (not a valid type code)'
      return
c 'spare8'
140   write(6,*)'spare8 just a placeholder (not a valid type code)'
      return
c
c 'marker': MAD marker
141   continue
cryne 5/4/2006      write(6,*)'ignoring input for MAD marker (not implemented)'
      write(12,*)'ignoring input for MAD marker (not implemented)'
      return
c
c=====================================================================
c                                DRIFT
c=====================================================================
c 'drift': MAD drift
142   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.      !isssflag
      pmenu(3+mpp(na))=1.      !slices
      cmenu(1+mppc(na))=' '     !wake
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)
      endif
c wake:
      call getparm(line,mmax,'wake=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))='wake'
        write(6,*)'found drift wake; name of wake =',cbuf
        call wake_init(cbuf)
      endif
c     write(6,*)'FOUND DRIFT; SLICES=',nint(bufr(1))
      return
c
c=====================================================================
c                                RBEND
c=====================================================================
c 'rbend': MAD rectangular bend
143   continue
c defaults:
      if(initparms.eq.1)then
cryne NOTE TO MYSELF:
c The phrases "entry angle" and "exit angle" are WRONG!!!!!!!!!!!
c They should be "entry pole face rotation angle" and "exit pole face
c rotation angle" (This is what E1 and E2 mean in MAD notation.)
c These are only comments in the code, but they will cause confusion.
c FIX LATER!!!!!!!!!
      angle=0.
      pmenu(1+mpp(na))=angle   !bend angle
      pmenu(2+mpp(na))=0.      !B field
      pmenu(3+mpp(na))=0.      !entry angle
      pmenu(4+mpp(na))=0.      !exit angle
      pmenu(5+mpp(na))=3.  !entry fringe on/off (default on)
      pmenu(6+mpp(na))=3.  !exit fringe on/off (default on)
      pmenu(7+mpp(na))=0.d0 !gap size for leading (entry) fringe field calc
      pmenu(8+mpp(na))=0.d0 !gap size for trailing (exit) fringe field calc
      pmenu(9+mpp(na))=0.d0  !fint for leading (entry) fringe
      pmenu(10+mpp(na))=0.d0 !fint size for trailing (exit) fringe
      pmenu(11+mpp(na))=0. ! tilt
      pmenu(12+mpp(na))=5. ! order
      pmenu(13+mpp(na))=1. ! slices
      endif
c values provided by user:
c angle:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        angdeg=bufr(1)
        pmenu(1+mpp(na))=angdeg
        angle=angdeg*pi180
      else
        call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          angle=bufr(1)
          pmenu(1+mpp(na))=angle/pi180
        else
          write(6,*)'warning (rbend): bend angle not found'
        endif
      endif
c
c b field:
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          if(angle.ne.0.)then
            pmenu(2+mpp(na))=brho/(0.5*bufr(1)/sin(0.5*angle))
          else
            write(6,*)'error (rbend): angle=0 in b field calculation'
            call myexit
          endif
        else
          write(6,*)'Error: rbend L or B not found'
          stop
        endif
      endif
c
c entry angle:
      call getparm(line,mmax,'e1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e1deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(3+mpp(na))=bufr(1)
        else
c         write(6,*)'warning (rbend): entry angle not found'
        endif
      endif
c
c exit angle:
      call getparm(line,mmax,'e2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)/pi180
      else
        call getparm(line,mmax,'e2deg=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(4+mpp(na))=bufr(1)
        else
c         write(6,*)'warning (rbend): exit angle not found'
        endif
      endif
c
c leading fringe, trailing fringe:
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c gap sizes:
      call getparm(line,mmax,'gap=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(7+mpp(na))=bufr(1)
        pmenu(8+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'hgap=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=2.d0*bufr(1)
        if(keypres.and.numpres)pmenu(8+mpp(na))=2.d0*bufr(1)
c
        call getparm(line,mmax,'gap1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
        call getparm(line,mmax,'gap2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
c
        call getparm(line,mmax,'hgap1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(7+mpp(na))=2.d0*bufr(1)
        call getparm(line,mmax,'hgap2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(8+mpp(na))=2.d0*bufr(1)
      endif
c     write(6,*)'rbend gap/hgap: found the following:'
c     write(6,*)pmenu(7+mpp(na)),pmenu(8+mpp(na))
c normalized field integrals:
      call getparm(line,mmax,'fint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(9+mpp(na))=bufr(1)
        pmenu(10+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'fint1=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
        call getparm(line,mmax,'fint2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
      endif
c check consistency of input parameters regarding fringe fields:
      if( (pmenu(7+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
       write(6,*)'warning (rbend): input specifies a nonzero gap size'
       write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(8+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
       write(6,*)'warning (rbend): input specifies a nonzero gap size'
       write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
      if( (pmenu(9+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
      write(6,*)'warning (rbend):input specifies nonzero field integral'
      write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(10+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
      write(6,*)'warning (rbend):input specifies nonzero field integral'
      write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
      if( (pmenu(3+mpp(na)).ne.0.) .and. (pmenu(5+mpp(na)).eq.0.) )then
      write(6,*)'warning (rbend):input specifies nonzero entrance angle'
      write(6,*)'but leading-edge fringe effects are turned off?'
      endif
      if( (pmenu(4+mpp(na)).ne.0.) .and. (pmenu(6+mpp(na)).eq.0.) )then
      write(6,*)'warning (rbend):input specifies nonzero exit angle'
      write(6,*)'but trailing-edge fringe effects are turned off?'
      endif
c
c
c tilt:
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
c order:
      call getparm(line,mmax,'order=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
c slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(13+mpp(na))=bufr(1)
c
c     write(6,*)'done reading RBEND parameters for element: ',lmntname
c     write(6,*)'bend angle=',pmenu(1+mpp(na))
c     write(6,*)'B field=',pmenu(2+mpp(na))
c     write(6,*)'entry angle=',pmenu(3+mpp(na))
c     write(6,*)'exit angle=',pmenu(4+mpp(na))
c     write(6,*)'lfrn=',pmenu(5+mpp(na))
c     write(6,*)'tfrn=',pmenu(6+mpp(na))
c     write(6,*)'gap1=',pmenu(7+mpp(na))
c     write(6,*)'gap2=',pmenu(8+mpp(na))
c     write(6,*)'fint1=',pmenu(9+mpp(na))
c     write(6,*)'fint2=',pmenu(10+mpp(na))
c     write(6,*)'tilt=',pmenu(11+mpp(na))
c     write(6,*)'order=',pmenu(12+mpp(na))
c     write(6,*)'slices=',pmenu(13+mpp(na))
      return
c
c 'sbend': MAD sector bend
c 144   continue
c see statement 1045
c
c 'gbend': MAD general bend
145   continue
      write(6,*)'MAD general bend not implemented'
      stop
c
c=====================================================================
c                                QUADRUPOLE
c=====================================================================
c 'quadrupole': MAD quadrupole
146   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.     !isssflag
      pmenu(6+mpp(na))=1.     !slices
      cmenu(1+mppc(na))=' '    !wake
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g1=',bufr,keypres1,numpres,0,cbuf)
      if(keypres1.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k1=',bufr,keypres2,numpres,0,cbuf)
        if(keypres2.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho
        endif
      endif
      if((.not.keypres1).and.(.not.keypres2))then
      write(6,*)lmntname,'WARNING: neither g1 nor k1 has been specified'
      endif
      call getparm(line,mmax,'lfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfrn=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(6+mpp(na))=bufr(1)
      endif
c wake:
      call getparm(line,mmax,'wake=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))='wake'
        write(6,*)'found quadrupole wake; name of wake =',cbuf
        call wake_init(cbuf)
      endif
      return
c
c
c=====================================================================
c                             SEXTUPOLE
c=====================================================================
c 'sextupol': MAD sextupole
147   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=1.     !slices
      cmenu(1+mppc(na))=' '    !wake
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k2=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho/2.d0
        endif
      endif
c
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c wake:
      call getparm(line,mmax,'wake=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))='wake'
        write(6,*)'found sextupole wake; name of wake =',cbuf
        call wake_init(cbuf)
      endif
      return
c
c=====================================================================
c                             OCTUPOLE
c=====================================================================
c 'octupole': MAD octupole
148   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=1.     !slices
      cmenu(1+mppc(na))=' '    !wake
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'g3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      else
        call getparm(line,mmax,'k3=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)*brho/6.d0
        endif
      endif
c
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c wake:
      call getparm(line,mmax,'wake=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))='wake'
        write(6,*)'found octupole wake; name of wake =',cbuf
        call wake_init(cbuf)
      endif
      return
c
c=====================================================================
c                                MULTIPOLE
c=====================================================================
c 'multipole': MAD general thin multipole
c NOTE WELL: no factors of brho here; dealt with in subroutine lmnt
149   continue
      if(multmsg.eq.1)then
      if(idproc.eq.0)write(6,*)'MAD multipole partially implemented'
      multmsg=0
      endif
c     if(idproc.eq.0)write(6,*)'reading parameters for MAD multipole'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=0.
      pmenu(6+mpp(na))=0.
      endif
c values provided by user:
c normal quadrupole (n=0):
        call getparm(line,mmax,'k1l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
c       if(idproc.eq.0)write(6,*)'found k1l'
         if(numpres)then
           pmenu(1+mpp(na))=bufr(1)
c          if(idproc.eq.0)write(6,*)'param1=',pmenu(1+mpp(na))
         endif
        endif
c skew quadrupole:
c fix later.
c
c normal sextupole:
      call getparm(line,mmax,'bsex=',bufr,keypres,numpres,0,cbuf)
        call getparm(line,mmax,'k2l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
c       if(idproc.eq.0)write(6,*)'found k2l'
         if(numpres)then
           pmenu(3+mpp(na))=bufr(1)
c          if(idproc.eq.0)write(6,*)'param3=',pmenu(3+mpp(na))
         endif
        endif
c skew sextupole:
c fix later.
c
c normal octupole:
        call getparm(line,mmax,'k3l=',bufr,keypres,numpres,0,cbuf)
        if(keypres)then
          if(numpres)pmenu(5+mpp(na))=bufr(1)
        endif
c skew sextupole:
c fix later.
c higher order stuff from MAD:
c fix later.
c printout for debugging:
      if(idproc.eq.0)then
        write(12,*)'(sif) multipole parameters read in, set to:'
        do ii=1,6
          write(12,*)pmenu(ii+mpp(na))
        enddo
      endif
      return
c
c=====================================================================
c
c 'solenoid': MAD solenoid
150   continue
c defaults:
      if (initparms.eq.1) then
        pmenu(1+mpp(na))=   0.  ! zstart
        pmenu(2+mpp(na))=   0.  ! zend
        pmenu(3+mpp(na))= 100.  ! steps
        pmenu(4+mpp(na))=   0.  ! iprofile
        pmenu(5+mpp(na))=  -1.  ! ipset (not used, but -1 ==> SIF-style)
        pmenu(6+mpp(na))=   0.  ! multipoles (not used)
        pmenu(7+mpp(na))=   0.  ! ldrift
        pmenu(8+mpp(na))=   0.  ! lbody
        pmenu(9+mpp(na))=   0.  ! clength
        pmenu(10+mpp(na))=  0.  ! b
        pmenu(11+mpp(na))=  0.  ! iecho
        pmenu(12+mpp(na))=  0.  ! iopt
        pmenu(13+mpp(na))=  1.  ! slices
      end if
c values provided by user:
c first make sure the user is not using a pset:
      call getparm(line,mmax,'ipset=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        if(idproc.eq.0)then
          write(6,*)'error: solenoid with MAD-style input does not'
          write(6,*)'require specification of ipset'
        endif
        stop
      endif
c
      call getparm(line,mmax,'zstart=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'zend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'steps=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'iprofile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
ccccc call getparm(line,mmax,'ipset=',bufr,keypres,numpres,0,cbuf)
ccccc if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'multipoles=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'ldrift=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'lbody=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      call getparm(line,mmax,'clength=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
      call getparm(line,mmax,'b=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
      call getparm(line,mmax,'iecho=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(13+mpp(na))=bufr(1)
      return
c
c 'hkicker': MAD horizontal kicker
151   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.      !isssflag
      pmenu(4+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'kick=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      return
c
c 'vkicker': MAD vertical kicker
152   continue
c defaults:
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.      !isssflag
      pmenu(4+mpp(na))=0.
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'kick=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      return
c
c 'kicker': MAD kicker
153   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.      !isssflag
      pmenu(5+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'hkick=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'vkick=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      return
c
c 'rfcavity': MAD rf cavity
154   continue
      write(6,*)'MAD RFCAVITY'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.   ! l
      pmenu(2+mpp(na))=0.   ! volt
      pmenu(3+mpp(na))=0.   ! lag
      pmenu(4+mpp(na))=0.   ! harmon
      pmenu(5+mpp(na))=0.   ! betrf
      pmenu(6+mpp(na))=0.   ! pg
      pmenu(7+mpp(na))=0.   ! shunt
      pmenu(8+mpp(na))=0.   ! tfill
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'volt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'lag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'harmon=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'betrf=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'pg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'shunt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'tfill=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      return
c
c 'elsepara': MAD electrostatic separator
155   continue
      write(6,*)'MAD ELSEPARATOR'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.      !isssflag
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'e=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'tilt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      return
c
c 'hmonitor': MAD horizontal monitor
156   continue
      write(12,*)'reading SIF input for MAD hmonitor'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.      !length
      pmenu(2+mpp(na))=0.      !isssflag
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'vmonitor': MAD vertical monitor
157   continue
      write(12,*)'reading SIF input for MAD vmonitor'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.      !length
      pmenu(2+mpp(na))=0.      !isssflag
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'monitor': MAD monitor
158   continue
      write(12,*)'reading SIF input for MAD monitor'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.      !length
      pmenu(2+mpp(na))=0.      !isssflag
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'instrume': MAD instrument
159   continue
      write(12,*)'reading SIF input for MAD instrument'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.      !length
      pmenu(2+mpp(na))=0.      !isssflag
      endif
c values provided by user:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'sssflag=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'sparem1'
160   write(6,*)'sparem1 just a placeholder (not a valid type code)'
      return
c
c 'rcollima': MAD rectangular collimator
161   continue
      write(6,*)'ignoring input for MAD rcollimator (not implemented)'
      return
c
c 'ecollima': MAD elliptical collimator
162   continue
      write(6,*)'ignoring input for MAD ecollimator (not implemented)'
      return
c
c 'yrot'
163   write(6,*)'YROT input: using MaryLie prot'
      write(6,*)'NOTE WELL: this will not work!'
      write(6,*)'YROT has one parameter; MaryLie prot has two.'
      write(6,*)'this needs to be fixed.'
cryne goto 105
      stop
c
c 'srot'
164   write(6,*)'SROT input: using MaryLie arot (check sign!)'
      goto 114
c
c=====================================================================
c                                PROT3
c=====================================================================
c
c 'prot3    ': 3rd order rotation of reference plane
165   continue
      if(idproc.eq.0)write(6,*)'reading prot3 data in sif.f'
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=0.
      endif
c values provided by user:
      call getparm(line,mmax,'angdeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'angle=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)/pi180
      call getparm(line,mmax,'kind=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c 'beambeam': MAD beam-beam element
166   write(6,*)'MAD beambeam not implemented'
      stop
c
c 'matrix': MAD matrix command
167   write(6,*)'MAD matrix not implemented'
      stop
c
c=====================================================================
c                               PROFILE1
c=====================================================================
c 'profile1d': 1D profile monitor
168   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.d0     !column for 1D histrogram
      pmenu(2+mpp(na))=128.d0   !number of bins
      pmenu(3+mpp(na))=0.d0     !sequencelength (=max # of files in sequence)
      pmenu(4+mpp(na))=5.d0     !precision
      pmenu(5+mpp(na))=0.d0 !assigned unit# for output file (or can be read in)
      pmenu(6+mpp(na))=0.d0     !assigned counter for sequence of files
      pmenu(7+mpp(na))=0.d0     !rwall
      cmenu(1+mppc(na))=' '     !file= : name of output file
      cmenu(2+mppc(na))='true'  !close= : flag to close output file
      cmenu(3+mppc(na))='false' !flush= : flag to flush output file
      endif
c values provided by user:
      call getparm(line,mmax,'column=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'bins=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'sequencelength=',bufr,keypres,numpres,      &
     &             0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,           &
     &             0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'unit=',bufr,keypres,numpres,                &
     &             0,cbuf)
      if(keypres.and.numpres)then
        pmenu(5+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'close=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'rwall=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
!     if(idproc.eq.0)then
!       write(6,*)'done reading data for profile1'
!       write(6,*)'column=',pmenu(1+mpp(na))
!       write(6,*)'bins=',pmenu(2+mpp(na))
!     endif
      return
c 'sparem4'
169   write(6,*)'sparem4 just a placeholder (not a valid type code)'
      return
c 'sparem5'
170   write(6,*)'sparem5 just a placeholder (not a valid type code)'
      return
c
c 'hkick  ': synonym for MAD horizontal kicker
171   continue
      goto 151
c
c 'vkick  ': synonym for MAD vertical kicker
172   continue
      goto 152
c
c 'kick  ': synonym for MAD kicker
173   continue
      goto 153
c
c 'sparem6'
174   write(6,*)'sparem6 just a placeholder (not a valid type code)'
      return
c
c=====================================================================
c                                NLRF
c=====================================================================
c 'nlrf': nonlinear rf cavity
175   continue
      if(initparms.eq.1)then
        pmenu(1+mpp(na))=0.d0  !zstart
        pmenu(2+mpp(na))=0.d0  !zend
        pmenu(3+mpp(na))=0.d0  !frequency
        pmenu(4+mpp(na))=0.d0  !phasedeg
        pmenu(5+mpp(na))=1.d0  !escale
        pmenu(6+mpp(na))=0.d0  !steps
        pmenu(7+mpp(na))=1.d0  !slices
        cmenu(1+mppc(na))='rfdata'  !E-field data file
        cmenu(2+mppc(na))='crz.dat'  !generalized gradients file
      endif
c values provided by user:
cp1
      call getparm(line,mmax,'zstart=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
cp2
      call getparm(line,mmax,'zend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
cp3
      call getparm(line,mmax,'frequency=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
cp4
      call getparm(line,mmax,'phasedeg=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
cp5
      call getparm(line,mmax,'escale=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
cp6
      call getparm(line,mmax,'nz=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
cp7
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
cc1
      call getparm(line,mmax,'efile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
cc2
      call getparm(line,mmax,'crzfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
c
      if(idproc.eq.0)then
        write(6,*) 'Read nlrf parameters:'
        write(6,*) '  zstart=',pmenu(1+mpp(na))
        write(6,*) '  zend=',pmenu(2+mpp(na))
        write(6,*) '  frequency=',pmenu(3+mpp(na))
        write(6,*) '  phasedeg=',pmenu(4+mpp(na))
        write(6,*) '  escale=',pmenu(5+mpp(na))
        write(6,*) '  steps=',pmenu(6+mpp(na))
        write(6,*) '  slices=',pmenu(7+mpp(na))
        write(6,*) '  efile=',cmenu(1+mppc(na))
        write(6,*) '  crzfile=',cmenu(2+mppc(na))
      endif
      return
c
c
c     2: user-supplied elements ************************************
c
c           'usr1    ','usr2    ','usr3    ','usr4    ','usr5    ',
12    go to (201,       202,       203,       204,       205,           &
c           'usr6    ','usr7    ','usr8    ','usr9    ','usr10   '
     &       206,       207,       208,       209,       210,           &
c           'usr11   ','usr12   ','usr13   ','usr14   ','usr15   ',
     &       211,       212,       213,       214,       215,           &
c           'usr16   ','usr17   ','usr18   ','usr19   ','usr20   ',
     &       216,       217,       218,       219,       220),kt2
c
201   continue
      return
202   continue
      return
203   continue
      return
204   continue
      return
205   continue
      return
206   continue
      return
207   continue
      return
208   continue
      return
209   continue
      return
210   continue
      return
211   continue
      return
212   continue
      return
213   continue
      return
214   continue
      return
215   continue
      return
216   continue
      return
217   continue
      return
218   continue
      return
219   continue
      return
220   continue
      return
c
c     3: parameter sets **********************************************
c
c           'ps1     ','ps2     ','ps3     ','ps4     ','ps5     ',
13    go to (301,       302,       303,       304,       305,           &
c           'ps6     ','ps7     ','ps8     ','ps9     '/
     &       306,       307,       308,       309),kt2
c
301   continue
      return
302   continue
      return
303   continue
      return
304   continue
      return
305   continue
      return
306   continue
      return
307   continue
      return
308   continue
      return
309   continue
      return
c
c     4-6: random elements
c
14    continue
15    continue
16    continue
      return
c
c     7: simple commands *************************************
c
c           'rt      ','sqr     ','symp    ','tmi     ','tmo     ',
17    go to (701,       702,       703,       704,       705,
c           'pmif    ','circ    ','stm     ','gtm     ','end     ',
     &       706,       707,       708,       709,       710,
c           'ptm     ','iden    ','whst    ','inv     ','tran    ',
     &       711,       712,       713,       714,       715,
c           'revf    ','rev     ','mask    ','num     ','rapt    ',
     &       716,       717,       718,       719,       720,
c           'eapt    ','of      ','cf      ','wnd     ','wnda    ',
     &       721,       722,       723,       724,       725,
c           'ftm     ','wps     ','time    ','cdf     ','bell    ',
     &       726,       727,       728,       729,       730,
c           'wmrt    ','wcl     ','paws    ','inf     ','dims    ',
     &       731,       732,       733,       734,       735,
c           'zer     ','sndwch  ','tpol    ','dpol    ','cbm     ',
     &       736,       737,       738,       739,       740,
c           'poisson ','preapply','midapply','autoapply','autoconc',
     &       741,       742,       743,       744,       745,
c           'rayscale','beam    ','units   ','autoslic','verbose ',
     &       746,       747,       748,       749,       750,
c           'mask6   ','arcreset','symbdef ','particledump','raytrace',
     &       751,       752,       753,       754,       755,
c           'autotrack','sckick','moments  ','maxsize','reftraj',
     &       756,       757,       758,       759,       760,
c           'initenv','envelopes','contractenv','setreftraj','setarclen',
     &       761,       762,       763,       764,       765,
c           'wakedefault','emittance','matchenv','fileinfo','egengrad',
     &       766,       767,       768,       769,       770,
c           'wrtmap','rdmap','sparec7','sparec8','sparec9'/
     &       771,       772,       773,       774,       775),kt2


c
c 'rt': ray trace
701   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=13.
      pmenu(2+mpp(na))=14.
      pmenu(3+mpp(na))=5.
      pmenu(4+mpp(na))=1.
      pmenu(5+mpp(na))=1.
      pmenu(6+mpp(na))=0.
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))=' '
      endif
c values provided by user:
      call getparm(line,mmax,'unit1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'order=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'ntrace=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'nwrite=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'ibrief=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'file1=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file2=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      return
c
c     square the existing map:
c
702   continue
      return
c
c     symplectify matrix in transfer map
c
703   continue
      return
c
c     input transfer map from an external file:
c
704   continue
      return
c
c     output transfer map to an external file (tmo):
c
705   continue
      return
c
c 'pmif':   print contents of file master input file:
706   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.
      pmenu(2+mpp(na))=12.
      pmenu(3+mpp(na))=3.
      cmenu(1+mppc(na))=' '
      endif
c values provided by user:
      call getparm(line,mmax,'itype=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      return
c
707   continue
      return
c
c     store the existing transfer map
c
708   continue
      return
c
c     get transfer map from storage
c
709   continue
      return
c
c     end of job:
c
710   continue
c defaults:
      cmenu(1+mppc(na))='false'  !timers=
c values provided by user:
      call getparm(line,mmax,'timers=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      return
c
c 'ptm':  print transfer map:
711   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=3.
      pmenu(2+mpp(na))=3.
      pmenu(3+mpp(na))=0.
      pmenu(4+mpp(na))=0.
      pmenu(5+mpp(na))=1.
      endif
c values provided by user:
      call getparm(line,mmax,'matrix=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'poly=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'t2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'u3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'basis=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      return
c
c     identity mapping:
c
712   continue
      return
c
c     write history of beam loss
c
713   continue
      return
c
c     inverse:
c
714   continue
      return
c
c     transpose:
c
715   continue
      return
c
c     reverse factorization:
c
716   continue
      return
c
c     Dragt's reversal
c
717   continue
      return
c
c     mask off selected portions of transfer map:
c
718   continue
c defaults:
      if(initparms.eq.1)then
c default it to not mask anything:
      pmenu(1+mpp(na))=1.   !f1
      pmenu(2+mpp(na))=1.   !matrix and f2
      pmenu(3+mpp(na))=1.   !f3
      pmenu(4+mpp(na))=1.   !f4
      pmenu(5+mpp(na))=1.   !f5
      pmenu(6+mpp(na))=1.   !f6
      endif
c values provided by user:
      call getparm(line,mmax,'f1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'f2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'f3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'f4=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'f5=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'f6=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c "matrix=" means the same thing as "f2="
      call getparm(line,mmax,'matrix=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      return
c
c     number lines in a file
c
719   continue
      return
c
c     aperture particle distribution
c
720   continue
      return
c
721   continue
      return
c
c     open files
c
722   continue
      return
c
c     close files
c
723   continue
      return
c
c     window particle distribution
c
724   continue
      return
725   continue
      return
c
c     filter transfer map
c
726   continue
      return
c
c     write parameter set
c
727   continue
      return
c
c     write time
c
728   continue
      return
c
c 'cdf': change output drop file
729   continue
c default:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=12.
      endif
c values provided by user:
      call getparm(line,mmax,'ifile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c     ring bell
c
730   continue
      return
c
c     write value of merit function
c
731   continue
      return
c
c     write contents of loop
c
732   continue
      return
c
c     pause (paws)
c
733   continue
      return
c
c     change or write out infinities (inf)
c
734   continue
      return
c
c get dimensions (dims)
c
735   continue
      return
c
c     change or write out values of zeroes (zer)
c
736   continue
      return
c
c sndwch
c
737   continue
      return
c
c     twiss polynomial (tpol)
c
738   continue
      return
c
c     dispersion polynomial (dpol)
c
739   continue
      return
c
c     change or write out beam parameters (cbm)
c
740   continue
      return
c
c     set parameters for poisson solver
c=====================================================================
c                               POISSON
c=====================================================================
741   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.  !nx
      pmenu(2+mpp(na))=0.  !ny
      pmenu(3+mpp(na))=0.  !nz
      pmenu(4+mpp(na))=0.  !xmin
      pmenu(5+mpp(na))=0.  !xmax
      pmenu(6+mpp(na))=0.  !ymin
      pmenu(7+mpp(na))=0.  !ymax
      pmenu(8+mpp(na))=0.  !zmin
      pmenu(9+mpp(na))=0.  !zmax
      pmenu(10+mpp(na))=0. !anag_patchsize
      pmenu(11+mpp(na))=0. !anag_refineratio
      cmenu(1+mppc(na))='fft'  !solver (fft or fft2 or chombo)
      cmenu(2+mppc(na))='rectangular'  !geometry [not currently used]
      cmenu(3+mppc(na))='variable' !gridsize (fixed or variable)
      cmenu(4+mppc(na))='fixed'    !gridpoints (fixed or variable)
      cmenu(5+mppc(na))='open'  !xboundary (open, dirichlet, or periodic)
      cmenu(6+mppc(na))='open'  !yboundary
      cmenu(7+mppc(na))='open'  !zboundary
      cmenu(8+mppc(na))='open'  !boundary
      cmenu(9+mppc(na))='E'  !solving_for (phi, E)
      cmenu(10+mppc(na))='delta' !densityfunction (delta or linear)
      cmenu(11+mppc(na))='undefined' !chombo_file
      cmenu(12+mppc(na))='none'      !anag_smooth
      cmenu(13+mppc(na))='undefined' !spare_13
      cmenu(14+mppc(na))='undefined' !spare_14
      cmenu(15+mppc(na))='undefined' !spare_15
      cmenu(16+mppc(na))='undefined' !spare_16
      cmenu(17+mppc(na))='undefined' !spare_17
      endif
c values provided by user:
c nx:
      call getparm(line,mmax,'nx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
c ny:
      call getparm(line,mmax,'ny=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
c nz:
      call getparm(line,mmax,'nz=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c xmin:
      call getparm(line,mmax,'xmin=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
c xmax:
      call getparm(line,mmax,'xmax=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c ymin:
      call getparm(line,mmax,'ymin=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c ymax:
      call getparm(line,mmax,'ymax=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
c zmin:
      call getparm(line,mmax,'zmin=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
c zmax:
      call getparm(line,mmax,'zmax=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
c anag_patchsize:
      call getparm(line,mmax,'anag_patchsize=',bufr,keypres,numpres,0
     &            ,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
c anag_refineratio:
      call getparm(line,mmax,'anag_refineratio=',bufr,keypres,numpres,0
     &            ,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
c solver:
      call getparm(line,mmax,'solver=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      endif
c geometry:
      call getparm(line,mmax,'geometry=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
c gridsize:
      call getparm(line,mmax,'gridsize=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(3+mppc(na))=cbuf
      endif
c gridpoints:
      call getparm(line,mmax,'gridpoints=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(4+mppc(na))=cbuf
      endif
c xboundary::
      call getparm(line,mmax,'xboundary=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(5+mppc(na))=cbuf
      endif
c yboundary::
      call getparm(line,mmax,'yboundary=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(6+mppc(na))=cbuf
      endif
c zboundary::
      call getparm(line,mmax,'zboundary=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(7+mppc(na))=cbuf
      endif
c boundary:
      call getparm(line,mmax,'boundary=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(8+mppc(na))=cbuf
      endif
c solving_for:
      call getparm(line,mmax,'solving_for=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(9+mppc(na))=cbuf
      endif
c densityfunction:
      call getparm(line,mmax,'densityfunction=',bufr,keypres,           &
     &             numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(10+mppc(na))=cbuf
      endif
c chombo_file:
      call getparm(line,mmax,'chombo_file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(11+mppc(na))=cbuf
      endif
c anag_smooth:
      call getparm(line,mmax,'anag_smooth=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(12+mppc(na))=cbuf
      endif
c spare_13:
      call getparm(line,mmax,'spare_13=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(13+mppc(na))=cbuf
      endif
c spare_14:
      call getparm(line,mmax,'spare_14=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(14+mppc(na))=cbuf
      endif
c spare_15:
      call getparm(line,mmax,'spare_15=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(15+mppc(na))=cbuf
      endif
c spare_16:
      call getparm(line,mmax,'spare_16=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(16+mppc(na))=cbuf
      endif
c spare_17:
      call getparm(line,mmax,'spare_17=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(17+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               PREAPPLY
c=====================================================================
c     preapply commands automatically
c
742   continue
      if(initparms.eq.1)then
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))='physical' !applyto=  'physical' or 'all' elements
      endif
      call getparm(line,mmax,'name=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'applyto=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               MIDAPPLY
c=====================================================================
c     midapply commands automatically
c
743   continue
      if(initparms.eq.1)then
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))='physical' !applyto=  'physical' or 'all' elements
      endif
      call getparm(line,mmax,'name=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'applyto=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               AUTOAPPLY
c=====================================================================
c     autoapply commands automatically
c
744   continue
      if(initparms.eq.1)then
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))='physical' !applyto=  'physical' or 'all' elements
      endif
      call getparm(line,mmax,'name=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'applyto=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
c
      if(idproc.eq.0 .and. iverbose.ge.1)then
       write(6,*)'(routine lmntparm) results from POST-APPLY command:'
       write(6,*)'cmenu(1+mppc(na))=',cmenu(1+mppc(na))
      endif
      return
c
c     autoconc ('auto-concatenate')
c=====================================================================
c                                AUTOCONC
c=====================================================================
c
745   continue
c default:
      if(initparms.eq.1)then
      cmenu(1+mppc(na))='true' !set
      endif
c
      call getparm(line,mmax,'set=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      else
        call getparm(line,mmax,'true',bufr,keypres,numpres,1,cbuf)
        if(keypres)cmenu(1+mppc(na))='true'
        call getparm(line,mmax,'false',bufr,keypres,numpres,1,cbuf)
        if(keypres)cmenu(1+mppc(na))='false'
      endif
      if(idproc.eq.0)then
        write(6,*)'(sif) result of autoconc input:'
        write(6,*)'set=',cmenu(1+mppc(na))
      endif
      return
c
c     rayscale ('scale zblock array')
c=====================================================================
c                                RAYSCALE
c=====================================================================
c
746   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.d0
      pmenu(2+mpp(na))=1.d0
      pmenu(3+mpp(na))=1.d0
      pmenu(4+mpp(na))=1.d0
      pmenu(5+mpp(na))=1.d0
      pmenu(6+mpp(na))=1.d0
      endif
c note: this code allows the use of a 'div' option instead of
c a multiplicative factor. However, when the div option is used
c the data are immediately converted to multiplicative factors,
c which is how they are stored in the pmenu array.
c
c values provided by user:
c X:
      call getparm(line,mmax,'xmult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'xdiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=1.d0/bufr(1)
      endif
c PX:
      call getparm(line,mmax,'pxmult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'pxdiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=1.d0/bufr(1)
      endif
c Y:
      call getparm(line,mmax,'ymult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'ydiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(3+mpp(na))=1.d0/bufr(1)
      endif
c PY:
      call getparm(line,mmax,'pymult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'pydiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(4+mpp(na))=1.d0/bufr(1)
      endif
c T:
      call getparm(line,mmax,'tmult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(5+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'tdiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(5+mpp(na))=1.d0/bufr(1)
      endif
c PT:
      call getparm(line,mmax,'ptmult=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(6+mpp(na))=bufr(1)
      endif
      call getparm(line,mmax,'ptdiv=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(6+mpp(na))=1.d0/bufr(1)
      endif
      return
c
c=====================================================================
c                                BEAM
c=====================================================================
c  beam: set parameters for an ensemble of particles
747   continue
c     write(6,*)'******************************here I am at beam'
c defaults:
      if(initparms.eq.1)then
      maxray=0
!     energy=1.d9
!     pmass=938.27200d6
!     gamma=energy/pmass
!     gamm1=gamma-1.d0
!     beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
!     brho=gamma*beta/clite*pmass
      energy=0.d0
      pmass=0.d0
      gamma=1.d0
      gamm1=0.d0
      beta=0.d0
      brho=0.d0
      achg=1.d0
      bcurr=0.d0
      bfreq=0.d0
      bomega=twopi*bfreq
!
      icompmass=0
!
      pmenu(1+mpp(na))=maxray
      pmenu(2+mpp(na))=brho
      pmenu(3+mpp(na))=gamma
      pmenu(4+mpp(na))=gamm1
      pmenu(5+mpp(na))=beta
      pmenu(6+mpp(na))=pmass
      pmenu(7+mpp(na))=achg
      pmenu(8+mpp(na))=bcurr
      pmenu(9+mpp(na))=bfreq
      pmenu(10+mpp(na))=bomega
      cmenu(1+mppc(na))='proton'
      endif
c values provided by user:
c maximum number of particles:
      call getparm(line,mmax,'maxray=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        maxray=nint(bufr(1))
        pmenu(1+mpp(na))=maxray
c       allocate space for the particle array, etc:
        if(.not.allocated(zblock))call new_particledata
      endif
c particle:
      call getparm(line,mmax,'particle=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
        if(cbuf.eq.'proton')then
          achg=1.d0
          pmass=938.27231d6
        endif
        if(cbuf.eq.'H+')then
          achg=1.d0
          pmass=939.29400d6
        endif
        if(cbuf.eq.'electron')then
          achg=-1.d0
          pmass=0.511d6
        endif
        if(cbuf.eq.'positron')then
          achg=1.d0
          pmass=0.511d6
        endif
      endif
c mass (read in as GeV, same as MAD):
      call getparm(line,mmax,'mass=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmass=bufr(1)*1.d9
        pmenu(6+mpp(na))=pmass
      endif
c charge:
      call getparm(line,mmax,'charge=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        achg=bufr(1)
        pmenu(7+mpp(na))=achg
      endif
c energy:
      call getparm(line,mmax,'energy=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        energy=bufr(1)*1.d9
        if(pmass.ne.0.d0)then
          gamma=energy/pmass
          gamm1=gamma-1.d0
          beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
          brho=gamma*beta/clite*pmass
          pc=sqrt(energy**2-pmass**2)
          ekinetic=energy-pmass
        endif
      endif
c pc:
      call getparm(line,mmax,'pc=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pc=bufr(1)*1.d9
        if(pmass.ne.0.d0)then
          energy=sqrt(pc**2+pmass**2)
          gamma=energy/pmass
          gamm1=gamma-1.d0
          beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
          brho=gamma*beta/clite*pmass
          ekinetic=energy-pmass
        endif
      endif
c gamma:
      call getparm(line,mmax,'gamma=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        gamma=bufr(1)
        gamm1=gamma-1.d0
        beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
        if(pmass.ne.0.d0)then
          brho=gamma*beta/clite*pmass
          energy=gamma*pmass
          pc=sqrt(energy**2-pmass**2)
          ekinetic=energy-pmass
        endif
      endif
c gamma-1:
      call getparm(line,mmax,'gamma1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        gamm1=bufr(1)
        gamma=gamm1+1.d0
        beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
        if(pmass.ne.0.d0)then
          brho=gamma*beta/clite*pmass
          energy=gamma*pmass
          pc=sqrt(energy**2-pmass**2)
          ekinetic=energy-pmass
        endif
      endif
c ekinetic:
      call getparm(line,mmax,'ekinetic=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        ekinetic=bufr(1)*1.d9
        if(pmass.ne.0.d0)then
          gamm1=ekinetic/pmass
          gamma= gamm1+1.d0
          beta=sqrt((gamma+1.d0)*(gamma-1.d0))/gamma
          brho=gamma*beta/clite*pmass
          energy=ekinetic+pmass
          pc=sqrt(energy**2-pmass**2)
        endif
      endif
c brho:
      call getparm(line,mmax,'brho=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        brho=bufr(1)
        if(gamma.ne.1.d0 .and. beta.ne.0.d0)then
          pmass=brho/(gamma*beta/clite)
          icompmass=1
!         if(idproc.eq.0)then
!         write(6,*)'computed pmass (from brho,gamma,beta) = ',pmass
!         endif
          energy=gamma*pmass
          pc=sqrt(energy**2-pmass**2)
          ekinetic=energy-pmass
        endif
      endif
c store brho,gamma,gamm1,beta in pmenu
      pmenu(2+mpp(na))=brho
      pmenu(3+mpp(na))=gamma
      pmenu(4+mpp(na))=gamm1
      pmenu(5+mpp(na))=beta
c beam current:
      call getparm(line,mmax,'bcurr=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        bcurr=bufr(1)
        pmenu(8+mpp(na))=bcurr
      endif
c beam frequency:
      call getparm(line,mmax,'bfreq=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        bfreq=bufr(1)
        pmenu(9+mpp(na))=bfreq
        bomega=twopi*bfreq
      endif
c beam angular frequency:
      call getparm(line,mmax,'bomega=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        bomega=bufr(1)
        pmenu(10+mpp(na))=bomega
        bfreq=bomega/twopi
      endif
c print results:
      if(idproc.eq.0)then
      write(6,*)' '
      write(6,*)'done reading BEAM parameters:'
      write(6,'(1x,''brho='',1pd23.16)')brho
      write(6,                                                              &
     &'(1x,''gamma, gamma-1='',1pd23.16,2x,1pd23.16)')gamma,gamm1
      write(6,'(1x,''beta='',1pd23.16)')beta
      write(6,                                                              &
     &'(1x,''energy, ekinetic='',1pe18.11,2x,1pe18.11)')energy,ekinetic
      write(6,'(1x,''pc='',1pe18.11)')pc
      write(6,'(1x,''mass, charge='',1pe18.11,2x,1pe18.11)')pmass,achg
      write(6,'(1x,''beam current='',1pe18.11)')bcurr
      write(6,'(1x,''beam frequency='',1pe18.11)')bfreq
c     write(6,*)'beam angular freq=',bomega
c     write(6,*)'array size=',maxray
c     write(6,*)'symbol(1)=',cmenu(1+mppc(na))
!     write(6,*)' '
      if(bcurr.ne.0.d0 .and. bfreq.eq.0.d0)then
      write(6,*)'**********************WARNING************************'
      write(6,*)                                                            &
     &'Beam current has been specfied, but the beam rf frequency that'
      write(6,*)                                                            &
     &'defines the total bunch charge, Q=I/freq, has not been specified'
      write(6,*)'*****************************************************'
      endif
!     write(6,*)' '
      endif
cryne 09/20/02
      constr(nconst+1)='brho'
      conval(nconst+1)=brho
      constr(nconst+2)='gamma'
      conval(nconst+2)=gamma
      constr(nconst+3)='charge'
      conval(nconst+3)=charge
      constr(nconst+4)='mass'
      conval(nconst+4)=pmass
      constr(nconst+5)='beta'
      conval(nconst+5)=beta
      constr(nconst+6)='energy'
      conval(nconst+6)=energy
      constr(nconst+7)='ekinetic'
      conval(nconst+7)=ekinetic
      constr(nconst+8)='pc'
      conval(nconst+8)=pc
      nconst=nconst+8
cryne
c---------
c also set the default units in case the user forgets to use the
c units command later:
!Jan23, 2003      sl=1.d0
!Jan23, 2003      ts=sl/clite
!Jan23, 2003      omegascl=1.d0/ts
!Jan23, 2003      p0sc=pc
c other stuff:
!Jan23, 2003      magunits=1
!Jan23, 2003      iverbose=0
!jan 23, 2003 set p0sc if it has not been set by the user already:
!jan 24, 2003      if(p0sc.eq.0.d0)p0sc=pc
!!!!!      if(p0sc.eq.0.d0)p0sc=pc/clite
!!!!!      if(idproc.eq.0)write(6,*)'p0sc=pc/clite=',p0sc
!!!!!      if(idproc.eq.0)write(6,*)'note also that pc=',pc
!     if(idproc.eq.0 )then
!     write(6,*)'default scaling variables follow. (Note:'
!     write(6,*)'default scaling can be changed using the units command)'
!     write(6,*)'scale length, sl=',sl
!     write(6,*)'scale time, ts=',ts
!     write(6,*)'scale omega, omegascl=',omegascl
!     write(6,*)'scale momentum, p0sc=',pc
!     write(6,*)'end of reporting of data for BEAM command'
!     endif
c---------
c this does not make sense if the default value for maxray is zero,
c so I have commented it out.
c allocate space for the particle array, etc:
c     if(.not.allocated(zblock))call new_particledata
c
cryne March 18, 2004
c store ref particle data in the default location (location #1)
c in case the user wants to do multiple runs from the same initial values:
c [note: initial reftraj data, refsave(1:6), is set/stored in subroutine tran]
      brhosav(1)=brho
      gamsav(1)=gamma
      gam1sav(1)=gamm1
      betasav(1)=beta

! Global defautls for wakefields
      call readin(line,leof)
      if(line(1:9) .ne. 'wakedflt:')then
!cryne         write(6,*)'\nWakefield global defaults are not specified'
!cryne         write(6,*)'You may specify wakefield global defaults'
!cryne         write(6,*)'in the following format'
!cryne         write(6,*)'wakedflt: type=  nturns=  nmodes=  r=  conduct= '
!cryne         write(6,*)'nx=  ny=  nz=  ndx=  ndy=   \n'
         backspace lf
      else
         call wake_defaults(line)
      endif

      return
c
c
c=====================================================================
c                                UNITS
c=====================================================================
c  units: set scaling variables
748   continue
c defaults:
!Jan23, 2003: now these are set in routine dumpin
!Jan23, 2003      sl=1.d0
!Jan23, 2003      ts=sl/clite
!Jan23, 2003      omegascl=1.d0/ts
!Jan23, 2003      p0sc=pc
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0    !scale length
      pmenu(2+mpp(na))=0.d0    !scale momentum
      pmenu(3+mpp(na))=0.d0    !scale time (sec)
      pmenu(4+mpp(na))=0.d0    !scale ang freq (rad/sec)
      pmenu(5+mpp(na))=0.d0    !scale freq (Hz)
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))=' '
      endif
      lflagmagu=.true.
      lflagdynu=.false.
      itscaleset=0
      ilscaleset=0
      if(idproc.eq.0)then
      write(6,*)' '
      write(6,*)'reading UNITS parameters'
      endif
c values provided by user:
c type (magnetic, static, or general):
      call getparm(line,mmax,'type=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        if( (cbuf.eq.'magnetic') .or. (cbuf.eq.'static') )then
          p0sc=gamma*beta*pmass/clite
          pmenu(2+mpp(na))=p0sc
          cmenu(1+mppc(na))='static'
          cmenu(2+mppc(na))='p0'
          lflagmagu=.true.
          lflagdynu=.false.
          if(idproc.eq.0)then
          write(6,*)'STATIC UNITS WILL BE USED. p0sc=',p0sc
          endif
        elseif(cbuf.eq.'dynamic')then
          p0sc=pmass/clite
          pmenu(2+mpp(na))=p0sc
          cmenu(1+mppc(na))='dynamic'
          cmenu(2+mppc(na))='mc'
          lflagmagu=.false.
          lflagdynu=.true.
          if(idproc.eq.0)then
          write(6,*)'DYNAMIC UNITS WILL BE USED. p0sc=',p0sc
          endif
        else
          cmenu(1+mppc(na))='general'
          cmenu(2+mppc(na))='none'
          lflagmagu=.false.
          lflagdynu=.false.
          if(idproc.eq.0)write(6,*)'GENERAL UNITS WILL BE USED.'
        endif
      endif
!
! momentum:
      if(cmenu(1+mppc(na)).ne.'static' .and.                               &
     &   cmenu(1+mppc(na)).ne.'magnetic' .and.                             &
     &   cmenu(1+mppc(na)).ne.'dynamic')then
        call getparm(line,mmax,'p=',bufr,keypres,numpres,0,cbuf)
        if(keypres.and.numpres)then
          pmenu(2+mpp(na))=bufr(1)
          p0sc=bufr(1)
        endif
        call getparm(line,mmax,'psymbolic=',bufr,keypres,numpres,1,cbuf)
        if(keypres.and.numpres)then
          cmenu(2+mppc(na))=cbuf
          if(trim(cbuf).eq.'p0')then
            pmenu(2+mpp(na))=pc
            p0sc=pc
          endif
          if(trim(cbuf).eq.'mc')then
            pmenu(2+mpp(na))=pmass/clite
            p0sc=pmass/clite
          endif
        endif
        if(pmenu(2+mpp(na)).eq.0.d0)then
          if(idproc.eq.0)then
            write(6,*)'ERROR: general units are being used but the'
            write(6,*)'scale momentum has not been specified; halting'
            call myexit
          endif
        endif
      endif

! length:
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=bufr(1)
        sl=bufr(1)
        ilscaleset=1
      else
        if(idproc.eq.0)then
        write(6,*)'WARNING: scale length has not been specified'
        endif
      endif
! time:
      call getparm(line,mmax,'w=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        omegascl=bufr(1)
        ts=1.d0/omegascl
        freqscl=omegascl/twopi
        itscaleset=1
      endif
      call getparm(line,mmax,'f=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        freqscl=bufr(1)
        omegascl=twopi*freqscl
        ts=1.d0/omegascl
        itscaleset=1
      endif
      call getparm(line,mmax,'t=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        ts=bufr(1)
        omegascl=1.d0/ts
        freqscl=omegascl/twopi
        itscaleset=1
      endif
      if(itscaleset.eq.0)then
        if(idproc.eq.0)then
          write(6,*)'WARNING: scale time has not been specified'
        endif
      endif
! choose some sensible defaults for data that have not been provided:
      if(itscaleset.eq.0 .and. ilscaleset.eq.0)then
!      if(idproc.eq.0)then
!      write(6,*)'neither the scale length nor the scale time have ',     &
!    &           'been specified'
!      endif
! static/magnetic:
       if(bfreq.eq.0.d0)then
         if(idproc.eq.0)write(6,*)'setting sl=1, omega=clite/sl'
         sl=1.d0
         omegascl=clite/sl
         ts=1.d0/omegascl
         freqscl=omegascl/twopi
       endif
! dynamic:
       if(bfreq.ne.0.d0)then
         if(idproc.eq.0)then
         write(6,*)'setting scalefreq=beamfreq, l=c/omega'
         endif
         freqscl=bfreq
         omegascl=twopi*freqscl
         ts=1.d0/omegascl
         sl=clite/omegascl
       endif
      elseif(ilscaleset.eq.0)then
        if(idproc.eq.0)then
        write(6,*)'The scale length has not been specified'
        write(6,*)'Setting sl=clite/omega=',clite/omegascl
        endif
        sl=clite/omegascl
      elseif(itscaleset.eq.0)then
        if(idproc.eq.0)then
        write(6,*)'The scale time has not been specified'
        write(6,*)'Setting omegascl=clite/sl=',clite/sl
        endif
        omegascl=clite/sl
        ts=1.d0/omegascl
        freqscl=omegascl/twopi
      endif
!
c print results:
      if(idproc.eq.0)then
      write(6,*)'done reading UNITS parameters:'
      write(6,'(1x,''scale length (m) ='',1pe18.11)')sl
      write(6,'(1x,''scale momentum ='',1pe18.11)')p0sc
      write(6,'(1x,''scale time (sec)='',1pe18.11)')ts
      write(6,'(1x,''scale freq (Hz)='',1pe18.11)')freqscl
      write(6,'(1x,''scale angular freq (rad/sec)='',1pe18.11)')omegascl
!     write(6,*)'symbol(1)=',cmenu(1+mppc(na))
!     write(6,*)'symbol(2)=',cmenu(2+mppc(na))
      write(6,*)' '
      endif
      return
c
c=====================================================================
c                               AUTOSLICE
c=====================================================================
c autoslice
c
749   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0       !slices=N
      pmenu(2+mpp(na))=0.d0       !l=
      cmenu(1+mppc(na))='local'    !control="local", "global", or "none"
      cmenu(2+mppc(na))='false'    !sckick="true" or "false"; determines
c                                 !if a space charge kick will be performed
c                                 !automatically in the middle of each slice
      cmenu(3+mppc(na))='false'    !includemlstyle="true" or "false"
      endif
c values provided by user:
!slices:
      call getparm(line,mmax,'slices=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=bufr(1)
      endif
!l (interval between slices):
      call getparm(line,mmax,'l=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(2+mpp(na))=bufr(1)
      endif
!control:
      localcontrol=0
      call getparm(line,mmax,'control=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
        if(cbuf.eq.'local')localcontrol=1
      endif
!sckick:
      call getparm(line,mmax,'sckick=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
!includemlstyle:
      inclmlstyle=0
      call getparm(line,mmax,'includemlstyle=',bufr,keypres,numpres,       &
     &             1,cbuf)
      if(keypres.and.numpres)then
        cmenu(3+mppc(na))=cbuf
        if(cbuf.eq.'true')inclmlstyle=1
      endif
!
! Dec 3, 2002:
! Augment parameter lists for old-style MaryLie input to include # of slices.
! Note that includemlstyle appears to be no longer needed; instead I just
! force it to be enabled if the user specifies local control.
! However, I have left the code in place in case it is needed later.
      if(localcontrol.eq.1 .or. inclmlstyle.eq.1)call sliceml
!
!     if(idproc.eq.0)then
!     write(6,*)'(routine lmntparm) results from AUTOSLICE command:'
!     write(6,*)'slices=',pmenu(1+mpp(na)),'l=',pmenu(2+mpp(na))
!     write(6,*)'control=',cmenu(1+mppc(na))
!!!!!!write(6,*)'sckick=',cmenu(2+mppc(na))
!     endif
      return
c
c=====================================================================
c                               VERBOSE
c=====================================================================
c autoslice
c
750   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0
      endif
c values provided by user:
      call getparm(line,mmax,'iverbose=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)then
        pmenu(1+mpp(na))=bufr(1)
        iverbose=nint(bufr(1))
      endif
      if(idproc.eq.0 .and. iverbose.ge.1)then
      write(6,*)'(routine lmntparm) results from VERBOSE command:'
      write(6,*)'iverbose=',iverbose
      endif
      return
c
c  mask6:
cryne 12/31/2004 mask6 now appears to be unecessary.
751   continue
      return
c
c  arcreset
752   continue
      return
c
c  symbdef
753   continue
      isymbdef=1
      if(idproc.eq.0)then
        write(6,*)'SYMBDEF: from now on symbolic names that appear in'
        write(6,*)'arithmetic expressions will have default value 0'
      endif
      return
c
c=====================================================================
c                               PARTICLEDUMP
c=====================================================================
c  particledump:
754   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0     !min= :  ...
      pmenu(2+mpp(na))=0.d0     !max= :  to print particles # min to max
      pmenu(3+mpp(na))=0.d0     !sequencelength (=max # of files in sequence)
      pmenu(4+mpp(na))=5.d0     !precision
      pmenu(5+mpp(na))=0.d0 !assigned unit# for output file (or can be read in)
      pmenu(6+mpp(na))=0.d0     !assigned counter for sequence of files
      pmenu(7+mpp(na))=1.d0  !nunits (=0 for dimensionless; =1 for physical)
      cmenu(1+mppc(na))=' '     !file= : name of output file
      cmenu(2+mppc(na))='true'  !close= : flag to close output file
      cmenu(3+mppc(na))='false' !flush= : flag to flush output file
      cmenu(4+mppc(na))='false' !printarc=:flag to print arc length in column 1
      endif
c values provided by user:
      call getparm(line,mmax,'min=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'max=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'sequencelength=',bufr,keypres,numpres,      &
     &             0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'close=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'printarc=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
!     if(idproc.eq.0)then
!     write(6,*)'(lmntparm) results from PARTICEDUMP sif input:'
!     write(6,*)'min=',pmenu(1+mpp(na))
!     write(6,*)'max=',pmenu(2+mpp(na))
!     write(6,*)'sequencelength=',pmenu(3+mpp(na))
!     write(6,*)'precision=',pmenu(4+mpp(na))
!     write(6,*)'unit=',pmenu(5+mpp(na))
!     write(6,*)'sequence counter=',pmenu(6+mpp(na))
!     write(6,*)'file=',cmenu(1+mppc(na))
!     write(6,*)'close=',cmenu(2+mppc(na))
!     write(6,*)'flush=',cmenu(3+mppc(na))
!     endif
      return
c
c
c=====================================================================
c                               RAYTRACE
c=====================================================================
c  raytrace:
755   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0   !min= :  ...
      pmenu(2+mpp(na))=0.d0   !max= :  to print particles # min to max
      pmenu(3+mpp(na))=5.     !norder
      pmenu(4+mpp(na))=1.     !ntrace
      pmenu(5+mpp(na))=1.     !nwrite
      pmenu(6+mpp(na))=0.     !ibrief
      pmenu(7+mpp(na))=0.     !sequencelength (=max # of files in sequence)
      pmenu(8+mpp(na))=6.     !precision
      pmenu(9+mpp(na))=0.     !assigned unit# for initial condition file
      pmenu(10+mpp(na))=0.    !assigned unit# for final condition file
      pmenu(11+mpp(na))=0.d0  !assigned counter for sequence of files
      pmenu(12+mpp(na))=0.d0  !nrays=# to read when reading from data file
      cmenu(1+mppc(na))=' '   !name of initial condition file
      cmenu(2+mppc(na))=' '   !name of final condition file
      cmenu(3+mppc(na))='false'   !flag to close final condition file
      cmenu(4+mppc(na))='false'   !flag to flush final condition file
      cmenu(5+mppc(na))='undefined' !alternate way to specify track type&order
      endif
c values provided by user:
      call getparm(line,mmax,'min=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'max=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'order=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'ibrief=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'sequencelength=',bufr,keypres,numpres,0,      &
     &             cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
      call getparm(line,mmax,'nrays=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'file1=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file2=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'close=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
c type:
      call getparm(line,mmax,'type=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(5+mppc(na))=cbuf
        if(cbuf.eq.'readonly')then
c       change defaults:
        pmenu(4+mpp(na))=0.     !ntrace
        pmenu(5+mpp(na))=0.     !nwrite
        endif
      endif
c ntrace, nwrite:
      call getparm(line,mmax,'ntrace=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'nwrite=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      return
c
c
c=====================================================================
c                               AUTOTRACK
c=====================================================================
c  autotrack
756   continue
c defaults:
      if(initparms.eq.1)then
      cmenu(1+mppc(na))='true'  !set
      cmenu(2+mppc(na))='undefined'  !type (taylorN, symplecticN, undefined)
      cmenu(3+mppc(na))='undefined'  !sckick (true,false)
      cmenu(4+mppc(na))='false'  !env (true,false)
      endif
c set:
      call getparm(line,mmax,'set=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      else
        call getparm(line,mmax,'true',bufr,keypres,numpres,1,cbuf)
        if(keypres)cmenu(1+mppc(na))='true'
        call getparm(line,mmax,'false',bufr,keypres,numpres,1,cbuf)
        if(keypres)cmenu(1+mppc(na))='false'
      endif
c type:
      call getparm(line,mmax,'type=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(2+mppc(na))=cbuf
      endif
c sckick:
      call getparm(line,mmax,'sckick=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(3+mppc(na))=cbuf
      endif
c env:
      call getparm(line,mmax,'env=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(4+mppc(na))=cbuf
      endif
c
c     if(idproc.eq.0)then
c       write(6,*)'(sif) results of autotrack input:'
c       write(6,*)'set= ',cmenu(1+mppc(na))
c       write(6,*)'type= ',cmenu(2+mppc(na))
c       write(6,*)'sckick= ',cmenu(3+mppc(na))
c       write(6,*)'env= ',cmenu(4+mppc(na))
c     endif
      return
c
c=====================================================================
c                               SCKICK
c=====================================================================
c  sckick
757   continue
      return
c
c=====================================================================
c                               MOMENTS
c=====================================================================
c moments:
758   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !assigned unit# for xfile
      pmenu(2+mpp(na))=0.d0  !assigned unit# for yfile
      pmenu(3+mpp(na))=0.d0  !assigned unit# for tfile
      pmenu(4+mpp(na))=5.d0  !precision
      pmenu(5+mpp(na))=1.d0  !nunits (=0 for dimensionless; =1 for physical)
      cmenu(1+mppc(na))='xrms.out'   !name of xfile
      cmenu(2+mppc(na))='yrms.out'   !name of yfile
      cmenu(3+mppc(na))='trms.out'   !name of tfile
      cmenu(4+mppc(na))='true'   !flag to flush xfile
      cmenu(5+mppc(na))='true'   !flag to flush yfile
      cmenu(6+mppc(na))='true'   !flag to flush tfile
      cmenu(7+mppc(na))='ratio'  !crossterm [ <xpx> or <xpx>/(xrms*pxrms)]
      cmenu(8+mppc(na))='false'  !set to true to divide emittances by pi
      cmenu(9+mppc(na))='keep'   !flag to keep/remove beam centroid
      endif
c values provided by user:
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'xfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'yfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'tfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'xflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'yflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(5+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'tflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(6+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'crossterm=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(7+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'includepi=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(8+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'centroid=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(9+mppc(na))=cbuf
      endif
      return
c
c
c
c=====================================================================
c                               MAXSIZE
c=====================================================================
c maxsize:
759   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !assigned unit# for output file 1
      pmenu(2+mpp(na))=0.d0  !assigned unit# for output file 2
      pmenu(3+mpp(na))=0.d0  !assigned unit# for output file 3
      pmenu(4+mpp(na))=0.d0  !assigned unit# for output file 4
      pmenu(5+mpp(na))=5.d0  !precision
      pmenu(6+mpp(na))=1.d0  !nunits (=0 for dimensionless; !=0 for physical)
      cmenu(1+mppc(na))=' '   !name of output file 1
      cmenu(2+mppc(na))=' '   !name of output file 2
      cmenu(3+mppc(na))=' '   !name of output file 3
      cmenu(4+mppc(na))=' '   !name of output file 4
      cmenu(5+mppc(na))='true'   !flag to flush output file 1
      cmenu(6+mppc(na))='true'   !flag to flush output file 2
      cmenu(7+mppc(na))='true'   !flag to flush output file 3
      cmenu(8+mppc(na))='true'   !flag to flush output file 4
      cmenu(9+mppc(na))='automatic'   !flag to use 'standard' file names
      endif
c values provided by user:
      call getparm(line,mmax,'unit1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit3=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit4=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'file1=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file2=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file3=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file4=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush1=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(5+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush2=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(6+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush3=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(7+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush4=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(8+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'files=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(9+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               REFTRAJ
c=====================================================================
c reftraj:
760   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !assigned unit# for output file
      pmenu(2+mpp(na))=5.d0  !precision
      pmenu(3+mpp(na))=1.d0  !nunits (=0 for dimensionless; !=0 for physical)
      cmenu(1+mppc(na))='reftraj.out'   !name of output file
      cmenu(2+mppc(na))='true'   !flag to flush output file
      endif
c values provided by user:
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      return
c=====================================================================
c                               INITENV
c=====================================================================
c initenv (initialize rms envelopes):
761   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  ! x=
      pmenu(2+mpp(na))=0.d0  ! px=
      pmenu(3+mpp(na))=0.d0  ! xpx=
      pmenu(4+mpp(na))=0.d0  ! xemit=
      pmenu(5+mpp(na))=0.d0  ! y=
      pmenu(6+mpp(na))=0.d0  ! py=
      pmenu(7+mpp(na))=0.d0  ! ypy=
      pmenu(8+mpp(na))=0.d0  ! yemit=
      pmenu(9+mpp(na))=0.d0  ! t=
      pmenu(10+mpp(na))=0.d0 ! pt=
      pmenu(11+mpp(na))=0.d0 ! tpt=
      pmenu(12+mpp(na))=0.d0 ! temit=
      pmenu(13+mpp(na))=0.d0 ! cpx=
      pmenu(14+mpp(na))=0.d0 ! cpy=
      pmenu(15+mpp(na))=0.d0 ! cpt=
      cmenu(1+mppc(na))='ratio'  !crossterm [ <xpx> or <xpx>/(xrms*pxrms)]
      endif
c values provided by user:
      call getparm(line,mmax,'x=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'px=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'xpx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'xemit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'y=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'py=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'ypy=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'yemit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      call getparm(line,mmax,'t=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
      call getparm(line,mmax,'pt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
      call getparm(line,mmax,'tpt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
      call getparm(line,mmax,'temit=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
      call getparm(line,mmax,'cpx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(13+mpp(na))=bufr(1)
      call getparm(line,mmax,'cpy=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(14+mpp(na))=bufr(1)
      call getparm(line,mmax,'cpt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(15+mpp(na))=bufr(1)
      call getparm(line,mmax,'crossterm=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               ENVELOPES
c=====================================================================
c envelopes:
762   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !assigned unit# for xfile
      pmenu(2+mpp(na))=0.d0  !assigned unit# for yfile
      pmenu(3+mpp(na))=0.d0  !assigned unit# for tfile
      pmenu(4+mpp(na))=5.d0  !precision
      pmenu(5+mpp(na))=1.d0  !nunits (=0 for dimensionless; =1 for physical)
      cmenu(1+mppc(na))='xenv.out'   !name of xfile
      cmenu(2+mppc(na))='yenv.out'   !name of yfile
      cmenu(3+mppc(na))='tenv.out'   !name of tfile
      cmenu(4+mppc(na))='true'   !flag to flush xfile
      cmenu(5+mppc(na))='true'   !flag to flush yfile
      cmenu(6+mppc(na))='true'   !flag to flush tfile
      cmenu(7+mppc(na))='ratio'  !crossterm [ <xpx> or <xpx>/(xrms*pxrms)]
      endif
c values provided by user:
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'xfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'yfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'tfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'xflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'yflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(5+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'tflush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(6+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'crossterm=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(7+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               CONTRACTENV
c=====================================================================
c contractenv (apply contraction map to rms envelopes):
763   continue
      return
c
c=====================================================================
c                               SETREFTRAJ
c=====================================================================
764   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na)) =0.d0  ! x=
      pmenu(2+mpp(na)) =0.d0  ! px=
      pmenu(3+mpp(na)) =0.d0  ! y=
      pmenu(4+mpp(na)) =0.d0  ! py=
      pmenu(5+mpp(na)) =0.d0  ! t=
      pmenu(6+mpp(na)) =0.d0  ! pt=
      pmenu(7+mpp(na)) =-9999.d0  ! s=
      pmenu(8+mpp(na)) =0.d0  ! sto=
      pmenu(9+mpp(na)) =0.d0  ! get=
      pmenu(10+mpp(na))=0.d0  ! write=
      cmenu(1+mppc(na))='false' !restart= true/false
!                                resets to values
!                                stored in loc 1, which contain values from
!                                start of run (unless overwritten by user)]
      cmenu(2+mppc(na))='true' !includearc=true/false [used w/ sto,get]
      endif
c values provided by user:
      call getparm(line,mmax,'x=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'px=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'y=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'py=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'t=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'pt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'arclen=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
c sto:
      call getparm(line,mmax,'sto=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      if(bufr(1).le.0 .or. bufr(1).ge.9)then
        if(idproc.eq.0)then
         write(6,*)'setreftraj error: sto must be in (1,...,8)'
        endif
        call myexit
      endif
c get:
      call getparm(line,mmax,'get=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
      if(bufr(1).le.0 .or. bufr(1).ge.9)then
        if(idproc.eq.0)then
         write(6,*)'setreftraj error: get must be in (1,...,8)'
        endif
        call myexit
      endif
c
      call getparm(line,mmax,'write=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'restart=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'includearc=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               SETARCLEN
c=====================================================================
765   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na)) =0.d0  ! s=
      pmenu(2+mpp(na))=0.d0  ! write=
      cmenu(1+mppc(na))='false' !restart= true/false
      endif
c values provided by user:
      call getparm(line,mmax,'s=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'write=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'restart=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                               WAKEDEFAULT
c=====================================================================
766   continue
      if(idproc.eq.0)then
       write(6,*)'CODE UNDER DEVELOPMENT. POSSIBLY PUT WAKEDEFAULT HERE'
      endif
      return
c
c=====================================================================
c                               EMITTANCE
c=====================================================================
c emittance:
767   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=0.d0  !assigned unit# for file
      pmenu(2+mpp(na))=5.d0  !precision
      pmenu(3+mpp(na))=1.d0  !nunits (=0 for dimensionless; =1 for physical)
      cmenu(1+mppc(na))='true'  !flag to compute 2-D emittances
      cmenu(2+mppc(na))='true'  !flag to compute 4-D transverse emittance
      cmenu(3+mppc(na))='false'  !flag to compute 6-D emittance
      cmenu(4+mppc(na))='keep'  !flag to keep/remove beam centroid
      cmenu(5+mppc(na))='emitrms.out'  !name of output file
      cmenu(6+mppc(na))='true'  !flag to flush output file
      endif
c values provided by user:
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'nunits=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'2d=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'4d=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'6d=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'centroid=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(5+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(6+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                           MATCHENV
c=====================================================================
768   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na)) =10.d0  ! iterations=
      pmenu(2+mpp(na)) =1.d-8  ! tolerance=
      cmenu(1+mppc(na))=' '    ! name=
      endif
      call getparm(line,mmax,'iterations=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'tolerance=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'name=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        cmenu(1+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                           FILEINFO
c       (controls printing of file open/close info)
c=====================================================================
769   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na)) =0.  ! info=
      endif
      call getparm(line,mmax,'info=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                               EGENGRAD
c=====================================================================
c egengrad:
770   continue
      if(initparms.eq.1)then
        pmenu(1+mpp(na))=0.d0  !assigned unit# for output file
        pmenu(2+mpp(na))=0.d0  !zstart
        pmenu(3+mpp(na))=0.d0  !zend
        pmenu(4+mpp(na))=1.d0  !nz (intervals)
        pmenu(5+mpp(na))=0.d0  !frequency
        pmenu(6+mpp(na))=0.d0  !radius
        pmenu(7+mpp(na))=3.d2  !kmax
        pmenu(8+mpp(na))=3.d3  !nk (intervals)
        pmenu(9+mpp(na))=0.d0  !infiles
        pmenu(10+mpp(na))=5.d0 !precision
        cmenu(1+mppc(na))='rfdata'  ! name of input datafile
        cmenu(2+mppc(na))='crz.dat' ! name of output datafile
        cmenu(3+mppc(na))='true'    ! flag to flush output file(s)
        cmenu(4+mppc(na))='e0.dat'  ! name of optional output datafile
        cmenu(5+mppc(na))='false'   ! flag to (not) write char. fns.
        cmenu(6+mppc(na))='t7'      ! E-field file type
        cmenu(7+mppc(na))=' '       ! name of optional E-field diagnos.
      endif
c values provided by user:
cp2
      call getparm(line,mmax,'zstart=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
cp3
      call getparm(line,mmax,'zend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
cp4
      call getparm(line,mmax,'nz=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
cp5
      call getparm(line,mmax,'frequency=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
cp6
      call getparm(line,mmax,'radius=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
cp7
      call getparm(line,mmax,'kmax=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
cp8
      call getparm(line,mmax,'nk=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
cp9
      call getparm(line,mmax,'infiles=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
cp10
      call getparm(line,mmax,'precision=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
cc1
      call getparm(line,mmax,'efile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
cc2
      call getparm(line,mmax,'crzfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
cc3
      call getparm(line,mmax,'flush=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
cc4
      call getparm(line,mmax,'charfile=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(4+mppc(na))=cbuf
      endif
cc5
      call getparm(line,mmax,'wrtchar=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(5+mppc(na))=cbuf
      endif
cc6
      call getparm(line,mmax,'ftype=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(6+mppc(na))=cbuf
      endif
cc7
      call getparm(line,mmax,'ediag=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(7+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                WRTMAP: Write Map to File (KMP: 8 Nov 2006)
c=====================================================================
c
771   continue
c
c Initialized Parameters to default values (in case they are not specified):
      if(initparms.eq.1)then
        cmenu(1+mppc(na))='map.dat'     ! name of file for map write
        cmenu(2+mppc(na))='lastslice'   ! kind of write: accumulated, lastslice
        cmenu(3+mppc(na))='append'      ! iostatus of write to file: overwrite, append
      endif
c
c Get file name for map write:
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
c
c Get kind of write to perform:
      call getparm(line,mmax,'kind=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
c
c Get iostatus of write:
      call getparm(line,mmax,'iostat=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(3+mppc(na))=cbuf
      endif
      return
c
c=====================================================================
c                RDMAP: Read a Map from File (KMP: 14 Dec 2006)
c=====================================================================
c
772   continue
c
c Initialized Parameters to default values (in case they are not specified):
      if(initparms.eq.1)then
        cmenu(1+mppc(na))='map.dat'     ! name of file for map read
        cmenu(2+mppc(na))='true'        ! name of file for map read
      endif
c
c Get file name for map write:
      call getparm(line,mmax,'file=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
c
c Get rewind information
      call getparm(line,mmax,'rewind=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      return

c=====================================================================
c                           SPARES TYPE CODES
c=====================================================================

773   continue
774   continue
775   continue
      return
c
c
c     8: advanced commands ********************************
c
c           'cod     ','amap    ','dia     ','dnor    ','exp     ',
18    go to (801,       802,       803,       804,       805,
c           'pdnf    ','psnf    ','radm    ','rasm    ','sia     ',
     &       806,       807,       808,       809,       810,
c           'snor    ','tadm    ','tasm    ','tbas    ','gbuf    ',
     &       811,       812,       813,       814,       815,
c           'trsa    ','trda    ','smul    ','padd    ','pmul    ',
     &       816,       817,       818,       819,       820,
c           'pb      ','pold    ','pval    ','fasm    ','fadm    ',
     &       821,       822,       823,       824,       825,
c           'sq      ','wsq     ','ctr     ','asni    ','pnlp    ',
     &       826,       827,       828,       829,       830,
c           'csym    ','psp    ','mn       ','bgen    ','tic     ',
     &       831,       832,      833,        834,       835,
c           'ppa     ','moma   ','geom     ','fwa     '/
     &       836,       837,      838,        839),kt2
c
c     off-momentum closed orbit analysis
c
801   continue
      return
c
c     apply map to a function or moments
c
802   continue
      return
c
c     dynamic invariant analysis
c
803   continue
      return
c
c     dynamic normal form analysis
c
804   continue
      return
c
c     compute exponential
c
805   continue
      return
c
c     compute power of dynamic normal form
c
806   continue
      return
c
c     compute power of static normal form
c
807   continue
      return
c
c     resonance analyze dynamic map
c
808   continue
      return
c
c     resonance analyze static map
c
809   continue
      return
c
c     static invariant ayalysis
c
810   continue
      return
c
c     static normal form analysis
c
811   continue
      return
c
c     twiss analyze dynamic map
c
812   continue
      return
c
c 'tasm': twiss analyze static map
813   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=2
      pmenu(2+mpp(na))=1.d-3
      pmenu(3+mpp(na))=1
      pmenu(4+mpp(na))=0
      pmenu(5+mpp(na))=3
      pmenu(6+mpp(na))=0
      endif
c values provided by user:
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'delta=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'idata=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'ipmaps=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'iwmaps=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      return
c
c     translate basis
c
814   continue
      return
c
c     get buffer contents
c
815   continue
      return
c
c     transport static (script) A
c
816   continue
      return
c
c    transport dynamic script A
c
817   continue
      return
c
c     multiply polynomial by a scalar
c
818   continue
      return
c
c     add two polynomials
c
819   continue
      return
c
c     multiply two polynomials
c
820   continue
      return
c
c     Poisson bracket two polynomials
c
821   continue
c defaults:
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.  !map1in
      pmenu(2+mpp(na))=2.  !map2in
      pmenu(3+mpp(na))=0.  !mapout
      endif
c values provided by user:
      call getparm(line,mmax,'map1in=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'map2in=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'mapout=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      return
c
c     polar decompose matrix portion of transfer map
c
822   continue
      return
c
c     evaluate a polynomial
c
823   continue
      return
c
c     fourier analyze static map
c
824   continue
      return
c
c     fourier analyze dynamic map
c
825   continue
      return
c
c     select quantities
c
826   continue
      return
c
c     write selected quantities
c
827   continue
      return
c
c     change tune ranges
c
828   continue
      return
c
c     apply script N inverse
c
829   continue
      return
c
c     compute power of nonlinear part
c
830   continue
      return
c
c     check for symplecticity
c
831   continue
      return
c
c     (psp) compute scalar product of two polynomials
c
832   continue
      return
c
c     (mn) compute matrix norm
c
833   continue
      return
c
c=====================================================================
c                               BGEN
c=====================================================================
c (bgen) generate a beam
834   continue
c     write(6,*)'HERE I AM AT BGEN; MMAX=',mmax, '   LINE(1:800)='
c     write(6,*)line(1:800)
      write(6,*)'getting bgen parameters in sif.f'
c defaults:
      if(initparms.eq.1)then
      job=6
      iopt=1
      nray=10000
      iseed=1234567
      isend=3
      ipset=0
      xx=1.d-3
      yy=1.d-3
      tt=1.d-3
      sigmax=5.d0
      unit1=0
      unit2=0
      pmenu(1+mpp(na))=job    !this is called 'dist' by the parser
      pmenu(2+mpp(na))=iopt
      pmenu(3+mpp(na))=nray   !this is called 'maxray' by the parser
      pmenu(4+mpp(na))=iseed
      pmenu(5+mpp(na))=isend
      pmenu(6+mpp(na))=ipset
      pmenu(7+mpp(na))=xx
      pmenu(8+mpp(na))=yy
      pmenu(9+mpp(na))=tt
      pmenu(10+mpp(na))=sigmax
      pmenu(11+mpp(na))=unit1
      pmenu(12+mpp(na))=unit2
      cmenu(1+mppc(na))=' '
      cmenu(2+mppc(na))=' '
      endif
c values provided by user:
      call getparm(line,mmax,'dist=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'maxray=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'iseed=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'ipset=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      call getparm(line,mmax,'xx=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(7+mpp(na))=bufr(1)
      call getparm(line,mmax,'yy=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(8+mpp(na))=bufr(1)
      call getparm(line,mmax,'tt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(9+mpp(na))=bufr(1)
      call getparm(line,mmax,'sigmax=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(10+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'unit1=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(11+mpp(na))=bufr(1)
      call getparm(line,mmax,'unit2=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(12+mpp(na))=bufr(1)
c
      call getparm(line,mmax,'file1=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(1+mppc(na))=cbuf
      endif
      call getparm(line,mmax,'file2=',bufr,keypres,numpres,1,cbuf)
      if(keypres.and.numpres)then
        call quotechk(cbuf,16)
        cmenu(2+mppc(na))=cbuf
      endif
      return
c
c     (tic) translate (move) initial conditions
c
835   continue
      return
c
c     (ppa) principal planes analysis
c
836   continue
      return
c
c     (moma) moment and map analysis
c
837   continue
      return
c
c     (geom) compute geometry of a loop
c
838   continue
      return
c
c     (fwa) copy file to working array
c
839   continue
      return
c
c     9: procedures and fitting and optimization *************************
c
c           'bip     ','bop     ','tip     ','top     ',
19    go to (901,       902,       903,       904,
c           'aim     ','vary    ','fit     ','opt     ',
     &       905,       906,       907,       908,
c           'con1    ','con2    ','con3    ','con4    ','con5    ',
     &       909,       910,       911,       912,       913,
c           'mrt0    ',
     &       914,
c           'mrt1    ','mrt2    ','mrt3    ','mrt4    ','mrt5    ',
     &       915,       916,       917,       918,       919,
c           'fps     ',
     &       920,
c           'cps1    ','cps2    ','cps3    ','cps4    ','cps5    ',
     &       921,       922,       923,       924,       925,
c           'cps6    ','cps7    ','cps8    ','cps9    ',
     &       926,       927,       928,       929,
c           'dapt    ','grad    ','rset    ','flag    ','scan    ',
     &       930,       931,       932,       933,       934,
c           'mss     ','spare1  ','spare2  '/
     &       935,       936,       937),kt2
c     begin procedures
c
c=====================================================================
c                                BIP
c=====================================================================
c 'bip    ': begin inner procedure
901   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=10.    !ntimes
      endif
      call getparm(line,mmax,'ntimes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                BOP
c=====================================================================
c 'bop    ': begin outer procedure
902   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=10.    !ntimes
      endif
      call getparm(line,mmax,'ntimes=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                TIP
c=====================================================================
c 'tip    ': terminate inner procedure
903   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.    !iopt
      endif
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                TOP
c=====================================================================
c 'tip    ': terminate outer procedure
904   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.    !iopt
      endif
      call getparm(line,mmax,'iopt=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      return
c
c     end procedures
c
c
c     specify aims
c
c=====================================================================
c                                AIM
c=====================================================================
c 'aim    ': specify quantities to be fit/optimized and set target values
905   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=2.    !job
      pmenu(2+mpp(na))=0.    !infile
      pmenu(3+mpp(na))=0.    !logfile
      pmenu(4+mpp(na))=0.    !iquiet
      pmenu(5+mpp(na))=0.    !istore
      pmenu(6+mpp(na))=1.    !isend
      endif
      call getparm(line,mmax,'job=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'infile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'logfile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'iquiet=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'istore=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      return
c
c     specify quantities to be varied
c
c=====================================================================
c                                VARY
c=====================================================================
c 'vary   ': specify quantities to be varied
906   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=1.    !job
      pmenu(2+mpp(na))=0.    !infile
      pmenu(3+mpp(na))=0.    !logfile
      pmenu(4+mpp(na))=0.    !isb (scaling and bounds)
      pmenu(5+mpp(na))=0.    !istore
      pmenu(6+mpp(na))=1.    !isend
      endif
      call getparm(line,mmax,'job=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'infile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'logfile=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'isb=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'istore=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      return
c
c=====================================================================
c                                FIT
c=====================================================================
c 'fit    ': fit to achieve aims
907   continue
      if(initparms.eq.1)then
      pmenu(1+mpp(na))=10.     !job
      pmenu(2+mpp(na))=1.      !aux
      pmenu(3+mpp(na))=1.d-8   !error
      pmenu(4+mpp(na))=1.d-3   !delta
      pmenu(5+mpp(na))=0.    !mprint
      pmenu(6+mpp(na))=1.    !isend
      endif
      call getparm(line,mmax,'job=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(1+mpp(na))=bufr(1)
      call getparm(line,mmax,'aux=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(2+mpp(na))=bufr(1)
      call getparm(line,mmax,'error=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(3+mpp(na))=bufr(1)
      call getparm(line,mmax,'delta=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(4+mpp(na))=bufr(1)
      call getparm(line,mmax,'mprint=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(5+mpp(na))=bufr(1)
      call getparm(line,mmax,'isend=',bufr,keypres,numpres,0,cbuf)
      if(keypres.and.numpres)pmenu(6+mpp(na))=bufr(1)
      return
c
c
c     optimize
c
908   continue
      return
c
c     constraints
c
909   continue
      return
910   continue
      return
911   continue
      return
912   continue
      return
913   continue
      return
c
c     merit functions
c
c     least squares merit function
c
914   continue
      return
c
c     user supplied merit functions
c
915   continue
      return
916   continue
      return
917   continue
      return
918   continue
      return
919   continue
      return
c
c     free parameter sets
c
920   continue
      return
c
c     capture parameter sets
c
921   continue
      return
922   continue
      return
923   continue
      return
924   continue
      return
925   continue
      return
926   continue
      return
927   continue
      return
928   continue
      return
929   continue
      return
c
c     compute dynamic aperture (dapt)
c
930   continue
      return
c
c     gradient (grad)
c
931   continue
      return
c
c     rset
c
932   continue
      return
c
c     flag
c
933   continue
      return
c
c     scan
c
934   continue
      return
c
c     mss
c
935   continue
      return
c
c     spare1

936   continue
      return
c
c     spare2

937   continue
      return
      return
      end
c
      subroutine getparm(line,mmax,kywrd,bufr,keypres,numpres,iopt,cbuf)
c "string" has to be long enough to hold a single number or expression
      use parallel, only : idproc
      use acceldata
      include 'impli.inc'
      include 'files.inc'
cryne 7/15/2002      parameter (nstrmax=80)
c nstrmax=80 should be sufficient for most purposes. The only reason
c to make it longer is if there is an arithmetic (symbolic) expression
c that extends over several lines. So, just in case, set nstrmax=800
      parameter (nstrmax=800)
      character (len=*) line
      character(len=nstrmax) string
      character(len=*) kywrd
      character(len=*) cbuf
      logical keypres,numpres
      dimension bufr(1)
      character*16 symb(50)
      integer istrtsym(50)
      common/showme/iverbose

      keylen=len(kywrd)
      if(kywrd(1:keylen).eq.'=')then
        n1=index(line,kywrd)
        goto 500
      endif
c--------------------------------
      if(iverbose.eq.2 .and. idproc.eq.0)then
      write(6,*)'[getparm] kywrd=',kywrd(1:keylen),', line(1:',mmax,')='
      write(6,*)line(1:mmax)
      endif
cccc  n1=index(line,kywrd)
      call getsymb(line,mmax,symb,istrtsym,nsymb)
      if(iverbose.eq.2 .and. idproc.eq.0)then
        write(6,*)'(getparm)returned from getsymb with nsymb=',nsymb
        do n=1,nsymb
        write(6,*)'n,symb(n)=',n,symb(n)
        enddo
      endif
      n1=0
      if(nsymb.gt.0)then
        do n=1,nsymb
c         write(6,*)n,' ',kywrd(1:keylen),' ',trim(symb(n))
          if(kywrd(1:keylen).eq.trim(symb(n))//'=')then
          n1=istrtsym(n)
c         write(6,*)'found a match with n, n1 = ',n,n1
          exit
          endif
        enddo
      endif
c--------------------------------
  500 continue
      if(n1.eq.0)then
        keypres=.false.
        numpres=.false.
c       if(jwarn.eq.1)then
c        write(6,*)'warning: could not find ',kywrd,' on input line:'
c        write(6,*)line(1:mmax)
c       endif
        return
      endif
c possibly this came from a multiline input.
c in that case, there may be '&' that need to be removed.
      do n=1,mmax
        if(line(n:n).eq.'&')line(n:n)=' '
      enddo
c obtain value for kywrd:
c     write(6,*)'obtaining value for keyword:',kywrd
c     write(6,*)'with line equal to'
c     write(6,*)line(1:mmax)
      keypres=.true.
      string(1:nstrmax)=' '
      m=1
c     write(6,*)'starting do loop with n=',n1+keylen,'  to',mmax
      do n=n1+keylen,mmax
cryne 7/9/2002          if(line(n:n).eq.' ')exit
        if(line(n:n).eq.',')exit
        if((line(n:n).eq.'!').or.(line(n:n).eq.';'))exit
        string(m:m)=line(n:n)
cryne can exit after storing another '=' (if any)
        if(line(n:n).eq.'=')exit
        m=m+1
        if(m.gt.nstrmax)then
          write(6,*)'TROUBLE: m > nstrmax'
        endif
      enddo
c     write(6,*)'finished do loop. now string='
c     write(6,*)string(1:nstrmax)
cccccccccccccccccccccc
c in case the user omitted the commas:
c find the next occurence of '=':
      mm1=index(string,'=')
c     write(6,*)'mm1=',mm1
      if(mm1.eq.0)goto 150
      string(mm1:nstrmax)=' '
      do n=mm1-1,1,-1
        if((string(n:n).eq.' ').or.(string(n:n).eq.','))exit
        string(n:n)=' '
      enddo
cccccccccccccccccccccc
  150 continue
c     write(6,*)'*string* = '
c     write(6,*)string(1:nstrmax)
c option to return string instead of value (i.e. if looking for a string):
c use numpres to indicate that this is normal return
      if(iopt.eq.1)then
        nbuf=len(cbuf)
        cbuf(1:nbuf)=string(1:nbuf)
        numpres=.true.
c       write(6,*)'returning from getparm w/ iopt=1'
c       write(6,*)'kywrd,cbuf=',kywrd(1:keylen),cbuf(1:nbuf)
        return
      endif
c
      numpres=.false.
c------------------------------------
c COMMENT THIS SECTION OUT WHEN FPARSER IS WORKING ON YOUR MAC:
c     call txtnum(string,nstrmax,1,nget,bufr)
c     read(string,*,iostat=numerr)value
c     if(numerr.eq.0)then
c       numpres=.true.
c       bufr(1)=value
c     endif
c-------------------------------------
      if(.not.numpres)then
        call strngfeval(string,nstrmax,value,numpres)
        if(numpres)bufr(1)=value
c       do j=1,nconst
c         if(string.eq.constr(j))then
c           write(6,*)'j,constr,conval=',j,constr(j),conval(j)
c           bufr(1)=conval(j)
c           numpres=.true.
c           exit
c         endif
c       enddo
        if(.not.numpres)write(6,*)'ERROR READING ELEMENT PARAMTERS:',   &
     &  'The string ',trim(string), ' has not been defined'
      endif
      if(numpres)then
c        write(6,*)'bufr(1)=',bufr(1)
      else
         write(6,*)'number not found for kywrd ',kywrd
      endif
      return
      end
c
      subroutine strngfeval(string,nstrmax,result,numpres)
      use parallel, only : idproc
      use acceldata
      include 'impli.inc'
c     parameter (nstrmax=80)
      character(len=nstrmax) string
      character*16 symb(50),newsymb
      integer istrtsym(50)
      logical numpres
      dimension val(50),bufr(1)
      common/showme/iverbose
      common/symbdef/isymbdef
      numpres=.false.
c     write(6,*)'evaluating string:',string
c     iverbose=2
      if(iverbose.eq.2 .and. idproc.eq.0)then
        write(6,*)'inside strngfeval; nstrmax=',nstrmax,' ;string='
        write(6,*)string(1:nstrmax)
      endif
      call getsymb(string,nstrmax,symb,istrtsym,nsymb)
      if(iverbose.eq.2 .and. idproc.eq.0)then
        write(6,*)'(strngfeval)returned from getsymb; nsymb=',nsymb
        do n=1,nsymb
        write(6,*)'n,symb(n)=',n,symb(n)
        enddo
      endif
c this is wrong; there can be an expression even though there are no symbols.
c so evaluate the symbolic expression regardless.
c     if(nsymb.eq.0)then
c       call txtnum(string,80,1,nget,bufr)
c       if(nget.ge.1)then
c         result=bufr(1)
c         numpres=.true.
c       else
c         write(6,*)'trouble(strngfeval):txtnum did not return a number'
c       endif
c       return
c     endif
c found symbols on this line; evaluate the symbolic expression:
      do n=1,nsymb
c this is the place to check for function names (e.g. sqrt) and skip:
        if(iverbose.eq.2 .and. idproc.eq.0)then
        write(6,*)'[strngfeval] n,symb(n),nconst =',n,symb(n),nconst
        endif
        if(trim(symb(n)).eq.'sqrt')cycle
        if(trim(symb(n)).eq.'sin')cycle
        if(trim(symb(n)).eq.'cos')cycle
        if(trim(symb(n)).eq.'tan')cycle
        if(trim(symb(n)).eq.'asin')cycle
        if(trim(symb(n)).eq.'acos')cycle
        if(trim(symb(n)).eq.'atan')cycle
        if(trim(symb(n)).eq.'abs')cycle
        do j=1,nconst
c         write(6,*)'j,constr(j)=',j,constr(j)
          if(symb(n).eq.constr(j))then
ccccc       write(6,*)'j,constr,conval=',j,constr(j),conval(j)
            val(n)=conval(j)
            goto 199
          endif
        enddo
c the symbol was not found in the constr array; check for [...]
c new code-------------------------
c if you find it, set val(n)=... and goto 199
      ir1=index(symb(n),'[')
      ir2=index(symb(n),']')
      if(ir1.eq.0 .or. ir2.eq.0)goto 198
c found [...] in the symbol. Now determine the value of xxxx[...]:
        write(6,*)'reached the [...] section of code'
        newsymb=symb(n)(ir1+1:ir2-1)
        nlength=ir2-ir1-1
c       write(6,*)'calling str2int w/ newsymb,nlength=',newsymb,nlength
        call str2int(newsymb,nlength,mval)
        write(6,*)'returned from str2int w/ mval=',mval
c now that the value of the symbol has been found, make sure that
c it corresponds to a menu item:
        newsymb=symb(n)(1:ir1-1)
        write(6,*)'looking up the following symbol: ',newsymb
        call lookup(newsymb,itype,indx)
        if(itype.ne.1)then
          write(6,*)'error: the symbol ',newsymb,' is not in the menu'
          call myexit
        endif
        write(6,*)'Found ',newsymb,' in the menu (element #',indx,')'
        write(6,*)'Parameter ',mval,' has the value ',                     &
     &            pmenu(mval+mpp(indx))
        val(n)=pmenu(mval+mpp(indx))
        goto 199
c ---------------------------------
  198   continue
cryne 9/26/02:
c the symbol has not been defined. deal with it:
        if(isymbdef.eq.1)then
         if(idproc.eq.0)then
         write(6,*)'error: symbol ',symb(n),' is used in expression:'
         write(6,*)string
         write(6,*)'but the symbol has not been defined. Setting to 0.'
         endif
         val(n)=0.d0
        else
         if(idproc.eq.0)then
         write(6,*)'error: symbol ',symb(n),' is used in expression:'
         write(6,*)string
         write(6,*)'but the symbol has not been defined. Halting.'
         endif
         call myexit
        endif
  199 continue
      enddo
      nfunc=1
!dbs -- fproutine expects an array so construct one on the fly
!dbs  call fproutine(string,nfunc,symb,val,nsymb,result,nget)
      call fproutine( (/ string /),nfunc,symb,val,nsymb,result,nget)
      if(nget.eq.1)numpres=.true.
      return
      end
c
      subroutine str2int(newsymb,nlength,mval)
      character*16 newsymb
      integer nlength,mval,n,ntmp
      if(nlength.le.0)then
        write(6,*)'error (str2int): nlength = ',nlength
        call myexit
      endif
      mval=0
      do n=nlength,1,-1
        if(newsymb(n:n).eq.'0')ntmp=0
        if(newsymb(n:n).eq.'1')ntmp=1
        if(newsymb(n:n).eq.'2')ntmp=2
        if(newsymb(n:n).eq.'3')ntmp=3
        if(newsymb(n:n).eq.'4')ntmp=4
        if(newsymb(n:n).eq.'5')ntmp=5
        if(newsymb(n:n).eq.'6')ntmp=6
        if(newsymb(n:n).eq.'7')ntmp=7
        if(newsymb(n:n).eq.'8')ntmp=8
        if(newsymb(n:n).eq.'9')ntmp=9
        npow=nlength-n
c       write(6,*)'n,ntmp,npow=',n,ntmp,npow
        mval=mval+ntmp*10**npow
      enddo
      return
      end
c
      subroutine getsymb(line,nstrmax,symb,istrtsym,nsymb)
cryne 09/21/02 modified to include ] and [
      use acceldata
      include 'impli.inc'
      character (len=*) line
      character*16 symb(*)
      integer istrtsym(*)
      character*1 str
      character*26 alpha
ccc   character*38 alphaug
      character*40 alphaug
      character*11 numerdot
      logical acheck,echeck,dcheck
      alpha='abcdefghijklmnopqrstuvwxyz'
ccc   alphaug='abcdefghijklmnopqrstuvwxyz0123456789_.'
      alphaug='abcdefghijklmnopqrstuvwxyz0123456789_.[]'
      numerdot='0123456789.'
      nsymb=0
      i=0
! search for a new symbol beginning with an alpha character:
  100 i=i+1
c     if(i.eq.81)goto 9999
      if(i.ge.nstrmax+1)goto 9999
      istart=i
c     if(line(i:i) .eq. an alpha)then
      str=line(i:i)
      echeck=.false.
      dcheck=.false.
      if(i.ge.2)then
        echeck=(str.eq.'e').and.(scan(line(i-1:i-1),numerdot).ne.0)
        dcheck=(str.eq.'d').and.(scan(line(i-1:i-1),numerdot).ne.0)
      endif
      acheck=.false.
      if(scan(str,alpha).ne.0)acheck=.true.
      if( (acheck) .and. (.not.echeck) .and. (.not.dcheck) )then
        nsymb=nsymb+1
  200   i=i+1
c       if(i.eq.81)goto 9999
        if(i.eq.nstrmax+1)goto 9999
cryne 5/4/2006     if( (line(i:i).eq.'!') .or. (line(i:i).eq.';') )goto 9999
        if( (line(i:i).eq.'!') .or. (line(i:i).eq.';') )goto 201
c       if(line(i:i).eq. an alphanumeric or _ or . or "'")goto 200
        str=line(i:i)
!dbs -- if the last char in the line is the last char in the symbol
!dbs    then process the symbol
        if( scan(str,alphaug).ne.0 )then
          if( i.lt.nstrmax ) goto 200
          !make it seem like the next char is a symbol terminator
          i = i + 1
        endif
!dbs        if(scan(str,alphaug).ne.0)goto 200
cryne 5/4/2006 ------
  201   continue
c--------------------
cryne 11/04/02
c store the symbol but first make sure it is not too long:
        if(i-istart.gt.16)then
         write(6,*)'(getsymb)warning: long character string on line='
         write(6,*)line
        endif
c
        imax=i-1
        if(i-istart.gt.16)imax=istart+16-1
c       symb(nsymb)=line(istart:i-1)
        symb(nsymb)=line(istart:imax)
        istrtsym(nsymb)=istart
        goto 100
      endif
      goto 100
 9999 continue
      return
      end
c
      subroutine getsymb60(line,nstrmax,symb,istrtsym,nsymb)
cryne 09/21/02 modified to include ] and [
      use acceldata
      include 'impli.inc'
      character (len=*) line
      character*60 symb(*)
      integer istrtsym(*)
      character*1 str
      character*26 alpha
ccc   character*38 alphaug
      character*40 alphaug
      character*11 numerdot
      logical acheck,echeck,dcheck
      alpha='abcdefghijklmnopqrstuvwxyz'
ccc   alphaug='abcdefghijklmnopqrstuvwxyz0123456789_.'
      alphaug='abcdefghijklmnopqrstuvwxyz0123456789_.[]'
      numerdot='0123456789.'
      nsymb=0
      i=0
! search for a new symbol beginning with an alpha character:
  100 i=i+1
c     if(i.eq.81)goto 9999
      if(i.eq.nstrmax+1)goto 9999
      istart=i
c     if(line(i:i) .eq. an alpha)then
      str=line(i:i)
      echeck=.false.
      dcheck=.false.
      if(i.ge.2)then
        echeck=(str.eq.'e').and.(scan(line(i-1:i-1),numerdot).ne.0)
        dcheck=(str.eq.'d').and.(scan(line(i-1:i-1),numerdot).ne.0)
      endif
      acheck=.false.
      if(scan(str,alpha).ne.0)acheck=.true.
      if( (acheck) .and. (.not.echeck) .and. (.not.dcheck) )then
        nsymb=nsymb+1
  200   i=i+1
c       if(i.eq.81)goto 9999
        if(i.eq.nstrmax+1)goto 9999
        if( (line(i:i).eq.'!') .or. (line(i:i).eq.';') )goto 9999
c       if(line(i:i).eq. an alphanumeric or _ or . or "'")goto 200
        str=line(i:i)
        if(scan(str,alphaug).ne.0)goto 200
cryne 11/04/02
c store the symbol but first make sure it is not too long:
        if(i-istart.gt.60)then
         write(6,*)'(getsymb)warning: long character string on line='
         write(6,*)line
        endif
c
        imax=i-1
        if(i-istart.gt.60)imax=istart+16-1
c       symb(nsymb)=line(istart:i-1)
        symb(nsymb)=line(istart:imax)
        istrtsym(nsymb)=istart
        goto 100
      endif
      goto 100
 9999 continue
      return
      end
c
c
c
      subroutine fproutine(func,nfunc,var,val,nvar,res,nget)
      USE parameters, ONLY: rn
      USE fparser,   ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
      IMPLICIT NONE
      INTEGER nfunc,nvar,i,nget
      CHARACTER (LEN=*), DIMENSION(nfunc) :: func
      CHARACTER (LEN=*), DIMENSION(nvar)  :: var
      REAL(rn),          DIMENSION(nvar)  :: val
      REAL(rn)                            :: res
!
      CALL initf (nfunc)     ! Initialize function parser for nfunc functions
      DO i=1,nfunc
       CALL parsef(i,func(i),var) ! Parse and bytecompile ith function string
      END DO
!     WRITE(*,*)'==> Bytecode evaluation:'
      DO i=1,nfunc
        res=evalf(i,val) ! Interprete bytecode representation of ith function
        IF(EvalErrType.gt.0)WRITE(6,*)'*** Error: ',EvalErrMsg ()
        nget=0
        IF(EvalErrType.eq.0)nget=1
!       WRITE(6,*)trim(func(i)),'=',res
      END DO
      return
      end
c
      subroutine quotechk(fname,lmax)
c check for unnecessary quotes and apostrophes
      character (len=*) fname
c
      if(fname.ne.' ')then
        if((fname(1:1).eq.'"').or.(fname(1:1).eq."'"))then
         write(6,*)'(pmif) quotes and apostrophes not needed in: ',fname
         write(6,*)'They will be stripped off.'
c        lmax=len(fname)
         do i=1,lmax-1
           fname(i:i)=fname(i+1:i+1)
         enddo
         fname(lmax:lmax)=' '
         n=index(fname,'"')
         if(n.ne.0)fname(n:n)=' '
         n=index(fname,"'")
         if(n.ne.0)fname(n:n)=' '
        endif
      endif
      return
      end
c---------------------------------------------------
