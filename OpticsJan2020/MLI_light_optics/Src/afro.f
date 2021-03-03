************************************************************************
*  header                      AFRONT                                  *
*   Front end of MARYLIE                                               *
************************************************************************
*
************************************************************************
*      MARYLIE and all its subroutines are copyrighted (1987) by       *
*                        Alex J. Dragt.                                *
*      All rights to MARYLIE and its subroutines are reserved.         *
************************************************************************
*
*  Version Date: 5/25/98
*
***********************************************************************
*  All of MARYLIE 3.0 is believed to conform to Fortran 77 standards.
*  In addition, all of MARYLIE 3.0 is self contained save for
*
*     1) the time functions called in the MARYLIE
*        subroutine cputim
*
*        and
*
*     2) the 'exit' subroutine called (and only called) in the
*        MARYLIE subroutine myexit.
*
*  The two subroutines cputim and myexit and the subroutine mytime
* (which is the only one that calls cputim) are kept separately
*  under the heading "XTRA".
*********************************************************************
c
*********************************************************************
c
*  Main program for MARYLIE 3.0
*  Written originally by Rob Ryne ca 1984 and revised by
*  Petra Schuett November 1987
*
      program mainmary
c
      use rays
      use acceldata
      use lieaparam, only : monoms
      use ml_timer
      include 'impli.inc'
      include 'files.inc'
      include 'time.inc'
      include 'ind.inc'
      include 'codes.inc'
c
      dimension p(6)
      common/xtrajunk/jticks1,iprinttimers
c
      call init_parallel
      call init_ml_timers
      call system_clock(count=jticks1)
c
c
c initiate starting time
c nint(p(1)).ne.0 resets a variable called tstart in common/mltime/
c which is contained in the file time.inc     rdr 08/25/2001
      p(1)=1.d0
      p(2)=12.d0
      p(3)=0.d0
      call mytime(p)
c
c write out copyright message at terminal
      call cpyrt
c
c initialize commons (other than those in block data's)
c initialize lie algebraic things
c
      call new_acceldata
!     call new_particledata  !do this for now for ML30 compat. RDR 8/10/02
c
      call setup
c
c read master input file
c
      call dumpin
c
c========================================================================
c print debug info on file 79
      write(79,*)'nconst=',nconst
      write(79,*)'   n    constr      conval'
      do n=1,nconst
      write(79,*)n,constr(n),conval(n)
c     call flush(79)
      enddo
c
      write(79,*)' '
      write(79,*)'na=',na
      write(79,*)'listing of n,nt1(n),nt2(n) for n=1,...,na'
      do n=1,na
      k1=nt1(n)
      k2=nt2(n)
      call lookup(lmnlbl(n),ntype,ith)
      nn0=mpp(ith)+1
      nn1=mppc(ith)+1
      write(79,1236)n,                                                     &
     & nt1(n),nt2(n),ltc(k1,k2),lmnlbl(n),nn0,nn1,nrp(k1,k2),ncp(k1,k2)
c     call flush(79)
 1236 format(3(i4,1x),2(a16,1x),2(i5,1x),6x,2(i4,1x))
      enddo
c
      write(79,*)' '
      write(79,*)'mpprpoi=',mpprpoi
      write(79,*)'listing of n,pmenu(n) for n=1,...,mpprpoi'
      do n=1,mpprpoi
c          the following loop makes this whole thing very inefficent,
c          but I am going to do it anyway in order to make this info
c          more readable:
!          do mmm=1,na
!          k1=nt1(mmm)
!          k2=nt2(mmm)
!          call lookup(lmnlbl(mmm),ntype,ith)
!          nn0=mpp(ith)+1
!          if(nn0.eq.n.and.nrp(k1,k2).ne.0)then
!            write(79,*)lmnlbl(mmm),ltc(k1,k2)
!          endif
!          enddo
        write(79,*)n,pmenu(n)
c       call flush(79)
      enddo
      write(79,*)' '
      write(79,*)'mppcpoi=',mppcpoi
      write(79,*)'listing of n,cmenu(n) for n=1,...,mppcpoi'
      do n=1,mppcpoi
        write(79,*)n,cmenu(n)
c       call flush(79)
      enddo
c========================================================================

c
c debug output
c      call dump(jodf)
c      call dump(jof)
c
c write out status at terminal
      if(idproc.eq.0)then
      write(jof,100)
  100 format(1h ,'Data input complete; going into #labor.')
      endif
c
c begin actual calculations
      call tran
c
c the end of the program is reached in MYEXIT
cryne ...but call myexit here in case the user forgot to use an "end" command
      call myexit
      end
c
***********************************************************************
      block data misc
c-----------------------------------------------------------------------
c initialize miscellaneous common variables
c Written by Petra Schuett, November 1987 based on earlier work of
c Rob Ryne and Liam Healy
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c--------------------
      include 'core.inc'
      include 'files.inc'
      include 'lims.inc'
      include 'maxcat.inc'
      include 'infin.inc'
      include 'zeroes.inc'
c--------------------
      data ordcat/6/,topcat/monoms/
      data ordlib/6/,toplib/monoms/
c--------------------
      data bottom/0,1, 7,28, 84,210,462, 924,1716,3003,5005, 8008,12376/
      data top   /0,6,27,83,209,461,923,1715,3002,5004,8007,12375,18563/
c--------------------
      data lf,icf,jfcf,mpi,mpo,jif,jof,jodf,ibrief,iquiet/              &
     &     11, 13,  14, 15, 16,  5,  6,  12,     2,     0/
c--------------------
      data inuse/maxlum*0/
c--------------------
      data xinf,     yinf,    tinf,    ginf/                            &
     &     1000.,    1000.,   1000.,   1000./
c--------------------
      data fzer,     detz/                                              &
     &     0.d0,     0.d0/
c
c--------------------
c
cryne 8/15/2001 initialize parameters for poisson calculation:
      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0
      data nxsave,nysave,nzsave,noresize,nadj0/0,0,0,0,0/
      common/poiblock1/nspchset,nsckick
      data nspchset,nsckick/0,-1/
      common/autotrk/lautotrk,ntrktype,ntrkorder
      data lautotrk,ntrktype,ntrkorder/0,2,5/
cryne 8/15/2001 parameters for automatic application of commands
      character*16 autostr
      common/autopnt/autostr(3)
      common/autolog/lautoap1,lautoap2,lautoap3,lrestrictauto
      common/autocon/lautocon
      data autostr                                                      &
     &/'xxxxxxxxxxxxxxxx','yyyyyyyyyyyyyyyy','zzzzzzzzzzzzzzzz'/
      data lautoap1,lautoap2,lautoap3/0,0,0/
      data lautocon/1/
c
      common/bfileinfo/ifileinfo
      data ifileinfo/0/
c
cryne 12 Nov 2003
      character*16 cparams
      common/scparams/rparams(60),cparams(60)
      data rparams/60*-9999999.d0/
      data cparams/60*'xxxxxxxxxxxxxxxx'/
c
      common/accum/at10,at21,at32,at43,at54,at65,at76,at70
      data at10,at21,at32,at43,at54,at65,at76,at70/0.,0.,0.,0.,0.,0.,0.,&
     &     0./
      common/ffttime/ft10,ft21,ft32,ft30
      data ft10,ft21,ft32,ft30/0.,0.,0.,0./
c
cryneneriwalstrom:
      parameter(maxcof=35)
      common /bicof/ bcoeff(maxcof,maxcof)
      data (bcoeff(k,1),k=1,2) /1.d0,1.d0/
      data (bcoeff(k,2),k=1,3) /1.d0,2.d0,1.d0/
      end
c
************************************************************************
c
      block data names
c-----------------------------------------------------------------------
c define type codes allowed in MARYLIE
c Written by Petra Schuett in November 1987 based on earlier
c work of Rob Ryne
c-----------------------------------------------------------------------
      include 'impli.inc'
c---------
c commons
c---------
      include 'codes.inc'
c-----------------------------------------------------------------------
c components of master input file
      data ling/'#comment',                                             &
     &          '#beam   ',                                             &
     &          '#menu   ',                                             &
     &          '#lines  ',                                             &
     &          '#lumps  ',                                             &
     &          '#loops  ',                                             &
     &          '#labor  ',                                             &
     &          '#call   ',                                             &
     &          '#const'/
c
c menu entries
c
c 1: simple elements
cryne 7/14/2002 added various elements to be compatible w/ MAD input.
cryne Also added some commonly used synonyms (e.g 'kick' for 'kicker')
cryne Need to work on the thin multipole.
c
      data (ltc(1,j),j=1,75)/
     &  'drft    ','nbnd    ','pbnd    ','gbnd    ','prot    ',         &
     &  'gbdy    ','frng    ','cfbd    ','quad    ','sext    ',         &
     &  'octm    ','octe    ','srfc    ','arot    ','twsm    ',         &
     &  'thlm    ','cplm    ','cfqd    ','dism    ','sol     ',         &
     &  'mark    ','jmap    ','dp      ','recm    ','spce    ',         &
     &  'cfrn    ','coil    ','intg    ','rmap    ','arc     ',         &
     &  'rfgap   ','confoc  ','transit ','interface','rootmap',         &
     &  'optirot ','spare5  ','spare6  ','spare7  ','spare8  ',         &
     &  'marker  ','drift   ','rbend   ','sbend   ','gbend   ',         &
     &  'quadrupole','sextupole','octupole','multipole','solenoid',     &
     &  'hkicker ','vkicker ','kicker  ','rfcavity','elsepararator',    &
     &  'hmonitor','vmonitor','monitor ','instrument','sparem1 ',       &
     &  'rcollimator','ecollimator','yrot    ','srot    ','prot3   ',   &
     &  'beambeam','matrix  ','profile1d','yprofile','tprofile',        &
     &  'hkick   ','vkick   ','kick    ','hpm     ','nlrf    '/
      data (nrp(1,j),j=1,75)/                                           &
     &    1,6,4,6,2,                                                    &
     &    4,5,6,4,2,                                                    &
     &    2,2,2,1,4,                                                    &
     &    6,5,4,5,6,                                                    &
     &    0,0,0,6,1,                                                    &
     &    4,11,6,6,0,                                                   &
     &    11,5,3,5,4,                                                   &
     &    2,0,0,0,0,                                                    &
     &    0,3,13,21,0,                                                  &
     &    6,3,3, 6,13,                                                  &
     &    4,4,5,8,4,                                                    &
     &    2,2,2,2,0,                                                    &
     &    3,3,1,1,2,                                                    &
     &    0,0,7,6,6,                                                    &
     &    4,4,5,3,7/
      data (ncp(1,j),j=1,75)/                                           &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    2,0,0,0,1,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,1,0,0,0,                                                    &
     &    1,1,1,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,0,3,1,1,                                                    &
     &    0,0,0,0,2/
      data (nrpold(1,j),j=1,75)/                                        &
     &    1,6,4,6,2,                                                    &
     &    4,5,6,4,2,                                                    &
     &    2,2,2,1,4,                                                    &
     &    6,5,4,5,6,                                                    &
     &    0,0,0,6,1,                                                    &
     &    4,11,6,6,0,                                                   &
     &    7,4,0,0,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    0,2,0,21,0,                                                   &
     &    5,3,3, 6,2,                                                   &
     &    3,3,4,8,3,                                                    &
     &    1,1,1,1,0,                                                    &
     &    3,3,1,1,0,                                                    &
     &    0,0,0,0,0,                                                    &
     &    3,3,4,3,0/
c
c 2: user-supplied elements
c
      data (ltc(2,j),j=1,20)/                                           &
     &  'usr1    ','usr2    ','usr3    ','usr4    ','usr5    ',         &
     &  'usr6    ','usr7    ','usr8    ','usr9    ','usr10   ',         &
     &  'usr11   ','usr12   ','usr13   ','usr14   ','usr15   ',         &
     &  'usr16   ','usr17   ','usr18   ','usr19   ','usr20   '/
      data (nrp(2,j),j=1,20)/20*6/
      data (ncp(2,j),j=1,20)/20*0/
      data (nrpold(2,j),j=1,20)/20*6/
c
c 3: parameter sets
c
      data (ltc(3,j),j=1,9)/                                            &
     &  'ps1     ','ps2     ','ps3     ','ps4     ','ps5     ',         &
     &  'ps6     ','ps7     ','ps8     ','ps9     '/
      data (nrp(3,i),i=1,9)/9*6/
      data (ncp(3,i),i=1,9)/9*0/
      data (nrpold(3,i),i=1,9)/9*6/
c
c 4: random elements
c    Note: The j indices for these entries must be aligned with
c          the j indices for their "simple element" counterparts.
c          This is done by using "dummy" entries.
c
      data (ltc(4,j),j=1,24)/                                           &
     &  'rdrft   ','rnbnd   ','rpbnd   ','rgbnd   ','rprot   ',         &
     &  'rgbdy   ','rfrng   ','rcfbd   ','rquad   ','rsext   ',         &
     &  'roctm   ','rocte   ','rsrfc   ','rarot   ','rtwsm   ',         &
     &  'rthlm   ','rcplm   ','rcfqd   ','rdism   ','rsol    ',         &
     &  'dummark ','dumjmap ','dumdp   ','rrecm   '/
      data (nrp(4,i),i=1,24)/24*2/
      data (ncp(4,i),i=1,24)/24*0/
      data (nrpold(4,i),i=1,24)/24*2/
c
c 5: random user-supplied elements
c
      data (ltc(5,j),j=1,9)/                                            &
     &  'rusr1   ','rusr2   ','rusr3   ','rusr4   ','rusr5   ',         &
     &  'rusr6   ','rusr7   ','rusr8   ','rusr9   '/
      data (nrp(5,i),i=1,9)/9*2/
      data (ncp(5,i),i=1,9)/9*0/
      data (nrpold(5,i),i=1,9)/9*2/
c
c 6: random parameter sets
c
      data (ltc(6,j),j=1,9)/                                            &
     &  'rps1    ','rps2    ','rps3    ','rps4    ','rps5    ',         &
     &  'rps6    ','rps7    ','rps8    ','rps9    '/
      data (nrp(6,i),i=1,9)/9*2/
      data (ncp(6,i),i=1,9)/9*0/
      data (nrpold(6,i),i=1,9)/9*2/
c
c 7: simple commands
cKMP: Replaced 'sparec5' with 'wrtmap'
cKMP: Replaced 'spacec6' with 'rdmap'
c
      data (ltc(7,j),j=1,75)/                                           &
     &  'rt      ','sqr     ','symp    ','tmi     ','tmo     ',         &
     &  'pmif    ','circ    ','stm     ','gtm     ','end     ',         &
     &  'ptm     ','iden    ','whst    ','inv     ','tran    ',         &
     &  'revf    ','rev     ','mask    ','num     ','rapt    ',         &
     &  'eapt    ','of      ','cf      ','wnd     ','wnda    ',         &
     &  'ftm     ','wps     ','time    ','cdf     ','bell    ',         &
     &  'wmrt    ','wcl     ','paws    ','inf     ','dims    ',         &
     &  'zer     ','sndwch  ','tpol    ','dpol    ','cbm     ',         &
     &  'poisson ','preapply','midapply','autoapply','autoconcat',      &
     &  'rayscale','beam    ','units   ','autoslice','verbose ',        &
     &  'mask6   ','arcreset','symbdef ','particledump','raytrace',     &
     &  'autotrack','sckick','moments  ','maxsize  ','reftraj',         &
     &  'initenv  ','envelopes','contractenv','setreftraj','setarclen', &
     &  'wakedefault','emittance','matchenv','fileinfo','egengrad',     &
     &  'wrtmap  ','rdmap   ','sparec7','sparec8','sparec9'/
      data (nrp(7,i),i=1,75)/                                           &
     &   6,0,2,4,1,                                                     &
     &   3,6,1,2,0,                                                     &
     &   5,0,2,0,0,                                                     &
     &   1,0,6,4,4,                                                     &
     &   3,6,6,5,6,                                                     &
     &   4,2,3,1,0,                                                     &
     &   2,3,0,5,6,                                                     &
     &   3,1,6,5,3,                                                     &
     &   12,0,0,0,0,                                                    &
     &   6,10,5,2,1,                                                    &
     &   6,0,0,7,12,                                                    &
     &   0,0,5,6,3,                                                     &
     &   15,5,0,10,2,                                                   &
     &   4,3,2,1,10,                                                    &
     &   0,0,0,0,0/
      data (ncp(7,i),i=1,75)/                                           &
     &   2,0,0,0,0,                                                     &
     &   1,0,0,0,1,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   17,2,2,2,1,                                                    &
     &   0,1,2,3,0,                                                     &
     &   0,0,0,4,5,                                                     &
     &   4,0,9,9,2,                                                     &
     &   1,7,0,2,1,                                                     &
     &   1,6,1,0,7,                                                     &
     &   3,2,0,0,0/
      data (nrpold(7,i),i=1,75)/                                        &
     &   6,0,2,4,1,                                                     &
     &   3,6,1,2,0,                                                     &
     &   5,0,2,0,0,                                                     &
     &   1,0,6,4,4,                                                     &
     &   3,6,6,5,6,                                                     &
     &   4,2,3,1,4,                                                     &
     &   2,3,0,5,6,                                                     &
     &   3,1,6,5,3,                                                     &
     &   6,0,0,0,1,                                                     &
     &   6,10,4,1,0,                                                    &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0/
c
c 8: advanced commands
c
      data (ltc(8,j),j=1,39)/                                           &
     &  'cod     ','amap    ','dia     ','dnor    ','exp     ',         &
     &  'pdnf    ','psnf    ','radm    ','rasm    ','sia     ',         &
     &  'snor    ','tadm    ','tasm    ','tbas    ','gbuf    ',         &
     &  'trsa    ','trda    ','smul    ','padd    ','pmul    ',         &
     &  'pb      ','pold    ','pval    ','fasm    ','fadm    ',         &
     &  'sq      ','wsq     ','ctr     ','asni    ','pnlp    ',         &
     &  'csym    ','psp     ','mn      ','bgen    ','tic     ',         &
     &  'ppa     ','moma    ','geom    ','fwa     '/
      data (nrp(8,i),i=1,39)/                                           &
     &   6,6,5,5,3,                                                     &
     &   5,5,5,5,5,                                                     &
     &   5,4,6,1,2,                                                     &
     &   6,6,6,3,3,                                                     &
     &   3,5,3,6,6,                                                     &
     &   4,6,4,6,4,                                                     &
     &   1,3,2,12,6,                                                     &
     &   6,6,6,6/
      data (ncp(8,i),i=1,39)/                                           &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,2,0,                                                     &
     &   0,0,0,0/
      data (nrpold(8,i),i=1,39)/                                        &
     &   6,6,5,5,3,                                                     &
     &   5,5,5,5,5,                                                     &
     &   5,4,6,1,2,                                                     &
     &   6,6,6,3,3,                                                     &
     &   3,5,3,6,6,                                                     &
     &   4,6,4,6,4,                                                     &
     &   1,3,2,6,6,                                                     &
     &   6,6,6,6/
c
c 9: procedures and fitting and optimization
c
      data (ltc(9,j),j=1,37)/                                           &
     &  'bip     ','bop     ','tip     ','top     ',                    &
     &  'aim     ','vary    ','fit     ','opt     ',                    &
     &  'con1    ','con2    ','con3    ','con4    ','con5    ',         &
     &  'mrt0    ',                                                     &
     &  'mrt1    ','mrt2    ','mrt3    ','mrt4    ','mrt5    ',         &
     &  'fps     ',                                                     &
     &  'cps1    ','cps2    ','cps3    ','cps4    ','cps5    ',         &
     &  'cps6    ','cps7    ','cps8    ','cps9    ',                    &
     &  'dapt    ','grad    ','rset    ','flag    ','scan    ',         &
     &  'mss     ','spare1  ','spare1  '/
      data (nrp(9,i),i=1,37)/                                           &
     &   1,1,1,1,                                                       &
     &   6,6,6,6,                                                       &
     &   6,6,6,6,6,                                                     &
     &   1,                                                             &
     &   6,6,6,6,6,                                                     &
     &   1,                                                             &
     &   6,6,6,6,6,                                                     &
     &   6,6,6,6,                                                       &
     &   4,6,5,6,6,                                                     &
     &   6,6,6/
      data (ncp(9,i),i=1,37)/                                           &
     &   0,1,1,1,                                                       &
     &   0,0,0,0,                                                       &
     &   0,0,0,0,0,                                                     &
     &   0,                                                             &
     &   0,0,0,0,0,                                                     &
     &   0,                                                             &
     &   0,0,0,0,0,                                                     &
     &   0,0,0,0,                                                       &
     &   0,0,0,0,0,                                                     &
     &   0,0,0/
      data (nrpold(9,i),i=1,37)/                                        &
     &   1,1,1,1,                                                       &
     &   6,6,6,6,                                                       &
     &   6,6,6,6,6,                                                     &
     &   1,                                                             &
     &   6,6,6,6,6,                                                     &
     &   1,                                                             &
     &   6,6,6,6,6,                                                     &
     &   6,6,6,6,                                                       &
     &   4,6,5,6,6,                                                     &
     &   6,6,6/
c
      end
c
c**********************************************************************
c
      subroutine buffin(th,tmh,thsave,tmhsav,lumpno)
c-----------------------------------------------------------------------
c     stores a map and polynomial coeffs into a buffer
c  Written by J. Howard, Fall 1986
c-----------------------------------------------------------------------
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'stack.inc'
c
      dimension th(monoms),tmh(6,6)
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c
      do 100 i = 1,6
      do 100 j = 1,6
      tmhsav(i,j,lumpno) = tmh(i,j)
100   continue
c
      do 200 i = 1,monoms
      thsave(i,lumpno) = th(i)
200   continue
      return
      end
************************************************************************
      subroutine bufout(thsave,tmhsav,th,tmh,lumpno)
c-----------------------------------------------------------------------
c     reads a map and polynomial coefficients out of a buffer
c  Written by J. Howard, Fall 1986
c-----------------------------------------------------------------------
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'stack.inc'
c
      dimension th(monoms),tmh(6,6)
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c
      do 100 i = 1,6
      do 100 j = 1,6
      tmh(i,j) = tmhsav(i,j,lumpno)
100   continue
c
      do 200 i = 1,monoms
      th(i) = thsave(i,lumpno)
200   continue
c
      return
      end
************************************************************************
      subroutine cnumb(string,num,lnum)
c-----------------------------------------------------------------------
c  This routine finds out, whether string codes an integer number.
c  if so, it converts it to num
c
c  Input: string character*16 input string
c  Output:num    integer      correspondent number
c         lnum   logical      =.true. if string is a number
c
c  Author: Petra Schuett
c          October 19, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
c
cryne 5/5/2006 declared slen an integer
      integer num,slen
cryne 8/5/2004      character string*16
      character (LEN=*) :: string
      logical lnum
c
      character*10 digits
cryne data digits /'0123456789'/
cryne save digits
      digits='0123456789'
cryne 8/5/2004:
      slen = LEN(string)
c-----------------
c The first char may be minus or digit
c      write(jodf,*)'cnumb: string= ',string
      if(index(digits,string(1:1)).eq.0 .and. string(1:1).ne.'-') then
        lnum=.false.
c        write(jodf,*)'first char is no digit'
        return
      endif
c All other characters must be digits...
cryne 8/5/2004      do 1 k=2,16
      do 1 k=2,slen
      if(index(digits,string(k:k)).eq.0) then
c ...or trailing blanks
cryne 8/5/2004        if(string(k:16).ne.' ') then
        if( string(k:k).ne.' ') then
c        write(jodf,*)'string(',k,':16) is not blank'
c        write(jodf,*)'string(',k,':slen) is not blank'
         lnum=.false.
         return
        else
         goto 11
        endif
      endif
   1  continue
c string is a number
  11  read(string,*,err=999) num
      lnum=.true.
      return
c-----------------
c error exit
 999  write(jof ,99) string
      write(jodf,99) string
 99   format(' ---> error in cnumb: string ',a16,' could not be',       &
     &       ' converted to number')
      call myexit
      end
c
************************************************************************
c
      subroutine comwtm(j,jrep)
c-----------------------------------------------------------------------
c  routine to combine the jth lump with the total map n times
c  Written by J. Howard, Fall 1986
c  Modified by Alex Dragt, 15 June 1988, to save storage
c-----------------------------------------------------------------------
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
      include 'map.inc'
      include 'core.inc'
c----------
c start
c----------
c make n positive in case lump has a negative repetition number
      n=iabs(jrep)
c
      do 100 i=1,n
      call concat(th,tmh,thl(1,j),tmhl(1,1,j),th,tmh)
  100 continue
      return
      end
c
*******************************************************************
c
      subroutine cqlate(icfile,norder,ntimes,nwrite,isend)
c
c-----------------------------------------------------------------------
c  circulate ntimes times; print every nwrite
c
c  input : icfile (integer) = filenumber for reading input rays
c          norder (integer) = NOT USED
c          ntimes (integer) = number of turns through actual line
c          nwrite (integer) = modulo for writing coordinates to jfcf
c          isend  (integer) = 1 output only to jof
c                           = 2 output only to jodf
c                           = 3 both
c  in common /files/: jfcf  < 0 full precision coordinates
c                           > 0 standard format coordinates
c
c  Written by Rob Ryne ca 1984
c  slightly changed by Petra Schuett (labels, if-then-else ...)
c                      October 30,1987
c------------------------------------------------------------------------
      use parallel
      use acceldata
      use rays
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'talk.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'map.inc'
      include 'deriv.inc'
      include 'core.inc'
      include 'nturn.inc'
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
      character*16 string(5),str,ubuf
      logical     ljof,ljodf
cryne:
cryne dimension zi(6),zf(6) this is defined in rays.inc
c-----------------------------------------------------------------------
c  start routine
c-----------------------------------------------------------------------
c  see if a loop exists
      if(nloop.le.0) then
      write(jof ,*) ' error from cqlate: no loop has been specified'
      write(jodf,*) ' error from cqlate: no loop has been specified'
      return
      endif
c
c  read initial conditions, if requested:
c     if(icfile.gt.0) call raysin(icfile)
      if(icfile.gt.0)then
      scaleleni=0.d0
      scalefrqi=0.d0
      scalemomi=0.d0
      call raysin(icfile,scaleleni,scalefrqi,scalemomi,ubuf,jerror)
      call rbcast(scaleleni)
      call rbcast(scalefrqi)
      call rbcast(scalemomi)
      endif
c check for initial conditions:
      if(nrays.eq.0)then
        if(idproc.eq.0)then
           write(6,*)'error(cqlate) : nrays=0'
           write(6,*)'forgot to generate or read initial conditions?'
        endif
      call myexit
      endif
c--------------------
c  write items in loop, if requested:
      ljof  = isend.eq.1 .or. isend.eq.3
      ljodf = isend.eq.2 .or. isend.eq.3
      if (ljof .or. ljodf) then
c  write loop name
      if(ljof ) write(jof ,510) ilbl(nloop)
      if(ljodf) write(jodf,510) ilbl(nloop)
  510 format(1h ,'circulating through ',a16,' :')
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
  511 format(' ',5(1x,a16))
   1  continue
      endif
c--------------------
c  circulation procedure
c
c  set up control indices
      ibrief=0
      kwrite=nwrite
      if(kwrite.eq.0)kwrite=1
c--------------------
c  circulation through loop
      do 1000 nturn=1,ntimes
c----------
c  handle each loop-element
      do 100 kk=1,joy
c  check to see if all particles lost
c  Go at least one turn...
      if(nturn.gt.1) then
      if(nlost.ge.nraysp) then
       write(jof,*) 'all particles lost on processor ',idproc
       return
      endif
      endif
c
c check to see if item is a lump, user routine, or element
      ip = mim(kk)
      if(ip.gt.5000)then
c procedure for a user routine
        nip=ip-5000
        if(nt2(nip).eq.1) call user1(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.2) call user2(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.3) call user3(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.4) call user4(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.5) call user5(pmenu(1+mpp(nip)))
c
      else if(ip.lt.0) then
c procedure for turtling thru elements (which may be in lines):
cryne 08/15/2001 autoslicing does not apply when cirulating thru a loop
cryne note: this has not been tested as of 8/15/2001:
        jslice=1
        nslices=1
        slfrac=1.
c       ntrk=0
c       if(nspchset.eq.1)ntrk=1
c-----------------------------------------------------
c       call isthick(kt1,kt2,ithick,thlen,0)
c       if(iautosl.lt.0.and.ithick.eq.1)nslices=pmenu(nrp(1,kt2))
c       if(iautosl.gt.0.and.ithick.eq.1)nslices=iautosl
c       slfrac=1./nslices
c       if(lautoap2.eq.1 .or. nspchset.eq.1)slfrac=slfrac*0.5d0
c-----------------------------------------------------
        narg1=1+mpp(-ip)
        narg2=1+mppc(-ip)
        call lmnt(nt1(-ip),nt2(-ip),pmenu(narg1),cmenu(narg2),          &
     &  1,jslice,nslices,slfrac,0)
c
      else
c procedure for a lump
       do 110 nn=1,nraysp
c check to see if particle has already been lost
       if (istat(nn).eq.0) then
         do 112 l=1,6
  112    zi(l)=zblock(l,nn)
c call symplectic tracker
         jwarn=0
         call evalsr_old(tmhl(1,1,ip),zi,zf,dfl(1,1,ip),                    &
     &               rdfl(1,1,ip),rrjacl(1,1,1,ip))
c check to see if particle was 'lost' by the symplectic ray tracer
         if (jwarn.ne.0) then
           istat(nn)=nturn
           nlost=nlost + 1
           ihist(1,nlost)=nturn
           ihist(2,nlost)=nn
         endif
c copy results of ray trace into storage array
         do 114 l=1,6
  114    zblock(l,nn)=zf(l)
       endif
  110  continue
      endif
  100 continue
c
c end of one turn
c----------
c check to see if results should be written out
      if(mod(nturn,kwrite).eq.0) then
      do 200 m1=1,nraysp
      if(istat(m1).eq.0) then
c  procedure for writing out results of ray trace
c
c  procedure for standard format
       if(jfcf.gt.0) then
        write(jfcf,520)(zblock(m2,m1),m2=1,6)
  520   format(6(1x,1pe12.5))
c
c  procedure for full precision
       else if(jfcf.lt.0) then
        do 525 m2=1,6
        write(-jfcf,*) zblock(m2,m1)
  525   continue
       endif
c
      endif
  200 continue
      endif
c
 1000 continue
c end of loop thru turns
c-----------------------------------------------------------------------
      return
      end
c
***********************************************************************
c
      subroutine cread(kbeg,msegm,line,string,lfound)
c-----------------------------------------------------------------------
c  This routine searches line(kbeg:) for the next string; strings are
c  delimited by ' ','*' or ','.
c
c  Input: line   character*(*)input line
c         kbeg   integer      line is only searched behind kbeg
c                             if a string is found, kbeg is set to
c                             possible start of next string
c         msegm  integer      segment number. for msegm=1 (comment
c                             section), the length of a string is not
c                             checked.
c  output:string character*(*)found string
c         lfound logical      =.true. if a string was found
c
c  Written by Rob Ryne ca 1984
c  Rewritten by Petra Schuett
c          October 19, 1987
cryne
c  Modified by Rob Ryne July 2002 to parse lines (with some exceptions)
c  in the Standard Input Format
c  Modified by David Serafini to allow lines longer than 80 chars
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
c
      integer kbeg,msegm
      character line*(*), string*(*)
      logical lfound
      integer llen,slen

      llen = LEN(line)
      slen = LEN(string)
c
c look for first character of string
cryne 3/6/2004 but first check for *(
cryne 5/4/2006 but don't do this if in #comment
      if(msegm.eq.1)goto 666
cryne 5/4/2006 and only do this if the word 'line' appears:
      if( index(line,'line').eq.0 )  goto 666
c
      do k=kbeg,llen-2
cryne   this version of the parser cannot handle lines defined as
cryne   repeated quantities between parenthesese
        if( line(k:k+1).eq.'*(' .or.
     &      line(k:k+2).eq.'* (' )then
        write(6,*)'error in cread:'
        write(6,*)'this version cannot parse a line of the form:'
        write(6,*)'N*(...)  or N* (...)'
        call myexit
        endif
      enddo
  666 continue
c
      do 1 k=kbeg,llen
        if(     (line(k:k).ne.' ')
     &     .and.(line(k:k).ne.'*')
cryne-
     &     .and.(line(k:k).ne.';')
     &     .and.(line(k:k).ne.':')
     &     .and.(line(k:k).ne.'=')
     &     .and.(line(k:k).ne.'(')
     &     .and.(line(k:k).ne.')')
cryne-
     &     .and.(line(k:k).ne.',')) then
c found:
          lfound=.true.
c dont run off the end of either char variable
          lmax  = min(k+slen,llen)
c look for delimiter
          do 2 l=k,lmax
            if((    line(l:l).eq.' ')
     &         .or.(line(l:l).eq.'*')
cryne-
     &         .or.(line(l:l).eq.':')
     &         .or.(line(l:l).eq.';')
     &         .or.(line(l:l).eq.'=')
     &         .or.(line(l:l).eq.'(')
     &         .or.(line(l:l).eq.')')
cryne-
     &         .or.(line(l:l).eq.',')) then
c if found, string is known
              string=line(k:l-1)
              kbeg  =l+1
c done.
              return
            endif
    2     continue
c no delimiter behind string found...
          string=line(k:lmax)
c ...end of line
          if(lmax.eq.llen) then
            kbeg=llen
c ...or string too long (ignored for comment section)
          else if (msegm.ne.1) then
            write(jof,99) kbeg
   99       format(' ---> warning from cread:'/,
     &            '      the following line has a very long string at ',
     &            'position ',i2,' :')
            write(jof,*) line
            kbeg=lmax+1
          endif
c ...anyway, its done.
          return
        endif
    1 continue
c No string was found after all.
      lfound=.false.
      return
      end
c
************************************************************************
      subroutine dump(iu)
c-----------------------------------------------------------------------
c  this subroutine dumps the status of the common-blocks
c  as they are filled by dumpin
c  input: iu   output-file
c
c  written by Petra Schuett
c             October 26, 1987
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'codes.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c start routine
c-----------------------------------------------------------------------
c comments
      if(np.ne.0) then
        write(iu,500) ling(1)
        write(iu,*) (mline(i),i=1,npcom)
      endif
c--------------------
c  beam
      write(iu,500) ling(2)
      write(iu,*) brho
      write(iu,*) gamm1
      write(iu,*) achg
      write(iu,*) sl
c--------------------
c  menu
      write(iu,500) ling(3)
      do 10 k=1,na
       write(iu,520) lmnlbl(k),ltc(nt1(k),nt2(k))
       imax=nrp(nt1(k),nt2(k))
       if(imax.ne.0) then
        write(iu,522)(pmenu(i+mpp(k)),i=1,imax)
       endif
   10 continue
c--------------------
c  lines,lumps,loops
      if(nb.ne.0) then
       do 40 ii=2,4
       write(iu,500) ling(ii+2)
        do 40 k=1,nb
        if(ityp(k).eq.ii) then
         write(iu,530) ilbl(k)
         write(iu,532)(irep(l,k),icon(l,k),l=1,ilen(k))
        endif
   40  continue
      endif
c--------------------
c labor
      if(noble.ne.0) then
      write(iu,500) ling(7)
      do 100 j=1,noble
       write(iu,540) num(j),latt(j)
  100 continue
      endif
      return
c-----------------------------------------------------------------------
c format
c-----------------------------------------------------------------------
  500 format(1h ,a8)
!  510 format(1h ,a80)
  520 format(1h ,1x,a16,1x,a16)
  522 format((1h ,3(1x,1pg22.15)))
  530 format(1h ,1x,a8)
  532 format((1h ,1x,5(i5,'*',a8),1x,:'&'))
  540 format(1h ,1x,i4,'*',a8)
      end
************************************************************************
      subroutine initst(string,nrept)
c-----------------------------------------------------------------------
c initialize stack, putting element 'string' with rep. factor nrept
c on top
c
c input: string character*(*) initial stack element (loop,line or lump!)
c        nrept  integer     rep. factor
c
c  Petra Schuett, October 30,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      use parallel, only : idproc
      include 'impli.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      character string*(*)
c-------
c start
c-------
      np = 1
      lstac(1) = string
      loop (1) = nrept
      call lookup(string,ntype,ith)
      if(ntype.eq.1 .or. ntype.eq.5) then
       if(idproc.eq.0)write(jof,510) string
  510  format(' >>>>warning from initst: stack element',a  ,' is an',   &
     &        ' element or an unused label')
      else
       nslot(1) = newsl(1,ith)
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
c-----------------------------------------------------------------------
      return
      end
************************************************************************
c
      subroutine lmnt(nt1,nt2,prms,crms,ntrk,jslice,jsltot,slfrac,ihalf)
c-----------------------------------------------------------------------
c this routine actually switches to the routines that handle single
c elements, commands, etc.
c
c input: nt1  = group type code of element
c        nt2  = index of element in its group
c        prms = array of numerical (double) params for this element
c        crms = array of character params for this element
c        ntrk = 0 if accumulated transfer map is to be constructed
c             = 1 if tracking
c Modified by Rob Ryne 08/15/2001 to handle sliced elements
c jslice=present slice number; jsltot=total number of slices,
c slfrac=fraction by which the element length is multiplied for slicing
c ihalf = 1 (first half) or 2 (second half) if a slice is split in half,
c         as when performing space-charge kicks or "midapply" operations
c ihalf = 0 if slices are not split in half.
c
cryne: this routine should be cleaned up in order to (among other things)
c have some consistency regarding the units of angles. Most of the
c inconsistency is historical: MaryLie uses degrees in the argument lists
c for most (all?) bending magnet routines, but the arot is in radians.
c MAD uses radians for everything, but, rather than change the MaryLie
c routines, the normal MAD input is converted to degrees.
c
c  Written by Rob Ryne ca 1984
c  modified by J. Howard 7-87
c              P. Schuett 11-87
c              A. Dragt 6/15/88
c-----------------------------------------------------------------------
      use beamdata
      use acceldata, only : slicetype,sliceprecedence,slicevalue
      use rays
      use lieaparam, only : monoms
      use spchdata
      use e_gengrad
      use multitrack
      include 'impli.inc'
c--------
c commons
c--------------------
      include 'map.inc'
      include 'parset.inc'
      include 'files.inc'
      include 'codes.inc'
      include 'pie.inc'
      include 'infin.inc'
      include 'zeroes.inc'
      include 'frnt.inc'
      include 'previous.inc'
      include 'setref.inc'
      include 'fitbuf.inc'
      common/rfxtra/zedge,gaplen,thetastore
      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0
      common/poiblock1/nspchset,nsckick
      common/gxyzsave/xmin0,xmax0,ymin0,ymax0,zmin0,zmax0,kfixbdy,madegr
      integer           idirectfieldcalc,idensityfunction,isolve
     &                 ,anagpatchsize,anagrefratio,anagsmooth
      common/newpoisson/idirectfieldcalc,idensityfunction,isolve
     &                 ,anagpatchsize,anagrefratio,anagsmooth
      character*256     chombofilename
      common/chombochar/chombofilename
      common/dirichlet/idirich
      common/autotrk/lautotrk,ntrktype,ntrkorder
      common/envstuff/nenvtrk
      character*16 fname,fname1,fname2,fname3,fname4
      character*16 estrng,autostr,icname,fcname
      character*5 anum
      common/autopnt/autostr(3)
      common/autolog/lautoap1,lautoap2,lautoap3,lrestrictauto
      common/autocon/lautocon
      common/showme/iverbose
      real*8 mhprev
      common/prevmap/hprev(monoms),mhprev(6,6)
      common/bfileinfo/ifileinfo
      common/xtrajunk/jticks1,iprinttimers
cryne
      character*3 kynd
c--------------------------------
c parameter types and dimensions
c--------------------------------
      integer nt1,nt2,ntrk
      character*16 crms,cparams
      character*5 aseq
      dimension prms(*),crms(*)
c-----------------
c local variables
c-----------------
c maps
      double precision h(monoms),ht1(monoms),ht2(monoms),ht3(monoms)
      double precision mh(6,6),mht1(6,6),mht2(6,6),mht3(6,6)
      double precision dreftraj(6)
c
c leading .. and trailing fringe field: lfrn (tfrn) .ne.0, if it is
c to be taken into account
      integer lfrn,tfrn
c Input args nt1,nt2,prms must not be changed. work with kt1,kt2,p.
c     dimension p(6),ptmp(6)
      dimension p(60),ptmp(6)
cryne 7/13/2002
      dimension coefpset(6)
cryne 12 Nov 2003 added common/scparams/
      common/scparams/rparams(60),cparams(60)
      equivalence (p1,p(1)),(p2,p(2)),(p3,p(3)),(p4,p(4)),(p5,p(5)),    &
     &            (p6,p(6))
      data imsg1/0/
      save imsg1
      logical, save :: first_entry = .true.
cKMP: 8 Noc 2006 - added clinestrng
      character*132 clinestrng
c-----------------------------------------------------------------------
c start
c----------------
      if(iverbose.eq.2.and.idproc.eq.0)                                 &
     &write(6,*)'inside LMNT with nt1,nt2=',nt1,nt2
c     write(6,*)'(LMNT) nt1,nt2,ntrk,jslice,jsltot,slfrac,ihalf='
c     write(6,*)        nt1,nt2,ntrk,jslice,jsltot,slfrac,ihalf
c     write(6,*)'reftraj(5),reftraj(6)=',reftraj(5),reftraj(6)
c-----------------------------------------------------------------------
      multitrac=0   !cryne May 21, 2006
c
c first-entry initializations
c
      if (first_entry) then
        first_entry = .false.
        nullify(eggrdata%zvals)
        nullify(eggrdata%G1c); nullify(eggrdata%G1s)
        nullify(eggrdata%G2c); nullify(eggrdata%G2s)
        nullify(eggrdata%G3c); nullify(eggrdata%G3s)
      endif
c
c Preliminary procedure for random elements or commands
c
      if(nt1.eq.4 .or. nt1.eq.5 .or. nt1.eq.6) then
        nfile = nint(prms(1))
        iecho = nint(prms(2))
        kt1 = nt1 - 3
        kt2 = nt2
        call randin(nfile,kt1,kt2,p)
c
c Echo back if requested
      if( iecho.eq.1 .or. iecho.eq.3) then
        write(jof, 7005) nfile,kt1,kt2
        write(jof ,*)(p(i),i=1,nrp(kt1,kt2))
      endif
      if( iecho.eq.2 .or. iecho.eq.3) then
        write(jodf,7005) nfile,kt1,kt2
        write(jodf,*)(p(i),i=1,nrp(kt1,kt2))
      endif
 7005   format(' random lmnt check:nfile,kt1,kt2=',3i5/
     &         ' parameters found:')
c
c
      else
c
c Preliminary procedure for all other type codes
c
c
        kt1 = nt1
        kt2 = nt2
        nparams=nrp(nt1,nt2)
        ncparams=ncp(nt1,nt2)
c       write(6,*)'nparams=',nparams
        if(nparams.gt.0)then
          p(1:nparams)=prms(1:nparams)
        endif
      endif
c
c Select appropriate action depending on value of kt1,kt2
c
      go to (11,12,13,14,15,16,17,18,19), kt1
c
c     1: simple elements *************************************
c
11    continue
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
c          'rfgap   ','confoc  ','transit','interface','rootmap',
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
c          'hkick   ','vkick   ','kick    ','hpm     ','nlrf    '/
     &      171,       172,       173,       174,       175),kt2
c
c 'drft    ': drift
c
101   continue
      if(slfrac.ne.1.0)p1=p1*slfrac
c     write(6,*)'DRFT:jslice,jsltot,slfrac=',jslice,jsltot,slfrac
      call drift3(p1,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
c note: under the assumption that omega*sl/c=1,  it follows that omega=c/sl
c therefore, multiplying by c/sl is the correct thing to do
c to make reftraj(5) dimensionless under this assumption.
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'nbnd    ': normal entry bend
c the map for this element is of the form
c gfrngg*nbend*gfrngg with the leading and trailing fringe field
c maps optional
c
102   continue
c     angdeg=p1
      gap=p2
      xk1=p3
      rho=brho/p4
      lfrn=nint(p5)
      tfrn=nint(p6)
c compute nbend
      if(slfrac.ne.1.0)p1=p1*slfrac
      call nbend(rho,p1,h,mh)
c put on leading fringe field
cryne
      azero=0.
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
        if(lfrn.ne.0)then
cryne   call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
        call gfrngg(azero,rho,1,ht1,mht1,gap,xk1)
        call concat(ht1,mht1,h,mh,h,mh)
        endif
      endif
c put on trailing fringe field
      if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
        if(tfrn.ne.0)then
cryne   call gfrngg(0.d0,rho,2,ht1,mht1,gap,xk1)
        call gfrngg(azero,rho,2,ht1,mht1,gap,xk1)
        call concat(h,mh,ht1,mht1,h,mh)
        endif
      endif
c==============
      call set_pscale_mc(h,mh)
      arclen=arclen+rho*p1*(pi/180.)
      refprev=reftraj
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c  'pbnd    ':
c  parallel faced bend, including leading and trailing
c  pole face rotations and fringe fields.
c  Symmetric bends only are permitted under this option.
c  For parallel-faced magnet with asymmetric entry and exit,
c  use the general bending magnet.
c  The map for the parallel faced bend is of the form
c  prot*gfrng*pbend*gfrng*prot
c
103   continue
      psideg=p1/2.d0
      gap=p2
      xk1=p3
      rho=brho/p4
c
cryne Jan 5, 2005      if(jsltot.eq.1)then
      if(jsltot.eq.1 .and. ihalf.eq.0)then
c procedure for a complete pbnd (no slices): [Jan 5, 2005 and not split in half]
c
cryne 3/17/2004
c but it still could be cut in half due to space charge or autoapply. So...
cryne Jan 5, 2005        if(slfrac.ne.1.0)p1=p1*slfrac
c
        call pbend(rho,p1,h,mh)
c put on the leading fringe field
        call gfrngg(psideg,rho,1,ht1,mht1,gap,xk1)
        call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
        call prot(psideg,1,ht1,mht1)
        call concat(ht1,mht1,h,mh,h,mh)
c put on trailing fringe field
        call gfrngg(psideg,rho,2,ht1,mht1,gap,xk1)
        call concat(h,mh,ht1,mht1,h,mh)
c put on trailing prot
        call prot(psideg,2,ht1,mht1)
        call concat(h,mh,ht1,mht1,h,mh)
      else
c
c procedure for a pbnd in slices.
c This follows the approach due to Dragt that is described
c in Exhibit 6.6.3 and section 10.9 of the MaryLie manual.
        if(slfrac.ne.1.0)p1=p1*slfrac
        azero=0.
        call nbend(rho,p1,h,mh)
        if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
c         put on the leading gbody
          call gbend(rho,azero,psideg,azero,ht1,mht1)
          call concat(ht1,mht1,h,mh,h,mh)
c         put on the leading fringe field
          call gfrngg(psideg,rho,1,ht1,mht1,gap,xk1)
          call concat(ht1,mht1,h,mh,h,mh)
c         put on leading prot
          call prot(psideg,1,ht1,mht1)
          call concat(ht1,mht1,h,mh,h,mh)
        endif
        if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
c         put on trailing gbody
          call gbend(rho,azero,azero,psideg,ht1,mht1)
          call concat(h,mh,ht1,mht1,h,mh)
c         put on trailing fringe field
          call gfrngg(psideg,rho,2,ht1,mht1,gap,xk1)
          call concat(h,mh,ht1,mht1,h,mh)
c         put on trailing prot
          call prot(psideg,2,ht1,mht1)
          call concat(h,mh,ht1,mht1,h,mh)
        endif
      endif
c==============
cryne 08/15/2001
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c 'gbnd    ': general bending magnet
c
c  The map for the general bend is of the form
c  prot*gfrng*gbend*gfrng*prot
c
104   continue
      if(jsltot.ne.1)then
        write(6,*)'error: gbnd not yet available with slices!'
        stop
      endif
      gap=p4
      xk1=p5
      if(p6.eq.0.d0)then
        write(6,*)'Error (gbnd): B field equal zero not allowed'
        stop
      endif
      rho=brho/p6
c compute gbend
      call gbend(rho,p1,p2,p3,h,mh)
c
c put on the leading fringe field
      call gfrngg(p2,rho,1,ht1,mht1,gap,xk1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
      call prot(p2,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on trailing fringe field
      call gfrngg(p3,rho,2,ht1,mht1,gap,xk1)
      call concat(h,mh,ht1,mht1,h,mh)
c put on trailing prot
      call prot(p3,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c 'sbend   ': MAD sector bend
c  similar to MaryLie cfbd, but with arbitrary entrance/exit angles
c
c  The map for this element is treated as the following:
c  prot*cgfrngg*cgbend*cgfrngg*prot
c where:
c  prot is a pole face rotation
c  cgfrngg is a fringe field (analgous to that used in cfbd)
c  cgbend is a gbend but w/ multipole coefficients
c  note that, in analogy with the fact that gbend=hpf1*nbend*hpf2,
c  in this case we set cgbend=hpf1*(cfbd w/ no fringe fields)*hpf2,
c  i.e.                cgbend=hpf1*cfbend*hpf2
c
1045  continue
      bndang=p1
      rho=brho/p2
      psideg=p(3)
      phideg=p(4)
      lfrn=nint(p5)
      tfrn=nint(p6)
      gap1=prms(7)
      gap2=prms(8)
      fint1=prms(9)
      fint2=prms(10)
      iopt=nint(prms(11))
      ipset=nint(prms(12))
c
      do j=1,6
      coefpset(j)=0.d0
      enddo
      if( (ipset.gt.0).and.(ipset.le.maxpst))then
        do j=1,6
        coefpset(j)=pst(j,ipset)
c       write(6,*)j,coefpset(j)
        enddo
      else
c       write(6,*)'In lmnt at sbend; j,coefpset(j)='
        do j=1,6
        if(prms(j+12).ne.0.d0)coefpset(j)=prms(j+12)
c       write(6,*)j,coefpset(j)
        enddo
      endif
c
      tiltang=prms(19)
      if(tiltang.ne.0. .and. idproc.eq.0)then
        write(6,*)'WARNING(sbend): tilt keyword being tested'
      endif
      myorder=nint(prms(20))
      if(idproc.eq.0)then
        if(myorder.ne.5)write(6,*)'(sbend) order = ',myorder
      endif
      azero=0.d0
c----------------
c compute gbend
c     call cgbend(rho,bndang,psideg,phideg,iopt,coefpset,h,mh)
      aldeg=bndang-psideg-phideg
      ptmp(1)=aldeg+psideg+phideg    ! = bndang
      if(slfrac.ne.1.0)ptmp(1)=ptmp(1)*slfrac
      ptmp(2)=brho/rho
      ptmp(3)=0.
      ptmp(4)=0.
      ptmp(5)=iopt
cryne 1 August 2004 ptmp(6)=0.    !not used
      ptmp(6)=myorder
c
c     ptmp(1)=psideg
c     call cfbend(ptmp,coefpset,ht2,mht2)
c     ptmp(1)=aldeg
c     call cfbend(ptmp,coefpset,ht3,mht3)
c     call concat(ht2,mht2,ht3,mht3,ht3,mht3)
c     ptmp(1)=phideg
c     call cfbend(ptmp,coefpset,ht2,mht2)
c     call concat(ht2,mht2,ht3,mht3,ht3,mht3)
c     if(idproc.eq.0)then
c     write(6,*)'call cfbend:aldeg,slfrac,ptmp(1)=',aldeg,slfrac,ptmp(1)
c     endif
c***
cryne 12/22/2004  mods per Marco Venturini bypass the slow cfbend routines
c if there are no multipole coefficients (note that, at the moment, the
c cfbend routine is only third order)
      scpset=0.d0
      do ii=1,6
        scpset=scpset+coefpset(ii)
      enddo
      if(scpset.eq.0.d0)then    ! use nbend if field is pure dipole
        slang=bndang
        if(slfrac.ne.1.0)slang=slang*slfrac
        call nbend(rho,slang,h,mh)
      else
        call cfbend(ptmp,coefpset,h,mh)
      endif
c***
c----------------
c
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
c put on the leading gbend
      call gbend(rho,azero,psideg,azero,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on the leading fringe field
      call cgfrngg(1,psideg,rho,lfrn,gap1,fint1,iopt,coefpset,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
      call prot(psideg,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on leading arot
      call arot(tiltang,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      endif
c
      if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
c put on the trailing gbend
      call gbend(rho,azero,azero,phideg,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c put on trailing fringe field
      call cgfrngg(2,phideg,rho,tfrn,gap2,fint2,iopt,coefpset,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c put on trailing prot
      call prot(phideg,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c put on trailing arot
      tiltang=-tiltang
      call arot(tiltang,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      endif
c==============
cryne*** 3/16/2004
      if(slfrac.ne.1.0)p1=p1*slfrac
cryne***
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c
c 'prot    ': rotation of reference plane
c
105   continue
cryne Aug 5, 2003: now this is the 5th order version:
      ijkind=nint(p2)
      call prot(p1,ijkind,h,mh)
      goto 2000
cryne 9/25/2007      return
c
c 'gbdy    ': body of a general bending magnet
c
106   continue
      rho=brho/p4
      call gbend(rho,p1,p2,p3,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c 'frng    ': hard edge dipole fringe fields
c
107   continue
      rho=brho/p4
      iedge=nint(p5)
      gap = p2
      xk1 = p3
      call gfrngg(p1,rho,iedge,h,mh,gap,xk1)
      goto 2000
c
c 'cfbdy    ': body of a combined function bending magnet
c
1075  continue
      write(6,*)'cfbdy not implemented'
      return
c
c 'cfbd    ': combined function bend (normal entry and exit)
c
c The map for the combined function bend is of the form
c [(arot*frquad*arotinv)frquad*gfrngg]*
c cfbend*
c [gfrngg*frquad*(arot*frquad*arotinv)]
c with the leading and trailing fringe fields optional.
c The gap size cfbgap and normalized leading and trailing
c field integrals cfblk1 and cfbtk1 for the dipole
c are taken from the table in block common frnt.
c These entries are initialized to 0, .5, .5, respectively
c at the beginning of a Marylie run, and can be changed using
c the type code cfrn.
c The factors of the form (arot*frquad*arotinv) give skew
c quad fringe fields.
c
108   continue
      if(jsltot.ne.1)then
        write(6,*)'error: cfbd not yet available with slices!'
        stop
      endif
      rho=brho/p2
      lfrn=nint(p3)
      tfrn=nint(p4)
      ijopt=nint(p5)
      iopt=mod(ijopt,10)
c      write(6,*) 'iopt=',iopt
      if ((iopt.lt.1) .or. (iopt.gt.3)) then
      write(jof,*) 'WARNING: parameter iopt outside allowed range in',  &
     & ' element with type code cfbd'
      endif
      ipset=nint(p6)
      do j=1,6
      coefpset(j)=0.d0
      enddo
      if( (ipset.gt.0).and.(ipset.le.maxpst))then
        do j=1,6
        coefpset(j)=pst(j,ipset)
        enddo
      endif
c compute cfbend
      call cfbend(p,coefpset,h,mh)
c put on leading dipole and quad fringe fields
      if(lfrn.ne.0) then
c compute and put on leading dipole fringe field
      gap=cfbgap
      xk1=cfblk1
      call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
c      call nfrng(rho,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c compute and put on leading quad fringe fields
      if((ipset.gt.0) .and. (ipset.le.maxpst)) then
      if(iopt.eq.1) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.2) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.3) then
      bqd=brho*coefpset(1)
      aqd=brho*coefpset(2)
      endif
c compute and put on normal quad fringe field
      call frquad(bqd,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c compute and put on skew quad fringe field
      angr=-pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      call frquad(aqd,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      endif
      endif
c put on trailing dipole and quad fringe fields
      if(tfrn.ne.0) then
c compute and put on trailing dipole fringe field
      gap=cfbgap
      xk1=cfbtk1
      call gfrngg(0.d0,rho,2,ht1,mht1,gap,xk1)
cryne 7/11/2002      call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
      write(6,*)'7/11/2002 fixed 2nd call to gfrngg for cfbd.  RDR'
c      call nfrng(rho,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c compute and put on trailing quad fringe fields
      if((ipset.gt.0) .and. (ipset.le.maxpst)) then
      if(iopt.eq.1) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.2) then
      bqd=coefpset(1)
      aqd=coefpset(2)
      endif
      if(iopt.eq.3) then
      bqd=brho*coefpset(1)
      aqd=brho*coefpset(2)
      endif
c compute and put on normal quad fringe field
      call frquad(bqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c compute and put on skew quad fringe field
      angr=pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      call frquad(aqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      endif
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
c 'quad    ': quadrupole
c
c the map for this element is of the form
c frquad*quad*frquad
c with the leading and trailing fringe field maps optional
c
109   continue
c     write(6,*)'QUAD:jslice,jsltot,slfrac=',jslice,jsltot,slfrac
      lfrn=nint(p3)
      tfrn=nint(p4)
      if(slfrac.ne.1.0)p1=p1*slfrac
c compute the map for the quad
      if(p2.lt.0.d0) call dquad3(p1,-p2,h,mh)
      if(p2.eq.0.d0) call drift3(p1,h,mh)
      if(p2.gt.0.d0) call fquad3(p1,p2,h,mh)
c put on leading fringe field
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
        if(lfrn.ne.0)then
        call frquad(p2,-1,ht1,mht1)
        call concat(ht1,mht1,h,mh,h,mh)
        endif
      endif
c put on trailing fringe field
      if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
        if(tfrn.ne.0)then
        call frquad(p2,1,ht1,mht1)
        call concat(h,mh,ht1,mht1,h,mh)
        endif
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'sext    ': sextupole
c
110   continue
      if(slfrac.ne.1.0)p1=p1*slfrac
      call sext3(p1,p2,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'octm    ': mag. octupole
c
111   continue
      if(slfrac.ne.1.0)p1=p1*slfrac
      call octm3(p1,p2,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'octe    ': elec. octupole
c
112   continue
      if(jsltot.ne.1)then
        write(6,*)'error: octe not yet available with slices!'
        stop
      endif
      call octe(p1,p2,h,mh)
      goto 2000
c
c 'srfc    ': short rf cavity
c
113   omega=twopi*p2
      call srfc(p1,omega,h,mh)
      goto 2000
c
c 'rfgap  ': rf gap
c
1135  continue
c     write(6,*)'RFGAP:jslice,jsltot,slfrac=',jslice,jsltot,slfrac
c----------------------------------------------
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
        nseq=prms(5)
        if(nseq.eq.0 .and. crms(1).eq.' ')then
          if(idproc.eq.0)write(6,*)'ERROR(RFGAP):no data file specified'
          call myexit
        endif
        if(crms(1).eq.' ')then
          fname1='rfdata'
          ndigits=nseq/10+1
          call num2string(nseq,aseq,ndigits)
          j=len_trim(fname1)
          fname=fname1(1:j)//aseq(1:ndigits)
        else
          fname=crms(1)
        endif
        nunit=0
        estrng='rfgap'
        call fnamechk(fname,nunit,ierr,estrng)
        if(ierr.eq.1)then
          write(6,*)'(rfgap)leaving lmnt due to problem w/ fname'
          return
        endif
        call read_RFdata(nunit,numdata,gaplen)
c       Set zedge, which defines where the field is located in space:
        zedge=arclen
cryne 5/2/2006:
c determine theta for this cavity now, while jslice=1
        if(prms(7).eq.0.d0)then
          thetastore=p(4)*asin(1.0d0)/90.
          if(idproc.eq.0)then
          write(6,*)'found rfgap absolute phase (deg) =',prms(4)
          write(6,*)'found rfgap absolute phase (rad) =',thetastore
          endif
        else
          thetastore=prms(7)*asin(1.0d0)/90.d0 - reftraj(5)
          if(idproc.eq.0)then
         write(6,*)'found rfgap phase entrance argument (deg) =',prms(7)
          write(6,*)'=',prms(7)*asin(1.0d0)/90.d0,' rad'
          write(6,*)'computed rfgap absolute phase (rad)=',thetastore
          write(6,*)'=',thetastore*90./asin(1.0d0),' deg'
          endif
        endif
      endif
c----------------------------------------------
cryne 8/12/2001 did this in trlmnt      zedge=arclen
cryne 8/12/2001      zmap=p(1)
cryne 8/12/2001 also added the following 2 lines:
      zedgg=zedge
      if(p1.le.0)zmap=gaplen*slfrac
      if(p1.gt.0)then
       if(p1.lt.gaplen)zmap=p1*slfrac
       if(p1.gt.gaplen)write(6,*)'error:integration length > gap length'
       if(p1.gt.gaplen)stop
      endif
      rffreq=p(2)*4.0*asin(1.0d0)
      escale=p(3)
cryne 5/2/2006      theta0=p(4)*asin(1.0d0)/90.
      theta0=thetastore
      itfile=nint(p(5))
      nstep=nint(p(6))
cryne 8/11/2001 itype not used for now
      itype=0
cryne-abell Mon Nov 24 18:23:05 PST 2003
      if(lflagmagu)then
       if(idproc.eq.0)then
        write(6,*)'ERROR: The code is being run with static units, but'
        write(6,*)'an rf cavity has been encountered. Simulations that'
        write(6,*)'include rf cavities must utilize dynamic units.'
        write(6,*)'Re-run using dynamic units.'
       endif
       call myexit
      endif
cryne 5/1/2006 useful constants:
      gscal=-omegascl*sl*p0sc/pmass
      tscal= omegascl*sl/c
cryne 5/1/2006 for future mods, save reference phase and energy at entry:
      refentry5=reftraj(5)      !common block variable used in trac
      refentry6=reftraj(6)      !common block variable used in trac
cryne==========================
cryne 5/5/2006 new code for rfgap multiple reference trajectories
      multitrac=0
      imap5=prms(9)
      jmap6=prms(10)
      if(imap5.ne.1 .or. jmap6.ne.1)multitrac=1 !flag for multitrack
      if(multitrac.eq.1)then
        call multirfsetup(imap5,jmap6)
        do i=1,imap5
        do j=1,jmap6
          t00=tlistin(i)
          gam00=ptlistin(j)*gscal
          sltemp=sl
          sl=c/omegascl
          call rfgap(zedgg,zmap,nstep,rffreq,escale,theta0,t00,gam00,         &
     &               itype,h,mh)
          sl=sltemp
c Ensure rf cavity map is in correct units for ML/I. Then store it.
          call set_rfscale(h,mh)
          tmhlist(1:6,1:6,i,j)=mh(1:6,1:6)
          hlist(1:923,i,j)=h(1:923)  ! h irrelevant for rfgap (rfgap is linear)
c Store the final values of t and pt:
         tlistfin(i,j)=t00
         ptlistfin(i,j)=gam00/gscal
        enddo
        enddo
      endif
cryne==========================
c If multiple maps are being used, that have now all been computed.
c But, regardless, ML/I needs to compute the evolution of the reference traj:
c
cryne 5/1/2006 rfgap code was originally written assuming omegascl*sl/c=1.
c              This is not necessarily the case for ML/I.
c              Temporarily set sl to satisfy this, then reset upon return.
      t00=refentry5
      gam00=refentry6*gscal
      sltemp=sl
      sl=c/omegascl
      call rfgap(zedgg,zmap,nstep,rffreq,escale,theta0,t00,gam00,itype, &
     &h,mh)
      sl=sltemp
c--
cryne 5/1/2006 moved these from rfgap.f ; The philosophy is that rfgap does not
c               make global changes, it just computes a map and the ref traj
      gamma=gam00
      gamm1=gamma-1.0
      beta=sqrt((gamma+1.0)*(gamma-1.0))/gamma
      brho=pmass*gamma*beta/c
c--
      call set_rfscale(h,mh)
      refprev=reftraj
      reftraj(5)=t00
      reftraj(6)=gam00/gscal
      arclen=arclen+zmap
      prevlen=zmap
      nt2prev=nt2
      goto 2000
c
c 'confoc  ': "constant focusing" element   rdr 08/29/2001
c
1136  continue
c p1 is the length, p2,p3,p4 are the 3 focusing constants
c     if(idproc.eq.0)write(6,*)'at confoc with slfrac=',slfrac
      if(slfrac.ne.1.0)p1=p1*slfrac
      call confoc(p1,p2,p3,p4,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
cryne 5/1/2006      reftraj(5)=reftraj(5)+p1
c==============
      goto 2000
c
c 'arot    ': axial rotation
c
114   p1rad=pi180*p1
      call arot(p1rad,h,mh)
      goto 2000
c
c 'twsm    ': linear matrix via twiss parameters
c
115   iplane=nint(p1)
      call twsm(iplane,p2,p3,p4,h,mh)
      goto 2000
c
c 'thlm    ': thin lens low order multipole
c
116   continue
      call thnl(p1,p2,p3,p4,p5,p6,h,mh)
      goto 2000
c
c 'cplm    ': "compressed" low order multipole
c
117   call cplm (p,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'cfqd    ': combined function quadrupole
c
118   continue
      if(jsltot.ne.1)then
        write(6,*)'error: cfqd not yet available with slices!'
        stop
      endif
c     call cfdrvr(p,h,mh)
c replaced cfdrvr with cfqd in new (this) afro. rdr and ctm 5/22/02
c*****
      call cfqd(p,h,mh)
c quad fringe fields:
      eangle=0.d0  !not needed
      rho=1.d0     !not needed
      gap=0.d0     !not needed
      fint=0.d0    !not needed
      iopt=1  !standard MaryLie interpretation of multipole coeffs
      kfrn=2  !no dipole fringe fields (1=dipole, 2=quad, 3=both)
      ipset=nint(p(2))
      coefpset(1:6)=pst(1:6,ipset)
c leading fringe:
      ilfrn=nint(p(3))
      if(ilfrn.ne.0)then
c      if(idproc.eq.0)write(6,*)'adding leading fringe to cfqd'
       iw=1
       call cgfrngg(iw,eangle,rho,kfrn,gap,fint,iopt,coefpset,ht1,mht1)
       call concat(ht1,mht1,h,mh,h,mh)
      endif
c trailing fringe:
      itfrn=nint(p(4))
      if(itfrn.ne.0)then
c      if(idproc.eq.0)write(6,*)'adding trailing fringe to cfqd'
       iw=2
       call cgfrngg(iw,eangle,rho,kfrn,gap,fint,iopt,coefpset,ht1,mht1)
       call concat(h,mh,ht1,mht1,h,mh)
      endif
c*****
c
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
c dispersion matrix (dism)
c
119   call dism(p,h,mh)
      go to 2000
c
c 'sol     ': solenoid
c
120   continue
      write(6,*) ' == afro::lmnt::sol =='
      write(6,*) ' jslice=',jslice
      write(6,*) ' jsltot=',jsltot
      write(6,*) ' slfrac=',slfrac
      if(jsltot.ne.1)then
        write(6,*)'error: sol not yet available with slices!'
        stop
      endif
      call gensol(p,h,mh,jslice,jsltot,slfrac)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p2-p1
      prevlen=p2-p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+(p2-p1)/(beta*c)*(omegascl)
c==============
      go to 2000
c
c 'mark    ': MaryLie marker
c
121   continue
      return
c
c 'jmap    ': j mapping
c
122   call jmap(h,mh)
      go to 2000
c
c 'dp      ': data point
c
123   continue
      return
c
c  REC multiplet (recm)
c
124   continue
      if(jsltot.ne.1)then
        write(6,*)'error: recm not yet available with slices!'
        stop
      endif
      call gnrec3(p,h,mh)
      go to 2000
c
c 'spce    ': space
c
125   continue
      return
c
c 'cfrn    ': change or write out fringe field parameters
c for combined function dipole
c
126   continue
      mode=nint(p1)
      if (mode .eq. 0) then
c change fringe field parameters
      cfbgap=p2
      cfblk1=p3
      cfbtk1=p4
      return
      endif
c write out values of fringe field parameters
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of fringe field parameters gap, ennfi, exnfi are:'
      write(jof,*) cfbgap,cfblk1,cfbtk1
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of fringe field parameters gap, ennfi, exnfi are:'
      write(jodf,*) cfbgap,cfblk1,cfbtk1
      endif
      return
c
c 'coil'
c
127   continue
      if(jsltot.ne.1)then
        write(6,*)'error: coil not yet available with slices!'
        stop
      endif
      call coil(prms)
      return
c
c 'intg'
c
128   continue
      if(jsltot.ne.1)then
        write(6,*)'error: intg not yet available with slices!'
        stop
      endif
      call integ(p,h,mh)
      go to 2000
c
c 'rmap'
c
129   continue
      call rmap(p,h,mh)
      go to 2000
c
c 'arc'
c
130   continue
      write(6,*)'(lmnt) arc: not implemented'
      return
c
c 'transit'
133   continue
      write(6,*)'(lmnt) transit'
      write(6,*)'p1,p2=',p1,p2
      if(slfrac.ne.1.0)p1=p1*slfrac
      call transit(p1,p2,h,mh)
      arclen=arclen+p1
      goto 2000
c 'interface'
134   continue
      write(6,*)'(lmnt) interface'
      write(6,*)'p1,p2,p3,p4,p5=',p1,p2,p3,p4,p5
      call interface(p1,p2,p3,p4,p5,h,mh)
      goto 2000
c 'rootmap'
135   continue
      write(6,*)'(lmnt) rootmap'
      write(6,*)'p1,p2,p3,p4=',p1,p2,p3,p4
      write(6,*)'crms(1)=',crms(1)
      call rootmap(p1,p2,p3,p4,h,mh)
      if(crms(1).eq.'true')then
        write(6,*)'computing inverse of rootmap'
        call inv(h,mh)
      endif
      goto 2000
c
c 'optiprot ': rotation of reference plane for optics calculations
136   continue
      write(6,*)'(lmnt) optirot'
      ijkind=nint(p2)
      call optirot(p1,ijkind,h,mh)
      goto 2000
c
c 'spare5'
137   continue
      write(6,*)'(lmnt) spare5'
      return
c 'spare6'
138   continue
      write(6,*)'(lmnt) spare6'
      return
c 'spare7'
139   continue
      write(6,*)'(lmnt) spare7'
      return
c 'spare8'
140   continue
      write(6,*)'(lmnt) spare8'
      return
c
c
141   continue
c 'marker': MAD marker
      write(6,*)'MAD marker not implemented'
      return
c
142   continue
c 'drift': MAD drift
      if(slfrac.ne.1.0)p1=p1*slfrac
c     write(6,*)'DRIFT:jslice,jsltot,slfrac=',jslice,jsltot,slfrac
      isssflag=nint(p2)
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
      goto 2000
c
143   continue
c 'rbend': MAD rectangular bend
c     write(6,*)'computing a MAD rbend'
c     write(6,*)'computing a MAD rbend; jslice,jsltot=',jslice,jsltot
c based on a MaryLie gbnd:
c     if(jsltot.ne.1 .and. idproc.eq.0)then
c       write(6,*)'WARNING: rbend slices being tested'
c     endif
      if(p2.eq.0.d0)then
        write(6,*)'Error (gbnd): B field equal zero not allowed'
        stop
      endif
      rho=brho/p2
!
cryne Jan 3, 2005
c p3 and p4 and the pole face rotation angles in MAD notation, e1 and e2:
      e1=p(3)
      e2=p(4)
cryne Jan 3, 2005
c The arguments to gbend, gfrngg, and prot are angles between the
c pole face and the reference trajectory , *not* angles between the
c pole face and reference plane (vertical line for rbends, radial line
c for sbends). The latter is the convention used in MAD.
c These angles need to be checked!!! At this point the RBEND
c implemented here is highly dubious!!!
      bndang=p1
      psideg=0.5d0*bndang+e1
      phideg=0.5d0*bndang+e2
      if( abs(psideg-phideg).gt.1.d-6 )then
        write(6,*)'ERROR (RBEND): The present RBEND impelmentation'
        write(6,*)'requires equal entry and exit angles.'
        call myexit
      endif
c
      lfrn=nint(p5)
      tfrn=nint(p6)
      gap1=prms(7)
      gap2=prms(8)
      fint1=prms(9)
      fint2=prms(10)
      tiltang=prms(11)
      if(tiltang.ne.0. .and. idproc.eq.0)then
        write(6,*)'WARNING(rbend): tilt keyword being tested'
      endif
c***********************************************************
c***********************************************************
      if(jsltot.eq.1 .and. ihalf.eq.0)then
c procedure for a complete pbnd (no slices, no splitting)
        call pbend(rho,p1,h,mh)
c put on the leading fringe field
        call gfrngg(psideg,rho,1,ht1,mht1,gap1,fint1)
        call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
        call prot(psideg,1,ht1,mht1)
        call concat(ht1,mht1,h,mh,h,mh)
c       put on leading arot
        call arot(tiltang,ht1,mht1)
        call concat(ht1,mht1,h,mh,h,mh)
c put on trailing fringe field
        call gfrngg(psideg,rho,2,ht1,mht1,gap2,fint2)
        call concat(h,mh,ht1,mht1,h,mh)
c put on trailing prot
        call prot(psideg,2,ht1,mht1)
        call concat(h,mh,ht1,mht1,h,mh)
c       put on trailing arot
        tiltang=-tiltang
        call arot(tiltang,ht1,mht1)
        call concat(h,mh,ht1,mht1,h,mh)
      else
c
c procedure for a sliced and/or split pbnd
c This follows the approach due to Dragt that is described
c in Exhibit 6.6.3 and section 10.9 of the MaryLie manual.
        if(slfrac.ne.1.0)p1=p1*slfrac
        azero=0.
        call nbend(rho,p1,h,mh)
        if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
c         put on the leading gbody
          call gbend(rho,azero,psideg,azero,ht1,mht1)
          call concat(ht1,mht1,h,mh,h,mh)
c         put on the leading fringe field
          call gfrngg(psideg,rho,1,ht1,mht1,gap1,fint1)
          call concat(ht1,mht1,h,mh,h,mh)
c         put on leading prot
          call prot(psideg,1,ht1,mht1)
          call concat(ht1,mht1,h,mh,h,mh)
c         put on leading arot
          call arot(tiltang,ht1,mht1)
          call concat(ht1,mht1,h,mh,h,mh)
        endif
        if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
c         put on trailing gbody
          call gbend(rho,azero,azero,psideg,ht1,mht1)
          call concat(h,mh,ht1,mht1,h,mh)
c         put on trailing fringe field
          call gfrngg(psideg,rho,2,ht1,mht1,gap2,fint2)
          call concat(h,mh,ht1,mht1,h,mh)
c         put on trailing prot
          call prot(psideg,2,ht1,mht1)
          call concat(h,mh,ht1,mht1,h,mh)
c         put on trailing arot
          tiltang=-tiltang
          call arot(tiltang,ht1,mht1)
          call concat(h,mh,ht1,mht1,h,mh)
        endif
      endif
c***********************************************************
c***********************************************************
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+rho*p1*(pi/180.)
      reftraj(5)=reftraj(5)+p1*rho*(pi/180.)/(beta*c)*(omegascl)
      prevlen=rho*p1*(pi/180.)
      nt2prev=nt2
      rhoprev=rho
c==============
      goto 2000
c
145   continue
c 'gbend ': MAD general bend
      write(6,*)'MAD gbend not implemented'
      return
c
146   continue
c 'quadrupo': MAD quadrupole
c 7/14/2002: need to modify 'quadrupo' to handle tilt
c     write(6,*)'QUADRUPOLE:jslice,jsltot,slfrac=',jslice,jsltot,slfrac
      lfrn=nint(p3)
      tfrn=nint(p4)
      if(slfrac.ne.1.0)p1=p1*slfrac
c compute the map for the quad
      isssflag=nint(p5)
!     write(6,*)'QUAD: isssflag=',isssflag
      if(isssflag.eq.1)then
        if(p2.lt.0.d0) call dquad(p1,-p2,h,mh)
        if(p2.eq.0.d0) call drift(p1,h,mh)
        if(p2.gt.0.d0) call fquad(p1,p2,h,mh)
      else
        if(p2.lt.0.d0) call dquad3(p1,-p2,h,mh)
        if(p2.eq.0.d0) call drift3(p1,h,mh)
        if(p2.gt.0.d0) call fquad3(p1,p2,h,mh)
      endif
c put on leading fringe field
      if(jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
        if(lfrn.ne.0)then
        call frquad(p2,-1,ht1,mht1)
        call concat(ht1,mht1,h,mh,h,mh)
        endif
      endif
c put on trailing fringe field
      if(jslice.eq.jsltot.and.(ihalf.eq.0.or.ihalf.eq.2))then
        if(tfrn.ne.0)then
        call frquad(p2,1,ht1,mht1)
        call concat(h,mh,ht1,mht1,h,mh)
        endif
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
147   continue
c 'sextupol': MAD sextupole
c 7/14/2002: need to modify 'sextupol' to handle tilt
      goto 110
c
148   continue
c 'octupole': MAD octupole
c 7/14/2002: need to modify 'octupole' to handle tilt
      goto 111
c
149   continue
c 'multipol': MAD general thin multipole
      if(imsg1.eq.0)then
      imsg1=1
c     write(6,*)'MULTIPOLE MULTIPLIER CORRECT?'
      endif
      p1b=p1*brho
      p2b=p2*brho
      p3b=p3*brho/2.d0
      p4b=p4*brho/2.d0
      p5b=p5*brho/6.d0
      p6b=p6*brho/6.d0
      call thnl(p1b,p2b,p3b,p4b,p5b,p6b,h,mh)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+0.d0
      prevlen=0.d0
      nt2prev=nt2
      reftraj(5)=reftraj(5)+0.d0/(beta*c)*(omegascl)
c==============
      goto 2000
c
c 'solenoid': MAD solenoid
c
150   continue
      write(6,*) ' == afro::lmnt::solenoid =='
      write(6,*) ' jslice=',jslice
      write(6,*) ' jsltot=',jsltot
      write(6,*) ' slfrac=',slfrac
      call gensol(p,h,mh,jslice,jsltot,slfrac)
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      prevlen=(p2-p1)*slfrac
      arclen=arclen+prevlen
      nt2prev=nt2
      reftraj(5)=reftraj(5)+prevlen/(beta*c)*(omegascl)
c==============
      go to 2000
      return
c
151   continue
c 'hkicker': MAD hkicker
      if(p2.ne.0.d0)then
       write(6,*)'error: hkicker strength, kick = ',p2
       write(6,*)'but hkicker currently implemented only w/ kick=0.'
       write(6,*)'kick will be ignored (i.e. treat like a drift)'
      endif
c
      isssflag=nint(p3)
!     write(6,*)'HKICKER: isssflag=',isssflag
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
152   continue
c 'vkicker': MAD vkicker
      if(p2.ne.0.d0)then
       write(6,*)'error: vkicker strength, kick = ',p2
       write(6,*)'but vkicker currently implemented only w/ kick=0.'
       write(6,*)'kick will be ignored (i.e. treat like a drift)'
      endif
c
      isssflag=nint(p3)
!     write(6,*)'VKICKER: isssflag=',isssflag
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
153   continue
c 'kicker': MAD kicker
      if(p2.ne.0.d0)then
       write(6,*)'error: kicker strength hkick = ',p2
       write(6,*)'but kicker currently implemented only w/ hkick=0.'
       write(6,*)'hkick will be ignored (i.e. treat like a drift)'
      endif
      if(p3.ne.0.d0)then
       write(6,*)'error: kicker strength vkick = ',p3
       write(6,*)'but kicker currently implemented only w/ vkick=0.'
       write(6,*)'vkick will be ignored (i.e. treat like a drift)'
      endif
      isssflag=nint(p4)
!     write(6,*)'KICKER: isssflag=',isssflag
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
154   continue
c 'rfcavity': MAD RF cavity
      write(6,*)'MAD RF cavity commented out'
      zlen=p1
      volt=p2
      phlagrad=p3
      nharm=nint(p4)
!     call rfcavmad(zlen,volt,phlagrad,nharm,h,mh)
      call ident(h,mh)
!     call set_pscale_mc(h,mh)
c==============
!     call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
c
155   continue
c 'elsepara': MAD electrostatic separator
      zlen=p1
      efield=p2
c tilt is parameter 3; not used for now.
      isssflag=nint(p4)
      if(efield.ne.0.d0)then
        write(6,*)'MAD elec separator: zlen,efield=',zlen,efield
        call elsepmad(zlen,efield,h,mh)
      else
        if(isssflag.eq.1)then
          write(6,*)'MAD elsep has zero field; treat as 5th order drift'
          call drift(p1,h,mh)
        else
          write(6,*)'MAD elsep has zero field; treat as 3rd order drift'
          call drift3(p1,h,mh)
        endif
      endif
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
cryne 12/27/02      return
c
156   continue
c 'hmonitor': MAD hmonitor
      zlen=p1
      if(zlen.eq.0.d0)return
      isssflag=nint(p2)
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
      call hmonitor
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
cryne 12/27/02      return
c
157   continue
c 'vmonitor': MAD vmonitor
      zlen=p1
      if(zlen.eq.0.d0)return
      isssflag=nint(p2)
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
      call vmonitor
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
cryne 12/27/02      return
c
158   continue
c 'monitor': MAD monitor
      zlen=p1
      if(zlen.eq.0.d0)return
      isssflag=nint(p2)
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
      call monitor
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
cryne 12/27/02      return
c
159   continue
c 'instrument': MAD instrument
      zlen=p1
      if(zlen.eq.0.d0)return
      isssflag=nint(p2)
      if(isssflag.eq.1)then
        call drift(p1,h,mh)
      else
        call drift3(p1,h,mh)
      endif
      call instrument
c==============
      call set_pscale_mc(h,mh)
      refprev=reftraj
      arclen=arclen+p1
      prevlen=p1
      nt2prev=nt2
      reftraj(5)=reftraj(5)+p1/(beta*c)*(omegascl)
c==============
      goto 2000
cryne 12/27/02      return
c
160   continue
c 'sparem1'
      write(6,*)'(lmnt) sparem1'
      return
c
161   continue
c 'rcollima': MAD rectangular collimator
      write(6,*)'MAD rectangular collimator not implemented'
      return
c
162   continue
c 'ecollima': elliptic collimator
      write(6,*)'MAD elliptic collimator not implemented'
      return
c
163   continue
c 'yrot': MAD yrot
cryne 7/14/2002 same as MaryLie prot? Check this. (sign, etc.)
      goto 105
c
164   continue
c 'srot': MAD srot
cryne 7/14/2002 same as MaryLie arot? Check this. (sign, etc.)
      goto 114
c
165   continue
c 'prot3   ': (original 3rd order version of) rotation of reference plane
      ijkind=nint(p2)
      call prot3(p1,ijkind,h,mh)
      goto 2000
c
166   continue
c 'beambeam:' MAD beam-beam kick
      write(6,*)'MAD beambeam element not implemented'
      return
c
167   continue
c 'matrix': input a matrix using the MAD syntax
      write(6,*)'matrix input (MAD syntax) element not implemented'
      return
c
168   continue
c 'profile1d'
      ncol=nint(p1)
      nbins=nint(p2)
!p3 is the max number of file names in a sequence of output files
      nprecision=nint(p4)-1
      nunit=nint(p5)
      rwall=prms(7)
      if(idproc.eq.0)write(6,*)'(lmnt) profile1d, rwall=',rwall
!p6 is a counter that gets incremented and appended to a sequence of file names
      fname1=crms(1)
      fname=fname1
      if(p3.ne.0)then
        xdigits=log(p3)/log(10.d0)
        ndigits=1+int(xdigits)
        prms(6)=prms(6)+1.d0
        nseq=nint(prms(6)) !nseq is the counter that gets converted to a string
        call num2string(nseq,aseq,ndigits)
        j=len_trim(fname1)
        fname=fname1(1:j)//aseq(1:ndigits)
      endif
      if(nunit.eq.0)then
        estrng='profile1d'
        call fnamechk(fname,nunit,ierr,estrng)
        if(ierr.eq.1)then
          write(6,*)'(profile1d)leaving lmnt due to problem w/ fname'
          return
        endif
      endif
      call ibcast(nunit)
      if(nunit.gt.0)prms(5)=nunit
!     call pwritez(nunit,idmin,idmax,nprecision)
      if(ncol.gt.0)then
        call profile1d(ncol,nbins,rwall,arclen,nunit,fname,nprecision)
      else
        call profilerad(ncol,nbins,rwall,arclen,nunit,fname,nprecision)
      endif
      if(crms(2).eq.'true' .or. p3.ne.0)then
        if(idproc.eq.0)write(6,*)'closing file connected to unit ',nunit
        close(nunit)
        prms(5)=0.d0
      endif
      if(crms(3).eq.'true')then
!       if(idproc.eq.0)write(6,*)'flushing file unit ',nunit
        call myflush(nunit)
      endif
      if(idproc.eq.0)then
        write(6,*)'s=',arclen,' ;profile1d output to file ',fname
      endif
      return
c
169   continue
c 'yprofile'
      if(idproc.eq.0)write(6,*)'(lmnt) yprofile'
      return
c
170   continue
c 'tprofile'
      if(idproc.eq.0)write(6,*)'(lmnt) tprofile'
      return
c
171   continue
c 'hkick': synonym for MAD hkicker
      goto 151
c
172   continue
c 'vkick': synonym for MAD vkicker
      goto 152
c
173   continue
c 'kick': synonym for MAD kicker
      goto 153
c
174   continue
c 'hpm'
      write(6,*)'(hpm) half parallel faced magnet'
      write(6,*)'(hpm) b,phideg,iwhich=',p1,p2,p3
      rho=brho/p1
      phideg=p2
      iwhich=p3
      call hpf(rho,iwhich,phideg,h,mh)
      goto 2000
c
175   continue
c 'nlrf': non-linear RF cavity
      if (idproc.eq.0) then
        write(6,*) '(lmnt) nlrf: ntrk jslice jsltot slfrac ihalf = ',   &
     &              ntrk,jslice,jsltot,slfrac,ihalf
        write(6,*)'  ...under construction...'
      endif
      if (jslice.eq.1.and.(ihalf.eq.0.or.ihalf.eq.1))then
        ! for the first slice, read or compute generalized gradients
        ! (for now assume pre-computed generalized gradients)
        nunit=0
        estrng='nlrf'
        call fnamechk(crms(2),nunit,ierr,estrng)
        if (ierr.eq.1) then
          if (idproc.eq.0) then
            write(6,*) '<*** ERROR ***> (lmnt) nlrf: '
            write(6,*) 'leaving because of problems with file ',crms(2)
          endif
          call myexit()
        endif
        call read_egengrads(nunit)
        ! also record element length and z at entrance of element
        gaplen=eggrdata%zmax-eggrdata%zmin
        zedge=arclen
      endif
      if (idproc.eq.0) then
        write(6,*) 'zedge, gaplen =',zedge,gaplen; call myflush(6)
      end if
      ! for all slices:
      !   compute and check slice length
      zedgg=zedge
      zlen=p2-p1
      if (zlen.le.0) then
        zmap=gaplen*slfrac
      else
        if (zlen.le.gaplen) then
          zmap=zlen*slfrac
        else
          if (idproc.eq.0) then
            write(6,*) '<*** ERROR ***> (lmnt) nlrf:'
            write(6,*) '  Integration length exceeds gap length!'
          endif
          call myexit()
        endif
      endif
      !   note parameters
      zlc=arclen-zedge
      rffreq=twopi*p(3)
      rf_phase=p(4)*pi180
      rf_escale=p(5)
      nstep=nint(p(6))
      !nstep=eggrdata%nz_intrvl
      nslices=nint(p(7))
      if(lflagmagu)then
        if(idproc.eq.0)then
          write(6,*) '<*** ERROR ***> (lmnt) nlrf:'
          write(6,*) '  Simulations that contain RF cavities must use'
          write(6,*) '  dynamic units!  Re-run using dynamic units.'
        endif
        call myexit
      endif
      !   note initial time-of flight and initial gamma
      t00=reftraj(5)
      gam00=(-reftraj(6)*omegascl*sl*p0sc)/pmass
      if (idproc.eq.0) then
        write(6,*) 'calling subroutine nlrfcav():'
        write(6,*) '  zedge =',zedge
        write(6,*) '  zedgg, zlc, zmap, nstep =',zedgg,zlc,zmap,nstep
        write(6,*) '  rffreq =',rffreq
        write(6,*) '  t00, gam00 =',t00,gam00
        write(6,*) '  reftraj(5), reftraj(6) =',reftraj(5),reftraj(6)
      endif
      call nlrfcav(zedgg,zlc,zmap,nstep,t00,gam00,h,mh)
      if (idproc.eq.0) then
        write(6,*) 'returned from subroutine nlrfcav():'
        write(6,*) '  zedge =',zedge
        write(6,*) '  zedgg, zlc, zmap, nstep =',zedgg,zlc,zmap,nstep
        write(6,*) '  rffreq =',rffreq
        write(6,*) '  t00, gam00 =',t00,gam00
      endif
      refprev=reftraj
      !   rescale map
      !call set_rfscale(h,mh)
      !   update reftraj(5:6) and current arc length,
      reftraj(5)=t00
      reftraj(6)=-gam00*pmass/(omegascl*sl*p0sc)
      if (idproc.eq.0) then
        write(6,*) '  reftraj(5), reftraj(6) =',reftraj(5),reftraj(6)
      endif
      arclen=arclen+zmap
      !   and update length and index code of previous element
      prevlen=zmap
      nt2prev=nt2
      goto 2000
c
      return
c
c
c     2: user-supplied elements ************************************
c
c           'usr1    ','usr2    ','usr3    ','usr4    ','usr5    ',
12    go to (201,       202,       203,       204,       205,
c           'usr6    ','usr7    ','usr8    ','usr9    ','usr10   ',
     &       206,       207,       208,       209,       210,
c           'usr11   ','usr12   ','usr13   ','usr14   ','usr15   ',
     &       211,       212,       213,       214,       215,
c           'usr16   ','usr17   ','usr18   ','usr19   ','usr20   '
     &       216,       217,       218,       219,       220),kt2
c
201   call user1(p)
      return
202   call user2(p)
      return
203   call user3(p)
      return
204   call user4(p)
      return
205   call user5(p)
      return
206   call user6(p,th,tmh)
      return
207   call user7(p,th,tmh)
      return
208   call user8(p,th,tmh)
      return
209   call user9(p,th,tmh)
      return
210   call user10(p,th,tmh)
      return
211   call user11(p,th,tmh)
      return
212   call user12(p,th,tmh)
      return
213   call user13(p,th,tmh)
      return
214   call user14(p,th,tmh)
      return
215   call user15(p,th,tmh)
      return
216   call user16(p,th,tmh)
      return
217   call user17(p,th,tmh)
      return
218   call user18(p,th,tmh)
      return
219   call user19(p,th,tmh)
      return
220   call user20(p,th,tmh)
      return
c
c
c     3: parameter sets **********************************************
c
c           'ps1     ','ps2     ','ps3     ','ps4     ','ps5     ',
13    go to (301,       302,       303,       304,       305,           &
c           'ps6     ','ps7     ','ps8     ','ps9     '/
     &       306,       307,       308,       309),kt2
c
301   call pset(p,1)
      return
302   call pset(p,2)
      return
303   call pset(p,3)
      return
304   call pset(p,4)
      return
305   call pset(p,5)
      return
306   call pset(p,6)
      return
307   call pset(p,7)
      return
308   call pset(p,8)
      return
309   call pset(p,9)
      return
c
c     4-6: random elements
c
14    continue
15    continue
16    continue
      write(jof ,9014)
      write(jodf,9014)
 9014 format(' error in lmnt: random element reached')
      call myexit
c
c     7: simple commands *************************************
c
17    continue
c           'rt      ','sqr     ','symp    ','tmi     ','tmo     ',
      go to (701,       702,       703,       704,       705,
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
c           'poisson ','preapply','midapply','autoapply','autoconc'
     &       741,       742,       743,       744,       745,
c           'rayscale','beam    ','units   ','autoslic','verbose '
     &       746,       747,       748,       749,       750,
c           'mask6   ','arcreset','symbdef ','particledump','raytrace'
     &       751,       752,       753,       754,       755,
c           'autotrack','sckick','moments  ','maxsize','reftraj',
     &       756,       757,       758,       759,       760,
c           'initenv','envelopes','contractenv','setreftraj','setarclen',
     &       761,       762,       763,       764,       765,
c           'wakedefault','emittance','matchenv','fileinfo','egengrad',
     &       766,       767,       768,       769,       770,
c           'wrtmap  ','rdmap   ','sparec7','sparec8','sparec9'/
     &       771,       772,       773,       774,       775),kt2
c
c 'rt      ': ray trace
c
701   continue
      icfile=nint(p1)
      nfcfle=nint(p2)
      norder=nint(p3)
      ntrace=nint(p4)
      nwrite=nint(p5)
      fname1=crms(1)
      fname2=crms(2)
      iotemp=jof
      jfctmp=jfcf
      jfcf=nfcfle
cryne 3/27/04 ntaysym=1 for taylor, =2 for symplectic
cryne note: originally MaryLie code did symplectic if norder=5
      ntaysym=1
      if(norder.ge.5)ntaysym=2
      if(p6.lt.0.)jof=jodf
      ibrief=iabs(nint(p6))
      estrng='rt'
      if(icfile.gt.0)then
        if(fname1.ne.' ')then
          call fnamechk(fname1,icfile,ierr,estrng)
          if(ierr.eq.1)then
            write(6,*)'(rt) leaving lmnt due to problem w/ fname1'
            return
          endif
        endif
      endif
!bug fix on jan 3       if(jfcfile.gt.0)then
      if(nfcfle.gt.0)then
        if(fname2.ne.' ')then
          call fnamechk(fname2,nfcfle,ierr,estrng)
          if(ierr.eq.1)then
            write(6,*)'(rt) leaving lmnt due to problem w/ fname2'
            return
          endif
        endif
      endif
      call trace(icfile,jfcf,ntaysym,norder,ntrace,nwrite,0,0,5,0,         &
     &           th,tmh)
      jof=iotemp
      jfcf=jfctmp
      return
c
c     square the existing map:
c
702   call concat(th,tmh,th,tmh,th,tmh)
      write(jof,5702)
 5702 format(1h ,'existing map concatenated with itself')
      if(ntrk.eq.1)write(jof,9702)
 9702 format(1h ,'warning: map squared in turtle mode')
      return
c
c     symplectify matrix in transfer map
c
703   continue
      iopt=nint(p1)
      ijkind=nint(p2)
      call sympl(iopt,ijkind,th,tmh)
      return
c
c     input transfer map from an external file:
c
704   continue
      iopt=nint(p1)
      ifile=nint(p2)
      nopt=nint(p3)
      nskp=nint(p4)
c    rewind only option
      if (nopt.eq.1 .and. nskp.eq.-1) then
       rewind ifile
       write(jof,5704) ifile
 5704  format(1x,'file unit ',i3,' rewound')
       return
      endif
c
c    other options
c
      mpitmp=mpi
      mpi=ifile
      call mapin(nopt,nskp,h,mh)
      mpi=mpitmp
c  option when tracking (procedure at end)
      if(ntrk.eq.1) goto 2000
c  options when not tracking
      if(iopt.eq.1) call concat(th,tmh,h,mh,th,tmh)
      if(iopt.eq.2) call mapmap(h,mh,th,tmh)
      return
c
c     output transfer map to an external file (tmo):
c
705   ifile=nint(p1)
      mpotmp=mpo
      mpo=ifile
      call mapout(0,th,tmh)
      mpo=mpotmp
      return
c
c     print contents of file master input file:
c
706   continue
      itype=nint(p1)
      ifile=nint(p2)
      isend=nint(p3)
      fname=crms(1)
      if((ifile.eq.5).or.(ifile.eq.11).or.(ifile.eq.13))then
        write(6,*)'Error: cannot write to unit 5, 11, or 13 w/ PMIF.'
        write(6,*)'PMIF command will be ignored'
        return
      endif
      estrng='PMIF'
      if((ifile.ne.6).and.(fname.ne.' '))                               &
     &   call fnamechk(fname,ifile,ierr,estrng)
      if(ierr.eq.1)return
c
      jtmp=jodf
      jodf=ifile
      if(isend.eq.1.or.isend.eq.3)call pmif(jof,itype,fname)
      if(isend.eq.2.or.isend.eq.3)call pmif(jodf,itype,fname)
      jodf=jtmp
      return
c
c     Note: the program should never get here.
c     (cqlate is called directly from tran)
c
707   write(jof ,9707)
      write(jodf,9707)
 9707 format(' error: reached element "circ" in routine lmnt')
      call myexit
c
c     store the existing transfer map
c
708   continue
      kynd='stm'
      nmap=nint(p1)
      if ((nmap.gt.20).or.(nmap.lt.1)) then
        write(jof,9708) nmap
 9708   format(1x,'nmap=',i3,1x,'trouble with stm:nmap < 1 or > 20')
        call myexit
      else
        call strget(kynd,nmap,th,tmh)
        return
      endif
c
c     get transfer map from storage
c
709   continue
      kynd='gtm'
      iopt=nint(p1)
      nmap=nint(p2)
c
      if ((nmap.gt.20).or.(nmap.lt.1)) then
        write(jof,9709) nmap
 9709   format(1x,'nmap=',i3,1x,'trouble with gtm:nmap < 1 or > 20')
        call myexit
        return
      endif
c
      call strget(kynd,nmap,h,mh)
cryne--- 08/21/2001
      call mapmap(h,mh,hprev,mhprev)
c option when tracking(procedure at end)
      if(ntrk .eq. 1) goto 2000
c options when not tracking
      if(iopt.eq.1) call concat(th,tmh,h,mh,th,tmh)
      if(iopt.eq.2) call mapmap(h,mh,th,tmh)
      return
c
c     end of job:
c
710   continue
      if(idproc.eq.0)then
        write(jof,5710)
 5710   format(/1x,'end of MARYLIE run')
!       write(6,*)'total length=',arclen
      endif
      iprinttimers=0
      if(crms(1).eq.'true')iprinttimers=1
      call myexit
      return
c
c     print transfer map:
c
711   continue
      n1=nint(p1)
      n2=nint(p2)
      n3=nint(p3)
      n4=nint(p4)
      n5=nint(p5)
      if(n5.eq.1)then
        if(lautotrk.eq.0)then
          call pcmap(n1,n2,n3,n4,th,tmh)
        else
          call pcmap(n1,n2,n3,n4,hprev,mhprev)
        endif
      endif
      if(n5.eq.2)then
        if(lautotrk.eq.0)then
          call psrmap(n1,n2,th,tmh)
        else
          call psrmap(n1,n2,hprev,mhprev)
        endif
      endif
      if(n5.eq.3)then
        if(lautotrk.eq.0)then
          call pdrmap(n1,n2,th,tmh)
        else
          call pdrmap(n1,n2,hprev,mhprev)
        endif
      endif
      return
c
c     identity mapping:
c
712   continue
      call ident(th,tmh)
cryne 08/16/2001
      call ident(hprev,mhprev)
      return
c
c     write history of beam loss
c
713   call whst(p)
      goto 2000
c
c     inverse:
c
714   continue
      call inv(th,tmh)
      return
c
c     transpose:
c
715   call mtran(tmh)
      return
c
c     reverse factorization:
c
716   iord=nint(p1)
      call revf(iord,th,tmh)
      return
c
c     Dragt's reversal
c
717   call rev(th,tmh)
      return
c
c     mask off selected portions of transfer map:
c
718   call mask(p,th,tmh)
      return
c
c     number lines in a file
c
719   call numfile(p)
      return
c
c     aperture particle distribution
c
720   call rapt(p)
      return
c
721   continue
      mode=nint(p1)
      if (mode .eq. 1) call eapt(p)
      return
c
c     open files
c
722   call of(p)
      return
c
c     close files
c
723   call cf(p)
      return
c
c     window particle distribution
c
724   call wnd(p)
      return
725   call wnda(p)
      return
c
c     filter transfer map
c
726   call ftm(p,th,tmh)
      return
c
c     write parameter set
c
727   ipset = nint(p(1))
      isend = nint(p(2))
      call wps(ipset,isend)
      return
c
c     write time
c
728   call mytime(p)
      return
c
c     change output drop file
c
729   jodf = nint(p(1))
      return
c
c     ring bell
c
730   continue
cryne call bell
      if(idproc.eq.0)write(6,*)'BELL BELL BELL BELL BELL BELL BELL '
      return
c
c     write value of merit function
c
731   ifn = nint(p(1))
      isend = nint(p(2))
      call wmrt(ifn,isend)
      return
c
c     write contents of loop
c
732   call wcl(p)
      return
c
c     pause (paws)
c
733   continue
      write(jof,*) ' press return to continue'
      read(5,7330) iwxyz
7330  format(a1)
      return
c
c     change or write out infinities (inf)
c
734   continue
      mode=nint(p1)
      if (mode .eq. 0) then
c change infinities
      xinf=p2
      yinf=p3
      tinf=p4
      ginf=p5
      return
      endif
c write out values of infinities
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of infinities xinf, yinf, tinf, ginf are:',
     & xinf,yinf,tinf,ginf
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of infinities xinf, yinf, tinf, ginf are:',
     & xinf,yinf,tinf,ginf
      endif
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
      mode=nint(p1)
      if (mode .eq. 0) then
c change values of zeroes
      fzer=p2
      detz=p3
      return
      endif
c write out values of zeroes
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of zeroes fzero, detzero are:',
     & fzer,detz
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of zeroes fzero, detzero are:',
     & fzer,detz
      endif
      return
c
c sndwch
c
737   continue
      n1=nint(p1)
      if(n1.eq.0)then
      call sndwch(hprev,mhprev,th,tmh,th,tmh)
      else
      call sndwchi(hprev,mhprev,th,tmh,th,tmh)
      endif
      return
c
c     twiss polynomial (tpol)
c
738   continue
      call tpol(p,th,tmh)
      return
c
c     dispersion polynomial (dpol)
c
739   continue
      call dpol(p,th,tmh)
      return
c
c     change or write out beam parameters (cbm)
c
740   continue
      job=nint(p1)
      if (job .eq. -1) then
      prms(1)=0.d0
      prms(2)=brho
      prms(3)=gamm1
      endif
      if (job .eq. 0) then
      brho=p2
      gamm1=p3
c  recomputation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      endif
      if ((job .eq. 1) .or. (job .eq. 3)) then
      write(jof,*) ' beam parameters are: ', brho,gamm1
      endif
      if ((job .eq. 2) .or. (job .eq. 3)) then
      write(jodf,*) ' beam parameters are: ', brho,gamm1
      endif
      return
c
c     POISSON
c     set parameters for poisson solver
c
741   continue
! Here are the real (array p) and character (array cparams) parameters:
! array p and cparam are dimensioned p(60) and cparams(60), but only
! the first few elements (12 and 17, resp) are used in this case.
! Somehow these need to be made available to spch3d.
! For now they are stored in common, then passed to spch3d through
! its argument list:
      rparams(1:nparams)=prms(1:nparams)
      cparams(1:ncparams)=crms(1:ncparams)
!
! this needs a lot of cleaning up. RDR 12 Nov 2003
!
!  p(1-3)=nx,ny,nz
!  p(4-5)=xmin,xmax
!  p(6-7)=ymin,ymax
!  p(8-9)=zmin,zmax
!  p(10-11)=anag_patchsize,anag_refineratio
!  cparam(1) : solver [i.e. solver type]
!  cparam(2) : geometry [not currently used]
!  cparam(3) : gridsize [fixed or variable]
!  cparam(4) : [determines whether # of gridpoints is fixed or variable]
!  cparam(5-7) : xboundary,yboundary,zboundary [open,dirichlet,periodic]
!  cparam(8) : boundary [not currently used; will specify x,y,z together]
!  cparam(9) : solving_for [determines whether solving for phi or E]
!  cparam(10): densityfunction ["delta" (old method); "linear" (new)]
!  cparam(11): chombo_input_file [default: 'undefined']
!  cparam(12): anag_smooth [default: 'none']
!  cparam(13-17)= 5 spares
!
      nxsave=nint(prms(1))
      nysave=nint(prms(2))
      nzsave=nint(prms(3))
c
      noresize=1
      if(crms(4).eq.'variable')noresize=0
      nadj0=0
      if(crms(7).eq.'periodic')nadj0=1
cryne--- Dec 29, 2002
      if(crms(5).eq.'dirichlet' .and.                                   &
     &   crms(6).eq.'dirichlet' .and.                                   &
     &   crms(7).eq.'dirichlet')then
          idirich=1
      else
          idirich=0
      endif
cryne---
      kfixbdy=0
      if(crms(3).eq.'fixed')kfixbdy=1
      xmin0=prms(4)
      xmax0=prms(5)
      ymin0=prms(6)
      ymax0=prms(7)
      zmin0=prms(8)
      zmax0=prms(9)
c Dec 29, 2002
c initialize the green function to "not made"
      madegr=0
c---------------------------------------
cryne Nov 12, 2003 Removed code dealing with ntrkorder and ntrktype
cryne because that has to do with dynamics, and this type code
cryne has to do with the Poisson solver.
cryne The dynamics info should be specified in the autotrack command
c---------------------------------------
cryne April 23, 2003
      if (crms(9).eq.'e'.or.crms(9).eq.'E') then
        idirectfieldcalc=1
      else if (crms(9).eq.'phi'.or.crms(9).eq.'Phi'                     &
     &         .or.crms(9).eq.'PHI') then
        idirectfieldcalc=0
      else
        if (idproc.eq.0) then
          write(6,*) ' <*** ERROR ***> unrecognized value ',crms(9)
          write(6,*) '   for parameter \"solving_for\" in POISSSON.'
        end if
        call myexit()
      end if
      if (crms(10).eq.'delta'.or.crms(9).eq.'Delta') then
        idensityfunction=0
      else if (crms(10).eq.'linear'.or.crms(10).eq.'Linear') then
        idensityfunction=1
      else
        if (idproc.eq.0) then
          write(6,*) ' <*** ERROR ***> unrecognized value ',crms(9)
          write(6,*) '   for parameter \"densityfunction\" in POISSSON.'
        end if
        call myexit()
      end if
      if (idproc.eq.0) then
        write(6,*) 'idirectfieldcalc=',idirectfieldcalc
        if (idirectfieldcalc.eq.0) then
          write(6,*) ' solving for the scalar potential phi'
        else
          write(6,*) ' solving for the electric field'
        end if
        write(6,*) 'idensityfunction=',idensityfunction
        if (idensityfunction.eq.0) then
          write(6,*) ' using delta \"interpolation\" for charges'
        else
          write(6,*) ' using linear interpolation charge distribution,'
          write(6,*) ' (i.e. using integrated Green function technique)'
        end if
      endif
c
c---------------------------------------
!dbs Nov03
      if(crms(1).eq.'fft')then
        isolve=1
      elseif(crms(1).eq.'fft2')then
        if(idirich.eq.0)then
          isolve=10  !alternate (ANAG) infinite domain solver
        else
          isolve=20  !alternate (ANAG) homogenous Dirichlet solver
        endif
      elseif(crms(1).eq.'chombo')then
        if(idirich.eq.0)then
          isolve=30  !Chombo infinite domain solver
        else
          isolve=40  !Chombo homogenous Dirichlet solver
        endif
      else
        if(idproc.eq.0)
     &    write(6,*) '(poisson) error: solver=' ,TRIM(crms(1))
     &              ,' is invalid.  Use fft, fft2 or chombo'
        call myexit
      endif

      ! set extra paramters for ANAG Poisson solver (James algorithm)
      if( isolve .NE. 1 )then
        ! parameters for James algorithm
        if( crms(5).eq.'open' .AND. crms(6).eq.'open' .AND.
     &      crms(7).eq.'open' )then
          ! anag_patchsize: size of grid blocks in FFT2 solver (must be multiple of 4)
          ! anag_refineratio: to create MLC global coarse grid from MLI grid
          anagpatchsize = prms(10)
          anagrefratio = prms(11)
        endif

        !mehrstellen smoothing on phi after MLC solve
        !anag_smooth=<anything except 'none' or 'off'>
        if( crms(12).eq.' ' .OR. crms(12).eq.'none' .OR.
     &      crms(12).eq.'off' )then
          anagsmooth = 0
        else
          anagsmooth = 1
        endif
      endif

      !chombo_file=<chombo parameters file> (see spch3d_chombo.f::CH_READINFILE())
      if( crms(11).eq.'undefined')then
        chombofilename = ' '
      else
        if( LEN_TRIM(crms(11)) .GT. LEN( chombofilename ) )then
          if(idproc.eq.0)
     &       write(6,*)'(poisson) error: chombo_input_file is too long'
          call myexit
        endif
        chombofilename = crms(11)
      endif

      if( idproc.eq.0 .and. iverbose.gt.3 )then
        write(6,*) 'lmnt: (poisson) input solver ',TRIM(crms(1))
     &            ,', isolve = ',isolve
      endif
!dbs
c---------------------------------------
c
c set the flag to indicate that s.c. parameters have been set:
cryne Nov 12, 2003: this was commented out Dec 2, 2002, but
cryne it belongs here so I am putting it back.
      nspchset=1
      if(idproc.eq.0)then
        write(6,*)'setting poisson parameters:'
        write(6,*)'nx,ny,nz=',nxsave,nysave,nzsave
        write(6,*)'noresize,nadj=',noresize,nadj0
cryne Nov 12, 2003: this really *is* wrong, so I am commenting out:
c       if(nspchset.eq.1)then
c         write(6,*)'autotracking w/ space charge from this point on'
c         write(6,*)'order=',ntrkorder
c       endif
      endif
cryne Dec 29, 2002: added module spchdata
        nx=nxsave
        ny=nysave
        nz=nzsave
        n1=2*nx
        n2=2*ny
        n3=2*nz
        nadj=nadj0
        n3a=2*nz-nadj*nz
cryne Dec 29, 2002: for Dirichlet case, double grid in all 3 dimensions
        if(idirich.eq.1)n3a=2*nz
        call new_spchdata(nx,ny,nz,n1,n2,n3a,isolve)
        if(idproc.eq.0)then
          write(6,*)'DONE ALLOCATING SPACE CHARGE ARRAYS'
        endif
      return
c
c     preapply commands automatically
c
742   continue
c     n1=nint(p1)
c     autostr(1)=cmenu(n1)
      autostr(1)=crms(1)
      if(iverbose.ge.1)write(6,*)'(lmnt)autostr(1)=',autostr(1)
      lautoap1=1
      if(autostr(1).eq.'off')then
        if(iverbose.ge.1)write(6,*)'turning off pre-autoapply'
        lautoap1=0
      endif
c restrict automatic application unless "applyto=all" :
      lrestrictauto=1
      if(crms(2).eq.'all')lrestrictauto=0
      return
c
c     midapply commands automatically
c
743   continue
c     n1=nint(p1)
c     autostr(2)=cmenu(n1)
      autostr(2)=crms(1)
      if(iverbose.ge.1)write(6,*)'(lmnt)autostr(2)=',autostr(2)
      lautoap2=1
      if(autostr(2).eq.'off')then
        if(iverbose.ge.1)write(6,*)'turning off mid-autoapply'
        lautoap2=0
      endif
c restrict automatic application unless "applyto=all" :
      lrestrictauto=1
      if(crms(2).eq.'all')lrestrictauto=0
      return
c
c     autoapply (actually, this is 'post-apply') commands automatically
c
744   continue
c     n1=nint(p1)
c     autostr(3)=cmenu(n1)
      autostr(3)=crms(1)
      if(iverbose.ge.1 .and. idproc.eq.0)                                &
     &   write(6,*)'(lmnt/posapply)autostr(3)=',autostr(3)
      lautoap3=1
      if(autostr(3).eq.'off')then
        if(iverbose.ge.1 .and. idproc.eq.0)                              &
     &  write(6,*)'turning off post-autoapply'
        lautoap3=0
      endif
c restrict automatic application unless "applyto=all" :
      lrestrictauto=1
      if(crms(2).eq.'all')lrestrictauto=0
      return
c
c     autoconc ('auto-concatenate')
c
745   continue
      if(crms(1).eq.'true')then
        lautocon=1
        lautotrk=0
        if(idproc.eq.0)write(6,*)'lautocon=1, lautotrk=0'
      elseif(crms(1).eq.'false')then
        lautocon=0
        lautotrk=1
        if(idproc.eq.0)write(6,*)'lautocon=0, lautotrk=1'
      elseif(crms(1).eq.'sandwich')then
        lautocon=2
        lautotrk=0
        if(idproc.eq.0)write(6,*)'lautocon=2, lautotrk=0'
        if(idproc.eq.0)write(6,*)'using sndwch instead of concat'
      elseif(crms(1).eq.'revsandwich')then
        lautocon=-2
        lautotrk=0
        if(idproc.eq.0)write(6,*)'lautocon=-2, lautotrk=0'
        if(idproc.eq.0)write(6,*)'using sndwchi instead of concat'
      else
        if(idproc.eq.0)then
        write(6,*)'error(autoconc): do not understand parameter'
        endif
        call myexit
      endif
      return
c
c     rayscale ('scale zblock array')
c
746   continue
c     if(idproc.eq.0.and.iverbose.gt.0)write(6,*)'scaling zblock array'
      if(idproc.eq.0)then
        write(6,*)'scaling zblock array; p(1)-p(6)='
        write(6,*)p(1),p(2)
        write(6,*)p(3),p(4)
        write(6,*)p(5),p(6)
      endif
      do 7462 j=1,nraysp
      do 7461 i=1,6
      zblock(i,j)=zblock(i,j)*p(i)
 7461 continue
 7462 continue

      return
c
c     spare1c
c
747   continue
      write(6,*)'(lmnt) BEAM type code codes here'
      return
c
c     spare2c
c
748   continue
      write(6,*)'(lmnt) UNITS type code codes here'
      return
c
c autoslice:
749   continue
      if(crms(1).ne.'local' .and. crms(1).ne.'global' .and.               &
     &   crms(1).ne.'none')then
        if(idproc.eq.0)then
        write(6,*)'Error(autoslice):'
        write(6,*)'Invalid slice control (local/global/none):',crms(1)
        write(6,*)'This autoslice command will be ignored'
        endif
        return
      endif
      if(prms(1).ne.0.d0 .and. prms(2).ne.0.d0)then
        if(idproc.eq.0)then
        write(6,*)'Error(autoslice): Cannot specify both'
        write(6,*)'# of slices and length between slices.'
        write(6,*)'Current values are:'
        write(6,*)'# of slices= ',prms(1)
        write(6,*)'length between slices= ',prms(2)
        write(6,*)'This autoslice command will be ignored'
        endif
        return
      endif
! input values are all valid; set autoslice parameters:
! default:
      slicetype='slices'
! set parameters:
      if(crms(1).eq.'none')then
        slicetype='none'
        slicevalue=0.d0
      endif
      if(prms(1).gt.0.d0)then
        slicetype='slices'
        slicevalue=prms(1)
      endif
      if(prms(2).gt.0.d0)then
        slicetype='interval'
        slicevalue=prms(2)
      endif
      sliceprecedence=crms(1)
c space charge kick:
c---------------------------------------
cryne Nov 14, 2003
cryne should probably get rid of this code eventually.
cryne does not make sense to specify whether or not to sckick
cryne when specifying the autoslice command.
cryne sckick info should be specified in the autotrack command,
cryne and in fact it should be enabled by default if autotracking
cryne is enabled, unless the user explicitly prevents it.
c     if(crms(2).eq.'true')nsckick=1
c inefficient, but works. fix later. rdr dec 9 2002
c     n1=index(crms(2),'1')
c     n2=index(crms(2),'2')
c     n3=index(crms(2),'3')
c     n4=index(crms(2),'4')
c     n5=index(crms(2),'5')
c     n6=index(crms(2),'6')
c     n7=index(crms(2),'7')
c     n8=index(crms(2),'8')
c     n9=index(crms(2),'9')
c     if(n1.ne.0)ntrkorder=1
c     if(n2.ne.0)ntrkorder=2
c     if(n3.ne.0)ntrkorder=3
c     if(n4.ne.0)ntrkorder=4
c     if(n5.ne.0)ntrkorder=5
c for now 6, 7, 8, or 9 are all the same as order 5:
c     if(n6.ne.0)ntrkorder=5
c     if(n7.ne.0)ntrkorder=5
c     if(n8.ne.0)ntrkorder=5
c     if(n9.ne.0)ntrkorder=5
c     ntay=index(crms(2),'tay')
c     nsym=index(crms(2),'sym')
c     if(ntay.ne.0)ntrktype=1
c     if(nsym.ne.0)ntrktype=2
c     if(n1+n2+n3+n4+n5+n6+n7+n8+n9+ntay+nsym.ne.0)nspchset=1
c---------------------------------------
      if(idproc.eq.0)then
      write(6,*)'autoslice;'
      write(6,*)'crms(1)=',crms(1)
      write(6,*)'crms(2)=',crms(2)
      write(6,*)'crms(3)=',crms(3)
      write(6,*)'slicetype=',slicetype
      write(6,*)'sliceprecedence=',sliceprecedence
      write(6,*)'slicevalue=',slicevalue
c     write(6,*)'nsckick=',nsckick
c     if(nspchset.ne.0)write(6,*)'ntrktype=',ntrktype
c     if(nspchset.ne.0)write(6,*)'ntrkorder=',ntrkorder
      endif
      return
c
c verbose:
750   continue
      iverbose=p(1)
      return
c
c mask6:
751   continue
      write(6,*)'(LMNT) Performing MASK6 command'
      call mask6(p,th,tmh)
      return
c
c arcreset:
752   continue
      if(idproc.eq.0)write(6,*)'resetting arc length to zero'
      arclen=0.d0
      return
c symbdef (="symbolic default"):
c this does not have any meaning in this routine; it only affects
c how the parser interprets in arithmentic expressions
c symbolic names that have not been define. Simply return from here.
753   continue
      return
c
c particledump:
754   continue
      idmin=nint(p1)
      idmax=nint(p2)
      if(idmax.eq.0)idmax=maxray
      nphysunits=nint(prms(7))
!p3 is the max number of file names in a sequence of output files
      nprecision=nint(p4)-1
      nunit=nint(p5)
!p6 is a counter that gets incremented and appended to a sequence of file names
      fname1=crms(1)
      fname=fname1
      if(p3.ne.0)then
        xdigits=log(p3)/log(10.d0)
        ndigits=1+int(xdigits)
        prms(6)=prms(6)+1.d0
        nseq=nint(prms(6)) !nseq is the counter that gets converted to a string
        call num2string(nseq,aseq,ndigits)
        j=len_trim(fname1)
        fname=fname1(1:j)//aseq(1:ndigits)
      endif
      estrng='particledump'
      if(nunit.eq.0)then
        call fnamechk(fname,nunit,ierr,estrng)
        if(ierr.eq.1)then
          write(6,*)'(particledump)leaving lmnt due to problem w/ fname'
          return
        endif
        call ibcast(nunit)
        prms(5)=nunit
      endif
      iprintarc=0
      if(crms(4).eq.'true')iprintarc=1
      call pwritez(nunit,idmin,idmax,nprecision,nphysunits,iprintarc)
      if(crms(2).eq.'true' .or. p3.ne.0)then
        if(idproc.eq.0)write(6,*)'closing file connected to unit ',nunit
        close(nunit)
        prms(5)=0.d0
      endif
      if(crms(3).eq.'true')then
!       if(idproc.eq.0)write(6,*)'flushing file unit ',nunit
        call myflush(nunit)
      endif
      if(idproc.eq.0)then
        write(6,*)'s=',arclen,' ;particledump to file ',fname
      endif
      return
c
c raytrace: (has more parameters than original rt command)
755   continue
      idmin=nint(p1)
      idmax=nint(p2)
      if(idmax.eq.0)idmax=maxray
c
      if(nint(p3).eq.0 .and. crms(5).eq.'undefined')then
c     Since the default value of p(3)=5, the user must have set it to zero.
c     In this case, he/she is probably just reading in some rays.
       if(idproc.eq.0)then
       write(6,*)'Error: you are using old style to specify raytrace'
       write(6,*)'order=0. Instead, if you only want to read particles,'
       write(6,*)'use type=readonly'
       endif
       call myexit
      endif
      if(crms(5).eq.'undefined')then
        if(idproc.eq.0)then
        write(6,*)'Error: you are using old style method to specify'
        write(6,*)'raytrace order.'
        write(6,*)'Instead use type=taylorN or type=symplecticN'
        endif
        call myexit
      endif
c
      if(crms(5).eq.'readonly' .or. crms(5).eq.'readwrite')then
        if(crms(5).eq.'readonly')  norder=-1
        if(crms(5).eq.'readwrite') norder=0
        ntaysym=1 !irrelevent in this case
c       note: when the user specifies readonly or readwrite, the default
c       should be ntrace=0,nwrite=0; this is taken care of in sif.f
      else
        call getordtyp(crms(5),ntaysym,norder)
      endif
      ntrace=nint(p4)
      nwrite=nint(p5)
      ibrief=iabs(nint(p6))
!prms(7) is the max number of file names in a sequence of output files
      nprecision=nint(prms(8))-1
      icunit=nint(prms(9))  !assigned unit number for input file
      nunit=nint(prms(10))  !assigned unit number for output file
      nraysinp=nint(prms(12)) !# of rays to read from data file (optional)
!p11 is a counter that gets incremented and appended to a sequence of file names
      icname=crms(1)
      fcname=crms(2)
      fname=fcname
c
      if(prms(7).ne.0)then
        xdigits=log(prms(7))/log(10.d0)
        ndigits=1+int(xdigits)
        prms(11)=prms(11)+1.d0
        nseq=nint(prms(11)) !nseq is the counter that gets converted to a string
        call num2string(nseq,aseq,ndigits)
        j=len_trim(fcname)
        fname=fcname(1:j)//aseq(1:ndigits)
      endif
c
      estrng='raytrace'
      if(icunit.eq.0 .and. icname.ne.' ')then
        call fnamechk(icname,icunit,ierr,estrng)
        if(ierr.eq.1)then
          if(idproc.eq.0)then
          write(6,*)'(rt) leaving lmnt due to problem w/ icunit'
          write(6,*)'icunit=',icunit
          write(6,*)'icname=',icname
          endif
          return
        endif
        call ibcast(icunit)
        prms(9)=icunit
      endif
!
      if(nunit.eq.0 .and. fcname.ne.' ')then
        call fnamechk(fname,nunit,ierr,estrng)
        if(ierr.eq.1)then
          write(6,*)'(rt) leaving lmnt due to problem w/ fcname'
          return
        endif
        call ibcast(nunit)
        prms(10)=nunit
      endif
!
      call trace(icunit,nunit,ntaysym,norder,ntrace,nwrite,idmin,idmax, &
     &           nprecision,nraysinp,th,tmh)
!close?
      if(icunit.ne.0)then
        if(idproc.eq.0)then
       write(6,*)'closing particle initial condition file connected to',&
     &  ' unit ',icunit
        endif
        close(icunit)
        prms(9)=0.d0
      endif
      if(crms(3).eq.'true' .or. prms(7).ne.0)then
        write(6,*)'closing particle final condition file connected to', &
     &  ' unit ',nunit
        close(nunit)
        prms(10)=0.d0
      endif
!flush?
      if(crms(4).eq.'true')then
        if(idproc.eq.0)write(6,*)'flushing file unit ',nunit
        call myflush(nunit)
      endif
      return
c
c autotrack:
756   continue
c---------------------------------------
cryne 4/14/04 crms(1) is no longer used
c crms(1) is 'set='
c     if(crms(1).eq.'true')then
c       lautotrk=1
c       lautocon=0
c       if(idproc.eq.0)write(6,*)'lautotrk=1, lautocon=0'
c     elseif(crms(1).eq.'false')then
c       lautotrk=0
c       lautocon=1
c       if(idproc.eq.0)write(6,*)'lautotrk=0, lautocon=1'
c     else
c       if(idproc.eq.0)then
c       write(6,*)'error(autotrack): do not understand parameter'
c       endif
c       call myexit
c     endif
      if(crms(2).ne.'undefined' .or. crms(4).eq.'true')then
        lautotrk=1
        lautocon=0
      else
        if(idproc.eq.0)then
        write(6,*)'error: autotrack has been invoked, but neither'
        write(6,*)'a particle tracking algorithm nor the envelope'
        write(6,*)'option have been specified'
        endif
        call myexit
      endif
c
c crms(2) is 'type=' [specifies particle tracking algorithm]
      ntrktype=0
      ntrkorder=0
      if(crms(2).ne.'undefined')then
        call getordtyp(crms(2),ntrktype,ntrkorder)
      endif
c
c crms(4) is 'env='
      if(crms(4).eq.'true')then
        nenvtrk=1
      else
        nenvtrk=0
      endif
c
c crms(3) is 'sckick='
c  1 means it has been turned on; 0 means it has been turned off;
c -1 means it has not been set. done this way for future use.
      if(crms(3).eq.'true')then
        nsckick=1
      elseif(crms(3).eq.'false')then
        nsckick=0
      else
        nsckick=-1
      endif
c write results of parsing this command for debug:
      if(idproc.eq.0)then
        if(ntrktype.ne.0)then
          write(6,*)'autotracking particles'
          if(ntrktype.eq.1)write(6,*)'ntrktype=1 (taylor)'
          if(ntrktype.eq.2)write(6,*)'ntrktype=2 (symplectic)'
          write(6,*)'ntrkorder=',ntrkorder
          write(6,*)'nsckick=',nsckick
        endif
        if(nenvtrk.eq.1)write(6,*)'autotracking envelopes'
      endif
      return
c
c sckick:
757   continue
      if(idproc.eq.0)write(6,*)'sckick command not implemented yet'
      return
c
c moments (write some 2nd order moments):
758   continue
      nfile1=nint(prms(1))  !assigned unit number for output xfile
      nfile2=nint(prms(2))  !assigned unit number for output yfile
      nfile3=nint(prms(3))  !assigned unit number for output tfile
      nprecision=nint(prms(4))-1
      nunits=nint(prms(5))
      fname1=crms(1)
      fname2=crms(2)
      fname3=crms(3)
      estrng='moments'
c1:
      if(nfile1.eq.0)then
        call fnamechk(fname1,nfile1,ierr,estrng)
        if(ierr.eq.1)write(6,*)'moments error, problem w/ xfile'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile1)
      if(nfile1.gt.0)prms(1)=nfile1
c2:
      if(nfile2.eq.0)then
        call fnamechk(fname2,nfile2,ierr,estrng)
        if(ierr.eq.1)write(6,*)'moments error, problem w/ yfile'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile2)
      if(nfile2.gt.0)prms(2)=nfile2
c3:
      if(nfile3.eq.0)then
        call fnamechk(fname3,nfile3,ierr,estrng)
        if(ierr.eq.1)write(6,*)'moments error, problem w/ tfile'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile3)
      if(nfile3.gt.0)prms(3)=nfile3
c7:
      if(crms(7).eq.'ratio')then
        ncorr=1
      else
        ncorr=0
      endif
c8:
      if(crms(8).eq.'true')then
        includepi=1
      else
        includepi=0
      endif
c9:
      if(crms(9).eq.'remove')then
        ncent=0
      elseif(crms(9).eq.'keep')then
        ncent=1
      else
        if(idproc.eq.0)
     &    write(6,*) 'error(emittance): centroid=',TRIM(crms(4)),
     &               ' is invalid.  Use keep or remove.'
        call myexit()
      endif
c===
      call writemom2d(arclen,nfile1,nfile2,nfile3,nprecision,nunits,    &
     &ncorr,includepi,ncent)
c===
      if(crms(4).eq.'true')call myflush(nfile1)
      if(crms(5).eq.'true')call myflush(nfile2)
      if(crms(6).eq.'true')call myflush(nfile3)
      return
c
c
c maxsize (write max beam sizes)
759   continue
      nfile1=nint(prms(1))  !assigned unit number for output file 1
      nfile2=nint(prms(2))  !assigned unit number for output file 2
      nfile3=nint(prms(3))  !assigned unit number for output file 3
      nfile4=nint(prms(4))  !assigned unit number for output file 4
      nprecision=nint(prms(5))-1
      nunits=nint(prms(6))
      if(crms(9).eq.'automatic'.or.crms(9).eq.'auto')then
c       if(idproc.eq.0.and.nfile1.eq.0)                                   &
c    &    write(6,*)'(maxsize) automatic file names'
        fname1='xmax.out'
        fname2='ymax.out'
        fname3='tmax.out'
        fname4='zmax.out'
      else
        fname1=crms(1)
        fname2=crms(2)
        fname3=crms(3)
        fname4=crms(4)
      endif
      estrng='maxsize'
c1:
      if(nfile1.eq.0)then
        call fnamechk(fname1,nfile1,ierr,estrng)
        if(ierr.eq.1)write(6,*)'maxsize error, problem w/ fname1'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile1)
      if(nfile1.gt.0)prms(1)=nfile1
c2:
      if(nfile2.eq.0)then
        call fnamechk(fname2,nfile2,ierr,estrng)
        if(ierr.eq.1)write(6,*)'maxsize error, problem w/ fname2'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile2)
      if(nfile2.gt.0)prms(2)=nfile2
c3:
      if(nfile3.eq.0)then
        call fnamechk(fname3,nfile3,ierr,estrng)
        if(ierr.eq.1)write(6,*)'maxsize error, problem w/ fname3'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile3)
      if(nfile3.gt.0)prms(3)=nfile3
c4:
      if(nfile4.eq.0)then
        call fnamechk(fname4,nfile4,ierr,estrng)
        if(ierr.eq.1)write(6,*)'maxsize error, problem w/ fname4'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile4)
      if(nfile4.gt.0)prms(4)=nfile4
c===
!XXX -- nunits arg is not in subroutine WRITEMAXSIZE()
      call writemaxsize(arclen,nfile1,nfile2,nfile3,nfile4,nprecision)   !XXX,     &
!XXX     &nunits)
c===
      if(crms(5).eq.'true')call myflush(nfile1)
      if(crms(6).eq.'true')call myflush(nfile2)
      if(crms(7).eq.'true')call myflush(nfile3)
      if(crms(8).eq.'true')call myflush(nfile4)
c     if(idproc.eq.0)then
c       write(6,*)'s=',arclen,'; maximum beam sizes written'
c     endif
      return
c
c reftraj:
760   continue
      nfile1=nint(prms(1))  !assigned unit number for output file
      nprecision=nint(prms(2))-1
      nunits=nint(prms(3))
      fname1=crms(1)
      estrng='reftraj'
      if(nfile1.eq.0)then
        call fnamechk(fname1,nfile1,ierr,estrng)
        if(ierr.eq.1)write(6,*)'(reftraj) problem w/ fname1:',fname1
        if(ierr.eq.1)call myexit
      endif
      call ibcast(nfile1)
      prms(1)=nfile1
c===
      call writereftraj(arclen,nfile1,nprecision,nunits)
c===
      if(crms(2).eq.'true')call myflush(nfile1)
      return
c
c initenv:
761   continue
      if(crms(1).eq.'ratio')then
        ncorr=1
      else
        ncorr=0
      endif
      call initenv(prms,ncorr)
      return
c
c envelopes: (write the envelopes stored in the env array)
762   continue
c     write(6,*)'here I am at envelopes'
      nfile1=nint(prms(1))  !assigned unit number for output file 1
      nfile2=nint(prms(2))  !assigned unit number for output file 2
      nfile3=nint(prms(3))  !assigned unit number for output file 3
c     write(6,*)'nfile1,nfile2,nfile3=',nfile1,nfile2,nfile3
      nprecision=nint(prms(4))-1
      nunits=nint(prms(5))
      if(crms(7).eq.'ratio')then
        ncorr=1
      else
        ncorr=0
      endif
      fname1=crms(1)
      fname2=crms(2)
      fname3=crms(3)
c     write(6,*)'fname1,fname2,fname3=',fname1,fname2,fname3
      estrng='envelopes'
c1:
      if(nfile1.eq.0)then
        call fnamechk(fname1,nfile1,ierr,estrng)
        if(ierr.eq.1)write(6,*)'envelopes error, problem w/ fname1'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile1)
      if(nfile1.gt.0)prms(1)=nfile1
c2:
      if(nfile2.eq.0)then
        call fnamechk(fname2,nfile2,ierr,estrng)
        if(ierr.eq.1)write(6,*)'envelopes error, problem w/ fname2'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile2)
      if(nfile2.gt.0)prms(2)=nfile2
c3:
      if(nfile3.eq.0)then
        call fnamechk(fname3,nfile3,ierr,estrng)
        if(ierr.eq.1)write(6,*)'envelopes error, problem w/ fname3'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile3)
      if(nfile3.gt.0)prms(3)=nfile3
c     write(6,*)'done checking file names'
c===
      call writeenv2d(arclen,nfile1,nfile2,nfile3,nprecision,nunits,    &
     &ncorr)
c     write(6,*)'returned from writeenv2d'
c===
c===
      if(crms(4).eq.'true')call myflush(nfile1)
      if(crms(5).eq.'true')call myflush(nfile2)
      if(crms(6).eq.'true')call myflush(nfile3)
      return
c contractenv (apply contraction map to envelopes):
763   continue
      if(idproc.eq.0)then
      write(6,*)'(lmnt) warning: user is calling contractenv directly,'
      write(6,*)'but this is normally called by ML/I when the user'
      write(6,*)'specifies a matchenv command. Make sure you know'
      write(6,*)'what you are doing!'
      endif
      call contractenv(delta)
      return
c
c setreftraj:
764   continue
      nsto=nint(p(8))
      nget=nint(p(9))
      nwrt=nint(p(10))
c reset from initial values? :
c crms(1) corresponds to "restart=' which refers to data in location 9
      if(crms(1).eq.'true')nget=9
c include arc length? :
      includearc=1     ! default is to include arc length
      if(crms(2).eq.'false')includearc=0
c
      if(nsto.gt.0. or. nget.gt.0)then
        if(nsto.gt.0)then
          if(idproc.eq.0)then
           if(nwrt.gt.0)then
           write(6,*)'storing ref particle data in location',nsto
           endif
          endif
          refsave(nsto,1:6)=reftraj(1:6)
          if(includearc.eq.1)arcsave(nsto)=arclen
          brhosav(nsto)=brho
          gamsav(nsto)=gamma
          gam1sav(nsto)=gamm1
          betasav(nsto)=beta
        else
          if(idproc.eq.0)then
           if(nwrt.gt.0)then
             write(6,*)'getting ref particle data from location ',nsto
           endif
          endif
          reftraj(1:6)=refsave(nget,1:6)
          if(includearc.eq.1)arclen=arcsave(nget)
          brho=brhosav(nget)
          gamma=gamsav(nget)
          gamma1=gam1sav(nget)
          beta=betasav(nget)
        endif
      else
        if(idproc.eq.0 .and. nwrt.gt.0)then
         write(6,*)'setting ref particle data'
         write(6,*)'x,px=',p(1),p(2)
         write(6,*)'y,py=',p(3),p(4)
         write(6,*)'t,pt=',p(5),p(6)
        endif
        reftraj(1)=p(1)
        reftraj(2)=p(2)
        reftraj(3)=p(3)
        reftraj(4)=p(4)
        reftraj(5)=p(5)
        reftraj(6)=p(6)
c only change arc length if the user has provided it (set to -9999 in sif.f):
        if(p(7).ne.-9999.d0)arclen=p(7)
        gamma=-reftraj(6)/pmass*omegascl*sl*p0sc
        gamm1=gamma-1.d0
        gbet=sqrt(gamm1*(gamma+1.d0))
        clite=299792458.d0
        brho=gbet/clite*pmass
        beta=gbet/gamma
      endif
      return
c
c setarclen:
765   continue
c crms(1) corresponds to "restart='
      if(crms(1).eq.'true')then
        arclen=0.d0
      else
        arclen=p(1)
      endif
      nwrt=nint(p(2))
      if(nwrt.gt.0 .and. idproc.eq.0)then
        write(6,*)'setting arc length to ',arclen
      endif
      return
c
c wakedefault: (actually, not needed here; dealt with in sif.f)
766   continue
      return
c
c emittance (write 2D, 4D, and/or 6D 2nd order moments):
767   continue
      nfile=nint(prms(1))  !assigned unit number for output file
      nprecision=nint(prms(2))-1
      nunits=nint(prms(3))
      fname=crms(5)
      estrng='emittance'
c
      if(nfile.eq.0)then
        call fnamechk(fname,nfile,ierr,estrng)
        if(ierr.eq.1)write(6,*)'error(emittance): problem w/ file'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile)
      if(nfile.gt.0)prms(1)=nfile
      ne2=1
      if(crms(1).eq.'false')ne2=0
      ne4=1
      if(crms(2).eq.'false')ne4=0
      ne6=1
      if(crms(3).eq.'false')then
        ne6=0
      else
        write(6,*) 'warning(emittance): 6d=true not yet implemented!'
        ne6=0
      endif
      if(crms(4).eq.'remove')then
        ncent=0
      elseif(crms(4).eq.'keep')then
        ncent=1
      else
        if(idproc.eq.0)
     &    write(6,*) 'error(emittance): centroid=',TRIM(crms(4)),
     &               ' is invalid.  Use keep or remove.'
        call myexit()
      endif
c===
      if(ne2.eq.1.or.ne4.eq.1.or.ne6.eq.1) then
        call writeemit(arclen,nfile,nprecision,nunits,ne2,ne4,ne6,ncent)
      endif
c===
      if(crms(6).eq.'true')call myflush(nfile)
      return
c
c matchenv
768   continue
      if(idproc.eq.0)then
        write(6,*)'error (lmnt): at matchenv in subroutine lmnt'
        write(6,*)'but the code should not get here.'
        write(6,*)                                                          &
     &  'a matchenv command can ONLY be placed under #labor in mli.in'
      endif
      call myexit
c     stop
c
c fileinfo:
769   continue
      ifileinfo=nint(p1)
      return
c
c egengrad: compute generalized gradients from E-field surface data
770   continue
      nfile=nint(prms(1))  !assigned unit number for output file
      zst=prms(2)
      zen=prms(3)
      nz=nint(prms(4))
      f=prms(5)
      r=prms(6)
      fkmx=prms(7)
      nk=nint(prms(8))
      infiles=nint(prms(9))
      nprec=nint(prms(10))-1
      !fnin=crms(1)       ! filename for input data (E-field)
      !fnout=crms(2)      ! filename for output data (genlzd grads)
      estrng='egengrad'
      kfile=0
      jfile=0
c
      if(nfile.eq.0)then
        call fnamechk(crms(2),nfile,ierr,estrng)
        if(ierr.eq.1)write(6,*)'error(egengrad): problem w/ input file'
        if(ierr.eq.1)return
      endif
      call ibcast(nfile)
      if(nfile.gt.0)prms(1)=nfile
      if(crms(5).eq.'true')then
        call fnamechk(crms(4),kfile,ierr,estrng)
        if(ierr.eq.1)write(6,*)'error(egengrad): problem w/ char file'
        if(ierr.eq.1)return
      endif
      call ibcast(kfile)
      if(crms(7).ne.' ')then
        call fnamechk(crms(7),jfile,ierr,estrng)
        if(ierr.eq.1)write(6,*)'error(egengrad): problem w/ diagn file'
        if(ierr.eq.1)return
      endif
      call ibcast(jfile)
      call cegengrad(crms(1),crms(6),infiles,nfile,kfile,jfile,         &
     &               zst,zen,nz,f,r,fkmx,nk,nprec)
      if(crms(3).eq.'true')then
        call myflush(nfile)
        if(kfile.ne.0) call myflush(kfile)
        if(jfile.ne.0) call myflush(jfile)
      end if
      return
c
c wrtmap:  Write a map to file
771   continue
c
c     crms(1) : filename for write
c     crms(2) : kind of map to write ('accumulated', 'lastslice')
c     crms(3) : i/o status for write ('overwrite','append')
c
      ifile=0
      estrng='wrtmap'
      call fnamechk(crms(1),ifile,ierr,estrng)
      if(ierr.eq.1)then
         if(idproc.eq.0)write(6,*)'error(wrtmap): problem w/ file'
         return
      endif
c
c If overwriting file, rewind...
c KMP: This is a hack, and this should be passed to the fnamechk routine
c      in order to be implemented properly.
      if(crms(3).eq.'overwrite')then
         close(ifile)
         open(unit=ifile, file=crms(1), status='unknown', err=7710)
         goto 7711
7710     if(idproc.eq.0)write(6,*)'error(wrtmap): problem overwriting'
         ierr=1
         return
7711     continue
      endif
c
c Only write file from processor 0
      if(idproc.eq.0)then
         if(crms(2).eq.'accumulated')then
!             call wrtmap(ifile,refsave(9,1),reftraj,arclen,th,tmh)
            call wrtmap(ifile,refprev,reftraj,arclen,th,tmh)
         elseif(crms(2).eq.'lastslice')then
            call wrtmap(ifile,refprev,reftraj,prevlen,hprev,mhprev)
         else
            ierr=1
            write(6,*)'error(wrtmap): invalid kind of map [',crms(2),']'
         endif
      endif
      return
c rdmap:
772   continue
c
c     crms(1) : filename to read from
c     crms(2) : rewind file or not
c
      ifile=0
      estrng='rdmap'
      call fnamechk(crms(1),ifile,ierr,estrng)
      if(ierr.eq.1)then
         if(idproc.eq.0)write(6,*)'error(rdmap): problem w/ file'
         return
      endif
      if(idproc.eq.0)write(6,*)'rdmap: reading map from ',crms(1)
      if(crms(2).eq.'true')then
         if(idproc.eq.0)write(6,*)'       rewinding file before use'
         rewind ifile
      endif

      mpitmp=mpi
      mpi=ifile
      call rdmap(ifile,darclen,dreftraj,h,mh)
      mpi=mpitmp

      refprev=reftraj
      arclen=arclen+darclen
      prevlen=darclen
      nt2prev=nt2
      reftraj=reftraj+dreftraj
      goto 2000

c sparec7:
773   continue
      return
c sparec8:
774   continue
      return
c sparec9:
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
801   call cod(p,th,tmh)
      return
c
c     apply map to a function or moments
c
802   call amap(p,th,tmh)
      return
c
c     dynamic invariant analysis
c
803   call dia(p,th,tmh)
      return
c
c     dynamic normal form analysis
c
804   call dnor(p,th,tmh)
      return
c
c     compute exponential
c
805   call cex(p,th,tmh)
      return
c
c     compute power of dynamic normal form
c
806   call pdnf(p,th,tmh)
      return
c
c     compute power of static normal form
c
807   call psnf(p,th,tmh)
      return
c
c     resonance analyze dynamic map
c
808   call radm(p,th,tmh)
      return
c
c     resonance analyze static map
c
809   call rasm(p,th,tmh)
      return
c
c     static invariant ayalysis
c
810   call sia(p,th,tmh)
      return
c
c     static normal form analysis
c
811   call snor(p,th,tmh)
      return
c
c     twiss analyze dynamic map
c
812   call tadm(p,th,tmh)
      return
c
c     twiss analyze static map
c
813   call tasm(p,th,tmh)
      return
c
c     translate basis
c
814   call tbas(p,th,tmh)
      return
c
c     get buffer contents
c
815   continue
c
c test control parameters
c
      nmap=nint(p2)
      if (nmap.gt.5 .or. nmap.lt.1) then
        write(jof,9815) nmap
 9815   format(1x,'nmap=',i3,1x,'trouble with gbuf:nmap < 1 or > 5')
        call myexit
      endif
c option when tracking (procedure at end)
      if(ntrk .eq. 1) then
c       if(idproc.eq.0)write(6,*)'at gbuf w/ ntrk=',ntrk
        p1temp=p(1)
        p(1)=2
        call gbuf(p,h,mh)
        p(1)=p1temp
        goto 2000
      endif
c options when not tracking
c     if(idproc.eq.0)write(6,*)'at gbuf w/ ntrk.ne.1, ntrk=',ntrk
      call gbuf(p,th,tmh)
        return
c
c     transport static (script) A
c
816   call trsa(p,th,tmh)
      return
c
c    transport dynamic script A
c
817   call trda(p,th,tmh)
      return
c
c     multiply polynomial by a scalar
c
818   call smul(p,th,tmh)
      return
c
c     add two polynomials
c
819   call padd(p,th,tmh)
      return
c
c     multiply two polynomials
c
820   call pmul(p,th,tmh)
      return
c
c     Poisson bracket two polynomials
c
821   call pbpol(p,th,tmh)
      return
c
c     polar decompose matrix portion of transfer map
c
822   call pold(p,th,tmh)
      return
c
c     evaluate a polynomial
c
823   call pval(p,th,tmh)
      return
c
c     fourier analyze static map
c
824   call fasm(p,th,tmh)
      return
c
c     fourier analyze dynamic map
c
825   call fadm(p,th,tmh)
      return
c
c     select quantities
c
826   call sq(p)
      return
c
c     write selected quantities
c
827   call wsq(p)
      return
c
c     change tune ranges
c
828   continue
      call subctr(p)
      return
c
c     apply script N inverse
c
829   continue
      call asni(p)
      return
c
c     compute power of nonlinear part
c
830   continue
      call pnlp(p,th,tmh)
      return
c
c     check for symplecticity
c
831   continue
      isend=nint(p1)
      call csym(isend,tmh,ans)
      return
c
c     (psp) compute scalar product of two polynomials
c
832   call psp(p,th,tmh)
      return
c
c     (mn) compute matrix norm
c
833   call submn(p,th,tmh)
      return
c
c     (bgen) generate a beam
c
834   call bgen(p,th,tmh)
      return
c
c     (tic) translate (move) initial conditions
c
835   call tic(prms)
      return
c
c     (ppa) principal planes analysis
c
836   call ppa(p,th,tmh)
      return
c
c     (moma) moment and map analysis
c
837   call moma(p)
      return
c
c     (geom) compute geometry of a loop
c
838   call geom(p)
      return
c
c     (fwa) copy file to working array
c
839   call fwa(p)
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
901   continue
cryne make sure that  the code knows it has not converged yet
cryne (supresses printing until convergence is achieved)
      kfit=0
      call bip(p)
      return
902   call bop(p)
      return
c
c     end procedures
c
903   continue
      call subtip(p)
      return
904   call subtop(p)
      return
c
c     specify aims
c
905   call aim(p)
      return
c
c     specify quantities to be varied
c
906   call vary(p)
      return
c
c     fit to achieve aims
c
907   continue
cryne
ctm   write(6,*)'(lmnt)calling fit; kfit=',kfit
      call fit(p)
ctm   write(6,*)'(lmnt)back from fit; kfit=',kfit
      return
c
c     optimize
c
908   call opt(p)
      return
c
c     constraints
c
909   call con1(p)
      return
910   call con2(p)
      return
911   call con3(p)
      return
912   call con4(p)
      return
913   call con5(p)
      return
c
c     merit functions
c
c     least squares merit function
c
914   call mrt0
      return
c
c     user supplied merit functions
c
915   call mrt1(p)
      return
916   call mrt2(p)
      return
917   call mrt3(p)
      return
918   call mrt4(p)
      return
919   call mrt5(p)
      return
c
c     free parameter sets
c
920   ipset=nint(p1)
      call fps(ipset)
      return
c
c     capture parameter sets
c
921   ipset=1
      call cps(prms,p,ipset)
      return
922   ipset=2
      call cps(prms,p,ipset)
      return
923   ipset=3
      call cps(prms,p,ipset)
      return
924   ipset=4
      call cps(prms,p,ipset)
      return
925   ipset=5
      call cps(prms,p,ipset)
      return
926   ipset=6
      call cps(prms,p,ipset)
      return
927   ipset=7
      call cps(prms,p,ipset)
      return
928   ipset=8
      call cps(prms,p,ipset)
      return
929   ipset=9
      call cps(prms,p,ipset)
      return
c
c     compute dynamic aperture (dapt)
c
930   continue
      call dapt(p)
      return
c
c     gradient (grad)
c
931   continue
      call grad(p)
      return
c
c     rset
c
932   continue
      call rset(p)
      return
c
c     flag
c
933   continue
      call flag(p)
      return
c
c     scan
c
934   continue
      call scan(p)
      return
c
c     mss
c
935   continue
      call mss(p)
      return
c
c     spare1
c
936   continue
      write(6,*) 'spare1 not yet available'
      return
c
c     spare2
c
937   continue
      write(6,*) 'spare2 not yet available'
      return
c
c  ......... concatenate or track before exiting ........
c
 2000 continue
cryne--- 08/16/2001
cryne This statement has been added because h and mh disappear after
cryne leaving the routine, and there may be times when we would like
cryne to make use of the last map that was constructed.
      call mapmap(h,mh,hprev,mhprev)
cryne---
cryne 08/16/2001 note well: now the user can turn off auto-concatenation
cryne Use with care!!!
      if(ntrk.eq.0)then
        if(lautocon.eq.0)write(6,*)'WARNING: not concatenating!'
        if(lautocon.eq.1)call concat(th,tmh,h,mh,th,tmh)
        if(lautocon.eq.2)call sndwchi(h,mh,th,tmh,th,tmh)
        if(lautocon.eq.-2)call sndwch(h,mh,th,tmh,th,tmh)
        return
      else
cryne 3/27/04 ntaysym=1 for taylor, =2 for symplectic
        if(ntrktype.ne.0)then
          ntaysym=ntrktype
          call trace(0,jfcf,ntaysym,ntrkorder,1,0,0,0,5,0,h,mh)
        endif
        if(nenvtrk.eq.1)then
c         write(6,*)'calling envtrace(mh) with mh(1:6,1:6)='
c         write(6,51)mh(1:6,1:6)
   51     format(6(1x,1pe12.5))
          call envtrace(mh)
        endif
c       write(6,*)'returning from bottom of lmnt'
        return
      endif
      end
c
************************************************************************
c
      subroutine lookup(string,itype,index)
c-----------------------------------------------------------------------
c  this subroutine determines whether the input string 'string'
c  is an element,line,lump,loop, or unused label.
c
c  input:  string  character*(*) item name
c  output: itype   integer       =1, if string is a menu entry
c                                =2, .... line
c                                =3, .... lump
c                                =4, .... loop
c                                =5, .... unknown label
c          index   integer       index of 'string' in its array
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c             October 30, 1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c parameter types
c-----------------------------------------------------------------------
      character string*(*)
      integer itype,index
c-----------------------------------------------------------------------
c start routine
c-----------------------------------------------------------------------
c element?
      do 10 n=1,na
      if(string.eq.lmnlbl(n)) then
       itype = 1
       index = n
       return
      endif
   10 continue
c item?
      do 20 n=1,nb
      if(string.eq.ilbl(n)) then
       itype = ityp(n)
       index = n
       return
      endif
   20 continue
c not found:
      itype=5
      return
      end
c
***********************************************************************
c
      subroutine low(line)
c  Converts all uppercase characters to lowercase.
c  Written by Liam Healy, Feb. 28, 1985.
c
      character line*(*)
c
      character*26 lower,upper
      character*1 blank
      data lower/'abcdefghijklmnopqrstuvwxyz'/
      data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data blank/' '/
      save lower,upper,blank
c
      lnbc=1
      do 100 i=1,LEN(line)
        loc=index(upper,line(i:i))
        if(loc.gt.0) then
          line(i:i)=lower(loc:loc)
        endif
        if (line(i:i).ne.blank) lnbc=i
  100 continue
      return
      end
c
***********************************************************************
      subroutine lumpit(mth)
c-----------------------------------------------------------------------
c  translate map of mth lump to special format and save it in core
c
c  input mth (integer) : index of lump in /items/
c
c  Written by Rob Ryne ca 1984
c  changed by Petra Schuett
c             October 30,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'map.inc'
      include 'deriv.inc'
      include 'core.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
c inew keeps track of the lump to be deleted next, if core is full
      integer inew
      save inew
      data inew /0/
c-----------------------------------------------------------------------
c  start routine
c-----------------------------------------------------------------------
c find vacant space in core
      icore = 0
      do 1 i=maxlum,1,-1
       if(inuse(i).eq.0) icore=i
   1  continue
c if core is full, destroy oldest lump
      if(icore.eq.0) then
       inew=inew+1
       if(inew.eq.maxlum)inew=1
       icore=inew
       lmade(inuse(icore))=0
       write(jof,510) ilbl(inuse(icore))
 510   format(1h ,'lump ',a16,' deleted;')
      endif
c--------------------
c calculate rjac,df,..
      call canx(tmh,th,5)
c ..rrjac and rdf
      call rearr
c--------------------
c now store all the info
      do 10 n1=1,6
!!!!! do 10 n2=1,83
      do 10 n2=1,monoms
      dfl(n1,n2,icore)=df(n1,n2)
   10 continue
      do 20 n1=1,3
!!!!! do 20 n2=1,84
!!!!! do 20 n2=1,monom1+1
      do 20 n2=1,monoms
   20 rdfl(n1,n2,icore)=rdf(n1,n2)
      do 30 n1=1,3
      do 30 n2=1,3
!!!!! do 30 n3=1,28
!!!!! do 30 n3=1,monom2+1
      do 30 n3=1,monoms
   30 rrjacl(n1,n2,n3,icore)=rrjac(n1,n2,n3)
      do 40 n1=1,6
      do 40 n2=1,6
   40 tmhl(n1,n2,icore)=tmh(n1,n2)
      do 50 n1=1,monoms
   50 thl(n1,icore)=th(n1)
c--------------------
c set pointers
      inuse(icore) = mth
      lmade(mth)   = icore
c--------------------
      write(jof,520) ilbl(mth),icore
      write(jodf,520)  ilbl(mth),icore
  520 format(1h ,'lump ',a16,' constructed and stored.','(',i2,')')
c--------------------
!     ifile=50+icore
!     write(ifile,*)'DFL:'
!     do n1=1,6
!     do n2=1,monoms
!     if(dfl(n1,n2,icore).eq.0.d0)cycle
!     write(ifile,*)n1,n2,dfl(n1,n2,icore)
!     enddo
!     enddo
!     write(ifile,*)'RDFL:'
!     do n1=1,6
!     do n2=1,monoms
!     if(rdfl(n1,n2,icore).eq.0.d0)cycle
!     write(ifile,*)n1,n2,rdfl(n1,n2,icore)
!     enddo
!     enddo
!     write(ifile,*)'RRJACL:'
!     do n1=1,3
!     do n2=1,3
!     do n3=1,monoms
!     if(rrjacl(n1,n2,n3,icore).eq.0.d0)cycle
!     write(ifile,*)n1,n2,n3,rrjacl(n1,n2,n3,icore)
!     enddo
!     enddo
!     enddo
!     write(ifile,*)'TMHL:'
!     do n1=1,6
!     do n2=1,6
!     if(tmhl(n1,n2,icore).eq.0.d0)cycle
!     write(ifile,*)n1,n2,tmhl(n1,n2,icore)
!     enddo
!     enddo
!     write(ifile,*)'THL:'
!     do n1=1,monoms
!     if(thl(n1,icore).eq.0.d0)cycle
!     write(ifile,*)n1,thl(n1,icore)
!     enddo
c--------------------
      return
      end
c
************************************************************************
c
      function newsl(mp,kth)
c-----------------------------------------------------------------------
c newsl is the first icon to be treated. Depending on the sign of the
c repetition factor loop(mp), it is either the first or the last one.
c
c Petra Schuett, November 6,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
c--------
c commons
c--------
      include 'stack.inc'
c-------
c start
c-------
      if(loop(mp).ge.0) then
        newsl = 1
      else
        newsl = ilen(kth)
      endif
      return
      end
************************************************************************
      function npm1(number)
c Petra Schuett, November 6,1987
c
      if (number .ge. 0) then
        npm1 = 1
      else
        npm1 = -1
      endif
      return
      end
************************************************************************
      subroutine pop(lempty)
c-----------------------------------------------------------------------
c  pop stack
c
c  output: lempty = .true. if stack is empty
c
c  Petra Schuett  October 30,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      logical lempty
c-------
c start
c-------
      np = np - 1
      if(np .le. 0) then
       lempty = .true.
      else
       lempty = .false.
       call lookup(lstac(np),ntype,ith)
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
      return
      end
************************************************************************
      subroutine push
c-----------------------------------------------------------------------
c  push stack: add actual icon on top of it
c
c  Petra Schuett  October 30,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'files.inc'
c-------
c start
c-------
      np = np + 1
      if(np .gt. mstack) then
       write(jof ,510)
       write(jodf,510)
 510   format(' error in push: stack overflow')
       call myexit
      else
       lstac(np)  =icon(nslot(np-1),ith)
       loop(np)   =irep(nslot(np-1),ith) * npm1(loop(np-1))
       nslot(np)  =newsl(np,jth)
       nslot(np-1)=nslot(np-1) + npm1(loop(np-1))
       ith         = jth
       ntype       = mtype
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
      return
      end
************************************************************************
      subroutine readin(line,leof)
c-----------------------------------------------------------------------
c  Reads a line from file lf, puts it in the character variable 'line'.
c  Designed so that filters may put on, such as the routine to convert
c  all uppercase characters to lower case.
c  Written by Liam Healy, May 1, 1985.
c
c  Changed by Petra Schuett, October 19, 1987 :
c    use logical: leof=.true. , if end of file is encountered
c-----------------------------------------------------------------------
      include 'files.inc'
      character line*(*)
      logical leof
c----Routine----
      read(lf,800,end=100,err=200) line
  800 format(a)
cryne 9/11/2002 uncommented the following:
       call low(line)
      return
  100 continue
      leof=.true.
      return
  200 continue
cryne 08/26/2001 error exit added by ryne
      write(6,*)'error reading data in routine readin'
      stop
      end
c
************************************************************************
c
      subroutine rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
c-----------------------------------------------------------------------
c  This routine reads a new line and interprets it partly:
c   - if a line starts with '!' (comment line), it is skipped and
c                  another line is read
c   - if a #... code is read, it sets msegm, indicating the start of
c                  a new input component (segment).
c   - lines in the comment component (segment) are not interpreted.
c   - all other lines are cut into stings and numbers assuming a form:
c        n1*str1 n2*str2 ... n3*str3 n4*str4,n5*str5  (etc)
c        the numbers n1,n2,... are optional
c        the strings and numbers are separated by blanks and/or commas
c   - a line ending with & is continued on the next line (lcont=true)
c   - anything after a '!' is ignored
c
c  Input: msegm      integer      input component (segment) currently
c                                 being read
c
c  Output:line       character*(*)line read
c         leof       logical      =.true. if end of file is encountered
c         msegm      integer      future component (segment) to be read
c         strarr(40) character*(*)array of strings str1,str2,...
c         narr(40)   integer      array of numbers n1,n2,...
c         itot       integer      number of strings found
c         lcont      logical      =.true. if line to be continued
c
c  Author: Petra Schuett
c          October 19, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
c
      include 'files.inc'
      include 'codes.inc'
c arguments:
      integer msegm,narr(40),itot
      character strarr(40)*(*)
      character line*(*)
      logical leof,lcont,npound
c local variables:
      character*16 string
      logical lnum,lfound
      logical defcon,defelem,defline
      external defcon,defelem,defline
cryne 5/4/2006
      logical MAD8,MADX
      common/mad8orx/MAD8,MADX
c-----------------
c init
      lcont=.false.
      npound=.false.
c
ctm 9/01  skip read if coming from #include
c
      if(itot.eq.-1) go to 3
c read new line
  1   call readin(line,leof)
c check for end of file
      if (leof) return
c-----------------------------------------
c check for blank lines or comment lines
c find the first nonblank character:
      do n=1,LEN(line)
        nfirst=n
        if(line(n:n).ne.' ')goto 100
      enddo
c blank line:
      goto 1
  100 continue
c is the first character a '!' ?
      if(line(nfirst:nfirst).eq.'!')goto 1

c-----------------------------------------
c     write(6,*)'(rearec got)',line
c skip this check if entering the routine in #comment or #beam
!     if((msegm.eq.1).or.(msegm.eq.2))goto 123
cryne+++ August 5, 2004 : some time ago I commented out the previous line.
cryne    However, it should be present for the case msegm.eq.1, so I am
cryne    putting it back now. The reason is that the parser could find
cryne    symbolic expressions in the comments, and it should not interpret
cryne    these as symbolic constants whose values need to be calculated.
cryne    Due to the following if test, the only way to get out of
cryne    the #comments section is by finding #something (normally #beam)
      if(msegm.eq.1)goto 123
cryne+++
c
      if(msegm.ne.9)then
        if(defcon(line))msegm=9
      endif
      if(msegm.ne.3)then
        if(defelem(line))msegm=3
      endif
      if(msegm.ne.4)then
        if(defline(line))msegm=4
      endif
  123 continue
c
cryne 5/4/2006
cryne I am struggling w/ MADX input (originally wrote for MAD8).
cryne Deal with this by deciding here if the line is to be continued,
cryne First blank out everything to the right of and including '!'
      lenline=len(line)
      n123=index(line,'!')
      if(n123.ne.0)line(n123:lenline)=' '
cryne in MADX style, there could be space between ; and !, so deal with it:
      if(MADX)then
        if(n123.ne.0)iqtop=n123-1
        if(n123.eq.0)iqtop=lenline
        do iq=iqtop,1,-1
          ilast=iq
          if(line(ilast:ilast).ne.' ')exit
        enddo
c       write(6,*)'ilast,line(ilast:ilast)=',ilast,line(ilast:ilast)
        if(line(ilast:ilast) .NE. ';' )lcont=.true.
      endif
cryne
c
c
  3   continue
      itot = 0
c convert to lower case, except for comment component (segment)
      if (msegm.ne.1) call low(line)
c find first string
      kbeg=1
      call cread(kbeg,msegm,line,string,lfound)
c empty line is ignored, except in comment component (segment)
      if(.not.lfound) then
        if(msegm.eq.1) then
         return
        else
cryne 7/7/2002         write(6,*) msegm,'= segment with blank line'
         goto 1
        endif
      endif
c new component (segment) starts...
      if(string(1:1).eq.'#')then
c in case previous component (segment) was comment, conv to lower case
        call low(line)
c ...which one?
cryne   do 2 k=1,8 !July 7, 2002
        do 2 k=1,nintypes
          if(string(1:8).eq.ling(k)) then
            msegm = k
            npound=.true.
c           write(6,*)'found a ',ling(k)
ctm            if(msegm.eq.8)then
ctm              call cread(kbeg,msegm,line,string,lfound)
ctm              if(.not.lfound)then
ctm not any more-- write(6,*)'error(rearec): file name must follow #include'
ctm              endif
ctm              strarr(1)=string(1:8)
ctm            return
ctm            endif
c            write(jodf,*) msegm
cryne 5/4/2006 a ";" on the line with #comment means "use MADX style"
            if(msegm.eq.1)then
              if(index(line,';').ne.0)then
                MAD8=.false.
                MADX=.true.
              endif
            endif
c------------------
            return
          endif
  2     continue
c ... no match:
c default is #beam after #comment
        if(msegm.eq.1) then
          msegm = 2
          return
        else
c but in all other cases, this should not happen!
          write(jof,99) string
  99      format(' ---> warning from rearec:'/                          &
     &           '      user name ',a,' begins with #')
        endif
      endif
c if #comment line is read, no interpretation
      if(msegm.eq.1) return
c comment line
      if (string(1:1).eq.'!') goto 1
c.........................................................
c now interpret line
c first init itot
      itot = 0
      do 10 i=1,40
c end of line
      if((.not.lfound).or.(string(1:1).eq.'!')) return
c line to be continued
      if((MAD8) .and. string(1:1).eq.'&') then
        lcont = .true.
        return
      endif
ccc   if((MADX).and.len_trim(string).eq.1.and.string(1:1).ne.';') then
ccc     lcont = .true.
ccc     return
ccc   endif
c so we found another string
      itot=itot+1
c is it a number?
      call cnumb(string,num,lnum)
c      write(jodf,*)string,'=',num,'lnum=',lnum
      if((.not.lnum).and.(string(1:1).ne.'-')) then
c.. this must be "string"
        narr(i)=1
        strarr(i)=string(1:16)
      else if(.not.lnum) then
c.. it is "-string"
        narr(i)=-1
cryne rhs should say string(2:29)???
        strarr(i)=string(2:16)
      else
c.. it is the number of "n*string"
        narr(i)=num
        call cread(kbeg,msegm,line,string,lfound)
c       write(jodf,*)kbeg,string,lfound
cryne 08/14/2001 normally this would indicate a problem, but due to
cryne changes in the #beam component, it could be ok. So skip warning.
        if (.not.lfound .and. msegm.eq.2)return
c.. error
        if (.not.lfound) then
          write(jof ,98) line
          write(jodf,98)
  98      format(' ---> warning from rearec:'/,                         &
     &           '      the following line contains a number which is', &
     &           ' not followed by a string:')
          write(jodf,*) '      ',line
          call myexit
        endif
c.. normal way
        strarr(i)=string(1:16)
      endif
cryne July 4, 2002
cryne if using the Standard Input Format to read menu items,
cryne it is only necessary to see if there are at least 3 strings
cryne on this record
      if(msegm.eq.3 .and. itot.eq.3)return
cryne July 14/2002
cryne if reading the definition of a constant (msegm=9), only need one.
      if(msegm.eq.9 .and. itot.eq.1)return
c find next string and start over again
      call cread(kbeg,msegm,line,string,lfound)
c      write(jodf,*)kbeg,string,lfound
  10  continue
c
c more than 40 strings, which are separated by single characters,
c cannot occur in a line of 80 characters.
      end
c
c***********************************************************************
c
      subroutine tran
c-----------------------------------------------------------------------
c organizes translation of input into work to be done
c
c Written by Rob Ryne ca 1984
c adapted to new version 9 Nov 87 by
c Petra Schuett
c Modified 31 Aug 88 by Alex Dragt
c Modified 20 Aug 98 to use a home-made do loop.  AJD
c 4/1/00: Home-made loop re-inserted into Mottershead version. RDR
c-----------------------------------------------------------------------
      use parallel
      use acceldata
      use beamdata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'map.inc'
      include 'core.inc'
      include 'loop.inc'
      include 'files.inc'
      include 'labpnt.inc'
      include 'fitbuf.inc'
      include 'setref.inc'   ! added by RDR, April 18, 2004
      real*8 mhprev
      character*16 strng
      dimension reftmp(6)
      common/prevmap/hprev(monoms),mhprev(6,6)
      common/envdata/env(6),envold(6),emap(6,6)
      common/envstuff/nenvtrk
      common/autotrk/lautotrk,ntrktype,ntrkorder
c
cryne 5/4/2006 added this to allowing printing of messages from #labor:
      integer, parameter :: lmaxmsg=1000
      character*256 lattmsg(lmaxmsg)
      common/lattmsgb/lattmsg,lmsgpoi
      character*80 linestar
c
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
      dimension rh(monoms),rmh(6,6)
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c initialize
      call ident(th,tmh)
cryne 08/17/201:
      call ident(hprev,mhprev)
ctm      reftraj(1:5)=0.
      do k = 1,5
        reftraj(k)=0.
      enddo
c
cryne 5/4/2006
      linestar='********************************************************&
     &*************************'
c
c     if(idproc.eq.0)write(6,*)'setting reference trajectory'
c     if(idproc.eq.0)write(6,*)'sl=',sl
c     if(idproc.eq.0)write(6,*)'omegascl=',omegascl
c     if(idproc.eq.0)write(6,*)'p0sc=',p0sc
c     if(idproc.eq.0)write(6,*)'gamma=',gamma
c     if(idproc.eq.0)write(6,*)'pmass=',pmass
cryne fixed 4/18/04 : set reftraj(6) correctly for all choices of units:
      reftraj(6)=-gamma*pmass/(omegascl*sl*p0sc)
cryne reftraj(6)=-gamma
c     if(idproc.eq.0)write(6,*)'reftraj(6)=',reftraj(6)
      arclen=0.
c
c save the initial reftraj data in location #9:
      refsave(9,1:6)=reftraj(1:6)
      arcsave(9)=arclen
cryne
cryne initialize this here because of a test on kfit in eintrp (gapmap)
      kfit=0
cryne
      jicnt=0
      jocnt=0
c main loop through latt (labor)
c
c Program changes made 20 August 98 by AJD.
c Change Fortran do loop into a home-made do loop because
c various other routines change the control index lp.
c
c      do 5 lp=1,noble
      lp=1
   4  continue
c     write(6,*)'AT START OF MAIN LABOR LOOP; lp=',lp
c For remaining changes, see end of this subroutine
c
cryne 5/4/2006:
        strng=latt(lp)
        if(strng(1:1).eq.'>')then
ccc       msglen=len_trim(lattmsg(num(lp)))
ccc       write(6,*)linestar(1:min(msglen,80))
ccc       write(6,*)trim(lattmsg(num(lp)))
ccc       write(6,*)linestar(1:min(msglen,80))
c
c  msgout =1 (>), =2 (>>), or =3 (>>>)
          msgout=1
          if(strng(1:2).eq.'>>')msgout=2
          if(strng(1:3).eq.'>>>')msgout=3
          msglen=len_trim(lattmsg(num(lp)))-msgout !length w/out > or >> or >>>
          if(msgout.eq.1 .or. msgout.eq.3)then
c write to terminal:
           if(idproc.eq.0)then
           write(jof,*)linestar(1:min(msglen,80))
           write(jof,*)lattmsg(num(lp))(msgout+1:msglen+msgout)
           write(jof,*)linestar(1:min(msglen,80))
           endif
          endif
          if(msgout.eq.2 .or. msgout.eq.3)then
c write to output file:
           if(idproc.eq.0)then
           write(jodf,*)linestar(1:min(msglen,80))
           write(jodf,*)lattmsg(num(lp))(msgout+1:msglen+msgout)
           write(jodf,*)linestar(1:min(msglen,80))
           endif
          endif
          lp=lp+1
          goto 4
        endif
c
        call lookup(latt(lp),ntype,ith)
        if (ntype.eq.1) then
c  command 'circ' ...
          if(nt1(ith).eq.7 .and. nt2(ith).eq.7 ) then
            icfile=nint(pmenu(1+mpp(ith)))
            nfcfle=nint(pmenu(2+mpp(ith)))
            norder=nint(pmenu(3+mpp(ith)))
            ntimes=nint(pmenu(4+mpp(ith)))
            nwrite=nint(pmenu(5+mpp(ith)))
            isend =nint(pmenu(6+mpp(ith)))
            jfctmp=jfcf
            jfcf  =nfcfle
            call cqlate(icfile,norder,ntimes,nwrite,isend)
            jfcj=jfctmp
c  command 'contractenv'
          elseif(nt1(ith).eq.7 .and. nt2(ith).eq.68) then
            write(6,*)'**********************************************'
            write(6,*)'**********************************************'
            write(6,*)'*********IN TRAN; found matchenv**************'
            write(6,*)'**********************************************'
            write(6,*)'**********************************************'
            niter=nint(pmenu(1+mpp(ith)))
            tolerance=pmenu(2+mpp(ith))
            strng=cmenu(1+mppc(ith))
c           write(6,*)'niter,tolerance,strng=',niter,tolerance,strng
c           note: this assumes that env(:) has been set. should check!
            envold(:)=env(:)
c           check that tracking variables nenvtrk,lautotrk,ntrktype make sense:
            nenvtrkold=nenvtrk
            lautotrkold=lautotrk
            ntrktypeold=ntrktype
            if(nenvtrk.ne.1)then
              nenvtrk=1
              if(idproc.eq.0)then
              write(6,*)'Envelope tracking turned on for rms matching.'
              write(6,*)'It will be turned off when finished matching.'
              endif
            endif
            if(lautotrk.ne.1)then
              lautotrk=1  !confusing; this really means "don't concatenate"
              if(idproc.eq.0)then
              write(6,*)'autoconcatenation turned off for rms matching.'
              write(6,*)'it will be turned on when finished matching.'
              endif
            endif
            if(ntrktype.ne.0)then
              ntrktype=0  !this, plus lautotrk=1,  means "don't call trace"
              if(idproc.eq.0)then
              write(6,*)'particle tracking turned off for rms matching.'
              write(6,*)'it will be turned on when finished matching.'
              endif
            endif
c start of do loop for rms matching:
            do iii=1,niter
c             store reference energy, etc.
              reftmp(1:6)=reftraj(1:6)
              arctmp=arclen
              brhotmp=brho
              gammatmp=gamma
              gamm1tmp=gamm1
              betatmp=beta
              call trobj(strng,1,0)
              call contractenv(delta)
c             check for convergence; if converged, break.
              if(idproc.eq.0)write(6,*)                                    &
     &        'rms matching: iteration,delta=',iii,delta
              if(delta.lt.tolerance)then
                if(idproc.eq.0)write(6,*)'SEARCH CONVERGED'
                exit
              endif
c             if not converged, restore reference energy, etc. and keep looping
              reftraj(1:6)=reftmp(1:6)
              arclen=arctmp
              brho=brhotmp
              gamma=gammatmp
              gamm1=gamm1tmp
              beta=betatmp
            enddo
            nenvtrk=nenvtrkold
            lautotrk=lautotrkold
            ntrktype=ntrktypeold
          else
c  ... or other element/command
            call trlmnt(ith,num(lp))
          endif
        else if (ntype.eq.2) then
c  line
          call trobj(latt(lp),num(lp),0)
        else if (ntype.eq.3) then
c  lump
c
c  ignore or destroy lump with zero repetition number
        if( num(lp).eq.0) then
c  if lump is unmade, ignore it
        if( lmade(ith).eq.0 ) then
            write(jof, 510) ilbl(ith)
            write(jodf,510) ilbl(ith)
 510        format(1x,'unmade lump ',a16,                               &
     &      ' in #labor with rep no. = 0 ignored')
        else
c  if lump is made, destroy it
          do 30 k = 1,maxlum
           if(inuse(k).eq.ith) then
            inuse(k) = 0
            write(jof, 520) ilbl(ith)
            write(jodf,520) ilbl(ith)
 520        format(1x,'lump ',a16,                                      &
     &      ' in #labor with rep n. = 0 destroyed')
           endif
 30       continue
          lmade(ith) = 0
        endif
        endif
c
c   otherwise
        if( num(lp).ne.0) then
          call trobj(latt(lp),num(lp),0)
        endif
c
        else if (ntype.eq.4) then
c  loop
          if( num(lp).ne.0) nloop=ith
          call trloop(ilbl(ith),num(lp))
c  store unmade lumps and replace lump-number by number in core
          do 40 i=1,joy
            if(mim(i).ge.0.and.mim(i).le.5000) then
              jth=mim(i)
              if(lmade(jth).eq.0) then
                call mapmap(th,tmh,rh,rmh)
                call ident(th,tmh)
                call trobj(ilbl(jth),1,0)
c               call lumpit(jth)
                call mapmap(rh,rmh,th,tmh)
              endif
              mim(i)=lmade(jth)
            endif
  40      continue
        else
c unused label
          if(idproc.eq.0)then
          write(jof,610) latt(lp)
  610     format(1h ,'warning from tran: ',a16,' not found.')
          endif
        endif
c
c Further modifications required to produce home-made do loop.
c
c   5  continue
c
      lp=lp+1
c     write(6,*)'AT BOTTOM OF MAIN LABOR LOOP; lp=',lp
      if(lp .le. noble) go to 4
c
c End of modifications. AJD 20 August 98
c     write(6,*)'returning from tran'
c
      return
      end
c
c***********************************************************************
c
      subroutine trlmnt(ith,mrep)
c-----------------------------------------------------------------------
c handle single item (element/command) in the menu
c
c input: ith  index of the item in menu
c        mrep repetition factor
c
c Petra Schuett, Nov.9,1987
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      use beamdata, only : bcurr
      use parallel, only : idproc
      use ml_timer
      include 'impli.inc'
      include 'map.inc'
      include 'codes.inc'
      include 'previous.inc'
      logical ldoit
c
      common/rfxtra/zedge,gaplen,thetastore
      common/showme/iverbose
      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0
      common/poiblock1/nspchset,nsckick
      common/autotrk/lautotrk,ntrktype,ntrkorder
      common/envstuff/nenvtrk
      common/autolog/lautoap1,lautoap2,lautoap3,lrestrictauto
c
c     iverbose=2
cryne 8/31/2001
      kt1=nt1(ith)
      kt2=nt2(ith)
cryne 1/11/2005 new code to only pre/post apply after "physical" elements.
c     the following is just a first attempt at this:
      ldoit=.false.
      if(lrestrictauto.eq.1 .and.                                       &
     &(kt1.eq.1 .or. kt1.eq.2 .or. kt1.eq.4 .or. kt1.eq.5))ldoit=.true.
c
c------ Put ith in inmenu common variable ------------
cryne 8/31/2001 I don't know who put this statement here or
cryne what it is good for:
      inmenu = ith
c----------
cryne 8/31/2001
c check to see whether or not this is a thick element:
c     write(6,*)' === afro::trlmnt(ith,mrep) ==='
c     write(6,*)' ith =',ith,'  mrep =',mrep
c     write(6,*)' calling isthick from trlmnt with kt1,kt2=',kt1,kt2
      call isthick(kt1,kt2,ithick,thlen,ith)
c     write(6,*)' done: ithick, thlen=', ithick,thlen
      if(iverbose.eq.2)then
        if(idproc.eq.0)write(6,*)'lmnlbl(ith)=',lmnlbl(ith)
        if(idproc.eq.0)then
          write(6,*)'kt1,kt2,ithick,thlen=',kt1,kt2,ithick,thlen
        endif
      endif
c
cryne 08/25/2001
c autotrack if the user has so specified (w/ autotrack command)
c also, for now the code is set to autotrack even if the user has NOT
c specified, but if the user HAS in fact set the poisson parameters.
c This is consistent with the code further down in this routine,
c where slices are split in half on the basis of whether or not
c nspchset has been set to 1.
c It would perhaps be better to do on the basis of nsckick, not nspchset.
      ntrk=0
      if(lautotrk.eq.1)ntrk=1
      if(nspchset.eq.1 .and. nsckick.ne.0)ntrk=1  !???we should omit this. rdr
cryne 11/15/2003 don't let the user mistakenly autotrack if the
cryne current is nonzero but the poisson params have not been set:
cryne this only applies if the code is about to track through a
cryne "thick" beamline element
      if(ntrk.eq.1 .and. ithick.eq.1 .and. ntrktype.ne.0)then
        if(bcurr.ne.0 .and. nspchset.eq.0)then
          if(idproc.eq.0)then
            write(6,*)'error: the code is about to track particles,'
            write(6,*)'but the current is nonzero and the Poisson'
            write(6,*)'solver has not been specified.'
            write(6,*)'Use the poisson type code and re-run'
          endif
          call myexit
        endif
      endif
c-----------------------------------------------------
      nslices=1
      slfrac=1.
c if autoslicing, get the number of slices:
c THIS ASSUMES THAT, IF THE CODE IS RUNNING IN THE MODE WHERE THE # OF
c SLICES IF EXPLICITLY SPECIFIED FOR EACH THICK ELEMENT, THEN IT IS
c ALWAYS THE LAST PARAMETER IN THE LIST:
c     if(idproc.eq.0)write(6,*)'sliceprecedence=',sliceprecedence
c     if(idproc.eq.0)write(6,*)'slicetype=',slicetype
c     if(idproc.eq.0)write(6,*)'slicevalue=',slicevalue
      if(ithick.eq.1 .and. slicetype.ne.'none')then
        if(sliceprecedence.eq.'local')then
c         if(idproc.eq.0)write(6,*)'determining slices locally'
! determine the number of slices locally:
          if(slicetype.eq.'slices')then
            nn0=mpp(ith)+nrp(1,kt2)
            nslices=pmenu(nn0)
c           if(idproc.eq.0)write(6,*)'thlen,nslices=',thlen,nslices
            if(nslices.le.0)then
              write(6,*)'error: nslices .le. 0; nslices=',nslices
              call myexit
            endif
          else
            nn0=mpp(ith)+nrp(1,kt2)
            if(pmenu(nn0).le.0)then
              write(6,*)'something wrong; slices/meter .le. 0'
              write(6,*)'pmenu(nn0)=',pmenu(nn0)
              call myexit
            endif
!!!!!!!!!!! elementlength=pmenu(mpp(ith)+1)
            elementlength=thlen
            nslices=nint(elementlength/pmenu(nn0))
c           if(idproc.eq.0)write(6,*)'thlen,nslices=',thlen,nslices
            if(nslices.eq.0)then
              write(6,*)'computed value of nslices =0; resetting to 1'
              nslices=1
            endif
ccc         write(6,*)'elementlength,slices/m,nslices=',                  &
ccc  &      elementlength,pmenu(nn0),nslices
          endif
        else
! determine the number of slices from the globally set value:
c         if(idproc.eq.0)write(6,*)'determining slices globally'
          if(slicetype.eq.'slices')then
            nslices=nint(slicevalue)
ccc         write(6,*)'global_nslices=',nslices
          else
!!!!!!!!!!! elementlength=pmenu(mpp(ith)+1)
            elementlength=thlen
            nslices=nint(elementlength/slicevalue)
            if(nslices.eq.0)then
              write(6,*)'computed value of nslices =0; resetting to 1'
              nslices=1
            endif
ccc         write(6,*)'elementlength,global_slices/m,nslices=',           &
ccc  &      elementlength,slicevalue,nslices
          endif
        endif
      endif
      kspchset=0
      if(nspchset.eq.1 .or. nenvtrk.eq.1)kspchset=1
      slfrac=1.d0/nslices
      if(lautoap2.eq.1 .or. kspchset.eq.1)slfrac=slfrac*0.5d0
c-----------------------------------------------------
c this is a repeat counter that is usually equal to one:
      do 20 i = 1,iabs(mrep)
c      write(6,991)lmnlbl(ith),nslices,slfrac,ntrk,arclen
       if(iverbose.eq.1)write(6,991)lmnlbl(ith),nslices,ntrk,arclen
       if(iverbose.eq.2)write(6,992)                                    &
     & lmnlbl(ith),(pmenu(m),m=1+mpp(ith),nrp(kt1,kt2)+mpp(ith))
c main loop: "tlmt,spch,lmnt" analogous to "map1,map2,map1"
c------------------
       nn1=1+mpp(ith)
       nn2=1+mppc(ith)
       do 15 j=1,nslices
       call init_step_timer
       call increment_timer('step',0)
cryne   if(lautoap1.eq.1)call autoapp(1)
        if(lautoap1.eq.1 .and. (kt1.ne.7.or.kt2.ne.42) .and. ldoit)       &
     &     call autoapp(1)
        if((ithick.eq.1).and.((lautoap2.eq.1).or.(kspchset.eq.1)))then
        call lmnt(kt1,kt2,pmenu(nn1),cmenu(nn2),ntrk,j,nslices,slfrac,1)
cryne8/21 if(lautoap2.eq.1)call autoapp(2)
          if(lautoap2.eq.1 .and. (kt1.ne.7.or.kt2.ne.43) .and. ldoit)     &
     &       call autoapp(2)
          tau=2.d0*prevlen
          if(kspchset.eq.1)call autospch(tau,ntrk)
        call lmnt(kt1,kt2,pmenu(nn1),cmenu(nn2),ntrk,j,nslices,slfrac,2)
        else
        call lmnt(kt1,kt2,pmenu(nn1),cmenu(nn2),ntrk,j,nslices,slfrac,0)
        endif
cryne821if(lautoap3.eq.1)call autoapp(3)
        if(lautoap3.eq.1 .and. (kt1.ne.7.or.kt2.ne.44) .and. ldoit)       &
     &     call autoapp(3)
! Short range wakefield forces
       if(cmenu(nn2).eq.'wake')then
         if(idproc.eq.0)then
         write(6,*)'short-range wakes temporarly commented out'
         write(6,*)'for debugging. RDR'
         endif
cryne    call wkfld_srange(nslices,lmnlbl(ith),tau)
       endif
       call increment_timer('step',1)
       call increment_timer('step',1)
       call step_timer(lmnlbl(ith),iverbose)
 15    continue
! Long range wakefield forces
       if(cmenu(nn2).eq.'wake')then
         if(idproc.eq.0)then
         write(6,*)'long-range wakes temporarly commented out'
         write(6,*)'for debugging. RDR'
         endif
cryne    call wkfld_lrange
       endif
c------------------
 20   continue
  991 format(a16,' nslices=',i5,' slfrac=',1pe14.7,' ntrk=',i1,         &
     &              ' arclen_in=',1pe14.7)
  992 format(a16,1x,9(1pe10.3,1x))
      return
      end
c
      subroutine isthick(kt1,kt2,ithick,thlen,ith)
c This routine determines if an element of type "kt1,kt2" is thick.
c If so, and if ith .ne.0, then it also returns its thickness.
c The routine is used in the autoslicing capability of MaryLie/IMPACT.
c It is necessary to return the length only if the user has requested
c that slicing be done N times per meter, i.e. the code needs to
c determine the number of slices based on the length.
c If the user specifies the number of slices, that all that is
c required is to determine whether the element is thick or not.
c Robert Ryne Nov. 30, 2002.
      use beamdata
      use acceldata
      integer kt1,kt2,ithick,ith
      integer itfile,numdata
      real*8 thlen,p1,p2,rho,gaplen
      character*16 estrng,fname,fname1
      character*5 aseq
      ithick=0
      thlen=0.
      if(kt1.ne.1)return
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
c          'rfgap   ','confoc  ','transit','interface','rootmap',
     &      1135,      1136,      133,       134,       135,
c          'optirot ','spare5  ','spare6  ','spare7  ','spare8  ',
     &      136,       137,       138,       139,       140,
c          'marker  ','drift   ','rbend   ','sbend   ','gbend   ',
     &      141,       142,       143,       1045,      145,
c          'quadrupo','sextupol','octupole','multipol','solenoid',
     &      146,       147,       148,       149,       150,
c          'hkicker ','vkicker ','kicker  ','rfcavity','elsepara',
     &      151,       152,       153,       154,       155,
c          'hmonitor','vmonitor','monitor ','instrume','sparem1 ',
     &      156,       157,       158,       159,       160,
c          'rcollima','ecollima','yrot    ','srot    ','sparem2 ',
     &      161,       162,       163,       164,       165,
c          'beambeam','matrix  ','profile1d','yprofile','tprofile',
     &      166,       167,       168,       169,       170,
c          'hkick   ','vkick   ','kick    ','hpm     ','nlrf    '/
     &      171,       172,       173,       174,       175),kt2

!drft 1:
  101 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
!     write(6,*)'isthick: drft length=',thlen
      return
!nbnd 2:
  102 continue
      ithick=1
      if(ith.eq.0)return
      rho=brho/pmenu(4+mpp(ith))
      p1=pmenu(1+mpp(ith))
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!pbnd 3:
  103 continue
      ithick=1
      if(ith.eq.0)return
      rho=brho/pmenu(4+mpp(ith))
      p1=pmenu(1+mpp(ith))
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!gbnd 4:
  104 continue
      ithick=1
      if(ith.eq.0)return
      rho=brho/pmenu(6+mpp(ith))
      p1=pmenu(1+mpp(ith))
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!prot 5:
  105 continue
      return
!gbdy 6:
  106 continue
      ithick=1
      if(ith.eq.0)return
      rho=brho/pmenu(4+mpp(ith))
      p1=pmenu(1+mpp(ith))
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!frng 7:
  107 continue
      return
!cfbd 8:
  108 continue
      ithick=1
      if(ith.eq.0)return
      rho=brho/pmenu(2+mpp(ith))
      p1=pmenu(1+mpp(ith))
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!quad 9:
  109 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
!     write(6,*)'isthick: quad length=',thlen
      return
!sext 10:
  110 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!octm 11:
  111 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!octe 12:
  112 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!srfc 13:
  113 continue
      return
!arot 14:
  114 continue
      return
!twsm 15:
  115 continue
      return
!thlm 16:
  116 continue
      return
!cplm 17:
  117 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!cfqd 18:
  118 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!dism 19:
  119 continue
      return
!sol  20:
  120 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      thlen=p2-p1
      return
!mark 21:
  121 continue
      return
!jmap 22:
  122 continue
      return
!dp   23:
  123 continue
      return
!recm 24:
  124 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      thlen=p2-p1
      return
!spce 25:
  125 continue
      return
!cfrn 26:
  126 continue
      return
!coil 27:
  127 continue
      return
!intg 28:
  128 continue
      return
!rmap 29:
  129 continue
      return
!arc  30:
  130 continue
      return
!rfgap 31:
 1135 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      if(p1.ge.0)then
        thlen=p1
        return
      endif
      if(p1.lt.0)then
        nseq=pmenu(5+mpp(ith))
        if(nseq.eq.0 .and. cmenu(1+mppc(ith)).eq.' ')then
          if(idproc.eq.0)write(6,*)'ISTHICK ERROR(rfgap):no data file'
          call myexit
        endif
        if(cmenu(1+mppc(ith)).eq.' ')then
          fname1='rfdata'
          ndigits=nseq/10+1
          call num2string(nseq,aseq,ndigits)
          j=len_trim(fname1)
          fname=fname1(1:j)//aseq(1:ndigits)
        else
          fname=cmenu(1+mppc(ith))
        endif
        nunit=0
        estrng='rfgap'
        call fnamechk(fname,nunit,ierr,estrng)
        if(ierr.eq.1)then
          write(6,*)'(isthick/rfgap) exiting due to problem w/ fname'
          call myexit
        endif
c       write(6,*)'calling read_RFdata from isthick'
        call read_RFdata(nunit,numdata,gaplen)
        thlen=gaplen
!       write(6,*)'isthick: rf gap length=',thlen
        return
      endif
!
!confoc 32:
 1136 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!transit:
  133 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!interface:
  134 continue
      return
!rootmap:
  135 continue
      return
!optirot:
  136 continue
      return
!spare5:
  137 continue
      return
!spare6:
  138 continue
      return
!spare7:
  139 continue
      return
!spare8:
  140 continue
      return
!marker:
  141 continue
      return
!
!drift 42:
  142 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
!     write(6,*)'isthick: drift length=',thlen
      return
!rbend 43:
  143 continue
cryne 12/21/2004
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      rho=brho/p2
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!sbend: 44
 1045 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      rho=brho/p2
      thlen=p1*rho*asin(1.d0)/90.d0
      return
!gbend: 45
  145 continue
      return
!quadrupo 46:
  146 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
!     write(6,*)'isthick: quadrupole length=',thlen
      return
!sextupol 47:
  147 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!octupole 48:
  148 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!multipol 49:
  149 continue
      return
!solenoid 50:
  150 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      thlen=p2-p1
      return
!hkicker  51:
  151 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!vkicker  52:
  152 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!kicker   53:
  153 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!rfcavity   54:
  154 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!elsepara   55:
  155 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!hmonitor 56
  156 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!vmonitor 57
  157 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!monitor  58
  158 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!instrume  59
  159 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!sparem1 60
  160 continue
      return
!rcollima 61
  161 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!ecollima 62
  162 continue
      ithick=1
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      thlen=p1
      return
!yrot
  163 continue
      return
!srot
  164 continue
      return
!sparem2
  165 continue
      return
!beambeam
  166 continue
      return
!matrix
  167 continue
      return
!profile1d 68
  168 continue
      return
!yprofile 69
  169 continue
      return
!tprofile 70
  170 continue
      return
!hkick 71
  171 continue
      goto 151
!vkick 72
  172 continue
      goto 152
!kick 73
  173 continue
      goto 153
!hpm 74
!fix hpm later
  174 continue
      return
!nlrf 75
  175 continue
      ithick=1
c     write(6,*) " === afro::isthick ==="
c     write(6,*) " ith    = ",ith
c     write(6,*) " ithick = ",ithick
      if(ith.eq.0)return
      p1=pmenu(1+mpp(ith))
      p2=pmenu(2+mpp(ith))
      thlen=p2-p1
      return
      end
c
c***********************************************************************
c
      subroutine trloop(string,nrept)
c-----------------------------------------------------------------------
c     interprets loops
c     Written by Rob Ryne ca 1984
c     adapted to use of strings and slightly simplified in logic
c     by Petra Schuett  11-10-87
c     Fixed by Filippo Neri and Alex Dragt 8-29-88
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'files.inc'
      include 'loop.inc'
c---------------
c parameter type
c---------------
      character string*(*)
c-----------------
c local variables
c-----------------
c lempty = .true. , when stack is empty
      logical lempty
      save lempty
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c???????????????????????????????????????????????????????????????????????
c       printout counters: prints out first nm icons
c
c      nm = 100
c      nc = 0
c      write(jof, 700)
c      write(jodf,700)
c 700  format(/' entering trloop:'/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c procedure when nrept=0
c
      if(nrept.eq.0) then
          write(jodf,508) string
          write(jof, 508) string
 508      format(1x,'loop ',a16,' with rep no. = 0 has been ignored')
      return
      endif
c
c procedure when nrept .ne. 0
c
      lnrept=npm1(nrept)
      joy=1
c  initialize the stacks
      call initst(string,lnrept)
c
c-----------------------------------------------------------------------
c here we start with a stack element, either a new one or one that
c just has been popped
c
 1000 continue
c
      if(ntype .eq. 1) then
c menu entry should never occur here ...
        write(jodf,910) lstac(np)
        write(jof ,910) lstac(np)
  910   format(1x,'error in trloop: menu entry ',a16,
     &            ' found at start of routine.')
        call myexit
      else if (ntype.eq.3) then
c ... neither should a lump
        write(jodf,920) lstac(np)
        write(jof ,920) lstac(np)
  920   format(1x,'error in trloop: lump ',a16,
     &            ' found at start of routine.')
        call myexit
      else if (ntype.eq.5) then
c unknown label is ignored
        write(jodf,930) lstac(np)
        write(jof, 930) lstac(np)
  930   format(1x,'warning from trloop: object ',a16,' not found.')
c
c main part is loops or lines (should be the only part used)
c
      else if (ntype.eq.2 .or. ntype.eq.4) then
c
c here we begin a computational loop, handling simple
c entries in lines or loops
c
 2000   continue
c
        if((ntype.eq.2 .or. ntype.eq.4) .and.
     &     (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a line or loop
          loop(np) = loop(np) - npm1(loop(np))
          if(loop(np).ne.0) then
            nslot(np) = newsl(np,ith)
            call lookup(icon(nslot(np),ith),mtype,jth)
            goto 1000
          endif
        else if (mtype.eq.2 .and. irep(nslot(np),ith).eq.0) then
c if  a line and rep no. = 0 ignore line
          write(jodf,509) icon(nslot(np),ith)
          write(jof, 509) icon(nslot(np),ith)
 509      format(1x,'line ',a16,' with rep no. = 0 has been ignored')
          nslot(np) = nslot(np) + npm1(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.4 .and. irep(nslot(np),ith).eq.0) then
c if a loop and rep no. = 0 ignore loop
          write(jodf,510) icon(nslot(np),ith)
          write(jof, 510) icon(nslot(np),ith)
 510      format(1x,'loop ',a16,' with rep no. = 0 has been ignored')
          nslot(np) = nslot(np) + npm1(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.5) then
c ignore unknown label
          write(jodf,930) icon(nslot(np),ith)
          write(jof, 930) icon(nslot(np),ith)
          nslot(np) = nslot(np) + npm1(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.1 .or. mtype.eq.3) then
c line or loop points to a menu entry or lump
          do 10 i=1,iabs(irep(nslot(np),ith))
            if(mtype.eq.1 .and. nt1(jth).eq.2) then
c user supplied element
              if(nt2(jth).gt.5) then
                write(jof ,940) icon(nslot(np),ith)
                write(jodf,940) icon(nslot(np),ith)
 940            format(/' warning from trloop: ',a16,' is a usern',
     &                  ' with n>5')
              endif
              mim(joy) = 5000+jth
            elseif (mtype.eq.1) then
c other menu entry
              mim(joy)=-jth
            elseif (mtype.eq.3) then
c lump
              mim(joy)=jth
            endif
            joy=joy+1
  10      continue
          if(joy.gt.joymax)then
            write(jodf,950) joymax
            write(jof ,950) joymax
 950        format(1x,'error: array length >= joymax (',                &
     &      i6,') in trloop')
            call myexit
          endif
c - next item:
          nslot(np) = nslot(np) + npm1(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.2 .or. mtype.eq.4) then
c line or loop points to line or loop
          call push
          goto 1000
        endif
c
c end of interpretation of a line/loop
c
      endif
c end of main way through stack-element
c-----------------------------------------------------------------------
c     pop stack; see what's there. this is the normal way to return
c
      call pop(lempty)
      if (lempty) then
        joy = joy-1
        return
      endif
      goto 1000
c
      end
c
c***********************************************************************
c
      subroutine trobj(string,nrept,ntrk)
c-----------------------------------------------------------------------
c     interprets lines and lumps
c  Written by Rob Ryne ca 1984
c        Modified to store unmade lumps (tro5)
c
c        New features:
c          1. Stores and retrieves unmade lumps in lines
c          2. Stores and retrieves nested lumps
c          3. Ignores lumps with nrept = 0
c          4. Ignores lines with nrept = 0
c
c        Jim Howard   CDG   7-7-87
c
c    adapted to use of strings and slightly simplified in logic
c        Petra Schuett  11-9-87
c
c    modified by Alex Dragt 31 August 1988
c-----------------------------------------------------------------------
      use acceldata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'core.inc'
      include 'map.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      character*16 string
c-----------------
c local variables
c-----------------
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c making(np) = 1  while this lump is under construction
      dimension making(mstack)
c lempty = .true. , when stack is empty
      logical lempty
      save thsave,tmhsav,making,lempty
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c  initialize arrays
c
      do 5 k = 1,mstack
        making(k) = 0
  5   continue
      call initst(string,nrept)
c???????????????????????????????????????????????????????????????????????
c       printout counters: prints out first nm calls to lmnt
c
c      nm = 100
c      nc = 0
c      write(jof, 700)
c      write(jodf,700)
c 700  format(/' entering trobj:'/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c ignore incoming line with zero repetition number
c
      if(ntype.eq.2.and.nrept.eq.0) then
        write(jof , 510) ilbl(ith)
        write(jodf, 510) ilbl(ith)
 510    format(1x,'line ',a16,' with rep no. = 0 has been ignored')
        return
      endif
c
c-----------------------------------------------------------------------
c here we start with a stack element, either a new one or one that
c just has been popped
c
 1000 continue
c
      if(ntype .eq. 1) then
c menu element should never occur here ...
        write(jodf,910) lstac(np)
        write(jof ,910) lstac(np)
  910   format(1x,'error in trobj: menu element ',a16,                  &
     & ' found at beginning of routine.')
        call myexit
      else if (ntype.eq.4) then
c ... neither should a loop
        write(jodf,920) lstac(np)
        write(jof ,920) lstac(np)
  920   format(1x,'error in trobj: loop ',a16,                          &
     & ' found at beginning of routine.')
        call myexit
      else if (ntype.eq.5) then
c unknown label is ignored
        write(jodf,930) lstac(np)
        write(jof, 930) lstac(np)
  930   format(1x,'warning from trobj: object ',a16,' not found.')
c
c main part is lumps or lines (should be the only part used)
c
      else if (ntype.eq.3 .and. loop(np).eq.0) then
c ignore lump with rep factor 0
          write(jof, 520) ilbl(ith)
          write(jodf,520) ilbl(ith)
 520      format(1x,'lump ',a16,' with rep no. = 0 has been ignored')
      else if (ntype.eq.3 .and. lmade(ith).ne.0) then
c a lump which is already in the core
c???????????????????????????????????????????????????????????????????????
c        write(jof, 710) ith
c        write(jodf,710) ith
c 710    format(/' picking up old lump, no.',i4)
c        if(nc.lt.nm) write(jodf,715) jth,lmade(jth),loop(np),np
c        if(nc.lt.nm) write(jof, 715) jth,lmade(jth),loop(np),np
c 715    format(/5x,'jth =',i4,'    lmade(jth) =',i4,
c     &          5x,'lumprep =',i3,'  np =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - combine previously made lump with the total map
        call comwtm(lmade(ith),loop(np))
c
      else if (ntype.eq.2 .or. ntype.eq.3) then
c all other lines and lumps
c
c first, lumps need special handling, if they are new
        if (ntype.eq.3 .and. making(np).ne.1) then
          making(np) = 1
c???????????????????????????????????????????????????????????????????????
c          write(jodf,720) ith,np,loop(np)
c          write(jof, 720) ith,np,loop(np)
c 720      format(/' starting lump no. ',i3,
c     1            '  np =',i3,'  loop(np) =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - save current map and initialize a new lump
          call buffin(th,tmh,thsave,tmhsav,np)
          call ident(th,tmh)
        endif
c
c here we begin a computational loop, handling simple elements
c in lines or lumps
c
 2000   continue
c
        if(ntype.eq.2 .and.                                             &
     &     (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a line
c????????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,730)
c          if(nc.lt.nm) write(jof, 730)
c  730     format(/' *********************** end of line')
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          loop(np) = loop(np) - npm1(loop(np))
          if(loop(np).ne.0) then
            nslot(np) = newsl(np,ith)
            call lookup(icon(nslot(np),ith),mtype,jth)
            goto 1000
          endif
        else if (ntype.eq.3 .and.                                       &
     &           (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a lump
c????????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,740)
c          if(nc.lt.nm) write(jof, 740)
c  740     format(/' *********************** end of lump')
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - store lump
          call lumpit(ith)
c - recall running total map
          call bufout(thsave,tmhsav,th,tmh,np)
c - combine new lump with total map:
c???????????????????????????????????????????????????????????????????????
c          write(jof, 750) loop(np)
c          write(jodf,750) loop(np)
c 750      format(/'  combine new lump with total map: lumprep =',i3/)
c          write(jof, 755) ith,lmade(ith)
c          write(jodf,755) ith,lmade(ith)
c 755      format(/' ith =',i3,'   lmade =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call comwtm(lmade(ith),loop(np))
          making(np) = 0
        else if (mtype.eq.4) then
c lump or line points to a loop; error
          write(jof, 940) lstac(np),icon(nslot(np),ith)
          write(jodf,940) lstac(np),icon(nslot(np),ith)
 940      format(1x,'error in trobj: ',a16,                             &
     &    ' contains a loop (',a16,').')
          call myexit
        else if (mtype.eq.2 .and. irep(nslot(np),ith).eq.0) then
c if rep no. = 0 ignore line
          write(jof ,510) icon(nslot(np),ith)
          write(jodf,510) icon(nslot(np),ith)
          nslot(np) = nslot(np) + npm1(loop(np))
          call lookup(icon(nslot(np),ith),mtype,jth)
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,760) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 760) ith,mtype,jth
c 760      format(/' 0*line: ith, mtype,   jth:'/5x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.5) then
c ignore unknown label
          write(jodf,930) icon(nslot(np),ith)
          write(jof, 930) icon(nslot(np),ith)
          nslot(np) = nslot(np) + npm1(loop(np))
          call lookup(icon(nslot(np),ith),mtype,jth)
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,762) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 762) ith,mtype,jth
c 762      format(/' unkn. : ith, mtype,   jth:'/5x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.1) then
c line or lump points to an element
c???????????????????????????????????????????????????????????????????????
c          nc = nc + 1
c          if(nc.lt.nm) write(jodf,770) np
c          if(nc.lt.nm) write(jof, 770) np
c 770      format(/' element in line or lump: np =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call trlmnt(jth,irep(nslot(np),ith))
c - next element:
          nslot(np) = nslot(np) + npm1(loop(np))
cryne Jan 6, 2005 added "if" test:
          if(nslot(np).ge.1)then
            call lookup(icon(nslot(np),ith),mtype,jth)
          endif
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,764) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 764) ith,mtype,jth
c 764      format(/' next icon: ith, mtype,   jth:'/7x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.2 .or. mtype.eq.3) then
c line or lump points to line or lump
          call push
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jof, 780) np
c          if(nc.lt.nm) write(jodf,780) np
c780      format(/' *********************** pushing stack: new np =',i3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 1000
        endif
c
c end of interpretation of a line/lump
c
      endif
c end of main way through stack-element
c-----------------------------------------------------------------------
c     pop stack; see what's there. this is the normal way to return
c
      call pop(lempty)
      if (lempty) return
c
c???????????????????????????????????????????????????????????????????????
c      write(jof, 790) np
c      write(jodf,790) np
c790  format(/' *************************** popping stack: new np =',i3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      goto 1000
      end
c
c end of file
c
      subroutine old_set_pscale_mc(h,mh)
c scale a map to use mc as the scale momentum (instead of p0)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      double precision h(monoms),mh(6,6)
c
c This routine is needed if we want to use "dynamic" units, which
c makes sense, e.g., if the beam is being accelerated.
c It converts maps from the "usual" units (in which MaryLie lets the user
c specify any scale length and assumes a scale momentum p0)
c to units where
c THE REST OF THIS COMMENT IS WRONG. FIX LATER. RDR
cthe new scale length is given by l=c/w (where w is a scale
c frequency = 2pi*the freq specified by the user in the new #beam component)
c and where the new scale momentum is m*c
c
c Skip this routine if we are using magnetostatic units (the "usual" units)
c     write(6,*)'****************************************************'
c     write(6,*)'**************HERE I AM IN SETPSCALEMC**************'
c     write(6,*)'****************************************************'
c
c exit if "magnetic" units (i.e. standard MaryLie units) are being used:
      if(lflagmagu)return
c see if "dynamic" units are being used; change only affects the momenta
c since the MaryLie library routines include the effect of sl:
      if(lflagdynu)then
c       write(6,*)'dynamic units.'
        x2ovrx1=1.
        p2ovrp1=1./(beta*gamma)
      else
c some other units are being used:
c       write(6,*)'general units.'
        x2ovrx1=1.
        p2ovrp1=p0sc/(beta*gamma*pmass/c)
      endif

      x1ovrx2=1./x2ovrx1
      p1ovrp2=1./p2ovrp1
cdebug
c     write(6,*)'inside set_pscale'

c     bg=beta*gamma
c     bgi=1./bg
c
ctm      mh(1,2:6:2)=mh(1,2:6:2)/bg
ctm      mh(3,2:6:2)=mh(3,2:6:2)/bg
ctm      mh(5,2:6:2)=mh(5,2:6:2)/bg
c
ctm      mh(2,1:5:2)=mh(2,1:5:2)*bg
ctm      mh(4,1:5:2)=mh(4,1:5:2)*bg
ctm      mh(6,1:5:2)=mh(6,1:5:2)*bg
ctm replace f90 array code with f77 do loops:
crdr replace beta*gamma with more general form:
         do 100 i = 2, 6, 2
            do 70 j = 1, 5, 2
               mh(j,i) = mh(j,i)*p2ovrp1*x1ovrx2
               mh(i,j) = mh(i,j)*p1ovrp2*x2ovrx1
   70       continue
  100    continue
c
cryne 08/26/2001
c now change the units of the nonlinear part of the map:
      do 200 i=28,209
c     h(i)=h(i)*bgi**(expon(2,i)+expon(4,i)+expon(6,i))
c     h(i)=h(i)*bgi**exp246(i)
      h(i)=h(i)*(p2ovrp1**nxp246(i))*(x2ovrx1**nxp135(i))
  200 continue
      return
      end
c
      subroutine set_pscale_mc(h,mh)
c scale a map to use mc as the scale momentum (instead of p0)
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      double precision h(monoms),mh(6,6)
c
c     write(6,*)'=== afro::set_pscale_mc() ==='
c     write(6,*) 'beta, gamma, pmass, brho = ',beta,gamma,pmass,brho
c     write(6,*) 'scale l, p, t, w, f= ',sl,p0sc,ts,omegascl,freqscl
c
cryne Jan 30,2003
c this routine does 2 things:
c (1) convert from static to dynamic units (if needed),
c (2) convert variables 5 and 6 to take into account the fact that the
c     scale varables omega and l are independent and do not need to
c     satisfy omega*l/c=1 (as is assumed in the MaryLie library routines)
c
      clite=299792458.d0
c dynamic units:
calsowrong      p2ovrp1=1.d0
      p2ovrp1=p0sc/(gamma*beta*pmass/clite)
c     write(6,*)'initializing p2ovrp1 to ',p2ovrp1
c
      if(lflagdynu)then
        p2ovrp1=1./(beta*gamma)
c       write(6,*)'lflagdynu is true; resetting p2ovrp1 to ',p2ovrp1
      endif
cwrong:      if(lflagdynu)p2ovrp1=(pmass/clite)/p0sc
cryne the reason this was wrong is that the MaryLie library
cryne routines use brho, the p0sc.In other words, everything
cryne is computed relative to the present momentum.
      p1ovrp2=1.d0/p2ovrp1
c omega*l/c:
      scal5=1.d0
      scal6=1.d0
      wlbyc=omegascl*sl/clite
c     write(6,*)'initializing wlbyc to ',wlbyc
      if(wlbyc.ne.1.d0)then
c       write(6,*)'scal5 being set to wlbyc, scal6 set to 1/wlbyc'
        scal5=wlbyc
        scal6=1.d0/wlbyc
      endif
c
      mh(1,2)=mh(1,2)*p2ovrp1
      mh(1,4)=mh(1,4)*p2ovrp1
      mh(1,5)=mh(1,5)/wlbyc
      mh(1,6)=mh(1,6)*p2ovrp1*wlbyc
c
      mh(2,1)=mh(2,1)*p1ovrp2
      mh(2,3)=mh(2,3)*p1ovrp2
      mh(2,5)=mh(2,5)*p1ovrp2/wlbyc
      mh(2,6)=mh(2,6)*wlbyc
c
      mh(3,2)=mh(3,2)*p2ovrp1
      mh(3,4)=mh(3,4)*p2ovrp1
      mh(3,5)=mh(3,5)/wlbyc
      mh(3,6)=mh(3,6)*p2ovrp1*wlbyc
c
      mh(4,1)=mh(4,1)*p1ovrp2
      mh(4,3)=mh(4,3)*p1ovrp2
      mh(4,5)=mh(4,5)*p1ovrp2/wlbyc
      mh(4,6)=mh(4,6)*wlbyc
c
      mh(5,1)=mh(5,1)*wlbyc
      mh(5,2)=mh(5,2)*p2ovrp1*wlbyc
      mh(5,3)=mh(5,3)*wlbyc
      mh(5,4)=mh(5,4)*p2ovrp1*wlbyc
      mh(5,6)=mh(5,6)*p2ovrp1*wlbyc**2
c
      mh(6,1)=mh(6,1)*p1ovrp2/wlbyc
      mh(6,2)=mh(6,2)/wlbyc
      mh(6,3)=mh(6,3)*p1ovrp2/wlbyc
      mh(6,4)=mh(6,4)/wlbyc
      mh(6,5)=mh(6,5)*p1ovrp2/wlbyc**2
c
c now change the units of the nonlinear part of the map:
cryne 12/21/2004 changed "i=28,209" to "i=28,monoms"  !!!!!!!!!!!!!!
c     write(6,*)'modified set_pscale_mc do loop to use monoms=',monoms
      do 200 i=28,monoms
ccc   h(i)=h(i)*(p2ovrp1**nxp246(i))*(x2ovrx1**nxp135(i))
!!!   h(i)=h(i)*(p2ovrp1**nxp24(i))*(scal5**nxp5(i))*(scal6**nxp6(i))
cryne 12/21/2004 note that this is a minus 1 only in nxp13 and nxp24
cryne note that we are assuming that the scale length has been done already
!     h(i)=h(i)*(p2ovrp1**(nxp24(i)+nxp6(i)))*(scal5**(nxp5(i)-nxp6(i)))
      h(i)=h(i)*(p2ovrp1**(nxp24(i)+nxp6(i)))*(scal5**(nxp6(i)-nxp5(i)))
  200 continue
      return
      end
c
      subroutine set_rfscale(h,mh)
c Nov 5, 2003
c scale the map produced by subroutine rfgap (if needed)
c note that the rf gap routine returns a map assuming
c dynamic units, with a scale length l=c/omega
      use beamdata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      double precision h(monoms),mh(6,6)
c
c switch to static units if needed:
c rf gap routine works in "dynamic" units, with p1=mc;
c If MLI is using "static" units, this routine transforms
c        to p2=the scale momentum p0sc
c therefore p2ovrp1=p0sc/mc if shifting to static mode
      clite=299792458.d0
      p2ovrp1=1.d0
      if(.not.lflagdynu)p2ovrp1=p0sc/(pmass/clite)
      p1ovrp2=1.d0/p2ovrp1
c
c Use this if the rfgap routine works internally assuming q1=c/omega,
c with the scale ang freq (internally) equal to the MLI scale ang freq:
c In this case, q2 is the MLI scale length, therefore q2ovrq1=sl/(c/omega)
cryne 5/1/2006 uncommented these:
      q2ovrq1=sl/(clite/omegascl)
      q1ovrq2=1.d0/q2ovrq1
      w2ovrw1=1.d0
      w1ovrw2=1.d0/w2ovrw1
c
c Use this if the rfgap routine works internally assuming a
c scale ang freq given by omega=c/sl, with the scale length
c (internally) equal to the MLI scale length: In this case
c w2 is the MLI scale ang freq, therefore w2ovrw1=omegascl/(c/sl)
cryne 5/1/2006 commented out the following:
c     w2ovrw1=omegascl/(clite/sl)
c     w1ovrw2=1.d0/w2ovrw1
c     q2ovrq1=1.d0
c     q1ovrq2=1.d0/q2ovrq1
c
c Scale linear part of the map.
      mh(1,2)=mh(1,2)*p2ovrp1*q1ovrq2
      mh(1,4)=mh(1,4)*p2ovrp1*q1ovrq2
      mh(1,5)=mh(1,5)*w1ovrw2*q1ovrq2
      mh(1,6)=mh(1,6)*p2ovrp1*w2ovrw1
c
      mh(2,1)=mh(2,1)*p1ovrp2*q2ovrq1
      mh(2,3)=mh(2,3)*p1ovrp2*q2ovrq1
      mh(2,5)=mh(2,5)*p1ovrp2*w1ovrw2
      mh(2,6)=mh(2,6)*w2ovrw1*q2ovrq1
c
      mh(3,2)=mh(3,2)*p2ovrp1*q1ovrq2
      mh(3,4)=mh(3,4)*p2ovrp1*q1ovrq2
      mh(3,5)=mh(3,5)*w1ovrw2*q1ovrq2
      mh(3,6)=mh(3,6)*p2ovrp1*w2ovrw1
c
      mh(4,1)=mh(4,1)*p1ovrp2*q2ovrq1
      mh(4,3)=mh(4,3)*p1ovrp2*q2ovrq1
      mh(4,5)=mh(4,5)*p1ovrp2*w1ovrw2
      mh(4,6)=mh(4,6)*w2ovrw1*q2ovrq1
c
      mh(5,1)=mh(5,1)*w2ovrw1*q2ovrq1
      mh(5,2)=mh(5,2)*w2ovrw1*p2ovrp1
      mh(5,3)=mh(5,3)*w2ovrw1*q2ovrq1
      mh(5,4)=mh(5,4)*w2ovrw1*p2ovrp1
      mh(5,6)=mh(5,6)*q2ovrq1*p2ovrp1*w2ovrw1**2
c
      mh(6,1)=mh(6,1)*w1ovrw2*p1ovrp2
      mh(6,2)=mh(6,2)*w1ovrw2*q1ovrq2
      mh(6,3)=mh(6,3)*w1ovrw2*p1ovrp2
      mh(6,4)=mh(6,4)*w1ovrw2*q1ovrq2
      mh(6,5)=mh(6,5)*q1ovrq2*p1ovrp2*w1ovrw2**2
c
cryne 08/26/2001
cryne 11/06/2003
c now change the units of the nonlinear part of the map:
      scal6=w2ovrw1*q2ovrq1*p2ovrp1
cryne 12/21/2004 changed "i=28,209" to "i=28,monoms"
cryne (this does not matter now, but will matter when we do nonlinear rf maps)
c     write(6,*)'(subroutine set_rfscale): need to fix scaling of the'
c     write(6,*)'nonlinear part of rfgap map; it is probably wrong!!'
cabell =Tue Dec 28 08:32:15 PST 2004= Think I've fixed scaling here.
cccryne 22 May 2006 -- Since this is only called by rfgap (i.e. only the linear
ccc                    map is produced), skip the scaling of the nonlinear part
ccc   write(6,*)'(subroutine set_rfscale): still need to verify scaling'
ccc   write(6,*)'of nonlinear part of rfcavity!! **********************'
ccc   do 200 i=28,monoms
ccc     h(i)=h(i)*(q2ovrq1**(nxp13(i)+nxp6(i)))                         &
ccc  &           *(p2ovrp1**(nxp24(i)+nxp6(i)))                         &
ccc  &           *(w2ovrw1**(nxp6(i)-nxp5(i)))
cabell this is equivalent to above
c       h(i)=h(i)*(p2ovrp1**nxp24(i))*(q2ovrq1**nxp13(i))               &
c    &           *(w1ovrw2**nxp5(i))*(scal6**nxp6(i))
cabell this is ok if scale lengths are equal
c       h(i)=h(i)*(p2ovrp1**(nxp24(i)+nxp6(i)))                         &
c    &           *(w2ovrw1**(nxp6(i)-nxp5(i)))
ccc200 continue
      return
      end
c

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c==============================================================================
c//////////////////////////////////////////////////////////////////////////////

      subroutine rescale_map(lsold,psold,fsold,h,mh,lsnew,psnew,fsnew)
c------------------------------------------------------------------------------
c
c  rescale_map:   Takes given values for the scale factors of a given map
c                 and rescales the map to the currently defined scale factors.
c                 Optionally, it can take a list of new scale factors, instead
c                 of assuming the currently defined scale factors.
c
c  input:
c
c     lsold -  given length scale of map
c     psold -  given momentum scale of map
c     fsold -  given frequency scale of map
c     h     -  non-linear coefficients of given map
c     mh    -  matrix part of given map
c
c  OPTIONAL input: (all must be present to work!!!)
c
c     lsnew -  length scale to scale map TO
c     psnew -  momentum scale to scale map TO
c     fsnew -  frequency scale to scale map TO
c
c  output:
c
c     h  -  rescaled map
c     mh -  rescaled map
c
c  KMP - 15 Dec 2006
c
c------------------------------------------------------------------------------
      use lieaparam, only : monoms
      use beamdata, only : sl,p0sc,freqscl
      implicit none
      include 'expon.inc'
      double precision, intent(in) :: lsold,psold,fsold ! scale factors of given map
      double precision, optional, intent(in) :: lsnew,psnew,fsnew ! new scale factors to use in scaling
      double precision, intent(inout), dimension(1:6,1:6) :: mh
      double precision, intent(inout), dimension(1:monoms) :: h
      double precision p2ovrp1,p1ovrp2
      double precision q2ovrq1,q1ovrq2
      double precision w2ovrw1,w1ovrw2
      integer i,j
c
      if(psold.eq.0.)then
         write(*,*) 'error (rescale_map): zero momentum scale given'
         call myexit
      elseif(lsold.eq.0.)then
         write(*,*) 'error (rescale_map): zero length scale given'
         call myexit
      elseif(fsold.eq.0.)then
         write(*,*) 'error (rescale_map): zero frequency scale given'
         call myexit
      endif

      if(present(lsnew).and.present(psnew).and.(present(fsnew)))then
         if(psnew.eq.0.)then
            write(*,*) 'error (rescale_map): zero momentum scale given'
            call myexit
         elseif(lsnew.eq.0.)then
            write(*,*) 'error (rescale_map): zero length scale given'
            call myexit
         elseif(fsnew.eq.0.)then
            write(*,*) 'error (rescale_map): zero frequency scale given'
            call myexit
         endif
         p1ovrp2 = psold / psnew
         p2ovrp1 = psnew / psold
         q1ovrq2 = lsold / lsnew
         q2ovrq1 = lsnew / lsold
         w1ovrw2 = fsold / fsnew
         w2ovrw1 = fsnew / fsold
      else
         p1ovrp2 = psold / p0sc
         p2ovrp1 = p0sc / psold
         q1ovrq2 = lsold / sl
         q2ovrq1 = sl / lsold
         w1ovrw2 = fsold / freqscl
         w2ovrw1 = freqscl / fsold
      endif
c
      mh(1,2)=mh(1,2)*p2ovrp1*q1ovrq2
      mh(1,4)=mh(1,4)*p2ovrp1*q1ovrq2
      mh(1,5)=mh(1,5)*w1ovrw2*q1ovrq2
      mh(1,6)=mh(1,6)*p2ovrp1*w2ovrw1
c
      mh(2,1)=mh(2,1)*p1ovrp2*q2ovrq1
      mh(2,3)=mh(2,3)*p1ovrp2*q2ovrq1
      mh(2,5)=mh(2,5)*p1ovrp2*w1ovrw2
      mh(2,6)=mh(2,6)*w2ovrw1*q2ovrq1
c
      mh(3,2)=mh(3,2)*p2ovrp1*q1ovrq2
      mh(3,4)=mh(3,4)*p2ovrp1*q1ovrq2
      mh(3,5)=mh(3,5)*w1ovrw2*q1ovrq2
      mh(3,6)=mh(3,6)*p2ovrp1*w2ovrw1
c
      mh(4,1)=mh(4,1)*p1ovrp2*q2ovrq1
      mh(4,3)=mh(4,3)*p1ovrp2*q2ovrq1
      mh(4,5)=mh(4,5)*p1ovrp2*w1ovrw2
      mh(4,6)=mh(4,6)*w2ovrw1*q2ovrq1
c
      mh(5,1)=mh(5,1)*w2ovrw1*q2ovrq1
      mh(5,2)=mh(5,2)*w2ovrw1*p2ovrp1
      mh(5,3)=mh(5,3)*w2ovrw1*q2ovrq1
      mh(5,4)=mh(5,4)*w2ovrw1*p2ovrp1
      mh(5,6)=mh(5,6)*q2ovrq1*p2ovrp1*w2ovrw1**2
c
      mh(6,1)=mh(6,1)*w1ovrw2*p1ovrp2
      mh(6,2)=mh(6,2)*w1ovrw2*q1ovrq2
      mh(6,3)=mh(6,3)*w1ovrw2*p1ovrp2
      mh(6,4)=mh(6,4)*w1ovrw2*q1ovrq2
      mh(6,5)=mh(6,5)*q1ovrq2*p1ovrp2*w1ovrw2**2
c
      do i=28,monoms
         h(i)=h(i)*(q2ovrq1**(nxp13(i)+nxp6(i)))
     &            *(p2ovrp1**(nxp24(i)+nxp6(i)))
     &            *(w2ovrw1**(nxp6(i)-nxp5(i)))
      enddo
      return
      end


************************************************************************
      subroutine mlfinc(incnam)
c-----------------------------------------------------------------------
      use beamdata
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'codes.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c arguments:
      character incnam*(*)
c local variables:
      character*16 strarr(40),string
      character*100 line
      integer narr(40)
      logical leof,lcont,npound
cryne 8/15/2001  new code to read in a character string instead of numbers:
c-----------------------------------------------------------------------
c start
      lfold=lf
      lf=98
cryne----- 15 Sept 2000 modified to check for input file:
      open(lf,file=incnam,status='old',err=357)
      itot = 0
      leof = .false.
      call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
ctm      write(6,333) line
ctm 333  format(' first line in include file:',/,a)
      if(leof)then
        write(6,*)'error: empty file:',incnam
        stop
      endif
      goto 4320
ctm  357 write(6,*)'error: empty file:', incnam
c
ctm     error opening include file
c
  357 continue
      write(jof,358) incnam
      write(jodf,358) incnam
  358 format(' ERROR: cannot open #include file: ',/,a)
      call myexit
cryne-----
c-----------------------------------------------------------------------
 4320 continue
      leof = .false.
      rewind lf
      msegm = -1
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800),msegm
c
c error exit:
      write(jof ,99) msegm
      write(jodf,99) msegm
  99  format(1x,i6,'=msegm problem at big goto of routine mlfinc')
      call myexit
c--------------------
c  #comment
c
  100 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
c
      npcom=npcom+1
      if(npcom.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'error in dumpin:',                                   &
     & ' too many comment lines (>= maxcmt = ',i6,') in #comment')
        call myexit
      endif
      mline(npcom)=line
      goto 100
c--------------------
c  #beam
c
  200 continue
      write(6,*)'error: cannot have #beam when using mlfinc'
      stop
c--------------------
c  #menu
c
  300 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 3) goto 1
c
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
     & a8,' is doubly defined'/                                         &
     &         ' the second definition will be ignored')
        na = na-1
        goto 300
      endif
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a8,' has no type code!')
        call myexit
      endif
      if ( itot .gt. 2 ) then
        write(jof ,1397) strarr(1)
        write(jodf ,1397) strarr(1)
 1397   format(' error in #menu: ',a8,' has more than one type code!')
        call myexit
      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
c string is name of element/command type, look up element indices
      string=strarr(2)
      do 325 m = 1,9
cryne do 325 n = 1,40
cryne 7/3/02   do 325 n = 1,45
      do 325 n = 1,nrpmax
      if(string.eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        if(imax.eq.0)goto 300
        mpp(na) = mpprpoi
cryne---8/15/2001
        if((m.eq.1).and.(n.eq.42.or.n.eq.43.or.n.eq.44))then
c       read a character string
cryne 08/24/2001
        mpp1=mppc(na)+1
cryne   read(lf,*,err=390)(cmenu(1+mppc(na)))
        read(lf,*,err=390)cmenu(mpp1)
        pmenu(1+mpp(na))=1+mpp(na)
        mpprpoi=mpprpoi+1
        goto 300
        else
cryne---
        read(lf,*,err=390)(pmenu(i+mpp(na)),i=1,imax)
        mpprpoi = mpprpoi + imax
        goto 300
        endif
c normal end of a menu element/command
      endif
  325 continue
c
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin error in ',a8,': type code ',a8,' not found.'/ &
     &       1h ,'this item will be ignored')
      na=na-1
      goto 300
c
c error in parameter input
  390 write(jof ,397)lmnlbl(na)
      write(jodf,397)lmnlbl(na)
  397 format(1h ,'dumpin data input error at item ',a8)
      call myexit
c--------------------
c  #lines,#lumps,#loops
c
  400 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
      if (leof) goto 1000
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
     &  a8,': name ',a8,' is doubly defined'/                           &
     &         ' the second definition will be ignored')
        na = na-1
        if(msegm .eq.4) then
          goto 400
        else if(msegm.eq.5) then
          goto 500
        else if(msegm.eq.6) then
          goto 600
        endif
      endif
c new item
      ilbl(nb)=strarr(1)
c read components of item
      imin=0
c--
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
  700 call rearec(line,leof,msegm,strarr,narr,itot,lcont,npound)
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
c  #include
  800 continue
      write(6,*)'error: cannot use #include with mlfinc'
      stop
c normal return at end of file
 1000 continue
      close(98)
      lf=lfold
      return
      end
c
c***********************************************************************
      subroutine autoapp(n)
      use parallel, only : idproc
      use acceldata
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'map.inc'
      character*16 string
      character*16 autostr
      real*8 thlen
      common/autopnt/autostr(3)
      common/autolog/lautoap1,lautoap2,lautoap3,lrestrictauto
      common/showme/iverbose
c---dummy variables needed in calls to lmnt:
      jslice=1
      jsltot=1
      slfrac=1.
      ntrk=0
c---
c string is the value of parameter "name" in the autoapply command
      string=autostr(n)
      if(iverbose.eq.2)write(6,*)'(inside autoapp) n,string=',n,string
c this version works only if the object is a either menu element or a line
      call lookup(string,itype,ith)
      if(iverbose.eq.2)write(6,*)'back from lookup; string,itype,ith='
      if(iverbose.eq.2)write(6,*)string,itype,ith
cryne Oct 30, 2003
      if(itype.eq.5)then
        if(idproc.eq.0)then
          write(6,*)'error: trying to perform autoapply, but the item'
          write(6,*)string(1:16)
          write(6,*)'has not been successfuly defined in the menu'
        endif
        call myexit
      endif
      mt1=nt1(ith)
      mt2=nt2(ith)
      l1=1+mpp(ith)
      l2=1+mppc(ith)
c element:
      if(itype.eq.1)then
      if(iverbose.eq.2)write(6,*)'string ',string,' is a menu item'
c       don't allow automatic application of thick elements
        call isthick(mt1,mt2,ithick,thlen,0)
        if(ithick.eq.1)write(6,*)'error: thick element found!:',string
        if(ithick.eq.1)return
      call lmnt(mt1,mt2,pmenu(l1),cmenu(l2),ntrk,jslice,jsltot,slfrac,0)
        if(iverbose.eq.2)write(6,*)'normal return from autoapp'
        return
      endif
c step through the line:
      if(itype.eq.2)then
       if(iverbose.eq.2)write(6,*)'string ',string,' is a line'
       if(iverbose.eq.2)                                                &
     & write(6,*)'ith,ilbl(ith),ilen(ith)=',ith,ilbl(ith),ilen(ith)
       do 100 k=1,ilen(ith)
        call lookup(icon(k,ith),ktype,kth)
        if(ktype.ne.1)then
          if(idproc.eq.0)then
            write(6,*)'error in autoapp:',string,icon(k,ith)
            if(ktype.eq.2)write(6,*)'cannot handle nested lines'
            if(ktype.eq.3)write(6,*)'cannot handle lumps'
            if(ktype.eq.4)write(6,*)'cannot handle loops'
            if(ktype.eq.5)write(6,*)'cannot find entry: ',icon(k,ith)
          endif
          call myexit
        endif
        kt1=nt1(kth)
        kt2=nt2(kth)
        l1=1+mpp(kth)
        l2=1+mppc(kth)
        call isthick(kt1,kt2,ithick,thlen,0)
        if(ithick.eq.1)write(6,*)'error:thick lmnt in line!:',icon(k,2)
        if(ithick.eq.1)return
        do 99 iii=1,iabs(irep(k,ith))
        if(iverbose.eq.2)write(6,*)'icon(k,ith)=',icon(k,ith)
      call lmnt(kt1,kt2,pmenu(l1),cmenu(l2),ntrk,jslice,jsltot,slfrac,0)
   99  continue
  100  continue
       if(iverbose.eq.1)write(6,*)'normal return from autoapp'
       return
      endif
c can autoapply only menu entries or simple lines.
c complain that you cannot autoapply lumps or loops
      if(itype.eq.3 .or. itype.eq.4)then
        write(6,*)'autoapp: error: lump or loop not allowed:',string
        return
      endif
      if(itype.eq.5)then
        write(6,*)'autoapp:string not found:',string
        return
      endif
      end
c
************************************************************************
      subroutine autospch(tau,ntrk)
      use beamdata
      use rays
      use lieaparam, only : monoms
      use ml_timer
      include 'impli.inc'
      include 'map.inc'
      include 'files.inc'
      logical msk,straight
      character*16 cparams
!3/6/03      dimension c4(4,maxray),c6(6,maxray),msk(maxray)
!3/6/03      dimension ex(maxray),ey(maxray),ez(maxray)
!3/6/03      dimension tmp(maxray),tmpr(maxray)
      dimension c4(4,maxrayp),c6(6,maxrayp),msk(maxrayp)
      dimension ex(maxrayp),ey(maxrayp),ez(maxrayp)
      dimension tmp(maxrayp),tmpr(maxrayp)
!hpf$ distribute (*,block) :: c6
!hpf$ align (*,:) with c6(*,:) :: c4
!hpf$ align (:) with c6(*,:) :: ex,ey,ez,tmp,tmpr,msk
      common/robhack/zlocate
      common/nxyzsave/nxsave,nysave,nzsave,noresize,nadj0
      common/poiblock1/nspchset,nsckick
      common/gridsz3d/xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz,hxi,hyi,hzi
      common/newpoisson/idirectfieldcalc,idensityfunction,idummy(4)
      common/scparams/rparams(60),cparams(60)
      common/showme/iverbose
      common/envstuff/nenvtrk
      common/autotrk/lautotrk,ntrktype,ntrkorder
c
c     write(6,*)'inside autospch'
cryne-April 8, 2004 new code for integrating envelope equations
      if(nenvtrk.eq.1)then
c       write(6,*)'calling envkick'
        call envkick(tau)
c       write(6,*)'returned from envkick'
      endif
      if(ntrktype.eq.0)return
c     write(6,*)'still inside autospch, performing space charge calc'
cryne-Now continue with the usual space-charge calculation
c
c check that space charge parameters have been set:
c note that the elements of rparams are initialized to -9999999.
c
      if(rparams(1).lt.-1.d6)then
      if(idproc.eq.0)then
      write(6,*)'ERROR: ATTEMPT TO TRACK WITH SPACE CHARGE,'
      write(6,*)'BUT POISSON SOLVER PARAMETERS HAVE NOT BEEN SET.'
      write(6,*)'RERUN WITH INPUT MODIFIED TO INCLUDE POISSON'
      endif
      call myexit
      endif
c
c
c     if(curr0.eq.-99999.)then
c      write(6,*)                                                       &
c    & 'error:entered space charge routine but current has not been set'
c      stop
c     endif
cryne 1/6/03      tau=2.d0*prevlen
      curr=bcurr
!     if(idproc.eq.0)write(12,*)'[autospch]curr,tau=',curr,tau
      if(curr.eq.0.)then
        if(iverbose.eq.2)write(6,*)'returning from sc routine (curr=0)'
        return
      endif
ccc
      straight=.true.
c     if(nt2prev.eq.2 .or. nt2prev.eq.4)straight=.false.

      nx0=nxsave
      ny0=nysave
      nz0=nzsave
c     write(6,*)'nx0,ny0,nz0=',nx0,ny0,nz0
c
      clite=299792458.d0
      fpei=clite*clite*1.d-7
      vz0=beta*clite
      perv=2.d0*curr*fpei/(brho*vz0*vz0*gamma*gamma)
c     write(6,*)'perv=',perv
c     write(6,*)'curr,brho,gamma=',curr,brho,gamma
c     write(6,*)'fpei,beta,vz0=',fpei,beta,vz0
cryne 11/19/02      if(nraysp.ne.maxrayp)then
cryne 11/19/02        write(6,*)'error:autospch works only w/ nrays=maxrayp'
cryne 11/19/02        write(6,*)'nraysp,maxrayp=',nraysp,maxrayp
cryne 11/19/02        stop
cryne 11/19/02      endif
      if(nz0.eq.0)then
        do j=1,nraysp
        c4(1,j)=sl*zblock(1,j)
        c4(2,j)=   zblock(2,j)
        c4(3,j)=sl*zblock(3,j)
        c4(4,j)=   zblock(4,j)
        enddo
!*******************************************************
        do j=1,nraysp
          if(istat(j).eq.0)then
            msk(j)=.true.
          else
            msk(j)=.false.
          endif
        enddo
        if(nraysp.lt.maxrayp)then
          do j=nraysp+1,maxrayp
            msk(j)=.false.
          enddo
        endif
!*******************************************************
        nx=nx0
        ny=ny0
        n1=2*nx
        n2=2*ny
!       if(idproc.eq.0)then
!       write(6,*)'Note: 2D space charge only works in static units'
!       endif
        if(idirectfieldcalc.eq.0)then
!         if(idproc.eq.0)write(6,*)'calling spch2dphi'
          call spch2dphi(c4,ex,ey,msk,nraysp,nrays,nx,ny,n1,n2)
        else
!         if(idproc.eq.0)write(6,*)'calling slic2d'
          call spch2ddirect(c4,ex,ey,msk,nraysp,nrays,nx,ny,n1,n2)
        endif
!4/28/03 put 1/nrays in spch2d
!4/28/03 also put 1/2pi in spch2d, which means need *2pi here
!4/28/03xycon=perv/nrays
        twopi=4.d0*asin(1.d0)
        xycon=twopi*perv
c 4/26/03 case of dynamic units not yet tested, but this should work:
        if(.not.lflagmagu)then
        xycon=xycon*gamma*beta
        endif
        do i=1,nraysp
        ex(i)=ex(i)*xycon
        ey(i)=ey(i)*xycon
        enddo
c
        do j=1,nraysp
        zblock(2,j)=zblock(2,j)+tau*ex(j)
        zblock(4,j)=zblock(4,j)+tau*ey(j)
        enddo
        zlocate=zlocate+tau
        return
      endif
c
      if(nz0.ne.0)then
        if(iverbose.eq.2.and.idproc.eq.0)write(6,*)'3D spchrg'
        do j=1,nraysp
        do i=1,6
        c6(i,j)=zblock(i,j)
        enddo
        enddo
c
c <convert to real units here>
cryne   bbyk=beta*sl !Jan 31, 2003
        bbyk=beta*clite/omegascl
c use ex as a temporary to hold beta of the nth particle:
cryne 3/4/2004        do j=1,nraysp
cryne 3/4/2004        ex(j)=(gamma-c6(6,j))**2-1.-c6(2,j)**2-c6(4,j)**2
cryne 3/4/2004        ex(j)=sqrt(ex(j))/(gamma-c6(6,j))
cryne 3/4/2004        enddo
c gamma (not gamma of the nth particle) takes us to the bunch frame:
c this converts phases to z-positions relative to the fiducial ptcl:
c Rob: need to check signs!!!
c Also, convert x and y to their values at the fixed time:
c       c6(5,1:nraysp)=c6(5,1:nraysp)*ex(1:nraysp)*sl*gamma
c       c6(1,1:nraysp)=(c6(1,1:nraysp)+c6(2,1:nraysp)*c6(5,1:nraysp))*sl
c       c6(3,1:nraysp)=(c6(3,1:nraysp)+c6(4,1:nraysp)*c6(5,1:nraysp))*sl
c same formulae as Impact code:
c 09/04/01 commented out code to treat space charge in bending magnets.
c will fix later. rdr and ctm
c       if(.not.straight)then
c         rho=rhoprev
c         theta0=arclen/rho+3.*asin(1.0d0)
c         write(6,*)'******************  NT2PREV,RHO=',nt2prev,rho
c 0th order approximation for theta:
c         tmp(1:nraysp)=-zblock(5,1:nraysp)*sl*beta/rho
c 1st order approximation for theta:
c         tmp(:)=tmp(:)-sl/rho*sin(tmp(:))*c6(1,:)
c         tmpr(:)=zblock(1,:)*cos(tmp(:))+
c    &    zblock(2,:)*(rho/sl)/(beta*gamma)*sin(tmp(:))-
c    &    zblock(6,:)*(rho/sl)/(beta*gamma)/beta*(1.-cos(tmp(:)))
c         c6(1,:)=(rho/sl+tmpr(:))*sl*sin(theta0+tmp(:))
c         c6(5,:)=(rho/sl+tmpr(:))*sl*cos(theta0+tmp(:))
c         c6(3,1:nraysp)=zblock(3,1:nraysp)*sl
c       else
          if(nadj0.ne.0)then
            toopi=4.d0*asin(1.d0)
            do j=1,nraysp
             c5j=c6(5,j)
             ac5j=abs(c5j)
             if(c5j.gt.0.)c6(5,j)=mod(c5j,toopi)
             if(c5j.lt.0.)c6(5,j)=toopi-mod(ac5j,toopi)
            enddo
c?!       endif
c  place the value which is in(-pi,pi) back in the zblock array:
          do j=1,nraysp
          zblock(5,j)=c6(5,j)-0.5d0*toopi
          enddo
          endif
c convert phase to z (bbyk) and transform to the beam frame (gamma)
c also multiply by scale length (sl) to convert x and y to units of meters
          do j=1,nraysp
            c6(5,j)=c6(5,j)*gamma*bbyk
            c6(1,j)=c6(1,j)*sl
            c6(3,j)=c6(3,j)*sl
          enddo
c       endif
c
cryne Feb 3, 2003:
c this is the "bunch length" in the beam frame when it extends from 0 to 2pi:
c       gblam=gamma*beta*4.*asin(1.0d0)*sl
        gblam=gamma*beta*299792458.d0/bfreq
c       if(idproc.eq.0)write(96,*)'gblam=',gblam
cryne 11/19/02  mask off lost particles
        do j=1,nraysp
          msk(j)=.true.
        enddo
        if(nraysp.lt.maxrayp)then
          do j=nraysp+1,maxrayp
            msk(j)=.false.
          enddo
        endif
c       if(idproc.eq.0)then
c       zmintmp=minval(c6(5,:),msk)
c       zmaxtmp=maxval(c6(5,:),msk)
c       write(6,*)'zmintmp,zmaxtmp=',zmintmp,zmaxtmp
c       endif
c       write(6,*)'calling setbound'
        call setbound(c6,msk,nraysp,gblam,gamma,nx0,ny0,nz0,noresize,   &
     &                nadj0)
        gaminv=1.d0/gamma
        if(idproc.eq.0)then
          write(91,891)arclen,xmin,xmax,hx,nx0
          write(92,891)arclen,ymin,ymax,hy,ny0
cryne 11/27/03     write(93,891)arclen,zmin*gaminv,zmax*gaminv,hz*gaminv,nz0
          write(93,892)arclen,zmin,zmax,hz,nz0
  891     format(4(1pe14.7,1x),i5)
  892     format(1x,'(bunch frame) ',4(1pe14.7,1x),i5)
        endif
c
        nxsave=nx0
        nysave=ny0
        nzsave=nz0
        nx=nx0
        ny=ny0
        nz=nz0
        n1=2*nx
        n2=2*ny
        n3=2*nz
        nadj=nadj0
        n3a=2*nz-nadj*nz
        if(perv.ne.0)then
          if (iverbose.ge.2) write(6,*) 'calling spch3d'
! rob: this will not work until the units are dealt with
          call increment_timer('spch3d',0)
          call spch3d(c6,ex,ey,ez,msk,nraysp,nrays,                     &
     &                nx,ny,nz,n1,n2,n3,n3a,nadj,rparams,cparams)
          call increment_timer('spch3d',1)
          if (iverbose.ge.2) write(6,*) 'returned from spch3d'
        else
cryne Due to code near the top, the following will never be executed:
          write(6,*)'code should not get here'
c         ex=0.
c         ez=0.
c         ey=0.
        endif
        do j=1,nraysp
          ex(j)=ex(j)*gamma
          ey(j)=ey(j)*gamma
        enddo
c       if(idproc.eq.0)then
c         write(96,*)'premult: ex(1),ey(1),ez(1)=',ex(1),ey(1),ez(1)
c       endif
c
        twopi=4.*asin(1.0d0)
c       xycon=0.5*perv*gamma*beta*(bbyk*twopi)
cc      xycon=0.5*perv*gamma*beta*(beta*sl*twopi)
ccc     wrf=twopi*bfreq
cccc    xycon=0.5*perv*gamma*beta*(twopi*beta*clite/wrf)
c (bfreq is in module beamdata)
c********
c 3/21/2003 factor of 1/(4pi) now appears in G, so *4pi is needed here:
      fourpi=8.d0*(asin(1.d0))
c********
c3/21/03xycon=0.5*perv*gamma*beta*(beta*clite/bfreq)
      xycon=fourpi*0.5*perv*gamma*beta*(beta*clite/bfreq)
cryne 11/6/03triedthis,no good:xycon=fourpi*0.5*perv*gamma*beta*(beta*clite/freqscl)
        tcon=beta*xycon*gamma**2
c----------
cryne Jan 28, 2003
        ratio3=omegascl*sl/299792458.d0
c       if(idproc.eq.0)write(6,*)'ratio3,ratio3inv=',ratio3,1.d0/ratio3
c       if(idproc.eq.0)write(96,*)'ratio3,r3inv=',ratio3,1.d0/ratio3
ccc     xycon=xycon*ratio3
        tcon=tcon/ratio3
c----------
c       write(6,*),'xycon,tcon,perv=',xycon,tcon,perv
c       if(idproc.eq.0)then
c       write(96,*)'post:*xytcon=',xycon*ex(1),xycon*ey(1),tcon*ez(1)
c       endif
        gbi=1./(gamma*beta)
cryne 12/24/2002
        if(lflagmagu)then
c       if(idproc.eq.0)write(6,*)'magnetic conversion for sc kick'
c       if(idproc.eq.0)write(96,*)'magnetic conversion for sc kick'
c       if(idproc.eq.0)write(6,*)'gamma*beta,gbi=',gamma*beta,gbi
c       if(idproc.eq.0)write(96,*)                                          &
c    &  'post:*con*gbi=',xycon*ex(1)*gbi,xycon*ey(1)*gbi,tcon*ez(1)*gbi
!       xycon=xycon*gamma*beta
!       tcon=  tcon*gamma*beta
        xycon=xycon*gbi
        tcon=tcon*gbi
!
!       if(idproc.eq.0)write(6,*) 'TCON SET TO ZERO!!!!!!!!!!!!!!!!!!!!'
!       if(idproc.eq.0)write(96,*)'TCON SET TO ZERO!!!!!!!!!!!!!!!!!!!!'
!       tcon=tcon*0.d0
        endif
c
cccc    zblock(2,1:nraysp)=zblock(2,1:nraysp)+tau*xycon*ex(1:nraysp)*gbi
cccc    zblock(4,1:nraysp)=zblock(4,1:nraysp)+tau*xycon*ey(1:nraysp)*gbi
cccc    zblock(6,1:nraysp)=zblock(6,1:nraysp)+tau* tcon*ez(1:nraysp)*gbi
c       if(straight)then
          do j=1,nraysp
          zblock(2,j)=zblock(2,j)+tau*xycon*ex(j)
          zblock(4,j)=zblock(4,j)+tau*xycon*ey(j)
          zblock(6,j)=zblock(6,j)+tau* tcon*ez(j)
          enddo
c       else
c        zblock(2,1:nraysp)=zblock(2,1:nraysp)+
c    &                     tau*xycon*sqrt(ex(1:nraysp)**2+ey(1:nraysp)**2)
c        zblock(4,1:nraysp)=zblock(4,1:nraysp)-tau*xycon*ey(1:nraysp)
c        zblock(6,1:nraysp)=zblock(6,1:nraysp)+tau* tcon*ez(1:nraysp)
c       endif
        zlocate=zlocate+tau
      endif
      if(iverbose.eq.2)write(6,*)'returning from space charge routine'
      return
      end
c
      subroutine confoc(l,kx,ky,kz,h,mh)
      use beamdata
      use lieaparam, only : monoms
      double precision kx,ky,kz,l
      double precision h,mh,gb
      dimension h(monoms),mh(6,6)
      call ident(h,mh)
      gb=gamma*beta
      mh(1,1)=cos(kx*l)
      if(kx.ne.0.d0)then
        mh(1,2)=1./kx/sl*sin(kx*l)
      else
        mh(1,2)=l/sl
      endif
      mh(2,1)=-kx*sl*sin(kx*l)
      mh(2,2)=cos(kx*l)
      mh(3,3)=cos(ky*l)
      if(ky.ne.0.d0)then
        mh(3,4)=1./ky/sl*sin(ky*l)
      else
        mh(3,4)=l/sl
      endif
      mh(4,3)=-ky*sl*sin(ky*l)
      mh(4,4)=cos(ky*l)
      mh(5,5)=cos(kz*l)
      if(kz.ne.0.d0)then
        mh(5,6)=1./kz/sl/gb**2*sin(kz*l)
      else
        mh(5,6)=l/sl
      endif
      mh(6,5)=-kz*sl*gb**2*sin(kz*l)
      mh(6,6)=cos(kz*l)
      return
      end
c
      logical function defcon(line)
cryne 09/21/02 modified to deal with the issue that a symbol may be
c defined using ":=" or "="
c The original version of this routine could not handle "="
c Now it can, and as a result MaryLie/IMPACT can parse both cases.
c (Note, however, that MaryLie/IMPACT treats them identically;
c  MAD parses both but treats them differently)
      character (len=*) line
      character*16 symb(50)
      integer istrtsym(50)
c     write(6,*)'INSIDE DEFCON; line ='
c     write(6,*)line
      defcon=.false.
c find the first nonblank character:
      do n=1,LEN(line)
        nfirst=n
        if(line(n:n).ne.' ')goto 100
      enddo
c blank line:
      return
  100 continue
c is the first character a "!" ?
      if(line(nfirst:nfirst).eq.'!')return
c is the first character a '#' ?
      if(line(nfirst:nfirst).eq.'#')return
c now inspect the contents of the line:
c check for ':' and '=' :
      n1=index(line,':')
      n2=index(line,'=')
c     if((n1.eq.0).or.(n2.eq.0))return
      if(n2.eq.0)return
c If both are present, there can be nothing in between:
      if(n1.ne.0)then
        if(n2.ne.n1+1)return
      endif
c Set ncheck equal to the starting point of "=" or ":="
      if(n1.ne.0)then
        ncheck=n1
      else
        ncheck=n2
      endif
c
c there must be at least one symbol, and ":="  or "=" must occur after it:
      call getsymb(line,LEN(line),symb,istrtsym,nsymb)
      if(nsymb.eq.0)return
      m1=istrtsym(1)
c     m11=index(line,symb(1))
c     write(6,*)'(defcon) m1,m11=',m1,m11
      m2=len_trim(symb(1))
c n0=location of last character of first string:
      n0=m1+m2-1
c     if(n0.lt.n1)defcon=.true.
      if(ncheck.eq.n0+1)defcon=.true.
      if(ncheck.gt.n0+1)then
        do mmm=n0+1,ncheck-1
          if(line(mmm:mmm).ne.' ')return
        enddo
        defcon=.true.
      endif
c     write(6,*)'**FOUND A DEFCON**'
c     if(defcon)then
c       if(ncheck.eq.n2)then
c         write(6,*)'*****FOUND A DEFCON defined with an equal sign in'
c       else
c         write(6,*)'*****FOUND A DEFCON defined with a :=  in'
c       endif
c       write(6,*)line
c     endif
      return
      end
c
      logical function defline(line)
      character (len=*) line
      character*16 symb(50)
      integer istrtsym(50)
c     write(6,*)'INSIDE DEFLINE; line ='
c     write(6,*)line
      defline=.false.
c find the first nonblank character:
      do n=1,LEN(line)
        nfirst=n
        if(line(n:n).ne.' ')goto 100
      enddo
c blank line:
      return
  100 continue
c is the first character a '!' ?
      if(line(nfirst:nfirst).eq.'!')return
c is the first character a '#' ?
      if(line(nfirst:nfirst).eq.'#')return
c now check the contents of the line:
c '=' must be present
      n1=index(line,'=')
      if(n1.eq.0)return
c     write(6,*)'n1 (location of equal sign) is',n1
c and 'line' must be present
      call getsymb(line,LEN(line),symb,istrtsym,nsymb)
c     write(6,*)'(defline) nsymb=',nsymb
c     if(nsymb.gt.0)then
c       do ijk=1,nsymb
c         write(6,*)'i,symb(i)=  ',ijk,symb(ijk)
c       enddo
c     endif
      if(nsymb.lt.2)return
      if((trim(symb(2)).ne.'line').and.(trim(symb(2)).ne.'LINE'))return
c and there can only be blanks (or nothing) between 'line' and '='
      m1=istrtsym(2)
c     m11=index(line,trim(symb(2)))
c     write(6,*)'(defline)m1,m11=',m1,m11
c     write(6,*)'m1 (start of second symbol)=',m1
c n0=location of last character of second string:
      n0=m1+3
c     write(6,*)'n0 (end of second symbol)=',n0
      if(n1.eq.n0+1)then
        defline=.true.
c       write(6,*)'*****(1)FOUND A DEFLINE*****,line='
c       write(6,*)line
        return
      endif
      do n=n0+1,n1-1
        if(line(n:n).ne.' ')return
      enddo
      defline=.true.
c     write(6,*)'*****(2)FOUND A DEFLINE*****'
c     write(6,*)line
      return
      end
c
      logical function defelem(line)
      include 'codes.inc'
      character (len=*) line
      character*16 symb(50)
      integer istrtsym(50)
      character*16 str
      dimension nstrmax(9)
      nstrmax(1)=75
      nstrmax(2)=10
      nstrmax(3)=9
      nstrmax(4)=24
      nstrmax(5)=9
      nstrmax(6)=9
      nstrmax(7)=50
      nstrmax(8)=39
      nstrmax(9)=37
c     write(6,*)'INSIDE DEFELEM; line ='
c     write(6,*)line
      defelem=.false.
c find the first nonblank character:
      do n=1,LEN(line)
        nfirst=n
        if(line(n:n).ne.' ')goto 100
      enddo
c blank line:
      return
  100 continue
c is the first character a '!' ?
      if(line(nfirst:nfirst).eq.'!')return
c is the first character a '#' ?
      if(line(nfirst:nfirst).eq.'#')return
c now check the contents of the line
c first obtain the symbols on the line:
      call getsymb(line,LEN(line),symb,istrtsym,nsymb)
c it must contain at least TWO symbols
      if(nsymb.lt.2)return
c     write(6,*)'symb(1),symb(2)  =',symb(1),'----',symb(2),'----'
c check for 'beam' or 'units'
      if((trim(symb(1)).eq.'beam').or.(trim(symb(1)).eq.'units').or.      &
     &   (trim(symb(1)).eq.'BEAM').or.(trim(symb(1)).eq.'UNITS'))then
c       write(6,*)'FOUND A BEAM OR UNITS DEFELEM'
        defelem=.true.
        return
      endif
c
c a ':' must occur between them
      m0=index(line,':')
c     write(6,*)'m0=',m0
      if(m0.eq.0)return
      ms1=istrtsym(1)
c     ms11=index(line,trim(symb(1)))
c     write(6,*)'(defelem)  ms1,ms11=',ms1,ms11
      ms2=istrtsym(2)
c     ms22=index(line,trim(symb(2)))
c     write(6,*)'(defelem)  ms2,ms22=',ms2,ms22
c     write(6,*)'m0,ms1,ms2 =',m0,ms1,ms2
      if((m0.lt.ms1).or.(m0.gt.ms2))return
c '=' must not occur between them
      k0=index(line,'=')
c     write(6,*)'k0=',k0
      if(k0.ne.0)then
        if((k0.gt.ms1).and.(k0.lt.ms2))return
      endif
c the second symbol must match a type code name:
      str=symb(2)(1:16)
c     write(6,*)'str=',str
      do i=1,9
      do j=1,nstrmax(i)
      if(str(1:16).eq.ltc(i,j))then
        defelem=.true.
c       write(6,*)'***FOUND A DEFELEM***; i,j,ltc(i,j)=',i,j,ltc(i,j)
      endif
      enddo
      enddo
      return
      end
c
      subroutine rfcavmad(zlen,volt,phlagrad,nharm,h,mh)
c MAD RF CAVITY routine
c note that the length,zlen, is not used.
c
c trevi is the inverse revolution frequency.
c how/where is this info obtained???
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision mh
      double precision omega
      dimension h(monoms),mh(6,6)
c--------
      write(6,*)'Inside routine for MAD RF CAVITY; do not know T_rev.'
      write(6,*)'This needs to be fixed. Setting 1./T_rev = 0'
      trevi=0.d0
c--------
      omega=4.d0*asin(1.d0)*nharm*trevi
      call ident(h,mh)
      mh(6,5)=-omega*volt/brho*cos(phlagrad)
      return
      end
c
      subroutine elsepmad(zlen,efield,h,mh)
c MAD ELECTROSTATIC SEPARATOR routine
      use lieaparam, only : monoms
      use beamdata
      include 'impli.inc'
      double precision mh
      double precision omega
      double precision xk
      dimension h(monoms),mh(6,6)
      xk=efield/brho
      co=cosh(xk*zlen)
      si=sinh(xk*zlen)
      write(6,*)'in elsepmad w/ zlen,efield,brho=',zlen,efield
      write(6,*)'and brho,beta=',brho,beta
      write(6,*)'and xk,co,si=',xk,co,si
      call ident(h,mh)
      mh(1,2)=zlen
      mh(3,3)=co-xk*zlen*si/beta**2
      mh(3,4)=si/xk
      mh(3,6)=(co-1.d0)/xk - zlen*si/beta**2
      mh(4,3)=xk*(si-xk*zlen*co/beta**2)
      mh(4,4)=co
      mh(4,6)=si-xk*zlen*co/beta**2
      mh(5,3)=-(si-xk*zlen*co/beta**2)
      mh(5,4)=-(co-1.d0)/xk
      mh(5,6)=-si/xk + zlen*co/beta**2
      return
      end
c
      subroutine hmonitor
      include 'impli.inc'
c     write(6,*)'MAD hmonitor not implemented'
      return
      end
c
      subroutine vmonitor
      include 'impli.inc'
c     write(6,*)'MAD vmonitor not implemented'
      return
      end
c
      subroutine monitor
      include 'impli.inc'
c     write(6,*)'MAD monitor not implemented'
      return
      end
c
      subroutine instrument
      include 'impli.inc'
c     write(6,*)'MAD instrument not implemented'
      return
      end
c
c***********************************************************************
c
      subroutine sliceml
c the routine augments the nrp ("number of required parameters") values
c of thick elements associated with original MaryLie type codes.
c By putting this in the menu first, users with original MaryLie
c input files can add a parameter to the parameter lists of the
c thick elements in order to do autoslicing.
      use parallel
      include 'codes.inc'
      nrp(1,1)=nrp(1,1)+1   !drft
      nrp(1,2)=nrp(1,2)+1   !nbnd
      nrp(1,3)=nrp(1,3)+1   !pbnd
      nrp(1,4)=nrp(1,4)+1   !gbnd
      nrp(1,6)=nrp(1,6)+1   !gbdy
      nrp(1,8)=nrp(1,8)+1   !cfbd
      nrp(1,9)=nrp(1,9)+1   !quad
      nrp(1,10)=nrp(1,10)+1   !sext
      nrp(1,11)=nrp(1,11)+1   !octm
      nrp(1,12)=nrp(1,12)+1   !octe
      nrp(1,18)=nrp(1,18)+1   !cfqd
      nrp(1,20)=nrp(1,20)+1   !sol
      nrp(1,24)=nrp(1,24)+1   !recm
      nrpold(1,1)=nrpold(1,1)+1   !drft
      nrpold(1,2)=nrpold(1,2)+1   !nbnd
      nrpold(1,3)=nrpold(1,3)+1   !pbnd
      nrpold(1,4)=nrpold(1,4)+1   !gbnd
      nrpold(1,6)=nrpold(1,6)+1   !gbdy
      nrpold(1,8)=nrpold(1,8)+1   !cfbd
      nrpold(1,9)=nrpold(1,9)+1   !quad
      nrpold(1,10)=nrpold(1,10)+1   !sext
      nrpold(1,11)=nrpold(1,11)+1   !octm
      nrpold(1,12)=nrpold(1,12)+1   !octe
      nrpold(1,18)=nrpold(1,18)+1   !cfqd
      nrpold(1,20)=nrpold(1,20)+1   !sol
      nrpold(1,24)=nrpold(1,24)+1   !recm
      if(idproc.eq.0)then
!     write(6,*)'Enabling slicing of original MaryLie thick elements.'
!     write(6,*)'For each thick element, the MaryLie parser will'
!     write(6,*)'expect 1 more parameter than is shown in the MaryLie'
!     write(6,*)'manual. The extra (last) parameter = the # of slices.'
      write(12,*)'Enabling slicing of original MaryLie thick elements.'
      write(12,*)'For each thick element, the MaryLie parser will'
      write(12,*)'expect 1 more parameter than is shown in the MaryLie'
      write(12,*)'manual. The extra (last) parameter = the # of slices.'
      endif
      return
      end
c
!***********************************************************************
!     subroutine myflush(nfile)
!     integer nfile
!     double precision x
!     close(nfile)
!     open(nfile,position='append')
c     open(nfile)
c 100 read(nfile,*,end=123,err=999)x
c     goto 100
c 123 continue
c     return
c 999 continue
c     if(idproc.eq.0)write(6,*)'trouble in routine myflush'
c     call myexit
!     return
!     end
!***********************************************************************
      subroutine num2string(num,a,ndigits)
c converts an integer "num" with at most "ndigits" digits
c into a character string "a"
      integer num,ndigits
      character a*(*)
      integer n,m,k
      m=num
      do n=1,ndigits
      k=m/10**(ndigits-n)
      if(k.eq.0)a(n:n)='0'
      if(k.eq.1)a(n:n)='1'
      if(k.eq.2)a(n:n)='2'
      if(k.eq.3)a(n:n)='3'
      if(k.eq.4)a(n:n)='4'
      if(k.eq.5)a(n:n)='5'
      if(k.eq.6)a(n:n)='6'
      if(k.eq.7)a(n:n)='7'
      if(k.eq.8)a(n:n)='8'
      if(k.eq.9)a(n:n)='9'
      m=m-k*10**(ndigits-n)
      enddo
      return
      end
c***********************************************************************
      subroutine ibcast(ival)
      use parallel
      integer ival
      if(idproc.eq.0)then
        do l=1,nvp-1
        call MPI_SEND(ival,1,mntgr,l,99,lworld,ierr)
!       write(6,*)'processor ',idproc,' sending ival=',ival,' to PE',l
        enddo
      else
        call MPI_RECV(ival,1,mntgr,0,99,lworld,mpistat,ierr)
c?      call MPI_GET_COUNT(mpistat,mreal,nraysp,ierr)
!       write(6,*)'processor ',idproc,' received ival=',ival
      endif
      return
      end
c***********************************************************************
      subroutine rbcast(rval)
      use parallel
      real*8 rval
      if(idproc.eq.0)then
        do l=1,nvp-1
        call MPI_SEND(rval,1,mreal,l,99,lworld,ierr)
!       write(6,*)'processor ',idproc,' sending rval=',rval,' to PE',l
        enddo
      else
        call MPI_RECV(rval,1,mreal,0,99,lworld,mpistat,ierr)
c?      call MPI_GET_COUNT(mpistat,mreal,nraysp,ierr)
!       write(6,*)'processor ',idproc,' received rval=',rval
      endif
      return
      end
***********************************************************************
c
      subroutine myexit
!     use parallel
      use acceldata
      use spchdata
      use ml_timer
      use rays
      common/xtrajunk/jticks1,iprinttimers
c
c  Written by Rob Ryne ca 1984
c
c Exit the program (may be replaced by STOP, if exit is not available).
c This is the only place in MaryLie where exit is called.
c
      call system_clock(count=jticks2)
      elapsed=(jticks2-jticks1)/hertz
      write(6,*)'ELAPSED TIME=',elapsed
      if(iprinttimers.eq.1)call end_ml_timers
c     write(6,*)'deallocating...'
      call del_acceldata
      if(allocated(zblock))call del_particledata
      if(allocated(grnxtr))call del_spchdata
      call end_parallel
      stop
      end
c
c
      subroutine getordtyp(cstring,ntaysym,norder)
      use parallel
      include 'impli.inc'
      character*16 cstring
      ntay=index(cstring,'tay')
      nsym=index(cstring,'sym')
      if(ntay.ne.0)ntaysym=1
      if(nsym.ne.0)ntaysym=2
c inefficient, but works. fix later. rdr dec 9 2002 & march 27 2004
      n0=index(cstring,'0')
      if(n0.ne.0)then
        norder=0
        return
      endif
      n1=index(cstring,'1')
      if(n1.ne.0)then
        norder=1
        return
      endif
      n2=index(cstring,'2')
      if(n2.ne.0)then
        norder=2
        return
      endif
      n3=index(cstring,'3')
      if(n3.ne.0)then
        norder=3
        return
      endif
      n4=index(cstring,'4')
      if(n4.ne.0)then
        norder=4
        return
      endif
      n5=index(cstring,'5')
      if(n5.ne.0)then
        norder=5
        return
      endif
      n6=index(cstring,'6')
      if(n6.ne.0)then
        norder=6
        return
      endif
      n7=index(cstring,'7')
      if(n7.ne.0)then
        norder=7
        return
      endif
      n8=index(cstring,'8')
      if(n8.ne.0)then
        norder=8
        return
      endif
      n9=index(cstring,'9')
      if(n9.ne.0)then
        norder=9
        return
      endif
      if(idproc.eq.0)then
       write(6,*)'(GETORDTYP)Unable to determine tracking type & order'
      endif
      call myexit
      end
c
