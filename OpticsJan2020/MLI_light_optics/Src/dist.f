*******************************************************************************
* header:    DIST                                                             *
*  A collection of programs for generating, modifying, and                    *
*  evaluating particle distributions                                          *
*******************************************************************************
      subroutine bgen(p,fa,fm)
c
c subroutine corresponding to the type code bgen
c written by Alex Dragt 6/14/91
c modified 10/14/98 AJD
c modified by RDR on 7/28/2002. The quantity "nray" which is a dummy
c argument to several routines could be eliminated, in which case
c it should be replaced by nraysp in the body of the routines.
c I have left it in place in order to minimize changes to these routines.
c Note that, except for this routine (bgen), the variable "nrays"
c should not appear anywhere, since it is the global # of particles,
c and the remaining routines all use the local number, nraysp, which is
c passed in the argument lists and corresponds to "nray" in the routines.
c 
c
      use lieaparam, only : monoms
      use rays
      include 'impli.inc'
      include 'buffer.inc'
      include 'files.inc'
      include 'parset.inc'
c
c calling arrays
      dimension p(*)
      dimension fa(*), fm(6,6)
c 
c working arrays
c use buf1 and buf2: put numerical moments in buf1 and analytic
c moments in buf2
c use buf3 as working space to store desired eigen moments
c
c set up control indices
c
      job=nint(p(1))
      iopt=nint(p(2))
      nrays=nint(p(3))
      iseed=nint(p(4))
      isend=nint(p(5))
      ipset=nint(p(6))
c
c read contents of pset ipset
c
      if(ipset.ne.0)then
        x2mom=pst(1,ipset) 
        y2mom=pst(2,ipset)
        t2mom=pst(3,ipset)
        if (x2mom .lt. 0.) then
c get eigenmoments from current map
        x2mom=fa(7)
        y2mom=fa(18)
        t2mom=fa(25)
        endif
        sigmax=pst(4,ipset)
        nof1=nint(pst(5,ipset))
        nof2=nint(pst(6,ipset))
      else
        x2mom=p(7)
        y2mom=p(8)
        t2mom=p(9)
        if (x2mom .lt. 0.) then
c get eigenmoments from current map
        x2mom=fa(7)
        y2mom=fa(18)
        t2mom=fa(25)
        endif
        sigmax=p(10)
        nof1=nint(p(11))
        nof2=nint(p(12))
      endif
c     write(6,*)'done setting bgen parameters'
c     write(6,*)'job,iopt,nrays=',job,iopt,nrays
c     write(6,*)'iseed,isend,ipset=',iseed,isend,ipset
c     write(6,*)'x2mom,y2mom,t2mom=',x2mom,y2mom,t2mom
c     write(6,*)'sigmax=',sigmax
c     write(6,*)'nof1,nof2=',nof1,nof2
c
c put eigen moments in buf3
c
      call clear(buf3a,buf3m)
      buf3a(7)=x2mom
      buf3a(13)=x2mom
      buf3a(18)=y2mom
      buf3a(22)=y2mom
      buf3a(25)=t2mom
      buf3a(27)=t2mom
c
c select job
c
c compute moments of particle distribution, and write them on
c an external file, the terminal and/or a drop file
      if (job .eq. 0) then
c
c	  if (iopt .ne. 6) then
c
c compute moments of particle distribution
      call cmom(buf1a,buf1m)
c write moments on file nof1
      mpt=mpo
      mpo=nof1
      if(mpo .gt. 0) call mapout(0,buf1a,buf1m)
      mpo=mpt
c
c	  else
c compute moments of particle distribution including <5> and <6>
c      call cmom5(buf1a)
c write moments on file nof1
c      mpt=mpo
c      mpo=nof1
c      call mapout5(0,buf1a,buf1m)
c      mpo=mpt
c	  endif
c
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) '  '
      write(jof,*)
     & ' numerically computed values of selected moments'
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) '  '
      write(jodf,*)
     & ' numerically computed values of selected moments'
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
c
c for other values of job generate various distributions 
c and/or moments
c
cryne 7/28/2002
c before particles can be generated, each processor needs to know
c how many it owns:
        nraysp=(nrays-1)/nvp + 1
 1221   continue
        if(nraysp*(nvp-1).ge.nrays)then
          nraysp=nraysp-1
          if(nraysp.eq.0)then
          write(6,*)'trouble: nraysp=0. something wrong. halting.'
          endif
          goto 1221
        endif
c
        if(idproc.eq.nvp-1)nraysp = nrays - nraysp*(nvp-1)
c has the zblock array been allocated?
        if(.not.allocated(zblock))then
        if(idproc.eq.0)write(6,*)'setting maxray=',nrays,' ;allocating'
        maxray=nrays
        call new_particledata
        endif
c also, need a different seed on every processor:
        seedfac=1.+float(idproc)/float(nvp)
        iseed=iseed*seedfac
c       write(6,*)'idproc,iseed=',idproc,iseed
c warm up the F90 random number generator because it is needed
c to generate big vectors of random numbers:
        call f90ranset(iseed)
c
c uniformly filled ellipse or ellipsoids
      if (job .eq. 1) then
c set up scaling factors
      sx=sqrt(4.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call re2d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmre2d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call re2d(nraysp,iseed,sx,sy,st)
      call cmre2d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call re2d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call re2d(nraysp,iseed,sx,sy,st)
      call cmre2d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 2) then
c set up scaling factors
      sx=sqrt(6.d0*x2mom)
      sy=sqrt(6.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call re4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmre4d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call re4d(nraysp,iseed,sx,sy,st)
      call cmre4d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call re4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call re4d(nraysp,iseed,sx,sy,st)
      call cmre4d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 3) then
c set up scaling factors
      sx=sqrt(8.d0*x2mom)
      sy=sqrt(8.d0*y2mom)
      st=sqrt(8.d0*t2mom)
      if(iopt .eq. 1) call re6d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmre6d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call re6d(nraysp,iseed,sx,sy,st)
      call cmre6d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call re6d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call re6d(nraysp,iseed,sx,sy,st)
      call cmre6d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
c gaussian distributions
      if (job .eq. 4) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call rg2d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) call cmrg2d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rg2d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg2d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rg2d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 123) then
      call rg2d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg2d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 5) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=sqrt(y2mom)
      st=0.d0
      if(iopt .eq. 1) call rg4d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) call cmrg4d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rg4d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg4d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rg4d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 123) then
      call rg4d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg4d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 6) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=sqrt(y2mom)
      st=sqrt(t2mom)
      if(iopt .eq. 1) call rg6d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) call cmrg6d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rg6d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg6d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rg6d(nraysp,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 123) then
      call rg6d(nraysp,iseed,sigmax,sx,sy,st)
      call cmrg6d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
c systematic uniform tori
      if (job .eq. 7) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call st2d(nraysp,sx,sy,st)
      if(iopt .eq. 2) call cmst2d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call st2d(nraysp,sx,sy,st)
      call cmst2d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call st2d(nraysp,sx,sy,st)
      if(iopt .eq. 123) then
      call st2d(nraysp,sx,sy,st)
      call cmst2d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 8) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call st4d(nraysp,sx,sy,st)
      if(iopt .eq. 2) call cmst4d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call st4d(nraysp,sx,sy,st)
      call cmst4d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call st4d(nraysp,sx,sy,st)
      if(iopt .eq. 123) then
      call st4d(nraysp,sx,sy,st)
      call cmst4d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 9) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=sqrt(2.d0*t2mom)
      if(iopt .eq. 1) call st6d(nraysp,sx,sy,st)
      if(iopt .eq. 2) call cmst6d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call st6d(nraysp,sx,sy,st)
      call cmst6d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call st6d(nraysp,sx,sy,st)
      if(iopt .eq. 123) then
      call st6d(nraysp,sx,sy,st)
      call cmst6d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
c random uniform tori
      if (job .eq. 10) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call rt2d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmst2d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rt2d(nraysp,iseed,sx,sy,st)
      call cmst2d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rt2d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call rt2d(nraysp,iseed,sx,sy,st)
      call cmst2d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 11) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call rt4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmst4d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rt4d(nraysp,iseed,sx,sy,st)
      call cmst4d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rt4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call rt4d(nraysp,iseed,sx,sy,st)
      call cmst4d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
      if (job .eq. 12) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=sqrt(2.d0*t2mom)
      if(iopt .eq. 1) call rt6d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmst6d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call rt6d(nraysp,iseed,sx,sy,st)
      call cmst6d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call rt6d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call rt6d(nraysp,iseed,sx,sy,st)
      call cmst6d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
c KV distribution in 4-D phase space
      if (job .eq. 13) then
c set up scaling factors
      sx=sqrt(4.d0*x2mom)
      sy=sqrt(4.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call kv4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 2) call cmkv4d(buf2a,buf2m,sx,sy,st)
      if(iopt .eq. 12) then
      call kv4d(nraysp,iseed,sx,sy,st)
      call cmkv4d(buf2a,buf2m,sx,sy,st)
      endif
      if(iopt .eq. 13) call kv4d(nraysp,iseed,sx,sy,st)
      if(iopt .eq. 123) then
      call kv4d(nraysp,iseed,sx,sy,st)
      call cmkv4d(buf2a,buf2m,sx,sy,st)
      endif
      endif
c
c write to files and write selected numerical/analytic moments 
c at terminal and/or drop file
c
      if (job .ne. 0) then
c
      if (iopt .eq. 2 .or. iopt .eq. 12) then
c write analytic moments to an external file and to 
c terminal and/or drop file 
      mpt=mpo
      mpo=nof2
      if(mpo .gt. 0) call mapout (0,buf2a,buf2m)
      mpo=mpt
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) '  '
      write(jof,*)
     & ' analytically computed values of selected moments'
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf2a(7),buf2a(8),buf2a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf2a(18),buf2a(19),buf2a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf2a(25),buf2a(26),buf2a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) '  '
      write(jodf,*)
     & ' analytically computed values of selected moments'
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf2a(7),buf2a(8),buf2a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf2a(18),buf2a(19),buf2a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf2a(25),buf2a(26),buf2a(27)
      endif
      endif
c
c compute numerical moments and write numerical moments to 
c an external file and to terminal and/or drop file 
      if (iopt .eq. 13) then
      call cmom(buf1a,buf1m)
      mpt=mpo
      mpo=nof1
      if(mpo .gt. 0) call mapout(0,buf1a,buf1m)
      mpo=mpt
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) '  '
      write(jof,*)
     & ' numerically computed values of selected moments'
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) '  '
      write(jodf,*)
     & ' numerically computed values of selected moments'
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
c
c compute numerical moments and write numerical and analytic
c moments to external files and to terminal and/or drop file 
      if (iopt .eq. 123) then
      call cmom(buf1a,buf1m)
      mpt=mpo
      mpo=nof1
      if(mpo .gt. 0) call mapout(0,buf1a,buf1m)
      mpo=nof2
      if(mpo .gt. 0) call mapout(0,buf2a,buf2m)
      mpo=mpt
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) '  '
      write(jof,*)
     & ' numerically computed values of selected moments'
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      write(jof,*) '  '
      write(jof,*)
     & ' analytically computed values of selected moments'
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf2a(7),buf2a(8),buf2a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf2a(18),buf2a(19),buf2a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf2a(25),buf2a(26),buf2a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) '  '
      write(jodf,*)
     & ' numerically computed values of selected moments'
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      write(jodf,*) '  '
      write(jodf,*)
     & ' analytically computed values of selected moments'
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf2a(7),buf2a(8),buf2a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf2a(18),buf2a(19),buf2a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf2a(25),buf2a(26),buf2a(27)
      endif
      endif
c
      endif
c
      return
      end
c
********************************************************************************
c
      subroutine cmre2d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 2D random
c  uniform distribution that fills a 2D ellipse in phase space
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms), fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/4.d0
      fa(13)=sx2/4.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/24.d0)
      fa(140)=sx2*sx2*(1.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmre4d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 4D random
c  uniform distribution that fills a 4D ellipsoid in phase space
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/6.d0
      fa(13)=sx2/6.d0
      fa(18)=sy2/6.d0
      fa(22)=sy2/6.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/16.d0)
      fa(90)=sx2*sx2*(1.d0/48.d0)
      fa(95)=sx2*sy2*(1.d0/48.d0)
      fa(99)=sx2*sy2*(1.d0/48.d0)
      fa(140)=sx2*sx2*(1.d0/16.d0)
      fa(145)=sx2*sy2*(1.d0/48.d0)
      fa(149)=sx2*sy2*(1.d0/48.d0)
      fa(175)=sy2*sy2*(1.d0/16.d0)
      fa(179)=sy2*sy2*(1.d0/48.d0)
      fa(195)=sy2*sy2*(1.d0/16.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmre6d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 6D random
c  uniform distribution that fills a 6D ellipsoid in phase space
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/8.d0
      fa(13)=sx2/8.d0
      fa(18)=sy2/8.d0
      fa(22)=sy2/8.d0
      fa(25)=st2/8.d0
      fa(27)=st2/8.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
      fm(5,5)=fa(25)
      fm(6,6)=fa(27)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/80.d0)
      fa(90)=sx2*sx2*(1.d0/80.d0)
      fa(95)=sx2*sy2*(1.d0/80.d0)
      fa(99)=sx2*sy2*(1.d0/80.d0)
      fa(102)=sx2*st2*(1.d0/80.d0)
      fa(104)=sx2*st2*(1.d0/80.d0)
      fa(140)=sx2*sx2*(3.d0/80.d0)
      fa(145)=sx2*sy2*(1.d0/80.d0)
      fa(149)=sx2*sy2*(1.d0/80.d0)
      fa(152)=sx2*st2*(1.d0/80.d0)
      fa(154)=sx2*st2*(1.d0/80.d0)
      fa(175)=sy2*sy2*(3.d0/80.d0)
      fa(179)=sy2*sy2*(1.d0/80.d0)
      fa(182)=sy2*st2*(1.d0/80.d0)
      fa(184)=sy2*st2*(1.d0/80.d0)
      fa(195)=sy2*sy2*(3.d0/80.d0)
      fa(198)=sy2*st2*(1.d0/80.d0)
      fa(200)=sy2*st2*(1.d0/80.d0)
      fa(205)=st2*st2*(3.d0/80.d0)
      fa(207)=st2*st2*(1.d0/80.d0)
      fa(209)=st2*st2*(3.d0/80.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg2d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 2D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(140)=sx2*sx2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg4d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 4D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
      fa(18)=sy2
      fa(22)=sy2
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(95)=sx2*sy2
      fa(99)=sx2*sy2
      fa(140)=sx2*sx2*(3.d0)
      fa(145)=sx2*sy2
      fa(149)=sx2*sy2
      fa(175)=sy2*sy2*(3.d0)
      fa(179)=sy2*sy2
      fa(195)=sy2*sy2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg6d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 6D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
      fa(18)=sy2
      fa(22)=sy2
      fa(25)=st2
      fa(27)=st2
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
      fm(5,5)=fa(25)
      fm(6,6)=fa(27)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(95)=sx2*sy2
      fa(99)=sx2*sy2
      fa(102)=sx2*st2
      fa(104)=sx2*st2
      fa(140)=sx2*sx2*(3.d0)
      fa(145)=sx2*sy2
      fa(149)=sx2*sy2
      fa(152)=sx2*st2
      fa(154)=sx2*st2
      fa(175)=sy2*sy2*(3.d0)
      fa(179)=sy2*sy2
      fa(182)=sy2*st2
      fa(184)=sy2*st2
      fa(195)=sy2*sy2*(3.d0)
      fa(198)=sy2*st2
      fa(200)=sy2*st2
      fa(205)=st2*st2*(3.d0)
      fa(207)=st2*st2
      fa(209)=st2*st2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst2d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 2D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst4d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 4D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
      fa(18)=sy2/2.d0
      fa(22)=sy2/2.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(95)=sx2*sy2*(1.d0/4.d0)
      fa(99)=sx2*sy2*(1.d0/4.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/4.d0)
      fa(149)=sx2*sy2*(1.d0/4.d0)
      fa(175)=sy2*sy2*(3.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/8.d0)
      fa(195)=sy2*sy2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst6d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 6D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
      fa(18)=sy2/2.d0
      fa(22)=sy2/2.d0
      fa(25)=st2/2.d0
      fa(27)=st2/2.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
      fm(5,5)=fa(25)
      fm(6,6)=fa(27)
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(95)=sx2*sy2*(1.d0/4.d0)
      fa(99)=sx2*sy2*(1.d0/4.d0)
      fa(102)=sx2*st2*(1.d0/4.d0)
      fa(104)=sx2*st2*(1.d0/4.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/4.d0)
      fa(149)=sx2*sy2*(1.d0/4.d0)
      fa(152)=sx2*st2*(1.d0/4.d0)
      fa(154)=sx2*st2*(1.d0/4.d0)
      fa(175)=sy2*sy2*(3.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/8.d0)
      fa(182)=sy2*st2*(1.d0/4.d0)
      fa(184)=sy2*st2*(1.d0/4.d0)
      fa(195)=sy2*sy2*(3.d0/8.d0)
      fa(198)=sy2*st2*(1.d0/4.d0)
      fa(200)=sy2*st2*(1.d0/4.d0)
      fa(205)=st2*st2*(3.d0/8.d0)
      fa(207)=st2*st2*(1.d0/8.d0)
      fa(209)=st2*st2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmkv4d(fa,fm,sx,sy,st)
c  subroutine to compute analytic moments of a 4-variable
c  KV distribution
c  written by Alex Dragt 7/15/91
c
      use lieaparam, only : monoms
      include 'impli.inc'
c
c  calling arrays
      dimension fa(monoms),fm(6,6)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/4.d0
      fa(13)=sx2/4.d0
      fa(18)=sy2/4.d0
      fa(22)=sy2/4.d0
c  matrix of moments
      call mclear(fm)
      fm(1,1)=fa(7)
      fm(2,2)=fa(13)
      fm(3,3)=fa(18)
      fm(4,4)=fa(22)
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/24.d0)
      fa(95)=sx2*sy2*(1.d0/24.d0)
      fa(99)=sx2*sy2*(1.d0/24.d0)
      fa(140)=sx2*sx2*(1.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/24.d0)
      fa(149)=sx2*sy2*(1.d0/24.d0)
      fa(175)=sy2*sy2*(1.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/24.d0)
      fa(195)=sy2*sy2*(1.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine re2d(nray,iseed,sx,sy,st)
c
c  subroutine to generate a 2D random uniform distribution that fills 
c  a 2D ellipse in phase space
c  written by Alex Dragt 7/6/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
c 
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray 
    5 sumsq=0.d0
      do 10 i=1,2
      call myrand(ktype,iseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=0.d0
      zblock(4,j)=0.d0
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by re2d on PE# ',idproc
c
      return
      end
c
********************************************************************************
c
      subroutine re4d(nray,iseed,sx,sy,st)
c
c  subroutine to generate a 4D random uniform distribution that fills 
c  a 4D ellipsoid in phase space
c  written by Alex Dragt 7/6/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray 
    5 sumsq=0.d0
      do 10 i=1,4
      call myrand(ktype,iseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by re4d on PE# ',idproc
c
      return
      end
c
********************************************************************************
c
      subroutine re6d(nray,iseed,sx,sy,st)
c
c  subroutine to generate a 6D random uniform distribution that fills 
c  a 6D ellipsoid in phase space
c  written by Alex Dragt 7/6/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
      integer, parameter :: ichunk=10000
c
c  working arrays
      dimension z(6),za(6*ichunk)
c
c  procedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      write(6,*) 'nray=',nray
      write(6,*) 'maxrayp=',maxrayp
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
c procedure to generate a small number (~few thousand) of rays:
cryne 1/11/2005      if(nray.gt.5000)goto 23
      if(nray.gt.500000)goto 23
      do 20 j=1,nray     
c generate a ray 
    5 sumsq=0.d0
      do 10 i=1,6
      call myrand(ktype,iseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=st*z(5)
      zblock(6,j)=st*z(6)
   20 continue
      goto 999
c
   23 continue
c procedure to generate a large number of rays:
      j=0
   24 continue
      call random_number(za)
      do i=1,6*ichunk-5,6
        sumsq=(2.d0*za(i)  -1.d0)**2+(2.d0*za(i+1)-1.d0)**2               &
     &       +(2.d0*za(i+2)-1.d0)**2+(2.d0*za(i+3)-1.d0)**2               &
     &       +(2.d0*za(i+4)-1.d0)**2+(2.d0*za(i+5)-1.d0)**2
        if(sumsq .gt. 1.d0)cycle
        j=j+1
        if(j.gt.nray)exit
c scale and store ray in zblock
        zblock(1,j)=sx*za(i)
        zblock(2,j)=sx*za(i+1)
        zblock(3,j)=sy*za(i+2)
        zblock(4,j)=sy*za(i+3)
        zblock(5,j)=st*za(i+4)
        zblock(6,j)=st*za(i+5)
      enddo
      if(j.lt.nray)goto 24
c
  999 continue
      write(6,*) nray, ' rays generated by re6d on PE# ',idproc
c
      return
      end
c
********************************************************************************
c
      subroutine rg2d(nray,iseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 2D gaussian distribution 
c  written by Alex Dragt 7/10/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      sigm2=sigmax**2
      do 20 j=1,nray     
c generate a ray 
   25 call normdv(z,1,ktype,iseed)
      ssq = z(1)**2 +z(2)**2
      if(ssq .gt. sigm2) go to 25
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=0.d0
      zblock(4,j)=0.d0
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by rg2d on PE# ',idproc
c
      return
      end
c
********************************************************************************
c
      subroutine rg4d(nray,iseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 4D gaussian distribution 
c  written by Alex Dragt 7/10/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      sigm2=sigmax**2
      do 20 j=1,nray     
c generate a ray 
   25 call normdv(z,2,ktype,iseed)
      ssq = z(1)**2 + z(2)**2 + z(3)**2 + z(4)**2
      if(ssq .gt. sigm2) go to 25
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by rg4d on PE# ',idproc
c
      return
      end
c
********************************************************************************
c
      subroutine rg6d(nray,iseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 6D gaussian distribution 
c  written by Alex Dragt 7/10/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
      integer, parameter :: ichunk=10000
c
c  working arrays
      dimension z(6),za(6*ichunk)
c
c  proceedure
c
c     write(6,*)'inside rg6d w/ sigmax,sx,sy,st='
c     write(6,*)sigmax,sx,sy,st
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      write(6,*) 'nray=',nray
      write(6,*) 'maxrayp=',maxrayp
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      sigm2=sigmax**2
c
c procedure to generate a small number (~few thousand) of rays:
cryne 1/11/2005      if(nray.gt.5000)goto 23
      if(nray.gt.500000)goto 23
      do 20 j=1,nray
c generate a ray
   21 call normdv(z,3,ktype,iseed)
      ssq=z(1)**2+z(2)**2+z(3)**2+z(4)**2+z(5)**2+z(6)**2
      if(ssq .gt. sigm2) go to 21
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=st*z(5)
      zblock(6,j)=st*z(6)
   20 continue
      goto 999
c
   23 continue
c procedure to generate a large number of rays:
      j=0
   24 continue
      call normdv(za,3*ichunk,ktype,iseed)
      do i=1,6*ichunk-5,6
        ssq=                                                             &
     &  za(i)**2+za(i+1)**2+za(i+2)**2+za(i+3)**2+za(i+4)**2+za(i+5)**2
        if(ssq .gt. sigm2)cycle
        j=j+1
        if(j.gt.nray)exit
c scale and store ray in zblock
        zblock(1,j)=sx*za(i)
        zblock(2,j)=sx*za(i+1)
        zblock(3,j)=sy*za(i+2)
        zblock(4,j)=sy*za(i+3)
        zblock(5,j)=st*za(i+4)
        zblock(6,j)=st*za(i+5)
      enddo
      if(j.lt.nray)goto 24
c
  999 continue
      write(6,*) nray, ' rays generated by rg6d on PE# ',idproc
      return
      end
c
******************************************************************************
c
      subroutine rt2d(nray,iseed,sx,sy,st)
c  subroutine to generate a 2D random uniform distribution on a 2-torus
c  written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray on the unit torus
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=0.d0
      zblock(4,j)=0.d0
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by rt2d on PE# ',idproc
c
      return
      end
c 
********************************************************************************
c
      subroutine rt4d(nray,iseed,sx,sy,st)
c  subroutine to generate a 4D random uniform distribution on a 4-torus
c  written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray on the unit torus
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(3) = cos(arg)
      z(4) = sin(arg)
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by rt4d on PE# ',idproc
c
      return
      end
c 
********************************************************************************
c
      subroutine rt6d(nray,iseed,sx,sy,st)
c  subroutine to generate a 6D random uniform distribution on a 6-torus
c  written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(6)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray on the unit torus
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(3) = cos(arg)
      z(4) = sin(arg)
      call myrand(ktype,iseed,ans)
      arg = twopi*ans
      z(5) = cos(arg)
      z(6) = sin(arg)
c scale and store ray in zblock
      zblock(1,j)=sx*z(1)
      zblock(2,j)=sx*z(2)
      zblock(3,j)=sy*z(3)
      zblock(4,j)=sy*z(4)
      zblock(5,j)=st*z(5)
      zblock(6,j)=st*z(6)
   20 continue
      write(6,*) nray, ' rays generated by rt6d on PE# ',idproc
c
      return
      end
c 
********************************************************************************
c
      subroutine st2d(nray,sx,sy,st)
c  subroutine to generate a 2D systematic uniform distribution on a 2-torus
c written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on circle
      npts=nray
c
      pidiv = twopi/float(npts)
      ntot = 0
      do 20 i=1,npts
c generate ray on unit 2-torus
      ai = float(i-1)
      arg= ai*pidiv
      z(1) = cos(arg)
      z(2) = sin(arg)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(1,ntot)=sx*z(1)
      zblock(2,ntot)=sx*z(2)
      zblock(3,ntot)=0.d0
      zblock(4,ntot)=0.d0
      zblock(5,ntot)=0.d0
      zblock(6,ntot)=0.d0
   20 continue
      nrays=ntot
      write(6,*) ntot, ' rays generated by st2d on PE# ',idproc
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c      
      return
      end
c
********************************************************************************
c
      subroutine st4d(nray,sx,sy,st)
c  subroutine to generate a 4D systematic uniform distribution on a 4-torus
c written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
c
c  working arrays
      dimension z(4)
      dimension c(100),s(100)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      if(maxrayp .gt. 10000) write(6,*)
     & 'error: maxrayp exceeds 10000, subroutines st* need rewriting'
      call myexit
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on each circle
      pow=1.d0/2.d0
      aray=float(nray)
      npts=int(aray**pow)
      if (npts.gt.100) then
      write(6,*)
     & 'error: maxrayp exceeds 10000, subroutines st* need rewriting'
      call myexit
      endif
c
c  compute needed sines and cosines    
c
      pidiv = twopi/float(npts)
      do 10 i=1,npts
      ai = float(i-1)
      arg= ai*pidiv
      c(i) = cos(arg)
      s(i) = sin(arg)
   10 continue
c
c generate rays on unit 4-torus
c
      ntot = 0
      do 20 i=1,npts
      z(1) = c(i)
      z(2) = s(i)
      do 30 j=1,npts
      z(3) = c(j)
      z(4) = s(j)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(1,ntot)=sx*z(1)
      zblock(2,ntot)=sx*z(2)
      zblock(3,ntot)=sy*z(3)
      zblock(4,ntot)=sy*z(4)
      zblock(5,ntot)=0.d0
      zblock(6,ntot)=0.d0
   30 continue
   20 continue
c
      nrays=ntot
      write(6,*) ntot, ' rays generated by st4d on PE# ',idproc
c      
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
      return
      end
c
********************************************************************************
c
      subroutine st6d(nray,sx,sy,st)
c  subroutine to generate a 6D systematic uniform distribution on a 6-torus
c written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
c
c  working arrays
      dimension z(6)
      dimension c(21),s(21)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      if(maxrayp .gt. 10000) write(6,*)
     & 'error: maxrayp exceeds 10000, subroutines st* need rewriting'
      call myexit
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on each circle
      pow=1.d0/3.d0
      aray=float(nray)
      npts=int(aray**pow)
      if (npts.gt.21) then
      write(6,*)
     & 'error: maxrayp exceeds 10000, subroutines st* need rewriting'
      call myexit
      endif
c
c  compute needed sines and cosines    
c
      pidiv = twopi/float(npts)
      do 10 i=1,npts
      ai = float(i-1)
      arg= ai*pidiv
      c(i) = cos(arg)
      s(i) = sin(arg)
   10 continue
c
c generate rays on unit 6-torus
c
      ntot = 0
      do 20 i=1,npts
      z(1) = c(i)
      z(2) = s(i)
      do 30 j=1,npts
      z(3) = c(j)
      z(4) = s(j)
      do 40 k=1,npts
      z(5) = c(k)
      z(6) = s(k)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(1,ntot)=sx*z(1)
      zblock(2,ntot)=sx*z(2)
      zblock(3,ntot)=sy*z(3)
      zblock(4,ntot)=sy*z(4)
      zblock(5,ntot)=st*z(5)
      zblock(6,ntot)=st*z(6)
   40 continue
   30 continue
   20 continue
c
      nrays=ntot
      write(6,*) ntot, ' rays generated by st6d on PE# ',idproc
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c      
      return
      end
c
***************************************************************************
c
      subroutine kv4d(nray,iseed,sx,sy,st)  
c subroutine to compute a KV distribution in 4-D phase space
c
c  written by Alex Dragt 7/15/91
c
      use rays
      include 'impli.inc'
      include 'files.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxrayp) then
      write(6,*) 'error: nray exceeds maxrayp'
      call myexit
      endif
cryne nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nray
      istat(k)=0
      ihist(1,k)=0
      ihist(2,k)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
c select and warmup random number generator
c
      if (iseed .lt. 0) ktype=1
      if (iseed .eq. 0) then
      write(jof,*) ' iseed = 0 in bgen call'
      call myexit
      return
      endif
      if (iseed .gt. 0) ktype=2
      iseed=-iabs(iseed)
      call myrand(ktype,iseed,ans)
c
      do 20 j=1,nray     
c generate a ray 
    5 sumsq=0.d0
      do 10 i=1,4
      call myrand(ktype,iseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
      if(sumsq .eq. 0.d0) goto 5
c scale and store ray in zblock
      rad=sqrt(sumsq)
      zblock(1,j)=sx*z(1)/rad
      zblock(2,j)=sx*z(2)/rad
      zblock(3,j)=sy*z(3)/rad
      zblock(4,j)=sy*z(4)/rad
      zblock(5,j)=0.d0
      zblock(6,j)=0.d0
   20 continue
      write(6,*) nray, ' rays generated by kv4d on PE# ',idproc
c
      return
      end
c
***********************************************************************
c
      subroutine myrand(ktype,iseed,ans)
c
c this subroutine generates random numbers using either myran1
c or myran2
c written by Alex Dragt, 7/28/91
c modified 7/3/98 AJD
c
      double precision ans
c
c      write(6,*) iseed
c
      if(ktype .eq. 1) call myran1(iseed,ans)
      if(ktype .eq. 2) call myran2(iseed,ans)
c
      return
      end
c
***********************************************************************
c
      subroutine myran1(idum,ans)
c
c This subroutine is a minor modification of the FUNCTION ran1(idum)
c given in the article by Press and Teukolosky in Computers in Physics,
c vol. 6, p. 522 (1992).  See also W. Press, S. Teukolsky, W. Vettering,
c and B. Flannery, Numerical Recipes in Fortran, Second edition, page 271,
c Cambridge University Press (1992).
c It requires seeds IDUM that are .lt. 0 to be reset.
c Written by Alex Dragt 7/3/98.
c
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL sans,arg1,AM,EPS,RNMX
      double precision ans
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
c
c warmup generator
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
c
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      arg1=AM*iy
      sans=amin1(arg1,RNMX)
      ans=dble(sans)
c
      return
      END
c
*************************************************************************
c
      subroutine myran2(idum,ans)
c
c This subroutine is a minor modification of the FUNCTION ran2(idum)
c given in the article by Press and Teukolosky in Computers in Physics,
c vol. 6, p. 522 (1992).  See also W. Press, S. Teukolsky, W. Vettering,
c and B. Flannery, Numerical Recipes in Fortran, Second edition, page 272,
c Cambridge University Press (1992).
c It requires seeds IDUM that are .lt. 0 to be reset.
c Written by Alex Dragt 7/3/98.
c
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL sans,arg1,AM,EPS,RNMX
      double precision ans
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     &NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
c
c warmup generator
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
c
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      arg1=AM*iy
      sans=amin1(arg1,RNMX)
      ans=dble(sans)
c
      return
      END
c
********************************************************************************
c
      subroutine normdv(x,n,ktype,iseed)
c  routine to generate 2n independent normal deviates
c  Polar method for normal deviates
c  written by Alex Dragt ca 1986, revised by Alex Dragt 7/15/91
c  modified by Rob Ryne July 28, 2002 to avoid a large number
c  of calls to myrand
c
c  References:
c  1)Knuth, The Art of Computer Programming (Vol. 2, page 117)
c  2)Irving Haber, NRL Memo Report #3705
c
      include 'impli.inc'
      integer, parameter :: ichunk=10000
c
c calling arrays
      dimension x(*)
      real*8, dimension(n) :: u1a,u2a
c
cryne 1/11/2005      if(n.gt.5000)goto 150
      if(n.gt.500000)goto 150
c procedure to generate a small number (~few thousand) of random numbers:
      do 100 k=1,n
   50 call myrand(ktype,iseed,ans)
      u1=ans
      call myrand(ktype,iseed,ans)
      u2=ans
      v1=2.d0*u1-1.d0
      v2=2.d0*u2-1.d0
      s=v1**2 + v2**2
      if(s .ge. 1.d0)goto 50
      arg=sqrt(-2.d0*log(s)/s)
      ans1=v1*arg
      ans2=v2*arg
      x(k)=ans1
      x(k+n)=ans2
  100 continue
      return
c procedure to generate a large number of random numbers:
  150 continue
c     write(6,*)'using packed version of normdv'
      j=0
  200 continue
      call random_number(u1a)
      call random_number(u2a)
      do i=1,ichunk
        v1=2.d0*u1a(i)-1.d0
        v2=2.d0*u2a(i)-1.d0
        s=v1**2 + v2**2
        if(s .gt. 1.d0)cycle
        j=j+1
        if(j.gt.n)exit
c store results:
        arg=sqrt(-2.d0*log(s)/s)
        ans1=v1*arg
        ans2=v2*arg
        x(j)=ans1
        x(j+n)=ans2
      enddo
      if(j.lt.n)goto 200
c
      return
      end
c
***********************************************************************
c
      subroutine tic(p)
c
c Translation of initial conditions.
c This subroutine produces a translation in 6-dimensional phase space.
c The parameters p(j) are used to specify translations deltaz(j)
c according to the relations deltaz(j)=p(j).
c The suffixes 'i' and 'f' refer to 'initial' and 'final' respectively.
c Written by Alex Dragt, Fall 1986
c
      use rays
      include 'impli.inc'
      dimension p(6)
c
      do 100 i=1,nraysp
c check to see if i'th ray has been lost
      if (istat(i).ne.0) goto 100
c if not, copy the i'th ray out of zblock
      do 110 j=1,6
  110 zi(j)=zblock(j,i)
c
c transform this ray
c
      do 10 j=1,6
   10 zf(j)=zi(j)+p(j)
c
c put the transformed ray back into zblock
      do 120 j=1,6
  120 zblock(j,i)=zf(j)
c
  100 continue
      return
      end
c
***********************************************************************
c routine to setup/warmup the F90 random number generator 
      subroutine f90ranset(iseed)
      implicit double precision(a-h,o-z)
      call random_seed(size=iseedsize)
c     write(6,*)'random number seed size equals ',iseedsize
      call f90ranwarm(iseedsize,iseed)
c     write(6,*)'here are a few numbers from random_number:'
c     do n=1,50
c       call random_number(x)
c       write(6,*)x
c     enddo
      return
      end
 
      subroutine f90ranwarm(iseedsize,iseed)
      implicit double precision(a-h,o-z)
      dimension iputarray(iseedsize)
      if(iseedsize.eq.1)return
      iputarray(1)=iseed
c multiply the seed by a random number
c note: the generator has not been initialized, so this is
c a bit flaky
      do n=2,iseedsize
        call random_number(x)
        iputarray(n)=iseed*x
      enddo
c set the seed:
      call random_seed(put=iputarray)
      return
      end

c
c end of file
