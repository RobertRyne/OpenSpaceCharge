c*********   USER subroutines pkg  *******
      subroutine allist(lun)
c
c       allist makes an alphabetical LaTeX table of all
c       type codes in the current vesion of marylie
c       T. Mottershead   LANL  AT-3  16 Jan 92
c-------------------------------------------------
      include 'codes.inc'
      character*8  tcode(360), word
      dimension nseq(360), npars(360), kind(360)
c
c     fill unsorted tcode array, stopping on null string for each group
c
      long = 40
      max = 0
      do 20 j = 1, 9
         do 10 n=1,40
         word = ltc(j,n)
         if(ichar(word(1:1)).eq.0) go to 20
         max = max + 1
         tcode(max) = word
         npars(max) = nrp(j,n)
         kind(max) = j
  10     continue
  20  continue
c
c      alphabetize and write Latex table
c
      call lexort(max,tcode,nseq)
      call headoc(lun,'MaryLie Type Codes')
      call tabdef(lun,5)
      write(lun,27)
  27  format(' No. & Code  & Kind & NRP  & \\')
      write(lun,*) '  \hline'
      do 50  k = 1,max
      jj = nseq(k)
      write(lun,37) k, tcode(jj), kind(jj), npars(jj)
  37  format(i6,' & ',a,' &',i3,' &',i3,' &  \\')
      if(mod(k,long).eq.0) then
        call tabend(lun)
        write(lun,*) ' \newpage'
        call tabdef(lun,5)
        write(lun,27)
        write(lun,*) '  \hline'
      endif
  50  continue
      call tabend(lun)
      call endoc(lun)
      return
      end
c---------------------------------------------------
      subroutine arcs(nseg,ktyp,path,bend,angle,delta,max,xx,yy,key,nn)
c
c        Generates (xx(i),yy(i),i=1,nn) for drawing a sequence of
c        connected arcs. The curve starts from x(1),y(1) at input angle,
c        and ends at x(nn),y(nn) headed in direction of final angle.
c        nseg = number of path segments k=1,nseg
c        path(k) = arclength of kth segment
c        bend(k) = bend angle in degrees of kth segment
c        ktyp(k) = ID or type code for kth segment
c        delta = degrees of bend per plot point (max allowed)
c        max = maximum number of output points allowed
c
c        C. T. Mottershead  24 Aug 99
c-------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension path(*),bend(*),ktyp(*),xx(max),yy(max),key(max)
      pi = 3.14159265358979d0
      radeg = pi/180.0d0
      bmin = 0.001
      one = 1.0d0
      xu = xx(1)
      yu = yy(1)
      key(1) = ktyp(1)
      lout = 44
      nn = 0
c
c       loop over all path segments
c
      do 200 k = 1,nseg
        theta = radeg*angle
        angle = angle + bend(k)
        xs = path(k)*dcos(theta)
        ys = path(k)*dsin(theta)
        write(lout,207) k,ktyp(k),path(k),bend(k),angle,xs,ys
  207   format(' Arc:',2i4,5f12.4)
c
c     straight segment
c
        if(dabs(bend(k)).lt.bmin) then
          nn = nn+1
          if(nn.gt.max) return
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          nn = nn+1
          if(nn.gt.max) return
          xu = xu + xs
          yu = yu + ys
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          np = 1
c          write(lout,97) nn,key(nn),xx(nn),yy(nn),xs,ys
          go to 190
        endif
c
c      curved segment (starts with duplicate of previous point)
c
        nn = nn+1
        if(nn.gt.max) return
        xx(nn) = xu
        yy(nn) = yu
        key(nn) = ktyp(k)
        phi = radeg*bend(k)
        radius = path(k)/phi
        np = nint(bend(k)/delta)
        if(np.lt.0) np = -np
        if(np.eq.0) np = 1
        fnp = float(np)
        alpha = phi/fnp
        xs = xs/fnp
        ys = ys/fnp
        cosa = dcos(alpha)
        sina = dsin(alpha)
        sfac = sina/alpha
        cfac = (one-cosa)/alpha
        x0 = xs
        y0 = ys
c
c     update plot points in arc
c
        do 100 j = 1,np
          dx = sfac*xs - cfac*ys
          dy = sfac*ys + cfac*xs
          xu = xu + dx
          yu = yu + dy
          nn = nn+1
          if(nn.gt.max) return
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          xs = cosa*x0 - sina*y0
          ys = cosa*y0 + sina*x0
          x0 = xs
          y0 = ys
  100   continue
  190 continue
  200 continue
      return
      end
c
ccccccccccccccccc     Blocks     cccccccccccccccccccccccccccccc
c
      subroutine blocks(nep,lun)
      use beamdata
      include 'impli.inc'
      include 'linbuf.inc'
      zero = 0.0d0
      ap = zero
      zbeg = zero
      nn = 0
      write(lun,12) nep,nrays,npos,nang,brho,beta,gamm1
      write(lun,13) zbeg, ap, strong(1), nn,nn, cname(1)
  12  format(4i5,3f16.7)
  13  format(3f16.7,2i5,2x,a8)
c
c         write envelope file top section
c
      do 20 nn = nbgn, nend
        ltyp = lintyp(nn)
        zbeg = total
        total = total + elong(nn)
        cmtot = 100.0*total
        cml = 100.0*elong(nn)
        glp = elong(nn)*strong(nn)
        gltot = gltot + abs(glp)
        rcm = 100.0*radius(nn)
        kk = nn
        if(ltyp.ne.0) then
            ap = radius(nn)
            write(lun,13)  zbeg, ap, strong(nn), nn, ltyp, cname(nn)
            write(lun,13) total, ap, strong(nn), nn, ltyp, cname(nn)
            kk = nn + 1
        endif
        ap = 0.0
        write(lun,13) total, ap, zero, nn, ltyp, cname(kk)
  20  continue
      return
      end
c--------------------------------------

      subroutine caedit(max, nl, text)
c
c     caedit is a simple character array editing routine.
c     text(i), i=1,nl is the array of character strings to be edited.
c          The length of the strings is determined internally.
c          All strings are the same length.
c     nl = actual number of character strings, which may be increased
c          or decreased by this routine.
c     max = maximum number of character strings allowed by the calling
c           program.
c
c           T. Mottershead/ LANL  AT-3 /  24 Jan 91
c----------------------------------------------------------
      character*(*) text(max)
      character*74  card
c
c         print array
c
  10  do 20 j=1, nl
      card = text(j)
      write(6,13) j,card
  13  format(i3,a)
  20  continue
c
c         command input
c
  30  write(6,32)
  32  format(' enter command: +N = insert new line N, 0 = quit,',
     & ' -N = delete line N')
      read(5,*) kl
      if(kl.eq.0) return
c
c     delete line N
c
      if(kl.lt.0) then
         nd = -kl
         if(nd.gt.nl) go to 30
         last = nd
         card = ' '
         write(6,37) last
  37     format(' enter last line to delete <',i3,'>:')
         read(5,64) card
         if(card.ne.' ') then
            read(card,*) last
         endif
  40     nl = nl-1
         do 50 n = nd,nl
         text(n) = text(n+1)
  50     continue
         if(last.le.nd) go to 10
         last = last-1
         go to 40
      endif
c
c       add new line N
c
      if(kl.gt.0) then
         if(kl.gt.nl) kl = nl+1
         write(6,*) ' Enter new line(s) <CR>=print'
  60     card = ' '
         read(5,64) card
  64     format(a)
         if(card.eq.' ') go to 10
         nl = nl+1
         do 90 n = nl-1, kl, -1
         text(n+1) = text(n)
  90     continue
         text(kl)=card
         kl = kl+1
         go to 60
      endif
      go to 10
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine endoc(lun)
c
c       writes LaTeX document end to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c--------------------------------------------------
      write(lun,*) ' \end{document}'
      return
      end
c
***********************************************************************
c
      subroutine filin
c
c       filin's job is to fill in /linbuf/ from the contents of
c       a loop. Derived from wcl, written by Alex Dragt, 23 August 1988
c       Based on the SRs cqlate and pmif. WCL modified by F. Ner
c       Dec 90 to fill LINBUF. That function of wcl save is this routine
c       April 91, by T. Mottershead.
c
      use beamdata
      use acceldata
      use lieaparam, only : monoms,monom1,monom2
      include 'impli.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
c
c common blocks
c
      include 'codes.inc'
      include 'parset.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
c
      dimension ltype(50)
      data (ltype(j),j=1,25)/
     &   0,2,2,2,2,
     &   2,3,2,1,0,
     &   0,0,4,5,7,
     &  -1,7,1,7,6,
     &   7,7,7,7,7/
      save ltype
c
c    set parameter set to use (old=9, now pass in from caller CTM 5/00)
c
c     npsu = 9
c
c  start routine by checking to see if a loop exists
c
      if(nloop.le.0) then
         write(jodf,*)  '  error in filin: no loop in labor'
         write(jof,*)  '  error in filin: no loop in labor'
         return
      endif
c
c  initialize linbuf pointers:
c
      nbgn = 1
      nend = 0
      npst = 0
c
c  contents of biglist
c
c      write(jodf,*) joy,' = joy in filin'
      do 100 jj=1,joy
      mmv = mim(jj)
c       element
        if(mmv.lt.0) then
          ktyp = 1
          mnu = -mmv
c       user supplied element
        else if(mmv.gt.5000) then
          ktyp = 1
          mnu = mmv-5000
c       lump
        else
          ktyp = 3
        endif
c       write(jodf,*) '  filin: ',jj,'=jj',mmv,'=mim(jj)',mnu,'=mnu',
c     * ktyp,'=ktyp'
c procedure for a menu item
      if(ktyp.ne.1) go to 100
         imax=nrp(nt1(mnu),nt2(mnu))
         if(imax.eq.0)goto 100
c
c   reset default radius from pset(npsu) (T. Mottershead 23 Jan 91)
c   keep all values from the same pset for other purposes, e.g. FOV calc
c   (CTM 12/97)
c   CTM 11/99: save all instances of npsu in linbuf for future use
c
      if((nt1(mnu).eq.3).and.(nt2(mnu).eq.npsu)) then
         aapar = pmenu(1+mpp(mnu))
         bbpar = pmenu(2+mpp(mnu))
         ccpar = pmenu(3+mpp(mnu))
         ddpar = pmenu(4+mpp(mnu))
         eepar = pmenu(5+mpp(mnu))
         ffpar = pmenu(6+mpp(mnu))
         npst = npst + 1
         aap(npst) = aapar
         bbp(npst) = bbpar
         ccp(npst) = ccpar
         ddp(npst) = ddpar
         eep(npst) = eepar
         ffp(npst) = ffpar
      endif
c
c  Put simple element in linbuf (F. Neri 12/18/1990)
c
      if (nt1(mnu).eq.1) then
        lt = ltype(nt2(mnu))
        if(lt.eq.1.or.lt.eq.0) then
          nend = nend+1
          cname(nend) = lmnlbl(mnu)
          lintyp(nend) = lt
          mltyp1(nend) = nt1(mnu)
          mltyp2(nend) = nt2(mnu)
          elong(nend) = pmenu(1+mpp(mnu))
          radius(nend) = aapar
          strong(nend) = 0.d0
c
c        normal quad
c
          if(nt2(mnu).eq.9) strong(nend) = pmenu(2+mpp(mnu))
c
c   quick & dirty fix for new cfqd (think again about when psets are exe
c
          if(nt2(mnu).eq.18) then
               ids = nint(pmenu(2+mpp(mnu)))
               strong(nend) = pst(1,ids)
               pssext(nend) = pst(3,ids)
               psoctu(nend) = pst(5,ids)
          endif
        endif
      endif
 100  continue
      return
      end
c
ccccccccccccccccc     RMQUAD    ccccccccccccccccccccccccccccccc
c
      subroutine fvscan(wx,wy,aper,theta,xfov,yfov,kx,ky)
      implicit double precision (a-h,o-z)
      parameter (maxz=2000)
      common/csrays/jp,mode,cx(maxz),sx(maxz),cy(maxz),sy(maxz)
     & ,z(maxz)
      write(6,*)'WARNING (fvscan): this routine has a different'
      write(6,*)'declaration for common/csrays/ than in sincos.inc'
      write(6,*)'Rob Ryne 7/23/2002'
c
      eps = 1.0d-12
      xfov = 1.0d6
      yfov = xfov
      do 100 k = 1,jp
      denx = cx(k)+wx*sx(k)
      deny = cy(k)+wy*sy(k)
      if(dabs(denx).lt.eps) denx = eps
      if(dabs(deny).lt.eps) deny = eps
c      type *,k,denx,deny
      biga = aper/denx
      bigb = aper/deny
      xbar = sx(k)*theta/denx
      ybar = sy(k)*theta/deny
      xsiz = dabs(biga+xbar)
      if(xsiz.lt.xfov) then
         xfov = xsiz
         kx = k
      endif
      xsiz = dabs(biga-xbar)
      if(xsiz.lt.xfov) then
         xfov = xsiz
         kx = k
      endif
      ysiz = dabs(bigb+ybar)
      if(ysiz.lt.yfov) then
         yfov = ysiz
         ky = k
      endif
      ysiz = dabs(bigb-ybar)
      if(ysiz.lt.yfov) then
         yfov = ysiz
         ky = k
      endif
  100 continue
      xfov = 200.0d0*xfov
      yfov = 200.0d0*yfov
      return
      end
c--------------------------------------------------------------
      subroutine genray(wx,wy,lunz,npos,dx,dy,nang,dxp,dyp)
c
c        generates matched ray initial conditions for a map
c        and a set of scattered rays
c        C. T. Mottershead  LANL AOT-1  16 April 1996
c-------------------------------------------------------
      use rays
      implicit double precision (a-h,o-z)
      logical ltbl
      include 'parset.inc'
      dimension pp(6)
      write(6,*) igen,' = igen,  wx,wy=', wx,wy
      zero = 0.0d0
      ltbl = .false.
      if(lunz.gt.0) ltbl = .true.
      do 10 j = 1,6
      pp(j) = zero
  10  continue
      nn = 0
      do 100 i = 1, npos
         fi = npos - i + 1
         x = fi*dx
         xpref = x*wx
         y = fi*dy
         ypref = y*wy
         pp(1) = x
         pp(3) = y
c
c        positive scattering angles
c
         if(nang.ge.1) then
         do 50 j = 1, nang
             fj = nang - j + 1
             pp(2) = xpref + fj*dxp
             pp(4) = ypref + fj*dyp
             if(ltbl) write(lunz,57) pp
             nn = nn + 1
             do 40 k = 1, 6
                zblock(nn,k) = pp(k)
  40         continue
  50     continue
         endif
c
c       central matched ray
c
      pp(2) = xpref
      pp(4) = ypref
      if(ltbl) write(lunz,57) pp
          nn = nn + 1
          do 60 k = 1, 6
             zblock(nn,k) = pp(k)
  60      continue
  57  format(1x,4(1pe13.5),2(1pe12.4))
c
c        negative scattering angles
c
         if(nang.ge.1) then
         do 80 j = 1, nang
             fj = j
             pp(2) = xpref - fj*dxp
             pp(4) = ypref - fj*dyp
             if(ltbl) write(lunz,57) pp
             nn = nn + 1
             do 70 k = 1, 6
                zblock(nn,k) = pp(k)
  70         continue
  80     continue
         endif
 100  continue
c
c      on axis rays
c
      pp(1) = zero
      pp(3) = zero
      if(nang.ge.1) then
      do 120 j = 1, nang
          fj = nang - j + 1
          pp(2) = fj*dxp
          pp(4) = fj*dyp
          if(ltbl) write(lunz,57) pp
             nn = nn + 1
             do 110 k = 1, 6
                zblock(nn,k) = pp(k)
 110         continue
 120  continue
      endif
      pp(2) = zero
      pp(4) = zero
      if(ltbl) write(lunz,57) pp
         nn = nn + 1
         do 140 k = 1, 6
            zblock(nn,k) = pp(k)
 140     continue
      nrays = nn
      write(6,147) npos,nang,nrays
 147  format(' #GENRAY:',i6,' positions',i6,' angles',i6,' rays')
      return
      end
c
***********************************************************************
c
      subroutine headoc(lun,title)
c
c       writes LaTeX document header to LUN from fortran programs.
c
c      T. Mottershead  AT-3   16 Jan 91
c--------------------------------------------------
      character title*(*)
c
c           Write Document opening
c
      write(lun,*) ' \documentstyle[12pt]{article} '
      write(lun,*) ' \setlength{\textwidth}{6.5in}'
      write(lun,*) ' \setlength{\oddsidemargin}{0in} '
      write(lun,*) ' \setlength{\topmargin}{-0.5in}'
      write(lun,*) ' \setlength{\textheight}{9in} \pagestyle{empty}'
      write(lun,*) ' \begin{document}'
      write(lun,*) ' '
      write(lun,*) ' {\large '
      write(lun,14) title
  14  format(' \begin{center} {\bf ',a,' }  \\
     & \today    \end{center} }')
      return
      end
c***************************************************************
      subroutine kinema(bmass,ener,pmom)
      use beamdata
      implicit double precision (a-h,o-z)
cryne 7/23/2002      common/parm/brho,c,gamma,gamm1,beta,charge,sl,ts,pbeam
      two = 2.0d0
      light = 299792458
      fmev = 1.0d-6*float(light)
      pmom = fmev*brho
      gambet = dsqrt(gamm1*(gamm1+two))
      bmass = pmom/gambet
      ener = gamm1*bmass
      return
      end
c-----------------------------------------------------------------
      subroutine listyp(lun)
c
c       listyp  lists the type codes in the current vesion of marylie
c       T. Mottershead   LANL  AT-3  17 Jan 91
c-------------------------------------------------
      include 'codes.inc'
      character*8  tcode(9), word
      dimension numcmd(9), nord(9), jlo(2),jhi(2)
      data nord /1,4,7,8,9,2,3,5,6/
      save nord   !cryne 7/23/2002
      jlo(1)=1
      jhi(1)=5
      jlo(2)=6
      jhi(2)=9
c
c     scan for null string to indicate number of commands
c
      ntot = 0
      do 20 j = 1, 9
         max = 0
         do 10 n=1,40
         word = ltc(j,n)
         if(ichar(word(1:1)).eq.0) go to 15
         max = max + 1
  10     continue
  15     numcmd(j) = max
         ntot = ntot + max
  20  continue
c
c        report total
c
      write(lun,*) ' Marylie Type Codes '
      write(lun,*)
      write(lun,*) ntot,' type codes in this version'
c
c        writeout
c
      do 90  kk = 1,2
      jmin = jlo(kk)
      jmax = jhi(kk)
      nmax = 0
      do 30 j=jmin, jmax
      mm = numcmd(nord(j))
      if(mm.gt.nmax) nmax = mm
  30  continue
      write(lun,*)
      write(lun,33) (nord(k),k=jmin,jmax)
  33  format(' index ',2x,5('nrp (G',i1,')',4x))
      do 50 n = 1 ,nmax
         do 40 j=jmin, jmax
         tcode(j) = '        '
         jju = nord(j)
         if(n.le.numcmd(jju)) tcode(j) = ltc(jju,n)
  40  continue
      write(lun,43) n,(nrp(nord(j),n), tcode(j), j=jmin,jmax)
  43  format(i4,5x,5(i2,2x,a8))
  50  continue
      write(lun,31) (numcmd(nord(j)), j=jmin, jmax)
  31  format(' total:',i8,4i12)
  90  continue
      return
      end
c
c*******************************************************
      subroutine mlfov(wx,wy,aper,phiref,philim,lunf,ufov,kfov)
      implicit double precision (a-h,o-z)
      dimension ufov(*),kfov(*)
      parameter (maxz=2000)
      common/csrays/jp,mode,cx(maxz),sx(maxz),cy(maxz),sy(maxz)
     & ,z(maxz)
      logical ltbl, last
c
      write(6,*)'WARNING (mlfov): this routine has a different'
      write(6,*)'declaration for common/csrays/ than in sincos.inc'
      write(6,*)'Rob Ryne 7/23/2002'
c
c     write(99,*) '  MLFOV call:',jp,aper,philim,phiref,lunf
      ltbl = .false.
      last = .false.
      if(lunf.gt.0) ltbl = .true.
c
c    scan given sin/cosine rays for extrema
c
c     type *, jp,'=jp in mlfov (length of S/C buffer)'
      sxmax = 0.0d0
      cxlo = 1.0d30
      cxhi = -cxlo
      symax = 0.0d0
      cylo = 1.0d30
      cyhi = -cylo
      do 50 k = 1,jp
      cc = cx(k)
      ss = sx(k)
      xmat = cc + wx*ss
      if(xmat.gt.cxhi) cxhi = xmat
      if(xmat.lt.cxlo) cxlo = xmat
      sab = dabs(ss)
      if(sab.gt.sxmax) sxmax = sab
      cc = cy(k)
      ss = sy(k)
      ymat = cc + wy*ss
      if(ymat.gt.cyhi) cyhi = ymat
      if(ymat.lt.cylo) cylo = ymat
      sab = dabs(ss)
      if(sab.gt.symax) symax = sab
  50  continue
c      type *, '  C+wS:',cxlo,cxhi,cylo,cyhi
c
c       solve for envelope max
c
      pxm = 1000.0d0*aper/sxmax
      pym = 1000.0d0*aper/symax
      if(ltbl) write(lunf,73) sxmax,symax, pxm, pym
  73  format(' Max S: ',2f13.5,'    Max-phi(mR):',2f10.3)
      phi = 0.0d0
      phimax = pxm
      if(ltbl) write(lunf,506)
  506 format(' phi(mR)  xfov     yfov     kx   ky',
     & 6x,'ysize+/-',10x,'xsize+/-')
c
c      scan for FOV at the reference angle
c
      xrefov = 0.0d0
      yrefov = 0.0d0
      if(phiref.lt.phimax) then
         theta = 1.0d-3*phiref
         call fvscan(wx,wy,aper,theta,xrefov,yrefov,kx,ky)
      endif
c
c      initialize search over all angles
c
      kf = 0
      kx = 0
      ky = 0
      dphi = 0.05
      theta = 0.0d0
      eps = 1.0d-12
  510 continue
      kf = kf+1
      kxo = kx
      kyo = ky
      call fvscan(wx,wy,aper,theta,xfov,yfov,kx,ky)
c
c       check size of other axis ellipse at limit points (kx,ky)
c
      dyx = cy(kx) + wy*sy(kx)
      if(dabs(dyx).lt.eps) dyx = eps
      ayxm = 100.0d0*(aper - sy(kx)*theta)/dyx
      ayxp = 100.0d0*(aper + sy(kx)*theta)/dyx
      dxy = cx(ky) + wx*sx(ky)
      if(dabs(dxy).lt.eps) dxy = eps
      axym = 100.0d0*(aper - sx(ky)*theta)/dxy
      axyp = 100.0d0*(aper + sx(ky)*theta)/dxy
      phi = 1000.0d0*theta
      ango = ang
      xapo = xap
      ang = phi
      xap = xfov
      if(ltbl) write(lunf,107) phi,xfov,yfov,kx,ky
     & ,ayxm,ayxp,axym,axyp
 107  format(f7.2,2f9.4,2i5,4f9.3)
      if(kf.eq.1) xzero = xfov
      if((iabs(kxo-kx).gt.80).and.(iabs(ky-kyo).gt.80)) then
         xfb = xapo
         phib = ango
         kxb = kx
         kyb = ky
      endif
      if(last) then
        radphi = phib/philim
        bratio = zrad(radphi)
        rmax = phimax/philim
        radmax = zrad(rmax)
c        type *, 'break ratios: phi,z ',phi, radphi, bratio
c        type *, 'endpoint ratios: phi,z ',phimax,rmax,radmax
        pipe = 200.0d0*aper
        ufov(1) = pipe
        ufov(2) = xzero
        ufov(3) = xfb
        ufov(4) = phib
        ufov(5) = phimax
        ufov(6) = bratio
        ufov(7) = radmax
        ufov(8) = sxmax
        ufov(9) = symax
        ufov(10)= pym
        ufov(11)= xrefov
        ufov(12)= yrefov
        kfov(1) = kxb
        kfov(2) = kyb
      endif
      phi = phi + dphi
      theta = 1.0d-3*phi
      if(last) return
      if(phi.lt.phimax) go to 510
      phi = phimax
      theta = 1.0d-3*phi
      last = .true.
      go to 510
      end
c--------------------------------------------------------------
      subroutine rmquad(zlen,grad,rm)
c
c     computes matrix rm for a quadrupole of length zlen meters
c     and with field gradient of grad tesla/meter
c     if grad>0, the quad is horizontally focusing, vertically defocusin
c     if grad<0, the quad is vertically focusing, horizontally defocusin
c     If grad=0, it is simply a drift.
c       T. Mottershead LANL AT-3, 10/89, copied from MARYLIE routines
c       FQUAD and DQUAD by D. Douglas, 1982.
c------------------------------------------------------------------
      use beamdata
      implicit double precision (a-h,o-z)
c
c   beam and kinematic constants:
cryne 7/23/2002      common/parm/brho,c,gamma,gamm1,beta,achg,sl,ts
      double precision rm(6,6)
c
c     initialize rm to the identitiy matrix
c
      do 40 i=1,6
      do 20 j=1,6
      rm(i,j)=0.0
  20  continue
      rm(i,i)=+1.0d0
  40  continue
c
c     add drift terms to rm
c
      rm(1,2)= zlen
      rm(3,4)= zlen
      rm(5,6)= zlen/(gamma*beta)**2
c
c      drift return
c
      if(grad.eq.0.0) return
c
c    positive quad (horizontally focusing)
c
      square = grad/brho
      if(square.gt.0.0) then
         ifp = 1
         jfp = 2
         idp = 3
         jdp = 4
      endif
c
c    negative quad (horizontally defocusing)
c
      if(square.lt.0.0) then
         idp = 1
         jdp = 2
         ifp = 3
         jfp = 4
         square = -square
      endif
c
c     coefficients for transverse rm
c
      wavek=dsqrt(square)
      zk=zlen*wavek
      coshkz=(dexp(zk)+dexp(-zk))/(2.0d0)
      sinhkz=(dexp(zk)-dexp(-zk))/(2.0d0)
      coskz=dcos(zk)
      sinkz=dsin(zk)
c
c      matrix in the focusing plane
c
      rm(ifp,ifp) = coskz
      rm(ifp,jfp) = sinkz/wavek
      rm(jfp,jfp) = coskz
      rm(jfp,ifp) = -wavek*sinkz
c
c      matrix in the defocusing plane
c
      rm(idp,idp) = coshkz
      rm(idp,jdp) = sinhkz/wavek
      rm(jdp,jdp) = coshkz
      rm(jdp,idp) = wavek*sinhkz
      return
      end
c
***********************************************************************
c
      subroutine setcor(kor,wx,wy)
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'taylor.inc'
      include 'map.inc'
      include 'parset.inc'
      include 'usrdat.inc'
      zero = 0.0d0
      two = 2.0d0
c
c     No corr: Straight in
c
      if(kor.eq.0) then
         wx = zero
         wy = zero
      endif
c
c     Fourier corr
c
      if(kor.eq.1) then
         wx = -tmh(1,1)/tmh(1,2)
         wy = -tmh(3,3)/tmh(3,4)
      endif
c
c     Chromatic corr for identity lens
c
      if(kor.eq.2) then
         wx = -th(38)/(two*th(53))
         wy = -th(70)/(two*th(76))
      endif
c
c     Chromatic corr in general from Taylor map
c
      if(kor.eq.3) then
         wx = -tumat(12,1)/tumat(17,1)
         wy = -tumat(21,3)/tumat(24,3)
      endif
c
c     general correlation from ucalc(N) and N+1
c
      if(kor.lt.0) then
        nc = -kor
        if(nc.gt.250) return
        wx = ucalc(nc)
        wy = ucalc(nc+1)
      endif
      return
      end
c==========================================
      subroutine tabdef(lun,ncol)
c
c       writes LaTeX table headers to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c--------------------------------------------------
      write(lun,*) ' \begin{center}'
      write(lun,*) ' \begin{tabular}{',('|c',j=1,ncol),'|}'
      write(lun,*) ' \hline'
      return
      end
c***************************************************************
      subroutine tabend(lun)
c
c       writes LaTeX table end to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c--------------------------------------------------
      write(lun,*) ' \hline \end{tabular} \end{center}'
      return
      end
c***************************************************************
      subroutine tugen(fa,fm)
c  pcmap routine to print m,f3,f4 and t,u.
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
ctm modified to just generate the taylor map in 'taylor.inc'
c
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'files.inc'
      include 'expon.inc'
      include 'pbkh.inc'
      include 'taylor.inc'
c
      dimension fa(monoms),fm(6,6),fsav(monoms)
      dimension t(monoms),u(monoms),u2(monoms)
      write(6,*)'inside routine tugen'
c
c  prepare for higher order matrix generation
c
      do 10 jj = 1, monoms
      fsav(jj) = fa(jj)
  10  continue
      call brkts(fa)
c
c  procedure for generating t-matrix
c
        do 35 i=1,6
        call xform(pbh(1,i),2,fm,i-1,t)
        do 36 n=7,27
        tumat(n,i) = t(n)
   36   continue
   35   continue
c
c  procedure for generating U-matrix
c
         do  44 i=1,6
         call xform(pbh(1,i),3,fm,i-1,u)
         call xform(pbh(1,i+6),3,fm,1,u2)
         do  45 n=28,83
         u(n)=u(n)+u2(n)/2.d0
         tumat(n,i) = u(n)
   45    continue
   44    continue
      return
      end
c
c------------------------------------------------------
c
      subroutine u8out(lun,level)
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'usrdat.inc'
      include 'sigbuf.inc'
      character*3 tag
      write(lun,21)
  21  format(/' Output beam parameters:',/,'  x=u(1)    ',
     &'   alfx=u(2)    px=u(3)      y=u(4)       alfy=u(5) ',
     &'   py=u(6)   ')
      write(lun,25) (ucalc(j), j=1,6)
  25  format(6(1pe13.5))
      write(lun,22)
  22  format(6x,'beta',15x,'emittance',10x,'focus',14x,'cross-focus')
      tag = ' X:'
      write(lun,24) tag, (j,ucalc(j),j=7,10)
  24  format(a3,i3,'>',1pe14.7,3(i4,'>',1pe14.7))
      tag = ' Y:'
      write(lun,24) tag, (j,ucalc(j),j=11,14)
      if(level.lt.1) return
c
c      also write aberration coefficients for level 1
c
      write(lun,27)
  27  format(/,' Chromatic and RMS Geometric Aberration Coefficients:',
     &/,'     dFx/dp=u(15)      dFy/dp=u(16)      K2(R/m^2)=u(17)',
     &'   K3(R/m^3)=u(18)')
      write(lun,29) (ucalc(j),j=15,18)
  29  format(4(1pg18.7))
      if(level.lt.2) return
c
c     also write dsigma/de for level 2
c
      write(lun,23)
  23  format(/' Energy dependence of sigma matrix:',/,'  ds11=u(21)',
     &'   ds12=u(22)   ds22=u(23)   ds33=u(24)   ds34=u(25)',
     &'   ds44=u(26)')
      write(lun,25) (ucalc(j),j=21,26)
c      write(lun,*) ' Sigma matrix:'
c      do 201 i=1,4
c      write(lun,26) (sigf(i,j),j=1,4)
c  201 continue
c 26  format(3x,4(1pe15.7))
c         write(lun,*) ' Energy derivative of Sigma matrix:'
c      do 202 i=1,4
c         write(lun,26) (dsig(i,j),j=1,4)
c  202 continue
      return
      end
c
***********************************************************************
c
      subroutine zmono(nopt,zz,zvec)
      use lieaparam, only : monoms
      include 'impli.inc'
      include 'expon.inc'
      include 'vblist.inc'
      dimension zz(6)
cryne August 4, 2004
cryne setting the following to 209 to reproduce old results.
cryne Extend later to 5th order
c     dimension zvec(monoms)
      dimension zvec(209)
      write(6,*)'inside routine zmono'
c
c      do 20 n = 1, monoms
c      j = 1
c      kode = expon(j,n)
c  18  continue
c      j = j + 1
c      kode = 10*kode + expon(j,n)
c      if(j.lt.6) go to 18
c      write(89,19) n,kode,(expon(j,n),j=1,6),(vblist(k,n),k=1,4)
c  19  format(2i8,2x,6i2,'  vb:',4i2)
c  20  continue
      one = 1.0d0
      nord = 2
cryne August 4, 2004
cryne setting the following to 209 to reproduce old results.
cryne Extend later to 5th order
c     do 100 n = 1, monoms
      do 100 n = 1, 209
      if(n.eq.28) nord = 3
      if(n.eq.84) nord = 4
      prod = one
        do 50 j = 1,nord
           prod = prod*zz(vblist(j,n))
  50    continue
      zvec(n) = prod
 100  continue
      return
      end
      double precision function zrad(radphi)
      implicit double precision (a-h,o-z)
      common/iters/it
      zero = 0.0d0
      if(radphi.le.zero) then
         zrad = zero
         return
      endif
      tol = 1.0d-9
      derr = 1.0d9
      one = 1.0d0
      two = 2.0d0
      bb = 9.0d0*dlog(10.0d0)
      aa = one/bb
      bb = bb/two
c      type *, ' aa=',aa,'  bb=',bb
      maxit=50
      y = radphi**2
      zest = y
      do 80 j=1,maxit
      it = j
      vv = one+aa*dlog(zest)
      yc = zest*(vv)**2
      err = yc-y
c      type *,j,zest,err
      errold = derr
      derr = dabs(err)
      if(derr.lt.tol) then
        if(derr.ge.errold) go to 100
      endif
      zest = (vv*zest+bb*y)/(vv*(one+bb*vv))
  80  continue
 100  zrad = zest
      return
      end
c-----------------------------------------------
      subroutine trcdmp(ltran,ltrace)
      use beamdata
      include 'impli.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
      include 'sigbuf.inc'
      include 'usrdat.inc'
      ntrc = ltrace-1
      write(ltran,*)  '$data'
      do 500 nn = nbgn, nend
      nnu = nn + ntrc
      tlong = 1000.0*elong(nn)
c
c     drift
c
      if(lintyp(nn).eq.0)  then
         ityp = 1
         write(ltran,488) nnu,cname(nn),nnu,ityp,nnu,tlong
      endif
c
c     quad
c
      if(lintyp(nn).eq.1)  then
         ityp = 3
         write(ltran,488) nnu,cname(nn),nnu,ityp,nnu,strong(nn),tlong
      endif
 488  format(1x,'cmt(',i3,')=''',a,''' nt(',i3,')=',i3,', a(1,',i3
     & ,')=',2(1pg17.9,','))
 500  continue
      write(ltran,*)  '$end'
      return
      end
ccccccccccccccc++++++++++++++++++++++++++++c
      subroutine adump(ldump)
      use beamdata
      include 'impli.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
      include 'sigbuf.inc'
c      write(ldump,*) '  Beam:', beta, gamm1, brho
c      write(ldump,*) '  Xsig:', xxin, axin, pxin
c      write(ldump,*) '  Ysig:', yyin, ayin, pyin
c      write(ldump,*) '  Zsig:', zzin, azin, pzin
      write(ldump,*) beta, gamm1, brho
      write(ldump,*) xxin, axin, pxin
      write(ldump,*) yyin, ayin, pyin
      write(ldump,*) zzin, azin, pzin
      do 700 nn = nbgn, nend
      write(ldump,617) cname(nn),nn,lintyp(nn),mltyp1(nn),mltyp2(nn),
     & radius(nn), elong(nn), strong(nn)
 617  format(2x,a8,i4,3i3,2f12.7,1pg24.16)
 700  continue
      return
      end
c--------------------------
      subroutine trnspt(ltran)
      return
      end
      subroutine madump(ltran)
      return
      end
