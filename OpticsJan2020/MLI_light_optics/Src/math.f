***********************************************************************
* header                 MATH UTILITIES                               *
*  Routines for handling maps and general mathematical utilities      *
***********************************************************************
c
      subroutine dcdiv(a,b,c,d,e,f)
c   computes the complex division
c     a + ib = (c + id)/(e + if)
c  very slow, but tries to be as accurate as
c  possible by changing the order of the
c  operations, so to avoid under(over)flow
c  problems.
c  Written by F. Neri Feb. 12 1986
c
c      implicit none
      double precision a,b,c,d,e,f
      double precision s,t
      double precision cc,dd,ee,ff
      double precision temp
      integer flip
      flip = 0
      cc = c
      dd = d
      ee = e
      ff = f
      if( dabs(f).ge.dabs(e) ) then
        ee = f
        ff = e
        cc = d
        dd = c
        flip = 1
      endif
      s = 1.d0/ee
      t = 1.d0/(ee+ ff*(ff*s))
      if ( dabs(ff) .ge. dabs(s) ) then
        temp = ff
        ff = s
        s = temp
      endif
      if( dabs(dd) .ge. dabs(s) ) then
        a = t*(cc + s*(dd*ff))
      else if ( dabs(dd) .ge. dabs(ff) ) then
        a = t*(cc + dd*(s*ff))
      else
        a = t*(cc + ff*(s*dd))
      endif
      if ( dabs(cc) .ge. dabs(s)) then
        b = t*(dd - s*(cc*ff))
      else if ( dabs(cc) .ge. dabs(ff)) then
        b = t*(dd - cc*(s*ff))
      else
        b = t*(dd - ff*(s*cc))
      endif
      if (flip.ne.0 ) then
        b = -b
      endif
      return
      end
c
c     ******************************************************************
c
      subroutine dhqr2(nm,n,low,igh,h,wr,wi,z,ierr)
c
c      implicit none
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
     &        igh,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,machep
c     double precision dsqrt,dabs,dsign
c     integer min0
      logical notlas
c     complex z3
      double precision z3r,z3i
c     complex cmplx
c     double precision real,aimag
c
c
c
c     this subroutine is a translation of the algol procedure hqr2,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a real upper hessenberg matrix by the qr method.  the
c     eigenvectors of a real general matrix can also be found
c     if  elmhes  and  eltran  or  orthes  and  ortran  have
c     been used to reduce this general matrix to hessenberg form
c     and to accumulate the similarity transformations.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n,
c
c        h contains the upper hessenberg matrix,
c
c        z contains the transformation matrix produced by  eltran
c          after the reduction by  elmhes, or by  ortran  after the
c          reduction by  orthes, if performed.  if the eigenvectors
c          of the hessenberg matrix are desired, z must contain the
c          identity matrix.
c
c     on output-
c
c        h has been destroyed,
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n,
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if the i-th eigenvalue is real, the i-th column of z
c          contains its eigenvector.  if the i-th eigenvalue is complex
c          with positive imaginary part, the i-th and (i+1)-th
c          columns of z contain the real and imaginary parts of its
c          eigenvector.  the eigenvectors are unnormalized.  if an
c          error exit is made, none of the eigenvectors has been found,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     arithmetic is double precision. complex division
c     is simulated by routin dcdiv.
c
c     fortran routine by b. s. garbow.
c     modified by f. neri.
c
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
      machep = 1.d-15
c     machep = r1mach(4)
c
      ierr = 0
      norm = 0.0
      k = 1
c     ********** store roots isolated by balanc
c                and compute matrix norm **********
      do 50 i = 1, n
c
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0
   50 continue
c
      en = igh
      t = 0.0
c     ********** search for next eigenvalues **********
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
c     ********** look for single small sub-diagonal element
c                for l=en step -1 until low do -- **********
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0) s = norm
         if (dabs(h(l,l-1)) .le. machep * s) go to 100
   80 continue
c     ********** form shift **********
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (its .eq. 30) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     ********** form exceptional shift **********
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
  130 its = its + 1
c     ********** look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- **********
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         if (dabs(h(m,m-1)) * (dabs(q) + dabs(r)) .le. machep * dabs(p)
     &    * (dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0
  160 continue
c     ********** double qr step involving rows l to en and
c                columns m to en **********
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
c     ********** row modification **********
         do 210 j = k, n
            p = h(k,j) + q * h(k+1,j)
            if (.not. notlas) go to 200
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * zz
  200       h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
  210    continue
c
         j = min0(en,k+3)
c     ********** column modification **********
         do 230 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            if (.not. notlas) go to 220
            p = p + zz * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
  220       h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
  230    continue
c     ********** accumulate transformations **********
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            if (.not. notlas) go to 240
            p = p + zz * z(i,k+2)
            z(i,k+2) = z(i,k+2) - p * r
  240       z(i,k+1) = z(i,k+1) - p * q
            z(i,k) = z(i,k) - p
  250    continue
c
  260 continue
c
      go to 70
c     ********** one root found **********
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0
      en = na
      go to 60
c     ********** two roots found **********
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = dsqrt(dabs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0) go to 320
c     ********** real pair **********
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0) wr(en) = x - w / zz
      wi(na) = 0.0
      wi(en) = 0.0
      x = h(en,na)
      s = dabs(x) + dabs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
      p = p / r
      q = q / r
c     ********** row modification **********
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
c     ********** column modification **********
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
c     ********** accumulate transformations **********
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
c
      go to 330
c     ********** complex pair **********
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     ********** all roots found.  backsubstitute to find
c                vectors of upper triangular form **********
  340 if (norm .eq. 0.0) go to 1001
c     ********** for en=n step -1 until 1 do -- **********
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
c     ********** real vector **********
  600    m = en
         h(en,en) = 1.0
         if (na .eq. 0) go to 800
c     ********** for i=en-1 step -1 until 1 do -- **********
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = h(i,en)
            if (m .gt. na) go to 620
c
            do 610 j = m, na
  610       r = r + h(i,j) * h(j,en)
c
  620       if (wi(i) .ge. 0.0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0) go to 640
            t = w
            if (w .eq. 0.0) t = machep * norm
            h(i,en) = -r / t
            go to 700
c     ********** solve real equations **********
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 700
  650       h(i+1,en) = (-s - y * t) / zz
  700    continue
c     ********** end real vector **********
         go to 800
c     ********** complex vector **********
  710    m = na
c     ********** last vector component chosen imaginary so that
c                eigenvector matrix is triangular **********
         if (dabs(h(en,na)) .le. dabs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
c 720    z3 = cmplx(0.0,-h(na,en)) / cmplx(h(na,na)-p,q)
c        h(na,na) = real(z3)
c        h(na,en) = aimag(z3)
  720    call dcdiv(z3r,z3i,0.d0,-h(na,en),h(na,na)-p,q)
         h(na,na) = z3r
         h(na,en) = z3i
  730    h(en,na) = 0.0
         h(en,en) = 1.0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
c     ********** for i=en-2 step -1 until 1 do -- **********
         do 790 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0
            sa = h(i,en)
c
            do 760 j = m, na
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
c
            if (wi(i) .ge. 0.0) go to 770
            zz = w
            r = ra
            s = sa
            go to 790
  770       m = i
            if (wi(i) .ne. 0.0) go to 780
c           z3 = cmplx(-ra,-sa) / cmplx(w,q)
c           h(i,na) = real(z3)
c           h(i,en) = aimag(z3)
            call dcdiv(z3r,z3i,-ra,-sa,w,q)
            h(i,na) = z3r
            h(i,en) = z3i
            go to 790
c     ********** solve complex equations **********
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0 * q
            if (vr .eq. 0.0 .and. vi .eq. 0.0) vr = machep * norm
     &       * (dabs(w) + dabs(q) + dabs(x) + dabs(y) + dabs(zz))
c           z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
c           h(i,na) = real(z3)
c           h(i,en) = aimag(z3)
            call dcdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi)
            h(i,na) = z3r
            h(i,en) = z3i
            if (dabs(x) .le. dabs(zz) + dabs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
c 785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
c           h(i+1,na) = real(z3)
c           h(i+1,en) = aimag(z3)
  785       call dcdiv(z3r,z3i,-r-y*h(i,na),-s-y*h(i,en),zz,q)
            h(i+1,na) = z3r
            h(i+1,en) = z3i
  790    continue
c     ********** end complex vector **********
  800 continue
c     ********** end back substitution.
c                vectors of isolated roots **********
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = i, n
  820    z(i,j) = h(i,j)
c
  840 continue
c     ********** multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- **********
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zz = 0.0
c
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
c
            z(i,j) = zz
  880 continue
c
      go to 1001
c     ********** set error -- no convergence to an
c                eigenvalue after 30 iterations **********
 1000 ierr = en
 1001 return
c     ********** last card of dhqr2 **********
      end
c
c     ******************************************************************
c
      subroutine dorhes(nm,n,low,igh,a,ort)
c      implicit none
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision a(nm,n),ort(igh)
      double precision f,g,h,scale
c     double precision dsqrt,dabs,dsign
c
c     this subroutine is a translation of the algol procedure orthes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     orthogonal similarity transformations.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n,
c
c        a contains the input matrix.
c
c     on output-
c
c        a contains the hessenberg matrix.  information about
c          the orthogonal transformations used in the reduction
c          is stored in the remaining triangle under the
c          hessenberg matrix,
c
c        ort contains further information about the transformations.
c          only elements low through igh are used.
c
c     fortran routine by b. s. garbow
c     modified by filippo neri.
c
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0
         ort(m) = 0.0
         scale = 0.0
c     ********** scale column (algol tol then not needed) **********
         do 90 i = m, igh
   90    scale = scale + dabs(a(i,m-1))
c
         if (scale .eq. 0.0) go to 180
         mp = m + igh
c     ********** for i=igh step -1 until m do -- **********
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
c
         g = -dsign(dsqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
c     ********** form (i-(u*ut)/h) * a **********
         do 130 j = m, n
            f = 0.0
c     ********** for i=igh step -1 until m do -- **********
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
c
            f = f / h
c
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
c
  130    continue
c     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
         do 160 i = 1, igh
            f = 0.0
c     ********** for j=igh step -1 until m do -- **********
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
c
            f = f / h
c
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
c
  160    continue
c
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
c
  200 return
c     ********** last card of dorhes **********
      end
c
c     ******************************************************************
c
      subroutine dorttr(nm,n,low,igh,a,ort,z)
c
c      implicit none
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),ort(igh),z(nm,n)
      double precision g
c
c     this subroutine is a translation of the algol procedure ortrans,
c     num. math. 16, 181-204(1970) by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c
c     this subroutine accumulates the orthogonal similarity
c     transformations used in the reduction of a real general
c     matrix to upper hessenberg form by  dorhes.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n,
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  orthes
c          in its strict lower triangle,
c
c          ort contains further information about the trans-
c          formations used in the reduction by  dorhes.
c          only elements low through igh are used.
c
c     on output-
c
c        z contains the transformation matrix produced in the
c          reduction by  dorhes,
c
c        ort has been altered.
c
c     fortran routine by b. s. garbow.
c     modified by f. neri.
c
c
c     ********** initialize z to identity matrix **********
      do 80 i = 1, n
c
         do 60 j = 1, n
   60    z(i,j) = 0.0
c
         z(i,i) = 1.0
   80 continue
c
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
c     ********** for mp=igh-1 step -1 until low+1 do -- **********
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0) go to 140
         mp1 = mp + 1
c
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
c
         do 130 j = mp, igh
            g = 0.0
c
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
c     ********** divisor below is negative of h formed in orthes.
c                double division avoids possible underflow **********
            g = (g / ort(mp)) / a(mp,mp-1)
c
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
c
  130    continue
c
  140 continue
c
  200 return
c     ********** last card of dorttr **********
      end
c
************************************************************************
c
      subroutine drphse(fm)
c this subroutine readjusts the phases of the eigenvectors making up the
c matrix computed by the subroutine da2 in such a way that fm(1,2)=0 and
c fm(1,1).ge.0, etc.
c Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension fm(6,6)
      dimension t1m(6,6),t2m(6,6)
c
c clear the matrix t1m
      call mclear(t1m)
c compute required phases
      arg1=-fm(1,2)
      arg2=fm(1,1)
      wx=atan2(arg1,arg2)
      arg1=-fm(3,4)
      arg2=fm(3,3)
      wy=atan2(arg1,arg2)
      arg1=-fm(5,6)
      arg2=fm(5,5)
      wt=atan2(arg1,arg2)
c compute rephasing matrix
      cwx=cos(wx)
      swx=sin(wx)
      cwy=cos(wy)
      swy=sin(wy)
      cwt=cos(wt)
      swt=sin(wt)
   10 continue
      t1m(1,1)=cwx
      t1m(1,2)=swx
      t1m(2,1)=-swx
      t1m(2,2)=cwx
      t1m(3,3)=cwy
      t1m(3,4)=swy
      t1m(4,3)=-swy
      t1m(4,4)=cwy
      t1m(5,5)=cwt
      t1m(5,6)=swt
      t1m(6,5)=-swt
      t1m(6,6)=cwt
c rephase fm
      call mmult(fm,t1m,t2m)
c check that t2m(1,1).ge.0 etc.
      iflag=0
      if (t2m(1,1).lt.0.d0) then
      iflag=1
      cwx=-cwx
      swx=-swx
      endif
      if (t2m(3,3).lt.0.d0) then
      iflag=1
      cwy=-cwy
      swy=-swy
      endif
      if (t2m(5,5).lt.0.d0) then
      iflag=1
      cwt=-cwt
      swt=-swt
      endif
      if (iflag.eq.1) goto 10
      call matmat(t2m,fm)
c
      return
      end
c
************************************************************************
c
      subroutine dvsort(revec,aievec)
c this is a subroutine that sorts eigenvectors for a dynamic map.
c it also standardizes the magnitudes and signs of the eigenvctors.
c Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension revec(6,6),aievec(6,6)
      dimension rtv(6,6),aitv(6,6)
      dimension r2v(2),ai2v(2)
      dimension c(6),kpnt(6)
c store vectors in temporary array
      do 10 i=1,5,2
      do 10 j=1,6
      rtv(i,j)=revec(i,j)
      aitv(i,j)=aievec(i,j)
  10  continue
c find vertical components of vectors
      do 20 i=1,5,2
      r2v(1)=rtv(i,3)
      ai2v(1)=aitv(i,3)
      r2v(2)=rtv(i,4)
      ai2v(2)=aitv(i,4)
      sqnrm=r2v(1)**2+ai2v(1)**2+r2v(2)**2+ai2v(2)**2
      c(i)=sqnrm
  20  continue
c find vector with the largest vertical component
      kv=1
      big=c(1)
      if(c(3).gt.big) big=c(3)
      if(c(5).gt.big) big=c(5)
      if(big.eq.c(3)) kv=3
      if(big.eq.c(5)) kv=5
c set pointer
      kpnt(3)=kv
c find  horizontal components of vectors
      do 30 i=1,5,2
      c(i)=0.
      if(i.eq.kv) go to 30
      r2v(1)=rtv(i,1)
      ai2v(1)=aitv(i,1)
      r2v(2)=rtv(i,2)
      ai2v(2)=aitv(i,2)
      sqnrm=r2v(1)**2+ai2v(1)**2+r2v(2)**2+ai2v(2)**2
      c(i)=sqnrm
  30  continue
c find vector with the largest horizontal component
      kh=1
      big=c(1)
      if(c(3).gt.big) big=c(3)
      if(c(5).gt.big) big=c(5)
      if(big.eq.c(3)) kh=3
      if(big.eq.c(5)) kh=5
c set pointer
      kpnt(1)=kh
c find the remaining vector
      do 40 i=1,5,2
      if (i.eq.kv) go to 40
      if (i.eq.kh) go to 40
      kt=i
   40 continue
c set pointer
      kpnt(5)=kt
c reorder vectors
      do 50 i=1,5,2
      k=kpnt(i)
      do 60 j=1,6
      revec(i,j)=rtv(k,j)
      aievec(i,j)=aitv(k,j)
   60 continue
   50 continue
c standardize magnitudes and signs
c the goal is to maximize revec(1,1), revec(3,3), and revec(5,5)
      do 70 i=1,5,2
c see if u and v should be interchanged
      amaxu=abs(revec(i,i))
      amaxv=abs(aievec(i,i))
      if(amaxv.gt.amaxu) then
      do 80 j=1,6
      unew=aievec(i,j)
      vnew=-revec(i,j)
      revec(i,j)=unew
      aievec(i,j)=vnew
   80 continue
      endif
c see if signs of u and v should be changed
      if( revec(i,i) .lt. 0.d0) then
      do 90 j=1,6
      revec(i,j)=-revec(i,j)
      aievec(i,j)=-aievec(i,j)
   90 continue
      endif
   70 continue
      return
      end
c
***********************************************************************
c
      subroutine egphse(fm)
c this subroutine readjusts the phases of the eigenvectors making up the
c matrix computed by the subroutine da2 in such a way that fm(2,1)=0 and
c fm(1,1).ge.0, etc.  Consequently, when fm is transposed, the 1,2 entry
c will be zero, and the 1,1 entry will be .ge. 0, etc.
c Written by Alex Dragt, 5 June 1991
      include 'impli.inc'
      dimension fm(6,6)
      dimension t1m(6,6),t2m(6,6)
c
c clear the matrix t1m
      call mclear(t1m)
c compute required phases
      arg1=fm(2,1)
      arg2=fm(2,2)
      wx=0.d0
      if(arg1 .ne. 0.d0) wx=atan2(arg1,arg2)
      arg1=fm(4,3)
      arg2=fm(4,4)
      wy=0.d0
      if(arg1 .ne. 0.d0) wy=atan2(arg1,arg2)
      arg1=fm(6,5)
      arg2=fm(6,6)
      wt=0.d0
      if(arg1 .ne. 0.d0) wt=atan2(arg1,arg2)
c compute rephasing matrix
      cwx=cos(wx)
      swx=sin(wx)
      cwy=cos(wy)
      swy=sin(wy)
      cwt=cos(wt)
      swt=sin(wt)
   10 continue
      t1m(1,1)=cwx
      t1m(1,2)=swx
      t1m(2,1)=-swx
      t1m(2,2)=cwx
      t1m(3,3)=cwy
      t1m(3,4)=swy
      t1m(4,3)=-swy
      t1m(4,4)=cwy
      t1m(5,5)=cwt
      t1m(5,6)=swt
      t1m(6,5)=-swt
      t1m(6,6)=cwt
c rephase fm
      call mmult(fm,t1m,t2m)
c check that t2m(1,1).ge.0 etc.
      iflag=0
      if (t2m(1,1).lt.0.d0) then
      iflag=1
      cwx=-cwx
      swx=-swx
      endif
      if (t2m(3,3).lt.0.d0) then
      iflag=1
      cwy=-cwy
      swy=-swy
      endif
      if (t2m(5,5).lt.0.d0) then
      iflag=1
      cwt=-cwt
      swt=-swt
      endif
      if (iflag.eq.1) goto 10
      call matmat(t2m,fm)
c
      return
      end

c
************************************************************************
c
      subroutine eig4(fm,reval,aieval,revec,aievec)
c  this routine finds the eigenvalues and eigenvectors
c  of the 4X4 upper left part of fm.
c  the eigenvectors are normalized so that the real and
c  imaginary part of vectors 1 and 3 have +1 antisymmetric
c  product:
c      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1.
c  the eigenvector 2 and 4 have the opposite normalization.
c  written by F. Neri, Feb 26 1986.
c
c      implicit none
      integer i,j,ilo,ihi,nn,mdim,info
      double precision fm(6,6),reval(6),aieval(6)
      double precision revec(6,6),aievec(6,6),pbkt(6)
      double precision vv(6,6),ort(6),aa(6,6)
c clear vector arrays
      do 10 i=1,6
      do 20 j=1,6
      revec(j,i)=0.
      aievec(j,i)=0.
      vv(j,i) = 0.
   20 continue
   10 continue
      ilo = 1
      ihi = 4
      mdim = 6
      nn = 4
c copy matrix to temporary storage ( the matrix aa is destroyed).
      do 100 i=1,6
        do 200 j=1,6
          aa(j,i) = fm(j,i)
  200   continue
  100 continue
c compute eigenvalues and vectors using double precision
c Eispack routines:
      call dorhes(mdim,nn,ilo,ihi,aa,ort)
      call dorttr(mdim,nn,ilo,ihi,aa,ort,vv)
      call dhqr2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
      if ( info .ne. 0 ) then
        write(6,*) '   Something wrong in eig4'
        return
      endif
      vv(5,5) = 1.d0
      vv(6,6) = 1.d0
      call neigv(vv,pbkt)
      reval(5) = 1.d0
      reval(6) = 1.d0
      aieval(5) = 0.d0
      aieval(6) = 0.d0
      do 300 j=1,4
        revec(1,j) = vv(j,1)
        revec(2,j) = vv(j,1)
        revec(3,j) = vv(j,3)
        revec(4,j) = vv(j,3)
        aievec(1,j) = vv(j,2)
        aievec(2,j) = -vv(j,2)
        aievec(3,j) = vv(j,4)
        aievec(4,j) = -vv(j,4)
  300 continue
c  if poisson bracket is negative,change sign of imaginary part of
c  eigenvalue
      do 400 i=2,4,2
        if( pbkt(i) .lt. 0.d0 ) then
          aieval(i-1) = -aieval(i-1)
          aieval(i)   = -aieval(i)
        endif
  400 continue
      revec(5,5) = 1.d0
      revec(6,6) = 1.d0
c  if eigenvalues are off the unit circle, print warning message:
      do 600 i=1,4
        if(dabs(reval(i)**2+aieval(i)**2 - 1.d0).gt.1.d-10) then
           write(6,*) ' Eig4: Eigenvalues off the unit circle!'
           write(6,*) ' delta=',dabs(reval(i)**2+aieval(i)**2 - 1.d0)
           return
        endif
  600 continue
      return
      end
c
********************************************************************
c
      subroutine eig6(fm,reval,aieval,revec,aievec)
c  this routine finds the eigenvalues and eigenvectors
c  of the full matrix fm.
c  the eigenvectors are normalized so that the real and
c  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
c  product:
c      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
c      revec5 J aivec5 = 1.
c  the eigenvectors 2 ,4, and 6 have the opposite normalization.
c  written by F. Neri, Feb 26 1986.
c
c      implicit none
      integer nn
      integer nnr,ilo,ihi,mdim,info
      double precision reval(6),aieval(6),revec(6,6),aievec(6,6)
      double precision fm(6,6),aa(6,6)
      integer i,i1,j
      double precision ort(6),vv(6,6)
      double precision pbkt(6)
c  copy matrix to temporary storage (the matrix aa is destroyed)
      do 600 i=1,6
        do 600 i1=1,6
          aa(i1,i) = fm(i1,i)
  600 continue
      ilo = 1
      ihi = 6
      mdim = 6
      nn = 6
c  compute eigenvalues and eigenvectors using double
c  precision Eispack routines:
      call dorhes(mdim,nn,ilo,ihi,aa,ort)
      call dorttr(mdim,nn,ilo,ihi,aa,ort,vv)
      call dhqr2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
      if ( info .ne. 0 ) then
        write(6,*) '  ERROR IN EIG6'
        return
      endif
      call neigv(vv,pbkt)
      do 700 j=1,6
        revec(1,j) = vv(j,1)
        revec(2,j) = vv(j,1)
        revec(3,j) = vv(j,3)
        revec(4,j) = vv(j,3)
        revec(5,j) = vv(j,5)
        revec(6,j) = vv(j,5)
        aievec(1,j) = vv(j,2)
        aievec(2,j) = -vv(j,2)
        aievec(3,j) = vv(j,4)
        aievec(4,j) = -vv(j,4)
        aievec(5,j) = vv(j,6)
        aievec(6,j) = -vv(j,6)
  700 continue
c  if poisson bracket is negative, change sign of imginary part of
c  eigenvalue
      do 800 i=2,6,2
        if( pbkt(i) .lt. 0.d0) then
          aieval(i-1) = -aieval(i-1)
          aieval(i  ) = -aieval(i  )
        endif
  800 continue
c  if eigenvalues are off unit circle, print warning message:
      do 900 i=1,6
        if(dabs(reval(i)**2+aieval(i)**2 -1.d0).gt.1.d-10) then
          write(6,*) ' EIG6: Eigenvalues off the unit circle!'
          return
        endif
  900 continue
      return
      end
c
************************************************************************
c 
      subroutine eigemt(idim,ha)
c This subroutine computes eigen emittances, etc.
c The array ha contains the incoming moments.  
c As described below, rusults are put in buffers 1 through 5.
c This subroutine eventually needs to be improved because at present it
c does not deal properly with the case of degenerate eigenvalues.
c Written by Alex Dragt, 22 May 1991.
c Modified by Alex Dragt, 25 May 1998 to deal with
c phase spaces of various dimensions.
c Again modified 10/14/98 by AJD to change contents of buffers 1 and 2.
c
      use lieaparam
      include 'impli.inc'
      include 'buffer.inc'
c
      dimension ha(monoms)
c
c local arrays
      dimension t1m(6,6), t2m(6,6)
c
c Case of 2-dimensional phase space
c
      if (idim .eq. 2) then
c compute "diagonalizing" matrix and put it in the matrix part of buffer 1.
      z11=ha(7)
      if(z11 .le. 0.d0) then
      write(6,*) 'error: <xx> is negative or zero'
      return
      endif
      z12=ha(8)
      z22=ha(13)
      if(z22 .le. 0.d0) then
      write(6,*) 'error: <pxpx> is negative or zero'
      return
      endif
c remove possibly offensive moments in 2-dimensional case.
      do 2 i=7,27
    2 ha(i)=0.d0
      ha(7)=z11
      ha(8)=z12
      ha(13)=z22
      emit2=z11*z22-z12**2
      if(emit2 .le. 0.d0) then
      write(6,*) 'error: x emittance is zero or imaginary'
      return
      endif
      emit=dsqrt(emit2)
      beta=z11/emit
      rbeta=dsqrt(beta)
      alpha=-z12/emit
      call mclear(buf1m)
      buf1m(1,1)=1.d0/rbeta
      buf1m(1,2)=alpha/rbeta
      buf1m(2,1)=0.d0
      buf1m(2,2)=rbeta
      buf1m(3,3)=1.d0
      buf1m(4,4)=1.d0
      buf1m(5,5)=1.d0
      buf1m(6,6)=1.d0
      go to 100
      endif
c
c Case of 4- and 6-dimensional phase space
c
      if ((idim .eq. 4) .or. (idim .eq. 6)) then  
c compute "diagonalizing" matrix and put it in the matrix part of buffer 1.
c
c Remove possibly offensive moments in 4-dimensional case.
c
      if(idim .eq. 4) then
      ha(11)=0.d0
      ha(12)=0.d0
      ha(16)=0.d0
      ha(17)=0.d0
      ha(20)=0.d0
      ha(21)=0.d0
      ha(23)=0.d0
      ha(24)=0.d0
      ha(25)=0.d0
      ha(26)=0.d0
      ha(27)=0.d0
      endif
c
c let Z denote the matrix with entries Zij=<zizj>.
c compute the polynomial -(1/2)*(z,Zz) and temporarily
c put result in buf2a
      do 5 i=7,27
    5 buf2a(i)=-ha(i)
      buf2a(7)=buf2a(7)/2.d0
      buf2a(13)=buf2a(13)/2.d0
      buf2a(18)=buf2a(18)/2.d0
      buf2a(22)=buf2a(22)/2.d0
      buf2a(25)=buf2a(25)/2.d0
      buf2a(27)=buf2a(27)/2.d0
c
c matify buf2a.  this should produce the matrix JZ.
      call matify(buf1m,buf2a)
c compute norm of buf1m=JZ.
      call mnorm(buf1m,ans)
c compute taylor series result buf2m=exp(scale*buf1m)
      scale=.1d0/ans
      call smmult(scale,buf1m,buf1m)
      call exptay(buf1m,buf2m)
c compute normal form transformation
      if (idim .eq. 4) call sa2(buf2m,buf2a,buf1m)
      if (idim .eq. 6) call da2(buf2m,buf2a,buf1m)
c
      endif    !cryne 8/9/2002
  100 continue
c
c rephase and transpose result
      call egphse(buf1m)
      call mtran(buf1m)
cryne 8/9/2002 moved endif to before "100 continue"      endif
c
c act on moments with "diagonalizing" matrix and put result in buf2a
c letting the map (matrix) act on moments amounts to computing buf2a = D*ha
c
c clear arrays (the array buf2a has already been cleared by da2)
c temporarily use buf3a, buf4a, and buf5a.
      do 10 i=1,monoms
      buf3a(i)=0.d0
   10 buf4a(i)=0.d0
c
c compute "diagonalized" moments only through 2'nd moments
      imax=27
c
c perform calculation
      do 20 i=7,imax
      buf3a(i)=1.d0
      call fxform(buf4a,buf1m,buf3a,buf5a)
      buf3a(i)=0.d0
      do 30 j=1,imax
   30 buf2a(i)=buf2a(i)+buf5a(j)*ha(j)
   20 continue
c
c put result in buf1a
      do 40 i=7,monoms
   40 buf1a(i)=buf2a(i)
c
c at this stage buf1m contains the "diagonalizing" matrix.
c put this matrix in buf2m, and put the inverse of the
c diagonalizing matrix in buf1m.
      call matmat(buf1m,buf2m)
      call inv(buf3a,buf1m)
c
c put original moments in buf2a
      do 50 i=7,monoms
   50 buf2a(i)=ha(i)
c
c compute results for buffers 3 through 5
c need to compute buf1m*diagj*(buf1m transpose) for j=X,Y,Tau,
c and then multiply the result by eigemitj.
c
      call matmat(buf1m,t1m)
      call mtran(t1m)
c
      call mclear(t2m)
      t2m(1,1)=1.d0
      t2m(2,2)=1.d0
      call mmult(buf1m,t2m,buf3m)
      call mmult(buf3m,t1m,buf3m)
      scale=buf1a(7)
      call smtof(scale,buf3m,buf3a)       
c
      call mclear(t2m)
      t2m(3,3)=1.d0
      t2m(4,4)=1.d0
      call mmult(buf1m,t2m,buf4m)
      call mmult(buf4m,t1m,buf4m)
      scale=buf1a(18)
      call smtof(scale,buf4m,buf4a)       
c
      call mclear(t2m)
      t2m(5,5)=1.d0
      t2m(6,6)=1.d0
      call mmult(buf1m,t2m,buf5m)
      call mmult(buf5m,t1m,buf5m)
      scale=buf1a(25)
      call smtof(scale,buf5m,buf5a)       
c
      return
      end
c
************************************************************************
c
      subroutine exptay(em,fm)
c this subroutine computes fm=exp(em) using a finite taylor series
c Written by Alex Dragt, Fall 1986, based on work of Liam Healy
      include 'impli.inc'
      include 'id.inc'
      dimension em(6,6),fm(6,6),tm1(6,6),tm2(6,6)
c
c initialize fm and tm1 to be the identity matrix
      call matmat(ident,fm)
      call matmat(ident,tm1)
c
c begin calculation
      kmax=10
      do 10 k=1,kmax
      rk=1.d0/float(k)
      call mmult(em,tm1,tm2)
      call smmult(rk,tm2,tm1)
      call madd(fm,tm1,fm)
   10 continue
c
      return
      end
c
***********************************************************************
c
      subroutine leshs(soln, dim,matrix,rhs,augm,det)
c  Linear Equation Solver HandShaking routine.
c  Interfaces Etienne's need for linear equation solutions with my
c  'solver'.
c  Written by Liam Healy, May 9, 1985.
c
c----Variables----
c  dim = size of linear space
      integer dim
c  matrix = matrix supplied (dim x dim)
c  rhs = right hand side of equation, supplied (dim)
c  augm = augmented martrix (dim x dim+1), set up with arbitrary
c         contents by calling routine
c  det = determinant of original matrix
      double precision matrix(dim,dim),rhs(dim),soln(dim),det
      double precision augm(dim,dim+1)
c
c----Routine----
      do 120 i=1,dim
      do 100 j=1,dim
  100   augm(i,j)=matrix(i,j)
  120   augm(i,dim+1)=rhs(i)
      call solver(augm,dim,det)
      do 140 i=1,dim
  140   soln(i)=augm(i,dim+1)
      return
      end
c
      subroutine madd(a,b,c)
c This is a subroutine for adding two matrices.
c Written by Alex Dragt on 11 Nov 1985.
      double precision a(6,6),b(6,6),c(6,6)
      do 10 j=1,6
      do 10 i=1,6
   10 c(i,j)=a(i,j)+b(i,j)
      return
      end
c
***********************************************************************
c
      subroutine matmat(a,b)
c This is a subroutine for copying a matrix.
c Written by Alex Dragt on 11 Nov 1985.
      double precision a(6,6),b(6,6)
c---Routine---
      do 10 j=1,6
      do 10 i=1,6
   10 b(i,j)=a(i,j)
      return
      end
c
***********************************************************************
c
      subroutine mmult(a,b,c)
c This is a subroutine for matrix multiplication.
c Written by Alex Dragt on Friday, 13 Sept 1985.
      double precision a(6,6),b(6,6),c(6,6),ct(6,6),sum
c----Routine----
      do 10 k=1,6
        do 20 i=1,6
          sum = 0.d0
          do 30 j=1,6
            sum = sum + a(i,j)*b(j,k)
   30     continue
          ct(i,k) = sum
   20   continue
   10 continue
      call matmat(ct,c)
c
      return
      end
c
***********************************************************************
c
      subroutine mnorm(fm,res)
c  Computes the norm (maximum column sum norm) for the matrix fm.
c  Reference: L. Collatz, Functional Analysis & Numerical Mathematics,
c             p.177
c  Written by Alex Dragt, Fall 1986, based on work of Liam Healy
      include 'impli.inc'
      dimension fm(6,6),sum(6)
c
c initialize variables
      res=0.d0
      do 100 j=1,6
  100   sum(j)=0.d0
c
c perform calculation
      do 120 j=1,6
        do 110 i=1,6
  110  	  sum(j)=sum(j)+abs(fm(i,j))
  120   res=max(res,sum(j))
c
      return
      end
c
***********************************************************************
c
       subroutine neigv(m,pbkt)
c  this subroutine normalizes the eigenvectors of
c  a stable ( all roots on the unit circle ) symplectic
c  matrix. on entry m has on the odd numbered colums
c  the real part of the eigenvectors, on the even ones
c  the imaginary parts. only one real (imaginary ) vector
c  is included of each complex conjugate pair.
c  on exit the vectors are rescaled so that the poisson
c  brackets of the real x imaginary parts are equal to 1.
c  also, the real and imaginary parts of each eigenvector
c  are made orthogonal.
c  the resulting matrix is symplectic.  when used in sa2 or da2,
c  it tranforms
c  the original matrix to block diagonal form, with the
c  blocks being 2 dimensional rotations.
c  WARNING: this only works if the eigenvalues are on
c  the unit circle.
c  written by F. Neri Feb 10 1986.
c       implicit none
       integer n
       parameter ( n = 3 )
       double precision m(2*n,2*n),pbkt(2*n)
       double precision pb,s
       double precision usq,vsq,udotv,a,b,phi(6),sn,cn,unew,vnew
       integer k,iq,ip,i,j
c  rescale u and v
       do 10 k=1,5,2
         pb = 0.d0
         do 20 ip = 2,2*n,2
           iq = ip-1
           pb = pb + m(iq,k)*m(ip,k+1) - m(ip,k)*m(iq,k+1)
   20    continue
         s = dsqrt(dabs(pb))
c         write(6,*) ' PB=',pb
         pbkt(k) = pb
         pbkt(k+1)=pb
         do 30 i=1,2*n
         m(i,k) = m(i,k)/s
         m(i,k+1) = m(i,k+1)*(s/pb)
   30    continue
   10  continue
c  orthogonalize u and v
c  compute required phase
       do 40 k=1,5,2
       usq=0.d0
       vsq=0.d0
       udotv=0.d0
       do 50 j=1,6
       usq=usq+m(j,k)**2
       vsq=vsq+m(j,k+1)**2
       udotv=udotv+m(j,k)*m(j,k+1)
   50  continue
       phi(k)=0.d0
       a=udotv
       if(a.eq.0.d0) goto 40
       b=(usq-vsq)/2.d0
       phi(k)=(1/2.d0)*atan2(-a,b)
   40  continue
c  transform u and v
       do 60 k=1,5,2
       sn=sin(phi(k))
       cn=cos(phi(k))
       do 70 j=1,6
       unew=cn*m(j,k)-sn*m(j,k+1)
       vnew=sn*m(j,k)+cn*m(j,k+1)
       m(j,k)=unew
       m(j,k+1)=vnew
   70  continue
   60  continue
       return
       end
c end of file
      subroutine pmadd(f,n,coeff,h)
      implicit double precision (a-h,o-z)
      include 'len.inc'
      dimension f(923),h(923)
      if(n.eq.1)then
        istart=1
      else
        istart=len(n-1)+1
      endif
c
      if(coeff.eq.1.d0) goto 20
cryne do 10 i=len(n-1)+1,len(n)
      do 10 i=istart,len(n)
        h(i) = h(i) + f(i)*coeff
 10   continue
      return
 20   continue
cryne do 30 i = len(n-1)+1, len(n)
      do 30 i = istart, len(n)
        h(i) = h(i) + f(i)
 30   continue
      return
      end
c
      subroutine product(a,na,b,nb,c)
      use lieaparam, only : monoms
      implicit double precision(a-h,o-z)
      include 'len.inc'
      include 'expon.inc'
      include 'vblist.inc'
cryne 7/23/2002      common/expon/expon
cryne 7/23/2002      common/vblist/vblist
cryne 7/23/2002       common /len/ len(16)
cryne 7/23/2002       integer expon(6,0:923),vblist(6,0:923)
       dimension a(923),b(923),c(923),l(6)
      if(na.eq.1) then
        ia1 = 1
      else
        ia1 = len(na-1)+1
      endif
      if(nb.eq.1) then
        ib1 = 1
      else
        ib1 = len(nb-1)+1
      endif
       do 200 ia=ia1,len(na)
           if(a(ia).eq.0.d0) goto 200
           do 20 ib = ib1,len(nb)
               if(b(ib).eq.0.d0) goto 20
               do 2 m=1,6
                   l(m) = expon(m,ia) +  expon(m,ib)
   2            continue
                n = ndex(l)
                c(n) = c(n) + a(ia)*b(ib)
  20        continue
 200   continue
       return
       end
c
********************************************************************
c
       subroutine pttodp(f,ft)
c
c   This subroutine converts a power series in the variable Ptau
c   (which is the negative of the normalized energy deviation) to a power
c   series in the normalized momentum deviation [deltap=(delta p)/p]. 
c   This routine works through third order.
c   Reference: A. Dragt and M. Venturini, Relation between Expansions
c   in Energy and Momentum Deviation Variables, U.MD. Physics Dept.
c   Technical Report (1995).
c   Written in Sept. 1995 by A. Dragt and M. Venturini
c
c   f is an aray containing the coefficients of the power
c   series in terms of the variable Ptau.
c   ft is a transformed aray that contains the coefficients
c   of the transformed power series in terms of the variable deltap.
c   The transformation uses the matrix a(i,j).
c
       use beamdata
       include 'impli.inc'
       dimension f(*), ft(*)
       dimension a(3,3)
c
c Set up transforming coefficients
c
       beta2=beta*beta
       beta3=beta2*beta
       beta4=beta2*beta2
       beta5=beta2*beta3
c
       a(1,1)= -beta
       a(2,1)=(-beta + beta3)/(2.d0)
       a(2,2)=  beta2
       a(3,1)=( beta3 - beta5)/(2.d0)
       a(3,2)=  beta2 - beta4
       a(3,3)= -beta3
c
c Make transformation
c 
       ft(1)=a(1,1)*f(1)
       ft(2)=a(2,1)*f(1) + a(2,2)*f(2) 
       ft(3)=a(3,1)*f(1) + a(3,2)*f(2) + a(3,3)*f(3)
c
       return
       end
c
*******************************************************************
c
      subroutine seig6(fm,reval,revec)
c  this routine finds the eigenvalues and eigenvectors
c  of a real symmetric 6x6 matrix
c  the eigenvectors are normalized so that related pairs (those
c  with eigenvalues eval and 1/eval) have antisymmetric product:
c      revec1 J revec2 = 1 ; revec3 J revec4 = 1 ;
c      revec5 J revec6 = 1.
c  further details about normalization and ordering
c  conventions may be gleaned from examination of the code
c  written by A. Dragt 20 March 1987
c
      include 'impli.inc'
      dimension fm(6,6),reval(6),revec(6,6)
      dimension aieval(6),ort(6),pbkt(6)
      dimension tm(6,6),vv(6,6)
c
c  begin computation
c
c  copy matrix to temporary storage (the matrix tm is destroyed)
      call matmat(fm,tm)
c
c  set up control indices
      ilo = 1
      ihi = 6
      mdim = 6
      nn = 6
c
c  compute eigenvalues and eigenvectors using double
c  precision Eispack routines:
      call dorhes(mdim,nn,ilo,ihi,tm,ort)
      call dorttr(mdim,nn,ilo,ihi,tm,ort,vv)
      call dhqr2(mdim,nn,ilo,ihi,tm,reval,aieval,vv,info)
      if ( info .ne. 0 ) then
        write(6,*) '  Error in seig6 from Eispack routines '
        return
      endif
c
c  sort eigenvalues and eigenvectors into related pairs
c
c  normalize eigenvectors to have unit Euclidean norm
      do 10 i=1,6
      vsq=0.d0
      do 20 j=1,6
      vsq=vsq+vv(j,i)**2
   20 continue
      rv=1.d0/sqrt(vsq)
      do 30 k=1,6
      revec(k,i)=rv*vv(k,i)
   30 continue
   10 continue
c  find largest eigenvalue
      big=0.d0
      do 40 i=1,6
      if (reval(i).gt.big) then
      big=reval(i)
      imax1=i
      endif
   40 continue
c  find its reciprocal pair based on symplectic 2-form J
      big=0.d0
      do 50 i=1,6
      if (i.eq.imax1) goto 50
      call s2f(revec,imax1,i,val)
      aval=abs(val)
      if (aval.gt.big) then
      big=aval
      imax1r=i
      endif
   50 continue
c  find second largest eigenvalue
      big=0.d0
      do 60 i=1,6
      if (i.eq.imax1 .or. i.eq.imax1r) goto 60
      if (reval(i).gt.big) then
      big=reval(i)
      imax2=i
      endif
   60 continue
c  find its reciprocal pair based on symplectic 2-form J
      big=0.d0
      do 70 i=1,6
      if (i.eq.imax1 .or. i.eq.imax1r .or. i.eq.imax2) goto 70
      call s2f(revec,imax2,i,val)
      aval=abs(val)
      if (aval.gt.big) then
      big=aval
      imax2r=i
      endif
   70 continue
c  find third largest eigenvalue
      big=0.d0
      do 80 i=1,6
      if (i.eq.imax1 .or. i.eq.imax1r) goto 80
      if (i.eq.imax2 .or. i.eq.imax2r) goto 80
      if (reval(i).gt.big) then
      big=reval(i)
      imax3=i
      endif
   80 continue
c  find its reciprocal pair by elimination
      do 90 i=1,6
      if (i.eq.imax1 .or. i.eq.imax1r) goto 90
      if (i.eq.imax2 .or. i.eq.imax2r) goto 90
      if (i.eq.imax3) goto 90
      imax3r=i
   90 continue
c
c  renormalize eigenvectors
      call s2f(revec,imax1,imax1r,prd)
      aprod=abs(prd)
      sign=1.d0
      if (prd.lt.0.d0) sign=-1.d0
      rfact=1.d0/sqrt(aprod)
      do 100 i=1,6
      vv(i,imax1)=rfact*revec(i,imax1)
      vv(i,imax1r)=sign*rfact*revec(i,imax1r)
  100 continue
      call s2f(revec,imax2,imax2r,prd)
      aprod=abs(prd)
      sign=1.d0
      if (prd.lt.0.d0) sign=-1.d0
      rfact=1.d0/sqrt(aprod)
      do 110 i=1,6
      vv(i,imax2)=rfact*revec(i,imax2)
      vv(i,imax2r)=sign*rfact*revec(i,imax2r)
  110 continue
      call s2f(revec,imax3,imax3r,prd)
      aprod=abs(prd)
      sign=1.d0
      if (prd.lt.0.d0) sign=-1.d0
      rfact=1.d0/sqrt(aprod)
      do 120 i=1,6
      vv(i,imax3)=rfact*revec(i,imax3)
      vv(i,imax3r)=sign*rfact*revec(i,imax3r)
  120 continue
c
c  rearrange eigenvectors
      do 130 i=1,6
      revec(i,1)=vv(i,imax1)
      revec(i,2)=vv(i,imax1r)
      revec(i,3)=vv(i,imax2)
      revec(i,4)=vv(i,imax2r)
      revec(i,5)=vv(i,imax3)
      revec(i,6)=vv(i,imax3r)
  130 continue
c
c  rearrange eigenvalues
      do 140 i=1,6
      aieval(i)=reval(i)
  140 continue
      reval(1)=aieval(imax1)
      reval(2)=aieval(imax1r)
      reval(3)=aieval(imax2)
      reval(4)=aieval(imax2r)
      reval(5)=aieval(imax3)
      reval(6)=aieval(imax3r)
c
      return
      end
c
***********************************************************************
c
      subroutine smmult(s,a,b)
c This is a subroutine for multiplying a matrix by a scalar.
c Written by Alex Dragt on 11 Nov 1985.
      double precision s,a(6,6),b(6,6)
c---Routine---
      do 10 j=1,6
      do 10 i=1,6
   10 b(i,j)=s*a(i,j)
      return
      end
c
***********************************************************************
c
      subroutine solver(augmat,dim,det)
c  Solves the linear equation
c       m*a = b
c  where m is  n by n
c  On entry :   augmat = m augmented by b
c               dim = n
c  On return:   matrix = identity augmented by a
c               det = determinant of m
c  Coded from a routine originally written for the HP41C calculator
c  (Dearing, p.46).  Written by Liam Healy, February 14, 1985.
c
c----Variables----
c  dim= dimension of matrix
      integer dim
c  matrix= input matrix m, vector = input and output vector
      double precision augmat(dim,dim+1)
c  det = determinant returned
      double precision det
c  row,col,r,c,roff,rs = row and column numbers, row and column indices,
c               row offset, row number for finding max
      integer row,col,r,c,roff,rs
c  nrow,ncol = total number of rows and columns
      integer nrow,ncol
c  me, mer, h = matrix element and its row number, held value of mat elt
c  const = constant used in multiplication
      double precision me,h,const
      integer mer
c
c----Routine----
      det=1.
      nrow=dim
      ncol=dim+1
      col=0
      do 100 row=1,nrow
 300    col=col+1
        if (col.le.ncol) then
c find max of abs of mat elts in this col, and its row number
          me=0.
          mer=0
          do 120 rs=row,nrow
            if (abs(augmat(rs,col)).ge.abs(me)) then
              me=augmat(rs,col)
              mer=rs
            endif
  120     continue
          det=det*me
          if (me.eq.0.) goto 300
          do 140 c=1,ncol
            augmat(mer,c)=augmat(mer,c)/me
  140     continue
          r=0
          if (mer.ne.row) then
c               swap the rows
            do 160 c=1,ncol
              h=augmat(mer,c)
              augmat(mer,c)=augmat(row,c)
              augmat(row,c)=h
  160       continue
            det=-det
          endif
  320     r=r+1
          if (r.le.nrow.and.row.lt.nrow) then
            if (r.eq.row) goto 320
c               multiply row row by const & subtract from row r
            const=augmat(r,col)
            do 180 c=1,ncol
              augmat(r,c)=augmat(r,c)-augmat(row,c)*const
  180       continue
            goto 320
          endif
        endif
  100 continue
c
c  Matrix is now in upper triangular form.
c  To solve equation, must get it to the identity.
      do 200 roff=nrow,1,-1
        do 200 row=roff-1,1,-1
          const=augmat(row,roff)
          do 200 c=row,ncol
            augmat(row,c)=augmat(row,c)-const*augmat(roff,c)
  200 continue
      return
      end
c
************************************************************************
c
      subroutine srphse(fm)
c this subroutine readjusts the phases of the eigenvectors making up the
c matrix computed by the subroutine sa2 in such a way that fm(1,2)=0 and
c fm(1,1).ge.0, etc.
c Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension fm(6,6)
      dimension t1m(6,6),t2m(6,6)
c
c clear the matrix t1m
      call mclear(t1m)
c compute required phases
      arg1=-fm(1,2)
      arg2=fm(1,1)
      wx=atan2(arg1,arg2)
      arg1=-fm(3,4)
      arg2=fm(3,3)
      wy=atan2(arg1,arg2)
c compute rephasing matrix
      cwx=cos(wx)
      swx=sin(wx)
      cwy=cos(wy)
      swy=sin(wy)
   10 continue
      t1m(1,1)=cwx
      t1m(1,2)=swx
      t1m(2,1)=-swx
      t1m(2,2)=cwx
      t1m(3,3)=cwy
      t1m(3,4)=swy
      t1m(4,3)=-swy
      t1m(4,4)=cwy
      t1m(5,5)=1.d0
      t1m(6,6)=1.d0
c rephase fm
      call mmult(fm,t1m,t2m)
c check that t2m(1,1).ge.0 etc.
      iflag=0
      if (t2m(1,1).lt.0.d0) then
      iflag=1
      cwx=-cwx
      swx=-swx
      endif
      if (t2m(3,3).lt.0.d0) then
      iflag=1
      cwy=-cwy
      swy=-swy
      endif
      if (iflag.eq.1) goto 10
      call matmat(t2m,fm)
c
      return
      end
c
***********************************************************************
c
      subroutine svsort(revec,aievec)
c this is a subroutine that sorts eigenvectors for a static map.
c it also standardizes the magnitudes and signs of the eigenvectors.
c Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension revec(6,6),aievec(6,6)
      dimension rtv(6,6),aitv(6,6)
      dimension r2v(2),ai2v(2)
      dimension c(6),kpnt(6)
c store vectors in temporary array
      do 10 i=1,3,2
      do 10 j=1,6
      rtv(i,j)=revec(i,j)
      aitv(i,j)=aievec(i,j)
  10  continue
c find vertical components of vectors
      do 20 i=1,3,2
      r2v(1)=rtv(i,3)
      ai2v(1)=aitv(i,3)
      r2v(2)=rtv(i,4)
      ai2v(2)=aitv(i,4)
      sqnrm=r2v(1)**2+ai2v(1)**2+r2v(2)**2+ai2v(2)**2
      c(i)=sqnrm
  20  continue
c find vector with the largest vertical component
      kv=1
      big=c(1)
      if(c(3).gt.big) kv=3
c set pointer
      kpnt(3)=kv
c find the remaining vector
      do 30 i=1,3,2
      if (i.eq.kv) go to 30
      kh=i
   30 continue
c set pointer
      kpnt(1)=kh
c reorder vectors
      do 40 i=1,3,2
      k=kpnt(i)
      do 50 j=1,6
      revec(i,j)=rtv(k,j)
      aievec(i,j)=aitv(k,j)
   50 continue
   40 continue
c standardize magnitudes and signs
c the goal is to maximize revec(1,1) and revec(3,3)
      do 60 i=1,3,2
c see if u and v should be interchanged
      amaxu=abs(revec(i,i))
      amaxv=abs(aievec(i,i))
      if(amaxv.gt.amaxu) then
      do 70 j=1,6
      unew=aievec(i,j)
      vnew=-revec(i,j)
      revec(i,j)=unew
      aievec(i,j)=vnew
   70 continue
      endif
c see if signs of u and v should be changed
      if( revec(i,i) .lt. 0.d0) then
      do 80 j=1,6
      revec(i,j)=-revec(i,j)
      aievec(i,j)=-aievec(i,j)
   80 continue
      endif
   60 continue
      return
      end
c
************************************************************************
c
      subroutine sympl1(fm)
      include 'impli.inc'
      dimension fm(6,6)
      return
      end
c
************************************************************************
c
      subroutine sympl2(fm)
      include 'impli.inc'
      dimension fm(6,6)
      return
      end
c
***********************************************************
c
c    SYMPL3
c
c**********************************************************
c
c  Written by F. Neri  Feb 7 1986
c
      subroutine sympl3(m)
c      implicit none
      integer n
      parameter ( n = 3 )
      double precision m(2*n,2*n)
c
c   On return ,the matrix m(*,*), supposed to be almost
c   symplectic on entry is made exactly symplectic by
c   using a non iterative, constructive method.
c
      double precision qq,pq,qp,pp
      integer kp,kq,lp,lq,jp,jq,i
c
      do 100 kp=2,2*n,2
        kq = kp-1
        do 200 lp=2,kp-2,2
          lq = lp-1
          qq = 0.d0
          pq = 0.d0
          qp = 0.d0
          pp = 0.d0
          do 300 jp=2,2*n,2
            jq = jp-1
            qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq)
            pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq)
            qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq)
            pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq)
  300     continue
c         write(6,*) qq,pq,qp,pp
          do 400 i=1,2*n
            m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i)
            m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i)
  400     continue
  200   continue
        qp = 0.d0
        do 500 jp=2,2*n,2
          jq = jp-1
          qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq)
  500   continue
c       write(6,*) qp
        do 600 i=1,2*n
          m(kp,i) = m(kp,i)/qp
  600   continue
c
c  Maybe the following is a better idea ( uses sqrt and is slower )
c       sign = 1.d0
c       if ( qp.lt.0.d0 ) sign = -1.d0
c  OR, BETTER:
c       if ( qp.le.0.d0 ) then complain
c       qp = dabs(qp)
c       qp = dsqrt(qp)
c       do 600 i=1,2*n
c         m(kq,i) = m(kq,i)/qp
c         m(kp,i) = sign*m(kp,i)/qp
c 600   continue
  100 continue
      return
      end
