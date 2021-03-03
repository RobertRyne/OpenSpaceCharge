c
c***************************************
c
c   OPTI is an alphabetized collection of the routines in
c   NLS, QSO, and SHARED, by Kurt Overley, LANL 89,
c   as modified by Tom Mottershead, March 91
c--------------------------------------------------------
      subroutine condn( maxf,nx,fjac,info)
c
c subroutine condn returns the condition number of the jacobian.
c----------------------------------------------------------------------
      include 'impli.inc'
      parameter (epsmach  = 1.0d-16)
      dimension fjac(maxf,nx)
c
      dmax = dabs(fjac(1,1))
      dmin = dmax
      do 10 i = 2,nx
         diag = dabs(fjac(i,i))
         if (diag .gt. dmax) then
            dmax = diag
         else if (diag .lt. dmin ) then
            dmin = diag
         endif
   10 continue
      rcond = dmin / dmax
      if (rcond .lt. epsmach) info = 5
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function enorm(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function enorm
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     &                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      save one,zero,rdwarf,rgiant  !cryne 7/23/2002
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
c              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
c              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
c           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
c     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     &         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     &         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine findhi( dmuhi, ld, n, h, hmu,
     &                   g, step, region, info)
c
c  subroutine findhi finds the upper value of the levenberg-marquardt
c  parameter that brackets the root of the nonlinear equation,
c  phi = region - stepnm.
c--------------------------------------------------------
      include 'impli.inc'
      dimension h(ld,n), hmu(ld,n)
      dimension g(n), step(n)
c
   10 continue
      do 20 i = 1,n
         call vecopy( n, h(1,i), hmu(1,i))
         hmu(i,i) = h(i,i) + dmuhi
   20 continue
      call newton( ld, n, hmu, g, step, info)
      if ( info .gt. 0) return
      stepnm = enorm(n, step)
      if (stepnm .gt. region) then
         dmuhi = 2.0d0 * dmuhi
         go to 10
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine geta (maxv, maxq, nx, m, a, xu)
c
c  subroutine geta forms the linear system in the center of mass
c  coordinates whose solution yields the quadratic coefficients.
c---------------------------------------------------------------
      include 'impli.inc'
      dimension a(maxq,m), xu(maxv, m)
c
      do 10 k = 1, m
         a(k,1) = 1.0
         do 20 i = 1, nx
            a(k,i+1) = xu(i,k)
            if ( (i+1+nx) .gt. m) go to 20
            a(k,i+1+nx) = xu(i,k)**2
   20    continue
   25    indexa = 2 * nx + 1
         do 30 i = 1, nx-1
            do 40 j = i+1,nx
               indexa = indexa + 1
               if (indexa .gt. m) go to 10
               a(k,indexa) = xu(i,k) * xu(j,k)
   40       continue
   30    continue
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gradck( n, sgrad, grel, sigma, fv, tol, info)
c
c  subroutine gradck provides the main stopping criterion for
c  oracle.
c------------------------------------------------------------------
      include 'impli.inc'
      dimension sgrad(n), sigma(n)
c
      grel = 0
      fmag = dmax1(dabs(fv), 1.0d0)
      do 10 i = 1,n
        grcomp = dabs( sgrad(i)*dmax1(sigma(i), 1.0d0))/fmag
        grel = dmax1( grel, grcomp)
   10 continue
      if (grel .lt. tol) info = 9
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine hesian( ld, nx, h, g, nq, q)
c
c  subroutine hesian forms the hessian matrix and the gradient vector
c  from the vector of the quadratic coefficients.
c------------------------------------------------------
      include 'impli.inc'
      dimension h(ld,nx), g(nx), q(nq)
c
      do 10 i = 1, nx
         g(i) = - q(i+1)
         h(i,i) = 2.0d0 * q(i+1+nx)
   10 continue
      indexq = 2*nx + 1
      do 20 i = 1, nx-1
         do 30 j = i+1, nx
            indexq = indexq + 1
            h(i,j) = q(indexq)
            h(j,i) = h(i,j)
   30    continue
   20 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine initq( m, nx, q)
c
c  subroutine initq sets the initial vector of quadratic coefficients
c  so that the hessian will start as the identity.
c---------------------------------------------------
      include 'impli.inc'
      dimension q(m)
c
      do 10 i=1,m
         if ( (i .gt. (1+nx)) .and. (i .le. (1+2*nx))) then
            q(i) = 0.5
         else
            q(i) = 0.0
         endif
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine optisort( iter, dist, nseq, nq)
c
c  sorts the first nq entries of the nseq to the
c  list of distances from nearest to farthest.
cryne 8/17/02 name changed from isort to optisort because isort
cryne is a name found in a widely used library.
c---------------------------------------------------------------------
      include 'impli.inc'
      dimension dist(iter), nseq(iter)
c
c     initialize nseq array in reverse order to lessen sorting.
c
      do 10 i=1,iter
         nseq(i) = iter + 1 - i
   10 continue
c
c     make nq passes through the nseq list, sorting on distance.
c
      do 20 j=1,nq
         do 30 k= j+1, iter
            if ( dist(nseq(k)) .lt. dist(nseq(j)) ) then
c
c              switch indices.
c
               itemp = nseq(j)
               nseq(j) = nseq(k)
               nseq(k) = itemp
            endif
   30    continue
   20 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine jtjmul( ldj,ldjtj,nc,fjac,fjtj)
c
c subroutine jtjmul multiplies the jacobian transpose times the jacobian
c-----------------------------------------------------------------------
      include 'impli.inc'
      dimension fjac(ldj,nc), fjtj(ldjtj,nc)
      parameter (maxv = 10)
      dimension temp(maxv), tempt(maxv)
c
      do 10 i= 1,nc
         do 20 j = i,nc
            do 30 k=1,nc
               temp(k) = 0.0
               tempt(k) = 0.0
   30       continue
            call vecopy(i,fjac(1,i),tempt)
            call vecopy(j,fjac(1,j),temp)
            fjtj(i,j) = ddot(nc, tempt, 1, temp, 1)
            fjtj(j,i) = fjtj(i,j)
   20    continue
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lform(maxv,maxf,nx,nf,nq,xu,fu,fjac,info,rcond)
c
c  subroutine lform makes a linear approximation to the list of
c  xv points and corresponding fv values.
c-------------------------------------------------------------------
      include 'impli.inc'
      dimension xu(maxv, nq), fu(maxf,nq)
      dimension fjac(maxf, nx)
c
      parameter (nxdim = 10, maxq = nxdim+1, epsmach = 1.0d-16)
      dimension a(nxdim, nxdim), ascale(nxdim)
      dimension y(nxdim), ipvt(nxdim), w(nxdim)
c
      call nlsa(maxv, nx, nq, a, xu)
      call scale(nxdim, nx, ascale, a, info)
      if (info .ne. 0) return
      call dgeco(a,nxdim,nx,ipvt,rcond,w)
      if (rcond .lt. epsmach) then
         info = 4
         return
      endif
      do 10 i = 1,nf
         call dcopy(nx,fu(i,2),maxf,y,1)
         call dgesl(a,nxdim,nx,ipvt,y,0)
         do 20 j= 1,nx
            fjac(i,j) = y(j) * ascale(j)
   20    continue
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lpick(nxmax,nfmax,nx,nf,nvec,iter,xu,fcen,fmu,fu)
c
c  subroutine lpick selects the xu from among all xv points by
c  choosing the best point, xcen, and the nvec - 1 nearest points
c  to xcen.
c--------------------------------------------------------------------
      include 'impli.inc'
cryne 3/15/06      include 'nlsvar.inc'
cryne 3/15/06 here is nlsvar.inc:

      parameter (maxdat = 1000)
      parameter (maxv = 10, maxf = 15, maxq = maxv+1 )
      common /nlsvar/ modu, ntot, nomin, minpt, maxpt, nlsi, fmin,
     * fmax, regold, region, seed, oldmu, dmu, xdat(maxv,maxdat),
     * fdat(maxf,maxdat), fmdat(maxdat), fvmin(maxf), xcen(maxv)
cryne 3/15/06      save /nlsvar/

      dimension xu(nxmax,nvec), fu(nfmax,nvec), fmu(nvec), fcen(nf)
      dimension  dist(maxdat), nseq(maxdat), xdif(maxv)

cryne 3/15/06
      save
c
c     write(6,*)'HEEEEEEERRRRRREEEEE I am in LPICK'
      do 10 k = 1,iter
         call vecdif(nx, xdat(1,k), xcen, xdif)
         dist(k) = enorm(nx, xdif)
   10 continue
         call optisort( iter, dist, nseq, nvec)
c
c      check that the current point is included in the list
c
         do 20 k=1,nvec
         if(nseq(k).eq.iter) go to 30
   20     continue
c
c     didn't find current point (iter) in the list, so force it
c     by replacing worst of the points already chosen:
c
         nseq(nvec) = iter
c
c      copy the selected points into the working vectors xu. The
c      origin is taken as xcen.
c
   30     do 40 k = 1, nvec
            nu = nseq(k)
            call vecdif( nx, xdat(1,nu), xcen, xu(1,k))
            call vecdif( nf, fdat(1,nu), fcen, fu(1,k))
            fmu(k) = fmdat(nu)
   40    continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine message(info,ipr)
c
c  subroutine message prints out termination messages in the file
c  qout.dat.
c---------------------------------------------------------------------
      include 'impli.inc'
      parameter (epsmach = 1.0d-16, epsqrt = 1.0d-8)
c
      write(ipr,*)
      write(ipr,*) 'done'
      if (info .eq. 1) then
         write(ipr,*) 'the relative difference between '
         write(ipr,*) 'predicted and actual f values agree within ',
     &                 epsqrt
      else if (info .eq. 2) then
         write(ipr,*) 'the relative norm of the difference between'
         write(ipr,*) 'the last two x vectors lies within stol.'
      else if (info .eq. 3) then
         write(ipr,*) 'the relative norm of the difference between'
         write(ipr,*) 'the last two q vectors lies within stol.'
      else if (info .eq. 4) then
         write(ipr,*)'exit due to ill-conditioning of the q-coefficient'
         write(ipr,*) 'matrix, a. '
      else if (info .eq. 5) then
         write(ipr,*) 'exit due to ill-conditioning of the hessian.'
      else if (info .eq. 6) then
         write(ipr,*)'exit due to attempted step of enormous magnitude.'
         write(ipr,*)'If using mss, try increasing parameter 4'
      else if (info .eq. 7) then
         write(ipr,*) 'absolute f error lies within ',epsqrt
      else if (info .eq. 8) then
         write(ipr,*) 'absolute x error lies within stol.'
      else if (info .eq. 9) then
         write(ipr,*) 'relative gradient is less than gtol.'
      else if (info .eq. 10) then
         write(ipr,*) 'maximum number of iterations exceeded.'
      else if (info .eq. 11) then
         write(ipr,*) 'norm of residual lies within ftol.'
      else if (info .eq. 12) then
         write(ipr,*) 'successive steps have an identical component.'
      else if (info .eq. 13) then
         write(ipr,*) 'one of user-specified weights is zero.'
      else if (info .eq. 14) then
         write(ipr,*) '# of equations is less than # of unkowns.'
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mvmult(lda,nr,nc,a,x,result)
c
c subroutine mvmult multiplies a vector, x, by a matrix, a.
c-----------------------------------------------------------------------
      include 'impli.inc'
      dimension a(lda,nc), x(nc), result(nr)
      parameter (maxv = 10)
      dimension xtemp(maxv)
c
      call vecopy( nc, x, xtemp)
      do 10 i = 1, nr
         result(i) = ddot(nc, a(i,1), lda, xtemp, 1)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine newton( ld, nv, h, g, step, info)
c
c  subroutine newton calculates the newton step from the best point, xce
c  which is now the origin in the transformed frame of reference.
c-----------------------------------------------------------------------
      include 'impli.inc'
      parameter (maxv = 10, epsmach = 1.0d-16)
      dimension h(ld,nv), g(nv), step(nv)
      dimension w(maxv), ipvt(maxv), ht(maxv,maxv)
c
      do 10 i= 1,nv
         call vecopy( nv, h(1,i), ht(1,i))
   10 continue
      call dsico( ht, maxv, nv, ipvt, rcond, w)
      if (rcond .lt. epsmach) then
         info = 5
         return
      endif
      call vecopy(nv, g, step)
      call dsisl( ht, maxv, nv, ipvt, step)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine nls(mode,iter,nx,xv,nf,fv,iscale,wts,tol,info)
c
c subroutine nls solves nonlinear equations and least squares problems.
c
c written by tom mottershead and kurt overley of los alamos national
c laboratory in august 1989.
c
c
c     variables:
c
c     --> mode: running mode.
c               = 0 for normal unconstrained minimization
c               = 1 for constrained minimization in the box defined by
c                   (xlo(i),xhi(i)), i=1,nx
c               = 2 to reload internal data commons and return without
c                   computing new point.
c               = 3 for inquiry. The best point so far is returned in
c                   (xv,fv), its index is returned in iter, and the tota
c                   points stored is returned as info = - ntot
c
c     --> iter: the iteration counter. should be incremented before ever
c               call in normal running.
c                iter = +n  means (xv,fv) is stored and used as nth data
c                iter = 0 or 1 causes initialization, clearing stored da
c                iter = -n means return nth stored point in (xv,fv)
c
c     --> nx: the number of variables.
c
c     <--> xv(nx): the current point on input, the new point on output.
c
c     <--> fv: the residual vector of function values at the current poi
c              on input.
c
c     --> iscale: controls x and f scaling.
c                  iscale = 0: no scaling.  this is the recommended
c                              default value. if unsatisfactory performa
c                              results, the problem may be poorly scaled
c                              try the various scalings.
c                  iscale = 1: auto x scaling
c                  iscale = 2: auto f scaling
c                  iscale = 3: auto x and f scaling
c                  iscale = 4: scale f with wts, no x scaling
c                  iscale = 5: scale f with wts, auto x scaling
c
c     --> wts(nf): a vector of weights to scale f.  f(j) is replaced
c                  by w(j)*f(j).
c
c     --> tol(j), j=1,4
c               tol(1) = gtol: if the relative gradient lies within the
c                        tolerance, gtol, nls exits.
c               tol(2) = stol: nls exits if the absolute or relative
c                        difference between the the suggested minimizer,
c                        xnew, and the current point, xc, is less
c                        than stol.
c               tol(3) = ftol: if the magnitude of the residual f vector
c                        is less than ftol, qadnls exits.
c               tol(4) = initial step fraction stpfrc (default = 0.1)
c
c     <-- info: on return info contains information as to the
c               termination status of nls.
c               info = 0: no termination yet.
c                    = 1: the relative difference between the predicted
c                         and actual values agrees within ftol.
c                    = 2: the relative difference between the last two
c                         xv vectors lies within stol.
c                    = 3: the relative difference between the last two
c                         q vectors lies within qtol.
c                    = 4: exit due to ill-conditioning of the q-coeffici
c                         matrix, a.
c                    = 5: exit due to ill-conditioning of the hessian.
c                    = 6: exit due to attempted step of enormous magnitu
c                    = 7: absolute fv error lies within ftol.
c                    = 8: absolute xv error lies within stol.
c                    = 9: relative gradient is less than gtol.
c                    = 10: maximum iterations exceeded (this must be
c                          determined by the program calling nls).
c                    = 11: absolute norm of residual function vector lie
c                          within ftol.
c                    = 12: successive steps have an identical x componen
c                    = 13: one of the user-specified weights is zero.
c                    = 14: number of equations is less than number of un
c
c  parameters:
c
c     maxv: the maximum number of variables allowed.
c
c     maxf: the maximum number of equations allowed.
c
c     maxq: the maximum number of points needed for fitting to a
c           quadratic form.
c
c     maxdat: the maximum number of points that can be stored.
c
c     epsmach: the machine precision.
c
c     gtol: the relative gradient tolerance.  if the current relative
c           gradient falls below gtol, nls exits.
c
c     stol: the step tolerance.  if two successive points lie
c           within stol, nls exits.
c
c     ftol: the residual tolerance.  if the norm of the current residual
c           falls below ftol, nls exits.
c
c     stpfrc: the initial random steps are of length stpfrc * the
c             norm of the current best point.
c
c     trfrac: the initial trust region is set to trustfrac times
c               the 2 norm of a vector of unit components, representativ
c               of an "average" shifted and rms scaled xv vector.
c
c  internal variables:
c
c     xcen(maxv): the best point so far.  it is the center of the linear
c                 approximation to the nonlinear system.
c
c     xu(maxv,maxq): contains a list of the last nv + 1 points.
c
c     fu(maxf,maxq): contains a list of the last nv + 1 function vectors
c
c     fmu(maxq): the list of previous norms of residual function values.
c
c     xdat(maxv,maxdat): the entire list of all points.
c
c     fdat(maxf,maxdat): the entire list of all function vectors.
c
c     fvmin(maxf): the best function vector so far.
c
c     fcen(maxf): forms the centering function vector for the linear
c                 approximation to the nonlinear system.  it is the same
c                 as fvmin.
c
c     fmdat(maxdat): the entire list of norms of residual function vecto
c
c     minpt: location of the minimum fv and xv in the data pool.
c
c     maxpt: a pointer to the maximum fv and xv in the list.
c
c     a(maxq,maxq): the quadratic coefficient matrix.
c
c     ascale(maxq): a scaling vector that normalizes the columns
c                   of a to unit euclidean norm.
c
c     y(maxq): the solution of the scaled a matrix.
c
c     q(maxq): the vector of coefficients to the quadratic form.
c
c     g(maxv): the gradient vector.
c
c     fjac(maxv,maxv): the jacobian of the nonlinear system of equations
c                      or for nonlinear least squares is the jacobian of
c                      the norm of the function residual vector.
c
c     step(maxv): the trust region newton step (nlstep).
c
c     nq = nx+1 : the number of points needed to
c                 make a linear approximation to the jacobian of the
c                 norm of the function vector.
c
c     nomin: the number of iterations since a new minimum point
c            has been found.
c
c     region: the region or radius about the current best point
c               in which the quadratic model is trusted to accurately
c               model the function.
c
c     regold: the old trust region.
c
c     dmu: the levenberg-marquardt parameter that adjusts the hessian
c          to solve the constrained problem of minimizing the quadratic
c          form subject to lying in a trust region about the best point.
c
c     oldmu: the old levenberg-marquardt parameter.
c
c     sgrad: the gradient of the function with respect to a step
c            from the minimum point.  note that this gradient is
c            different from the gradient of fv with respect to xv.
c
c     sigma: a vector containing the diagonal elements of a
c             scaling matrix used to normalize each of the components of
c             xcen to one.
c
c     fsigma: a similar scaling vector for the function vectors.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'impli.inc'
cryne 3/15/06      include 'nlsvar.inc'
      dimension xv(nx), fv(nf), tol(4), wts(nf)
      parameter (epsmach = 1.0d-16)
      parameter ( trfrac = 2.0d0)

cryne 3/15/06 here is nlsvar.inc:
      parameter (maxdat = 1000)
      parameter (maxv = 10, maxf = 15, maxq = maxv+1 )
      common /nlsvar/ modu, ntot, nomin, minpt, maxpt, nlsi, fmin,
     * fmax, regold, region, seed, oldmu, dmu, xdat(maxv,maxdat),
     * fdat(maxf,maxdat), fmdat(maxdat), fvmin(maxf), xcen(maxv)
cryne 3/15/06      save /nlsvar/

      dimension xu(maxv, maxq), xnew(maxv)
      dimension fu(maxf, maxq), fmu(maxq), fcen(maxf)
      dimension fjac(maxf, maxv), g(maxv)
      dimension sigma(maxv), fsigma(maxf), step(maxv)
cryne 3/15/06
      save
c
      if (info .gt. 0) return
      if (nf .lt. nx) then
         info = 14
         return
      endif
      modu = mode
c
c       nq is number of points needed for full linear form.
c       nvec is the number available this iteration.
c
      nq = nx + 1
      nvec = min0(nq,iter)
      if(iter.gt.ntot) ntot=iter
      gtol = tol(1)
      stol = tol(2)
      ftol = tol(3)
      stpfrc = 0.1
      if(tol(4).gt.0.0) stpfrc = tol(4)
c     type *, ' nls:',nq,'=nq',iter,'=iter',nvec,'=nvec'
      if(iter.lt.0) return
      if (iter .gt. 1) go to 60
c
c      initialize varibles on first iteration
c
      seed = 1.0d0
      grel = 0.0
      fmin = 1.e30
      fmax = -fmin
      minpt = 1
      maxpt = 1
      ntot = 0
      nomin = 0
      regold = 0.0
      region = 0.0
      oldmu = 0.0
      dmu = 0.0
      do 10 i = 1,nx
         xnew(i) = 0.0
         xcen(i) = 0.0
         g(i) = 0.0
         sigma(i) = 0.0
         step(i) = 0.0
   10 continue
      do 15 i = 1,nf
         fcen(i) = 0.0
         fvmin(i) = 0.0
         fsigma(i) = 0.0
   15 continue
      do 20 j = 1,nq
         fmu(j) = 0.0
         do 25 i = 1,nf
            fu(i,j) = 0.0
   25    continue
         do 30 i = 1,nx
            xu(i,j) = 0.0
   30    continue
   20 continue
      do 40 j = 1,maxdat
         fmdat(j) = 0.0
         do 45 i = 1,nf
            fdat(i,j) = 0.0
   45    continue
         do 50 i = 1,nx
            xdat(i,j) = 0.0
   50    continue
   40 continue
c
c   add new xv & fv to the data pool.
c
  60  continue
      call vecopy( nx, xv, xdat(1,iter))
      call vecopy( nf, fv, fdat(1,iter))
      call vecopy( nf, fvmin, fcen)
      fmerit = 0.5d0 * (enorm( nf,fv) ** 2.0d0)
      fmdat(iter) = fmerit
c
c        check for new minimum point
c
      if(fmerit.lt.fmin) then
         fmin = fmerit
         minpt = iter
         call vecopy(nx,xv,xcen)
         call vecopy(nf,fv,fcen)
         call vecopy(nf,fv,fvmin)
      endif
c
c        check for new maximum point
c
      if(fmerit.gt.fmax) then
         fmax = fmerit
         maxpt = iter
      endif
c
c      random steps about the minimum the first nx times:
c
c     write(6,*)'INSIDE NLS with iter, nx=',iter,nx
      if (iter .le. nx) then
         call ranstep(iter, nx, seed, xcen, xv, stpfrc)
         return
      endif
c-------------------------------------------------------------
c    Normal running: fit and solve the linear approximation when
c    there are enough points:
c-------------------------------------------------------------
c
c      lpick the best points - xu - to fit to a quadratic form:
c
      call lpick(maxv,maxf,nx,nf,nq,iter,xu,fcen,fmu,fu)
c
c        normalize the working points - xu
c
      if ((iscale.eq.1) .or. (iscale.eq.3) .or. (iscale.eq.5)) then
         call nlscal( maxv, nx, nq, nvec, sigma, xu)
      endif
      if (iscale.ge.2) then
         call scalfu( maxf,nf,nq,fsigma,fu,fcen,iscale,wts,info)
      endif
c
c       determine quadratic form
c
      call lform(maxv,maxf,nx,nf,nq,xu,fu,fjac,info,rcond)
      if (info .gt. 0) go to 99
c
c      reset trust region
c
      call truset(nx, nq, fmu, fpred, fmerit,
     & trfrac, iter, nomin, minpt, nvec, regold, region)
c
c       compute step on surface of trust region
c
      call nlstep(maxf, nx, nf, fjac, fcen, g, step, trfrac,
     & nq, iter, regold, region, oldmu, dmu, info)
      if (info .gt. 0) go to 99
c
c       compute xnew by adding rescaled step to xcen
c
      do 80 i = 1,nx
         if ((iscale.eq.1).or.(iscale.eq.3).or.(iscale.eq.5)) then
            xnew(i) = xcen(i) + step(i)*sigma(i)
         else
            xnew(i) = xcen(i) + step(i)
         endif
   80 continue
c
c      check whether xnew has moved enough to continue the search
c
      call xstop(nx, xv, xnew, stol, info)
c
c     put xnew in the xv vector for the return
      call vecopy(nx, xnew, xv)
      if (info .gt. 0) go to 99
c
c       check gradient for convergence
c
cryne not sure I did this right, but here is how I changed this on 29 April 2008:
cryne call gradck( nx, g, grel, sigma, fv, gtol, info)
      fvmaxabs=maxval(abs(fv))
      call gradck( nx, g, grel, sigma, fvmaxabs, gtol, info)
      if (fmerit .lt. ftol) info = 11
   99 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine nlsa (maxv, nx, nq, a, xu)
c
c  subroutine nlsa forms the linear system in the center of mass
c  coordinates whose solution yields the quadratic coefficients.
c---------------------------------------------------------------
      include 'impli.inc'
      dimension a(maxv,nx), xu(maxv, nq)
c
      do 10 k = 1, nx
         call dcopy(nx,xu(1,k+1),1,a(k,1),maxv)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine nlscal( nxmax, nx, nq, nvec, sigma, xu)
c
c  subroutine nlscal finds the rms values of each of the components
c  in the xu.  then each of the xu vector components is scaled
c  by the inverse of the rms values.
c------------------------------------------------------------------
      include 'impli.inc'
      dimension sigma(nx), xu(nxmax, nq)
c
      tol = 1.d-15
      denom = nvec - 1
c
c       loop over all components
c
      do 40 i = 1, nx
      sigma(i) = 0.0
c
c      collect standard deviation of all components
c
         do 20 j = 2, nvec
         sigma(i) = sigma(i) + xu(i,j)**2
   20    continue
      sigma(i) = dsqrt(sigma(i)/ denom)
c
c       renormalize all vectors
c
      xdenom = dabs(xu(i,2))
      if(xdenom.le.tol) then
         sigma(i) = 1.0d0
         go to 40
      else
         sigma(i) = xdenom/ sigma(i)
         do 30 j = 2, nvec
           xu(i,j) = xu(i,j)/sigma(i)
   30    continue
      endif
   40 continue
      return
      end
c
c*********************************************************************
c
      subroutine nlstep( mf, nx, nf, fjac, fcen, g, step, trfrac,
     &   nq, iter, regold, region, oldmu, dmu, info)
c
c  nlstep (trust region step) calculates the new xv point that solves th
c  constrained problem of minimizing a quadratic functional subject to
c  lying in trust region about the best xv point.
c------------------------------------------------------------
      include 'impli.inc'
      dimension fjac(mf, nx), g(nx), step(nx), fcen(nf)
c
      parameter (maxv = 10, maxf = 15)
      dimension fjact(maxv, maxf),fjtj(maxv,maxv),fjtjmu(maxv,maxv)
      dimension qraux(maxv), rsd(maxf)
c
c c   first try the newton step.
c
      factor = -1.0d0
      call dscal(nf, factor, fcen, 1)
      call trnsps(maxf,maxv,nf,nx,fjac,fjact)
      call mvmult(maxv,nx,nf,fjact,fcen,g)
      call dqrdc(fjac,maxf,nf,nx,qraux,jdum,dum,0)
      call condn(maxf,nx,fjac,info)
      if (info .gt. 0) return
      if (nf .eq. nx) then
         call dqrsl(fjac,maxf,nf,nx,qraux,fcen,dum,rsd,step,rsd,
     &              dum,110,info)
      else
         call vecopy( nx,g,step)
         call dposl( fjac,maxf,nx,step)
      endif
      if ( info .gt. 0) return
      stepnm = enorm(nx, step)
      regmax = 1.5*region
      regmin = 0.75*region
      if (stepnm .lt. regmax) then
         oldmu = 0.0
         region = stepnm
         return
      endif
c
c c   newton step has failed, so calculate a dmu such that step(dmu)
c c   is fairly close to the region.
c
      call jtjmul(maxf,maxv,nx,fjac,fjtj)
      dmulow = 0.0
      gnorm = enorm(nx, g)
      dmuhi = gnorm/ region
      call findhi( dmuhi, maxv, nx, fjtj, fjtjmu,
     &             g, step, region, info)
      if (info .gt. 0) return
      if (oldmu .eq. 0.0) then
         dmu = dmuhi/ 2.0d0
      else
         dmu = (regold/ region) * oldmu
      endif
      if ( (dmu .lt. dmulow) .or. (dmuhi .lt. dmu) ) then
         dmu = dmulow + (dmuhi - dmulow)/ 2.0d0
      endif
c
c c   at 40 the mu loop begins by calculating a new step.
c
   40 continue
      do 30 i = 1,nx
         call vecopy( nx, fjtj(1,i), fjtjmu(1,i))
         fjtjmu(i,i) = fjtj(i,i) + dmu
   30 continue
      call newton( maxv, nx, fjtjmu, g, step, info)
      if (info .gt. 0) return
c      call dpofa(fjtjmu,maxv,nx,info)
c      if (info .gt. 0) then
c         info = 5
c         return
c      endif
c      call vecopy(nx, g, step)
c      call dposl(fjtjmu,maxv,nx,step)
      stepnm = enorm(nx, step)
      if ( (regmin .lt. stepnm) .and. (stepnm .lt. regmax) ) then
         go to 99
      endif
c
c c   update the mu bounds and get a new mu.
c
      phi = stepnm - region
      if (phi .gt. 0.0) then
         dmulow = dmu
      else
         dmuhi = dmu
      endif
      if (dmuhi .le. dmulow) then
         go to 99
      endif
       dmu = dmulow + (dmuhi - dmulow)/ 2.0d0
      go to 40
c
c c   at 99 do exit chores.
c
   99 continue
      oldmu = dmu
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine qform(maxv, nx, nq, nvec, xu, fu, h, g, info, rcond)
c
c  subroutine qform fits a quadratic form (finds the hessian matrix, h,
c  and the gradient, g) to the list of xv points and corresponding
c  fv values.
c-------------------------------------------------------------------
      include 'impli.inc'
      dimension xu(maxv, nq), fu(nq)
      dimension h(maxv, nx), g(nx)
c
      parameter (nxdim = 10, maxq = (nxdim+1)*(nxdim+2)/2)
      dimension a(maxq, maxq), ascale(maxq)
      dimension y(maxq), q(maxq)
c
c     type *, ' qform:',maxv,'=maxv',nx,'=nx',nvec,'=nvec',nq,'=nq'
      if (nvec .lt. nq) call initq( nq, nx, q)
      call geta(maxv, maxq, nx, nvec, a, xu)
cryne 3/15/06      call scale(maxq, nvec, ascale, a)
      call scale(maxq, nvec, ascale, a, info)
      call solve(maxq, nvec, a, y, fu, info, rcond)
      if (info .gt. 0) return
      do 10 i = 1,nvec
         q(i) = y(i) * ascale(i)
   10 continue
      call hesian(maxv, nx, h, g, nq, q)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine qso(mode,iter,nx,xv,fv,iscale,tol,info)
c
c subroutine qadmin is an unconstrained minimization routine that employ
c quadratic form fitting with a trust region step control stategy.
c
c written by tom mottershead and kurt overley of los alamos national
c laboratory in august 1988.
c   Renamed  22 Mar 91  to QFM = Quadratic Fit Minimizer (at first
c   installation in Marylie. Still not happy with name.
c
c     variables:
c
c     --> mode: running mode.
c               = 0 for normal unconstrained minimization
c               = 1 for constrained minimization in the box defined by
c                   (xlo(i),xhi(i)), i=1,nx
c               = 2 to reload internal data commons and return without
c                   computing new point.
c               = 3 for inquiry. The best point so far is returned in
c                   (xv,fv), its index is returned in iter, and the tota
c                   points stored is returned as info = - ntot
c
c     --> iter: the iteration counter. should be incremented before ever
c               call in normal running.
c                iter = +n  means (xv,fv) is stored and used as nth data
c                iter = 0 or 1 causes initialization, clearing stored da
c                iter = -n means return nth stored point in (xv,fv)
c
c     --> nx: the number of variables.
c
c     <--> xv(nx): the current point on input, the new point on output.
c
c     <--> fv: the function value at the current point on input, the new
c            predicted value of the minimun on output
c
c     --> iscale: controls x.
c                  iscale = 0: no scaling.  this is the recommended
c                              default value. if unsatisfactory performa
c                              results, the problem may be poorly scaled
c                              try scaling.
c                  iscale = 1: auto x scaling.
c
c     --> tol(j), j=1,4
c               tol(1) = gtol: if the relative gradient lies within the
c                        tolerance, gtol, qadmin exits.
c               tol(2) = stol: qadmin exits if the absolute or relative
c                        difference between the the suggested minimizer,
c                        xnew, and the current point, xc, is less
c                        than stol.
c
c               tol(3) = ftol: not used by qadmin, but by nls.
c               tol(4) = initial step fraction stpfrc (default = 0.1)
c
c     <-- info: on return info contains information as to the
c               termination status of qadmin.
c               info = 0: no termination yet.
c                    = 1: the relative difference between the predicted
c                         and actual values agrees within ftol.
c                    = 2: the relative difference between the last two
c                         xv vectors lies within stol.
c                    = 3: the relative difference between the last two
c                         q vectors lies within qtol.
c                    = 4: exit due to ill-conditioning of the q-coeffici
c                         matrix, a.
c                    = 5: exit due to ill-conditioning of the hessian.
c                    = 6: exit due to attempted step of enormous magnitu
c                    = 7: absolute fv error lies within ftol.
c                    = 8: absolute xv error lies within stol.
c                    = 9: relative gradient is less than gtol.
c
c  parameters:
c
c     maxv: then maximum number of variables allowed.
c
c     maxq: the maximum number of points needed for fitting to a
c           quadratic form.
c
c     maxdat: the maximum number of points that can be stored.
c
c     epsmach: the machine precision.
c
c     gtol: the relative gradient tolerance.  if the current relative
c           gradient falls below gtol, qadmin exits.
c
c     stol: the step tolerance.  if two successive points lie
c           within stol, qadmin exits.
c
c     ftol: the residual tolerance.  if the norm of the current residual
c           falls below ftol, qadmin exits.
c
c     stpfrc: the initial random steps are of length stpfrc * the
c             norm of the current best point.
c
c     trfrac: the initial trust region is set to trustfrac times
c               the 2 norm of a vector of unit components, representativ
c               of an "average" shifted and rms scaled xv vector.
c
c  internal variables:
c
c     xu(maxv,maxq): contains a list of the last (n+1)(n+2)/2 points.
c
c     fu(maxq): the list of previous function values.
c
c     xdat(maxv,maxdat): the entire list of all points.
c
c     fdat(maxdat): the entire list of function values.
c
c     minpt: location of the minimum fv and xv in the data pool.
c
c     maxpt: a pointer to the maximum fv and xv in the list.
c
c     a(maxq,maxq): the quadratic coefficient matrix.
c
c     ascale(maxq): a scaling vector that normalizes the columns
c                   of a to unit euclidean norm.
c
c     y(maxq): the solution of the scaled a matrix.
c
c     q(maxq): the vector of coefficients to the quadratic form.
c
c     h(maxv,maxv): the second derivative hessian matrix.
c
c     g(maxv): the gradient vector.
c
c     step(maxv): the trust region newton step (trstep).
c
c     nq = (nx+1)*(nx+2)/2 : the number of points needed to
c                            determine a quadratic form.
c
c     nomin: the number of iterations since a new minimum point
c            has been found.
c
c     region: the region or radius about the current best point
c               in which the quadratic model is trusted to accurately
c               model the function.
c
c     regold: the old trust region.
c
c     dmu: the levenberg-marquardt parameter that adjusts the hessian
c          to solve the constrained problem of minimizing the quadratic
c          form subject to lying in a trust region about the best point.
c
c     oldmu: the old levenberg-marquardt parameter.
c
c     sgrad: the gradient of the function with respect to a step
c            from the minimum point.  note that this gradient is
c            different from the gradient of fv with respect to xv.
c
c     sigma: a vector containing the diagonal elements of a
c             scaling matrix used to normalize each of the components of
c             xcen to one.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      include 'impli.inc'
cryne--- 23 March 2006      include 'minvar.inc'
      parameter ( maxdat = 1000)
      parameter (maxv = 10, maxq = (maxv+1)*(maxv+2)/2 )
      common /minvar/ jtyp, modu, ntot, nomin, minpt, maxpt, fmin, fmax,
     * regold, region, seed, oldmu, dmu, xdat(maxv,maxdat), fdat(maxdat)
cryne      save /minvar/
cryne---
      dimension xv(nx), tol(4)
      parameter (epsmach = 1.0d-16)
      parameter ( trfrac = 2.0d0)
      dimension xu(maxv, maxq), fu(maxq), xcen(maxv), xnew(maxv)
      dimension h(maxv, maxv), g(maxv)
      dimension sigma(maxv), step(maxv)
cryne--- 23 march 2006:
      save
c
      write(6,*)'------****-----INSIDE QSO with iter=',iter
c
      if (info .gt. 0) return
      modu = mode
c
c       nq is number of points needed for full quadratic form.
c       nvec is the number available this iteration.
c
      nq = (nx + 1)*(nx + 2)/2
      nvec = min0(nq,iter)
      if(iter.gt.ntot) ntot=iter
      stol = tol(1)
      gtol = tol(2)
      stpfrc = 0.1
      if(tol(4).gt.0.0) stpfrc = tol(4)
c     type *, ' qadmin:',nq,'=nq',iter,'=iter',nvec,'=nvec'
      if(iter.lt.0) return
      if (iter .gt. 1) go to 60
c
c      initialize varibles on first iteration
c
      seed = 1.0d0
      grel = 0.0
      fmin = 1.e30
      fmax = -fmin
      minpt = 1
      maxpt = 1
      ntot = 0
      nomin = 0
      regold = 0.0
      region = 0.0
      oldmu = 0.0
      dmu = 0.0
      do 10 i = 1,nx
         xnew(i) = 0.0
         xcen(i) = 0.0
         g(i) = 0.0
         sigma(i) = 0.0
         step(i) = 0.0
   10 continue
      do 20 j = 1,nq
         fu(j) = 0.0
         do 30 i = 1,nx
            xu(i,j) = 0.0
   30    continue
   20 continue
      do 40 j = 1,maxdat
         fdat(j) = 0.0
         do 50 i = 1,nx
            xdat(i,j) = 0.0
   50    continue
   40 continue
c
c   add new xv & fv to the data pool.
c
  60  continue
      call vecopy( nx, xv, xdat(1,iter))
      fdat(iter) = fv
c
c        check for new minimum point
c
      if(fv.lt.fmin) then
         fmin = fv
         minpt = iter
         call vecopy(nx,xv,xcen)
      endif
c
c        check for new maximum point
c
      if(fv.gt.fmax) then
         fmax = fv
         maxpt = iter
      endif
c
c      random steps about the minimum the first nx times:
c
      if (iter .le. nx) then
         call ranstep(iter, nx, seed, xcen, xv, stpfrc)
         return
      endif
c-------------------------------------------------------------
c    Normal running: fit and solve the quadratic form when
c    there are enough points:
c-------------------------------------------------------------
c
c      select the best points - xu - to fit to a quadratic form:
c
      call bestpt(maxv,nx,nvec,iter,xcen,xu,fu)
c
c        normalize the working points - xu
c
      if (iscale .eq. 1) then
         call scalxu( maxv, nx, nq, nvec, sigma, xu)
      endif
c
c       determine quadratic form
c
      call qform(maxv, nx, nq, nvec, xu, fu, h, g, info, rcond)
      if (info .gt. 0) go to 99
c
c      reset trust region
c
      call truset(nx, nq, fu, fpred, fv,
     & trfrac, iter, nomin, minpt, nvec, regold, region)
c
c       compute step on surface of trust region
c
      call trstep(maxv, nx, h, g, step, trfrac,
     & nq, iter, regold, region, oldmu, dmu, info)
      if (info .gt. 0) go to 99
c
c       compute xnew by adding rescaled step to xcen
c
      do 80 i = 1,nx
         if (iscale .eq. 1) then
            xnew(i) = xcen(i) + step(i)*sigma(i)
         else
            xnew(i) = xcen(i) + step(i)
         endif
   80 continue
c
c      check whether xnew has moved enough to continue the search
c
      call xstop(nx, xv, xnew, stol, info)
c
c     put xnew in the xv vector for the return
      call vecopy(nx, xnew, xv)
      if (info .gt. 0) go to 99
c
c       check gradient for convergence
c
      call gradck( nx, g, grel, sigma, fv, gtol, info)
   99 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ranstep( iter, nv, seed, xcen, xnew, stpfrc)
c
c  subroutine ranstep generates the first n+1 points by taking random
c  steps from the current best point.  note there is no protection
c  for components of xnew to remain positive.
c--------------------------------------------------------------------
      include 'impli.inc'
      dimension xcen(nv), xnew(nv)
c
c
      call vecopy( nv, xcen, xnew)
      do 10 i=1,nv
c        seed = ran( idint( seed*1.0d8 + 1.0d0) )
c        xnew(i) = seed - 0.5d0
      call random_number(x)
c     write(6,*)'i,x=',i,x
      xnew(i) = x - 0.5d0
   10 continue
      cnorm = enorm(nv, xnew)
      xnorm = enorm(nv, xcen)
      if (xnorm .eq. 0.) then
         dnv = nv
         xnorm = dsqrt(dnv)
      endif
      factor = stpfrc * xnorm/cnorm
c     write(6,*)'cnorm,xnorm=',cnorm,xnorm
c     write(6,*)'stpfrc,factor=',stpfrc,factor
      do 20 i=1,nv
         xnew(i) = xcen(i) + factor * xnew(i)
   20 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scale(ld, n, vscale, dmat, info)
c
c  subroutine scale scales a matrix so its columns have unit norm.
c----------------------------------------------------------------
      include 'impli.inc'
      dimension vscale(n), dmat(ld,ld)
c
      do 10 i = 1,n
         denom = enorm(n, dmat(1,i))
         if (denom .eq. 0.) then
            info = 12
            return
         endif
         vscale(i) = 1.0d0/denom
   10 continue
      do 30 i = 1,n
         do 40 j = 1,n
            dmat(j,i) = dmat(j,i) * vscale(i)
   40    continue
   30 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scalfu( nfmax, nf, nq, fsigma, fu, fcen,
     &                   iscale, wts, info)
c
c subroutine scalfu scales the list of function vectors using the same
c automatic scaling as the x points or else the user-specified weight
c vector.
c----------------------------------------------------------------------
      include 'impli.inc'
      dimension fsigma(nf), fu(nfmax, nq), fcen(nf), wts(nf)
c
      if (iscale.ge.4) then
         do 20 i=1,nf
            if (wts(i) .eq. 0.0) then
               info = 13
               return
            endif
            fsigma(i) = 1.0d0/ wts(i)
            do 30 j = 2, nq
               fu(i,j) = fu(i,j)/ fsigma(i)
   30       continue
   20    continue
      else
         call nlscal( nfmax, nf, nq, nq, fsigma, fu)
      endif
      do 10 i = 1,nf
         fcen(i) = fcen(i)/ fsigma(i)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scalxu( nxmax, nx, nq, nvec, sigma, xu)
c
c  subroutine scalxu finds the rms values of each of the components
c  in the xu.  then each of the xu vector components is scaled
c  by the inverse of the rms values - i.e. normalized to unit variance.
c------------------------------------------------------------------
      include 'impli.inc'
      dimension sigma(nx), xu(nxmax, nq)
c
      tol = 1.d-15
      denom = nvec - 1
c
c       loop over all components
c
      do 40 i = 1, nx
      sigma(i) = 0.0
c
c      collect standard deviation of all components
c
         do 20 j = 2, nvec
         sigma(i) = sigma(i) + xu(i,j)**2
   20    continue
      sigma(i) = dsqrt(sigma(i)/ denom)
c
c       renormalize all vectors to unit variance
c
      if(sigma(i).le.tol) go to 40
        do 30 j = 2, nvec
        xu(i,j) = xu(i,j)/sigma(i)
   30   continue
   40 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine bestpt(nxmax,nx,nvec,iter,xcen,xu,fu)
c
c  subroutine bestpt selects the xu from among all xv points by
c  choosing the best point, xcen, and the nvec - 1 nearest points
c  to xcen.
c--------------------------------------------------------------------
      include 'impli.inc'
cryne--- 23 March 2006      include 'minvar.inc'
      parameter ( maxdat = 1000)
      parameter (maxv = 10, maxq = (maxv+1)*(maxv+2)/2 )
      common /minvar/ jtyp, modu, ntot, nomin, minpt, maxpt, fmin, fmax,
     * regold, region, seed, oldmu, dmu, xdat(maxv,maxdat), fdat(maxdat)
cryne      save /minvar/
cryne---
      dimension xu(nxmax,nvec), fu(nvec), xcen(nx)
      dimension  dist(maxdat), nseq(maxdat), xdif(maxv)
cryne--- 23 march 2006:
      save
c
      do 10 k = 1,iter
         call vecdif(nx, xdat(1,k), xcen, xdif)
         dist(k) = enorm(nx, xdif)
   10 continue
         call optisort( iter, dist, nseq, nvec)
c
c      check that the current point is included in the list
c
         do 20 k=1,nvec
         if(nseq(k).eq.iter) go to 30
   20     continue
c
c     didn't find current point (iter) in the list, so force it
c     by replacing worst of the points already chosen:
c
         nseq(nvec) = iter
c
c      copy the selected points into the working vectors xu. The
c      origin is taken as xcen.
c
   30     do 40 k = 1, nvec
            nu = nseq(k)
            call vecdif( nx, xdat(1,nu), xcen, xu(1,k))
            fu(k) = fdat(nu)
   40    continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine solve(ld, n, dmat, dx, data, info, rcond)
c
c  subroutine solve uses the linpack subroutine dgeco and dgesl to
c  solve a general linear system.
c----------------------------------------------
      include 'impli.inc'
      dimension dmat(ld,ld), dx(n), data(n)
      parameter (maxv = 10, maxq = (maxv+1)*(maxv+2)/2 )
      parameter (epsmach = 1.0d-16)
      dimension w(maxq), ipvt(maxq), tempd(maxq,maxq)
c
      do 10 i = 1,n
         call vecopy( n, dmat(1,i), tempd(1,i))
   10 continue
      call dgeco( tempd, maxq, n, ipvt, rcond, w)
      if (rcond .lt. epsmach) then
         info = 4
         return
      endif
      call vecopy( n, data, dx)
      call dgesl( tempd, maxq, n, ipvt, dx, 0)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine trnsps( lda, ldat, nr, nc, a, at)
c
c subroutine trnsps transposes a matrix.
c-----------------------------------------------------------------------
      include 'impli.inc'
      dimension a(lda,nc), at(ldat,nr)
c
      do 10 i = 1,nc
         call dcopy(nr, a(1,i),1,at(i,1),ldat)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine trstep( ld, nx, h, g, step, trfrac,
     &   nq, iter, regold, region, oldmu, dmu, info)
c
c  trstep (trust region step) calculates the new xv point that solves th
c  constrained problem of minimizing a quadratic functional subject to
c  lying in trust region about the best xv point.
c------------------------------------------------------------
      include 'impli.inc'
      dimension h(ld,nx), g(nx), step(nx)
c
      parameter (maxv = 10)
      dimension hmu(maxv,maxv)
c
c c   first try the newton step.
c
      call newton( ld, nx, h, g, step, info)
      if ( info .gt. 0) return
      stepnm = enorm(nx, step)
c      if (iter .eq. nq) then
c         region = stepnm
c      endif
      regmax = 1.5*region
      regmin = 0.75*region
      if (stepnm .lt. regmax) then
         oldmu = 0.0
         region = stepnm
         return
      endif
c
c c   newton step has failed, so calculate a dmu such that step(dmu)
c c   is fairly close to the region.
c
      dmulow = 0.0
      gnorm = enorm(nx, g)
      dmuhi = gnorm/ region
      call findhi( dmuhi, maxv, nx, h, hmu,
     &             g, step, region, info)
      if (info .gt. 0) return
      if (oldmu .eq. 0.0) then
         dmu = dmuhi/ 2.0d0
      else
         dmu = (regold/ region) * oldmu
      endif
      if ( (dmu .lt. dmulow) .or. (dmuhi .lt. dmu) ) then
         dmu = dmulow + (dmuhi - dmulow)/ 2.0d0
      endif
c
c c   at 40 the mu loop begins by calculating a new step.
c
   40 continue
      do 30 i = 1,nx
         call vecopy( nx, h(1,i), hmu(1,i))
         hmu(i,i) = h(i,i) + dmu
   30 continue
      call newton( ld, nx, hmu, g, step, info)
      if ( info .gt. 0) return
      stepnm = enorm(nx, step)
      if ( (regmin .lt. stepnm) .and. (stepnm .lt. regmax) ) then
         go to 99
      endif
c
c c   update the mu bounds and get a new mu.
c
      phi = stepnm - region
      if (phi .gt. 0.0) then
         dmulow = dmu
      else
         dmuhi = dmu
      endif
c      call vecopy( nx, step, temp)
c      call dsisl( hmu, maxv, nx, ipvt, temp)
c      phip = ddot( nx, step, 1, temp, 1) / stepnm
      if (dmuhi .le. dmulow) then
         go to 99
      endif
c      if (phip .ne. 0.) then
c         dmu = dmu - (stepnm/ region) * (phi/ phip)
c      endif
       dmu = dmulow + (dmuhi - dmulow)/ 2.0d0
c       if ( (dmu .lt. dmulow) .or. (dmuhi .lt. dmu) ) then
c         dmu = dmulow + (dmuhi - dmulow)/ 2.0d0
c      endif
      go to 40
c
c c   at 99 do exit chores.
c
   99 continue
      oldmu = dmu
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine truset(nv, nq, fmu, fpred, fmerit,
     & trfrac, iter, nomin, minpt, maxpt, regold, region)
c
c  subroutine truset updates the size of the trust region.  if the
c  current fv value is less than the previous best, than the trust regio
c  doubles.  the trust region will also double if the current function
c  value is within 10% of the predicted function value.  if the fv value
c  worse than the previous worst, the trust region halves.
c
c  function values which lie between the best and worst function values
c  cause no change in the size of the trust region for n iterations.
c  thereafter, each in-between function value halves the trust region.
c
c  the reason n iterations are allowed without effect is to keep at leas
c  n function evaluations within the current trust region, to avoid the
c  problem of allowing the trust region to shrink so far that all the
c  other points besides xcen lie outside.  the latter situation is
c  undesirable because the model that we "trust" is determined by points
c  that lie outside the region in which we trust the model, a glaring
c  contradiction.
c---------------------------------------------------------------------
      include 'impli.inc'
      dimension fmu(nq)
c
      regold = region
      if ( iter .eq. (nv+1)) then
         dnv = nv
         region = trfrac * dsqrt(dnv)
         return
      else if ( fmerit .le. fmu(1) ) then
         region = 2.0d0 * regold
         nomin = 0
      else
         nomin = nomin + 1
c         if ( relerr(1, fv, fpred) .lt. 1.0d-1) then
c            region = 2.0d0 * regold
         if (nomin .ge. nv) then
            region = regold / 2.0d0
         else if ( fmerit .ge. fmu(maxpt) ) then
            region = regold / 2.0d0
         endif
      endif
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vecdif(n, v1, v2, dif)
c
c  subroutine vecdif subtracts two vectors.
c-------------------------------------------------
      include 'impli.inc'
      dimension v1(n), v2(n), dif(n)
c
      do 10 i = 1,n
         dif(i) = v1(i) - v2(i)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vecopy( n, xin, xout)
      include 'impli.inc'
      dimension xin(n), xout(n)
c
      do 10 i=1,n
         xout(i) = xin(i)
   10 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine xstop( n, xv, xnew, tol, info)
c
c  subroutine xstop determines whether current xv and the new xv lie
c  within the specified tolerance, either absolutely or relatively.
c------------------------------------------------------------------
      include 'impli.inc'
      parameter (maxv = 10, maxq = (maxv+1)*(maxv+2)/2 )
      dimension xv(n), xnew(n), xdif(maxv)
      call vecdif(n, xv, xnew, xdif)
c
c        check absolute difference
c
      absdif = enorm(n, xdif)
      if (absdif .lt. tol) then
         info = 8
         return
      endif
c
c         check for grossly oversized step
c
      xnorm = enorm(n, xv)
      xnewnm = enorm(n, xnew)
      if (xnorm .lt. tol*xnewnm) then
         info = 6
         return
      endif
c
c         check relative step size
c
      xnmax = dmax1( xnewnm, xnorm)
      reldif = absdif/ xnmax
      if (reldif .lt. tol) then
         info = 2
      endif
      return
      end
