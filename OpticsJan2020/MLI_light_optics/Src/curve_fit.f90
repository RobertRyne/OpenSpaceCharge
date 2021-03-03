!***********************************************************************
!
! curve_fit: module containing functions for fitting curves to data
!
! Description:
!     This module implements subroutines for fitting curves to data.
!
! Version: 0.1
! Author: D.T.Abell, Tech-X Corp., Apr.2005
!
! Comments
!   
!
!***********************************************************************
!
      module curve_fit
        use parallel, only : idproc
        implicit none
!
! functions and subroutines
!
      contains
!
!***********************************************************************
      subroutine cf_polyInterp(xx,yy,x,y,ix,iorder,iorigin)
!
!  Interpolate at abcissa x the ordinate value y determined by the
!  polynomial of degree iorder that fits the (iorder+1) data points
!  (xx,yy) nearest x.  The interpolation uses Neville's algorithm, and
!  iorder has a default value of three.
!  NB: The values in xx MUST increase (or decrease) monotonically.
      implicit none
      double precision, dimension(:), intent(in) :: xx,yy  ! input data
      double precision, intent(in) :: x  ! interpolate at this abcissa
      double precision, intent(out) :: y  ! interpolated ordinate value
      integer, intent(inout) :: ix  ! index of xx nearest x
      integer, optional, intent(in) :: iorder  ! polynomial order
      integer, optional, intent(in) :: iorigin  ! index origin of xx,yy
!-----!----------------------------------------------------------------!
      integer :: ioff,iord
      integer :: i,ii,il,m
      double precision, dimension(:), allocatable :: p

      ! use iorder to set iord
      iord=3
      if(present(iorder)) iord=iorder
      ! use index origin to set index offset
      ioff=0
      if(present(iorigin)) ioff=iorigin-1

      ! allocate temporary array p(0:iord)
      allocate(p(0:iord))

      ! determine lower boundary of interpolation range [il,il+iord]
      ! (within this subroutine, xx and yy have index origin one)
      ix=ix-ioff
      call cf_searchNrst(xx,x,ix)
      il=ix-iord/2
      if (il.lt.1) then
        il=1
      else if (il.gt.size(xx)-iord) then
        il=size(xx)-iord
      end if

      ! use array value or do interpolation using Neville's algorithm
      if (x.eq.xx(ix)) then
        y=yy(ix)
      else
        do i=0,iord
          p(i)=yy(il+i)
        end do
        do m=1,iord
          do i=0,iord-m
            ii=il+i
            p(i)=((x-xx(ii))*p(i+1)+(xx(ii+m)-x)*p(i))/(xx(ii+m)-xx(ii))
          end do
        end do
        y=p(0)
      end if

      ! offset index for external value of index origin
      ix=ix+ioff

      ! deallocate temporary array p(0:iord)
      deallocate(p)

      return
      end subroutine cf_polyInterp
!
!***********************************************************************
      subroutine cf_polyInterp3(xx,yy,x,y,ix,iorigin)
!
!  Interpolate at abcissa x the ordinate value y determined by the
!  third-order polynomial that fits the four data points (xx,yy)
!  nearest x.  The interpolation uses Neville's algorithm.
!  NB: The values in xx MUST increase (or decrease) monotonically.
      implicit none
      double precision, dimension(:), intent(in) :: xx,yy  ! input data
      double precision, intent(in) :: x  ! interpolate at this abcissa
      double precision, intent(out) :: y  ! interpolated ordinate value
      integer, intent(inout) :: ix  ! index of xx nearest x
      integer, optional, intent(in) :: iorigin  ! index origin of xx,yy
!-----!----------------------------------------------------------------!
      integer :: ioff,iord
      integer :: i,ii,il,m
      double precision, dimension(0:3) :: p

      ! use cubic interpolation
      iord=3
      ! use index origin to set index offset
      ioff=0
      if(present(iorigin)) ioff=iorigin-1

      ! determine lower boundary of interpolation range [il,il+iord]
      ! (within this subroutine, xx and yy have index origin one)
      ix=ix-ioff
      call cf_search(xx,x,ix)
      il=ix-1  ! il=ix-iorder/2=ix-3/2
      if (il.lt.1) then
        il=1
      else if (il.gt.size(xx)-3) then
        il=size(xx)-3
      end if

      ! use array value or do interpolation using Neville's algorithm
      if (x.eq.xx(ix)) then
        y=yy(ix)
      else
        do i=0,iord
          p(i)=yy(il+i)
        end do
        !write(21,*) x,xx(il:il+iord)
        !write(21,*) x,p(0:iord)
        do m=1,iord
          do i=0,iord-m
            ii=il+i
            p(i)=((x-xx(ii))*p(i+1)+(xx(ii+m)-x)*p(i))/(xx(ii+m)-xx(ii))
          end do
          !write(21,*) x,p(0:iord-m)
        end do
        y=p(0)
      end if

      ! offset index for external value of index origin
      ix=ix+ioff

      return
      end subroutine cf_polyInterp3
!
!***********************************************************************
      subroutine cf_bivarInterp(ss,tt,ff,s,t,f,is,it,sorder,torder,     &
     &                          iorigin)
!
!  The two-dimensional array ff records a set of data values over the
!  points in the direct product of the one-dimensional arrays ss and tt.
!  This subroutine performs bivariate polynomial interpolation of order
!  sorder in s and torder in t to interpolate a value f = f(s,t) using
!  the (sorder+1)*(torder+1) points nearest (s,t) and the corresponding
!  data in ff.  The interpolation uses Neville's algorithm in each
!  variable, and sorder and torder both have default value three.  The
!  arrays are assumed to have an index origin of iorigin, which has a
!  default value of one.
!  NB: The values in ss and tt MUST increase (or decrease) monotonically.
      implicit none
      double precision, dimension(:), intent(in) :: ss,tt  ! arg arrays
      double precision, dimension(:,:), intent(in) :: ff  ! ff(ss,tt)
      double precision, intent(in) :: s,t  ! interpolate at this point
      double precision, intent(out) :: f  ! interpolated value
      integer, intent(inout) :: is,it  ! indices of ss, tt nearest s,t
      integer, optional, intent(in) :: sorder,torder  ! polynomial orders
      integer, optional, intent(in) :: iorigin  ! index origin of arrays
!-----!----------------------------------------------------------------!
      integer :: ioff,sord,tord
      integer :: i,ii,isl
      double precision, dimension(:), allocatable :: p,v

      ! use sorder and torder to set sord and tord
      sord=3; tord=3
      if(present(sorder)) sord=sorder
      if(present(torder)) tord=torder

      ! use index origin to set index offset
      ioff=0
      if(present(iorigin)) ioff=iorigin-1
      ! adjust is and it for index origin (on entrance)
      ! (within this subroutine, ss tt, and ff have index origin one)
      is=is-ioff; it=it-ioff

      ! determine lower boundary of s interpolation range [isl,isl+sord]
      call cf_search(ss,s,is)
      isl=is-sord/2
      if (isl.lt.1) then
        isl=1
      else if (isl.gt.size(ss)-sord) then
        isl=size(ss)-sord
      end if

      ! allocate temporary arrays p(0:sord) and v(0:sord)
      allocate(p(0:sord))
      allocate(v(0:sord))

      ! do polynomial interpolation in t for nearby s values
      do i = 0,sord
        v(i)=ss(isl+i)
        call cf_polyInterp(tt,ff(isl+i,:),t,p(i),it,tord)
      end do
      ! do polynomial interpolation in s
      ii=is-isl
      call cf_polyInterp(v,p,s,f,ii,sord,0)

      ! adjust is and it for index origin (on exit)
      is=is+ioff; it=it+ioff

      ! deallocate temporary arrays p and v
      deallocate(p)
      deallocate(v)

      return
      end subroutine cf_bivarInterp
!
!***********************************************************************
      function cf_interp(xx,a,b,c,x)
!
!  Interpolate at x the value of a function described by the n-element
!  coefficient arrays a, b, and c.  These arrays hold the quadratic
!  fit parameters determined at locations xx by subroutine cf_parfit
!  such that f = a + b*x + c*x**2.
!  NB: The values in xx MUST increase (or decrease) monotonically.
      implicit none
      double precision :: cf_interp
      double precision, dimension(:) :: a,b,c,xx
      double precision :: x
!-----!----------------------------------------------------------------!
      integer :: j,n

      n=size(xx)
      call cf_search(xx,x,j)
      if(j.lt.1) then
        j=1
      else if(j.gt.(n-1))  then
        j=n-1
      endif
      cf_interp=a(j)+x*(b(j)+x*c(j))
      return
      end function cf_interp
!
!***********************************************************************
      subroutine cf_interp2(xx,yy,x,y)
!
!  Quadratic interpolation at x based on data pairs (xx,yy).
!  NB: The data in xx MUST increase (or decrease) monotonically.
      implicit none
      double precision, dimension(0:2), intent(in) :: xx,yy
      double precision, intent(in) :: x
      double precision, intent(out) :: y
!-----!----------------------------------------------------------------!
      double precision :: c1,c2

      c1=(yy(1)-yy(0))/(xx(1)-xx(0))
      c2=((yy(2)-yy(0))/(xx(2)-xx(0))-c1)/(xx(2)-xx(1))
      y=yy(0)+(x-xx(0))*(c1+(x-xx(1))*c2)
      return
      end subroutine cf_interp2
!
!***********************************************************************
      subroutine cf_parfit(xi,yi,ai,bi,ci)
!
!  This subroutine takes data pairs (x(i),y(i)), i=1..npts, fits a
!  quadratic form a+b*x+c*x**2 to each triple of adjacent points, and
!  returns the lists of coefficients a(i), b(i), and c(i).
!  Except at the ends, the coefficients reported for the i-th interval
!  are computed as the average of those calculated using the points
!  (i-1,i,i+1) and those calculated using the points (i,i+1,i+2).
      implicit none
      double precision, dimension(:), intent(in) :: xi,yi
      double precision, dimension(:), intent(out) :: ai,bi,ci
!-----!----------------------------------------------------------------!
      integer :: i,ip2,npts,npm1,npm2
      double precision :: x1,x2,x3,y1,y2,y3
      double precision :: d1,d2,d3,h1,h2,h3
      double precision :: a1,b1,c1,a2,b2,c2
!
      npts=size(xi)
      npm1=npts-1
      npm2=npts-2
!
! initialize the fitting
      x1=xi(1)
      x2=xi(2)
      x3=xi(3)
      y1=yi(1)
      y2=yi(2)
      y3=yi(3)
      h1=x2-x1
      h2=x3-x2
      h3=x3-x1
      d1= y1/(h3*h1)
      d2=-y2/(h1*h2)
      d3= y3/(h2*h3)
! and fit a parabola to the first three points
      a1=  d1*x2*x3   + d2*x3*x1   + d3*x1*x2
      b1= -d1*(x2+x3) - d2*(x3+x1) - d3*(x1+x2)
      c1=  d1         + d2         + d3
      ai(1)=a1
      bi(1)=b1
      ci(1)=c1
!
! compute the next set of quadratic coefficients,
! and average with the previous set
      do i=2,npm2
        ip2=i+2
        x1=x2
        x2=x3
        x3=xi(ip2)
        y1=y2
        y2=y3
        y3=yi(ip2)
        h1=h2
        h2=x3-x2
        h3=x3-x1
        d1= y1/(h1*h3)
        d2=-y2/(h1*h2)
        d3= y3/(h2*h3)
        a2=  d1*x2*x3   + d2*x3*x1   + d3*x1*x2
        b2= -d1*(x2+x3) - d2*(x3+x1) - d3*(x1+x2)
        c2=  d1         + d2         + d3
        ai(i)=0.5d0*(a1+a2)
        bi(i)=0.5d0*(b1+b2)
        ci(i)=0.5d0*(c1+c2)
        a1=a2
        b1=b2
        c1=c2
      end do
!
! rightmost interval
      ai(npm1)=a2
      bi(npm1)=b2
      ci(npm1)=c2
!
      return
      end subroutine cf_parfit
!
!***********************************************************************
      subroutine cf_search(xx,x,j,iorigin)
!
!  Search array xx(1:n) and return index j such that x lies in range
!  [xx(j),x(j+1)) or, if j=n-1, [xx(n-1),xx(n)].  If the value returned
!  in j is either 0 or n=size(xx), then x lies outside the range of xx.
!  Use the value of j on input as an initial guess for searching.
!  If the array submitted to this routine has index origin other than 1,
!  then one may specify it using the optional argument iorigin.  In this
!  case, of course, the values of j that indicate an out-of-range
!  condition are (iorigin-1) and (size(xx)+iorigin-1).
!  NB: The array xx must be monotonic, either increasing or decreasing.
!  Cf. Numerical Recipes, 2nd ed. (1992), p.121.
      implicit none
      double precision, dimension(:), intent(in) :: xx
      double precision, intent(in) :: x
      integer, intent(inout) :: j
      integer, optional, intent(in) :: iorigin
!-----!----------------------------------------------------------------!
      integer jl,jm,ju,incr,n
      logical ascend
      double precision, parameter :: epsilon=1.d-13

      ! get size of xx, and adjust for index origin
      n=size(xx)
      if(present(iorigin)) j=j-(iorigin-1)
      ! do array values ascend or descend
      ascend=(xx(n).gt.xx(1))
      ! use initial value of j to bracket binary search
      if (j.lt.1.or.j.gt.n) then
        ! bracket entire array
        jl=0; ju=n+1
      else
        ! hunt for smaller bracket
        jl=j; ju=j
        incr=1
        if ((x.gt.xx(jl).eqv.ascend).or.x.eq.xx(jl)) then
          ! move bracket -=>
 1        ju=jl+incr
          if (ju.gt.n) then
            ju=n+1
          else if ((x.gt.xx(ju).eqv.ascend).or.x.eq.xx(ju)) then
            jl=ju
            incr=incr+incr
            goto 1
          end if
        else
          ! move bracket <=-
 2        jl=ju-incr
          if (jl.lt.1) then
            jl=0
          else if ((x.lt.xx(jl).eqv.ascend).and.x.ne.xx(jl)) then
            ju=jl
            incr=incr+incr
            goto 2
          end if
        end if
      end if
      ! now do binary search
 3    if (ju-jl.gt.1) then
        jm=(jl+ju)/2
        if(((x.gt.xx(jm)).eqv.ascend).or.(x.eq.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        goto 3
      end if
      j=jl
      ! be nice if within epsilon of edges
      if (j.eq.0) then
        if (abs(x-xx(1)).le.epsilon) j=1
      else if (j.eq.n) then
        if (abs(x-xx(n)).le.epsilon) j=n-1
      end if
      ! adjust for index origin
      if(present(iorigin)) j=j+(iorigin-1)
!
      return
      end subroutine cf_search
!
!***********************************************************************
      subroutine cf_searchNrst(xx,x,j,iorigin)
!
!  Search array xx(1:n) and return index j such that x lies in range
!  [xx(j),x(j+1)) or, if j=n-1, [xx(n-1),xx(n)].  If the value returned
!  in j is either 0 or n=size(xx), then x lies outside the range of xx.
!  Use the value of j on input as an initial guess for searching.
!  If the array submitted to this routine has index origin other than 1,
!  then one may specify it using the optional argument iorigin.  In this
!  case, of course, the values of j that indicate an out-of-range
!  condition are (iorigin-1) and (size(xx)+iorigin-1).
!  NB: The array xx must be monotonic, either increasing or decreasing.
!  Cf. Numerical Recipes, 2nd ed. (1992), p.121.
      implicit none
      double precision, dimension(:), intent(in) :: xx
      double precision, intent(in) :: x
      integer, intent(inout) :: j
      integer, optional, intent(in) :: iorigin
!-----!----------------------------------------------------------------!
      integer jl,jm,ju,incr,n
      logical ascend
      double precision, parameter :: epsilon=1.d-13

      ! get size of xx, and adjust for index origin
      n=size(xx)
      if(present(iorigin)) j=j-(iorigin-1)
      ! do array values ascend or descend
      ascend=(xx(n).gt.xx(1))
      ! use initial value of j to bracket binary search
      if (j.lt.1.or.j.gt.n) then
        ! bracket entire array
        jl=0; ju=n+1
      else
        ! hunt for smaller bracket
        jl=j; ju=j
        incr=1
        if ((x.gt.xx(jl).eqv.ascend).or.x.eq.xx(jl)) then
          ! move bracket -=>
 1        ju=jl+incr
          if (ju.gt.n) then
            ju=n+1
          else if ((x.gt.xx(ju).eqv.ascend).or.x.eq.xx(ju)) then
            jl=ju
            incr=incr+incr
            goto 1
          end if
        else
          ! move bracket <=-
 2        jl=ju-incr
          if (jl.lt.1) then
            jl=0
          else if ((x.lt.xx(jl).eqv.ascend).and.x.ne.xx(jl)) then
            ju=jl
            incr=incr+incr
            goto 2
          end if
        end if
      end if
      ! now do binary search
 3    if (ju-jl.gt.1) then
        jm=(jl+ju)/2
        if(((x.gt.xx(jm)).eqv.ascend).or.(x.eq.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        goto 3
      end if
      j=jl
      ! be nice if within epsilon of edges,
      ! and be sure we return Nearest index
      if (j.eq.0) then
        if (abs(x-xx(1)).le.epsilon) j=1
      else if (j.eq.n) then
        if (abs(x-xx(n)).le.epsilon) j=n
      else
        if (abs(x-xx(j)).gt.abs(x-xx(j+1))) j=j+1
      end if
      ! adjust for index origin
      if(present(iorigin)) j=j+(iorigin-1)
!
      return
      end subroutine cf_searchNrst
!
!***********************************************************************
      subroutine cf_locate(xx,x,j,iorigin)
!
!  Search array xx(:) and return index j such that x lies in range
!  [xx(j),x(j+1)) or, if j=n-1, [xx(n-1),xx(n)].  If the value returned
!  in j is either 0 or n=size(xx), then x lies outside the range of xx.
!  If the array submitted to this routine has index origin other than 1,
!  then one may specify it using the optional argument iorigin.  In this
!  case, of course, the values of j that indicate an out-of-range
!  condition are (iorigin-1) and (size(xx)+iorigin-1).
!  NB: The array xx must be monotonic, either increasing or decreasing.
!  Cf. Numerical Recipes, 2nd ed. (1992), p.121.
      implicit none
      double precision, dimension(:), intent(in) :: xx
      double precision, intent(in) :: x
      integer, intent(inout) :: j
      integer, optional, intent(in) :: iorigin
!-----!----------------------------------------------------------------!
      integer jl,jm,ju,incr,n
      logical ascend
      double precision, parameter :: epsilon=1.d-13

      ! get size of xx
      n=size(xx)
      ! do array values ascend or descend?
      ascend=(xx(n).gt.xx(1))
      ! now initialize and perform binary search
      jl=0; ju=n+1
 3    if (ju-jl.gt.1) then
        jm=(jl+ju)/2
        if(((x.gt.xx(jm)).eqv.ascend).or.(x.eq.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        goto 3
      endif
      j=jl
      ! be nice if within epsilon of edges
      if (j.eq.0) then
        if (abs(x-xx(1)).le.epsilon) j=1
      else if (j.eq.n) then
        if (abs(x-xx(n)).le.epsilon) j=n-1
      end if
      ! adjust for index origin
      if(present(iorigin)) j=j+(iorigin-1)
!
      return
      end subroutine cf_locate
!
!***********************************************************************
      subroutine cf_locate_nr(xx,x,j,iorigin)
!
!  Search array xx(:) and return index j such that x lies in range
!  [xx(j),x(j+1)).  If the value j returned is either 0 or n=size(xx),
!  then x lies outside the range of xx.
!  If the array submitted to this routine has index origin other than 1,
!  then one may specify it using the optional argument iorigin.  In this
!  case, of course, the values of j that indicate an out-of-range
!  condition are (iorigin-1) and (size(xx)+iorigin-1).
!  NB: The array xx must be monotonic, either increasing or decreasing.
!  Cf. Numerical Recipes, 2nd ed. (1992), p.121.
      implicit none
      double precision, dimension(:), intent(in) :: xx
      double precision, intent(in) :: x
      integer, intent(inout) :: j
      integer, optional, intent(in) :: iorigin
!-----!----------------------------------------------------------------!
      integer jl,jm,ju,incr,n
      logical ascend
      double precision, parameter :: eps=1.d-13

      ! get size of xx
      n=size(xx)
      ! do array values ascend or descend?
      ascend=(xx(n).gt.xx(1))
      ! now initialize and perform binary search
      jl=0; ju=n+1
 3    if (ju-jl.gt.1) then
        jm=(jl+ju)/2
        if(((x.gt.xx(jm)).eqv.ascend).or.(abs(x-xx(jm)).lt.eps)) then
          jl=jm
        else
          ju=jm
        endif
        goto 3
      endif
      !if (jl.eq.n.and.x.eq.xx(n)) jl=n-1
      j=jl
      ! adjust for index origin
      if(present(iorigin)) j=j+(iorigin-1)
!
      return
      end subroutine cf_locate_nr
!
      end module curve_fit

