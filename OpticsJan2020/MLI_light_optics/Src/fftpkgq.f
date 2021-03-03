! 3D FFT package
      subroutine fft3dhpf (n1, n2, n3, ksign, scale, icpy, nadj, x, x3)
      use parallel  !for idproc
      use ml_timer
      implicit none
      integer n1, n2, n3, ksign, icpy, nadj, icpy2d, i,j,k
      real*8 scale,scale2d,ezero
      integer klen,k1,kn1,kn3,ierror
      complex*16 x,x2,x3
      complex*16 x2t,x3t
      dimension x(n1,n2,n3)
      dimension x3(n3,n2,n1)
      dimension x3t(n3,n2,n1)
      dimension x2(n2,n1,n3)
      real*8 alnum,rem
      dimension x2t(n2,n1,n3)

      call increment_timer('fft',0)

!!XXX <dbs> Feb 2003 -- moved here from afro.f, because it applies only to this version of FFT 
!---Jan 27, 2003
! make sure that these are powers of 2
      alnum=log(1.d0*n1)/log(2.d0)
      rem=abs(alnum-nint(alnum))
      if(rem.gt.1.d-8)then
        if(idproc.eq.0)then
        write(6,*)'(fft3dhpf) error: n1 has an input value of ',n1
        write(6,*)'           but must be a power of 2. Halting.'
        endif
        call myexit
      endif
      alnum=log(1.d0*n2)/log(2.d0)
      rem=abs(alnum-nint(alnum))
      if(rem.gt.1.d-8)then
        if(idproc.eq.0)then
        write(6,*)'(fft3dhpf) error: n2 has an input value of ',n2
        write(6,*)'           but must be a power of 2. Halting.'
        endif
        call myexit
      endif
      if(n3.ne.0)then
        alnum=log(1.d0*n3)/log(2.d0)
        rem=abs(alnum-nint(alnum))
        if(rem.gt.1.d-8)then
          if(idproc.eq.0)then
            write(6,*)'(fft3dhpf) error: n3 has an input value of ',n3
            write(6,*)'           but must be a power of 2. Halting.'
          endif
          call myexit
        endif
      endif

!---
!     write(6,*)'starting fft'
! FFTs along all dimensions except the last:
!
      scale2d=1.0
      icpy2d=0
!
      klen=n3/nvp
!logic goes here to decide whether to do serial or parallel fft.
!for example:
!     if(klen.le.1000000)then   !always do serial
!     if(klen.le.256)then !do parallel if every proc has a lot of work
      if(klen.eq.0)then   !do parallel if every proc has some work to do
        call increment_timer('timer1',0)
        do k=1,n3
        call fft2dhpf(n1,n2,ksign,nadj,scale2d,icpy2d,                  &
     &                x(1,1,k),x2(1,1,k))
        enddo
        call increment_timer('timer1',1)
      else
! note to myself: check to see if this will work if n3, nvp are not
! both powers of 2
        k1=idproc*klen+1
        kn3=k1+klen-1
        call increment_timer('timer1',0)
        do k=k1,kn3
        call fft2dhpf(n1,n2,ksign,nadj,scale2d,icpy2d,                  &
     &                x(1,1,k),x2t(1,1,k))
        enddo
        call increment_timer('timer1',1)
!
        call increment_timer('timer2',0)
        call MPI_ALLGATHER(x2t(1,1,k1),n1*n2*klen,mcplx,x2(1,1,1),      &
     &  n1*n2*klen,mcplx,lworld,ierror)
        call increment_timer('timer2',1)
      endif
!
! The 3D and higher, the previous call transposes the subarray of all but
! the last index; Now move the last index from the right-most position to
! the left-most position in preparation for final, in-place transformation:
! In 2D this is simply a tranpose of the entire array.
!     write(6,*)'starting triple do to move x2 into x3'
      call increment_timer('transp',0)
!2D   x3=transpose(x)
            do i = 1, n1
      do k = 1, n3
         do j = 1, n2
               x3(k,j,i) = x2(j,i,k)
            end do
         end do
      end do
!     x3(1:n3,1:n2,1:n1)=x2(1:n2,1:n1,1:n3)
!     write(6,*)'finished triple-do'
      call increment_timer('transp',1)
!
! Final transform along the left-most direction
!
      klen=n1/nvp
      if(klen.eq.0)then
        call increment_timer('timer3',0)
        call mfft_local13d(n3,n2,n1,ksign,nadj,1,n1,x3)
        call increment_timer('timer3',1)
      else
        k1=idproc*klen+1
        kn1=k1+klen-1
        call increment_timer('timer3',0)
        call mfft_local13d (n3,n2,n1,ksign,nadj,k1,kn1,x3)
        call increment_timer('timer3',1)
        x3t(:,:,k1:kn1)=x3(:,:,k1:kn1)
        call increment_timer('timer4',0)
        call MPI_ALLGATHER(x3t(1,1,k1),n3*n2*klen,mcplx,x3(1,1,1),      &
     &  n3*n2*klen,mcplx,lworld,ierror)
        call increment_timer('timer4',1)
      endif
!
!
!
      if(scale.ne.1.0)then
        do k=1,n1
        do j=1,n2
        do i=1,n3
!       x3(:,:,:)=x3(:,:,:)*scale
        x3(i,j,k)=x3(i,j,k)*scale
        enddo
        enddo
        enddo
      endif
! if icpy.eq.1, store the final result back in x
! if icpy.ne.1, store with zeros (ensure users know it's filled with junk)
      if(icpy.eq.1)then
      write(6,*)'error: should not get here! (icpy=1)'
!2D     x=transpose(x3)
        do k = 1, n3
           do j = 1, n2
              do i = 1, n1
                 x(i,j,k) = x3(k,j,i)
              end do
           end do
        end do
      else
        ezero=0.
        do k = 1, n3
           do j = 1, n2
              do i = 1, n1
                 x(i,j,k) = cmplx(ezero,ezero)
              end do
           end do
        end do
!       write(6,*)'...finished fft (without copying)'
      endif
      call increment_timer('fft',1)
      return
      end !subroutine fft3dhpf

!=======================================================================
!
! HPF_Local subroutine to perform multiple FFTs on the inner dimension
!
      subroutine mfft_local13d(n1,n2,n3,ksign,nadj,k1,kn3,x)
      implicit none
      integer n1,n2,n3,ksign,nadj,i,j,k,n3top,n2top
      integer k1,kn3
      complex*16 x,tempi
      dimension x(n1,n2,n3),tempi(n1)
!
! Perform multiple FFTs without scaling:
      n3top=n3
      n2top=n2
c this works for open bc's because it is not necessary to do the
c inverse transform over the complete (octupled) domain.
c Of the N^2 1D fft's performed over the ixj plane, only 1/4 are
c in the physical domain. However, for the periodic case, the
c full z-dimension has to be covered. At the time of the inverse
c transformation (when this routine is called with ksign=-1),
c the n2 and n3 dimensions are in fact y and z, respectively.
c So the full fft is needed along the third dimension in that case.
      if(ksign.eq.-1 .and. nadj.eq.0)then
!dec28        n3top=n3/2
!dec28        n2top=n2/2
      endif
      if(ksign.eq.-1 .and. nadj.eq.1)then
!       write(12,*)'8-29-01:testing needed for this part of fftpkgq.f'
!dec28      n2top=n2/2
      endif
cryne Nov 9, 2003      do k=1,n3top
      do k=k1,kn3
      do j=1,n2top
cryne        tempi(1:n1)=x(1:n1,j,k)
cryne        call ccfftnr(tempi,n1,ksign)
cryne        x(1:n1,j,k) = tempi(1:n1)
      call ccfftnr(x(1,j,k),n1,ksign)
      end do
      end do
      return
      end
!     end subroutine mfft_local13d
!=======================================================================
!=======================================================================
! 2D FFT package
      subroutine fft2dhpf (n1,n2,ksign,nadj,scale,icpy,x,x3)
         use parallel, only : idproc
         implicit none
         integer n1, n2, ksign, nadj, icpy, i,j,ihalf,itop
         real*8 scale,ezero
         complex*16 x,x3
         dimension x(n1,n2),x3(n2,n1)
!
! FFTs along all dimensions except the last:
! Note: In 3D and higher, this is accomplished with a routine mfft_local2
! But in 2D this isn't necessary--just used mfft_local1
!
         ihalf=0
         call mfft_local1 (ksign,x,n1,n2,ihalf)
!
! The 3D and higher, the previous call transposes the subarray of all but
! the last index; Now move the last index from the right-most position to
! the left-most position in preparation for final, in-place transformation:
! In 2D this is simply a tranpose of the entire array. 
cryne    x3=transpose(x)
         itop=n1
!dec28         if(ksign.eq.-1.and.nadj.eq.0)itop=n1/2
         do i=1,itop
         do j=1,n2
         x3(j,i)=x(i,j)
         enddo
         enddo
!
! Final transform along the left-most direction
!
         if(ksign.eq.-1.and.nadj.eq.0)ihalf=1
!        if(ihalf.eq.1)write(6,*)'ihalf=1'
         call mfft_local1 (ksign,x3,n2,n1,ihalf)
!
!        if(scale.ne.1.0)x3=x3*scale
         if(scale.ne.1.0)then
         do j=1,n1
         do i=1,n2
         x3(i,j)=x3(i,j)*scale
         enddo
         enddo
         endif
! if icpy.eq.1, store the final result back in x
! if icpy.ne.1, store with zeros (ensure users know it's filled with junk)
         if(icpy.eq.1)then
           write(6,*)'error:should not get here!(icpy.eq.1 in fft2dhpf)'
!          x=transpose(x3)
           do j=1,n2
           do i=1,n1
           x(i,j)=x3(j,i)
           enddo
           enddo
         endif
!        if(icpy.ne.1)x=(0.,0.)
         ezero=0.
         if(icpy.ne.1)then
           do j=1,n2
           do i=1,n1
           x(i,j)=cmplx(ezero,ezero)
           enddo
           enddo
         endif
         return
      end !subroutine fft2dhpf
!
!=======================================================================
!
! HPF_Local subroutine to perform multiple FFTs on the inner dimension
!
      subroutine mfft_local1 (ksign,x,m1,m2,ihalf)
         implicit none
         integer ksign,j,m1,m2,m2top,ihalf
         complex*16 x,tempi
         dimension x(m1,m2),tempi(m1)
!        complex tempo,table,work
!        dimension tempo(m1),table(m1),work(2*m1)
!
! Initialize:
!!!!!        call ccfft (0,m1,1.0,tempi,tempo,table,work,0)
! Perform multiple FFTs without scaling:
         m2top=m2
!dec28         if(ihalf.eq.1)m2top=m2/2
         do j = 1, m2top
c            tempi = x(:,j)
!!!!!           call ccfft(ksign,m1,1.0,tempi,tempo,table,work,0)
!!!!!           x(:,j) = tempo
cryne        call ccfftnr(tempi,m1,ksign)
            call ccfftnr(x(1,j),m1,ksign)
c            x(:,j) = tempi
         end do
         return
      end
!     end subroutine mfft_local1
!
      subroutine ccfftnr(cdata,nn,isign)
      use ml_timer
      implicit real*8(a-h,o-z)
      complex*16 cdata
      dimension cdata(nn),data(2*nn)
      call increment_timer('ccfftnr',0)
      do i=1,nn
        data(2*i-1)=real(cdata(i))
        data(2*i) =aimag(cdata(i))
      enddo
! bit reversal:
      n=2*nn
      j=1
      do i=1,n,2
       if(j.gt.i)then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
       endif
       m=n/2
    1  if((m.ge.2).and.(j.gt.m))then
        j=j-m
        m=m/2
        goto 1
       endif
       j=j+m
      enddo
! Danielson-Lanczos:
      twopi=4.0*asin(1.0d0)
      mmax=2
    2 if(n.gt.mmax)then
       istep=2*mmax
       theta=twopi/(isign*mmax)
       wpr=-2.*sin(0.5*theta)**2
       wpi=sin(theta)
       wr=1.0
       wi=0.0
       do m=1,mmax,2
        do i=m,n,istep
         j=i+mmax
         tempr=wr*data(j)-wi*data(j+1)
         tempi=wr*data(j+1)+wi*data(j)
         data(j)=data(i)-tempr
         data(j+1)=data(i+1)-tempi
         data(i)=data(i)+tempr
         data(i+1)=data(i+1)+tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
       enddo
       mmax=istep
       goto 2
      endif
c     ezero=0.
c     eunit=1.
      do i=1,nn
        cdata(i)=data(2*i-1)+(0.,1.)*data(2*i)
c       cdata(i)=data(2*i-1)+cmplx(ezero,eunit)*data(2*i)
      enddo
      call increment_timer('ccfftnr',1)
      return
      end
!=======================================================================
