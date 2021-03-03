! -*- Mode: Fortran; Modified: "Fri 14 Nov 2003 17:34:29 by dbs"; -*- 

! Computes 3D sine transform on a cubic grid using FFT from IBM ESSL

      subroutine sfft3d(f,N)
      implicit none
!Arguments
      integer N
      real*8 f(0:N,0:N,0:N)
!Locals
      integer i,j,k,ll ,naux1,naux2
      real*8 in1d(0:N-1)
      real*8 sincoef(0:N-1)
      real*8 out1d(0:N)
      real*8 aux1(22000+N+N),aux2(25000)
      real*8 Pi

!Initialization
      Pi = ATAN(1.0d0) * 4
      naux1 = 22000+N+N
      naux2 = 25000
!Execution

      ! N must be an even number.
      if ((N/2)*2 .ne. N) then
        write(6,*) "sfft3d: N is odd, N = ",N
        call abort
      endif

      ! Set up preliminary stuff.
      do ll = 0,N-1
        sincoef(ll) = sin(ll*Pi/N)
      enddo

      ! assuming cubic grids, initialize internal data just once
      call DRCFT( 1 ,in1d ,0 ,out1d,0 ,N,1 ,+1,1.0d0
     &    ,aux1,naux1 ,aux2,naux2 )

      ! Sine transform in the x direction.
      do k = 0,N-1
        do j = 0,N-1
          in1d(0) = 0.d0
          do i = 1,N-1
            in1d(i) = sincoef(i)*(f(i,j,k) + f(N-i,j,k)) +
     &          (f(i,j,k) - f(N-i,j,k))/2
          enddo

          call DRCFT( 0 ,in1d ,0 ,out1d,0 ,N,1 ,+1,1.0d0
     &        ,aux1,naux1 ,aux2,naux2 )
          
          f(1,j,k) = out1d(0)/2
          f(0,j,k) = 0.d0

          ! ESSL returns the complex result, but
          ! this is the same as the real result just
          ! reordered so the imaginary part of out(1:N/2-1)
          ! is the same as the real part of out(N/2:N-1)
          do i = 1,N/2-1
            f(2*i,j,k) = -out1d(2*i+1) !out(N - i)
            f(2*i+1,j,k) = out1d(2*i) + f(2*i-1,j,k) !out(i)
          enddo

        enddo
      enddo


      ! Sine transform in the y direction.
      do k = 0,N-1
        do j = 0,N-1
          in1d(0) = 0.d0
          do i = 1,N-1
            in1d(i) = sincoef(i)*(f(j,i,k) + f(j,N-i,k)) +
     &          (f(j,i,k) - f(j,N-i,k))/2
          enddo

          call DRCFT( 0 ,in1d ,0 ,out1d,0 ,N,1 ,+1,1.0d0
     &        ,aux1,naux1 ,aux2,naux2 )
          
          f(j,1,k) = out1d(0)/2
          f(j,0,k) = 0.d0
          
          ![NOTE: see comment above.]
          do i = 1,N/2-1
            f(j,2*i,k) = -out1d(2*i+1) !out(N - i)
            f(j,2*i+1,k) = out1d(2*i) + f(j,2*i-1,k) !out(i)
          enddo

        enddo
      enddo


      ! Sine transform in the z direction.
      do k = 0,N-1
        do j = 0,N-1
          in1D(0) = 0.d0

          do i = 1,N-1
            in1D(i) = sincoef(i)*(f(j,k,i) + f(j,k,N-i)) +
     &          (f(j,k,i) - f(j,k,N-i))/2
          enddo

          call DRCFT( 0 ,in1d ,0 ,out1d,0 ,N,1 ,+1,1.0d0
     &        ,aux1,naux1 ,aux2,naux2 )
          
          f(j,k,0) = 0.d0
          f(j,k,1) = out1d(0)/2
          
          ![NOTE: see comment above.]
          do i = 1,N/2-1
            f(j,k,2*i) = -out1d(2*i+1) !out(N - i)
            f(j,k,2*i+1) = out1d(2*i) + f(j,k,2*i-1) !out(i)
          enddo

        enddo
      enddo

!done
      return
      end
