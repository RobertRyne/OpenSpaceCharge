! for calling ESSL using new FFT*HPF interface


! implements FFTs using ESSL instead of Ryne's Num.Recipes versions

      subroutine fft3dhpf (N1, N2, N3, KSign, Scale, ICpy, NAdj, x, x3)
      implicit none
!Constants
!Args
      integer N1, N2, N3, KSign, ICpy, NAdj
      real*8 Scale
      complex*16 x,x3
      dimension x(N1,N2,N3)
      dimension x3(N3,N2,N1)
!Locals
      complex*16 aux  !just a placeholder, not really used
      integer i,j,k

!ValidateArgs

      if( KSign .EQ. 0 )then
          write(6,*) 'ERROR: FFT3DHPF/ESSL with zero ksign'
          stop 'FFT3err1'
      endif
      if( NAdj .NE. 0 )then
          write(6,*) 'ERROR: FFT3DHPF/ESSL with non-zero nadj'
          stop 'FFT3err2'
      endif
      if( ICpy .NE. 0 )then
          write(6,*) 'ERROR: FFT3DHPF/ESSL with non-zero icpy'
          stop 'FFT3err3'
      endif

      if( Scale .EQ. 0.0 )then
          write(6,*) 'ERROR: FFT3DHPF/ESSL with zero scale'
          stop 'FFT3err4'
      endif
      
!Execute
      ! the ESSL routine doesn't leave the output transposed, so 
      ! do the FFT in-place and transpose into the output variable.

      call DCFT3( x, N1,N1*N2, x, N1,N1*N2, N1,N2,N3
     &           ,KSign,Scale ,aux,0 )

      do k = 1 ,N3
          do j = 1,N2
              do i = 1,N1
                  x3(k,j,i) = x(i,j,k)
              enddo
          enddo
      enddo

!Done

      return
      end

!----------------------------------------------------------------

      subroutine fft2dhpf()
      stop 'FFT2Derr'
      return
      end

