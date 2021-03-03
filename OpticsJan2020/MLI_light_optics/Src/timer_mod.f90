module ml_timer
      implicit none
      integer, parameter :: ntimers=16
      integer :: ihertz
      real*8 :: hertz
      real*8, dimension(6*ntimers) :: time_vec
! each 5-vec contains, for a given operation:
! 1: start time for an operation during a step
! 2: elapsed time for an operation during a step
! 3: elasped time (minimum over PEs) during a step
! 4: elasped time (maximum over PEs) during a step
! 5: elasped time (average over PEs) during a step
! 6: elasped time (maximum over PEs) accumulated since start of run
!
      character*16 :: opname(ntimers)
      save
contains

      subroutine init_ml_timers
      use parallel
      implicit none
      integer n
      opname(1)='step'
      opname(2)='raysin'
      opname(3)='spch3d'
      opname(4)='rhoslo3d'
      opname(5)='ntrslo3d'
      opname(6)='greenf3d'
      opname(7)='fft'
      opname(8)='transp'
      opname(9)='eval'
      opname(10)='moments'
      opname(11)='ccfftnr'
      opname(12)='profile1d'
      opname(13)='timer1'
      opname(14)='timer2'
      opname(15)='timer3'
      opname(16)='timer4'
      call system_clock(count_rate=ihertz)
      hertz=ihertz
      time_vec(1:6*ntimers)=0.d0
      if(idproc.eq.0)then
      write(71,100)(opname(n)(1:10),n=1,ntimers)
      write(72,100)(opname(n)(1:10),n=1,ntimers)
      write(73,100)(opname(n)(1:10),n=1,ntimers)
  100 format(12x,20(a10,1x))
      endif
      end subroutine init_ml_timers
!
      subroutine init_step_timer
      implicit none
      integer n
      do n=1,ntimers
        time_vec(6*n-4)=0.d0
      enddo
      end subroutine init_step_timer
!
      subroutine increment_timer(str,init)
! stores the times at which various operations have been completed
      implicit none
      character(len=*) :: str
      integer :: init,iticks,n,n1,n2
! time_vec(n1) is the oldest stored time, time_vec(n2) is the most recent
      do n=1,ntimers
        n1=6*n-5
        n2=6*n-4
        if(str.ne.trim(opname(n)))cycle
          call system_clock(count=iticks)
          if(init.eq.0)then
            time_vec(n1)=iticks/hertz
          else
            time_vec(n2)=time_vec(n2)+(iticks/hertz-time_vec(n1))
          endif
        return
      enddo
      write(6,*)'timer error: string ',str,' not found.'
      end subroutine increment_timer
!
      subroutine step_timer(str,ioscreen)
! to be called at the end of each ste. computes, for all operations:
! the min/max times/proc in a step, the average time/proc in a step,
! and the total accumulated time 
      use parallel
      implicit none
      character*16 :: str
      integer ioscreen
      integer :: n,ierror
      real*8, dimension(2*ntimers) :: tmp,gtmp
!----------------
      do n=1,ntimers
!    "t_elapsed=time_vec(2)-time_vec(1)"
      tmp(2*n-1)=-time_vec(6*n-4)
      tmp(2*n  )=+time_vec(6*n-4)
      enddo
      call MPI_ALLREDUCE(tmp,gtmp,2*ntimers,mreal,mpimax,lworld,ierror)
      do n=1,ntimers
      time_vec(6*n-3)=-gtmp(2*n-1)
      time_vec(6*n-2)=+gtmp(2*n)
      enddo
!----------------
      do n=1,ntimers
      tmp(n)=time_vec(6*n-4)
      enddo
      call MPI_ALLREDUCE(tmp,gtmp,ntimers,mreal,mpisum,lworld,ierror)
      do n=1,ntimers
      time_vec(6*n-1)=gtmp(n)/nvp
      time_vec(6*n)=time_vec(6*n)+time_vec(6*n-2)
      enddo
!----------------
      if(idproc.eq.0)then
! minimum times:
      write(71,101)str(1:8),(time_vec(6*n-3),n=1,ntimers)
! maximum times:
      write(72,101)str(1:8),(time_vec(6*n-2),n=1,ntimers)
      if(ioscreen.eq.-1)write(6,101)str(1:8),(time_vec(6*n-2),n=1,ntimers)
! average times:
      write(73,101)str(1:8),(time_vec(6*n-1),n=1,ntimers)
  101 format(a8,1x,12(1pe10.3,1x))
      endif
!
      do n=1,ntimers
        time_vec(6*n-4)=0.d0
      enddo
      end subroutine step_timer
!
!
      subroutine end_ml_timers
      use parallel
      implicit none
      integer n
      if(idproc.eq.0)then
      do n=1,ntimers
        write(6,100)opname(n),time_vec(6*n)
      enddo
  100 format(a16,'::',1pe13.6)
      endif
      end subroutine end_ml_timers
end module ml_timer
