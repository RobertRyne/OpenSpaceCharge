      blockdata mpifbd
      logical           mpif_initialized
      common /MPILCOMN/ mpif_initialized
      data              mpif_initialized/.FALSE./
      end

      subroutine MPI_SEND( )
      include 'mpif.h'
      print*,'in MPI_SEND: not implemented'
      stop 'MPISEND'
      end

      subroutine MPI_RECV( )
      include 'mpif.h'
      print*,'in MPI_RECV: not implemented'
      stop 'MPRECV'
      end

      subroutine MPI_GET_COUNT( )
      include 'mpif.h'
      print*,'in MPI_GET_COUNT: not implemented'
      stop 'MPGETC'
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine MPI_ALLGATHER( Sendbuf,SendCount,SendType
     &                         ,recvbuf,RecvCount,RecvType
     &                         ,Comm,ierror )
      include 'mpif.h'
      integer Sendbuf(*) ,SendCount,SendType
     &       ,recvbuf(*) ,RecvCount,RecvType ,Comm ,ierror
      ! ignore possible differences in type and count
      call  MPI_ALLREDUCE( Sendbuf,recvbuf,RecvCount,RecvType
     &                    ,MPI_SUM,Comm,ierror )
      if( ierror .NE. 0 )then
        print*,'error: MPI_ALLGATHER: error in MPI_ALLREDUCE()'
      else
        ! check for count and type differences
        if( SendCount .NE. RecvCount )then
          ierror = 99
        elseif( SendType .NE. RecvType )then
          ierror = 98
        endif
      endif
      return
      end


      subroutine MPI_ALLREDUCE( Sendbuf,recvbuf,Count,Datatype,Op,Comm
     &                         ,ierror )
      include 'mpif.h'
      integer Sendbuf(*) ,recvbuf(*) ,Count,Datatype,Op,Comm ,ierror
      integer dtype ,dtypemult

      dtype = MOD( Datatype,100 ) / 10
      dtypemult = MOD( Datatype ,10 )

!XXX      if( Op .LT. MPII_OP_FIRST .OR. Op .GT. MPII_OP_LAST )then
!XXX        print*,'error: MPI_ALLREDUCE: unknown reduce operator ',Op
!XXX        print*,'info: op must be in ',MPII_OP_FIRST,' to ',MPII_OP_LAST
!XXX        ierror = 1
!XXX        return
!XXX      endif
      if( dtype .LT. 0 .OR. dtype .GT. MPII_MAX_TYPE )then
        print*,'error: MPI_ALLREDUCE: unknown datatype ',Datatype
        ierror = 2
      endif
      if( dtypemult .LE. 0 .OR. dtypemult .GT. MPII_MAX_TYPE_MULT )then
        print*,'error: MPI_ALLREDUCE: unknown datatype ',Datatype
        ierror = 3
        return
      endif

      call MPII_COPY( recvbuf,Sendbuf,count,dtypemult )
      ierror = 0
      ! warn about invalid args
      if( Op .LT. MPII_OP_FIRST .OR. Op .GT. MPII_OP_LAST )then
        ierror = 99
      elseif( Comm .NE. MPI_COMM_WORLD )then
        ierror = 98
      endif
      return
      end

      subroutine MPI_BCAST( Bdata,DataLen,DataType,Root,Comm,ierror )
      integer Bdata ,DataLen,DataType,Root,Comm,ierror
      ierror = 0
      return
      end


      subroutine MPI_INITIALIZED( flag ,ierror )
      include 'mpif.h'
      integer flag ,ierror
      if( mpif_initialized )then
        flag = 1
      else
        flag = 0
      endif
      ierror = 0
      return
      end

      subroutine MPI_INIT( ierror )
      include 'mpif.h'
      integer ierror
      mpif_initialized = .TRUE.
      ierror = 0
      return
      end

      subroutine MPI_COMM_SIZE( Communicator ,num_p ,status )
      include 'mpif.h'
      integer Communicator ,num_p ,status(MPI_STATUS_SIZE)
      num_p = 1
      status(1) = 0
      if( Communicator .NE. MPI_COMM_WORLD )then
        status(1) = 99
      endif
      return
      end

      subroutine MPI_COMM_RANK( Communicator ,rank_p ,ierror )
      include 'mpif.h'
      integer Communicator ,rank_p ,ierror
      rank_p = 0
      ierror = 0
      if( Communicator .NE. MPI_COMM_WORLD )then
        ierror = 99
      endif
      return
      end

      subroutine MPI_BARRIER( Communicator ,ierror )
      include 'mpif.h'
      integer Communicator ,ierror
      ierror = 0
      if( Communicator .NE. MPI_COMM_WORLD )then
        ierror = 99
      endif
      return
      end

      subroutine MPI_FINALIZE( ierror )
      include 'mpif.h'
      integer ierror
      mpif_initialized = .FALSE.
      ierror = 0
      return
      end

      subroutine MPII_COPY( recvbuf,Sendbuf,Count,Dtypemult )
      integer  recvbuf(*),Sendbuf(*) ,Count,Dtypemult
      integer i
      do i = 1,Count*Dtypemult
        recvbuf(i) = Sendbuf(i)
      enddo
      return
      end
