      implicit none

      integer    MPI_STATUS_SIZE   ,MPI_COMM_WORLD
      parameter( MPI_STATUS_SIZE=1 ,MPI_COMM_WORLD=2 )
      ! values of type parameters are 100*<count>+10*<type>+<mult>,
      ! where <count> is the number of words in the type,
      !       <type> indicates the base type,
      !   and <mult> is the multiple of the length of this type times
      !              int length
      integer    MPI_DOUBLE_PRECISION     ,MPI_DOUBLE_COMPLEX
      parameter( MPI_DOUBLE_PRECISION=112 ,MPI_DOUBLE_COMPLEX=124 )
      integer    MPI_2DOUBLE_PRECISION
      parameter( MPI_2DOUBLE_PRECISION=214 )
      integer    MPI_INTEGER     ,MPI_2INTEGER
      parameter( MPI_INTEGER=101 ,MPI_2INTEGER=202 )
      ![NOTE: have to update these when types change!!! -dbs]
      integer    MPII_TYPE_INT   ,MPII_TYPE_REAL   ,MPII_TYPE_CPLX
      parameter( MPII_TYPE_INT=0 ,MPII_TYPE_REAL=1 ,MPII_TYPE_CPLX=2 )
      integer    MPII_MAX_TYPE   ,MPII_MAX_TYPE_MULT
      parameter( MPII_MAX_TYPE=2 ,MPII_MAX_TYPE_MULT=4 )

      integer    MPI_SUM      ,MPI_MAX      ,MPI_MAXLOC
      parameter( MPI_SUM=1001 ,MPI_MAX=1002 ,MPI_MAXLOC=1003 )
      integer    MPII_OP_FIRST         ,MPII_OP_LAST
      parameter( MPII_OP_FIRST=MPI_SUM ,MPII_OP_LAST=MPI_MAXLOC )

      logical           mpif_initialized
      common /MPILCOMN/ mpif_initialized
      save   /MPILCOMN/
