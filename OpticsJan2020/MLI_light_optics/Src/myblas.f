      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     &                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     &   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END

      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     &   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     &       INTA.GE.145 .AND. INTA.LE.153 .OR.
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
*
**************************************************************************
****BEGIN PROLOGUE  DCOPY
****PURPOSE  Copy a vector.
****LIBRARY   SLATEC (BLAS)
****CATEGORY  D1A5
****TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
****KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
****AUTHOR  Lawson, C. L., (JPL)
*           Hanson, R. J., (SNLA)
*           Kincaid, D. R., (U. of Texas)
*           Krogh, F. T., (JPL)
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*       DY  double precision vector with N elements
*     INCY  storage spacing between elements of DY
*
*     --Output--
*       DY  copy of vector DX (unchanged if N .LE. 0)
*
*     Copy double precision DX to double precision DY.
*     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
*     defined in a similar way using INCY.
*
****REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
*                 Krogh, Basic linear algebra subprograms for Fortran
*                 usage, Algorithm No. 539, Transactions on Mathematical
*                 Software 5, 3 (September 1979), pp. 308-323.
****ROUTINES CALLED  (NONE)
****REVISION HISTORY  (YYMMDD)
*   791001  DATE WRITTEN
*   890831  Modified array declarations.  (WRB)
*   890831  REVISION DATE from Version 3.2
*   891214  Prologue converted to Version 4.0 format.  (BAB)
*   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
*   920501  Reformatted the REFERENCES section.  (WRB)
****END PROLOGUE  DCOPY

      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end

**************************************************************************

      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
****BEGIN PROLOGUE  DSWAP
****DATE WRITTEN   791001   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D1A5
****KEYWORDS  BLAS,DOUBLE PRECISION,INTERCHANGE,LINEAR ALGEBRA,VECTOR
****AUTHOR  LAWSON, C. L., (JPL)
*           HANSON, R. J., (SNLA)
*           KINCAID, D. R., (U. OF TEXAS)
*           KROGH, F. T., (JPL)
****PURPOSE  Interchange d.p. vectors
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*       DY  double precision vector with N elements
*     INCY  storage spacing between elements of DY
*
*     --Output--
*       DX  input vector DY (unchanged if N .LE. 0)
*       DY  input vector DX (unchanged if N .LE. 0)
*
*     Interchange double precision DX and double precision DY.
*     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
*     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
*     defined in a similar way using INCY.
****REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
****ROUTINES CALLED  (NONE)
****END PROLOGUE  DSWAP
*
      DOUBLE PRECISION DX(1),DY(1),DTEMP1,DTEMP2,DTEMP3
****FIRST EXECUTABLE STATEMENT  DSWAP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
*
*       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
*
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       CODE FOR BOTH INCREMENTS EQUAL TO 1
*
*
*       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
*
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
*
*     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
*
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
****BEGIN PROLOGUE  IDAMAX
****DATE WRITTEN   791001   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D1A2
****KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAXIMUM COMPONENT,
*             VECTOR
****AUTHOR  LAWSON, C. L., (JPL)
*           HANSON, R. J., (SNLA)
*           KINCAID, D. R., (U. OF TEXAS)
*           KROGH, F. T., (JPL)
****PURPOSE  Find largest component of d.p. vector
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*
*     --Output--
*   IDAMAX  smallest index (zero if N .LE. 0)
*
*     Find smallest index of maximum magnitude of double precision DX.
*     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
****REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
****ROUTINES CALLED  (NONE)
****END PROLOGUE  IDAMAX
*
      DOUBLE PRECISION DX(1),DMAX,XMAG
****FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
*
*        CODE FOR INCREMENTS NOT EQUAL TO 1.
*
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
*
*        CODE FOR INCREMENTS EQUAL TO 1.
*
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
****************************************************************
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
****BEGIN PROLOGUE  DASUM
****DATE WRITTEN   791001   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D1A3A
****KEYWORDS  ADD,BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAGNITUDE,SUM,
*             VECTOR
****AUTHOR  LAWSON, C. L., (JPL)
*           HANSON, R. J., (SNLA)
*           KINCAID, D. R., (U. OF TEXAS)
*           KROGH, F. T., (JPL)
****PURPOSE  Sum of magnitudes of d.p. vector components
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*
*     --Output--
*    DASUM  double precision result (zero if N .LE. 0)
*
*     Returns sum of magnitudes of double precision DX.
*     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
****REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
****ROUTINES CALLED  (NONE)
****END PROLOGUE  DASUM
*
      DOUBLE PRECISION DX(1)
****FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
*
*        CODE FOR INCREMENTS NOT EQUAL TO 1.
*
      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN
*
*        CODE FOR INCREMENTS EQUAL TO 1.
*
*
*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
*
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     &   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
****************************************************************
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
****BEGIN PROLOGUE  DAXPY
****DATE WRITTEN   791001   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D1A7
****KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,TRIAD,VECTOR
****AUTHOR  LAWSON, C. L., (JPL)
*           HANSON, R. J., (SNLA)
*           KINCAID, D. R., (U. OF TEXAS)
*           KROGH, F. T., (JPL)
****PURPOSE  D.P computation y = a*x + y
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DA  double precision scalar multiplier
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*       DY  double precision vector with N elements
*     INCY  storage spacing between elements of DY
*
*     --Output--
*       DY  double precision result (unchanged if N .LE. 0)
*
*     Overwrite double precision DY with double precision DA*DX + DY.
*     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
*       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
*       and LY is defined in a similar way using INCY.
****REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
****ROUTINES CALLED  (NONE)
****END PROLOGUE  DAXPY
*
      DOUBLE PRECISION DX(1),DY(1),DA
****FIRST EXECUTABLE STATEMENT  DAXPY
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
*
*        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
*
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        CODE FOR BOTH INCREMENTS EQUAL TO 1
*
*
*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
*
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
*
*        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
*
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
****************************************************************
      SUBROUTINE DSCAL(N,DA,DX,INCX)
****BEGIN PROLOGUE  DSCAL
****DATE WRITTEN   791001   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D1A6
****KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
****AUTHOR  LAWSON, C. L., (JPL)
*           HANSON, R. J., (SNLA)
*           KINCAID, D. R., (U. OF TEXAS)
*           KROGH, F. T., (JPL)
****PURPOSE  D.P. vector scale x = a*x
****DESCRIPTION
*
*                B L A S  Subprogram
*    Description of Parameters
*
*     --Input--
*        N  number of elements in input vector(s)
*       DA  double precision scale factor
*       DX  double precision vector with N elements
*     INCX  storage spacing between elements of DX
*
*     --Output--
*       DX  double precision result (unchanged if N.LE.0)
*
*     Replace double precision DX by double precision DA*DX.
*     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
****REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
*                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
*                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
*                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
****ROUTINES CALLED  (NONE)
****END PROLOGUE  DSCAL
*
      DOUBLE PRECISION DA,DX(1)
****FIRST EXECUTABLE STATEMENT  DSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
*
*        CODE FOR INCREMENTS NOT EQUAL TO 1.
*
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
*
*        CODE FOR INCREMENTS EQUAL TO 1.
*
*
*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
*
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
****************************************************************
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      implicit double precision (a-h, o-z)	
*     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
*     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)
*     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
*     DEFINED IN A SIMILAR WAY USING INCY.
      DOUBLE PRECISION DX(1),DY(1)
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
*         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*        CODE FOR BOTH INCREMENTS EQUAL TO 1.
*        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     &   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
*         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
