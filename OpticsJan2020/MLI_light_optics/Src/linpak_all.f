      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  DGECO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination
C            and estimates the condition of the matrix.
C***DESCRIPTION
C
C     DGECO factors a double precision matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGEFA is slightly faster.
C     To solve  A*X = B , follow DGECO by DGESL.
C     To compute  INVERSE(A)*C , follow DGECO by DGESL.
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
C     To compute  INVERSE(A) , follow DGECO by DGEDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an INTEGER vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
C***END PROLOGUE  DGECO
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
C***BEGIN PROLOGUE  DGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
C            using the factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)
C***BEGIN PROLOGUE  DSICO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2B1A
C***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX,SYMMETRIC
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a d.p. SYMMETRIC matrix by elimination with symmet-
C            ric pivoting and estimates the condition of the matrix.
C***DESCRIPTION
C
C     DSICO factors a double precision symmetric matrix by elimination
C     with symmetric pivoting and estimates the condition of the
C     matrix.
C
C     If  RCOND  is not needed, DSIFA is slightly faster.
C     To solve  A*X = B , follow DSICO by DSISL.
C     To compute  INVERSE(A)*C , follow DSICO by DSISL.
C     To compute  INVERSE(A) , follow DSICO by DSIDI.
C     To compute  DETERMINANT(A) , follow DSICO by DSIDI.
C     To compute  INERTIA(A), follow DSICO by DSIDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the symmetric matrix to be factored.
C                Only the diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     Output
C
C        A       a block diagonal matrix and the multipliers which
C                were used to obtain it.
C                The factorization can be written  A = U*D*TRANS(U)
C                where  U  is a product of permutation and unit
C                upper triangular matrices, TRANS(U) is the
C                transpose of  U , and  D  is block diagonal
C                with 1 by 1 and 2 by 2 blocks.
C
C        KPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     LINPACK DSIFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,IABS,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DSCAL,DSIFA
C***END PROLOGUE  DSICO
      INTEGER LDA,N,KPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
      DOUBLE PRECISION ANORM,S,DASUM,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
C***FIRST EXECUTABLE STATEMENT  DSICO
      DO 30 J = 1, N
         Z(J) = DASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DSIFA(A,LDA,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = 1.0D0
      DO 50 J = 1, N
         Z(J) = 0.0D0
   50 CONTINUE
      K = N
   60 IF (K .EQ. 0) GO TO 120
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 90
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 110
  100    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
      GO TO 60
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     &         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE U*D*V = Y
C
      K = N
  170 IF (K .EQ. 0) GO TO 230
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 200
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 220
  210    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
      GO TO 170
  230 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     &         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DSISL(A,LDA,N,KPVT,B)
C***BEGIN PROLOGUE  DSISL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2B1A
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE,
C             SYMMETRIC
C***AUTHOR  BUNCH, J., (UCSD)
C***PURPOSE  Solves the double precision SYMMETRIC system
C            A*X=B  using the factors computed by DSIFA.
C***DESCRIPTION
C
C     DSISL solves the double precision symmetric system
C     A * X = B
C     using the factors computed by DSIFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA,N)
C                the output from DSIFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        KPVT    INTEGER(N)
C                the pivot vector from DSIFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero may occur if  DSICO  has set RCOND .EQ. 0.0
C        or  DSIFA  has set INFO .NE. 0  .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DSIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C     Fortran IABS
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DSISL
      INTEGER LDA,N,KPVT(1)
      DOUBLE PRECISION A(LDA,1),B(1)
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
C***FIRST EXECUTABLE STATEMENT  DSISL
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
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
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
****BEGIN PROLOGUE  DGEFA
****DATE WRITTEN   780814   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D2A1
****KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
****AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
****PURPOSE  Factors a double precision matrix by Gaussian elimination.
****DESCRIPTION
*
*     DGEFA factors a double precision matrix by Gaussian elimination.
*
*     DGEFA is usually called by DGECO, but it can be called
*     directly with a saving in time if  RCOND  is not needed.
*     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
*
*     On Entry
*
*        A       DOUBLE PRECISION(LDA, N)
*                the matrix to be factored.
*
*        LDA     INTEGER
*                the leading dimension of the array  A .
*
*        N       INTEGER
*                the order of the matrix  A .
*
*     On Return
*
*        A       an upper triangular matrix and the multipliers
*                which were used to obtain it.
*                The factorization can be written  A = L*U  where
*                L  is a product of permutation and unit lower
*                triangular matrices and  U  is upper triangular.
*
*        IPVT    INTEGER(N)
*                an integer vector of pivot indices.
*
*        INFO    INTEGER
*                = 0  normal value.
*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
*                     condition for this subroutine, but it does
*                     indicate that DGESL or DGEDI will divide by zero
*                     if called.  Use  RCOND  in DGECO for a reliable
*                     indication of singularity.
*
*     LINPACK.  This version dated 08/14/78 .
*     Cleve Moler, University of New Mexico, Argonne National Lab.
*
*     Subroutines and Functions
*
*     BLAS DAXPY,DSCAL,IDAMAX
****REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
*                 *LINPACK USERS  GUIDE*, SIAM, 1979.
****ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
****END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
*
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
*
*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
*
****FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
*
*        FIND L = PIVOT INDEX
*
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
*
*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
*
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
*
*           INTERCHANGE IF NECESSARY
*
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
*
*           COMPUTE MULTIPLIERS
*
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
*
*           ROW ELIMINATION WITH COLUMN INDEXING
*
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)
****BEGIN PROLOGUE  DSIFA
****DATE WRITTEN   780814   (YYMMDD)
****REVISION DATE  820801   (YYMMDD)
****CATEGORY NO.  D2B1A
****KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,
*             SYMMETRIC
****AUTHOR  BUNCH, J., (UCSD)
****PURPOSE  Factors a d.p. SYMMETRIC matrix by elimination with
*            symmetric pivoting
****DESCRIPTION
*
*     DSIFA factors a double precision symmetric matrix by elimination
*     with symmetric pivoting.
*
*     To solve  A*X = B , follow DSIFA by DSISL.
*     To compute  INVERSE(A)*C , follow DSIFA by DSISL.
*     To compute  DETERMINANT(A) , follow DSIFA by DSIDI.
*     To compute  INERTIA(A) , follow DSIFA by DSIDI.
*     To compute  INVERSE(A) , follow DSIFA by DSIDI.
*
*     On Entry
*
*        A       DOUBLE PRECISION(LDA,N)
*                the symmetric matrix to be factored.
*                Only the diagonal and upper triangle are used.
*
*        LDA     INTEGER
*                the leading dimension of the array  A .
*
*        N       INTEGER
*                the order of the matrix  A .
*
*     On Return
*
*        A       a block diagonal matrix and the multipliers which
*                were used to obtain it.
*                The factorization can be written  A = U*D*TRANS(U)
*                where  U  is a product of permutation and unit
*                upper triangular matrices, TRANS(U) is the
*                transpose of  U , and  D  is block diagonal
*                with 1 by 1 and 2 by 2 blocks.
*
*        KPVT    INTEGER(N)
*                an integer vector of pivot indices.
*
*        INFO    INTEGER
*                = 0  normal value.
*                = K  if the K-th pivot block is singular.  This is
*                     not an error condition for this subroutine,
*                     but it does indicate that DSISL or DSIDI may
*                     divide by zero if called.
*
*     LINPACK.  This version dated 08/14/78 .
*     James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
*
*     Subroutines and Functions
*
*     BLAS DAXPY,DSWAP,IDAMAX
*     Fortran DABS,DMAX1,DSQRT
****REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
*                 *LINPACK USERS  GUIDE*, SIAM, 1979.
****ROUTINES CALLED  DAXPY,DSWAP,IDAMAX
****END PROLOGUE  DSIFA
      INTEGER LDA,N,KPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
*
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,IDAMAX
      LOGICAL SWAP
*
*     INITIALIZE
*
*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
****FIRST EXECUTABLE STATEMENT  DSIFA
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
*
      INFO = 0
*
*     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
*
      K = N
   10 CONTINUE
*
*        LEAVE THE LOOP IF K=0 OR K=1.
*
*     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0D0) INFO = 1
*     ......EXIT
            GO TO 200
   20    CONTINUE
*
*        THIS SECTION OF CODE DETERMINES THE KIND OF
*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
*        REQUIRED.
*
         KM1 = K - 1
         ABSAKK = DABS(A(K,K))
*
*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
*        COLUMN K.
*
         IMAX = IDAMAX(K-1,A(1,K),1)
         COLMAX = DABS(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
*
*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
*           ROW IMAX.
*
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))
   50       CONTINUE
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
*
*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
*
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
*
*           1 X 1 PIVOT BLOCK.
*
            IF (.NOT.SWAP) GO TO 120
*
*              PERFORM AN INTERCHANGE.
*
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
*
*           PERFORM THE ELIMINATION.
*
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
*
*           SET THE PIVOT ARRAY.
*
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
*
*           2 X 2 PIVOT BLOCK.
*
            IF (.NOT.SWAP) GO TO 160
*
*              PERFORM AN INTERCHANGE.
*
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
*
*           PERFORM THE ELIMINATION.
*
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0D0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
*
*           SET THE PIVOT ARRAY.
*
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
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
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      implicit double precision (a-h, o-z)	
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
C     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
      DOUBLE PRECISION DX(1),DY(1)
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
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
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
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
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      implicit double precision (a-h, o-z)	
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C           C.L.LAWSON, 1978 JAN 08
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C                        PHASE 1.  SUM IS ZERO
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C                                PREPARE FOR PHASE 4.
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C                  PREPARE FOR PHASE 3.
   75 SUM = (SUM * XMAX) * XMAX
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
   85 HITEST = CUTHI/FLOAT( N )
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C              END OF MAIN LOOP.
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
