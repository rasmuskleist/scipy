*> \brief \b DTRSV
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRSM(UPLO,TRANS,M,P,U,LDU,VT,LDVT,D,IDD,B,LDB,WORK,LWORK)
*
*       .. Scalar Arguments ..
*       INTEGER INCD,INCX,LDU,LDVT,M
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION D(*),U(LDU,*),VT(LDVT,*),X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix on extended generator form
*>
*>    A = tril(U*V**T,-1) + diag(D) or A = triu(U*V**T,1) + diag(D),
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the equations to be solved as
*>           follows:
*>
*>              TRANS = 'M' or 'n'   A*x = b.
*>
*>              TRANS = 'T' or 't'   A**T*x = b.
*>
*>              TRANS = 'C' or 'c'   A**T*x = b.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>           On entry, P specifies the semiseperability rank of the matrix A.
*>           P must be at least zero.
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is DOUBLE PRECISION array, dimension ( LDU, P )
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>           On entry, LDU specifies the first dimension of U as declared
*>           in the calling (sub) program. LDU must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          U is DOUBLE PRECISION array, dimension ( LDV, M )
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDVT
*> \verbatim
*>          LDVT is INTEGER
*>           On entry, LDVT specifies the first dimension of VT as declared
*>           in the calling (sub) program. LDVT must be at least
*>           max( 1, p ).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array D must contain the n
*>           element right-hand side vector b. On exit, X is overwritten
*>           with the solution vector x.
*> \endverbatim
*>
*> \param[in] INCD
*> \verbatim
*>          INCD is INTEGER
*>           On entry, INCD specifies the increment for the elements of
*>           D. INCD must not be zero.
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension ( LDB, N )
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in,out] WORK
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element right-hand side vector b. On exit, X is overwritten
*>           with the solution vector x.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Rasmus Kleist Hørlyck Sørensen
*
*> \date May 2024
*
*> \ingroup double_egrss_level3
*
*  =====================================================================
      SUBROUTINE DTRMM(UPLO,TRANS,M,N,P,U,LDU,VT,LDVT,D,INCD,B,LDB,WORK,
     + LWORK)
*
*     .. Scalar Arguments ..
      INTEGER INCD,LDB,LDU,LDVT,LWORK,M,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION B(LDB,*), D(*), U(LDU,*), VT(LDVT,*), WORK(*)
*     ..
*ª
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,ID,J,K,KD
      DOUBLE PRECISION TEMP
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
         INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. 
     + .NOT.LSAME(TRANS,'T') .AND.
     + .NOT.LSAME(TRANS,'C')) THEN
         INFO = 2
      ELSE IF (LWORK.LT.MAX(1,P)) THEN
         INFO = 14
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0 .OR. P.EQ.0) RETURN
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCD.LE.0) THEN
          KD = 1 - (N-1)*INCD
      ELSE
          KD = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS, 'N')) THEN
*
*        Form  x := inv( A )*x.
*
         IF (LSAME(UPLO, 'U')) THEN
            DO J = 1, N
               ID = KD + (M-1)*INCD
               DO K = 1, P
                  WORK(K) = B(M,J) * VT(K,M)
               ENDDO
               B(M,J) = B(M,J) * D(ID)
               DO I = M-1, 1, -1
                  ID = ID - INCD
                  TEMP = B(I,J) * D(ID)
                  DO K = 1, P
                     TEMP = TEMP + U(I,K) * WORK(K)
                  ENDDO
                  DO K = 1, P
                     WORK(K) = WORK(K) + B(I,J) * VT(K,I)
                  ENDDO
                  B(I,J) = TEMP
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO K = 1, P
                  WORK(K) = B(1,J) * VT(K,1)
               ENDDO
               ID = KD
               B(1,J) = B(1,J) * D(ID)
               DO I = 2, M
                  ID = ID + INCD
                  TEMP = B(I,J) * D(ID)
                  DO K = 1, P
                     TEMP = TEMP + U(I,K) * WORK(K)
                  ENDDO
                  DO K = 1, P
                     WORK(K) = WORK(K) + B(I,J) * VT(K,I)
                  ENDDO
                  B(I,J) = TEMP
               ENDDO
            ENDDO
         END IF
      ELSE
*
*        Form  x := inv( A**T )*x.
*
         IF (LSAME(UPLO, 'U')) THEN
            DO J = 1, N
               DO K = 1, P
                  WORK(K) = B(1,J) * U(1,K)
               ENDDO
               ID = KD
               B(1,J) = B(1,J) * D(ID)
               DO I = 2, M
                  ID = ID + INCD
                  TEMP = B(I,J) * D(ID)
                  DO K = 1, P
                     TEMP = TEMP + VT(K,I) * WORK(K)
                  ENDDO
                  DO K = 1, P
                     WORK(K) = WORK(K) + B(I,J) * U(I,K)
                  ENDDO
                  B(I,J) = TEMP
               ENDDO
            ENDDO
         ELSE
            DO J = 1, N
               DO K = 1, P
                  WORK(K) = B(M,J) * U(M,K)
               ENDDO
               ID = KD + (M-1)*INCD
               B(M,J) = B(M,J) * D(ID)
               DO I = M-1, 1, -1
                  ID = ID - INCD
                  TEMP = B(I,J) * D(ID)
                  DO K = 1, P
                     TEMP = TEMP + VT(K,I) * WORK(K)
                  ENDDO
                  DO K = 1, P
                     WORK(K) = WORK(K) + B(I,J) * U(I,K)
                  ENDDO
                  B(I,J) = TEMP
               ENDDO
            ENDDO
         END IF
      END IF
*
*     End of DTRMM .
*
      END SUBROUTINE DTRMM