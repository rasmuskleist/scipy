*> \brief \b DTRMV
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRMV(UPLO,TRANS,N,P,U,LDU,VT,LDVT,D,INCD,X,INCX,WORK,LWORK)
*
*       .. Scalar Arguments ..
*       INTEGER INCD,INCX,LDU,LDVT,N
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
*> DTRMV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix on extended generator form.
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
*>              TRANS = 'N' or 'n'   A*x = b.
*>
*>              TRANS = 'T' or 't'   A**T*x = b.
*>
*>              TRANS = 'C' or 'c'   A**T*x = b.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>           On entry, P specifies the order of the matrix A.
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
*> \param[in,out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element right-hand side vector b. On exit, X is overwritten
*>           with the solution vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*
*  Authors:
*  ========
*
*> \author Rasmus Kleist Hørlyck Sørensen
*
*> \date May 2024
*
*> \ingroup double_blas_level1
*
*  =====================================================================
      SUBROUTINE DTRMV(UPLO,TRANS,N,P,U,LDU,VT,LDVT,D,INCD,X,INCX,WORK,
     + LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCD,INFO,LDU,LDVT,LWORK,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*), U(LDU,*), VT(LDVT,*), X(*), WORK(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,J,JX,JD,KX,KD
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
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     + .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (N.LT.0 .OR. P.LT.0) THEN
          INFO = 4
      ELSE IF (LDU.LT.MAX(1,N) .OR. LDVT.LT.MAX(1,P)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0 .OR. INCD.EQ.0) THEN
          INFO = 8
      ELSE IF (LWORK.LT.MAX(1,P)) THEN
          INFO = 12
      END IF

      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0 .OR. P.EQ.0) RETURN
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 ) * INCX   too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE
          KX = 1
      END IF
*
*     Set up the start point in D if the increment is not unity. This
*     will be  ( N - 1 ) * INCD  too small for descending loops.
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
*        Form  y := A*x.
*
          IF (LSAME(UPLO, 'U')) THEN
              JX = KX + (N-1)*INCX
              JD = KD + (N-1)*INCD
              DO I = 1, P
                  WORK(I) = X(JX) * VT(I, N)
              ENDDO
              X(JX) = X(JX) * D(JD)
              DO J = N - 1, 1, -1
                  JX = JX - INCX
                  JD = JD - INCD
                  TEMP = X(JX) * D(JD)
                  DO I = 1, P
                      TEMP = TEMP + WORK(I) * U(J, I)
                  ENDDO
                  DO I = 1, P
                      WORK(I) = WORK(I) + X(JX) * VT(I, J)
                  ENDDO
                  X(JX) = TEMP
              ENDDO
          ELSE
              JX = KX
              JD = KD
              DO I = 1, P
                  WORK(I) = X(JX) * VT(I, 1)
              ENDDO
              X(JX) = X(JX) * D(JD)
              DO J = 2, N
                  JX = JX + INCX
                  JD = JD + INCD
                  TEMP = X(JX) * D(JD)
                  DO I = 1, P
                      TEMP = TEMP + WORK(I) * U(J, I)
                  ENDDO
                  DO I = 1, P
                      WORK(I) = WORK(I) + X(JX) * VT(I, J)
                  ENDDO
                  X(JX) = TEMP
              ENDDO
          END IF
      ELSE
*
*        Form  y := ( A**T )*x.
*
          IF (LSAME(UPLO, 'U')) THEN
              JX = KX
              JD = KD
              DO I = 1, P
                  WORK(I) = X(JX) * U(1, I)
              ENDDO
              X(JX) = X(JX) * D(JD)
              DO J = 2, N
                  JX = JX + INCX
                  JD = JD + INCD
                  TEMP = X(JX) * D(JD)
                  DO I = 1, P
                      TEMP = TEMP + WORK(I) * VT(I, J)
                  ENDDO
                  DO I = 1, P
                      WORK(I) = WORK(I) + X(JX) * U(J, I)
                  ENDDO
                  X(JX) = TEMP
              ENDDO
          ELSE
              JX = KX + (N-1)*INCX
              JD = KD + (N-1)*INCD
              DO I = 1, P
                  WORK(I) = X(JX) * U(N, I)
              ENDDO
              X(JX) = X(JX) * D(JD)
              DO J = N - 1, 1, -1
                  JX = JX - INCX
                  JD = JD - INCD
                  TEMP = X(JX) * D(JD)
                  DO I = 1, P
                      TEMP = TEMP + WORK(I) * VT(I, J)
                  ENDDO
                  DO I = 1, P
                      WORK(I) = WORK(I) + X(JX) * U(J, I)
                  ENDDO
                  X(JX) = TEMP
              ENDDO
          END IF
      END IF
*
*     End of DTRMV .
*
      END SUBROUTINE DTRMV