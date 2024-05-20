*> \brief \b DTRSV
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRSV(UPLO,TRANS,N,P,U,LDU,VT,LDVT,D,INCD,X,INCX,WORK,LWORK)
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
*> DTRSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix.
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
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N )
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
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
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
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
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
      SUBROUTINE DTRSV(UPLO,TRANS,N,P,U,LDU,VT,LDVT,D,INCD,B,INCB,WORK,
     + LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INCB,INCD,INFO,LDU,LDVT,LWORK,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*), U(LDU,*), VT(LDVT,*), B(*), WORK(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,J,JB,JD,KB,KD
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
      ELSE IF (INCB.EQ.0 .OR. INCD.EQ.0) THEN
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
*     Set up the start point in B if the increment is not unity. This
*     will be  ( N - 1 ) * INCB   too small for descending loops.
*
      IF (INCB.LE.0) THEN
        KB = 1 - (N-1)*INCB
      ELSE IF (INCX.NE.1) THEN
        KB = 1
      END IF
*
*     Set up the start point in D if the increment is not unity. This
*     will be  ( N - 1 ) * INCD  too small for descending loops.
*
      IF (INCD.LE.0) THEN
        KD = 1 - (N-1)*INCD
      ELSE IF (INCX.NE.1) THEN
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
              IF (INCB.EQ.1 .AND. INCD.EQ.1) THEN
                  B(N) = B(N) / D(N)
                  DO I = 1, P
                      WORK(I) = B(N) * VT(I, N)
                  ENDDO
                  DO J = N-1, 1, -1
                      DO I = 1, P
                          B(J) = B(J) - U(J, I) * WORK(I)
                      ENDDO
                      B(J) = B(J) / D(J)
                      DO I = 1, P
                          WORK(I) = WORK(I) + B(J) * VT(I, J)
                      ENDDO
                  ENDDO
              ELSE
              END IF
          ELSE
              IF (INCB.EQ.1 .AND. INCD.EQ.1) THEN
                  B(1) = B(1) / D(1)
                  DO I = 1, P
                      WORK(I) = B(1) * VT(I, 1)
                  ENDDO
                  DO J = 2, N
                      DO I = 1, P
                          B(J) = B(J) - U(J, I) * WORK(I)
                      ENDDO
                      B(J) = B(J) / D(J)
                      DO I = 1, P
                          WORK(I) = WORK(I) + B(J) * VT(I, J)
                      ENDDO
                  ENDDO
              ELSE
              END IF
          END IF
      ELSE
*
*        Form  x := inv( A**T )*x.
*
          IF (LSAME(UPLO, 'U')) THEN
              IF (INCB.EQ.1 .AND. INCD.EQ.1) THEN
                  B(1) = B(1) / D(1)
                  DO I = 1, P
                      WORK(I) = B(1) * U(1, I)
                  ENDDO
                  DO J = 2, N
                      DO I = 1, P
                          B(J) = B(J) - VT(I, J) * WORK(I)
                      ENDDO
                      B(J) = B(J) / D(J)
                      DO I = 1, P
                          WORK(I) = WORK(I) + B(J) * U(J, I)
                      ENDDO
                  ENDDO
              ELSE
              END IF
          ELSE
              IF (INCB.EQ.1 .AND. INCD.EQ.1) THEN
                  B(N) = B(N) / D(N)
                  DO I = 1, P
                      WORK(I) = B(N) * U(N, I)
                  ENDDO
                  DO J = N - 1, 1, -1
                      DO I = 1, P
                          B(J) = B(J) - VT(I, J) * WORK(I)
                      ENDDO
                      B(J) = B(J) / D(J)
                      DO I = 1, P
                          WORK(I) = WORK(I) + B(J) * U(J, I)
                      ENDDO
                  ENDDO
              ELSE
              END IF
          END IF
      END IF
*
*     End of DTRSV .
*
      END SUBROUTINE DTRSV