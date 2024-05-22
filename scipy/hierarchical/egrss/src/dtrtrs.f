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
      SUBROUTINE DTRTRS(UPLO,TRANS,N,P,NRHS,U,LDU,VT,LDVT,D,INCD,B,LDB,
     $                  WORK,LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INCD,INFO,LDB,LDU,LDVT,LWORK,N,NRHS,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*), U(LDU,*), VT(LDVT,*), B(*), WORK(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D+0 )
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Executable Statements ..
*     ..
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
         INFO = -1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. 
     + .NOT.LSAME(TRANS,'T') .AND.
     + .NOT.LSAME(TRANS,'C')) THEN
         INFO = -2
      ELSE IF (N.LT.0 .OR. P.LT.0) THEN
         INFO = -4
      ELSE IF (NRHS.LT.0) THEN
         INFO = -5
      ELSE IF (LDU.LT.MAX(1,N) .OR. LDVT.LT.MAX(1,P)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -9
      ELSE IF (LWORK.LT.MAX(1,P)) THEN
         INFO = -10
      END IF

      IF (INFO.NE.0) THEN
         CALL XERBLA('DTRSV ', -INFO)
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF (NHRS.EQ.0 .OR. N.EQ.0 .OR. P.EQ.0) RETURN
*
*     Check for singularity.
*
      DO INFO = 1, N
         IF (D(INFO).EQ.ZERO) RETURN
      END DO
*
*     Solve A * x = b  or  A**T * x = b.
*
      CALL DTRSM(UPLO,TRANS,N,NRHS,P,U,LDU,VT,LDVT,D,INCD,B,LDB,
     $ WORK, LWORK)
*
*     End of DTRTRS .
*
      END SUBROUTINE DTRTRS