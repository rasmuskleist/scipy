*> \brief \b DPOTRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DPOTRF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DPOTRF computes the Cholesky factorization of a real symmetric
*> positive definite matrix A.
*>
*> The factorization has the form
*>    A = U**T * U,  if UPLO = 'U', or
*>    A = L  * L**T,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the block version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization A = U**T*U or A = L*L**T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading principal minor of order i
*>                is not positive, and the factorization could not be
*>                completed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup potrf
*
*  =====================================================================
      SUBROUTINE DPOTRF(UPLO,N,P,U,LDU,VT,LDVT,D,INCD,WORK,LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INCD,LDU,LDVT,LWORK,N,P
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*), U(LDU,*), VT(LDVT,*), WORK(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL UPPER
      INTEGER I,ID,J,K,KD
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
      INTRINSIC SQRT
*     ..
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME(UPLO,'U')
      IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN
         INFO = -1
      ELSE IF (N.LT.0 .OR. P.LT.0) THEN
         INFO = -2
      ELSE IF (LDU.LT.MAX(1,N) .OR. LDVT.LT.MAX(1,P)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.MAX(1,P * P)) THEN
         INFO = 14
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DPOTRF ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0 .OR. P.EQ.0) RETURN
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
      IF (UPPER) THEN
      ELSE
         ID = KD
         DO I = 1, N
            DO J = 1, P
               DO K = 1, P
                  VT(J,I) = VT(J,I) - WORK(P * J + K) * U(I,K)
               ENDDO
            ENDDO
            DO K = 1, P
               D(ID) = D(ID) + U(I,K) * VT(I,K)
            ENDDO
            D(ID) = SQRT(D(ID))
            DO J = 1, P
               VT(J,I) = VT(J,I) / D(ID)
            ENDDO
            DO J = 1, P
               DO K = 1, P
                  WORK(P * J + K) = WORK(P * J + K) + VT(J,I) * VT(K,I)
               ENDDO
            ENDDO
            ID = ID + INCD
         ENDDO
      ENDIF
*
*     End of DPOTRF .
*
      END SUBROUTINE DPOTRF