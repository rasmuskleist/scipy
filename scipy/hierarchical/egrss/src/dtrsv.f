      SUBROUTINE DTRSV(UPLO,TRANS,N,P,U,VT,D,B,WORK,LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INFO,LWORK,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(N), U(N,P), VT(P,N), B(N), WORK(LWORK)
*     ..
*
*  =====================================================================
*     ..
*     .. External Functions ..
      EXTERNAL DDOT
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY
      EXTERNAL DCOPY
      EXTERNAL XERBLA
*     ..

*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     + .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 4
      END IF

      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN

*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      INFO = 0
      CALL DCOPY(P, 0.0D0, 0, WORK, 1)
      IF (LSAME(TRANS, 'N')) THEN
          IF (LSAME(UPLO, 'U')) THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(WORK, U(J,:))) / D(J)
                  CALL DAXPY(P, B(J), VT(1, J), 1, WORK, 1)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(WORK, U(J,:))) / D(J)
                  CALL DAXPY(P, B(J), VT(1, J), 1, WORK, 1)
              ENDDO
          END IF
      ELSE
          IF (LSAME(UPLO, 'U')) THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(WORK, VT(J,:))) / D(J)
                  CALL DAXPY(P, B(J), U(J, 1), P, WORK, 1)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(WORK, VT(J,:))) / D(J)
                  CALL DAXPY(P, B(J), U(J, 1), P, WORK, 1)
              ENDDO
          END IF
      END IF
*
*     End of DTRSV .
*
      END SUBROUTINE DTRSV