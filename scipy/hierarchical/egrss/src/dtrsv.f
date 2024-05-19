      SUBROUTINE DTRSV(UPLO,TRANS,N,P,U,VT,D,B,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INFO,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(N), U(N,P), VT(P,N), B(N)
*     ..
*
*  =====================================================================
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION Z(P)
*     ..
*     .. External Functions ..
      EXTERNAL DDOT
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL DAXPY
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
      Z = 0.00
      INFO = 0
      CALL DCOPY(P, 0.0D0, 0, Z, 1)
      IF (LSAME(TRANS, 'N')) THEN
          IF (LSAME(UPLO, 'U')) THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(Z, U(J,:))) / D(J)
                  CALL DAXPY(P, B(J), VT(:, J), 1, Z, 1)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(Z, U(J,:))) / D(J)
                  CALL DAXPY(P, B(J), VT(:, J), 1, Z, 1)
              ENDDO
          END IF
      ELSE
          IF (LSAME(UPLO, 'U')) THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(Z, VT(J,:))) / D(J)
                  CALL DAXPY(P, B(J), U(J, :), P, Z, 1)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(Z, VT(J,:))) / D(J)
                  CALL DAXPY(P, B(J), U(J, :), P, Z, 1)
              ENDDO
          END IF
      END IF
*
*     End of DTRSV .
*
      END SUBROUTINE DTRSV