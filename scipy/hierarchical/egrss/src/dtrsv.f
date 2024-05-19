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
      Z = 0.0
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      INFO = 0
      IF (TRANS == 'N') THEN
          IF (UPLO == 'U') THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(Z, U(J,:))) / D(J)
                  Z = Z + VT(:,J) * B(J)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(Z, U(J,:))) / D(J)
                  Z = Z + VT(:,J) * B(J)
              ENDDO
          END IF
      ELSE
          IF (UPLO == 'U') THEN
              DO J = N, 1, -1
                  B(J) = (B(J) - DOT_PRODUCT(Z, VT(J,:))) / D(J)
                  Z = Z + U(:,J) * B(J)
              ENDDO
          ELSE
              DO J = 1, N
                  B(J) = (B(J) - DOT_PRODUCT(Z, VT(J,:))) / D(J)
                  Z = Z + U(:,J) * B(J)
              ENDDO
          END IF
      END IF
*
*     End of DTRSV .
*
      END SUBROUTINE DTRSV