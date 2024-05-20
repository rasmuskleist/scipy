      SUBROUTINE DTRSV(UPLO,TRANS,N,P,U,LDU,VT,LDVT,D,B,WORK,LWORK,INFO)
*
*     .. Scalar Arguments ..
      INTEGER INFO,LDU,LDVT,LWORK,N,P
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION D(*), U(LDU,*), VT(LDVT,*), B(*), WORK(*)
*     ..
*
*  =====================================================================
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
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
      ELSE IF (LDU.LT.MAX(1,N) .OR. LDVT.LT.MAX(1,P)) THEN
          INFO = 6
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
      CALL DCOPY(P, 0.0D0, 0, WORK, 1)
      IF (LSAME(TRANS, 'N')) THEN
          IF (LSAME(UPLO, 'U')) THEN
              DO J = N, 1, -1
                  DO I = 1, P
                      B(J) = B(J) - U(J, I) * WORK(I)
                  ENDDO
                  B(J) = B(J) / D(J)
                  DO I = 1, P
                      WORK(I) = WORK(I) + B(J) * VT(I, J)
                  ENDDO
              ENDDO
          ELSE
              DO J = 1, N
                  DO I = 1, P
                      B(J) = B(J) - U(J, I) * WORK(I)
                  ENDDO
                  B(J) = B(J) / D(J)
                  DO I = 1, P
                      WORK(I) = WORK(I) + B(J) * VT(I, J)
                  ENDDO
              ENDDO
          END IF
      ELSE
          IF (LSAME(UPLO, 'U')) THEN
              DO J = 1, N
                  DO I = 1, P
                      B(J) = B(J) - VT(I, J) * WORK(I)
                  ENDDO
                  B(J) = B(J) / D(J)
                  DO I = 1, P
                      WORK(I) = WORK(I) + B(J) * U(J, I)
                  ENDDO
              ENDDO
          ELSE
              DO J = N, 1, -1
                  DO I = 1, P
                      B(J) = B(J) - VT(I, J) * WORK(I)
                  ENDDO
                  B(J) = B(J) / D(J)
                  DO I = 1, P
                      WORK(I) = WORK(I) + B(J) * U(J, I)
                  ENDDO
              ENDDO
          END IF
      END IF
*
*     End of DTRSV .
*
      END SUBROUTINE DTRSV