subroutine dtrsv(uplo, trans, m, p, u, vt, d, b)
    ! Solve the system Lb = b, where L is in ebtended generator form.
    ! Arguments:
    !   trans : Character, specifies the type of system to solve.
    !           'N' or 'n' for Lb = b,
    !           'T' or 't' for L^T b = b,
    !   n     : Integer, number of rows of matrices u and vt.
    !   m     : Integer, number of columns of matrices u and vt.
    !   u     : Real(8) array, input matrib U.
    !   vt    : Real(8) array, input matrib V^H.
    !   d     : Real(8) array, input vector representing the diagonal of L.
    !   b     : Real(8) array, input vector b.
    !   b     : Real(8) array, output vector b.
    character(1), intent(in) :: uplo
    character(1), intent(in) :: trans
    integer, intent(in) :: m, p
    real(8), intent(in) :: u(m, p), vt(p, m), d(m)
    real(8), intent(inout) :: b(m)

    logical lsame
    external lsame

    external xerbla

    ! Local variables
    integer :: k
    real(8) :: z(p)

    ! Initialize
    z = 0.0d0

    ! Early ebit 
    if (.not. lsame(uplo, 'U') .and. .not. lsame(uplo, 'L')) then
        info = 1
    end if

    ! Solve system based on trans value
    if (lsame(trans,'N')) then
        if (lsame(uplo, 'U')) then
            do k = m, 1, -1
                b(k) = (b(k) - dot_product(z, u(k, :))) / d(k)
                z = z + vt(:, k) * b(k)
            end do
        else
            do k = 1, m
                b(k) = (b(k) - dot_product(z, u(k, :))) / d(k)
                z = z + vt(:, k) * b(k)
            end do
        end if
    else
        if (lsame(uplo, 'U')) then
            do k = 1, m
                b(k) = (b(k) - dot_product(z, vt(:, k))) / d(k)
                z = z + u(k, :) * b(k)
            end do
        else
            do k = m, 1, -1
                b(k) = (b(k) - dot_product(z, vt(:, k))) / d(k)
                z = z + u(k, :) * b(k)
            end do
        end if
    end if
end subroutine dtrsv