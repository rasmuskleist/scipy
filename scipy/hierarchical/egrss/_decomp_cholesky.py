import numpy as np
from .egrss import get_egrss_func


def cholesky(u, vh, d):
    """
    Computes the Cholesky factorization of a symmetric matrix `A` in extended
    generator form.

    The symmetric matrix `A` is represented as an extended generator
    representable semiseparable matrix:

    .. math::
        A = \\text{tril}(U V^H, -1) + \\text{triu}(VU^H, 1) + \\text{diag}(d)

    Parameters
    ----------
    u : numpy.ndarray
        Generator matrix `U` for the symmetric matrix `A`.
    vh : numpy.ndarray
        Generator matrix `V^H` of the symmetric matrix `A`.
    d : numpy.ndarray
        The diagonal entries of the matrix A.

    Returns
    -------
    wh : numpy.ndarray
        Generator `W^H` for the Cholesky factor `L`
    c : numpy.ndarray, optional
        The diagonal entries of the Cholesky factor `L`.

    Notes
    -----
    The Cholesky factorization `A=LL^H` only exists when the input matrix `A`
    is symmetric positive definite.

    The Cholesky factorization is represented in the extended generator form:

    .. math::
        A = \\text{tril}(U W^H, -1) + \\text{diag}(c)

    Computational Complexity
    ------------------------
    The computational complexity of this operation is `O(p^2n)`, where `n` is
    the size of the matrix and `p` is the semiseparability rank.

    References
    ----------
    [1] M. S. Andersen, and T. Chen, "Smoothing Splines and Rank Structured
        Matrices: Revisiting the Spline Kernel", SIAM Journal on Matrix
        Analysis and Applications, 41 (2020), pp. 389-412.
    """

    potrf = get_egrss_func('potrf', (u, vh, d))

    u, wh, c, info = potrf('L', u, vh, d)
    if info == 0:
        return u, wh, c
    else:
        raise ValueError(
            f"EGRSS reported an illegal value in {-info}-th argument"
            'on entry to "POTRF".'
        )
