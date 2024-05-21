import numpy as np


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

    n, m = u.shape
    p = np.zeros((m, m))
    c = np.empty_like(d)
    wh = np.empty_like(vh)

    for k in range(n):
        wh[:, k] = vh[:, k] - np.dot(p, u[k, :])
        c[k] = np.sqrt(np.dot(u[k, :], wh[:, k]) + d[k])

        wh[:, k] = wh[:, k] / c[k]
        p = p + np.outer(wh[:, k], wh[:, k])

    return wh, c
