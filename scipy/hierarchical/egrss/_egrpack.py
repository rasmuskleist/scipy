__all__ = ["potrf", "trmv", "trsv", "trsm", "trsm"]

import numpy as np


def trmv(u, vh, d, x, trans=0, dtype=np.float64):
    """
    Evaluates the matrix-vector product of an extended generator form of a
    lower triangular matrix `L` and a vector `x`.

    The lower triangular matrix `L` is represented in the extended generator
    form:

    .. math::
        L = \\text{tril}(UV^H, -1) + \\text{diag}(d)

    where `U` and `V` are generator matrices, and `d` is a vector representing
    the diagonal.

    Parameters
    ----------
    u : numpy.ndarray
        Generator matrix `U` for the extended generator form of `L`.
    vh : numpy.ndarray
        Generator matrix `V` for the extended generator form of `L`.
    d : numpy.ndarray
        The diagonal entries of `L`.
    x : numpy.ndarray
        Input vector `x` for the matrix-vector product.
    trans : {0, 1, 2, 'N', 'T', 'C'}, optional
        Type of product to form:

        ========  =========
        trans     system
        ========  =========
        0 or 'N'  Lx  = y
        1 or 'T'  L^T x = y
        2 or 'C'  L^H x = y
        ========  =========

    Returns
    -------
    y : numpy.ndarray
        The result of the matrix-vector product `Lx`.

    Computational Complexity
    ------------------------
    The computational complexity of this operation is `O(pn)`, where `n` is the
    size of the matrix and `p` is the semiseperability rank.

    References
    ----------
    [1] M. S. Andersen, and T. Chen, "Smoothing Splines and Rank Structured
        Matrices: Revisiting the Spline Kernel", SIAM Journal on Matrix
        Analysis and Applications, 41 (2020), pp. 389-412.
    """

    n, m = u.shape
    y = np.empty((n,), dtype=dtype)
    z = np.zeros((m,), dtype=dtype)

    if trans == 0:
        for k in range(n):
            y[k] = d[k] * x[k] + np.dot(z, u[k, :])
            z = z + vh[:, k] * x[k]
    elif trans == 1:
        for k in range(n - 1, -1, -1):
            y[k] = d[k] * x[k] + np.dot(z, vh[:, k])
            z = z + u[k, :] * x[k]

    return y


def trmm(u, vh, c, x, trans=0, dtype = np.float64):
    if x.ndim == 1:
        return trmv(u, vh, c, x, trans=trans)
    else:
        n, _ = u.shape
        _, m = x.shape

        y = np.empty((n, m), dtype=dtype)
        for k in range(m):
            y[:, k] = trmv(u, vh, c, x[:, k], trans=trans)

        return y
    

def trsv(u, vh, d, b, trans=0, dtype=np.float64):
    """
    Solves the system `Lx = b` where the lower triangular matrix `L` is in
    extended generator form.

    The lower triangular matrix `L` is represented in the extended generator
    form:

    .. math::
        L = \\text{tril}(UV^H, -1) + \\text{diag}(d)

    Parameters
    ----------
    u : numpy.ndarray
        Generator matrix `U` for the extended generator form of `L`.
    vh : numpy.ndarray
        Generator matrix `V^H` for the extended generator form of `L`.
    d : numpy.ndarray
        Diagonal entries of the triangular matrix.
    x : numpy.ndarray
        Input vector `x` for the system `Lx = b`.
    trans : {0, 1, 2, 'N', 'T', 'C'}, optional
        Type of system to solve:

        ========  =========
        trans     system
        ========  =========
        0 or 'N'  L x  = b
        1 or 'T'  L^T x = b
        2 or 'C'  L^H x = b
        ========  =========

    Returns
    -------
    x : numpy.ndarray
        The solution vector `x` for the system `Lx = b`.

    Computational Complexity
    ------------------------
    The computational complexity of solving the system `Lx = b` is `O(pn)`,
    where `n` is the size of the matrix and `p` is the semiseperability rank.

    References
    ----------
    [1] M. S. Andersen, and T. Chen, "Smoothing Splines and Rank Structured
        Matrices: Revisiting the Spline Kernel", SIAM Journal on Matrix
        Analysis and Applications, 41 (2020), pp. 389-412.
    """
    n, m = u.shape
    z = np.zeros((m,), dtype=dtype)
    x = np.empty_like(b, dtype=dtype)

    if trans == 0:
        for k in range(n):
            x[k] = (b[k] - np.dot(z, u[k, :])) / d[k]
            z = z + vh[:, k] * x[k]
    elif trans == 1:
        for k in range(n - 1, -1, -1):
            x[k] = (b[k] - np.dot(z, vh[:, k])) / d[k]
            z = z + u[k, :] * x[k]

    return x


def trsm(u, vh, c, b, trans=0, dtype=np.float64):
    if b.ndim == 1:
        return trsv(u, vh, c, b, trans=trans)
    else:
        n, _ = u.shape
        _, m = b.shape

        x = np.empty((n, m), dtype=dtype)
        for k in range(m):
            x[:, k] = trsv(u, vh, c, b[:, k], trans=trans)

        return x


def potrf(u, vh, d, dtype=np.float64):
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
    p = np.zeros((m, m), dtype=dtype)
    c = np.empty_like(d, dtype=dtype)
    wh = np.empty_like(vh, dtype=dtype)

    for k in range(n):
        wh[:, k] = vh[:, k] - np.dot(p, u[k, :])
        c[k] = np.sqrt(np.dot(u[k, :], wh[:, k]) + d[k])

        wh[:, k] = wh[:, k] / c[k]
        p = p + np.outer(wh[:, k], wh[:, k])

    return wh, c
