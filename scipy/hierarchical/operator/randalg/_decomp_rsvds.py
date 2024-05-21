import numpy as np
from scipy.linalg import norm, qr
from scipy.sparse.linalg import svds

from scipy.hierarchical.operator import LinearOperator


def _rsvds(a: LinearOperator, eps=1e-6, p=1, r=5):
    """
    Compute a low-rank factorization A = UV^T with a predefined accuracy

    This function computes a low-rank factorization A = UV^T using random
    sampling to approximate a basis matrix Q for the range of A. Likewise The
    approximation error ||(I-QQ^T)A|| is approximated using random sampling.
    The function uses Q to help compute the low-rank factorization A = UV^T.

    Parameters
    ----------
    A : LinearOperator
        Function to evaluate matrix-vector products with A. The calling
        signature is ``fun(x, transpose)`` where where ``x`` is an ndarray with
        shape (n,) and ``transpose`` is a boolean indicating wether the
        transpose should be computed. ``fun`` must return an ndarray of the
        same shape as ``x``. See `vectorized` for more information.
    eps : float
        Maximum allowed tolerance for the approximation ||(I-QQ^T)A||

    Returns
    -------
    u : (M, K) ndarray
        U in the factorization A = UV^T
    vh : (K, N) ndarray
        V^T in the factorization A = UV^T

    References
    ----------
    .. [1]  F. Woolfe, E. Liberty, V. Rokhlin, and M. Tygert, "A fast
        randomized algorithm for the approximation of matrices", Appl. Comput.
        Harmon. Anal., 25 (2008), pp. 335-366.

    .. [2]  N. Halko, P. G. Martinsson, and J. A. Tropp, "Finding Structure
        with Randomness: Probabilistic Algorithms For Constructing Approximate
        Matrix Decompositions", Siam Review, 2021

    """

    m, n = a.shape

    # TODO: Make r an argument
    r = 10
    w = np.random.randn(n, r)

    critical = eps / (10 * np.sqrt(2 / np.pi))

    y = a.matmat(w)
    ynorm = norm(y, ord=2, axis=0)

    j = 0
    q = np.zeros((m, 0))

    # TODO: Implement maxiter
    # TODO: Use block samples
    # Find orthonormal basis Q for the range of A satisfying the tolerance
    while np.max(ynorm[-r:]) > critical:
        y[:, j] -= q @ q.T @ y[:, j]
        qj = y[:, j] / norm(y[:, j], ord=2)
        q = np.column_stack([q, qj])

        w = np.random.randn(n, 1)
        aw = a.matmat(w)

        # TODO: Implement power iteration

        y = np.column_stack([y, aw - q @ q.T @ aw])
        ynorm = norm(y, ord=2, axis=0)

        # TODO: Should this be vectorized?
        for i in range(j + 1, j + r):
            y[:, i] -= qj * np.inner(qj, y[:, i])

        j += 1

    # Compute B = Q.T @ A which yields the low-rank factorization Q @ B
    btilde = a.rmatmat(q)

    # TODO: Verify that using svds is correct
    # Compute an SVD factorization of the small matrix B
    u, s, vh = svds(btilde.T.conj(), j - 1)
    u = q @ u

    return u, s, vh


def rsvds(a: LinearOperator, k=10, r=5, p=0):
    """
    Compute a low-rank SVD factorization using a random power iteration
    method

    This function computes a low-rank SVD factorization using a random power
    iteration method. Specifically, it computes an m by k orthonormal matrix
    Q, whose range approximates the range of an m by n matrix A. The function
    uses Q to help compute the low-rank factorization A = UV^T.

    Parameters
    ----------
    A : Linear Operator
        Function to evaluate matrix-vector products with A. The calling
        signature is ``fun(x, transpose)`` where where ``x`` is an ndarray with
        shape (n,) and ``transpose`` is a boolean indicating wether the
        transpose should be computed. ``fun`` must return an ndarray of the
        same shape as ``x``. See `vectorized` for more information.
    k : int
        The target rank of the input matrix A
    p : int, optional
        Power iteration parameter, which is used if the singular values of A
        decrease very slowly.
    r : int, optional
        Oversapling parameter.

    Returns
    -------
    u : (M, K) ndarray
        U in the factorization A = USV^T
    s : (K) ndarray
        S
    vh : (K, N) ndarray
        V^T in the factorization A = UV^T

    References
    ----------
    .. [1]  N. Halko, P. G. Martinsson, and J. A. Tropp, "Finding Structure
        with Randomness: Probabilistic Algorithms For Constructing Approximate
        Matrix Decompositions", Siam Review, 2021

    """
    # TODO: Look at docstring from svd
    _, n = a.shape

    # Draw an n by (k + r) Gaussian random matrix
    w = np.random.randn(n, k + r)

    # Compute initial range and QR factorization
    y = a.matmat(w)
    q, _ = qr(y, mode="economic")

    # Form the m by (k + r) matrix Y = (A @ A.T)^p @ Q
    for _ in range(p):
        ytilde = a.rmatmat(q)
        qtilde, _ = qr(ytilde, mode="economic")

        y = a.matmat(qtilde)
        q, _ = qr(y, mode="economic")

    # Compute B = Q.T @ A which yields the low-rank factorization Q @ B
    btilde = a.rmatmat(q)

    # Compute an SVD factorization of the small matrix B
    u, s, vh = svds(btilde.T.conj(), k)
    u = q @ u

    return u, s, vh
