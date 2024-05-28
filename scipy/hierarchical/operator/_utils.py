import numpy as np
from scipy.linalg.lapack import get_lapack_funcs

from scipy.hierarchical.operator import LinearOperator
from scipy.hierarchical import sif, Partition, ReverseLevelOrderIterator, empty
from scipy.hierarchical.operator.randalg import rsvds
from scipy.hierarchical.linalg import solve_triangular
from scipy.hierarchical.egrss import get_egrss_func


class _PreconditionedMatrix(LinearOperator):
    def __init__(self, a: LinearOperator, l1: sif, l2: sif) -> sif:
        self.a = a
        self.l1 = l1
        self.l2 = l2
        _, m1 = l1.shape
        n2, _ = l2.shape
        super().__init__((m1, n2), a.dtype)

    def matvec(self, x: np.ndarray) -> np.ndarray:
        xtilde = solve_triangular(self.l1, x, trans=1)
        xtilde = self.a.matvec(xtilde)
        return solve_triangular(self.l2, xtilde, trans=0)

    def rmatvec(self, x: np.ndarray) -> np.ndarray:
        xtilde = solve_triangular(self.l2, x, trans=1)
        xtilde = self.a.rmatvec(xtilde)
        return solve_triangular(self.l1, xtilde, trans=0)

    def matmat(self, x: np.ndarray) -> np.ndarray:
        xtilde = solve_triangular(self.l1, x, trans=1)
        xtilde = self.a.matmat(xtilde)
        return solve_triangular(self.l2, xtilde, trans=0)

    def rmatmat(self, x: np.ndarray) -> np.ndarray:
        xtilde = solve_triangular(self.l2, x, trans=1)
        xtilde = self.a.rmatmat(xtilde)
        return solve_triangular(self.l1, xtilde, trans=0)


def tosif(a: LinearOperator, t: Partition, k, r=5) -> sif:
    """
    Compute the Cholesky factorization of a symmetric positive definite matrix.

    Parameters
    ----------
    a : LinearOperator
        The matrix to factorize.
    t : Partition
        The partition of the matrix.
    k : int
        The number of random vectors to use in the incomplete factorization.
    r : int, optional
        The number of iterations in the incomplete factorization.
    delta : float, optional
        The tolerance in the incomplete factorization.

    Returns
    -------
    sif
        The Cholesky factorization of the matrix.

    Notes
    -----
    The incomplete factorization is computed using the randomized incomplete factorization method.
    """
    l = empty(t, order="F")

    # TODO: Perform array validaiton (check finite)
    for li, i, j in ReverseLevelOrderIterator(l):
        if li.isleaf():
            np.copyto(li.d, a[i, j])
            (potrf,) = get_lapack_funcs(("potrf",), arrays=(li.d,))
            _, info = potrf(li.d, lower=True, overwrite_a=True, clean=True)
            if info != 0:
                raise ValueError(
                    f"LAPACK reported an illegal value in {-info}-th argument"
                    'on entry to "POTRF".'
                )

        else:
            l1 = li.a11
            l2 = li.a22

            _, n1 = l1.shape
            m2, _ = l2.shape
            j11 = slice(j.start, j.start + n1)
            i22 = slice(i.stop - m2, i.stop)

            a21 = _PreconditionedMatrix(a[i22, j11], l1, l2)
            u, s, vh = rsvds(a21, k=k, r=r)

            u = u * s
            potrf = get_egrss_func("potrf", dtype=a.dtype)
            _, wh, c, info = potrf('L', -u, u.T.conj(), np.ones(m2))
            if info == 0:
                li.u = u
                li.vh = vh
                li.wh = wh
                li.c = c
            else:
                raise ValueError(
                    f"LAPACK reported an illegal value in {-info}-th argument"
                    'on entry to "POTRF".'
                )

    return l
