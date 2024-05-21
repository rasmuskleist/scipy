import numpy as np
import scipy as sp
from scipy.linalg.lapack import get_lapack_funcs

from scipy.hierarchical import sif, LevelOrderIterator, ReverseLevelOrderIterator
from scipy.hierarchical.egrss import get_egrss_func


def solve_triangular(a: sif, b: np.ndarray, trans: int = 0) -> np.ndarray:
    """
    Solve the system `ax = b` where `a` is a triangular matrix on the extended generator form.

    Parameters
    ----------
    a : sif
        Triangular matrix `a` in structured incomplete factorization form.
    b : numpy.ndarray
        Right-hand side vector `b`.

    Returns
    -------
    x : numpy.ndarray
        Solution vector `x` for the system `ax = b`.
    """

    b = np.atleast_2d(b)
    n, m = a.shape
    p, q = b.shape

    x = np.empty((m, q), order="F")  # dtype = a.dtype
    if trans == 0:
        for aii, i, j in ReverseLevelOrderIterator(a):
            if aii.isleaf():
                xii = x[j, :]
                bii = b[j, :]

                trtrs = get_lapack_funcs("trtrs", (aii.d, bii), dtype=a.dtype)
                xii[:, :], info = trtrs(aii.d, bii, lower=True, trans=trans)
            else:
                _, m1 = aii.a11.shape
                _, m2 = aii.a22.shape
                j1 = slice(j.start, j.start + m1)
                j2 = slice(j.stop - m2, j.stop)

                xii = x[j1, :]
                xjj = x[j2, :]

                xjjtilde = xjj - aii.u @ (aii.vh @ xii)
                trsm = get_egrss_func("trsm", (aii.u, aii.vh, aii.d, xjjtilde), dtype=a.dtype)

                egtrans = {0: "N", 1:"T", 2:"C"}.get(trans, trans)
                xjj[:, :] = trsm(-aii.u, aii.wh, aii.c, xjjtilde, trans=egtrans)

    elif trans == 1:
        btilde = np.copy(b)
        for aii, i, j in LevelOrderIterator(a):
            if aii.isleaf():
                xii = x[j, :]
                biitilde = btilde[j, :]

                trtrs = get_lapack_funcs("trtrs", (aii.d, biitilde), dtype=a.dtype)
                xii[:, :], info = trtrs(aii.d, biitilde, lower=True, trans=trans)
            else:
                _, m1 = aii.a11.shape
                _, m2 = aii.a22.shape
                j1 = slice(j.start, j.start + m1)
                j2 = slice(j.stop - m2, j.stop)

                biitilde = btilde[j1, :]
                bjjtilde = btilde[j2, :]

                egtrans = {0: "N", 1:"T", 2:"C"}.get(trans, trans)
                trsm = get_egrss_func("trsm", (aii.u, aii.vh, aii.d, bjjtilde), dtype=a.dtype)
                bjjtilde[:, :] = trsm(-aii.u, aii.wh, aii.c, bjjtilde, trans=egtrans)
                biitilde[:, :] = biitilde - aii.vh.T.conj() @ (
                    aii.u.T.conj() @ bjjtilde
                )
    else:
        raise ValueError("Invalid")

    return x
