import numpy as np
import scipy as sp
from scipy.hierarchical.egrss import get_egrss_func
from scipy.linalg.lapack import get_lapack_funcs

from scipy.hierarchical import sif, LevelOrderIterator, ReverseLevelOrderIterator


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
    n, m = a.shape
    p, q = b.shape

    x = np.empty((m, q), order="F")  # dtype = a.dtype
    if trans == 0:
        for aii, i, j in ReverseLevelOrderIterator(a):
            if aii.isleaf():
                xii = x[j, :]
                bii = b[j, :]
                trtrs = get_lapack_funcs("trtrs", dtype=a.dtype)
                xii[:, :], _ = trtrs(aii.d, bii, lower=True, trans=trans)
            else:
                _, m1 = aii.a11.shape
                _, m2 = aii.a22.shape
                j1 = slice(j.start, j.start + m1)
                j2 = slice(j.stop - m2, j.stop)

                xii = x[j1, :]
                xjj = x[j2, :]

                xjj = xjj - aii.u @ (aii.vh @ xii)

                trsm = get_egrss_func("trsm", dtype=a.dtype)
                xjj[:,:] = trsm("L", "N", -aii.u, aii.wh, aii.c, xjj)

    elif trans == 1:
        btilde = np.copy(b, order="F")
        for aii, i, j in LevelOrderIterator(a):
            if aii.isleaf():
                xii = x[j, :]
                biitilde = btilde[j, :]
                trtrs = get_lapack_funcs("trtrs", dtype=a.dtype)
                xii[:, :], _ = trtrs(aii.d, biitilde, lower=True)
            else:
                _, m1 = aii.a11.shape
                _, m2 = aii.a22.shape
                j1 = slice(j.start, j.start + m1)
                j2 = slice(j.stop - m2, j.stop)

                biitilde = btilde[j1, :]
                bjjtilde = btilde[j2, :]

                trsm = get_egrss_func("trsm", dtype=a.dtype)
                bjjtilde[:,:] = trsm("L", "T", -aii.u, aii.wh, aii.c, bjjtilde)

                biitilde[:,:] = biitilde - aii.vh.T.conj() @ (
                    aii.u.T.conj() @ bjjtilde
                )
    else:
        raise ValueError("Invalid")

    return x
