import numpy as np
from .egrss import get_egrss_func
from scipy.linalg import LinAlgError

def solve_triangular(u, vh, d, b, trans=0, lower=False, overwrite_b=False):
    """
    Solve a triangular system of equations.

    Parameters
    ----------
    u : array_like
        The matrix U in the QR decomposition.
    vh : array_like
        The matrix V^H in the QR decomposition.
    """

    uplo = 'L' if lower else 'U'
    trans = {0: 'N', 1: 'T', 2: 'C'}.get(trans, trans)

    # Get the function
    trsm = get_egrss_func('trsm', (u, vh, d, b))

    # Solve the triangular system
    x = trsm(uplo, trans, u, vh, d, b, overwrite_b=overwrite_b)
    return x

