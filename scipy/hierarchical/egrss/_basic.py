import numpy as np
from .egrss import get_egrss_func
from scipy.linalg import LinAlgError

def solve_triangular(u, vh, d, b):
    """
    Solve a triangular system of equations.

    Parameters
    ----------
    u : array_like
        The matrix U in the QR decomposition.
    vh : array_like
        The matrix V^H in the QR decomposition.
    """

    # Get the function
    trsv = get_egrss_func('trsv', (u, vh, d, b))

    x = np.copy(b, order='F')

    # Solve the triangular system
    info = trsv("L", "N", u, vh, d, x)

    if info == 0:
        return x
    if info > 0:
        raise LinAlgError("singular matrix: resolution failed at diagonal %d" %
                          (info-1))
    raise ValueError('illegal value in %dth argument of internal trtrs' %
                     (-info))

