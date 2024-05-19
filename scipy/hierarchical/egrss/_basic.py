import numpy as np
from .egrss import get_egrss_func

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
    trsv("L", "N", u, vh, d, x)

    return x

