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
    print(trsv.__doc__)

    # Solve the triangular system
    x = trsv("L", "N", u, vh, d, b, n = 4, m = 2)

    print(x)
    return x

