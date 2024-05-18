import numpy as np

from scipy.hierarchical.egrss import _fegrss
from scipy.linalg.blas import find_best_blas_type


def get_egrss_func(name, arrays=(), dtype=None):
    """
    Return available EGRSS function objects from the name

    Parameters
    ----------
    name : str
        The name of the function to retrieve
    arrays : tuple
        The input arrays to determine the data type from
    dtype : type, optional
        The data type to use

    Returns
    -------
    function
        The function object

    """

    prefix, dtype, _ = find_best_blas_type(arrays, dtype)

    func = getattr(_fegrss, f"{prefix}{name}")
    if func is None:
        raise ValueError(f"Function {name} not found in the EGRSS module")

    return func
