import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from numpy import full
from numpy.random import randn
from scipy.hierarchical.egrss import cholesky
from scipy.linalg._testutils import assert_no_overwrite


class TestCholesky:
    def test_simple(self):
        U = randn(4, 2)
        V = randn(4, 2)
        d = full(4, 8, dtype=np.float64)

        A  = np.tril(U @ V.T) + np.triu(V @ U.T, 1) + np.diag(d)

        U, Wh, c = cholesky(U, V.T, d)
        L = np.tril(U @ Wh, -1) + np.diag(c)
        assert_array_almost_equal(A, L @ L.T)

        U = np.triu(Wh.T @ U.T, 1) + np.diag(c)
        assert_array_almost_equal(A, U.T @ U)
