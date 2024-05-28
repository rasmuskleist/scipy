import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_equal
from pytest import raises as assert_raises

from numpy import array, transpose, dot, conjugate, zeros_like, empty, full
from numpy.random import randn
from scipy.hierarchical.egrss import cholesky
from scipy.hierarchical.egrss._egrss import dpotrf
from scipy.linalg._testutils import assert_no_overwrite


class TestCholesky:

    def test_simple(self):
        U = randn(4, 2)
        V = randn(4, 2)
        d = full(4, 4, dtype=np.float64)

        A  = np.tril(U @ V.T) + np.triu(V @ U.T, 1) + np.diag(d)

        Wh, c = cholesky(U, V.T, d)
        print(U)
        print(Wh)
        print(c)

        L = np.tril(U @ Wh, -1) + np.diag(c)
        assert_array_almost_equal(A, L @ L.T)

        U, Wh, c, info = dpotrf('L', U, V.T, d)
        print(U)
        print(Wh)
        print(c)

        assert_equal(0, 1)
