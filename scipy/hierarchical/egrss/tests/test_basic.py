import itertools
import warnings

import numpy as np
from numpy import (arange, array, diag, dot, zeros, identity, tril, conjugate, transpose,
                   float32)
from numpy.random import randn

from numpy.testing import (assert_equal, assert_almost_equal, assert_,
                           assert_array_almost_equal, assert_allclose,
                           assert_array_equal, suppress_warnings,)
import pytest
from pytest import raises as assert_raises

from scipy.hierarchical.egrss import solve_triangular

from scipy.linalg._testutils import assert_no_overwrite
from scipy._lib._testutils import check_free_memory, IS_MUSL
from scipy.linalg.blas import HAS_ILP64
from scipy.hierarchical.egrss._egrss import dtrmv, dtrmm, dtrsm

class TestSolveTriangular:
    def test_multiply(self):
        U = np.random.randn(16, 8)
        Vh = np.random.randn(8, 16)
        d = np.ones(16)
        x = np.random.randn(16)

        A = np.tril(U @ Vh, -1) + np.diag(d)
        y = dtrmm('L', 'N', U, Vh, d, x)
        yh = dtrmm('L', 'T', U, Vh, d, x)

        assert_array_almost_equal(A @ x, y)
        assert_array_almost_equal(A.T @ x, yh)

        A = np.triu(U @ Vh, 1) + np.diag(d)
        y = dtrmm('U', 'N', U, Vh, d, x)
        yh = dtrmm('U', 'T', U, Vh, d, x)

        assert_array_almost_equal(A @ x, y)
        assert_array_almost_equal(A.T @ x, yh)

    def test_simple(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """

        U = np.random.randn(16, 8)
        Vh = np.random.randn(8, 16)
        d = np.ones(16)
        b = np.random.randn(16)

        A = np.tril(U @ Vh, -1) + np.diag(d)
        x = solve_triangular(U, Vh, d, b, lower=True)
        xh = solve_triangular(U, Vh, d, b, lower=True, trans=1)

        assert_array_almost_equal(A @ x, b)
        assert_array_almost_equal(A.T @ xh, b)

        A = np.triu(U @ Vh, 1) + np.diag(d)
        x = solve_triangular(U, Vh, d, b, lower=False)
        xh = solve_triangular(U, Vh, d, b, lower=False, trans=1)

        assert_array_almost_equal(A @ x, b)
        assert_array_almost_equal(A.T @ xh, b)

    def test_dtrsm(self):
        U = np.random.randn(16, 8)
        Vh = np.random.randn(8, 16)
        d = np.ones(16)
        b = np.random.randn(16)

        A = np.tril(U @ Vh, -1) + np.diag(d)
        x = dtrsm('L','N',U, Vh, d, b)
        xh = dtrsm('L','T',U, Vh, d, b)

        assert_array_almost_equal(A @ x, b)
        assert_array_almost_equal(A.T @ xh, b)

        A = np.triu(U @ Vh, 1) + np.diag(d)
        x = dtrsm('U','N',U, Vh, d, b)
        xh = dtrsm('U','T',U, Vh, d, b)

        assert_array_almost_equal(A @ x, b)
        assert_array_almost_equal(A.T @ xh, b)


