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
from scipy.hierarchical.egrss._egrss import fib

class TestSolveTriangular:

    def test_simple(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """

        print(fib.__doc__)

        f = fib(8)
        print(f)

        assert_equal(f, [0, 1, 1, 2, 3, 5, 8, 13])

        #s = dsum(a, b)
        #print(s)
        #print(dsum.__doc__)

        #assert_almost_equal(s, a + b)

        #U = np.random.randn(4, 2)
        #Vh = np.random.randn(2, 4)
        #d = np.random.randn(4)
        #b = np.random.randn(4)
#
        #A = np.tril(U @ Vh, -1) + np.diag(d)
        #x = dtrsv("L", "N", U, Vh, d, b)
#
        #print(x)
        #print(A @ x)
        #print(b)
#
        #print(np.allclose(A @ x, b))
        #assert_array_almost_equal(A @ x, b)

        #U = np.asfortranarray(randn(4, 2), dtype=np.float64)
        #Vh = np.asfortranarray(randn(2, 4), dtype=np.float64)
        #print(U)
        #d = np.asfortranarray(randn(4), dtype=np.float64)
        #b = np.asfortranarray(randn(4), dtype=np.float64)
        #print(b)
#
        #A = tril(U @ Vh, -1) + diag(d)
        #x = solve_triangular(U, Vh, d, b)
#
        #print(A)
        #print(x)
        #print(b)
        #assert_array_almost_equal(A @ x, b)
