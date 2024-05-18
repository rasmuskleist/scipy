import itertools
import warnings

import numpy as np
from numpy import (arange, array, diag, dot, zeros, identity, tril, conjugate, transpose,
                   float32)
from numpy.random import random

from numpy.testing import (assert_equal, assert_almost_equal, assert_,
                           assert_array_almost_equal, assert_allclose,
                           assert_array_equal, suppress_warnings)
import pytest
from pytest import raises as assert_raises

from scipy.hierarchical.egrss import solve_triangular

from scipy.linalg._testutils import assert_no_overwrite
from scipy._lib._testutils import check_free_memory, IS_MUSL
from scipy.linalg.blas import HAS_ILP64


class TestSolveTriangular:

    def test_simple(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """
        U = random((16, 4))
        Vh = random((4, 16))
        d = random(16)
        b = random(16)

        A = tril(dot(U, Vh), -1) + diag(d)
        sol = solve_triangular(U, Vh, d, b)
        assert_array_almost_equal(dot(A, sol), b)
