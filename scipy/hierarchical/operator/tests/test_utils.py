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

import scipy as sp
from scipy.sparse import diags
from scipy.linalg._testutils import assert_no_overwrite
from scipy._lib._testutils import check_free_memory, IS_MUSL
from scipy.linalg.blas import HAS_ILP64
from scipy.hierarchical.operator import LinearOperator, MatrixOperator, tosif
from scipy.hierarchical import LeafIterator, cluster_tree, empty


class Operator(LinearOperator):
    def __init__(self, u, vh, d):
        super().__init__(d.shape, np.float64)
        self.u = u
        self.vh = vh
        self.d = d

    def matvec(self, x) -> np.ndarray:
        xtilde = np.dot(self.vh, x)
        return np.dot(self.u, xtilde) + self.d @ x

    def rmatvec(self, x) -> np.ndarray:
        xtilde = np.dot(self.u.T, x)
        return np.dot(self.vh.T, xtilde) + self.d @ x

    def matmat(self, x) -> np.ndarray:
        xtilde = np.matmul(self.vh, x)
        return np.matmul(self.u, xtilde) + self.d @ x

    def rmatmat(self, x) -> np.ndarray:
        xtilde = np.matmul(self.u.T, x)
        return np.matmul(self.vh.T, xtilde) + self.d @ x

    def __getitem__(self, key) -> LinearOperator:
        # TODO: Fix this indexing
        i, j = key
        return Operator(self.u[i, :], self.vh[:, j], self.d[i, j])

    def toarray(self) -> np.ndarray:
        return np.dot(self.u, self.vh) + self.d.toarray()


class TestSif:
    def test_cholesky(self):
        block_size = np.full(32, 4)

        n = np.sum(block_size)
        index_set = np.arange(n)
        t = cluster_tree(index_set, block_size)

        n = np.sum(block_size)
        k = 1

        c = np.full(n, n)
        u = np.random.randn(n, k)
        d = diags(c, shape=(n, n), format="csr")

        a = u @ u.T + np.diag(c)
        op = Operator(u, u.T, d)

        l = tosif(op, t, k)
        cho = sp.linalg.cholesky(a, lower=True)
        assert_almost_equal(l.toarray(), cho)

    def test_operator(self):
        n = 64
        k = 8

        c = np.full(n, n)
        u = np.random.randn(n, k)
        d = diags(c, shape=(n, n), format="csr")

        a = u @ u.T + np.diag(c)
        op = Operator(u, u.T, d)

        x = np.random.randn(n)
        assert_almost_equal(op.matvec(x), np.dot(a, x))

        aii = a[:32, :32]
        opii = op[:32, :32]

        assert_almost_equal(aii, opii.toarray())
