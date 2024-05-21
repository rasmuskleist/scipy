from unittest import TestCase

import numpy as np
import scipy as sp
from scipy.linalg import cholesky
from scipy.sparse import diags

from hmpy.hieraichal import LeafIterator, cluster_tree, empty
from hmpy.hieraichal.linalg import solve_triangular
from hmpy.operator import LinearOperator, MatrixOperator, tosif


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


class TestSif(TestCase):
    def test_cholesky(self):
        block_size = np.full(2, 4)

        n = np.sum(block_size)
        index_set = np.arange(n)
        t = cluster_tree(index_set, block_size)

        n = np.sum(block_size)
        k = 2

        c = np.full(n, n)
        u = np.random.randn(n, k)
        d = diags(c, shape=(n, n), format="csr")

        a = u @ u.T + np.diag(c)
        op = Operator(u, u.T, d)

        l = tosif(op, t, k)
        cho = cholesky(a, lower=True)

        print(l.toarray())
        print(cho)

        self.assertTrue(np.allclose(l.toarray(), cho))

    def test_operator(self):
        n = 64
        k = 8

        c = np.full(n, n)
        u = np.random.randn(n, k)
        d = diags(c, shape=(n, n), format="csr")

        a = u @ u.T + np.diag(c)
        op = Operator(u, u.T, d)

        x = np.random.randn(n)

        self.assertTrue(np.allclose(op.matvec(x), np.dot(a, x)))

        aii = a[:32, :32]
        opii = op[:32, :32]

        self.assertTrue(np.allclose(aii, opii.toarray()))

        # self.assertTrue(False)
