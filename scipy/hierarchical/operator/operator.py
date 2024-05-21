import numpy as np


class LinearOperator:
    """
    Base class for linear operators.

    Parameters
    ----------
    shape : tuple
        Shape of the operator.
    dtype : dtype
        Data type of the operator.

    Attributes
    ----------
    shape : tuple
        Shape of the operator.
    dtype : dtype
        Data type of the operator.
    """

    def __init__(self, shape, dtype):
        self.shape = shape
        self.dtype = dtype

    def transpose(self):
        raise NotImplementedError

    def adjoint(self):
        raise NotImplementedError

    def inverse(self):
        """Not every linear operator has an inverse."""
        raise NotImplementedError

    def inv(self):
        return self.inverse()

    def conj(self):
        raise NotImplementedError

    def toarray(self):
        raise NotImplementedError

    def matvec(self, x) -> np.ndarray:
        raise NotImplementedError

    def rmatvec(self, x) -> np.ndarray:
        raise NotImplementedError

    def matmat(self, x) -> np.ndarray:
        raise NotImplementedError

    def rmatmat(self, x) -> np.ndarray:
        raise NotImplementedError

    def __getitem__(self, key):
        raise NotImplementedError

    def __array__(self):
        return self.toarray()


class TransposeLinearOperator(LinearOperator):
    def __init__(self, A: LinearOperator):
        self.A = A

    def transpose(self):
        return self.A

    def adjoint(self):
        return self.A.conj()

    def matvec(self, x):
        return self.A.rmatvec(x)

    def rmatvec(self, x):
        return self.A.matvec(x)

    def matmat(self, x):
        return self.A.rmatmat(x)

    def rmatmat(self, x):
        return self.A.matmat(x)

    def toarray(self):
        return self.A.toarray().T

    def __getitem__(self, key) -> LinearOperator:
        return TransposeLinearOperator(self.A[key])


class MatrixOperator(LinearOperator):
    def __init__(self, a):
        a1 = np.asarray(a)
        a1 = np.atleast_2d(a1)
        self.A = a1

        super().__init__(a.shape, a.dtype)

    def matvec(self, x) -> np.ndarray:
        return np.dot(self.A, x)

    def rmatvec(self, x) -> np.ndarray:
        return np.dot(x, self.A)

    def matmat(self, x) -> np.ndarray:
        return np.matmul(self.A, x, order="F")

    def rmatmat(self, x) -> np.ndarray:
        return np.matmul(self.A.T, x, order="F")

    def __getitem__(self, key) -> LinearOperator:
        return MatrixOperator(self.A[key])

    def toarray(self) -> np.ndarray:
        return self.A

    def __array__(self):
        return self.A


class LowRankOperator(LinearOperator):
    def __init__(self, u: np.ndarray, vh: np.ndarray):
        self.u = u
        self.vh = vh
        super().__init__((u.shape[0], vh.shape[1]), u.dtype)

    def matvec(self, x) -> np.ndarray:
        xtilde = np.dot(self.vh, x)
        return np.dot(self.u, xtilde)

    def rmatvec(self, x) -> np.ndarray:
        xtilde = np.dot(self.u.T.conj(), x)
        return np.dot(self.vh.T.conj(), xtilde)

    def matmat(self, x) -> np.ndarray:
        xtilde = np.matmul(self.vh, x)
        return np.matmul(self.u, xtilde)

    def rmatmat(self, x) -> np.ndarray:
        xtilde = np.matmul(self.u.T.conj(), x)
        return np.matmul(self.vh.T.conj(), xtilde)

    def __getitem__(self, key) -> LinearOperator:
        # TODO: Fix this indexing
        i, j = key
        return LowRankOperator(self.u[i, :], self.vh[:, j])

    def toarray(self) -> np.ndarray:
        return np.dot(self.u, self.vh)

    def transpose(self) -> LinearOperator:
        return LowRankOperator(self.vh.T, self.u.T)

    def adjoint(self) -> LinearOperator:
        return LowRankOperator(self.vh.T.conj(), self.u.T.conj())
