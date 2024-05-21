import numpy as np
from collections import deque

from scipy.hierarchical.egrss import get_egrss_func

class sif:
    ndim = 2

    def __init__(self, *args, **kwargs):
        self.dtype = kwargs.get('dtype', np.float64)
        self.order = kwargs.get('order', 'F')

        # TODO: Always property initialize u, vh, wh, c as empty arrays of shape (_, 0)
        if len(args) == 0:
            self.d = None
            self.a11 = None
            self.a22 = None
        elif len(args) == 1:
            (d,) = args
            self.d = d
            self.a11 = None
            self.a22 = None
        elif len(args) == 6:
            (a11, a22, u, vh, wh, c) = args
            self.a11 = a11
            self.a22 = a22
            self.u = u
            self.vh = vh
            self.wh = wh
            self.c = c
            self.d = None
        else:
            raise ValueError("Invalid number of input arguments")

    @property
    def shape(self):
        # TODO: Infer from the off-diagonal blocks
        n = m = 0
        for aii in LeafIterator(self):
            ni, mi = aii.d.shape
            n += ni
            m += mi

        return (n, m)

    def isleaf(self):
        return not (self.a11 and self.a22)

    def astype(self, dtype=None):
        pass

    def copy(self):
        pass

    def toarray(self):
        a = np.empty(self.shape)

        for aii, i, j in ReverseLevelOrderIterator(self):
            if aii.isleaf():
                a[i, j] = aii.d
            else:
                n1, m1 = aii.a11.shape
                n2, m2 = aii.a22.shape
                i11 = slice(i.start, i.start + n1)
                j11 = slice(j.start, j.start + m1)
                i22 = slice(i.stop - n2, i.stop)
                j22 = slice(j.stop - m2, j.stop)

                a22tilde = a[i22, j22]
                a21tilde = np.matmul(a22tilde, aii.u)
                a21 = np.matmul(a21tilde, aii.vh)

                trmm = get_egrss_func('trmm', dtype=self.dtype)
                a22h = trmm('L', 'T', -aii.u, aii.wh, aii.c, a22tilde.T.conj())

                a[i11, j22] = 0
                a[i22, j11] = a21
                a[i22, j22] = np.asfortranarray(a22h.T.conj())

        return a


class LevelOrderIterator:
    def __init__(self, a: sif):
        n, m = a.shape
        self.queue = deque([(a, slice(0, n), slice(0, m))])

    def __iter__(self):
        return self

    def __next__(self) -> sif:
        if self.queue:
            a, i, j = self.queue.popleft()
            if not a.isleaf():
                n1, m1 = a.a11.shape
                n2, m2 = a.a22.shape
                i11 = slice(i.start, i.start + n1)
                j11 = slice(j.start, j.start + m1)
                i22 = slice(i.stop - n2, i.stop)
                j22 = slice(j.stop - m2, j.stop)

                self.queue.append((a.a11, i11, j11))
                self.queue.append((a.a22, i22, j22))

            return a, i, j
        else:
            raise StopIteration


class ReverseLevelOrderIterator:
    # TODO: Investigate if this can be implemented in a lazy fashion
    def __init__(self, a: sif):
        self.stack = deque()

        n, m = a.shape
        queue = deque([(a, slice(0, n), slice(0, m))])
        while queue:
            a, i, j = queue.popleft()
            self.stack.appendleft((a, i, j))
            if not a.isleaf():
                n1, m1 = a.a11.shape
                n2, m2 = a.a22.shape
                i11 = slice(i.start, i.start + n1)
                j11 = slice(j.start, j.start + m1)
                i22 = slice(i.stop - n2, i.stop)
                j22 = slice(j.stop - m2, j.stop)

                queue.append((a.a22, i22, j22))
                queue.append((a.a11, i11, j11))

    def __iter__(self):
        return self

    def __next__(self) -> sif:
        if self.stack:
            return self.stack.popleft()
        else:
            raise StopIteration


class LeafIterator:
    def __init__(self, a: sif):
        self.queue = deque([a])

    def __iter__(self):
        return self

    def __next__(self) -> sif:
        while self.queue:
            a = self.queue.popleft()
            if not a.isleaf():
                self.queue.append((a.a11))
                self.queue.append((a.a22))
            else:
                return a

        raise StopIteration
