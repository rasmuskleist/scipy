from collections import deque

import numpy as np

from ._partition import LevelOrderIterator, Partition
from ._sif import sif


def empty(t: Partition, dtype=None, order="F") -> sif:
    """
    Create an empty sif matrix from a partition.

    Parameters
    ----------
    t : Partition
        The partition to create the sif matrix from.
    dtype : type, optional
        The data type of the matrix.

    Returns
    -------
    sif
        An empty sif matrix.
    """
    iterator = LevelOrderIterator(t)
    root = sif(dtype=dtype, order=order)
    queue = deque([root])

    while queue:
        node = queue.popleft()
        subtree = next(iterator)

        if subtree.isleaf():
            i, j = subtree.index_set
            node.d = np.empty((j - i, j - i), dtype=dtype, order=order)
        else:
            node.a11 = sif(dtype=dtype, order=order)
            node.a22 = sif(dtype=dtype, order=order)
            queue.append(node.a11)
            queue.append(node.a22)

    return root
