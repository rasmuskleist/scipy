from collections import deque

import numpy as np


class Partition:
    def __init__(self, *args):
        if len(args) == 0:
            self.index_set = None
            self.left_child = self.right_child = None
        if len(args) == 1:
            (self.index_set,) = args
            self.left_child = self.right_child = None
        elif len(args) == 2:
            (self.left_child, self.right_child) = args
            self.index_set = np.concatenate(
                self.left_child.index_set, self.right_child.index_set
            )
        else:
            raise RuntimeError

    def isleaf(self):
        return not (self.left_child and self.right_child)


def split_indices(size: int, block_size: int | np.ndarray[int]):
    """
    Create a list of index pairs that split a range into blocks of a given size.

    Parameters
    ----------
    size : int
        The size of the range to split.
    block_size : int | np.ndarray[int]
        The size of the blocks.

    Returns
    -------
    np.ndarray
        An array of index pairs.
    """

    if np.isscalar(block_size):
        block_count, remainder = np.divmod(size, block_size)
        block_size = np.full(block_count, block_size, dtype=int)
        block_size[:remainder] += 1

    idx = np.cumsum(block_size, dtype=int)
    return np.column_stack((idx - block_size, idx))


def cluster_tree(index_set: np.ndarray, block_size: int | np.ndarray[int]) -> Partition:
    """
    Create a cluster tree from an index set and a block size.

    Parameters
    ----------
    index_set : np.ndarray
        The index set to partition.
    block_size : int | np.ndarray[int]
        The size of the leaf blocks in the partition.

    Returns
    -------
    Partition
        The root of the cluster tree.

    """
    splits = split_indices(index_set.size, block_size)

    root = Partition(splits)
    queue = deque([root])

    while queue:
        node = queue.popleft()
        left_split, right_split = np.array_split(node.index_set, 2)

        if np.size(left_split) != 0 and np.size(right_split) != 0:
            node.left_child = Partition(left_split)
            node.right_child = Partition(right_split)

            queue.append(node.left_child)
            queue.append(node.right_child)

        i = node.index_set[0, 0]
        j = node.index_set[-1, -1]
        node.index_set = [i, j]

    return root


from collections import deque

import numpy as np

from ._partition import Partition


class LevelOrderIterator:
    def __init__(self, cluster_tree: Partition = None):
        self.queue = deque([cluster_tree])

    def __iter__(self):
        return self

    def __next__(self) -> Partition:
        if self.queue:
            subtree = self.queue.popleft()
            if not subtree.isleaf():
                self.queue.append(subtree.left_child)
                self.queue.append(subtree.right_child)

            return subtree
        else:
            raise StopIteration


class ReverseLevelOrderIterator:
    # TODO: Investigate if there is a better way to implement this perhaps using a stack
    # NOTE: This is litterally reverse level order iterator so the left and right children are reversed
    def __init__(self, cluster_tree: Partition = None):
        self.queue = deque([t for t in LevelOrderIterator(cluster_tree)])

    def __iter__(self):
        return self

    def __next__(self) -> Partition:
        if self.queue:
            return self.queue.pop()
        else:
            raise StopIteration


class LeafIterator:
    def __init__(self, cluster_tree: Partition = None):
        self.queue = deque([cluster_tree])

    def __iter__(self):
        return self

    def __next__(self) -> Partition:
        if self.queue:
            subtree = self.queue.popleft()
            if not subtree.isleaf():
                self.queue.append(subtree.left_child)
                self.queue.append(subtree.right_child)
            else:
                while not subtree.isleaf():
                    subtree = self.queue.popleft()
                return subtree
        else:
            raise StopIteration
