"""
Test the parallelization utility.
"""

import multiprocessing
import platform
import time
import unittest

import pytest

from wholecell.utils import parallelization


def _square(i):  # Runs in inner_pool
    square = i ** 2
    time.sleep(i / 100)
    return square


def _sum_squares(i, j):  # Runs in outer_pool
    with parallelization.pool(num_processes=2, nestable=False) as inner_pool:
        squares = inner_pool.map(_square, (i, j))

    sum_squares = sum(squares)
    time.sleep(0.01)
    return sum_squares


class Test_parallelization(unittest.TestCase):

    def _check_multi_sum(self, processes, nestable):
        i = 2
        j = 3
        expected_sum = i*i + j*j

        with parallelization.pool(num_processes=processes, nestable=nestable) as outer_pool:
            actual_sum = outer_pool.apply_async(_sum_squares, (i, j)).get()
            assert actual_sum == expected_sum

    def test_cpus(self):
        assert 1 <= parallelization.cpus() <= multiprocessing.cpu_count()
        assert 1 <= parallelization.cpus(4) <= 4

    def test_is_macos(self):
        assert parallelization.is_macos() == (platform.system() == 'Darwin')

    def test_nested_pools(self):
        """Test nested process pools. A nested multiprocessing.Pool will raise
        an AssertionError that daemonic processes are not allowed to have
        children.
        """
        self._check_multi_sum(1, nestable=False)  # InlinePool (always nestable)

        self._check_multi_sum(2, nestable=True)  # NoDaemonPool

        with pytest.raises(AssertionError, match='daemonic'):
            self._check_multi_sum(2, nestable=False)  # mp.Pool
