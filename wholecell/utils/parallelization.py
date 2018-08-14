'''Parallelization utilities.'''

from __future__ import absolute_import

import multiprocessing as mp
import os


def cpus():
	"""The number of CPUs available for a multiprocessing Pool whether running
	in SLURM or not. This reads the environment variable 'SLURM_CPUS_PER_TASK'
	containing the number of CPUs requested per task but since that's only set
	if the --cpus-per-task option was specified, this falls back to
	'SLURM_JOB_CPUS_PER_NODE' containing the number of processors available to
	the job on this node. Failing that, it reads the multiprocessing.cpu_count().

	By default, srun sets:
		SLURM_CPUS_ON_NODE=1
		SLURM_JOB_CPUS_PER_NODE=1

	srun -p mcovert --cpus-per-task=2:
		SLURM_CPUS_PER_TASK=2
		SLURM_CPUS_ON_NODE=2
		SLURM_JOB_CPUS_PER_NODE=2

	srun --ntasks-per-node=3 --cpus-per-task=4:
		SLURM_CPUS_PER_TASK=4
		SLURM_CPUS_ON_NODE=12
		SLURM_JOB_CPUS_PER_NODE=12

	srun --ntasks-per-node=3:
		SLURM_CPUS_ON_NODE=3
		SLURM_JOB_CPUS_PER_NODE=3

	See https://slurm.schedmd.com/sbatch.html
	"""
	value = os.environ.get('SLURM_CPUS_PER_TASK',
		os.environ.get('SLURM_JOB_CPUS_PER_NODE',
			mp.cpu_count()))
	return int(value)


def pool(processes=None):
	"""Return a multiprocessing.Pool or an InlinePool, depending on the
	requested number of processes. See `InlinePool` for why this is important.

	Requesting the default number of processes checks the SLURM CPU count,
	unlike `multiprocessing.Pool()`. See `cpus()`.
	"""
	if processes is None:
		processes = cpus()

	return mp.Pool(processes=processes) if processes > 1 else InlinePool()


class InlinePool(object):
	"""
	A substitute for multiprocessing.Pool() that runs the work inline in the
	current process. This is important because (1) a regular Pool worker cannot
	construct a nested Pool (even with processes=1) since daemon processes are
	not allowed to have children, and (2) it's easier to debug code running in
	the main process.
	"""

	def map(self, func, iterable):
		"""Map the function over the iterable."""
		return map(func, iterable)

	# TODO(jerry): Implement apply_async() if needed.

	def close(self):
		pass

	def terminate(self):
		pass

	def join(self):
		pass
