'''Parallelization utilities.'''

from __future__ import absolute_import, division, print_function

import multiprocessing as mp
import os


def is_macos():
	'''Return True if this is running on macOS.'''
	return os.uname()[0].lower() == 'darwin'


def cpus(requested_num_processes=None, **kwargs):
	"""Return the usable number of worker processes via `fork` (e.g. with
	`multiprocessing.Pool`), up to `requested_num_processes` (default: the max
	as reported by `multiprocessing.cpu_count()`) considering macOS and SLURM
	limitations. `1` means do the work in-process rather than forking
	subprocesses.

	On macOS: This returns 1 due to problems where `fork` can segfault or fail
	to parallelize -- unless the caller overrides that safety check. See
	Issue #392.

	TODO(jerry): Test if Python 3's `multiprocessing` "spawn" mode fixes the
		`fork` problems.

	On SLURM: This reads the environment variable 'SLURM_CPUS_PER_TASK'
	containing the number of CPUs requested per task but since that's only set
	if the --cpus-per-task option was specified, this falls back to
	'SLURM_JOB_CPUS_PER_NODE' containing the number of processors available to
	the job on this node.

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

	Args:
		requested_num_processes (int): the requested number of worker
			processes; pass None or 0 to default to the max available
		kwargs (Dict[str]): go ahead and pass in `advice='mac override'` to
			override the safety check if you're confident that `fork`ed
			processes parallelize OK in this caller on macOS; otherwise this
			function will return 1 on macOS
	Returns:
		num_cpus (int): the usable number of worker processes via `fork` (e.g.
			with `multiprocessing.Pool`) as limited by the hardware, macOS,
			SLURM, and `requested_num_processes`.

			==> 1 means DO NOT `fork` PROCESSES. (Try `exec`?)

	See also `pool()`.

	See https://slurm.schedmd.com/sbatch.html

	See https://github.com/CovertLab/wcEcoli/issues/392
	"""
	if is_macos() and kwargs.get('advice') != 'mac override':
		os_cpus = 1
	else:
		os_cpus = mp.cpu_count()

	value = os.environ.get('SLURM_CPUS_PER_TASK',
		os.environ.get('SLURM_JOB_CPUS_PER_NODE',
		os_cpus))
	slurm_cpus = int(value)

	available = min(os_cpus, slurm_cpus)

	if requested_num_processes is not None and requested_num_processes > 0:
		if requested_num_processes > available:
			print('Warning: Request for {} worker processes got limited to {}'
				.format(requested_num_processes, available))
		elif requested_num_processes < available:
			available = requested_num_processes

	return available


def pool(num_processes=None):
	"""Return an `InlinePool` if `cpus(num_processes) == 1`, else a
	`multiprocessing.Pool(cpus(num_processes))`, as suitable for the current
	runtime environment. See `cpus()` on figuring the number of usable
	processes and `InlinePool` about why running in-process is important.
	"""
	usable = cpus(num_processes)

	return mp.Pool(processes=usable) if usable > 1 else InlinePool()


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
