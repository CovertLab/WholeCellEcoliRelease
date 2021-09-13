'''Parallelization utilities.'''

import multiprocessing as mp
import multiprocessing.pool
import os

from typing import Any, Callable, Dict, Iterable, List, Optional, Union


def is_macos():
	# type: () -> bool
	'''Return True if this is running on macOS.'''
	return os.uname()[0].lower() == 'darwin'


def cpus(requested_num_processes=None):
	# type: (Optional[int]) -> int
	"""Return the usable number of worker processes for a multiprocessing Pool,
	up to `requested_num_processes` (default = max available), considering SLURM
	and any other environment-specific limitations.
	`1` means do all work in-process rather than forking subprocesses.

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
		requested_num_processes: the requested number of worker
			processes; None or 0 means return the max usable number

	Returns:
		num_cpus: the usable number of worker processes for a Pool, as limited
			by the hardware, OS, SLURM, and `requested_num_processes`.

			==> 1 means DO NOT CREATE WORKER SUBPROCESSES.

	See also `pool()`.

	See https://slurm.schedmd.com/sbatch.html

	See https://github.com/CovertLab/wcEcoli/issues/392
	"""
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


def pool(num_processes=None, nestable=False):
	# type: (Optional[int], bool) -> Union[mp.pool.Pool, InlinePool]
	"""Return an `InlinePool` if `cpus(num_processes) == 1`, else a
	multiprocessing `Pool(cpus(num_processes))`, as suitable for the current
	runtime environment.

	This uses the 'spawn' process start method to create a fresh python
	interpreter process, avoiding threading problems and cross-platform
	inconsistencies.

	nestable can create a pool of non-daemon worker processes that can spawn
	nested processes and have have nested pools.

	See `cpus()` on figuring the number of usable processes.
	See `InlinePool` about why running in-process is important.
	"""
	usable = cpus(num_processes)

	if usable == 1:
		return InlinePool()
	elif nestable:
		return NoDaemonPool()
	else:
		return mp.get_context(method='spawn').Pool(processes=usable)


class InlinePool(object):
	"""
	A substitute for multiprocessing.Pool() that runs the work inline in the
	current process. This is important because (1) a regular Pool worker cannot
	construct a nested Pool (even with processes=1) since daemon processes are
	not allowed to have children, and (2) it's easier to debug code running in
	the main process.
	"""

	def map(self, func, iterable, chunksize=None):
		# type: (Callable[..., Any], Iterable[Any], Optional[int]) -> List[Any]
		"""Map the function over the iterable."""
		return list(map(func, iterable))

	def apply_async(self, func, args=(), kwds=None, callback=None):
		# type: (Callable[..., Any], Iterable[Any], Optional[Dict[str, Any]], Optional[Callable[..., None]]) -> ApplyResult
		"""
		Apply the function to the args serially (not asynchronously since
		only one process available).
		"""
		if kwds is None:
			kwds = {}
		result = func(*args, **kwds)
		if callback:
			callback(result)
		return ApplyResult(result)

	def close(self):
		pass

	def terminate(self):
		pass

	def join(self):
		pass

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_val, exc_tb):
		pass


class ApplyResult(object):
	"""
	A substitute for multiprocessing.ApplyResult() to return with apply_async.
	Will get created after a successful function call so ready() and
	successful() are always True.
	"""
	def __init__(self, result):
		self._result = result

	def ready(self):
		return True

	def successful(self):
		return True

	def wait(self, timeout=None):
		pass

	def get(self, timeout=None):
		return self._result


class NoDaemonProcess(mp.Process):
	@property  # type: ignore[override]
	def daemon(self):
		return False

	@daemon.setter
	def daemon(self, value):
		pass


class NoDaemonContext(type(mp.get_context())):  # type: ignore
	Process = NoDaemonProcess


class NoDaemonPool(mp.pool.Pool):
	"""
	A substitute for multiprocessing.Pool() that creates a pool that is not a
	daemonic process. This allows for nesting pool calls that would otherwise
	be prevented with an assertion error (AssertionError: daemonic processes
	are not allowed to have children).
	"""
	def __init__(self, *args, **kwargs):
		kwargs['context'] = NoDaemonContext()
		super().__init__(*args, **kwargs)
