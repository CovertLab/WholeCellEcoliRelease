'''Parallelization utilities.'''

from __future__ import absolute_import

import multiprocessing
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
			multiprocessing.cpu_count()))
	return int(value)

def plotter_cpus():
	"""Return the number of CPUs available to use within an analysis plot,
	under the assumption that the analysis Firetask checks the
	"WC_ANALYZE_FAST" environment variable to run multiple plots in parallel,
	in which case it's counter-productive to fork a bunch more processes within
	the plot code.
	"""
	value = cpus()
	if "WC_ANALYZE_FAST" in os.environ:
		value = 1
	return value
