"""Generic Sisyphus/Gaia/Google Cloud workflow builder."""

# TODO(jerry): Use different utilities than os.path functions to construct
#  paths to use inside the linux containers when the builder runs on Windows.

from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import os
import sys
if os.name == 'posix' and sys.version_info[0] < 3:
	import subprocess32 as subprocess
else:
	import subprocess
from typing import Any, Callable, Dict, Iterable, List, Optional

from gaia.client import Gaia
from requests import ConnectionError

from wholecell.utils import filepath as fp


# Config details to pass to Gaia.
# ASSUMES: gaia_host is reachable e.g. via an ssh tunnel set up by
# runscripts/cloud/ssh-tunnel.sh.
GAIA_CONFIG = {'gaia_host': 'localhost:24442'}

STDOUT_PATH = '>'  # special pathname that captures stdout + stderror

MAX_WORKERS = 500  # don't launch more than this many worker nodes at a time


def _rebase(path, internal_prefix, storage_prefix):
	# type: (str, str, str) -> str
	"""Return a path rebased from internal_prefix to storage_prefix."""
	new_path = os.path.join(storage_prefix, os.path.relpath(path, internal_prefix))

	assert '..' not in new_path, (
		'''Can't rebase path "{}" that doesn't start with internal_prefix "{}"'''
			.format(path, internal_prefix))
	return new_path

def _keyify(paths, fn=lambda path: path):
	# type: (Iterable[str], Callable[[str], str]) -> Dict[str, str]
	"""Map sequential keys to fn(path)."""
	return {str(i): fn(path) for i, path in enumerate(paths)}

def _copy_as_list(value):
	# type: (Iterable[str]) -> List[str]
	"""Copy an iterable of strings as a list. Fail fast on improper input."""
	assert isinstance(value, Iterable) and not isinstance(value, basestring), (
		'Expected a list, not {}'.format(repr(value)))
	result = list(value)
	for s in result:
		assert isinstance(s, basestring), 'Expected a string, not {}'.format(s)
	return result

def _copy_path_list(value):
	# type: (Iterable[str]) -> List[str]
	"""Copy an iterable of strings as a list and check that they're absolute paths
	to catch goofs like `outputs=plot_dir` or `outputs=['out/metadata.json']`.
	Sisyphus needs absolute paths to mount into the Docker container. Fail fast
	on improper input. Handle the the STDOUT prefix '>'.
	"""
	result = _copy_as_list(value)
	for path in result:
		path = path.lstrip('>')
		assert os.path.isabs(path), 'Expected an absolute path, not {}'.format(path)
	return result

def _launch_workers(worker_names):
	# type: (List[str]) -> None
	"""Launch Sisyphus worker nodes with the given names."""
	path = os.path.join(fp.ROOT_PATH, 'runscripts', 'cloud', 'launch-workers.sh')
	subprocess.call([path] + worker_names)


class Task(object):
	"""A workflow task builder."""

	def __init__(self, name='', image='', command=(),
			inputs=(), outputs=(), storage_prefix='', internal_prefix=''):
		# type: (str, str, Iterable[str], Iterable[str], Iterable[str], str, str) -> None
		"""Construct a Workflow Task.

		input and output paths are internal to the worker container and get
		rebased from internal_prefix to storage_prefix for the storage paths.

		An output path that starts with '>' will capture a log from stdout +
		stderr. The rest of the path will get rebased to a storage path.
		"""
		assert name, 'Every task needs a name'
		assert image, 'Every task needs a Docker image name'
		assert command, 'Every task needs a command list of tokens'
		assert storage_prefix, 'Every task needs a storage_prefix'
		assert internal_prefix, 'Every task needs an internal_prefix'

		self.name = name
		self.image = image
		self.command = _copy_as_list(command)
		self.inputs = _copy_path_list(inputs)
		self.outputs = _copy_path_list(outputs)
		self.storage_prefix = storage_prefix
		self.internal_prefix = internal_prefix

	def build_command(self):
		# type: () -> Dict[str, Any]
		"""Build a Gaia Command to run this Task."""
		def specialize(path):
			# type: (str) -> str
			return STDOUT_PATH if path.startswith('>') else path

		return dict(
			name=self.name,
			image=self.image,
			command=self.command,
			inputs=_keyify(self.inputs),
			outputs=_keyify(self.outputs, specialize),
			vars={})

	def build_step(self):
		# type: () -> Dict[str, Any]
		"""Build a Gaia Step to run this Task."""
		def rebase(path):
			# type: (str) -> str
			return _rebase(path.lstrip('>'), self.internal_prefix, self.storage_prefix)

		return dict(
			name=self.name,
			command=self.name,
			inputs=_keyify(self.inputs, rebase),
			outputs=_keyify(self.outputs, rebase))


class Workflow(object):
	"""A workflow builder."""

	def __init__(self, name, verbose_logging=True):
		# type: (str, bool) -> None
		self.name = name
		self.verbose_logging = verbose_logging
		self._tasks = OrderedDict()  # type: Dict[str, Task]

	def log_info(self, message):
		# type: (str) -> None
		if self.verbose_logging:
			print(message)

	def get(self, task_name, default=None):
		# type: (str, Optional[Task]) -> Optional[Task]
		return self._tasks.get(task_name, default)

	def __getitem__(self, task_name):
		# type: (str) -> Task
		return self._tasks[task_name]

	def add_task(self, task):
		# type: (Task) -> Task
		"""Add a Task object. It creates a workflow step. Return it for chaining."""
		if task.name in self._tasks:
			print('Warning: Replacing the task named "{}"'.format(task.name))

		self._tasks[task.name] = task
		self.log_info('    Added step: {}'.format(task.name))
		return task

	def build_commands(self):
		# type: () -> List[dict]
		"""Build this workflow's Commands."""
		return [task.build_command() for task in self._tasks.itervalues()]

	def build_steps(self):
		# type: () -> List[dict]
		"""Build this workflow's Steps."""
		return [task.build_step() for task in self._tasks.itervalues()]

	def write(self):
		# type: () -> None
		"""Build the workflow and write it as JSON files for debugging that can
		be manually sent to the Gaia server.
		"""
		commands = self.build_commands()
		steps = self.build_steps()

		fp.makedirs('out')
		commands_path = os.path.join('out', 'workflow-commands.json')
		steps_path = os.path.join('out', 'workflow-steps.json')

		self.log_info('\nWriting {} {}'.format(commands_path, steps_path))
		fp.write_json_file(commands_path, commands)
		fp.write_json_file(steps_path, steps)

	def launch_workers(self, count):
		# type: (int) -> None
		"""Launch the requested number of Sisyphus worker nodes (GCE VMs)."""
		if count <= 0:
			return
		count = min(count, MAX_WORKERS)

		self.log_info('\nLaunching {} worker node(s).'.format(count))
		user = os.environ['USER']
		names = ['sisyphus-{}-{}'.format(user, i) for i in range(count)]
		_launch_workers(names)

	def send(self, worker_count=4):
		# type: (int) -> None
		"""Build the workflow and send it to the Gaia server to start running."""
		self.launch_workers(worker_count)

		commands = self.build_commands()
		steps = self.build_steps()

		gaia = Gaia(GAIA_CONFIG)

		try:
			self.log_info('\nUploading {} steps to Gaia for workflow {}'.format(
				len(steps), self.name))
			gaia.command(self.name, commands)
			gaia.merge(self.name, steps)
		except ConnectionError as e:
			print('\n*** Did you set up port forwarding to the gaia host? See'
				  ' runscripts/cloud/ssh-tunnel.sh ***\n')
			raise e
