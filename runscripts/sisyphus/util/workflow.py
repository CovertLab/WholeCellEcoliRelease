"""Generic Sisyphus/Gaia workflow builder."""

from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import os
import sys
if os.name == 'posix' and sys.version_info[0] < 3:
	import subprocess32 as subprocess
else:
	import subprocess
from typing import Any, Dict, Iterable, List, Optional

from gaia.client import Gaia
from requests import ConnectionError

from wholecell.utils import filepath as fp


# Config details to pass to Gaia.
# ASSUMES: gaia_host is reachable e.g. via an ssh tunnel set up by
# runscripts/sisyphus/ssh-tunnel.sh.
GAIA_CONFIG = {'gaia_host': 'localhost:24442'}


def _rebase(path, old_prefix, new_prefix):
	# type: (str, str, str) -> str
	"""Return a path starting with new_prefix in place of old_prefix."""
	new_path = os.path.join(new_prefix, os.path.relpath(path, old_prefix))

	### Should new_prefix end with os.sep iff path does? relpath() drops the
	### trailing os.sep. This `if` statement will reattach it but disable it
	### until we figure out what to do.
	# if path.endswith(os.sep) and not new_path.endswith(os.sep):
	# 	new_path = os.path.join(new_path, '')

	assert '..' not in new_path, "Rebasing a path that doesn't start with old_prefix?"
	return new_path

def _keyify(paths):
	# type: (Iterable[str]) -> Dict[str, str]
	"""Map sequential keys to the given paths."""
	return {str(i): path for i, path in enumerate(paths)}

def _re_keyify(paths, old_prefix, new_prefix):
	# type: (Iterable[str], str, str) -> Dict[str, str]
	"""Map sequential keys to rebased paths."""
	return _keyify([_rebase(path, old_prefix, new_prefix) for path in paths])

def _copy_as_list(value):
	# type: (Iterable[str]) -> List[str]
	"""Copy an iterable of strings as a list."""
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
	Sisyphus needs absolute paths to mount into the Docker container.
	"""
	result = _copy_as_list(value)
	for path in result:
		assert os.path.isabs(path), 'Expected an absolute path, not {}'.format(path)
	return result

def _launch_workers(worker_names):
	# type: (List[str]) -> None
	"""Launch Sisyphus worker nodes with the given names."""
	path = os.path.join(fp.ROOT_PATH, 'runscripts', 'sisyphus', 'launch-workers.sh')
	subprocess.call([path] + worker_names)


class Task(object):
	"""A workflow task builder."""

	def __init__(self, upstream_tasks=(), name='', image='', command=(),
			inputs=(), outputs=(), storage_prefix='', internal_prefix=''):
		# type: (Iterable[Task], str, str, Iterable[str], Iterable[str], Iterable[str], str, str) -> None
		"""Construct a Workflow Task.
		upstream_tasks and the `>>` operator are convenient ways to add inputs.
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

		for task in upstream_tasks:
			task >> self

	def __rshift__(self, t2):
		# type: (Task) -> Task
		"""Set downstream: `t1 >> t2` adds `t1`'s outputs to `t2`'s inputs.
		Return `t2` for chaining.
		"""
		t2.inputs.extend(self.outputs)
		return t2

	def build_command(self):
		# type: () -> Dict[str, Any]
		"""Build a Gaia Command to run this Task."""
		return dict(
			name=self.name,
			image=self.image,
			command=self.command,
			inputs=_keyify(self.inputs),
			outputs=_keyify(self.outputs),
			vars={})

	def build_step(self):
		# type: () -> Dict[str, Any]
		"""Build a Gaia Step to run this Task."""
		return dict(
			name=self.name,
			command=self.name,
			inputs=_re_keyify(self.inputs, self.internal_prefix, self.storage_prefix),
			outputs=_re_keyify(self.outputs, self.internal_prefix, self.storage_prefix))


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
		"""Add a task object. Return it for chaining."""
		self._tasks[task.name] = task
		self.log_info('    Added task: {}'.format(task.name))
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
		if count <= 0:
			return

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
			print('\n*** Did you set up port forwarding for gaia-base and'
				  ' zookeeper-prime? See runscripts/sisyphus/ssh-tunnel.sh ***\n')
			raise e
