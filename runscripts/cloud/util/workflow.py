"""Generic Sisyphus/Gaia/Google Cloud workflow builder."""

# TODO(jerry): For Windows: This code uses posixpath to construct paths to use
#  on linux servers (and os.path for local file I/O).
#  To finish the job, either make callers do likewise or replace os.sep with
#  posixpath.sep in argument paths, deal with isabs(), and test out on Windows.

from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import os
import posixpath
import re
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

STDOUT_PATH = '>'    # special path that captures stdout + stderror
LOG_OUT_PATH = '>>'  # special path for a fuller log; written even on task failure

STORAGE_ROOT_ENV_VAR = 'WORKFLOW_STORAGE_ROOT'
MAX_WORKERS = 500  # don't launch more than this many worker nodes at a time


def _rebase(path, internal_prefix, storage_prefix):
	# type: (str, str, str) -> str
	"""Return a path rebased from internal_prefix to storage_prefix and switch
	to "bucket:path" format for Gaia/Sisyphus."""
	relpath = posixpath.relpath(path, internal_prefix)
	new_path = posixpath.join(storage_prefix, relpath).replace('/', ':', 1)

	# posixpath.relpath removes a trailing slash if it exists.
	if path.endswith(posixpath.sep):
		new_path = posixpath.join(new_path, '')

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
		assert posixpath.isabs(path), 'Expected an absolute path, not {}'.format(path)
	return result

def _launch_workers(worker_names, workflow=''):
	# type: (List[str], str) -> None
	"""Launch Sisyphus worker nodes with the given names."""
	path = os.path.join(fp.ROOT_PATH, 'runscripts', 'cloud', 'launch-workers.sh')
	subprocess.call([path] + worker_names,
		env=dict(os.environ, WORKFLOW=workflow))

def _to_create_bucket(prefix):
	# type: (str) -> str
	return ('{0}'
			' Create your own Google Cloud Storage bucket if needed. This'
			' supports usage tracking, cleanup, and ACLs.'
			' Pick a name like "sisyphus-{2}" [it has to be globally unique;'
			' BEWARE that it\'s publicly visible so don\'t include login IDs,'
			' email addresses, project names, project numbers, or personally'
			' identifiable information (PII)], the same Region used with'
			' Compute Engine (run `gcloud info` for info), Standard storage'
			' class, and default access control.'
			' Then store it in an environment variable in your shell profile'
			' and update your current shell, e.g. `export {1}="sisyphus-{2}"`.'
				.format(prefix, STORAGE_ROOT_ENV_VAR, os.environ['USER']))

def bucket_from_path(storage_path):
	# type: (str) -> str
	"""Extract the bucket name (the first component) from a Google Cloud
	Storage path such as 'sisyphus-crick', 'sisyphus-crick/', or
	'sisyphus/data/crick'.
	"""
	return storage_path.split('/', 1)[0]

def validate_gcs_bucket(bucket_name):
	# type: (str) -> None
	"""Raise an exception if the Google Cloud Storage bucket name is malformed,
	doesn't exist, or inaccessible.
	"""
	pattern = r'[a-z0-9][-_a-z0-9]{1,61}[a-z0-9]$'  # no uppercase; dots require DNS approval
	if not re.match(pattern, bucket_name):
		raise ValueError('Storage bucket name "{}" doesn\'t match the required'
						 ' regex pattern "{}".'.format(bucket_name, pattern))

	completion = subprocess.run(
		["gsutil", "ls", "-b", "gs://" + bucket_name],  # bucket list!
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		universal_newlines=True,
		timeout=60)
	if completion.returncode:
		message = _to_create_bucket(
			'Couldn\'t access the Google Cloud Storage bucket "{0}" [{1}].'
				.format(bucket_name, completion.stderr.strip()))
		raise ValueError(message)
	# else completion.stdout.strip() should == 'gs://' + bucket_name + '/'


class Task(object):
	"""A workflow task builder."""

	DEFAULT_TIMEOUT = 60 * 60  # in seconds

	def __init__(self, name='', image='', command=(),
			inputs=(), outputs=(), storage_prefix='', internal_prefix='',
			timeout=0, store_log=True):
		# type: (str, str, Iterable[str], Iterable[str], Iterable[str], str, str, int, bool) -> None
		"""Construct a Workflow Task.

		The inputs and outputs are absolute paths internal to the worker's
		Docker container.
		Task will rebase them from internal_prefix to storage_prefix to construct
		the corresponding storage paths. An internal path ending with '/' will
		upload or download a directory tree.

		An output path that starts with '>' will capture stdout + stderr (if the
		task completes normally). The rest of the path will get rebased to a
		storage path.

		(An output path that starts with '>>' will capture a log of stdout +
		stderr + other log messages like elapsed time and task exit code, for
		for debugging, even if the task fails. Just default store_log=True to
		save a log.)
		"""
		assert name, 'Every task needs a name'
		assert image, 'Every task needs a Docker image name'
		assert command, 'Every task needs a command list of tokens'
		assert storage_prefix, 'Every task needs a storage_prefix'
		assert not storage_prefix.startswith('/'), (
			'storage_prefix must not start with "/": {}'.format(storage_prefix))
		assert internal_prefix, 'Every task needs an internal_prefix'
		assert posixpath.isabs(internal_prefix), (
			'internal_prefix must be an absolute path, not {}'.format(internal_prefix))

		self.name = name
		self.image = image
		self.command = _copy_as_list(command)
		self.inputs = _copy_path_list(inputs)
		self.outputs = _copy_path_list(outputs)
		self.storage_prefix = storage_prefix
		self.internal_prefix = internal_prefix
		self.timeout = timeout if timeout > 0 else self.DEFAULT_TIMEOUT

		if store_log:
			self.outputs[0:0] = [
				LOG_OUT_PATH + posixpath.join(internal_prefix, 'logs', name + '.log')]

	def build_command(self):
		# type: () -> Dict[str, Any]
		"""Build a Gaia Command to run this Task."""
		def specialize(path):
			# type: (str) -> str
			return (LOG_OUT_PATH if path.startswith(LOG_OUT_PATH)
					else STDOUT_PATH if path.startswith(STDOUT_PATH)
					else path)

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
			outputs=_keyify(self.outputs, rebase),
			timeout=self.timeout)


class Workflow(object):
	"""A workflow builder."""

	def __init__(self, name, verbose_logging=True):
		# type: (str, bool) -> None
		self.name = name
		self.verbose_logging = verbose_logging
		self._tasks = OrderedDict()  # type: Dict[str, Task]

	@classmethod
	def storage_root(cls, cli_arg=None):
		# type: (Optional[str]) -> str
		"""Validate the workflow cloud storage root from the given CLI argument
		or else an environment variable that's usually configured via shell
		profile. If invalid, raise an exception with instructions to set it up.

		BEST PRACTICE is a storage bucket per user like 'sisyphus-crick'. That
		supports usage tracking, cleanup, and ACLs. But as a fallback, use a
		subdirectory of the 'sisyphus' bucket, e.g. 'sisyphus/data/crick'.
		"""
		try:
			root = cli_arg or os.environ[STORAGE_ROOT_ENV_VAR]
		except KeyError:
			message = _to_create_bucket(
				'Environment variable ${} not found.'.format(STORAGE_ROOT_ENV_VAR))
			raise KeyError(STORAGE_ROOT_ENV_VAR, message)

		bucket = bucket_from_path(root)
		assert bucket, _to_create_bucket(
			'${} "{}" doesn\'t begin with a storage bucket name.'.format(
				STORAGE_ROOT_ENV_VAR, root))
		validate_gcs_bucket(bucket)

		return root

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

		self.log_info('\nWorkflow {} writing {} {}'.format(
			self.name, commands_path, steps_path))
		fp.write_json_file(commands_path, commands)
		fp.write_json_file(steps_path, steps)

	def launch_workers(self, count):
		# type: (int) -> None
		"""Launch the requested number of Sisyphus worker nodes (GCE VMs)."""
		if count <= 0:
			return
		count = min(count, MAX_WORKERS)

		self.log_info('\nLaunching {} worker node(s).'.format(count))

		# Convert the workflow name to valid GCE VM names.
		sanitized = re.sub('[^-a-z0-9]', '-', self.name.lower()).replace('workflow', '')
		names = ['sisyphus-{}-{}'.format(sanitized, i) for i in range(count)]
		_launch_workers(names, workflow=self.name)

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
