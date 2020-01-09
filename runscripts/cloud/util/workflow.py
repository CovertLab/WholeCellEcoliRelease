"""Google Cloud workflow builder."""

# TODO(jerry): For Windows: This code uses posixpath to construct paths to use
#  on linux servers (and os.path for local file I/O).
#  To finish the job, either make callers do likewise or replace os.sep with
#  posixpath.sep in argument paths, deal with isabs(), and test out on Windows.

from __future__ import absolute_import, division, print_function

from collections import defaultdict, OrderedDict
import os
import posixpath
import re
import sys
if os.name == 'posix' and sys.version_info[0] < 3:
	import subprocess32 as subprocess
else:
	import subprocess
from typing import Any, Callable, Dict, Iterable, List, Optional, Set

from fireworks import FiretaskBase, Firework, LaunchPad
from fireworks import Workflow as FwWorkflow
from gaia.client import Gaia
from requests import ConnectionError
import yaml

from wholecell.utils import filepath as fp
from cloud.docker_task import DockerTask
from .. import gce_vms


# Config details to pass to Gaia.
# ASSUMES: gaia_host is reachable e.g. via an ssh tunnel set up by
# runscripts/cloud/ssh-tunnel.sh.
GAIA_CONFIG = {'gaia_host': 'localhost:24442'}

STDOUT_PATH = '>'    # special path that captures stdout + stderror
LOG_OUT_PATH = '>>'  # special path for a fuller log; written even on task failure

STORAGE_ROOT_ENV_VAR = 'WORKFLOW_STORAGE_ROOT'

LAUNCHPAD_FILENAME = 'my_launchpad.yaml'


def _rebase_for_gaia(path, internal_prefix, storage_prefix):
	# type: (str, str, str) -> str
	"""Return a path rebased from internal_prefix to storage_prefix and switch
	to "bucket:path" format for Gaia/Sisyphus."""
	relpath = posixpath.relpath(path, internal_prefix)
	new_path = posixpath.join(storage_prefix, relpath).replace('/', ':', 1)

	# posixpath.relpath() removes a trailing slash if it exists. Restore it.
	if path.endswith(posixpath.sep):
		new_path = posixpath.join(new_path, '')
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

def _copy_path_list(value, internal_prefix, is_output=False):
	# type: (Iterable[str], str, bool) -> List[str]
	"""Copy an iterable of strings as a list and check that they're absolute
	paths that start with `internal_prefix` to catch improper args like
	`outputs=plot_dir` or `outputs=['out/metadata.json']`.
	(The task worker needs absolute paths to mount into the Docker container.)
	Handle STDOUT_PATH and LOG_OUT_PATH prefixes.
	"""
	result = _copy_as_list(value)
	for path in result:
		assert is_output or not path.startswith('>'), (
			'''Input paths must not start with '>', not "{}"'''.format(path))

		path = path.lstrip('>')
		assert posixpath.isabs(path), (
			'Expected an absolute path, not "{}"'.format(path))

		relpath = posixpath.relpath(path, internal_prefix)
		assert '..' not in relpath, (
			'Expected a path that starts with the internal_prefix "{}", not "{}"'
				.format(internal_prefix, path))
	return result

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
		Docker container. The corresponding storage paths will get constructed
		by rebasing each path from internal_prefix to storage_prefix.

		Each path indicates a directory tree of files to fetch or store (rather
		than a single file) iff it ends with '/'.

		Outputs will get written to GCS if the task completes normally. An
		output path that starts with '>' will capture stdout + stderr (if the
		task completes normally), while the rest of the path gets rebased to
		provide the storage path.

		An output path that starts with '>>' will capture a log of stdout +
		stderr + other log messages like elapsed time and task exit code, even
		if the task fails. This is useful for debugging. Just set
		store_log=True (the default) to save this log.
		"""
		assert name, 'Every task needs a name'
		assert image, 'Every task needs a Docker image name'
		assert command, 'Every task needs a command list of tokens'
		assert storage_prefix, 'Every task needs a storage_prefix'
		assert not storage_prefix.startswith('/'), (
			'storage_prefix must not start with "/": "{}"'.format(storage_prefix))
		assert internal_prefix, 'Every task needs an internal_prefix'
		assert posixpath.isabs(internal_prefix), (
			'internal_prefix must be an absolute path, not "{}"'.format(internal_prefix))

		self.name = name
		self.image = image
		self.command = _copy_as_list(command)
		self.inputs = _copy_path_list(inputs, internal_prefix)
		self.outputs = _copy_path_list(outputs, internal_prefix, is_output=True)
		self.storage_prefix = storage_prefix
		self.internal_prefix = internal_prefix
		self.timeout = timeout if timeout > 0 else self.DEFAULT_TIMEOUT

		if store_log:
			self.outputs[0:0] = [
				LOG_OUT_PATH + posixpath.join(internal_prefix, 'logs', name + '.log')]

	def __repr__(self):
		return 'Task{}'.format(vars(self))

	def build_firetask(self):
		# type: () -> FiretaskBase
		"""Build a FireWorks Firetask to run this Task."""
		return DockerTask(
			name=self.name,
			image=self.image,
			command=self.command,
			inputs=self.inputs,
			outputs=self.outputs,
			storage_prefix=self.storage_prefix,
			internal_prefix=self.internal_prefix,
			timeout=self.timeout)

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
			return _rebase_for_gaia(path.lstrip('>'), self.internal_prefix, self.storage_prefix)

		return dict(
			name=self.name,
			command=self.name,
			inputs=_keyify(self.inputs, rebase),
			outputs=_keyify(self.outputs, rebase),
			timeout=self.timeout)


class Workflow(object):
	"""A workflow builder."""

	def __init__(self, name, owner_id='', verbose_logging=True):
		# type: (str, str, bool) -> None
		"""
		Construct a workflow builder, ready to add Tasks.

		Args:
			name: Each workflow needs a name that's unique among all workflows
				in the Gaia server. Standard practice is to construct a name in
				the form `OWNER_PROGRAM_DATETIME`, e.g.
				`crick_DemoWorkflow_20191209.133734` to aid sorting and
				filtering workflows.
				See filepath.timestamp() for the datetime format.

			owner_id: Sets the `owner` name for this workflow to enable
				filtering the workflow list by owner. E.g. the user name or a
				CI build name. Defaults to the $USER environment variable.

			verbose_logging: Enables verbose logging while building a workflow.
		"""
		self.owner_id = owner_id or os.environ['USER']
		self.name = name
		self.properties = {'owner': self.owner_id}
		self.verbose_logging = verbose_logging
		self._tasks = OrderedDict()  # type: Dict[str, Task]
		self._output_to_taskname = {}  # type: Dict[str, str]
		self._input_to_tasknames = defaultdict(set)  # type: Dict[str, Set[str]]

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

	def __repr__(self):
		return 'Workflow{}'.format(vars(self))

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

	def add_properties(self, **kwargs):
		# type: (Any) -> None
		"""Add more Gaia properties to the workflow, e.g. 'description'."""
		self.properties.update(kwargs)

	def add_task(self, task):
		# type: (Task) -> Task
		"""Add a Task object. It will create a workflow step, aka a FireWorks
		"firework".
		Return it for chaining.
		Raise ValueError if there's a conflicting task name or output path.
		"""
		task_name = task.name
		if task_name in self._tasks:
			raise ValueError('''There's already a task named "{}"'''.format(
				task_name))

		for output_path in task.outputs:
			real_path = output_path.lstrip('>')
			existing_task = self._output_to_taskname.get(real_path)
			if existing_task:
				raise ValueError(
					'Task "{}" wants to write to the same output "{}" that Task'
					' "{}" writes to'.format(task_name, real_path, existing_task))

		# Checks passed. Now add the new task.
		self._tasks[task_name] = task
		for output_path in task.outputs:
			self._output_to_taskname[output_path.lstrip('>')] = task_name
		for input_path in task.inputs:
			self._input_to_tasknames[input_path].add(task_name)

		self.log_info('    Added task: {}'.format(task_name))
		return task

	def task_dependencies(self, task):
		# type: (Task) -> Set[str]
		"""Given a Task, return a set of the task names it depends on, i.e.
		that write its inputs. Print a warning about any unfulfilled inputs.
		"""
		dependencies = set()
		unfulfilled = set()

		for input_path in task.inputs:
			dependency = self._output_to_taskname.get(input_path)
			if dependency:
				dependencies.add(dependency)
			else:
				unfulfilled.add(input_path)

		if unfulfilled:
			print('WARNING: Task "{}" has inputs unfulfilled by the Tasks in'
				  ' this workflow: {}'.format(task.name, sorted(unfulfilled)))
		return dependencies

	def task_dependents(self, task):
		# type: (Task) -> Set[str]
		"""Given a Task, return a set of the task names that depend on it, i.e.
		that read its outputs.
		"""
		dependents = set()

		for output_path in task.outputs:
			dependents.update(self._input_to_tasknames[output_path])

		return dependents

	def _build_firework(self, task, built):
		# type: (Task, Dict[str, Optional[Firework]]) -> Firework
		"""Build a Firework for the given Task or return the already-built one,
		recursively building its parent Fireworks (dependencies) as needed,
		collecting them in topo-sorted `built` to ensure they're 1:1 and acyclic.
		"""
		# TODO(jerry): Optional optimization. Given Task links A -> B -> C and
		#  A -> C, omitting lots of parent links like A -> C can speed up the
		#  workflow engine (assuming nobody deletes tasks). Optimize down to
		#  {the set of direct parents} - {transitive closure of their parents}?
		#  Or safer, would it suffice to sort dependencies to begin with the
		#  later ones in the pipeline?
		task_name = task.name
		already_built = built.get(task_name, False)
		if already_built:
			return already_built

		if already_built is None:
			raise ValueError('Cyclic dependency with Task "{}"'.format(task_name))
		built[task_name] = None

		parents = []  # type: List[Firework]
		dependency_names = self.task_dependencies(task)

		for dependency_name in dependency_names:
			parent_task = self[dependency_name]
			parent_firework = self._build_firework(parent_task, built)
			parents.append(parent_firework)

		firework = Firework(
			task.build_firetask(),
			name=task_name,
			parents=parents)  # TODO(jerry): spec=...?
		built[task_name] = firework

		return firework

	def build_fireworks(self):
		# type: () -> List[Firework]
		"""Build this workflow's Firework objects."""
		built = OrderedDict()

		for task in self._tasks.itervalues():
			self._build_firework(task, built)

		return built.values()

	def build_workflow(self):
		# type: () -> FwWorkflow
		"""Build this workflow for FireWorks."""
		fireworks = self.build_fireworks()
		wf = FwWorkflow(fireworks, name=self.name, metadata=self.properties)
		return wf

	def launch_fireworkers(self, count):
		# type: (int) -> None
		"""Launch the requested number of fireworker nodes (GCE VMs)."""
		prefix = 'fireworker-{}'.format(self.name)
		gce_vms.launch_fireworkers(prefix, count=count, workflow=self.name)

	def send_to_lpad(self, lpad_filename=LAUNCHPAD_FILENAME, worker_count=4):
		# type: (str, int) -> FwWorkflow
		"""Build this workflow for FireWorks, upload it to the given or
		default LaunchPad, launch workers, and return the built workflow.
		"""
		with open(lpad_filename) as f:
			lpad = LaunchPad(**yaml.safe_load(f))

		wf = self.build_workflow()
		lpad.add_wf(wf)

		# Launch the workers after the successful upload.
		self.launch_fireworkers(worker_count)

		return wf

	def build_gaia_commands(self):
		# type: () -> List[dict]
		"""Build this workflow's Gaia Commands."""
		return [task.build_command() for task in self._tasks.itervalues()]

	def build_gaia_steps(self):
		# type: () -> List[dict]
		"""Build this workflow's Gaia Steps."""
		return [task.build_step() for task in self._tasks.itervalues()]

	def write_for_gaia(self):
		# type: () -> None
		"""Build the Gaia workflow spec and write it as JSON files for
		debugging that can be manually sent to the Gaia server.
		"""
		commands = self.build_gaia_commands()
		steps = self.build_gaia_steps()

		fp.makedirs('out')
		commands_path = os.path.join('out', 'workflow-commands.json')
		steps_path = os.path.join('out', 'workflow-steps.json')

		self.log_info('\nWorkflow {} writing {} {}'.format(
			self.name, commands_path, steps_path))
		fp.write_json_file(commands_path, commands)
		fp.write_json_file(steps_path, steps)

		self.log_info('Associated properties = {}'.format(self.properties))

	def launch_sisyphus_workers(self, count):
		# type: (int) -> None
		"""Launch the requested number of sisyphus worker nodes (GCE VMs)."""
		prefix = 'sisyphus-{}'.format(self.name)
		gce_vms.launch_sisyphus_workers(prefix, count=count, workflow=self.name)

	def send_to_gaia(self, worker_count=4):
		# type: (int) -> None
		"""Build the Gaia workflow and send it to the server to start running."""
		commands = self.build_gaia_commands()
		steps = self.build_gaia_steps()
		self.properties['requested-worker-count'] = worker_count

		gaia = Gaia(GAIA_CONFIG)

		try:
			self.log_info('\nUploading workflow {} to Gaia with {} steps'.format(
				self.name, len(steps)))
			gaia.upload(
				workflow=self.name,
				properties=self.properties,
				commands=commands,
				steps=steps)
		except ConnectionError as e:
			print('\n*** Did you set up port forwarding to the gaia host? See'
				  ' runscripts/cloud/ssh-tunnel.sh ***\n')
			raise e

		# Launch the workers after the successful upload.
		self.launch_sisyphus_workers(worker_count)
