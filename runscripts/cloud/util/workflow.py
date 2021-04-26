"""Google Cloud workflow builder."""

# TODO(jerry): For Windows: This code uses posixpath to construct paths to use
#  on linux servers (and os.path for local file I/O).
#  To finish the job, either make callers do likewise or replace os.sep with
#  posixpath.sep in argument paths, deal with isabs(), and test out on Windows.

from collections import defaultdict, OrderedDict
import os
import posixpath
import re
import subprocess
from typing import Any, Callable, Dict, Iterable, List, Optional, Set
import urllib.parse

from borealis import gce
from borealis.util import gcp
from borealis.docker_task import DockerTask
from fireworks import FiretaskBase, Firework, LaunchPad
from fireworks import Workflow as FwWorkflow
from future.utils import raise_with_traceback
import ruamel.yaml as yaml

from wholecell.utils import filepath as fp
from wholecell.utils.py3 import ANY_STRING


STDOUT_PATH = '>'    # special path that captures stdout + stderror
LOG_OUT_PATH = '>>'  # special path for a fuller log; written even on task failure

STORAGE_ROOT_ENV_VAR = 'WORKFLOW_STORAGE_ROOT'

# The MongoDB service must be reachable at the host:port named in the LaunchPad
# file, either directly or via an ssh port forwarding tunnel like
# runscripts/cloud/mongo-ssh.sh
DEFAULT_LPAD_YAML = 'my_launchpad.yaml'
DEFAULT_FIREWORKS_DATABASE = 'fireworks'  # in LaunchPad's constructor but it's broken in uri_mode
DEFAULT_FIREWORKER_LAUNCHPAD_HOST = 'mongo2'


def _keyify(paths, fn=lambda path: path):
	# type: (Iterable[str], Callable[[str], str]) -> Dict[str, str]
	"""Map sequential keys to fn(path)."""
	return {str(i): fn(path) for i, path in enumerate(paths)}

def _copy_as_list(value):
	# type: (Iterable[str]) -> List[str]
	"""Copy an iterable of strings as a list. Fail fast on improper input."""

	assert isinstance(value, Iterable) and not isinstance(value, ANY_STRING), (
		'Expected a list-like Iterable, not {!r}'.format(value))
	result = list(value)
	for s in result:
		assert isinstance(s, ANY_STRING), 'Expected a string, not {!r}'.format(s)
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
		if path.startswith('>'):
			assert is_output, (
				'''Input paths must not start with '>': "{}"'''.format(path))

			assert os.path.basename(path), (
				'A capture path must name a file, not a directory: "{}"'.format(path))

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
	bucket_name = 'sisyphus-' + os.environ['USER']
	region = gcp.gcloud_get_config('compute/region')
	return (f'{prefix}'
			f' Create your own Google Cloud Storage bucket to aid usage'
			f' tracking, cleanup, and ACLs.'
			f' Pick a name like "{bucket_name}", THE SAME REGION ({region})'
			f' used with Compute Engine (run `gcloud info` for more info),'
			f' Standard storage'
			f' class, and default access control. The name has to be globally'
			f' unique across EVERY GCS PROJECT IN THE WORLD so we\'re using'
			f' "sisyphus-" as a distinguishing prefix.'
			f' BEWARE that the name is publicly visible so don\'t include login'
			f' IDs, email addresses, project names, project numbers, or'
			f' personally identifiable information (PII).'
			f' Then store the bucket name in an environment variable in your'
			f' shell profile and update your current shell, e.g.'
			f' `export {STORAGE_ROOT_ENV_VAR}="{bucket_name}"`.')

def _sanitize_uri(host: str, uri_mode: Optional[bool] = True) -> str:
	"""Sanitize password, params, query, and fragment out of the Launchpad 'host'."""
	if uri_mode:
		p1 = urllib.parse.urlparse(host)
		p2 = (p1.scheme, p1.netloc.split(':')[-1].split('@')[-1], p1.path, '', '', '')
		return urllib.parse.urlunparse(p2)

	return host

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
		encoding='utf-8',
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
		"""Construct a "DockerTask" Workflow Firetask.

		Fireworks can run other Firetasks but DockerTask is especially handy for
		Google Cloud workflows since it pulls a Docker Image containing the
		payload code, copies input and output files from/to Google Cloud Storage
		(GCS), implements timeouts, and sends console output to cloud logs.

		The Workflow builder uses `inputs` and `outputs` to derive Task-to-Task
		dependencies.

		Args:
			name: The task name. It must be unique within the workflow.
			image: The Docker Image name to load as a Container.
			command: The shell command to run in the Container.
			inputs: The absolute pathnames (internal to the Docker Container)
				for the input files and directories to fetch from GCS.

				Their corresponding GCS storage paths will get constructed by
				rebasing each path from internal_prefix to storage_prefix.

				Each path indicates a directory tree of files (rather than a
				single file) iff it ends with '/'.
			outputs: The absolute pathnames (internal to the Docker Container)
				for the output files and directories to store to GCS.

				Output pathnames get rebased and treated as directories like
				input pathnames.

				Outputs will get written to GCS if the task completes normally.

				'>': An output path that starts with '>' will capture stdout +
				stderr (if the task completes normally). The rest of the path
				gets rebased to provide the storage pathname.

				'>>': An output path that starts with '>>' will capture a log of
				stdout + stderr + other log messages like elapsed time and the
				task exit code, and it gets written to GCS even if the task
				fails. This is useful for debugging. See `store_log`.
			storage_prefix: The GCS path prefix for the inputs and outputs.
			internal_prefix: The path prefix within the container that
				corresponds to storage_prefix.
			timeout: How many seconds to let the Docker task run.
				0 => DEFAULT_TIMEOUT.
			store_log: If True, save a '>>' log file to the logs/ dir.
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


class Workflow(object):
	"""A workflow builder."""

	def __init__(self, name, owner_id='', verbose_logging=True, description=''):
		# type: (str, str, bool, str) -> None
		"""
		Construct a workflow builder, ready to add Tasks.

		Args:
			name: Standard practice is to construct a name in
				the form `OWNER_PROGRAM_DATETIME`, e.g.
				`crick_DemoWorkflow_20191209.133734` to aid sorting and
				filtering workflows.
				See filepath.timestamp() for the datetime format.

			owner_id: Sets the `owner` name for this workflow to enable
				filtering the workflow list by owner. E.g. the user name or a
				CI build name. Defaults to the $USER environment variable.

			verbose_logging: Enables verbose logging while building a workflow.

			description: Optional description for this workflow. This gets stored
				in a FireWorks property and also used in the timestamp part of
				the GCS storage path (spaces replaced with underscores).
		"""
		self.owner_id = owner_id or os.environ['USER']
		self.name = name
		self.properties = {'owner': self.owner_id}
		self.verbose_logging = verbose_logging
		self._tasks = OrderedDict()  # type: Dict[str, Task]
		self._output_to_taskname = {}  # type: Dict[str, str]
		self._input_to_tasknames = defaultdict(set)  # type: Dict[str, Set[str]]

		if description:
			self.add_properties(description=description)

	@classmethod
	def sanitize_description(cls, description):
		# type (str) -> str
		"""Sanitize the description and check that it's legal in a file path."""
		description = description.replace(' ', '_')

		pattern = r'[-.\w]*$'
		assert re.match(pattern, description), (
			"description {!r} doesn't match the regex pattern {!r} for a file path."
				.format(description, pattern))

		return description

	@classmethod
	def timestamped_description(cls, timestamp, description=''):
		# type (str, str) -> str
		"""Construct a timestamped workflow description subdirectory name."""
		return timestamp + (
			'__' + cls.sanitize_description(description) if description else '')

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
		"""Add more properties to the workflow, e.g. 'description'."""
		self.properties.update(kwargs)

	def add_task(self, task):
		# type: (Task) -> Task
		"""Add a Task object. It will create a FireWorks "Firework" that runs a
		single "Firetask". Return it for chaining.
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
			print(f'\nWARNING: Task {task.name} has inputs unfulfilled by the'
				  f' Tasks in this workflow:')
			for input_path in sorted(unfulfilled):
				print(f'    {input_path}')
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
		"""Build all the FireWorks `Firework` objects for the workflow."""
		built = OrderedDict()  # type: Dict[str, Firework]

		for task in self._tasks.values():
			self._build_firework(task, built)

		return list(built.values())

	def build_workflow(self):
		# type: () -> FwWorkflow
		"""Build the FireWorks `Workflow` object."""
		fireworks = self.build_fireworks()
		wf = FwWorkflow(fireworks, name=self.name, metadata=self.properties)
		return wf

	def launch_fireworkers(self, count, config, gce_options=None):
		# type: (int, Dict[str, Any], Optional[Dict[str, Any]]) -> None
		"""Launch the requested number of fireworker nodes (GCE VMs) with
		(host, uri_mode, db name, username, password) metadata from config and
		GCE VM options from gce_options.
		"""
		def copy_key(src, key, dest):
			# type: (dict, str, dict) -> None
			"""Copy a keyed value from src to dest dicts unless the value is
			absent or None.
			"""
			val = src.get(key)
			if val is not None:
				dest[key] = val

		db_name = config['name']
		prefix = 'fireworker-{}'.format(db_name)
		options = {
			'image-family': 'fireworker',
			'network-interface': 'no-address',  # no External IP
			'description': f'FireWorks worker VM for user/ID {self.name}'}
		options.update(gce_options or {})

		# Tell the Fireworkers how to contact the Launchpad DB.
		metadata = {'db': db_name}
		copy_key(config, 'host', metadata)
		copy_key(config, 'uri_mode', metadata)
		copy_key(config, 'username', metadata)
		copy_key(config, 'password', metadata)

		info = ''
		if "host" in metadata:
			db_host = _sanitize_uri(
				metadata['host'], uri_mode=metadata.get('uri_mode'))
			info = f' to connect to MongoDB at {db_host}'
		self.log_info(f'\nCreating {count} GCE Fireworker '
					  f'{"VM" if count == 1 else "VMs"}{info}')

		engine = gce.ComputeEngine(prefix)
		engine.create(count=count, command_options=options, **metadata)

	def send_to_lpad(self, worker_count=4, lpad_filename=DEFAULT_LPAD_YAML,
					 gce_options=None):
		# type: (int, str, Optional[Dict[str, Any]]) -> FwWorkflow
		"""Build this workflow for FireWorks, upload it to the given or
		default LaunchPad, launch workers with the given GCE VM options, and
		return the built workflow.

		ASSUMES: When uploading to Launchpad host=localhost (and not in
		uri_mode) and launching GCE Fireworkers, they'll need to know to
		connect to the host without the tunnel, so assume it's our GCE MongoDB
		host. This matters only to GCE Fireworkers, not when you'll run
		`fireworker` locally. If this assumption becomes an issue, use uri_mode
		in the YAML file or add a field for the GCE-relative host.
		"""
		with open(lpad_filename) as f:
			yml = yaml.YAML(typ='safe')
			config = yml.load(f)
			lpad = LaunchPad(**config)

		# (See the "ASSUMES" comment above about 'localhost'.)
		config['host'] = (DEFAULT_FIREWORKER_LAUNCHPAD_HOST
						  if lpad.host == 'localhost' else lpad.host)
		if lpad.uri_mode:
			config['name'] = config['host'].split('/')[-1].split('?')[0]
		else:
			config['name'] = lpad.name or DEFAULT_FIREWORKS_DATABASE

		wf = self.build_workflow()

		try:
			self.log_info(
				f'\nSending the workflow to the Launchpad DB {config["name"]},'
				f' connecting to MongoDB at'
				f' {_sanitize_uri(lpad.host, lpad.uri_mode)}')
			lpad.add_wf(wf)
		except ValueError as e:
			raise_with_traceback(RuntimeError(
				'You might need to set up port forwarding to the MongoDB'
				' server. See `runscripts/cloud/mongo-ssh.sh`.\n'
				f'Caused by: {e!r}'))

		# Launch the workers after the successful upload.
		try:
			self.launch_fireworkers(worker_count, config, gce_options=gce_options)
		except Exception:
			print('\nNOTE: The workflow was uploaded to the Launchpad DB but'
				  ' creating Fireworker VMs raised an exception. You could use'
				  ' the `gce` command to create Fireworkers, or maybe run'
				  ' `lpad reset` to erase the workflow and start over.\n'
				  '    |\n'
				  '    V')
			raise

		return wf

	def write(self):
		# type: () -> None
		"""Write this workflow as a YAML file for FireWorks."""
		fp.makedirs('out')
		filename = os.path.join('out', 'workflow-{}.yaml'.format(self.name))
		self.log_info('\nWriting workflow to {}'.format(filename))

		fw_wf = self.build_workflow()

		with open(filename, 'w') as f:
			yml = yaml.YAML(typ='safe')
			yml.dump(fw_wf.to_dict(), f)
