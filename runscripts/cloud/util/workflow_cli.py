"""
Create and upload a simple workflow for simple runs of the workflow software.
"""

import argparse
import os
import posixpath as pp
from typing import Iterable, Optional

import wholecell.utils.filepath as fp
from wholecell.utils import scriptBase
from runscripts.cloud.util.workflow import (DEFAULT_LPAD_YAML,
	STORAGE_ROOT_ENV_VAR, Task, Workflow)


class WorkflowCLI(scriptBase.ScriptBase):
	"""Abstract base class for a Command Line Interface to build a workflow."""

	# Subclasses can override these:
	DOCKER_IMAGE = 'python:3.8.5'
	DEFAULT_TIMEOUT = Task.DEFAULT_TIMEOUT  # in seconds

	def __init__(self, internal_prefix='/tmp'):
		super().__init__()
		self.storage_prefix = ''
		self.internal_prefix = pp.join(pp.sep, internal_prefix, '')
		self.wf = None  # type: Optional[Workflow]

	def internal(self, *path_elements):
		# type: (*str) -> str
		"""Construct a docker container internal file path."""
		return pp.join(self.internal_prefix, *path_elements)

	def add_task(self, name='', inputs=(), outputs=(), command=(), timeout=0):
		# type: (str, Iterable[str], Iterable[str], Iterable[str], int) -> Task
		assert self.wf, 'Need to construct the Workflow object first.'
		task = Task(
			name=name,
			image=self.DOCKER_IMAGE,
			inputs=inputs,
			outputs=outputs,
			storage_prefix=self.storage_prefix,
			internal_prefix=self.internal_prefix,
			command=command,
			timeout=timeout if timeout > 0 else self.DEFAULT_TIMEOUT)
		return self.wf.add_task(task)

	def define_parameters(self, parser):
		"""Define command line parameters used by all workflow builders."""
		self.define_option(parser, 'storage_root', str,
			help='The cloud storage root for the output files, usually a GCS'
				 ' bucket name like "sisyphus-crick". Default = ${}'
				 ' environment variable.'.format(STORAGE_ROOT_ENV_VAR))
		self.define_option(parser, 'launchpad_filename', str, flag='l',
			default=DEFAULT_LPAD_YAML,
			help='Launchpad config YAML filename for sending the workflow to a'
				 ' launchpad DB server.')
		self.define_option(parser, 'workers', int, default=1, flag='w',
			help='The number of worker nodes to launch.')
		parser.add_argument('--dump', action='store_true',
			help='Dump the built workflow to a YAML file for review *instead*'
				 ' of sending it to a launchpad DB server. This is useful'
				 ' for testing and debugging.')

		self.define_parameter_bool(parser, 'verbose', True,
			help='Verbose workflow builder logging.')
		super().define_parameters(parser)

	def define_wf_name_parameters(self, parser):
		"""Define additional command line parameters for wcEcoli-related
		workflow naming: description, (owner) id, timestamp.

		The CLI run() method will use these args to construct a base directory
		name like '20200924.121917__with_Causality' and a workflow name like
		'jerry_WCM_20200924.121917'.

		Call this in a subclass define_parameters() method if relevant.
		"""
		self.define_option(parser, 'description', str, '',
			help='A simulation description; part of the storage folder name.')
		self.define_option(parser, 'id', str, default=None,
			help='Workflow ID or owner ID such as a user name or a CI build'
				 ' name to combine with the timestamp to form the unique'
				 ' workflow name. Default = $WF_ID environment variable or'
				 ' else the $USER environment variable.')
		self.define_option(parser, 'timestamp', str, fp.timestamp(),
			help='Timestamp for this workflow. It gets combined with the'
				 ' Workflow ID to form the workflow name. Set this if you want'
				 ' to upload new steps for an existing workflow. Default ='
				 ' the current local date-time.')

	def build(self, args):
		# type: (argparse.Namespace) -> None
		"""Build the workflow's Firetasks and call add_task() to add each one to
		the workflow."""
		raise NotImplementedError("WorkflowCLI subclass must implement build()")

	def dumpOrRun(self, args):
		# type: (argparse.Namespace) -> None
		"""Dump or run the workflow."""
		assert self.wf, 'Need to construct the Workflow first.'
		if args.dump:
			self.wf.write()
		else:
			self.wf.send_to_lpad(
				worker_count=args.workers, lpad_filename=args.launchpad_filename)

	def run(self, args):
		# type: (argparse.Namespace) -> None
		"""ScriptBase: Run the CLI. Construct the Workflow object, build its
		tasks, then upload & run it or dump it to a yaml file."""
		wf_basename = getattr(self, 'WORKFLOW_BASENAME', type(self).__name__)

		owner_id = getattr(args, 'owner_id', os.environ['USER'])
		timestamp = getattr(args, 'timestamp', fp.timestamp())
		description = getattr(args, 'description', '')
		name = '{}_{}_{}'.format(owner_id, wf_basename, timestamp)

		storage_basename = getattr(args, 'storage_basename', wf_basename)
		subdir = Workflow.timestamped_description(timestamp, description)
		self.storage_prefix = pp.join(
			Workflow.storage_root(args.storage_root), storage_basename, subdir, '')

		self.wf = Workflow(
			name,
			owner_id=owner_id,
			verbose_logging=args.verbose,
			description=description)

		self.wf.log_info(
			f'\nWorkflow: {name}\n'
			f'Storage prefix: {self.storage_prefix}\n'
			f'Docker internal path prefix: {self.internal_prefix}\n')

		self.build(args)
		self.dumpOrRun(args)
