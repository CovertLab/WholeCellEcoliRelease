"""
Create and upload a simple workflow for simple runs of the workflow software.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import posixpath
from typing import Iterable

import wholecell.utils.filepath as fp
from wholecell.utils import scriptBase
from runscripts.cloud.util.workflow import STORAGE_ROOT_ENV_VAR, Task, Workflow


class WorkflowCLI(scriptBase.ScriptBase):
	"""Abstract base class for a Command Line Interface to build a workflow."""

	# Subclasses can override these:
	DOCKER_IMAGE = 'python:2.7.16'
	DEFAULT_TIMEOUT = Task.DEFAULT_TIMEOUT  # in seconds

	def __init__(self):
		self.storage_prefix = None
		self.wf = None

	def add_task(self, name='', inputs=(), outputs=(), command=(), timeout=0):
		# type: (str, Iterable[str], Iterable[str], Iterable[str], int) -> Task
		task = Task(
			name=name,
			image=self.DOCKER_IMAGE,
			inputs=inputs,
			outputs=outputs,
			storage_prefix=self.storage_prefix,
			internal_prefix='/tmp',
			command=command,
			timeout=timeout if timeout > 0 else self.DEFAULT_TIMEOUT)
		return self.wf.add_task(task)

	def define_parameters(self, parser):
		super(WorkflowCLI, self).define_parameters(parser)
		self.define_option(parser, 'storage_root', str,
			help='The cloud storage root for the output files, usually a GCS'
				 ' bucket name like "sisyphus-crick". Default = ${}'
				 ' environment variable.'.format(STORAGE_ROOT_ENV_VAR))
		parser.add_argument('-w', '--workers', type=int, default=1,
			help='number of worker nodes to launch; default = 1')
		parser.add_argument('--dump', action='store_true',
			help='Dump the built workflow to JSON files for your review *instead* of'
				 ' sending them to the Gaia workflow server. This is useful for'
				 ' testing and debugging. You can upload them manually or re-run'
				 ' this program without `--dump`.')

	def build(self, args):
		# type: (argparse.Namespace) -> None
		"""Build the workflow."""
		raise NotImplementedError("WorkflowCLI subclass must implement build()")

	def dumpOrRun(self, args):
		# type: (argparse.Namespace) -> None
		"""Dump or run the workflow."""
		if args.dump:
			self.wf.write()
		else:
			self.wf.send(worker_count=args.workers)

	def run(self, args):
		# type: (argparse.Namespace) -> None
		basename = type(self).__name__
		owner_id = os.environ['USER']
		timestamp = fp.timestamp()
		name = '{}_{}_{}'.format(owner_id, basename, timestamp)
		self.storage_prefix = posixpath.join(
			Workflow.storage_root(args.storage_root), basename, timestamp, '')
		self.wf = Workflow(name)

		self.wf.log_info('\nStorage prefix: {}'.format(self.storage_prefix))

		self.build(args)
		self.dumpOrRun(args)
