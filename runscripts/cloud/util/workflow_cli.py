"""
Create and upload a simple workflow for simple runs of the workflow software.
"""

from __future__ import absolute_import, division, print_function

import argparse
import os
import posixpath
from pprint import pprint
from typing import Iterable, Optional

import wholecell.utils.filepath as fp
from wholecell.utils import scriptBase
from runscripts.cloud.util.workflow import (DEFAULT_LPAD_YAML,
	STORAGE_ROOT_ENV_VAR, Task, Workflow)


USE_GAIA = False

class WorkflowCLI(scriptBase.ScriptBase):
	"""Abstract base class for a Command Line Interface to build a workflow."""

	# Subclasses can override these:
	DOCKER_IMAGE = 'python:2.7.16'
	DEFAULT_TIMEOUT = Task.DEFAULT_TIMEOUT  # in seconds

	def __init__(self):
		super(WorkflowCLI, self).__init__()
		self.storage_prefix = ''
		self.wf = None  # type: Optional[Workflow]

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
		parser.add_argument('-l', dest='launchpad_filename',
			default=DEFAULT_LPAD_YAML,
			help='Launchpad config YAML filename (default="{}").'.format(
				DEFAULT_LPAD_YAML))
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
		if USE_GAIA:
			if args.dump:
				self.wf.write_for_gaia()
			else:
				self.wf.send_to_gaia(worker_count=args.workers)
			return

		if args.dump:
			# TODO(jerry): Write a yaml spec file.
			fw_wf = self.wf.build_workflow()
			pprint(fw_wf)
		else:
			self.wf.send_to_lpad(
				worker_count=args.workers, lpad_filename=args.launchpad_filename)

	def run(self, args):
		# type: (argparse.Namespace) -> None
		basename = type(self).__name__
		owner_id = os.environ['USER']
		timestamp = fp.timestamp()
		name = '{}_{}_{}'.format(owner_id, basename, timestamp)
		self.storage_prefix = posixpath.join(
			Workflow.storage_root(args.storage_root), basename, timestamp, '')
		self.wf = Workflow(name, owner_id=owner_id)

		self.wf.log_info('\nStorage prefix: {}'.format(self.storage_prefix))

		self.build(args)
		self.dumpOrRun(args)
