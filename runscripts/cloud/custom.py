#!/usr/bin/env python

"""Run a custom analysis command line in a new Fireworks workflow on Google
Compute Engine (GCE), with access to the Google Cloud Storage (GCS) files from a
WCM workflow run.

Purpose: Use this as as starting point for exploratory analyses and other custom
code that you want to iterate on or use temporarily rather than adding it to the
normal wcm.py workflow. Tweak it to fetch the relevant simOut/ files, e.g.
specific tables.

Alternative: An easier approach is to use gcsfuse to mount the storage bucket on
your development computer then run the analysis there. gcsfuse will fetch files
on demand. This might run slower since gcsfuse takes a bunch of network requests
to enumerate directories, fetch files, and check if cached files are up to date.
Also, fetching files incurs cloud data egress costs at $0.12 per GB.

Alternative: Run lots of `gsutil cp -m ...` commands to download the needed
simOut/ files to local directories, then process them locally. This would save
download time once you write the download code.

Prerequisite: Run `cloud/build-wcm.sh` to build a Docker Image containing the
custom code you want the GCE workers to run. That code needn't be checked in to
git. It just needs to be in the wcEcoli/ directory (and not listed in
`.gcloudignore` or `.dockerignore`).

For the custom code to read input files from a WCM workflow run, the owner_id,
timestamp arg, and description CLI args must match those of the previous run,
and `add_analysis_task()` must list the needed input directories and files.

This CustomWorkflow builder runs on your development computer so you can change
it without rebuilding the Docker Image. When you run it again, if your
fireworkers are still running in GCE, pass `-w0` so it won't try to recreate
them with clashing names. Or don't sweat it -- this uploads a workflow and then
tries to create fireworkers; if the second part fails, your running fireworkers
will still run the newly uploaded Firetasks.

You can create additional fireworker VMs in GCE like this:

	gce --base 20 --count 3 fireworker-$USER

to create workers named "fireworker-$USER-20", "fireworker-$USER-21", and
"fireworker-$USER-22". Here, `--base 20` sets the base index to 20 so the new
names won't collide with running fireworkers "fireworker-$USER-0" through
"fireworker-$USER-19".

You can iterate on the custom code by rebuilding the Docker Image then
re-running this CustomWorkflow builder or (if the same Firetask definition is
fine), just use `lpad rerun_fws -i <N>` to rerun a Firetask. Use the `gce`
command to launch a fireworker VM if they all timed out.
"""

import os
import posixpath as pp

from borealis.util import gcp

from runscripts.cloud.util.workflow import Task
from runscripts.cloud.util.workflow_cli import WorkflowCLI


# SEEDS = sorted(set(range(100)) - {9, 10, 28, 30, 49, 57, 58, 74, 80, 94} | {128} | set(range(300, 310)))


class CustomWorkflow(WorkflowCLI):
	"""Build a workflow that runs a custom command line in Google Cloud."""

	WORKFLOW_BASENAME = 'ComplexCounts'
	DEFAULT_TIMEOUT = 10 * 60  # add_task() default timeout, in seconds

	def __init__(self):
		super().__init__(internal_prefix=pp.join(pp.sep, 'wcEcoli', 'out', 'wf'))
		self.DOCKER_IMAGE = ''  # set it after parsing CLI parameters

	def add_analysis_task(self, seed, num_gens):
		# type: (int, int) -> Task
		"""Add a workflow task that analyzes some of the WCM output data."""
		def in_sim_dir(*path_elements):
			return pp.join(sim_dir, *path_elements)

		base = self.internal('wildtype_000000')
		seed_key = format(seed, '06')  # in 6-digit format
		inputs = []

		for generation in range(num_gens):
			generation_key = 'generation_{:06d}'.format(generation)
			sim_dir = pp.join(
				base, seed_key, generation_key, '000000', 'simOut')
			inputs.append(in_sim_dir('BulkMolecules', 'counts'))
			inputs.append(in_sim_dir('BulkMolecules', 'attributes.json'))
			inputs.append(in_sim_dir('Main', 'time'))
			inputs.append(in_sim_dir('Main', 'attributes.json'))

		return self.add_task(
			name='pull_complex_counts_cloud{}'.format(seed),  # unique task name
			command=['python',
					 '-u',
					 # 'prototypes/subgenerational_analysis/pull_complex_counts_cloud.py',
					 'runscripts/cloud/util/multigen_analysis_example.py',
					 pp.join(base, seed_key)],
			inputs=inputs,
			outputs=[pp.join(base, 'count_out', seed_key, 'complex', '')])

	def build(self, args):
		"""Build the workflow."""
		seeds = [int(seed) for seed in args.seed_list.split()]
		for seed in seeds:
			self.add_analysis_task(seed, args.generations)

	def define_parameters(self, parser):
		self.define_option(parser, 'generations', int, 1, flag='g',
			help='The number of generations to analyze.')
		self.define_option(parser, 'seed_list', str, '0', flag='s',
			help='The list of cell sim seed numbers to analyze, e.g. "0 1 2".')

		self.define_wf_name_parameters(parser)

		super().define_parameters(parser)

	def run(self, args):
		owner_id = args.id or os.environ.get('WF_ID', os.environ['USER'])
		setattr(args, 'owner_id', owner_id)

		# Fetch from and store to a WCM workflow's storage directory, assuming
		# the owner_id, timestamp arg, and description arg match.
		setattr(args, 'storage_basename', 'WCM')

		self.DOCKER_IMAGE = f'gcr.io/{gcp.project()}/{owner_id}-wcm-code'
		super().run(args)


if __name__ == '__main__':
	CustomWorkflow().cli()
