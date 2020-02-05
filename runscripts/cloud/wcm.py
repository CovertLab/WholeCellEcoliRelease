#!/usr/bin/env python

"""
Build a workflow for the Whole Cell Model then send it to the Gaia server in
Google Cloud.
"""

from __future__ import absolute_import, division, print_function

import json
import os
import posixpath
from pprint import pprint
import re
from typing import Any, Dict, Iterable, Optional, Type

from fireworks import FiretaskBase

from wholecell.fireworks.firetasks import (
	ParcaTask,
	VariantSimDataTask,
	SimulationTask,
	SimulationDaughterTask,
	AnalysisVariantTask,
	AnalysisCohortTask,
	AnalysisSingleTask,
	AnalysisMultiGenTask,
	BuildCausalityNetworkTask,
	WriteJsonTask)
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp
from runscripts.manual.analysisBase import AnalysisBase
from runscripts.cloud.util.workflow import STORAGE_ROOT_ENV_VAR, Task, Workflow


# ':latest' -- "You keep using that word. I do not think it means what you think it means."
DOCKER_IMAGE = 'gcr.io/allen-discovery-center-mcovert/{}-wcm-code'
USE_GAIA = False


class WcmWorkflow(Workflow):
	"""A Workflow builder for the Whole Cell Model."""

	def __init__(self, owner_id, timestamp, verbose_logging=True,
			description='', cli_storage_root=None):
		# type: (str, str, bool, str, Optional[str]) -> None
		name = '{}_WCM_{}'.format(owner_id, timestamp)
		super(WcmWorkflow, self).__init__(
			name, owner_id=owner_id, verbose_logging=verbose_logging)

		self.timestamp = timestamp
		self.image = DOCKER_IMAGE.format(owner_id)

		subdir = self.timestamp + (
			'__' + _sanitize_description(description) if description else '')
		self.storage_prefix = posixpath.join(
			self.storage_root(cli_storage_root), 'WCM', subdir, '')
		self.internal_prefix = posixpath.join(posixpath.sep, 'wcEcoli', 'out', 'wf')

		self.log_info('\nStorage prefix: {}'.format(self.storage_prefix))

		if description:
			self.add_properties(description=description)

	def internal(self, *path_elements):
		# type: (*str) -> str
		"""Construct a file path that's internal to the task's container."""
		return posixpath.join(self.internal_prefix, *path_elements)

	def remote(self, *path_elements):
		# type: (*str) -> str
		"""Construct a remote GCS storage path within the bucket."""
		return posixpath.join(self.storage_prefix, *path_elements)

	def add_python_task(self, firetask, python_args, name='', inputs=(),
			outputs=(), timeout=0):
		# type: (Type[FiretaskBase], Dict[str, Any], str, Iterable[str], Iterable[str], int) -> Task
		"""Add a Python task to the workflow and return it. Store its
		stdout + stderr as storage_prefix/logs/name.log .
		Turn on Python '-u' so it doesn't buffer output for long time.
		"""
		# TODO(jerry): An option to emit a task that runs `firetask` directly
		#  rather than via a command line in a Docker image. That requires
		#  running the firetask with access to the wcEcoli code and the file
		#  system. You wouldn't do that on GCE unless there's an NFS mount for
		#  the wcEcoli code and data.
		return self.add_task(Task(
			name=name,
			image=self.image,
			command=['python', '-u', '-m', 'wholecell.fireworks.runTask',
					firetask.__name__, json.dumps(python_args)],
			inputs=inputs,
			outputs=outputs,
			storage_prefix=self.storage_prefix,
			internal_prefix=self.internal_prefix,
			timeout=timeout))

	def build(self, args):
		# type: (Dict[str, Any]) -> None

		# Joining with '' gets a path that ends with the path separator, which
		# tells Sisyphus to pull or push an entire directory tree.
		kb_dir = self.internal(ParcaTask.OUTPUT_SUBDIR, '')
		sim_data_file = posixpath.join(kb_dir, constants.SERIALIZED_SIM_DATA_FILENAME)
		validation_data_file = posixpath.join(kb_dir, constants.SERIALIZED_VALIDATION_DATA)

		variant_arg = args['variant']
		variant_spec = (variant_arg[0], int(variant_arg[1]), int(variant_arg[2]))
		variant_type = variant_spec[0]
		variant_count = variant_spec[2] + 1 - variant_spec[1]

		run_analysis = args['run_analysis'] and args['generations'] > 0

		if args['workers'] is None:
			args['workers'] = variant_count * args['init_sims']

		metadata_file = self.internal('metadata', constants.JSON_METADATA_FILE)
		metadata = data.select_keys(
			args,
			scriptBase.METADATA_KEYS,
			git_hash="$IMAGE_GIT_HASH",  # expanded by the Docker Container shell!
			git_branch="$IMAGE_GIT_BRANCH",
			workflow_git_hash=fp.run_cmdline("git rev-parse HEAD"),
			workflow_git_branch=fp.run_cmdline("git symbolic-ref --short HEAD"),
			description=args['description'] or 'WCM',
			time="$IMAGE_TIMESTAMP",
			workflow_time=self.timestamp,
			variant=variant_type,
			total_variants=str(variant_count),
			total_gens=args['generations'])

		python_args = dict(output_file=metadata_file, data=metadata)
		metadata_task = self.add_python_task(WriteJsonTask, python_args,
			name='write_metadata',
			inputs=[],
			outputs=[metadata_file],
			timeout=90)

		python_args = data.select_keys(
			args,
			scriptBase.PARCA_KEYS,
			debug=args['debug_parca'],
			output_directory=kb_dir)
		parca_task = self.add_python_task(ParcaTask, python_args,
			name='parca',
			outputs=[kb_dir])

		variant_analysis_inputs = [kb_dir]

		sim_args = data.select_keys(args, scriptBase.SIM_KEYS)

		for i, subdir in fp.iter_variants(*variant_spec):
			variant_sim_data_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_KB, '')
			variant_metadata_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_METADATA, '')
			variant_sim_data_modified_file = posixpath.join(
				variant_sim_data_dir, constants.SERIALIZED_SIM_DATA_MODIFIED)
			md_cohort = dict(metadata, variant_function=variant_type,
				variant_index=i)

			python_args = dict(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=sim_data_file,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_dir)
			variant_task = self.add_python_task(VariantSimDataTask, python_args,
				name='variant_{}_{}'.format(variant_type, i),
				inputs=[kb_dir],
				outputs=[variant_sim_data_dir, variant_metadata_dir],
				timeout=90)

			this_variant_cohort_analysis_inputs = [kb_dir, variant_sim_data_dir]
			variant_analysis_inputs.append(variant_sim_data_dir)
			arg_seed = args['seed']

			for j in xrange(arg_seed, arg_seed + args['init_sims']):  # init sim seeds
				seed_dir = self.internal(subdir, '{:06d}'.format(j))
				md_multigen = dict(md_cohort, seed=j)

				this_variant_this_seed_multigen_analysis_inputs = [kb_dir, variant_sim_data_dir]

				for k in xrange(args['generations']):
					gen_dir = posixpath.join(seed_dir, "generation_{:06d}".format(k))
					md_single = dict(md_multigen, gen=k)

					# l is the daughter number among all of this generation's cells.
					# l in [0] for single daughters; l in range(2**k) for dual daughters.
					for l in [0]:
						cell_dir = posixpath.join(gen_dir, '{:06d}'.format(l))
						cell_sim_out_dir = posixpath.join(cell_dir, 'simOut', '')

						python_args = dict(sim_args,
							input_sim_data=variant_sim_data_modified_file,
							output_directory=cell_sim_out_dir)
						inputs = [kb_dir, variant_sim_data_dir]

						if k == 0:
							python_args['seed'] = j
							firetask = SimulationTask
						else:
							firetask = SimulationDaughterTask
							parent_gen_dir = posixpath.join(
								seed_dir, 'generation_{:06d}'.format(k - 1))
							parent_cell_dir = posixpath.join(parent_gen_dir, '{:06d}'.format(l // 2))
							parent_cell_sim_out_dir = posixpath.join(parent_cell_dir, 'simOut', '')
							daughter_state_path = posixpath.join(
								parent_cell_sim_out_dir,
								constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
							python_args['inherited_state_path'] = daughter_state_path
							python_args['seed'] = (j + 1) * ((2 ** k - 1) + l)
							inputs += [parent_cell_sim_out_dir]

						cell_id = 'Var{:03d}_Seed{:03d}_Gen{:03d}_Cell{:03d}'.format(i, j, k, l)
						sim_task = self.add_python_task(firetask, python_args,
							name='simulation_' + cell_id,
							inputs=inputs,
							outputs=[cell_sim_out_dir])

						variant_analysis_inputs.append(cell_sim_out_dir)
						this_variant_cohort_analysis_inputs.append(cell_sim_out_dir)
						this_variant_this_seed_multigen_analysis_inputs.append(cell_sim_out_dir)

						if run_analysis:
							plot_dir = posixpath.join(cell_dir, AnalysisBase.OUTPUT_SUBDIR, '')
							python_args = data.select_keys(
								args, scriptBase.ANALYSIS_KEYS,
								input_results_directory=cell_sim_out_dir,
								input_sim_data=variant_sim_data_modified_file,
								input_validation_data=validation_data_file,
								output_plots_directory=plot_dir,
								metadata=md_single)
							analysis_single_task = self.add_python_task(AnalysisSingleTask,
								python_args,
								name='analysis_' + cell_id,
								inputs=[kb_dir, variant_sim_data_dir, cell_sim_out_dir],
								outputs=[plot_dir])

						if args['build_causality_network']:
							cell_series_out_dir = posixpath.join(cell_dir, 'seriesOut', '')
							# NOTE: This could reuse the Causality network over the variant. For
							# Sisyphus it'd take moving that work from BuildCausalityNetworkTask
							# to VariantSimDataTask, but it wouldn't save much space and time.
							python_args = dict(
								input_results_directory=cell_sim_out_dir,
								input_sim_data=variant_sim_data_modified_file,
								output_network_directory=cell_series_out_dir,
								output_dynamics_directory=cell_series_out_dir,
								metadata=md_single)
							causality_task = self.add_python_task(BuildCausalityNetworkTask,
								python_args,
								name='causality_' + cell_id,
								inputs=[cell_sim_out_dir, variant_sim_data_dir],
								outputs=[cell_series_out_dir])

				if run_analysis:
					multigen_plot_dir = posixpath.join(seed_dir, AnalysisBase.OUTPUT_SUBDIR, '')
					python_args = data.select_keys(
						args, scriptBase.ANALYSIS_KEYS,
						input_seed_directory=seed_dir,
						input_sim_data=variant_sim_data_modified_file,
						input_validation_data=validation_data_file,
						output_plots_directory=multigen_plot_dir,
						metadata=md_multigen)
					analysis_multigen_task = self.add_python_task(AnalysisMultiGenTask,
						python_args,
						name='analysis_multigen_Var{}_Seed{}'.format(i, j),
						inputs=this_variant_this_seed_multigen_analysis_inputs,
						outputs=[multigen_plot_dir])

			if run_analysis:
				cohort_plot_dir = self.internal(subdir, AnalysisBase.OUTPUT_SUBDIR, '')
				python_args = data.select_keys(
					args, scriptBase.ANALYSIS_KEYS,
					input_variant_directory=self.internal(subdir),
					input_sim_data=variant_sim_data_modified_file,
					input_validation_data=validation_data_file,
					output_plots_directory=cohort_plot_dir,
					metadata=md_cohort)
				analysis_cohort_task = self.add_python_task(AnalysisCohortTask,
					python_args,
					name='analysis_cohort_Var{}'.format(i),
					inputs=this_variant_cohort_analysis_inputs,
					outputs=[cohort_plot_dir])

		if run_analysis:
			variant_plot_dir = self.internal(AnalysisBase.OUTPUT_SUBDIR, '')
			python_args = data.select_keys(
				args, scriptBase.ANALYSIS_KEYS,
				input_directory=self.internal(''),
				input_sim_data=sim_data_file,
				input_validation_data=validation_data_file,
				output_plots_directory=variant_plot_dir,
				metadata=metadata)
			analysis_variant_task = self.add_python_task(AnalysisVariantTask,
				python_args,
				name='analysis_variant',
				inputs=variant_analysis_inputs,
				outputs=[variant_plot_dir])

def _sanitize_description(description):
	# type (str) -> str
	"""Sanitize the description and check that it's legal in a file path."""
	description = description.replace(' ', '_')

	pattern = r'[-.\w]*$'
	assert re.match(pattern, description), (
		"description {!r} doesn't match the regex pattern {!r} for a file path."
			.format(description, pattern))

	return description


def wc_ecoli_workflow(args):
	# type: (Dict[str, Any]) -> WcmWorkflow
	"""Build a workflow for wcEcoli."""
	owner_id = args['id'] or os.environ.get('WF_ID', os.environ['USER'])
	timestamp = args['timestamp']
	description = args['description']

	wf = WcmWorkflow(owner_id, timestamp, verbose_logging=args['verbose'],
		description=description, cli_storage_root=args['storage_root'])
	wf.build(args)
	return wf


class RunWcm(scriptBase.ScriptBase):
	"""Command line interpreter to run a WCM workflow."""

	def description(self):
		return 'E. coli Whole Cell Model workflow'

	def epilog(self):
		return '''Set the environment variable $WF_ID if you want to select the
			Docker container image named "$WF_ID-wcm-code" for the workers.
			$WF_ID defaults to $USER. This also determines the output storage
			path prefix. Run `cloud/build-wcm.sh $WF_ID` to build the WCM
			container image with that name.

			(The command line option names are long but you can use any
			unambiguous prefix.)'''

	def define_parameters(self, parser):
		self.define_option(parser, 'description', str, '',
			help='A simulation description to append to the output folder name.')
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
		self.define_option(parser, 'storage_root', str,
			help='The cloud storage root for the output files, usually a GCS'
				 ' bucket name like "sisyphus-crick". Default = ${}'
				 ' environment variable.'.format(STORAGE_ROOT_ENV_VAR))
		self.define_parameter_bool(parser, 'verbose', True,
			help='Verbose workflow builder logging')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use in the Parca and analysis'
				 ' steps. Default = 1.')
		self.define_parameter_bool(parser, 'dump', False,
			help='Dump the built workflow to JSON files for'
				 ' review *instead* of sending them to the Gaia workflow'
				 ' server. This is useful for testing and debugging. You can'
				 ' upload them manually or re-run this program without `--dump`.')
		parser.add_argument('-w', '--workers', type=int,
			help='The number of worker nodes to launch, with a smart default.')

		# Parca
		self.define_parca_options(parser)

		# Simulation
		self.define_sim_loop_options(parser)
		self.define_sim_options(parser)

		# For Parca and Sim. define_parca_options() and define_sim_options()
		# can't both call this since ArgumentParser would raise an error.
		self.define_elongation_options(parser)

		# Analyses
		self.define_parameter_bool(parser, 'run_analysis', True,
			help="Run the Variant, Cohort, Multigen, and Single analyses.")
		parser.add_argument('-p', '--plot', nargs='+', default=[],
			help='''Names the analysis plots to run, e.g. plot filenames
				like "aaCounts.py" or "aaCounts" and tags like "METABOLISM"
				as defined in the __init__.py files. If omitted, the default is
				"CORE", which names the plots recommended for everyday
				development. Use "ACTIVE" to run all active plots in this
				category. You can name specific analysis files but any
				analysis categories that don't have those files will print
				error messages.''')
		self.define_parameter_bool(parser, 'compile', False,
			'Compiles output images into one file (only for .png).')
		self.define_parameter_bool(parser, 'build_causality_network', False,
			help="Build the Causality network files for each sim generation.")

		super(RunWcm, self).define_parameters(parser)

	def parse_args(self):
		args = super(RunWcm, self).parse_args()
		args.cpus = max(args.cpus, 1)

		assert args.generations >= 0
		assert args.init_sims > 0
		assert args.length_sec > 0
		assert args.cpus > 0
		return args

	def run(self, args):
		wf = wc_ecoli_workflow(vars(args))

		if USE_GAIA:
			if args.dump:
				wf.write_for_gaia()
			else:
				wf.send_to_gaia(worker_count=args.workers)
			return

		if args.dump:
			# TODO(jerry): Write a yaml spec file.
			fw_wf = wf.build_workflow()
			pprint(fw_wf)
		else:
			wf.send_to_lpad(worker_count=args.workers)


if __name__ == '__main__':
	script = RunWcm()
	script.cli()
