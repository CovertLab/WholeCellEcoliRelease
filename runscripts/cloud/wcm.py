#!/usr/bin/env python

"""
Build a workflow for the Whole Cell Model then send it to the workflow server.
"""

import json
import os
import posixpath as pp
import sys
from typing import Any, Dict, Iterable, Optional, Type

from borealis.util import gcp
from fireworks import FiretaskBase

from wholecell.fireworks.firetasks import (
	ParcaTask,
	VariantSimDataTask,
	SimulationTask,
	SimulationDaughterTask,
	AnalysisParcaTask,
	AnalysisVariantTask,
	AnalysisCohortTask,
	AnalysisSingleTask,
	AnalysisMultiGenTask,
	BuildCausalityNetworkTask,
	WriteJsonTask)
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp
from runscripts.manual.analysisBase import AnalysisBase
from runscripts.cloud.util.workflow import (DEFAULT_LPAD_YAML,
	STORAGE_ROOT_ENV_VAR, Task, Workflow)


# ':latest' -- "You keep using that word. I do not think it means what you think it means."
DOCKER_IMAGE = 'gcr.io/{}/{}-wcm-code'

GCE_VM_MACHINE_TYPE = 'custom-1-5120'  # N1 VM with 1 CPU and 5 GB RAM


class WcmWorkflow(Workflow):
	"""A Workflow builder for the Whole Cell Model."""

	def __init__(self, owner_id, timestamp, verbose_logging=True,
			description='', cli_storage_root=None):
		# type: (str, str, bool, str, Optional[str]) -> None
		name = '{}_WCM_{}'.format(owner_id, timestamp)

		super().__init__(
			name,
			owner_id=owner_id,
			verbose_logging=verbose_logging,
			description=description)

		self.timestamp = timestamp
		self.image = DOCKER_IMAGE.format(gcp.project(), owner_id)

		subdir = Workflow.timestamped_description(timestamp, description)
		self.storage_prefix = pp.join(
			self.storage_root(cli_storage_root), 'WCM', subdir, '')
		self.internal_prefix = pp.join(pp.sep, 'wcEcoli', 'out', 'wf', '')

		self.log_info(
			f'\nWorkflow: {name}\n'
			f'Storage prefix: {self.storage_prefix}\n'
			f'Docker internal path prefix: {self.internal_prefix}\n')

	def internal(self, *path_elements):
		# type: (*str) -> str
		"""Construct a docker container internal file path."""
		return pp.join(self.internal_prefix, *path_elements)

	def remote(self, *path_elements):
		# type: (*str) -> str
		"""Construct a remote GCS storage path within the bucket."""
		return pp.join(self.storage_prefix, *path_elements)

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

	# noinspection PyUnusedLocal
	def build(self, args):
		# type: (Dict[str, Any]) -> None

		# Joining with '' gets a path that ends with the path separator, which
		# tells DockerTask to fetch or store an entire directory tree.
		kb_dir = self.internal(ParcaTask.OUTPUT_SUBDIR, '')
		sim_data_file = pp.join(kb_dir, constants.SERIALIZED_SIM_DATA_FILENAME)
		validation_data_file = pp.join(kb_dir, constants.SERIALIZED_VALIDATION_DATA)

		variant_arg = args['variant']
		variant_spec = (variant_arg[0], int(variant_arg[1]), int(variant_arg[2]))
		variant_type = variant_spec[0]
		variant_count = variant_spec[2] + 1 - variant_spec[1]

		run_analysis = args['run_analysis'] and args['generations'] > 0
		cell_series_out_dir = ''

		if args['workers'] is None:
			args['workers'] = variant_count * args['init_sims']

		# Collect metadata. Analysis Firetasks will expand the _keyed $VARs
		# from (Docker Image) environment variables to update the regular dict
		# entries. Outside the Docker Image, the initial dict entries are good.
		metadata_file = self.internal('metadata', constants.JSON_METADATA_FILE)
		git_hash = fp.run_cmdline("git rev-parse HEAD")
		git_branch = fp.run_cmdline("git symbolic-ref --short HEAD")
		metadata = data.select_keys(
			args,
			scriptBase.METADATA_KEYS,
			git_hash=git_hash,
			_git_hash="$IMAGE_GIT_HASH",
			workflow_git_hash=git_hash,
			git_branch=git_branch,
			_git_branch="$IMAGE_GIT_BRANCH",
			workflow_git_branch=git_branch,
			description=args['description'] or 'WCM',
			time=self.timestamp,
			_time="$IMAGE_TIMESTAMP",
			python=sys.version.splitlines()[0],
			variant=variant_type,
			total_variants=str(variant_count),
			total_gens=args['generations'])

		python_args = dict(output_file=metadata_file, data=metadata)  # type: Dict[str, Any]
		metadata_task = self.add_python_task(WriteJsonTask, python_args,
			name='write_metadata',
			inputs=[],
			outputs=[metadata_file],
			timeout=90)

		if args['run_parca']:
			python_args = data.select_keys(
				args,
				scriptBase.PARCA_KEYS,
				debug=args['debug_parca'],
				output_directory=kb_dir)
			parca_task = self.add_python_task(ParcaTask, python_args,
				name='parca',
				outputs=[kb_dir])
		else:
			print('    (Skipping the Parca step per the --no-run-parca option.)')

		if run_analysis:
			parca_plot_dir = self.internal(constants.KB_PLOT_OUTPUT_DIR, '')
			python_args = data.select_keys(
				args, scriptBase.ANALYSIS_KEYS,
				input_directory=kb_dir,
				input_sim_data=sim_data_file,
				input_validation_data=validation_data_file,
				output_plots_directory=parca_plot_dir,
				metadata=metadata)
			analysis_parca_task = self.add_python_task(AnalysisParcaTask,
				python_args,
				name='analysis_parca',
				inputs=[kb_dir],
				outputs=[parca_plot_dir])

		variant_analysis_inputs = [kb_dir]

		sim_args = data.select_keys(args, scriptBase.SIM_KEYS)

		for i, subdir in fp.iter_variants(*variant_spec):
			variant_sim_data_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_KB, '')
			variant_metadata_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_METADATA, '')
			variant_sim_data_modified_file = pp.join(
				variant_sim_data_dir, constants.SERIALIZED_SIM_DATA_MODIFIED)
			md_cohort = dict(metadata, variant_function=variant_type,
				variant_index=i)
			variant_analysis_inputs.append(variant_metadata_dir)

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

			for j in range(arg_seed, arg_seed + args['init_sims']):  # init sim seeds
				seed_dir = self.internal(subdir, '{:06d}'.format(j))
				md_multigen = dict(md_cohort, seed=j)

				this_variant_this_seed_multigen_analysis_inputs = [kb_dir, variant_sim_data_dir]

				for k in range(args['generations']):
					gen_dir = pp.join(seed_dir, "generation_{:06d}".format(k))
					md_single = dict(md_multigen, gen=k)

					# l is the daughter number among all of this generation's cells.
					# l in [0] for single daughters; l in range(2**k) for dual daughters.
					for l in [0]:
						cell_dir = pp.join(gen_dir, '{:06d}'.format(l))
						cell_sim_out_dir = pp.join(cell_dir, 'simOut', '')

						python_args = dict(sim_args,
							input_sim_data=variant_sim_data_modified_file,
							output_directory=cell_sim_out_dir)
						inputs = [kb_dir, variant_sim_data_dir]

						if k == 0:
							python_args['seed'] = j
							firetask = SimulationTask
						else:
							firetask = SimulationDaughterTask
							parent_gen_dir = pp.join(
								seed_dir, 'generation_{:06d}'.format(k - 1))
							parent_cell_dir = pp.join(parent_gen_dir, '{:06d}'.format(l // 2))
							parent_cell_sim_out_dir = pp.join(parent_cell_dir, 'simOut', '')
							daughter_state_path = pp.join(
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
							plot_dir = pp.join(cell_dir, AnalysisBase.OUTPUT_SUBDIR, '')
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

						if args['build_causality_network'] and not cell_series_out_dir:
							cell_series_out_dir = pp.join(cell_dir, 'seriesOut', '')
							python_args = dict(
								input_results_directory=cell_sim_out_dir,
								input_sim_data=variant_sim_data_modified_file,
								output_dynamics_directory=cell_series_out_dir,
								metadata=md_single)
							causality_task = self.add_python_task(BuildCausalityNetworkTask,
								python_args,
								name='causality_' + cell_id,
								inputs=[cell_sim_out_dir, variant_sim_data_dir],
								outputs=[cell_series_out_dir])

				if run_analysis:
					multigen_plot_dir = pp.join(seed_dir, AnalysisBase.OUTPUT_SUBDIR, '')
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
		return '''The command line option names are long but you can use any
			unambiguous prefix.

			Set the environment variable WF_ID if you want to select the
			Docker Container Image named "${WF_ID}-wcm-code" for the workers.
			$WF_ID defaults to $USER. It also determines the output storage
			path prefix. Run `cloud/build-wcm.sh $WF_ID` to build the WCM
			Container Image with that name.'''

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
			help='Dump the built workflow to a YAML file for review *instead*'
				 ' of sending it to the launchpad DB server. This is useful'
				 ' for testing and debugging.')
		parser.add_argument('-l', dest='launchpad_filename',
			default=DEFAULT_LPAD_YAML,
			help='Launchpad config YAML filename (default="{}").'.format(
				DEFAULT_LPAD_YAML))
		parser.add_argument('-w', '--workers', type=int,
			help='The number of worker nodes to launch, with a smart default.')
		self.define_option(parser, 'machine_type', str,
			default=GCE_VM_MACHINE_TYPE,
			help='GCE (Compute Engine) worker VM machine-type to launch.'
				 ' E.g. "n1-standard-1" for 1st generation machine with 1 CPU'
				 ' and 3.75 GB RAM, "n2-standard-2" for 2nd gen machine with'
				 ' 2 CPUs and 8 GB RAM, "custom-1-5120" for a custom 1st gen'
				 ' machine with 1 CPU and 5 GB RAM, or "n2-custom-2-5120" for'
				 ' 2nd gen machine with 2 CPUs and 5 GB RAM.'
				 ' See https://cloud.google.com/compute/docs/machine-types,'
				 ' https://cloud.google.com/compute/vm-instance-pricing, and'
				 ' https://cloud.google.com/compute/docs/instances/'
				 'creating-instance-with-custom-machine-type ')

		# Parca
		self.define_parca_options(parser, run_parca_option=True)

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
			help="Build the Causality network files for one sim generation.")

		super().define_parameters(parser)

	def parse_args(self):
		args = super().parse_args()
		args.cpus = max(args.cpus, 1)

		assert args.generations >= 0
		assert args.init_sims > 0
		assert args.length_sec > 0
		assert args.cpus > 0
		return args

	def run(self, args):
		wf = wc_ecoli_workflow(vars(args))

		if args.dump:
			wf.write()
		else:
			gce_options = {'machine-type': args.machine_type}
			wf.send_to_lpad(
				worker_count=args.workers,
				lpad_filename=args.launchpad_filename,
				gce_options=gce_options)


if __name__ == '__main__':
	script = RunWcm()
	script.cli()
