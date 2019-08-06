#!/usr/bin/env python

"""
Build a workflow for the Whole Cell Model then send it to the Gaia server in
Google Cloud.
"""

from __future__ import absolute_import, division, print_function

import json
import os
import re
from typing import Any, Dict, Iterable, Mapping

from wholecell.fireworks.firetasks import ParcaTask, VariantSimDataTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp
from runscripts.manual.analysisBase import AnalysisBase
from runscripts.cloud.util.workflow import Task, Workflow


DOCKER_IMAGE = 'gcr.io/allen-discovery-center-mcovert/{}-wcm-code:latest'
STORAGE_PREFIX_ROOT = 'sisyphus:data/'

DEFAULT_VARIANT = ['wildtype', '0', '0']


def select_keys(mapping, keys, **kwargs):
	# type: (Mapping[str, Any], Iterable[str], **Any) -> Dict[str, Any]
	"""Return a dict of the mapping entries with the given keys plus the kwargs."""
	result = {key: mapping[key] for key in keys if mapping[key] is not None}
	result.update(**kwargs)
	return result


class WcmWorkflow(Workflow):
	"""A Workflow builder for the Whole Cell Model."""

	def __init__(self, owner_id, timestamp, verbose_logging=True, description=''):
		# type: (str, str, bool, str) -> None
		name = 'WCM_{}_{}'.format(owner_id, timestamp)
		super(WcmWorkflow, self).__init__(name, verbose_logging=verbose_logging)

		self.owner_id = owner_id
		self.timestamp = timestamp
		self.image = DOCKER_IMAGE.format(self.owner_id)

		subdir = self.timestamp + ('__' + description if description else '')
		self.storage_prefix = os.path.join(
			STORAGE_PREFIX_ROOT, self.owner_id, subdir, '')
		self.internal_prefix = os.path.join(os.sep, 'wcEcoli', 'out', 'wf')

		self.log_info('\nStorage prefix: {}'.format(self.storage_prefix))

	def internal(self, *path_elements):
		# type: (*str) -> str
		"""Construct a file path that's internal to the task's container."""
		return os.path.join(self.internal_prefix, *path_elements)

	def remote(self, *path_elements):
		# type: (*str) -> str
		"""Construct a remote GCS storage path within the bucket."""
		return os.path.join(self.storage_prefix, *path_elements)

	def add_python_task(self, firetask, python_args, upstream_tasks=(), name='',
			inputs=(), outputs=()):
		# type: (str, Dict[str, Any], Iterable[Task], str, Iterable[str], Iterable[str]) -> Task
		"""Add a Python task to the workflow and return it."""
		return self.add_task(Task(
			upstream_tasks=upstream_tasks,
			name=name,
			image=self.image,
			command=['python', '-u', '-m', 'wholecell.fireworks.runTask',
					firetask, json.dumps(python_args)],
			inputs=inputs,
			outputs=outputs,
			storage_prefix=self.storage_prefix,
			internal_prefix=self.internal_prefix))

	def build(self, args):
		# type: (Dict[str, Any]) -> None

		# Joining with '' gets a path that ends with the path separator, which
		# tells Sisyphus to pull or push an entire directory tree.
		kb_dir = self.internal(ParcaTask.OUTPUT_SUBDIR, '')
		sim_data_file = os.path.join(kb_dir, constants.SERIALIZED_SIM_DATA_FILENAME)
		validation_data_file = os.path.join(kb_dir, constants.SERIALIZED_VALIDATION_DATA)

		variant_arg = args['variant']
		variant_spec = (variant_arg[0], int(variant_arg[1]), int(variant_arg[2]))
		variant_type = variant_spec[0]
		variant_count = variant_spec[2] + 1 - variant_spec[1]

		run_analysis = args['run_analysis'] and args['generations'] > 0

		if args['workers'] is None:
			args['workers'] = variant_count * args['init_sims']

		metadata_file = self.internal('metadata', constants.JSON_METADATA_FILE)
		metadata = select_keys(args,
			('generations', 'mass_distribution', 'growth_rate_noise',
			'd_period_division', 'translation_supply', 'trna_charging'),
			git_hash=fp.run_cmdline("git rev-parse HEAD"),
			git_branch=fp.run_cmdline("git symbolic-ref --short HEAD"),
			description=args['description'] or 'WCM',
			time=self.timestamp,
			variant=variant_type,
			total_variants=str(variant_count),
			total_gens=args['generations'])

		python_args = dict(output_file=metadata_file, data=metadata)
		metadata_task = self.add_python_task('write_json', python_args, (),
			name='write_metadata',
			inputs=[kb_dir],  # TODO(jerry): TEMPORARY workaround to delay this
				# task so its worker doesn't exit while the Parca runs.
			outputs=[metadata_file])

		python_args = select_keys(args,
			('ribosome_fitting', 'rnapoly_fitting', 'cpus'),
			debug=args['debug_parca'],
			output_directory=kb_dir)
		parca_task = self.add_python_task('parca', python_args, (),
			name='parca',
			outputs=[kb_dir])

		variant_analysis_inputs = [kb_dir]

		sim_args = select_keys(args,
			('timeline', 'length_sec', 'timestep_safety_frac', 'timestep_max',
			'timestep_update_freq', 'mass_distribution', 'growth_rate_noise',
			'd_period_division', 'translation_supply', 'trna_charging'))

		for i, subdir in fp.iter_variants(*variant_spec):
			variant_sim_data_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_KB, '')
			variant_metadata_dir = self.internal(subdir,
				VariantSimDataTask.OUTPUT_SUBDIR_METADATA, '')
			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_dir, constants.SERIALIZED_SIM_DATA_MODIFIED)
			md_cohort = dict(metadata, variant_function=variant_type,
				variant_index=i)

			python_args = dict(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=sim_data_file,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_dir)
			variant_task = self.add_python_task('variant_sim_data', python_args,
				(parca_task,),
				name='variant_{}_{}'.format(variant_type, i),
				outputs=[variant_sim_data_dir, variant_metadata_dir])

			this_variant_cohort_analysis_inputs = [kb_dir, variant_sim_data_dir]
			variant_analysis_inputs.append(variant_sim_data_dir)

			for j in xrange(args['init_sims']):  # seed
				seed_dir = self.internal(subdir, '{:06d}'.format(j))
				md_multigen = dict(md_cohort, seed=j)

				this_variant_this_seed_multigen_analysis_inputs = [kb_dir, variant_sim_data_dir]

				for k in xrange(args['generations']):
					gen_dir = os.path.join(seed_dir, "generation_{:06d}".format(k))
					md_single = dict(md_multigen, gen=k)

					# l is the daughter number among all of this generation's cells.
					# l in [0] for single daughters; l in range(2**k) for dual daughters.
					for l in [0]:
						cell_dir = os.path.join(gen_dir, '{:06d}'.format(l))
						cell_sim_out_dir = os.path.join(cell_dir, 'simOut', '')

						python_args = dict(sim_args,
							input_sim_data=variant_sim_data_modified_file,
							output_directory=cell_sim_out_dir,
							seed=j)

						if k == 0:
							firetask = 'simulation'
							inputs = []
						else:
							firetask = 'simulation_daughter'
							parent_gen_dir = os.path.join(
								seed_dir, 'generation_{:06d}'.format(k - 1))
							parent_cell_dir = os.path.join(parent_gen_dir, '{:06d}'.format(l // 2))
							parent_cell_sim_out_dir = os.path.join(parent_cell_dir, 'simOut', '')
							daughter_state_path = os.path.join(
								parent_cell_sim_out_dir,
								constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
							python_args['inherited_state_path'] = daughter_state_path
							inputs=[parent_cell_sim_out_dir]

						cell_id = 'Var{}_Seed{}_Gen{}_Cell{}'.format(i, j, k, l)
						sim_task = self.add_python_task(firetask, python_args,
							(parca_task, variant_task),
							name='simulation_' + cell_id,
							inputs=inputs,
							outputs=[cell_sim_out_dir])

						variant_analysis_inputs.append(cell_sim_out_dir)
						this_variant_cohort_analysis_inputs.append(cell_sim_out_dir)
						this_variant_this_seed_multigen_analysis_inputs.append(cell_sim_out_dir)

						if run_analysis:
							plot_dir = os.path.join(cell_dir, AnalysisBase.OUTPUT_SUBDIR, '')
							python_args = dict(
								input_results_directory=cell_sim_out_dir,
								input_sim_data=variant_sim_data_modified_file,
								input_validation_data=validation_data_file,
								output_plots_directory=plot_dir,
								metadata=md_single,
								plots_to_run=args['plot'],
								cpus=args['cpus'])
							analysis_single_task = self.add_python_task('analysis_single',
								python_args,
								(parca_task, variant_task, sim_task),
								name='analysis_' + cell_id,
								outputs=[plot_dir])

						if args['build_causality_network']:
							cell_series_out_dir = os.path.join(cell_dir, 'seriesOut', '')
							# NOTE: This could reuse the Causality network over the variant. For
							# Sisyphus it'd take moving that work from BuildCausalityNetworkTask
							# to VariantSimDataTask, but it wouldn't save much space and time.
							python_args = dict(
								input_results_directory=cell_sim_out_dir,
								input_sim_data=variant_sim_data_modified_file,
								output_network_directory=cell_series_out_dir,
								output_dynamics_directory=cell_series_out_dir,
								metadata=md_single)
							causality_task = self.add_python_task('build_causality_network',
								python_args, (),
								name='causality_' + cell_id,
								inputs=[cell_sim_out_dir, variant_sim_data_dir],
								outputs=[cell_series_out_dir])

				if run_analysis:
					multigen_plot_dir = os.path.join(seed_dir, AnalysisBase.OUTPUT_SUBDIR, '')
					python_args = dict(
						input_seed_directory=seed_dir,
						input_sim_data=variant_sim_data_modified_file,
						input_validation_data=validation_data_file,
						output_plots_directory=multigen_plot_dir,
						plots_to_run=args['plot'],
						cpus=args['cpus'],
						metadata=md_multigen)
					analysis_multigen_task = self.add_python_task('analysis_multigen',
						python_args, (),
						name='analysis_multigen_Var{}_Seed{}'.format(i, j),
						inputs=this_variant_this_seed_multigen_analysis_inputs,
						outputs=[multigen_plot_dir])

			if run_analysis:
				cohort_plot_dir = self.internal(subdir, AnalysisBase.OUTPUT_SUBDIR, '')
				python_args = dict(
					input_variant_directory=self.internal(subdir),
					input_sim_data=variant_sim_data_modified_file,
					input_validation_data=validation_data_file,
					output_plots_directory=cohort_plot_dir,
					plots_to_run=args['plot'],
					cpus=args['cpus'],
					metadata=md_cohort)
				analysis_cohort_task = self.add_python_task('analysis_cohort',
					python_args, (),
					name='analysis_cohort_Var{}'.format(i),
					inputs=this_variant_cohort_analysis_inputs,
					outputs=[cohort_plot_dir])

		if run_analysis:
			variant_plot_dir = self.internal(AnalysisBase.OUTPUT_SUBDIR, '')
			python_args = dict(
				input_directory=self.internal(''),
				input_sim_data=sim_data_file,
				input_validation_data=validation_data_file,
				output_plots_directory=variant_plot_dir,
				plots_to_run=args['plot'],
				cpus=args['cpus'],
				metadata=metadata)
			analysis_variant_task = self.add_python_task('analysis_variant',
				python_args, (),
				name='analysis_variant',
				inputs=variant_analysis_inputs,
				outputs=[variant_plot_dir])


def wc_ecoli_workflow(args):
	# type: (Dict[str, Any]) -> WcmWorkflow
	"""Build a workflow for wcEcoli."""
	owner_id = os.environ.get('WF_ID', os.environ['USER'])
	timestamp = fp.timestamp()
	description = args['description'].replace(' ', '_')

	pattern = r'[-.\w]*$'
	assert re.match(pattern, description), (
		"description {!r} doesn't match the regex pattern {!r}.".format(
			description, pattern))

	wf = WcmWorkflow(owner_id, timestamp, verbose_logging=args['verbose'],
		description=description)
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

			The command line option names are long but you can use any
			unambiguous prefix.'''

	def define_parameters(self, parser):
		def add_option(name, key, datatype, help):
			"""Add an option with the given name and datatype to the parser using
			DEFAULT_SIMULATION_KWARGS[key] for the default value.
			"""
			default = DEFAULT_SIMULATION_KWARGS[key]
			self.define_option(parser, name, datatype, default, help)
		def add_bool_option(name, key, help):
			"""Add a boolean option parameter with the given name to the parser
			using DEFAULT_SIMULATION_KWARGS[key] for the default value. The CLI
			input can be `--name`, `--no_name`, `--name true`, `--name false`,
			`--name 1`, `--name 0`, `--name=true`, etc.
			"""
			self.define_parameter_bool(
				parser, name, DEFAULT_SIMULATION_KWARGS[key], help)

		self.define_option(parser, 'description', str, '',
			help='A simulation description to append to the output folder name.')
		self.define_parameter_bool(parser, 'verbose', True,
			help='Verbose workflow builder logging')
		parser.add_argument('-c', '--cpus', type=int, default=2,
			help='The number of CPU processes to use in the Parca and analysis'
				 ' steps. Default = 2.')
		self.define_parameter_bool(parser, 'dump', False,
			help='Dump the built workflow to JSON files for'
				 ' review *instead* of sending them to the Gaia workflow'
				 ' server. This is useful for testing and debugging. You can'
				 ' upload them manually or re-run this program without `--dump`.')
		parser.add_argument('-w', '--workers', type=int,
			help='The number of worker nodes to launch, with a smart default.')

		# Parca
		self.define_parameter_bool(parser, 'ribosome_fitting', True,
			help="Fit ribosome expression to protein synthesis demands")
		self.define_parameter_bool(parser, 'rnapoly_fitting', True,
			help="Fit RNA polymerase expression to protein synthesis demands")
		self.define_parameter_bool(parser, 'debug_parca', False,
			help='Make Parca calculate only one arbitrarily-chosen transcription'
				 ' factor condition when adjusting gene expression levels, leaving'
				 ' the other TFs at their input levels for faster Parca debugging.'
				 ' DO NOT USE THIS FOR A MEANINGFUL SIMULATION.')

		# Variant
		parser.add_argument('-v', '--variant', nargs=3, default=DEFAULT_VARIANT,
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='''The variant type name, first index, and last index to make.
				See models/ecoli/sim/variants/__init__.py for the variant
				type choices and their supported index ranges, e.g.: wildtype,
				condition, meneParams, metabolism_kinetic_objective_weight,
				nutrientTimeSeries, and param_sensitivity.
				Default = wildtype 0 0''')

		# Simulation
		parser.add_argument('-g', '--generations', type=int, default=1,
			help='Number of cell generations to run. Set it to 0 to just run'
				 ' Parca and make-variants with no sim generations or analysis.'
				 ' Default = 1')
		parser.add_argument('-i', '--init_sims', type=int, default=1,
			help='(int; 1) Number of initial sims (seeds) per variant.'
				 ' Default = 1')
		parser.add_argument('-t', '--timeline', type=str, default='0 minimal',
			help='set timeline. Default = "0 minimal". See'
				 ' environment/condition/make_media.py, make_timeline() for'
				 ' timeline formatting details')
		add_option('length_sec', 'lengthSec', int,
			help='The maximum simulation time, in seconds. Useful for short'
				 ' simulations; not so useful for multiple generations.'
				 ' Default is 3 hours')
		add_option('timestep_safety_frac', 'timeStepSafetyFraction', float,
			help='Scale the time step by this factor if conditions are'
				 ' favorable, up the the limit of the max time step')
		add_option('timestep_max', 'maxTimeStep', float,
			help='the maximum time step, in seconds')
		add_option('timestep_update_freq', 'updateTimeStepFreq', int,
			help='frequency at which the time step is updated')
		add_bool_option('mass_distribution', 'massDistribution',
			help='If true, a mass coefficient is drawn from a normal distribution'
				 ' centered on 1; otherwise it is set equal to 1')
		add_bool_option('growth_rate_noise', 'growthRateNoise',
			help='If true, a growth rate coefficient is drawn from a normal'
				 ' distribution centered on 1; otherwise it is set equal to 1')
		add_bool_option('d_period_division', 'dPeriodDivision',
			help='If true, ends simulation once D period has occurred after'
				 ' chromosome termination; otherwise simulation terminates once'
				 ' a given mass has been added to the cell')
		add_bool_option('translation_supply', 'translationSupply',
			help='If true, the ribosome elongation rate is limited by the'
				 ' condition specific rate of amino acid supply; otherwise the'
				 ' elongation rate is set by condition')
		add_bool_option('trna_charging', 'trna_charging',
			help='if true, tRNA charging reactions are modeled and the ribosome'
				 ' elongation rate is set by the amount of charged tRNA	present.'
				 ' This option will override TRANSLATION_SUPPLY in the simulation.')

		# TODO(jerry): Single/dual daughters.

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
		if args.dump:
			wf.write()
		else:
			wf.send(args.workers)


if __name__ == '__main__':
	script = RunWcm()
	script.cli()
