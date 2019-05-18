"""
Run the Parameter Calculator (Parca).
The output goes into the named subdirectory of wcEcoli/out/, defaulting to "manual".

TODO: Share lots of code with fw_queue.py and AnalysisPaths.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os

from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.utils import constants, scriptBase
from wholecell.utils import filepath as fp


class RunParca(scriptBase.ScriptBase):
	"""Runs the Parca, aka simulation parameter calculator."""

	def define_parameters(self, parser):
		super(RunParca, self).define_parameters(parser)

		# NOTE: Don't name this arg sim_dir since that makes parse_args() look
		# for an existing sim_dir directory while here we aim to create one.
		parser.add_argument('sim_outdir', nargs='?', default='manual',
			help='The simulation "out/" subdirectory to write to.'
				 ' Default = "manual".'
			)

		parser.add_argument('--timestamp', action='store_true',
			help='Timestamp the given `sim_outdir`, transforming e.g.'
				 ' "Fast run" to "20190514.135600.123456__Fast_run".')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. Default = 1.'
			)
		parser.add_argument('--cached', action='store_true',
			help='Copy a cached "' + constants.SERIALIZED_SIM_DATA_FILENAME
				 + '" file instead of generating it.'
			)
		parser.add_argument('-d', '--debug', action='store_true',
			help="Enable Parca debugging: For a faster debug cycle, calculate"
				 " only one arbitrarily-chosen transcription factor condition"
				 " when adjusting gene expression levels, leaving the others at"
				 " their input levels. Do not use this for an actual simulation."
			)
		self.define_parameter_bool(parser, 'ribosome_fitting', True,
			help="Fit ribosome expression to protein synthesis demands"
			)
		self.define_parameter_bool(parser, 'rnapoly_fitting', True,
			help="Fit RNA polymerase expression to protein synthesis demands"
			)

	def parse_args(self):
		args = super(RunParca, self).parse_args()

		if args.timestamp:
			args.sim_outdir = fp.timestamp() + '__' + args.sim_outdir.replace(' ', '_')

		args.sim_path = fp.makedirs(fp.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		kb_directory = fp.makedirs(args.sim_path, "kb")
		raw_data_file = os.path.join(kb_directory, constants.SERIALIZED_RAW_DATA)
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		cached_sim_data_file = os.path.join(
			fp.ROOT_PATH, 'cached', constants.SERIALIZED_SIM_DATA_FILENAME)
		raw_validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_RAW_VALIDATION_DATA)
		validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_VALIDATION_DATA)

		if args.debug or args.cached:
			print("{}{}Parca".format(
				'DEBUG ' if args.debug else '',
				'CACHED ' if args.cached else ''))

		tasks = [
			InitRawDataTask(
				output=raw_data_file,
				),

			FitSimDataTask(
				input_data=raw_data_file,
				output_data=sim_data_file,
				cached=args.cached,  # bool
				cached_data=cached_sim_data_file,  # cached file to copy
				cpus=args.cpus,
				debug=args.debug,
				disable_ribosome_capacity_fitting=not args.ribosome_fitting,
				disable_rnapoly_capacity_fitting=not args.rnapoly_fitting
				),

			InitRawValidationDataTask(
				output=raw_validation_data_file,
				),

			InitValidationDataTask(
				validation_data_input=raw_validation_data_file,
				knowledge_base_raw=raw_data_file,
				output_data=validation_data_file,
				),
			]
		for task in tasks:
			task.run_task({})

		print('\n\t'.join(['Wrote', raw_data_file, sim_data_file,
			raw_validation_data_file, validation_data_file]))


if __name__ == '__main__':
	script = RunParca()
	script.cli()
