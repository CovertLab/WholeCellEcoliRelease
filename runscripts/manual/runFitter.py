"""
Run the Fitter. The output goes into the named subdirectory of wcEcoli/out/,
defaulting to "manual".

TODO: Default to a typical timestamped directory name instead of "manual"?

TODO: Share lots of code with fw_queue.py and AnalysisPaths.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import os

from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.fireworks.firetasks import SymlinkTask
from wholecell.utils import constants, scriptBase
from wholecell.utils import filepath


class RunFitter(scriptBase.ScriptBase):
	"""Runs the Fitter, aka simulation parameter calculator."""

	def define_parameters(self, parser):
		super(RunFitter, self).define_parameters(parser)

		# NOTE: Don't name this arg sim_dir since that makes parse_args() look
		# for an existing sim_dir directory.
		parser.add_argument('sim_outdir', nargs='?', default='manual',
			help='The simulation "out/" subdirectory to write to.'
				 ' Default = "manual".')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. Default = 1.')
		parser.add_argument('--cached', action='store_true',
			help='Copy a cached "' + constants.SERIALIZED_FIT1_FILENAME
				 + '" file instead of generating it.')
		parser.add_argument('-d', '--debug', action='store_true',
			help="Enable Fitter debugging: Fit only one arbitrarily-chosen"
				 " transcription factor for a faster debug cycle. Don't use it"
				 " for an actual simulation.")

	def parse_args(self):
		args = super(RunFitter, self).parse_args()
		args.sim_path = filepath.makedirs(
			scriptBase.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		kb_directory = filepath.makedirs(args.sim_path, "kb")
		raw_data_file = os.path.join(kb_directory, constants.SERIALIZED_RAW_DATA)
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_FIT1_FILENAME)
		cached_sim_data_file = os.path.join(
			scriptBase.ROOT_PATH, 'cached', constants.SERIALIZED_FIT1_FILENAME)
		most_fit_filename = os.path.join(
			kb_directory, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME)
		raw_validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_RAW_VALIDATION_DATA)
		validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_VALIDATION_DATA)

		task1 = InitRawDataTask(
			output=raw_data_file,
		)
		task1.run_task({})

		if args.debug or args.cached:
			print "{}{}Fitter".format(
				'DEBUG ' if args.debug else '',
				'CACHED ' if args.cached else '',
			)

		task2 = FitSimDataTask(
			fit_level=1,
			input_data=raw_data_file,
			output_data=sim_data_file,
			cached=args.cached,  # bool
			cached_data=cached_sim_data_file,  # cached file to copy
			cpus=args.cpus,
			debug=args.debug,
		)
		task2.run_task({})

		task3 = SymlinkTask(
			to=constants.SERIALIZED_FIT1_FILENAME,
			link=most_fit_filename,
			overwrite_if_exists=True,
		)
		task3.run_task({})

		task4 = InitRawValidationDataTask(
			output=raw_validation_data_file,
		)
		task4.run_task({})

		task5 = InitValidationDataTask(
			validation_data_input=raw_validation_data_file,
			knowledge_base_raw=raw_data_file,
			output_data=validation_data_file,
		)
		task5.run_task({})

		print '\n\t'.join(['Wrote', raw_data_file, sim_data_file,
			most_fit_filename, raw_validation_data_file, validation_data_file])


if __name__ == '__main__':
	script = RunFitter()
	script.cli()