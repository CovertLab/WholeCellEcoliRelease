"""
Run the Parameter Calculator (Parca).
The output goes into the named subdirectory of wcEcoli/out/, defaulting to "manual".

TODO: Share lots of code with fw_queue.py and AnalysisPaths.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os

from wholecell.fireworks.firetasks import ParcaTask
from wholecell.utils import parallelization, scriptBase
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
				 ' "Fast run" to "20190514.135600__Fast_run".')
		parser.add_argument('-c', '--cpus', type=int, default=1,
			help='The number of CPU processes to use. Default = 1.'
			)
		parser.add_argument('-d', '--debug', action='store_true',
			help="Enable Parca debugging: For a faster debug cycle, calculate"
				 " only one arbitrarily-chosen transcription factor condition"
				 " when adjusting gene expression levels, leaving the others at"
				 " their input levels. Do not use this for an actual simulation."
			)
		self.define_parameter_bool(parser, 'variable_elongation_transcription', False,
			help="Use a variable elongation rate for transcription")
		self.define_parameter_bool(parser, 'variable_elongation_translation', False,
			help="Use a variable elongation rate for translation")
		self.define_parameter_bool(parser, 'ribosome_fitting', True,
			help="Fit ribosome expression to protein synthesis demands")
		self.define_parameter_bool(parser, 'rnapoly_fitting', True,
			help="Fit RNA polymerase expression to protein synthesis demands")

	def parse_args(self):
		args = super(RunParca, self).parse_args()
		args.cpus = parallelization.cpus(args.cpus)

		if args.timestamp:
			args.sim_outdir = fp.timestamp() + '__' + args.sim_outdir.replace(' ', '_')

		args.sim_path = fp.makedirs(fp.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, ParcaTask.OUTPUT_SUBDIR)

		if args.debug:
			print("{}Parca".format('DEBUG ' if args.debug else ''))

		task = ParcaTask(
			output_directory=kb_directory,
			ribosome_fitting=args.ribosome_fitting,
			rnapoly_fitting=args.rnapoly_fitting,
			variable_elongation_transcription=args.variable_elongation_transcription,
			variable_elongation_translation=args.variable_elongation_translation,
			cpus=args.cpus,
			debug=args.debug)
		task.run_task({})


if __name__ == '__main__':
	script = RunParca()
	script.cli()
