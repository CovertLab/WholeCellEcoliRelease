"""
Make one or more sim_data variants via VariantSimDataTask, writing e.g.
`wildtype_000000/kb/simData_Modified.cPickle`, in preparation for running
cell simulations via `runSim --require_variants`. (Without --require_variants,
runSim will make the sim_data variants, which is handy until you want to launch
multiple first-gen runSim runs in parallel without collisions writing the same
`simData_Modified.cPickle` file.)

See models/ecoli/sim/variants/*.py for the variant choices.

Prerequisite: Run the parameter calculator (runParca.py).

TODO(jerry): Write the SIM_PATH/metadata.json file here instead of in runSim?
Write it in runParca, then update it here and in runSim?

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import os

from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp

DEFAULT_VARIANT = ['wildtype', '0', '0']


class MakeVariants(scriptBase.ScriptBase):
	"""Make sim_data variants like wildtype_000000/kb/simData_Modified.cPickle."""

	def help(self):
		"""Return help text for the Command Line Interface."""
		return '''Run {}. Given a variant type name and an index range,
			this makes the variant subdirectories and their
			simData_Modified.cPickle files.'''.format(self.description())

	def define_parameters(self, parser):
		super(MakeVariants, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		parser.add_argument('-v', '--variant', nargs=3, default=DEFAULT_VARIANT,
			metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'),
			help='''The variant type name, first index, and last index to make.
				See models/ecoli/sim/variants/__init__.py for the variant
				type choices and their supported index ranges, e.g.: wildtype,
				condition, meneParams, metabolism_kinetic_objective_weight,
				nutrientTimeSeries, and param_sensitivity.
				Default = wildtype 0 0''')

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, constants.KB_DIR)
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		fp.verify_file_exists(sim_data_file, 'Run runParca?')

		variant_arg = args.variant
		variant_spec = (variant_arg[0], int(variant_arg[1]), int(variant_arg[2]))
		variant_type = variant_spec[0]

		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i, subdir in fp.iter_variants(*variant_spec):
			variant_sim_data_directory = os.path.join(args.sim_path, subdir,
				constants.VKB_DIR)
			variant_metadata_directory = os.path.join(args.sim_path, subdir,
				constants.METADATA_DIR)

			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			task = VariantSimDataTask(
				variant_function=variant_type,
				variant_index=i,
				input_sim_data=sim_data_file,
				output_sim_data=variant_sim_data_modified_file,
				variant_metadata_directory=variant_metadata_directory,
				)
			task.run_task({})


if __name__ == '__main__':
	script = MakeVariants()
	script.cli()
