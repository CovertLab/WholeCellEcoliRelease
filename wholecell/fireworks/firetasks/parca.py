from __future__ import absolute_import, division, print_function

import os

from fireworks import FireTaskBase, explicit_serialize

from . import InitRawDataTask
from . import InitRawValidationDataTask
from . import InitValidationDataTask
from . import FitSimDataTask

from wholecell.utils import constants
from wholecell.utils import filepath as fp

from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

@explicit_serialize
class ParcaTask(FireTaskBase):
	"""A complete Parameter Calculator Firetask. It makes its output directory
	and writes everything into it, to fit into a Gaia/Sisyphus workflow.
	"""

	_fw_name = 'ParcaTask'
	required_params = [
		'output_directory',  # e.g. 'out/manual/kb'
		'ribosome_fitting',
		'rnapoly_fitting']
	optional_params = [
		'cpus',
		'debug',
		'variable_elongation_transcription',
		'variable_elongation_translation']

	OUTPUT_SUBDIR = 'kb'  # the task's recommended output directory

	def run_task(self, fw_spec):
		kb_directory = fp.makedirs(self['output_directory'])
		raw_data_file = os.path.join(kb_directory, constants.SERIALIZED_RAW_DATA)
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		raw_validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_RAW_VALIDATION_DATA)
		validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_VALIDATION_DATA)

		tasks = [
			InitRawDataTask(
				output=raw_data_file),

			FitSimDataTask(
				input_data=raw_data_file,
				output_data=sim_data_file,
				cached=False,
				cpus=self.get('cpus', 1),
				debug=self.get('debug', False),
				variable_elongation_transcription=self.get('variable_elongation_transcription', DEFAULT_SIMULATION_KWARGS['variable_elongation_transcription']),
				variable_elongation_translation=self.get('variable_elongation_translation', DEFAULT_SIMULATION_KWARGS['variable_elongation_translation']),
				disable_ribosome_capacity_fitting=not self['ribosome_fitting'],
				disable_rnapoly_capacity_fitting=not self['rnapoly_fitting']),

			InitRawValidationDataTask(
				output=raw_validation_data_file),

			InitValidationDataTask(
				validation_data_input=raw_validation_data_file,
				knowledge_base_raw=raw_data_file,
				output_data=validation_data_file),
			]
		for task in tasks:
			task.run_task(fw_spec)

		print('\n\t'.join(['Wrote', raw_data_file, sim_data_file,
			raw_validation_data_file, validation_data_file]))
