from __future__ import absolute_import, division, print_function

import os

from fireworks import FiretaskBase, explicit_serialize

from . import InitRawDataTask
from . import InitRawValidationDataTask
from . import InitValidationDataTask
from . import FitSimDataTask

from wholecell.utils import constants
from wholecell.utils import filepath as fp
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


@explicit_serialize
class ParcaTask(FiretaskBase):
	"""A complete Parameter Calculator Firetask. It makes its output directory
	and writes everything into it, to fit into a cloud workflow.
	"""

	_fw_name = 'ParcaTask'
	required_params = [
		'output_directory',  # e.g. 'out/manual/kb'
		'ribosome_fitting',
		'rnapoly_fitting',
		'optimize_trna_charging_kinetics',
		]
	optional_params = [
		'operons',
		'load_intermediate',
		'save_intermediates',
		'intermediates_directory',
		'cpus',
		'debug',
		'variable_elongation_transcription',
		'variable_elongation_translation',
		]

	def _get_default(self, key):
		return self.get(key, DEFAULT_SIMULATION_KWARGS[key])

	def run_task(self, fw_spec):
		kb_directory = fp.makedirs(self['output_directory'])
		raw_data_file = os.path.join(kb_directory, constants.SERIALIZED_RAW_DATA)
		sim_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		metrics_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_METRICS_DATA_FILENAME)
		raw_validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_RAW_VALIDATION_DATA)
		validation_data_file = os.path.join(
			kb_directory, constants.SERIALIZED_VALIDATION_DATA)
		intermediates_dir = self.get('intermediates_directory') or kb_directory

		tasks = [
			InitRawDataTask(
				operons=self.get('operons'),
				output=raw_data_file),

			FitSimDataTask(
				input_data=raw_data_file,
				output_data=sim_data_file,
				output_metrics_data=metrics_data_file,
				cached=False,
				load_intermediate=self.get('load_intermediate', None),
				save_intermediates=self.get('save_intermediates', False),
				intermediates_directory=intermediates_dir,
				cpus=self.get('cpus', 1),
				debug=self.get('debug', False),
				variable_elongation_transcription=self._get_default('variable_elongation_transcription'),
				variable_elongation_translation=self._get_default('variable_elongation_translation'),
				disable_ribosome_capacity_fitting=not self['ribosome_fitting'],
				disable_rnapoly_capacity_fitting=not self['rnapoly_fitting'],
				optimize_trna_charging_kinetics=self['optimize_trna_charging_kinetics']),

			InitRawValidationDataTask(
				output=raw_validation_data_file),

			InitValidationDataTask(
				validation_data_input=raw_validation_data_file,
				knowledge_base_raw=raw_data_file,
				output_data=validation_data_file),
			]
		for task in tasks:
			task.run_task(fw_spec)

		print(f'Wrote: {kb_directory}')
