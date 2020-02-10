from __future__ import absolute_import, division, print_function

import os

from fireworks import FiretaskBase, explicit_serialize

from wholecell.utils import filepath as fp
from wholecell.utils import data


@explicit_serialize
class WriteJsonTask(FiretaskBase):
	"""A Firetask to write a JSON file with environment $VARIABLE expansion."""

	_fw_name = 'WriteJsonTask'
	required_params = [
		'output_file',  # e.g. 'out/manual/metadata/metadata.json'
		'data']

	def run_task(self, fw_spec):
		output_file = self['output_file']
		data1 = data.expand_keyed_env_vars(self['data'])

		output_dir = os.path.dirname(output_file)
		fp.makedirs(output_dir)
		fp.write_json_file(output_file, data1)
