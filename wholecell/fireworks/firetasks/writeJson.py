from __future__ import absolute_import, division, print_function

import os

from fireworks import FireTaskBase, explicit_serialize

from wholecell.utils import filepath as fp


@explicit_serialize
class WriteJsonTask(FireTaskBase):
	"""A Firetask to write a JSON file."""

	_fw_name = 'JsonTask'
	required_params = [
		'output_file',  # e.g. 'out/manual/metadata/metadata.json'
		'data']

	def run_task(self, fw_spec):
		output_file = self['output_file']
		data = self['data']

		output_dir = os.path.basename(output_file)
		fp.makedirs(output_dir)
		fp.write_json_file(output_file, data)
