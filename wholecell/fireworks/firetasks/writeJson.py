from __future__ import absolute_import, division, print_function

import os

from fireworks import FiretaskBase, explicit_serialize
from six import string_types

from wholecell.utils import filepath as fp


@explicit_serialize
class WriteJsonTask(FiretaskBase):
	"""A Firetask to write a JSON file with environment $VARIABLE expansion."""

	_fw_name = 'JsonTask'
	required_params = [
		'output_file',  # e.g. 'out/manual/metadata/metadata.json'
		'data']

	def run_task(self, fw_spec):
		def expand(var):
			"""Expand an environment $VARIABLE, if it is one."""
			if isinstance(var, string_types) and var.startswith('$'):
				return os.environ.get(var[1:], var)
			return var

		output_file = self['output_file']
		data = self['data']

		if isinstance(data, dict):
			data = {k: expand(v) for k, v in data.items()}

		output_dir = os.path.dirname(output_file)
		fp.makedirs(output_dir)
		fp.write_json_file(output_file, data)
