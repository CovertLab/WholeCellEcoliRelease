"""Run a shell script in a Docker container."""

from __future__ import absolute_import, division, print_function

from fireworks import FiretaskBase, explicit_serialize


@explicit_serialize
class DockerTask(FiretaskBase):
	_fw_name = "DockerTask"

	required_params = [
		'name',
		'image',
		'command',
		'storage_prefix',
		'internal_prefix']
	optional_params = [
		'inputs',
		'outputs',
		'timeout']

	def run_task(self, fw_spec):
		raise NotImplementedError(
			'This is a placeholder. The implementation is in the'
			' borealis-fireworks repo, which should become a pip.')
