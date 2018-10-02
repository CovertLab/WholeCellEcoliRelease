from __future__ import absolute_import, division, print_function

import os

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.analysis.causality_network.build_network import BuildNetwork
from models.ecoli.analysis.causality_network.network_components import NODELIST_FILENAME


@explicit_serialize
class BuildCausalityNetworkTask(FireTaskBase):

	_fw_name = "BuildCausalNetworkTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"output_network_directory",
		"output_dynamics_directory",
		"metadata",
		]
	optional_params = [
		"output_filename_prefix",
		"check_sanity",
		"cpus",
		]

	def run_task(self, fw_spec):
		if not os.path.isfile(
				os.path.join(self["output_network_directory"],
					NODELIST_FILENAME)
				):
			causality_network = BuildNetwork(
				self["input_sim_data"], self["output_network_directory"],
				self["check_sanity"])
			causality_network.run()
