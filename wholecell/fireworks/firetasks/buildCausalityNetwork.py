"""
Run the analysis scripts that generate input files to the Causality Network
tool.
"""
from __future__ import absolute_import, division, print_function

import datetime
import time

from fireworks import FiretaskBase, explicit_serialize
from models.ecoli.analysis.causality_network import read_dynamics
from models.ecoli.analysis.causality_network.build_network import BuildNetwork
from wholecell.utils import filepath as fp
from wholecell.utils.py3 import monotonic_seconds


@explicit_serialize
class BuildCausalityNetworkTask(FiretaskBase):

	_fw_name = "BuildCausalNetworkTask"

	required_params = [
		"input_results_directory",
		"input_sim_data",
		"output_dynamics_directory",
		]
	optional_params = [
		"check_sanity",
		"output_network_directory",  # no longer used
		"metadata",  # no longer used
		]

	def run_task(self, fw_spec):
		start_real_sec = monotonic_seconds()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(), type(self).__name__))

		print("{}: Building the Causality network".format(time.ctime()))
		causality_network = BuildNetwork(
			self["input_sim_data"],
			self["output_dynamics_directory"],
			self.get("check_sanity", False))
		node_list, edge_list = causality_network.build_nodes_and_edges()

		fp.makedirs(self["output_dynamics_directory"])

		print("{}: Converting simulation results to a Causality series"
			.format(time.ctime()))
		read_dynamics.convert_dynamics(
			self["input_results_directory"],
			self["output_dynamics_directory"],
			self["input_sim_data"],
			node_list,
			edge_list)

		elapsed_real_sec = monotonic_seconds() - start_real_sec

		duration = datetime.timedelta(seconds=elapsed_real_sec)
		print("{}: Completed building the Causality network in {}".format(
			time.ctime(), duration))
