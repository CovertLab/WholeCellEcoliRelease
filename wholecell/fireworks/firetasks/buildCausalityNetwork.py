"""
Run the analysis scripts that generate input files to the Causality Network
tool.
"""
from __future__ import absolute_import, division, print_function

import datetime
import importlib
import os
import time

from fireworks import FiretaskBase, explicit_serialize
from models.ecoli.analysis.causality_network.build_network import BuildNetwork
from models.ecoli.analysis.causality_network.network_components import NODELIST_JSON, DYNAMICS_FILENAME
from wholecell.utils import filepath as fp
from wholecell.utils import data
from wholecell.utils.py3 import monotonic_seconds


@explicit_serialize
class BuildCausalityNetworkTask(FiretaskBase):

	_fw_name = "BuildCausalNetworkTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"output_network_directory",  # an output once per variant, else an input!
		"output_dynamics_directory",
		"metadata",
		]
	optional_params = [
		"check_sanity",
		"force_update",
		]

	READER_FILE_PATH = 'models.ecoli.analysis.causality_network.read_dynamics'

	def plotter_args(self):
		self["metadata"] = dict(self["metadata"], analysis_type="causality_network")

		return (
			self["input_results_directory"],
			self["output_dynamics_directory"],
			DYNAMICS_FILENAME,
			self["input_sim_data"],
			self["node_list_file"],
			self["metadata"],
			)

	def run_task(self, fw_spec):
		start_real_sec = monotonic_seconds()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(), type(self).__name__))

		self['metadata'] = data.expand_keyed_env_vars(self['metadata'])
		self['node_list_file'] = os.path.join(
			self["output_network_directory"], NODELIST_JSON)
			# self["output_network_directory"], NODELIST_FILENAME)

		self["check_sanity"] = self.get("check_sanity", False)

		if self.get("force_update", False) or not os.path.isfile(self['node_list_file']):
			print("{}: Building causality network".format(time.ctime()))

			fp.makedirs(self["output_network_directory"])
			causality_network = BuildNetwork(
				self["input_sim_data"], self["output_network_directory"],
				self["check_sanity"])
			causality_network.run()

		fp.makedirs(self["output_dynamics_directory"])

		mod = importlib.import_module(self.READER_FILE_PATH)
		args = self.plotter_args()

		print("{}: Reading simulation results for causality network"
			.format(time.ctime()))
		mod.Plot.main(*args)

		elapsed_real_sec = monotonic_seconds() - start_real_sec

		duration = datetime.timedelta(seconds=elapsed_real_sec)
		print("{}: Completed building causality network in {}".format(
			time.ctime(), duration)
			)
