"""
Run the analysis scripts that generate input files to the Causality Network
tool.
"""
from __future__ import absolute_import, division, print_function

import importlib
import os
import time
import traceback

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.analysis.causality_network.build_network import BuildNetwork
from models.ecoli.analysis.causality_network.network_components import NODELIST_FILENAME, DYNAMICS_FILENAME
import models.ecoli.analysis.causality_network


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
		]

	READER_FILE_PATH = 'models.ecoli.analysis.causality_network.read_dynamics'

	def plotter_args(self):
		self["metadata"] = dict(self["metadata"], analysis_type = "causality_network")

		return (
			self["input_results_directory"],
			self["output_dynamics_directory"],
			self['output_filename_prefix'] + DYNAMICS_FILENAME,
			self["input_sim_data"],
			self["node_list_file"],
			self["metadata"],
			)

	def run_task(self, fw_spec):
		startTime = time.time()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(startTime), type(self).__name__))

		self['node_list_file'] = os.path.join(
			self["output_network_directory"], NODELIST_FILENAME)

		self["check_sanity"] = self.get("check_sanity", False)

		if not os.path.isfile(self['node_list_file']):
			print("{}: Building causality network".format(time.ctime()))

			causality_network = BuildNetwork(
				self["input_sim_data"], self["output_network_directory"],
				self["check_sanity"])
			causality_network.run()

		self['output_filename_prefix'] = self.get('output_filename_prefix', '')

		mod = importlib.import_module(self.READER_FILE_PATH)
		args = self.plotter_args()

		print("{}: Reading simulation results for causality network"
			.format(time.ctime()))
		mod.Plot.main(*args)

		timeTotal = time.time() - startTime

		duration = time.strftime("%H:%M:%S", time.gmtime(timeTotal))
		print("{}: Completed building causality network in {}".format(
			time.ctime(), duration)
			)
