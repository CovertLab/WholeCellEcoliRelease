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
		"cpus",
		]

	MODULE_PATH = 'models.ecoli.analysis.causality_network'
	READER_FILE = models.ecoli.analysis.causality_network.READER

	def plotter_args(self):
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

		exceptionFileList = []
		self['node_list_file'] = os.path.join(
			self["output_network_directory"], NODELIST_FILENAME)

		if not os.path.isfile(self['node_list_file']):
			print("{}: Building causality network".format(time.ctime()))

			try:
				causality_network = BuildNetwork(
					self["input_sim_data"], self["output_network_directory"],
					self["check_sanity"])
				causality_network.run()
			except Exception:
				traceback.print_exc()
				exceptionFileList.append("build_network.py")

		self['output_filename_prefix'] = self.get('output_filename_prefix', '')

		mod = importlib.import_module(
			self.MODULE_PATH + "." + self.READER_FILE)
		args = self.plotter_args()

		print("{}: Reading simulation results for causality network"
			.format(time.ctime()))
		try:
			mod.Plot.main(*args)
		except Exception:
			traceback.print_exc()
			exceptionFileList.append("read_dynamics.py")

		timeTotal = time.time() - startTime

		duration = time.strftime("%H:%M:%S", time.gmtime(timeTotal))
		if exceptionFileList:
			print("Completed analysis in {} with an exception in:".format(duration))
			for file in exceptionFileList:
				print("\t{}".format(file))
			raise Exception("Error in analysis")
		else:
			print("Completed analysis in {}".format(duration))