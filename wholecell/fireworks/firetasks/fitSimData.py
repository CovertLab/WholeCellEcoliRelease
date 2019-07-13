from __future__ import absolute_import, division, print_function

import cPickle
import time
import os
import shutil
import sys

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1


@explicit_serialize
class FitSimDataTask(FireTaskBase):

	_fw_name = "FitSimDataTask"
	required_params = [
		"cached",
		"debug",
		"input_data", "output_data", "cpus",
		"disable_ribosome_capacity_fitting",
		"disable_rnapoly_capacity_fitting",
		]
	optional_params = [
		"cached_data",
		"sim_out_dir",
		]

	def run_task(self, fw_spec):
		print("{}: Calculating sim_data parameters".format(time.ctime()))

		if self["cached"]:
			try:
				shutil.copyfile(self["cached_data"], self["output_data"])
				mod_time = time.ctime(os.path.getctime(self["cached_data"]))
				print("Copied sim data from cache (last modified {})".format(mod_time))
				return
			except Exception as exc:
				print("Warning: Could not copy cached sim data due to"
					  " exception ({}). Running Parca.".format(exc))

		cpus = self["cpus"]

		with open(self["input_data"], "rb") as f:
			raw_data = cPickle.load(f)

		sim_data = fitSimData_1(
			raw_data, cpus=cpus, debug=self["debug"],
			disable_ribosome_capacity_fitting=self['disable_ribosome_capacity_fitting'],
			disable_rnapoly_capacity_fitting=self['disable_rnapoly_capacity_fitting'],
			)

		sys.setrecursionlimit(4000) #limit found manually
		with open(self["output_data"], "wb") as f:
			cPickle.dump(sim_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
