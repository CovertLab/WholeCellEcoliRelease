from __future__ import absolute_import
from __future__ import division

import cPickle
import time
import os
import shutil
import sys

from fireworks import FireTaskBase, explicit_serialize
from wholecell.utils import parallelization
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from reconstruction.ecoli.fit_sim_data_2 import fitSimData_2

@explicit_serialize
class FitSimDataTask(FireTaskBase):

	_fw_name = "FitSimDataTask"
	required_params = ["fit_level", "input_data", "output_data", "cpus"]
	optional_params = [
		"sim_out_dir",
		"disable_ribosome_capacity_fitting",
		"disable_rnapoly_capacity_fitting"
		]

	def run_task(self, fw_spec):

		print "%s: Creating/Fitting sim_data (Level %d)" % (time.ctime(), self["fit_level"])

		if self["fit_level"] == 1:
			if self["cached"]:
				try:
					shutil.copyfile(self["cached_data"], self["output_data"])
					mod_time = time.ctime(os.path.getctime(self["cached_data"]))
					print "Copied sim data from cache (modified %s)" % (mod_time,)
					return
				except Exception as exc:
					print ("Warning: could not copy cached sim data due to"
						   " exception (%s), running fitter") % (exc,)

			cpus = min(self["cpus"], parallelization.cpus())
			if cpus > 1:
				print ("Warning: running fitter in parallel with %i processes -"
					   " ensure there are enough cpus_per_task allocated" % (cpus,))

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

		# TODO: Get rid of this if not used
		if self["fit_level"] == 2:
			with open(self["input_data"], "rb") as f:
				sim_data = cPickle.load(f)

			fitSimData_2(sim_data, self["sim_out_dir"])

			with open(self["output_data"], "wb") as f:
				cPickle.dump(sim_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
