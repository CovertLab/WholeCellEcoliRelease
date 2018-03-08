import cPickle
import time
import os
import shutil

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from reconstruction.ecoli.fit_sim_data_2 import fitSimData_2

@explicit_serialize
class FitSimDataTask(FireTaskBase):

	_fw_name = "FitSimDataTask"
	required_params = ["fit_level", "input_data", "output_data"]
	optional_params = ["sim_out_dir"]

	def run_task(self, fw_spec):

		print "%s: Creating/Fitting sim_data (Level %d)" % (time.ctime(), self["fit_level"])

		if self["fit_level"] == 1:
			if self["cached"]:
				try:
					shutil.copyfile(self["cached_data"], self["output_data"])
					print "Copied sim data from cache (modified %s)" % time.ctime(os.path.getctime(self["cached_data"]))
					return
				except Exception as exc:
					print "Warning: could not copy cached sim data due to exception (%s), running fitter" % (exc)

			if self["cpus"] > 1:
				print "Warning: running fitter in parallel with %i processes - ensure there are enough cpus_per_task allocated" % self["cpus"]

			raw_data = cPickle.load(open(self["input_data"], "rb"))
			sim_data = fitSimData_1(raw_data, self["cpus"])
			import sys; sys.setrecursionlimit(4000) #limit found manually
			cPickle.dump(
				sim_data,
				open(self["output_data"], "wb"),
				protocol = cPickle.HIGHEST_PROTOCOL
				)

		# TODO: Get rid of this if not used
		if self["fit_level"] == 2:
			sim_data = cPickle.load(open(self["input_data"], "rb"))
			fitSimData_2(sim_data, self["sim_out_dir"])
			cPickle.dump(
				sim_data,
				open(self["output_data"], "wb"),
				protocol = cPickle.HIGHEST_PROTOCOL
				)