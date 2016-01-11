import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.single
import importlib

@explicit_serialize
class AnalysisSingleTask(FireTaskBase):

	_fw_name = "AnalysisSingleTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		startTime = time.time()
		print "%s: Running single simulation analysis" % time.ctime(startTime)

		directory = os.path.dirname(models.ecoli.analysis.single.__file__)

		# Run analysis scripts in order of modification, most recently edited first
		fileList = os.listdir(directory)
		fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.single." + f[:-3])
			mod.main(
				self["input_results_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_sim_data"],
				self["input_validation_data"],
				self["metadata"],
				)

		timeTotal = time.time() - startTime
		print "Completed single simulation analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))
