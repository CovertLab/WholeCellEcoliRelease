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
		"input_kb",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		print "%s: Running single simulation analysis" % time.ctime()

		directory = os.path.dirname(models.ecoli.analysis.single.__file__)

		fileList = sorted(os.listdir(directory))

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.single." + f[:-3])
			mod.main(
				self["input_results_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_kb"],
				self["metadata"],
				)