import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.multigen
import importlib

@explicit_serialize
class AnalysisMultiGenTask(FireTaskBase):

	_fw_name = "AnalysisMultiGenTask"
	required_params = [
		"input_seed_directory",
		"input_kb",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		print "%s: Running multiple generation analysis" % time.ctime()

		directory = os.path.dirname(models.ecoli.analysis.multigen.__file__)

		fileList = sorted(os.listdir(directory))

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.multigen." + f[:-3])
			mod.main(
				self["input_seed_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_kb"],
				self["metadata"]
				)