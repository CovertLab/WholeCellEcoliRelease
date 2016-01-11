"""
AnalysisVariantTask

Analyzes across variants. Has access to all cells in the entire simulation run.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/06/2015
"""

import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.variant
import importlib

@explicit_serialize
class AnalysisVariantTask(FireTaskBase):

	_fw_name = "AnalysisVariantTask"
	required_params = [
		"input_directory",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		startTime = time.time()
		print "%s: Running variant analysis" % time.ctime(startTime)

		directory = os.path.dirname(models.ecoli.analysis.variant.__file__)

		# Run analysis scripts in order of modification, most recently edited first
		fileList = os.listdir(directory)
		fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.variant." + f[:-3])
			mod.main(
				inputDir = self["input_directory"],
				plotOutDir = self["output_plots_directory"],
				plotOutFileName = f[:-3],
				validationDataFile = self['input_validation_data'],
				metadata = self["metadata"]
				)

		timeTotal = time.time() - startTime
		print "Completed variant analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))