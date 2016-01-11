"""
AnalysisCohortTask

Analyzes all cells, all seeds, all generations.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
"""

import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.cohort
import importlib

@explicit_serialize
class AnalysisCohortTask(FireTaskBase):

	_fw_name = "AnalysisCohortTask"
	required_params = [
		"input_variant_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		startTime = time.time()
		print "%s: Running cohort analysis" % time.ctime(startTime)

		directory = os.path.dirname(models.ecoli.analysis.cohort.__file__)

		# Run analysis scripts in order of modification, most recently edited first
		fileList = os.listdir(directory)
		fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.cohort." + f[:-3])
			mod.main(
				variantDir = self["input_variant_directory"],
				plotOutDir = self["output_plots_directory"],
				plotOutFileName = f[:-3],
				simDataFile = self["input_sim_data"],
				validationDataFile = self['input_validation_data'],
				metadata = self["metadata"]
				)

		timeTotal = time.time() - startTime
		print "Completed cohort analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))
