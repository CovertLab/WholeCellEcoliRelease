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
		"variant_directory",
		"input_kb",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		print "%s: Running cohort analysis" % time.ctime()

		directory = os.path.dirname(models.ecoli.analysis.cohort.__file__)

		fileList = sorted(os.listdir(directory))

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			print "%s: Running %s" % (time.ctime(), f)

			mod = importlib.import_module("models.ecoli.analysis.cohort." + f[:-3])
			mod.main(
				variantDir = self["variant_directory"],
				plotOutDir = self["output_plots_directory"],
				plotOutFileName = f[:-3],
				kbFile = self["input_kb"],
				metadata = self["metadata"]
				)

# def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):