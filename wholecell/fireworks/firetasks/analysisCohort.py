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
import multiprocessing as mp

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

		if "WC_ANALYZE_FAST" in os.environ:
			pool = mp.Pool(processes = 8)

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			mod = importlib.import_module("models.ecoli.analysis.cohort." + f[:-3])
			args = (
				self["input_variant_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_sim_data"],
				self['input_validation_data'],
				self["metadata"],
				)

			if "WC_ANALYZE_FAST" in os.environ:
				pool.apply_async(run_function, args = (mod.main, args, f))
			else:
				print "%s: Running %s" % (time.ctime(), f)
				mod.main(*args)

		if "WC_ANALYZE_FAST" in os.environ:
			pool.close()
			pool.join()
		timeTotal = time.time() - startTime
		print "Completed cohort analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)
