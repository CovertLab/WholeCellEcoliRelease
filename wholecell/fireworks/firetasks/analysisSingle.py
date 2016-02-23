import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.single
import importlib
import multiprocessing as mp

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

		# Run files in runlast.txt after others
		if "runlast.txt" in fileList:
			fileList = run_last(directory, fileList, "runlast.txt", verbose=False)
		
		self.analyze_fast = "WC_ANALYZE_FAST" in os.environ

		if self.analyze_fast:
			self.analyse_fast = True
			pool = mp.Pool(processes = 8)

		for f in fileList:
			if f == "__init__.py" or not f.endswith(".py"):
				continue

			mod = importlib.import_module("models.ecoli.analysis.single." + f[:-3])
			args = (
				self["input_results_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_sim_data"],
				self["input_validation_data"],
				self["metadata"],
				)

			if self.analyze_fast:
				pool.apply_async(run_function, args = (mod.main, args, f))
			else:
				print "%s: Running %s" % (time.ctime(), f)
				mod.main(*args)

		if self.analyze_fast:
			pool.close()
			pool.join()
		timeTotal = time.time() - startTime
		print "Completed single simulation analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))


def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)

def run_last(directory, fileList, fileName, verbose=False):
	"""
	Moves files mentioned in fileName to end of fileList.
	"""
	with open(os.path.join(directory, "runlast.txt"), 'r') as f:
		runlast = []
		for line in f:
			line = line.strip("\n")
			if line in fileList:
				runlast.append(line)
				idx = fileList.index(line)
				fileList.append(fileList.pop(idx))
			else:
				if verbose:
					print "\nUnknown file in runlast.txt: %s" % (line)
	if verbose:
		print "Running files in runlast.txt after others: \n%s\n" % (", ".join(runlast))
	return fileList