"""
Common code for analysis plots. The abstract base class AnalysisPlot defines a
plot() method for scripts to call.

TODO: Enable future warnings, esp. for matplotlib.

TODO: Move the run_plot() args to instance variables?

TODO: Other shared code to simplify the subclasses, e.g. instantiate an
AnalysisPaths (except for SingleAnalysisPlot subclasses), etc.
"""

from __future__ import annotations

import abc
import os
import pickle
from typing import Tuple

import matplotlib as mp
from matplotlib import pyplot as plt
import numpy as np

from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.utils import constants, memory_debug, parallelization
from wholecell.utils import filepath as fp


class AnalysisPlot(metaclass=abc.ABCMeta):
	"""Abstract Base Class for analysis plots.

	Each analysis class must override do_plot().

	Call main() to run an analysis plot for a Firetask.

	Use the environment variable 'DEBUG_GC' to enable memory leak debugging.

	Inputs:
		cpus: allotted number of CPU cores; default (0) => all available cores
	"""

	#: Whether to suppress NumPy "divide" and "invalid" warnings like Theano
	#: used to do. Override this in subclasses as needed.
	_suppress_numpy_warnings = False

	def __init__(self, cpus=0):
		self.cpus = parallelization.cpus(cpus)
		self._axeses = {}
		self.ap = None

	@staticmethod
	def read_sim_data_file(sim_path: str) -> SimulationDataEcoli:
		"""Return the sim_data object read from sim_path."""
		simDataFile = os.path.join(sim_path,
								   constants.KB_DIR,
								   constants.SERIALIZED_SIM_DATA_FILENAME)
		with open(simDataFile, 'rb') as f:
			return pickle.load(f)

	@staticmethod
	def read_validation_data_file(sim_path: str) -> ValidationDataEcoli:
		"""Return the validation_data object read from sim_path."""
		validationDataFile = os.path.join(sim_path,
										  constants.KB_DIR,
										  constants.SERIALIZED_VALIDATION_DATA)
		with open(validationDataFile, 'rb') as f:
			return pickle.load(f)

	def subplot(self, *args):
		"""
		Create a subplot or return the axes previously created with the same
		args tuple. Either way, make it the current axes.

		The caller must use consistent args for each subplot, e.g. `(2, 2, 1)`
		vs. `(221,)`. Don't intermix them.

		This does not yet support keyword args.

		Use this to work around this matplotlib 2.2.2 deprecation when you want
		to retrieve a previously-created axes (subplot):

			MatplotlibDeprecationWarning: Adding an axes using the same
			arguments as a previous axes currently reuses the earlier instance.
			In a future version, a new instance will always be created and
			returned.  Meanwhile, this warning can be suppressed, and the
			future behavior ensured, by passing a unique label to each axes
			instance.

		Args:
			*args (Tuple[int]): subplot grid and index args

		Returns:
			(axes): the created or retrieved axes object; now the current one
		"""
		axes = self._axeses.get(args)
		if not axes:
			self._axeses[args] = axes = plt.subplot(*args)
		plt.sca(axes)
		return axes

	@staticmethod
	def set_ylim(axes, ymin, ymax):
		"""
		Set the axes y-limits, avoiding the matplotlib warning:

			UserWarning: Attempting to set identical bottom==top results

		Args:
			axes (axes): the axes to modify
			ymin (float): the minimum y value
			ymax (float): the maximum y value

		Returns:
			ylimits (Tuple[float]): the new y-axis limits as (`bottom`, `top`)
		"""
		if ymin == ymax:
			ymin -= 0.001
			ymax += 0.001
		return axes.set_ylim(ymin, ymax)

	@staticmethod
	def remove_border(ax=None, bottom=False):
		if ax is None:
			ax = plt.gca()

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		if bottom:
			ax.spines['bottom'].set_visible(False)
			ax.set_xticks([])

	@abc.abstractmethod
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method that each analysis class must override."""
		raise NotImplementedError("AnalysisPlot subclass must implement do_plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Public method to set up, make a plot, and cleanup."""
		def do_plot():
			self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
						 validationDataFile, metadata)

		if not os.path.isdir(inputDir):
			raise RuntimeError('Input directory ({}) does not currently exist.'
				.format(inputDir))
		fp.makedirs(plotOutDir)

		with memory_debug.detect_leaks(), mp.rc_context():
			if self._suppress_numpy_warnings:
				with np.errstate(divide='ignore'), np.errstate(invalid='ignore'):
					do_plot()
			else:
				do_plot()

		self._axeses = {}

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None, cpus=0, analysis_paths=None):
		"""Run an analysis plot for a Firetask."""
		instance = cls(cpus)
		instance.ap = analysis_paths
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
