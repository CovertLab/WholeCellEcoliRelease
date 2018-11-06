"""
Run the plots for a given list of `models.ecoli.analysis.multigen` analyses, or
by default the ACTIVE plots listed in that package's `__init__.py`.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import
from __future__ import division

from fireworks import explicit_serialize

from wholecell.fireworks.firetasks.analysisBase import AnalysisBase
import models.ecoli.analysis.multigen


@explicit_serialize
class AnalysisMultiGenTask(AnalysisBase):

	_fw_name = "AnalysisMultiGenTask"
	required_params = [
		"input_seed_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]
	optional_params = [
		"plots_to_run",  # absent or empty => run all active analysis plots
		"output_filename_prefix",
		"cpus",
		]
	MODULE_PATH = 'models.ecoli.analysis.multigen'
	ACTIVE_MODULES = models.ecoli.analysis.multigen.ACTIVE

	def plotter_args(self, module_filename):
		self["metadata"]["analysis_type"] = "multigen"

		return (
			self["input_seed_directory"],
			self["output_plots_directory"],
			self['output_filename_prefix'] + module_filename[:-3],
			self["input_sim_data"],
			self["input_validation_data"],
			self["metadata"],
			)
