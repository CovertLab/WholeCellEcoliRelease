"""
Run the plots for a given list of `models.ecoli.analysis.variant` analyses, or
by default the ACTIVE plots listed in that package's `__init__.py`.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import
from __future__ import division

from fireworks import explicit_serialize

from wholecell.fireworks.firetasks.analysisBase import AnalysisBase
import models.ecoli.analysis.variant


@explicit_serialize
class AnalysisVariantTask(AnalysisBase):

	_fw_name = "AnalysisVariantTask"
	required_params = [
		"input_directory",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]
	optional_params = [
		"plots_to_run",  # absent or empty => run all active analysis plots
		"output_filename_prefix",
		"cpus",
		]
	MODULE_PATH = 'models.ecoli.analysis.variant'
	ACTIVE_MODULES = models.ecoli.analysis.variant.ACTIVE

	def plotter_args(self, module_filename):
		return (
			self["input_directory"],
			self["output_plots_directory"],
			self['output_filename_prefix'] + module_filename[:-3],
			'unused_sim_data_filename',
			self['input_validation_data'],
			self["metadata"],
			)
