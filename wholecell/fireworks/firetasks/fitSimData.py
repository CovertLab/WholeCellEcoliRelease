from __future__ import absolute_import
from __future__ import division

import cPickle
import time
import os
import shutil
import sys

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from reconstruction.ecoli.fit_sim_data_2 import fitSimData_2

@explicit_serialize
class FitSimDataTask(FireTaskBase):

	_fw_name = "FitSimDataTask"
	required_params = ["fit_level", "input_data", "output_data"]
	optional_params = [
		"sim_out_dir",
		"disable_ribosome_capacity_fitting",
		"disable_rnapoly_capacity_fitting",
		"variable_elongation_transcription",
		"variable_elongation_translation",
		"rnapoly_activity_fitting",
		"mrna_half_life_fitting",
		"max_rnap_activity",
		"flat_elongation_transcription",
		"flat_elongation_translation",
		"adjust_rna_and_protein_parameters",
		"adjust_rnase_expression",
		"disable_measured_protein_deg",
		"alternate_mass_fraction_protein",
		"alternate_mass_fraction_rna",
		"alternate_mass_fraction_mrna",
		"alternate_r_protein_degradation",
		"alternate_rna_seq",
		"alternate_rna_half_life",
		"alternate_translation_efficiency",
		"alternate_ribosome_activity",
		"disable_rnap_fraction_increase",
		"disable_ribosome_activity_fix",
		"save_cell_specs",
		"cell_specs_file",
		"write_translation_efficiencies"
		]

	def run_task(self, fw_spec):

		print "%s: Creating/Fitting sim_data (Level %d)" % (time.ctime(), self["fit_level"])

		if self["fit_level"] == 1:
			if self["cached"]:
				try:
					shutil.copyfile(self["cached_data"], self["output_data"])
					mod_time = time.ctime(os.path.getctime(self["cached_data"]))
					print "Copied sim data from cache (modified %s)" % (mod_time,)
					return
				except Exception as exc:
					print ("Warning: could not copy cached sim data due to"
						   " exception (%s), running fitter") % (exc,)

			if self["cpus"] > 1:
				print ("Warning: running fitter in parallel with %i processes -"
					   " ensure there are enough cpus_per_task allocated" % (self["cpus"],))

			with open(self["input_data"], "rb") as f:
				raw_data = cPickle.load(f)

			options = dict(
				cpus=self["cpus"],
				debug=self["debug"],
				disable_ribosome_capacity_fitting=self['disable_ribosome_capacity_fitting'],
				disable_rnapoly_capacity_fitting=self['disable_rnapoly_capacity_fitting'],
				flat_elongation_transcription=not self['variable_elongation_transcription'],
				flat_elongation_translation=not self['variable_elongation_translation'],
				rnapoly_activity_fitting=self['rnapoly_activity_fitting'],
				mrna_half_life_fitting=self['mrna_half_life_fitting'],
				max_rnap_activity=self['max_rnap_activity'],
				adjust_rna_and_protein_parameters=self['adjust_rna_and_protein_parameters'],
				adjust_rnase_expression=self['adjust_rnase_expression'],
				disable_measured_protein_deg=self['disable_measured_protein_deg'],
				alternate_mass_fraction_protein=self['alternate_mass_fraction_protein'],
				alternate_mass_fraction_rna=self['alternate_mass_fraction_rna'],
				alternate_mass_fraction_mrna=self['alternate_mass_fraction_mrna'],
				alternate_r_protein_degradation=self['alternate_r_protein_degradation'],
				alternate_rna_seq=self['alternate_rna_seq'],
				alternate_rna_half_life=self['alternate_rna_half_life'],
				alternate_translation_efficiency=self['alternate_translation_efficiency'],
				alternate_ribosome_activity=self['alternate_ribosome_activity'],
				disable_rnap_fraction_increase=self['disable_rnap_fraction_increase'],
				disable_ribosome_activity_fix=self['disable_ribosome_activity_fix'],
				write_translation_efficiencies=self['write_translation_efficiencies'])

			sim_data, cell_specs = fitSimData_1(
				raw_data,
				options)

			sys.setrecursionlimit(4000) #limit found manually
			with open(self["output_data"], "wb") as f:
				cPickle.dump(sim_data, f, protocol = cPickle.HIGHEST_PROTOCOL)

			if self["save_cell_specs"]:
				with open(self["cell_specs_file"], "wb") as f:
					cPickle.dump(cell_specs, f, protocol = cPickle.HIGHEST_PROTOCOL)

		# TODO: Get rid of this if not used
		if self["fit_level"] == 2:
			with open(self["input_data"], "rb") as f:
				sim_data = cPickle.load(f)

			fitSimData_2(sim_data, self["sim_out_dir"])

			with open(self["output_data"], "wb") as f:
				cPickle.dump(sim_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
