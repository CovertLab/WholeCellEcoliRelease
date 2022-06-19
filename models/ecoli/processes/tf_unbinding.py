"""
TfUnbinding

Unbind transcription factors from DNA to allow signaling processes before
binding back to DNA.
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils import units


class TfUnbinding(wholecell.processes.process.Process):
	""" TfUnbinding """

	_name = "TfUnbinding"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super().initialize(sim, sim_data)

		# Get IDs of transcription factors
		self.tf_ids = sim_data.process.transcription_regulation.tf_ids

		# Build views
		self.promoters = self.uniqueMoleculesView("promoter")
		self.active_tf_view = {}

		for tf in self.tf_ids:
			self.active_tf_view[tf] = self.bulkMoleculeView(tf + "[c]")

		# Build array of active TF masses
		bulk_molecule_ids = sim_data.internal_state.bulk_molecules.bulk_data["id"]
		tf_indexes = [np.where(bulk_molecule_ids == tf_id + "[c]")[0][0]
			for tf_id in self.tf_ids]
		self.active_tf_masses = (sim_data.internal_state.bulk_molecules.bulk_data[
			"mass"][tf_indexes] / sim_data.constants.n_avogadro).asNumber(units.fg)

	def calculateRequest(self):
		# Request edit access to promoter molecules
		self.promoters.request_access(self.EDIT_ACCESS)

	def evolveState(self):
		# If there are no promoters, return immediately
		if self.promoters.total_count() == 0:
			return

		# Get attributes of all promoters
		bound_TF = self.promoters.attr("bound_TF")

		# Calculate number of bound TFs for each TF prior to changes
		n_bound_TF = bound_TF.sum(axis=0)

		# Free all DNA-bound transcription factors into free active
		# transcription factors
		for tf_idx, tf_id in enumerate(self.tf_ids):
			self.active_tf_view[tf_id].countInc(n_bound_TF[tf_idx])

		# Reset bound_TF attribute of promoters
		self.promoters.attrIs(bound_TF=np.zeros_like(bound_TF))

		# Add mass_diffs array to promoter submass
		mass_diffs = bound_TF @ -self.active_tf_masses
		self.promoters.add_submass_by_array(mass_diffs)
