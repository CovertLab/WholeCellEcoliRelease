#!/usr/bin/env python

"""
TfBinding

Bind transcription factors to DNA

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/14/16
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units
from itertools import izip

class TfBinding(wholecell.processes.process.Process):
	""" TfBinding """

	_name = "TfBinding"

	# Constructor
	def __init__(self):
		super(TfBinding, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TfBinding, self).initialize(sim, sim_data)

		# Get IDs of transcription factors
		self.tf_ids = sim_data.process.transcription_regulation.tf_ids
		self.n_TF = len(self.tf_ids)

		# Build dict that maps TFs to transcription units they regulate
		delta_prob = sim_data.process.transcription_regulation.delta_prob
		self.TF_to_TU_idx = {}

		for i, tf in enumerate(self.tf_ids):
			self.TF_to_TU_idx[tf] = delta_prob['deltaI'][
				delta_prob['deltaJ'] == i]

		# Get total counts of transcription units
		self.n_TU = delta_prob['shape'][0]

		# Get constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		# Create dictionaries and method
		self.pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF
		self.tfToTfType = sim_data.process.transcription_regulation.tfToTfType

		# Get DNA polymerase elongation rate (used to mask out promoters that
		# are expected to be replicated in the current timestep)
		self.dnaPolyElngRate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s)))

		# Build views
		self.promoters = self.uniqueMoleculesView("promoter")
		self.active_replisomes = self.uniqueMoleculesView("active_replisome")
		self.active_tf_view = {}
		self.inactive_tf_view = {}

		for tf in self.tf_ids:
			self.active_tf_view[tf] = self.bulkMoleculeView(tf + "[c]")

			if self.tfToTfType[tf] == "1CS":
				if tf == sim_data.process.transcription_regulation.activeToBound[tf]:
					self.inactive_tf_view[tf] = self.bulkMoleculeView(
						sim_data.process.equilibrium.getUnbound(tf + "[c]"))
				else:
					self.inactive_tf_view[tf] = self.bulkMoleculeView(
						sim_data.process.transcription_regulation.activeToBound[tf] + "[c]")
			elif self.tfToTfType[tf] == "2CS":
				self.inactive_tf_view[tf] = self.bulkMoleculeView(
					sim_data.process.two_component_system.activeToInactiveTF[tf + "[c]"])

		# Build array of active TF masses
		bulk_molecule_ids = sim_data.internal_state.bulkMolecules.bulkData["id"]
		tf_indexes = [np.where(bulk_molecule_ids == tf_id + "[c]")[0][0]
			for tf_id in self.tf_ids]
		self.active_tf_masses = (sim_data.internal_state.bulkMolecules.bulkData[
			"mass"][tf_indexes]/self.nAvogadro).asNumber(units.fg)


	def calculateRequest(self):
		# Request all counts of active transcription factors
		for view in self.active_tf_view.itervalues():
			view.requestAll()

		# Request edit access to promoter molecules
		self.promoters.request_access(self.EDIT_ACCESS)


	def evolveState(self):
		# Get attributes of all promoters
		TU_index, coordinates_promoters, domain_index_promoters, bound_TF = self.promoters.attrs(
			"TU_index", "coordinates", "domain_index", "bound_TF")

		# If there are active replisomes, construct mask for promoters that are
		# expected to be replicated in the current timestep. Transcription
		# factors should not bind to these promoters in this timestep.
		# TODO (ggsun): This assumes that replisomes elongate at maximum rates.
		# 	Ideally this should be done in the reconciler.
		collision_mask = np.zeros_like(coordinates_promoters, dtype=np.bool)

		if self.active_replisomes.total_counts()[0] > 0:
			domain_index_replisome, right_replichore, coordinates_replisome = self.active_replisomes.attrs(
				"domain_index", "right_replichore", "coordinates")

			elongation_length = np.ceil(self.dnaPolyElngRate*self.timeStepSec())

			for domain_index, rr, coord in izip(domain_index_replisome,
					right_replichore, coordinates_replisome):
				if rr:
					coordinates_mask = np.logical_and(
						coordinates_promoters >= coord,
						coordinates_promoters <= coord + elongation_length)
				else:
					coordinates_mask = np.logical_and(
						coordinates_promoters <= coord,
						coordinates_promoters >= coord - elongation_length)

				mask = np.logical_and(
					domain_index_promoters == domain_index, coordinates_mask)
				collision_mask[mask] = True

		# Calculate number of bound TFs for each TF prior to changes
		n_bound_TF = bound_TF[~collision_mask, :].sum(axis=0)

		# Initialize new bound_TF array
		bound_TF_new = np.zeros_like(bound_TF, dtype=np.bool)
		bound_TF_new[collision_mask, :] = bound_TF[collision_mask, :]

		# Create vectors for storing values
		pPromotersBound = np.zeros(self.n_TF, dtype=np.float64)
		nPromotersBound = np.zeros(self.n_TF, dtype=np.float64)
		nActualBound = np.zeros(self.n_TF, dtype=np.float64)
		n_bound_TF_per_TU = np.zeros((self.n_TU, self.n_TF), dtype=np.int16)

		for tf_idx, tf_id in enumerate(self.tf_ids):
			# Get counts of transcription factors
			active_tf_counts = self.active_tf_view[tf_id].count()
			bound_tf_counts = n_bound_TF[tf_idx]

			# If there are no active transcription factors to work with,
			# continue to the next transcription factor
			if active_tf_counts + bound_tf_counts == 0:
				continue

			# Free all DNA-bound transcription factors into free active
			# transcription factors
			self.active_tf_view[tf_id].countInc(bound_tf_counts)
			active_tf_counts += bound_tf_counts

			# Compute probability of binding the promoter
			if self.tfToTfType[tf_id] == "0CS":
				pPromoterBound = 1.
			else:
				inactive_tf_counts = self.inactive_tf_view[tf_id].total_counts()
				pPromoterBound = self.pPromoterBoundTF(
					active_tf_counts, inactive_tf_counts)

			# Determine the number of available promoter sites
			available_promoters = np.logical_and(
				np.isin(TU_index, self.TF_to_TU_idx[tf_id]),
				~collision_mask)
			n_available_promoters = available_promoters.sum()

			# Calculate the number of promoters that should be bound
			n_to_bind = int(stochasticRound(
				self.randomState, n_available_promoters*pPromoterBound))

			if n_to_bind > 0:
				# Determine randomly which DNA targets to bind based on which of
				# the following is more limiting:
				# number of promoter sites to bind, or number of active
				# transcription factors
				bound_locs = np.zeros(n_available_promoters, dtype=np.bool)
				bound_locs[
					self.randomState.choice(
						n_available_promoters,
						size=np.min((n_to_bind, self.active_tf_view[tf_id].count())),
						replace=False)
					] = True

				# Update count of free transcription factors
				self.active_tf_view[tf_id].countDec(bound_locs.sum())

				# Update bound_TF array
				bound_TF_new[available_promoters, tf_idx] = bound_locs

			n_bound_TF_per_TU[:, tf_idx] = np.bincount(
				TU_index[bound_TF_new[:, tf_idx]],
				minlength=self.n_TU)

			# Record values
			pPromotersBound[tf_idx] = pPromoterBound
			nPromotersBound[tf_idx] = n_to_bind
			nActualBound[tf_idx] = bound_locs.sum()

		delta_TF = bound_TF_new.astype(np.int8) - bound_TF.astype(np.int8)
		mass_diffs = delta_TF.dot(self.active_tf_masses)

		# Reset bound_TF attribute of promoters
		self.promoters.attrIs(bound_TF=bound_TF_new)

		# Add mass_diffs array to promoter submass
		self.promoters.add_submass_by_array(mass_diffs)

		# Write values to listeners
		self.writeToListener("RnaSynthProb", "pPromoterBound", pPromotersBound)
		self.writeToListener("RnaSynthProb", "nPromoterBound", nPromotersBound)
		self.writeToListener("RnaSynthProb", "nActualBound", nActualBound)
		self.writeToListener(
			"RnaSynthProb", "n_bound_TF_per_TU", n_bound_TF_per_TU)
