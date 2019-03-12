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
		self.tfs = sim_data.process.transcription_regulation.tf_ids
		self.n_tfs = len(self.tfs)

		# Build dict that maps TFs to transcription units they regulate
		delta_prob = sim_data.process.transcription_regulation.delta_prob
		self.tf_to_trs_unit_idx = {}

		for i, tf in enumerate(self.tfs):
			self.tf_to_trs_unit_idx[tf] = delta_prob['deltaI'][
				delta_prob['deltaJ'] == i]

		# Get total counts of transcription units
		self.n_trs_units = delta_prob['shape'][0]

		# Get constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		# Create dictionaries and method
		self.pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF
		self.tfToTfType = sim_data.process.transcription_regulation.tfToTfType

		# Get DNA polymerase elongation rate (used to mask out promoters that
		# is expected to be replicated in the current timestep)
		self.dnaPolyElngRate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s))
			)

		# Build views
		self.promoters = self.uniqueMoleculesView("promoter")
		self.active_replisomes = self.uniqueMoleculesView("active_replisome")
		self.active_tf_view = {}
		self.inactive_tf_view = {}

		for tf in self.tfs:
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
			for tf_id in self.tfs]
		self.active_tf_masses = (sim_data.internal_state.bulkMolecules.bulkData[
			"mass"][tf_indexes]/self.nAvogadro).asNumber(units.fg)


	def calculateRequest(self):
		# Request all counts of active transcription factors
		for view in self.active_tf_view.itervalues():
			view.requestAll()


	def evolveState(self):
		# Get attributes of all promoters
		promoters = self.promoters.molecules_read_and_edit()
		trs_unit_index, coordinates_promoters, domain_index_promoters, bound_tfs = promoters.attrs(
			"trs_unit_index", "coordinates", "domain_index", "bound_tfs"
			)

		# Get attributes of replisomes
		replisomes = self.active_replisomes.molecules_read_only()

		# If there are active replisomes, construct mask for promoters that are
		# expected to be replicated in the current timestep. Transcription
		# factors should not bind to these promoters in this timestep.
		# TODO (ggsun): This assumes that replisomes elongate at maximum rates.
		# 	Ideally this should be done in the reconciler.
		collision_mask = np.zeros_like(coordinates_promoters, dtype=np.bool)

		if len(replisomes) > 0:
			domain_index_replisome, right_replichore, coordinates_replisome = replisomes.attrs(
				"domain_index", "right_replichore", "coordinates",
				)

			elongation_length = np.ceil(self.dnaPolyElngRate*self.timeStepSec())

			for domain_index, rr, coord in izip(
					domain_index_replisome, right_replichore, coordinates_replisome):
				if rr:
					collision_mask[np.logical_and(
						domain_index_promoters == domain_index,
						np.logical_and(
							coordinates_promoters >= coord,
							coordinates_promoters <= coord + elongation_length
							)
						)] = True
				else:
					collision_mask[np.logical_and(
						domain_index_promoters == domain_index,
						np.logical_and(
							coordinates_promoters <= coord,
							coordinates_promoters >= coord - elongation_length
							)
						)] = True

		# Calculate number of bound TFs for each TF prior to changes
		n_bound_tfs = bound_tfs[~collision_mask, :].sum(axis=0)

		# Initialize new bound_tfs array
		bound_tfs_new = np.zeros_like(bound_tfs, dtype=np.bool)
		bound_tfs_new[collision_mask, :] = bound_tfs[collision_mask, :]

		# Create vectors for storing values
		pPromotersBound = np.zeros(self.n_tfs, dtype=np.float64)
		nPromotersBound = np.zeros(self.n_tfs, dtype=np.float64)
		nActualBound = np.zeros(self.n_tfs, dtype=np.float64)
		n_bound_tfs_per_trs_unit = np.zeros(
			(self.n_trs_units, self.n_tfs), dtype=np.int16)

		for tf_idx, tf_id in enumerate(self.tfs):
			# Get counts of transcription factors
			active_tf_counts = self.active_tf_view[tf_id].count()
			bound_tf_counts = n_bound_tfs[tf_idx]

			# If there are no active transcription factors to work with,
			# continue to the next transcription factor
			if active_tf_counts + bound_tf_counts == 0:
				continue

			# Free all DNA-bound transcription factors into free active
			# transcription factors
			self.active_tf_view[tf_id].countInc(bound_tf_counts)

			# Compute probability of binding the promoter
			if self.tfToTfType[tf_id] == "0CS":
				pPromoterBound = 1.
			else:
				inactive_tf_counts = self.inactive_tf_view[tf_id].total_counts()
				pPromoterBound = self.pPromoterBoundTF(
					active_tf_counts, inactive_tf_counts)

			# Determine the number of available promoter sites
			available_promoters = np.logical_and(
				np.isin(trs_unit_index, self.tf_to_trs_unit_idx[tf_id]),
				~collision_mask
				)
			n_available_promoters = available_promoters.sum()

			# Calculate the number of promoters that should be bound
			n_to_bind = int(stochasticRound(
				self.randomState, n_available_promoters*pPromoterBound)
				)

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

				# Update bound_tfs array
				bound_tfs_new[available_promoters, tf_idx] = bound_locs

			n_bound_tfs_per_trs_unit[:, tf_idx] = np.bincount(
				trs_unit_index[bound_tfs_new[:, tf_idx]],
				minlength=self.n_trs_units
				)

			# Record values
			pPromotersBound[tf_idx] = pPromoterBound
			nPromotersBound[tf_idx] = n_to_bind
			nActualBound[tf_idx] = bound_locs.sum()

		delta_tfs = bound_tfs_new.astype(np.int8) - bound_tfs.astype(np.int8)
		mass_diffs = delta_tfs.dot(self.active_tf_masses)

		# Reset bound_tfs attribute of promoters
		promoters.attrIs(bound_tfs=bound_tfs_new)

		# Add mass_diffs array to promoter submass
		promoters.add_submass_by_array(mass_diffs)

		# Write values to listeners
		self.writeToListener("RnaSynthProb", "pPromoterBound", pPromotersBound)
		self.writeToListener("RnaSynthProb", "nPromoterBound", nPromotersBound)
		self.writeToListener("RnaSynthProb", "nActualBound", nActualBound)
		self.writeToListener(
			"RnaSynthProb", "n_bound_tfs_per_trs_unit", n_bound_tfs_per_trs_unit
			)
