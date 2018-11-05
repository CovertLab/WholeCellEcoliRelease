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


class TfBinding(wholecell.processes.process.Process):
	""" TfBinding """

	_name = "TfBinding"

	# Constructor
	def __init__(self):
		super(TfBinding, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TfBinding, self).initialize(sim, sim_data)

		# Get IDs of gene copy numbers and transcription factor bound targets
		bound_tf_names = sim_data.process.transcription_regulation.boundTFColNames
		gene_copy_names = sim_data.process.transcription_regulation.geneCopyNumberColNames
		self.tfs = sorted(set([x.split("__")[-1] for x in bound_tf_names]))

		tf_to_bound_tfs = {}
		self.tf_to_target_index = {}
		for tf in self.tfs:
			tf_to_bound_tfs[tf] = [x for x in bound_tf_names if x.endswith("__" + tf)]
			self.tf_to_target_index[tf] = [
				gene_copy_names.index(x.split("__")[0] + "__alpha")
				for x in tf_to_bound_tfs[tf]
				]

		self.n_tfs = len(self.tfs)

		# Get constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		# Create dictionaries and method
		self.pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF
		self.tfToTfType = sim_data.process.transcription_regulation.tfToTfType

		# Build views
		self.gene_copy_view = self.bulkMoleculesView(gene_copy_names)
		self.bound_tf_views = {}
		self.active_tf_view = {}
		self.inactive_tf_view = {}

		for tf in self.tfs:
			self.bound_tf_views[tf] = self.bulkMoleculesView(tf_to_bound_tfs[tf])
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


	def calculateRequest(self):
		# Get current copy numbers of each gene
		self.gene_copy_numbers = self.gene_copy_view.total()

		# Request all counts of active transcription factors
		for view in self.active_tf_view.itervalues():
			view.requestAll()

		# Request all counts of DNA bound transcription factors
		for view in self.bound_tf_views.itervalues():
			view.requestAll()


	def evolveState(self):
		# Create vectors for storing values
		pPromotersBound = np.zeros(self.n_tfs, np.float64)
		nPromotersBound = np.zeros(self.n_tfs, np.float64)
		nActualBound = np.zeros(self.n_tfs, np.float64)

		for i, tf in enumerate(self.tfs):
			# Get counts of transcription factors
			active_tf_counts = self.active_tf_view[tf].count()
			bound_tf_counts = self.bound_tf_views[tf].counts()
			target_gene_counts = self.gene_copy_numbers[
				self.tf_to_target_index[tf]]

			# Free all DNA-bound transcription factors into free active
			# transcription factors
			self.bound_tf_views[tf].countsIs(0)
			self.active_tf_view[tf].countInc(bound_tf_counts.sum())

			# If there are no active transcription factors, continue to the
			# next transcription factor
			if active_tf_counts == 0:
				continue

			# Compute probability of binding the promoter
			if self.tfToTfType[tf] == "0CS":
				if active_tf_counts + bound_tf_counts.sum() > 0:
					pPromoterBound = 1.
				else:
					pPromoterBound = 0.
			else:
				tfInactiveCounts = self.inactive_tf_view[tf].total()
				pPromoterBound = self.pPromoterBoundTF(
					active_tf_counts, tfInactiveCounts)

			# Determine the number of available promoter sites to bind
			n_to_bind = int(stochasticRound(
				self.randomState, target_gene_counts.sum()*pPromoterBound)
				)

			# If there are no promoter sites to bind, continue to the next
			# transcription factor
			if n_to_bind == 0:
				continue

			# Determine randomly which DNA targets to bind based on which of
			# the following is more limiting:
			# number of promoter sites to bind, or number of active
			# transcription factors
			bound_locs = np.zeros(target_gene_counts.sum(), dtype=np.int)
			bound_locs[
				self.randomState.choice(
					target_gene_counts.sum(),
					size=np.min((n_to_bind, self.active_tf_view[tf].count())),
					replace=False)
				] = 1

			# Update count of free transcription factors
			self.active_tf_view[tf].countDec(bound_locs.sum())

			# Update count of bound transcription factors
			bound_tf_counts_updated = np.ediff1d(
				bound_locs.cumsum()[np.cumsum(target_gene_counts) - 1],
				to_begin=bound_locs[:target_gene_counts[0]].sum())

			self.bound_tf_views[tf].countsIs(bound_tf_counts_updated)

			# Record values
			pPromotersBound[i] = pPromoterBound
			nPromotersBound[i] = pPromoterBound*target_gene_counts.sum()
			nActualBound[i] = bound_locs.sum()

		# Write values to listeners
		self.writeToListener("RnaSynthProb", "pPromoterBound", pPromotersBound)
		self.writeToListener("RnaSynthProb", "nPromoterBound", nPromotersBound)
		self.writeToListener("RnaSynthProb", "nActualBound", nActualBound)
