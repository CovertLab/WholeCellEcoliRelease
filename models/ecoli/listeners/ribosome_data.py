"""
RibosomeData

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/14
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RibosomeData(wholecell.listeners.listener.Listener):
	""" RibosomeData """

	_name = 'RibosomeData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RibosomeData, self).initialize(sim, sim_data)

		self.monomerIds = sim_data.process.translation.monomerData['id'].tolist()
		self.nMonomers = len(self.monomerIds)


	# Allocate memory
	def allocate(self):
		super(RibosomeData, self).allocate()

		# Attributes broadcast by the PolypeptideElongation process
		self.aaCountInSequence = np.zeros(21, np.int64)
		self.aaCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.actualElongationHist = np.zeros(22, np.int64)
		self.elongationsNonTerminatingHist = np.zeros(22, np.int64)
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0
		self.effectiveElongationRate = 0.
		self.rrn16S_produced = 0
		self.rrn23S_produced = 0
		self.rrn5S_produced = 0
		self.rrn16S_init_prob = 0.
		self.rrn23S_init_prob = 0.
		self.rrn5S_init_prob = 0.
		self.total_rna_init = 0
		self.processElongationRate = 0.
		self.translationSupply = np.zeros(21, np.float64)
		self.numTrpATerminated = 0.
		self.probTranslationPerTranscript = np.zeros(self.nMonomers, np.float64)


	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'probTranslationPerTranscript': 'monomerIds'}

		tableWriter.writeAttributes(
			monomerIds = self.monomerIds,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			aaCountInSequence = self.aaCountInSequence,
			aaCounts = self.aaCounts,
			actualElongations = self.actualElongations,
			actualElongationHist = self.actualElongationHist,
			elongationsNonTerminatingHist = self.elongationsNonTerminatingHist,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			terminationLoss = self.terminationLoss,
			effectiveElongationRate = self.effectiveElongationRate,
			rrn16S_produced = self.rrn16S_produced,
			rrn23S_produced = self.rrn23S_produced,
			rrn5S_produced = self.rrn5S_produced,
			rrn16S_init_prob = self.rrn16S_init_prob,
			rrn23S_init_prob = self.rrn23S_init_prob,
			rrn5S_init_prob = self.rrn5S_init_prob,
			total_rna_init = self.total_rna_init,
			processElongationRate = self.processElongationRate,
			translationSupply = self.translationSupply,
			numTrpATerminated = self.numTrpATerminated,
			probTranslationPerTranscript = self.probTranslationPerTranscript,
			)
