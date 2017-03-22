#!/usr/bin/env python

"""
ChromosomeFormation

Takes partial chromosomes and concatonates them into full chromosomes
This is a modeling approximation that cleans up the way replication elongation
is modeled.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/20/2015
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

class ChromosomeFormation(wholecell.processes.process.Process):
	""" ChromosomeFormation """

	_name = "ChromosomeFormation"

	def __init__(self):
		super(ChromosomeFormation, self).__init__()

	def initialize(self, sim, sim_data):
		super(ChromosomeFormation, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		self.D_period = sim_data.growthRateParameters.d_period.asNumber(units.s)

		# Create views on partial and full chromosomes
		# NOTE: Parial chromosomes are modeling artifact. They are single strands of half a chromosome.
		# this was done to simplify the chromosome replication polymerization process and make it analogous
		# to what was done in transcription and translation elongation. This process removes the artifact.
		self.partialChromosomes = self.bulkMoleculesView(sim_data.moleculeGroups.partialChromosome)
		self.fullChromosome = self.bulkMoleculeView(sim_data.moleculeGroups.fullChromosome[0])

		# Placeholder for cell division data
		self.fullChromosomeUnique = self.uniqueMoleculesView("fullChromosome")

	def calculateRequest(self):
		self.partialChromosomes.requestAll()

	def evolveState(self):
		partialChromosomes = self.partialChromosomes.counts()

		# If >1 of each partial chromosome exists, turn them into a standard full chromosome
		if partialChromosomes.min():
			fullUniqueChrom = self.fullChromosomeUnique.moleculesNew("fullChromosome", partialChromosomes.min())
			# Log the time that the C period finishes and set the corresponding division time to be D period time later in seconds
			fullUniqueChrom.attrIs(division_time = [self.time() + self.D_period] * partialChromosomes.min())

		# Decrement and increment counts
		self.fullChromosome.countInc(partialChromosomes.min())
		self.partialChromosomes.countsDec(partialChromosomes.min())