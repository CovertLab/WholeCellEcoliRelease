#!/usr/bin/env python

"""
ChromosomeFormation

Takes partial chromosomes and concatenates them into full chromosomes
This is a modeling approximation that cleans up the way replication elongation
is modeled.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/20/2015
"""

from __future__ import division

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
		# NOTE: Partial chromosomes are a modeling artifact. They are single
		# strands of half a chromosome. This was done to simplify the
		# chromosome replication polymerization process and make it analogous
		# to what was done in transcription and translation elongation. This
		# process removes the artifact.
		self.partialChromosomes = self.bulkMoleculesView(
			sim_data.moleculeGroups.partialChromosome)
		self.fullChromosome = self.bulkMoleculeView(
			sim_data.moleculeIds.fullChromosome)

		# Placeholder for cell division data
		self.fullChromosomeUnique = self.uniqueMoleculesView("fullChromosome")

	def calculateRequest(self):
		self.partialChromosomes.requestAll()

	def evolveState(self):
		# Get counts of each of the four partial chromosomes
		# Note: only the fully elongated partial chromosomes with half the
		# length of the full chromosome are counted as one.
		partialChromosomeCounts = self.partialChromosomes.counts()

		# If all four partial chromosomes exist, convert to full chromosome
		if partialChromosomeCounts.min() > 0:
			fullUniqueChrom = self.fullChromosomeUnique.moleculesNew(
				"fullChromosome", partialChromosomeCounts.min())

			# Log the current time (end of C period) and set the corresponding
			# division time to be D period time later in seconds
			fullUniqueChrom.attrIs(
				division_time = [self.time() + self.D_period]
							    * partialChromosomeCounts.min()
				)

		# Decrement and increment counts
		self.fullChromosome.countInc(partialChromosomeCounts.min())
		self.partialChromosomes.countsDec(partialChromosomeCounts.min())
