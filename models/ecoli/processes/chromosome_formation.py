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

	# Constructor
	def __init__(self):
		super(ChromosomeFormation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ChromosomeFormation, self).initialize(sim, kb)

		# Load constants
		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / units.mol)

		self.partialChromosomes = self.bulkMoleculesView(kb.moleculeGroups.partialChromosome)
		self.fullChromosome = self.bulkMoleculeView(kb.moleculeGroups.fullChromosome[0])

	def calculateRequest(self):
		self.partialChromosomes.requestAll()

	def evolveState(self):
		partialChromosomes = self.partialChromosomes.counts()
		self.fullChromosome.countInc(partialChromosomes.min())
		self.partialChromosomes.countsDec(partialChromosomes.min())