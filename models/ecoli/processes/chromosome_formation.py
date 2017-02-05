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
	def initialize(self, sim, sim_data):
		super(ChromosomeFormation, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)

		self.partialChromosomes = self.bulkMoleculesView(sim_data.moleculeGroups.partialChromosome)
		self.fullChromosome = self.bulkMoleculeView(sim_data.moleculeGroups.fullChromosome[0])

		self.fullChromosomeUnique = self.uniqueMoleculesView("fullChromosome")

	def calculateRequest(self):
		self.partialChromosomes.requestAll()

	def evolveState(self):
		partialChromosomes = self.partialChromosomes.counts()
		if partialChromosomes.min():
			fullUniqueChrom = self.fullChromosomeUnique.moleculesNew("fullChromosome", partialChromosomes.min())
			fullUniqueChrom.attrIs(division_time = [self.time() + 20. * 60.] * partialChromosomes.min())

			print "grep_marker replication termination - time: {}".format(self.time())
			print "grep_marker cell division occurs - relative time: {}".format(self.time() - self._sim.initialTime())
			print "grep_marker replication termination - time set for division: {}".format(self.time() + 20. * 60.)
			print "grep_marker replication termination - cell mass: {}".format(self.readFromListener("Mass", "cellMass"))
			print "grep_marker replication termination - partial chromosome counts: {}".format(partialChromosomes)

		self.fullChromosome.countInc(partialChromosomes.min())
		self.partialChromosomes.countsDec(partialChromosomes.min())