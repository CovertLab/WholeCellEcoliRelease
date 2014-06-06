#!/usr/bin/env python

"""
BulkChromosome.py

State which represents for a class of molecules or chromosome sites by their bulk copy numbers.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division
import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.states.bulk_molecules import BulkMolecules, BulkMoleculesView, BulkMoleculeView

class BulkChromosome(BulkMolecules):
	_name = 'BulkChromosome'

	def initialize(self, sim, kb):
		super(BulkChromosome, self).initialize(sim, kb)

		# Load constants
		self._moleculeIDs = kb.bulkChromosome['moleculeId']
		self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkChromosome['mass'].to('fg / mol').magnitude / kb.nAvogadro.to('1 / mole').magnitude

		self._typeIdxs = {
			'genes':kb.bulkChromosome['isGene'],
			# Eventually, when there are "bound" states, these will be non-empty
			# (the mass refactoring may happen before that, however)
			'metabolites':np.zeros(0, np.int64),
			'rnas':np.zeros(0, np.int64),
			'proteins':np.zeros(0, np.int64),
			'water':np.zeros(0, np.int64),
			}

		self._compIndexes = {
			compartmentKey:(kb.bulkChromosome['compartment'] == compartmentKey)
			for compartmentKey in kb.compartments['compartmentAbbreviation']
			}

		# Create the container for molecule counts
		self.container = BulkObjectsContainer(self._moleculeIDs)
		
		# TODO: restore this behavior or replace it with something bettter

		self._isRequestAbsolute = np.zeros(self._nProcesses, np.bool)


class BulkChromosomesView(BulkMoleculesView):
	_stateID = 'BulkChromosome'


class BulkChromosomeView(BulkMoleculeView):
	_stateID = 'BulkChromosome'
