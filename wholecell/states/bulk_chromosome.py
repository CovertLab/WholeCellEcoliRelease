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
from wholecell.utils import units

class BulkChromosome(BulkMolecules):
	_name = 'BulkChromosome'

	def initialize(self, sim, kb):
		super(BulkChromosome, self).initialize(sim, kb)

		# Load constants
		self._moleculeIDs = kb.bulkChromosome['moleculeId']
		self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkChromosome['mass'].asNumber(units.fg / units.mol) / kb.nAvogadro.asNumber(1 / units.mol)

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
