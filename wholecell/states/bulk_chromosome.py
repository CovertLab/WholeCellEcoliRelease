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

	def initialize(self, sim, sim_data):
		super(BulkChromosome, self).initialize(sim, sim_data)

		# Load constants
		self._moleculeIDs = sim_data.state.bulkChromosome.bulkData['id']
		self._compartmentIDs = sim_data.state.compartments['compartmentAbbreviation']

		self._moleculeMass = sim_data.state.bulkChromosome.bulkData['mass'].asNumber(units.fg / units.mol) / sim_data.constants.nAvogadro.asNumber(1 / units.mol)

		#self._compIndexes = {compartmentKey:(sim_data.state.bulkChromosome.bulkData['compartment'] == compartmentKey) for compartmentKey in sim_data.state.compartments['compartmentAbbreviation']}

		# Create the container for molecule counts
		self.container = BulkObjectsContainer(self._moleculeIDs)
		
		# TODO: restore this behavior or replace it with something bettter

		self._isRequestAbsolute = np.zeros(self._nProcesses, np.bool)


class BulkChromosomesView(BulkMoleculesView):
	_stateID = 'BulkChromosome'


class BulkChromosomeView(BulkMoleculeView):
	_stateID = 'BulkChromosome'
