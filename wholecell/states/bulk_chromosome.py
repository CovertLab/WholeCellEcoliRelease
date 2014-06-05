#!/usr/bin/env python

"""
BulkChromosome.py

State which represents for a class of molecules or chromosome sites by their bulk copy numbers.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

from wholecell.containers.bulk_objects_container import BulkObjectsContainer

class BulkChromosome(wholecell.states.state.bulk_molecules.BulkMolecules):
	_name = 'BulkChromosome'

	def initialize(self, sim, kb):
		super(BulkChromosome, self).initialize(sim, kb)

		# Load constants
		self._moleculeIDs = moleculeIds(kb)
		self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkChromosome['mass'].to('fg / mol').magnitude / kb.nAvogadro.to('1 / mole').magnitude

		self._typeIdxs = {
			'genes'	:	kb.bulkChromosome['isGene']
			}

		self._compIndexes = {
			compartmentKey:(kb.bulkChromosome['compartment'] == compartmentKey)
			for compartmentKey in kb.compartments['compartmentAbbreviation']
			}

		# Create the container for molecule counts
		self.container = bulkObjectsContainer(kb)
		
		# TODO: restore this behavior or replace it with something bettter

		self._isRequestAbsolute = np.zeros(self._nProcesses, np.bool)

def moleculeIds(kb):
	return kb.bulkChromosome['moleculeId']

def bulkObjectsContainer(kb, dtype = np.int64):
	return BulkObjectsContainer(moleculeIds(kb), dtype)

class BulkChromosomeViewBase(wholecell.states.bulk_moelcules.BulkMoleculesViewBase):
	_stateID = 'BulkChromosome'

class BulkChromosomeView(BulkChromosomeViewBase):
	pass

class BulkChromosomeView(BulkChromosomeViewBase):
	pass