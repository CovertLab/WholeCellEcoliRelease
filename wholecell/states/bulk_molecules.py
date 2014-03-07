#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

# TODO: break this into multiple files, it's becoming unbearably long

import re

import numpy as np
import tables

import wholecell.states.state
import wholecell.states.partition
import wholecell.utils.bulk_objects_container


class BulkMolecules(wholecell.states.state.State):
	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'BulkMolecules',
			'name':'Bulk Molecules',
			'dynamics':[],
			'units':{}
			}

		self.time = None
		self.partitionClass = BulkMoleculesPartition

		self._container = None

		super(BulkMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(BulkMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# self._container = wholecell.utils.bulk_objects_container.BulkObjectsContainer(...)

		# TODO


	def allocate(self):
		# TODO
		pass

	
	def calcInitialConditions(self):
		# TODO
		pass


	def partition(self):
		# TODO
		pass


	def countView(self, names):
		return self._container.countsView(names)


	def countView(self, name):
		return self._container.countView(name)

	# TODO: more container method wrapping


	def pytablesCreate(self, h5file, expectedRows):
		# TODO
		pass


	def pytablesAppend(self, h5file):
		# TODO
		pass


	def pytablesLoad(self, h5file, timePoint):
		# TODO
		pass


class BulkMoleculesPartition(wholecell.states.partition.Partition):
	def __init__(self, *args, **kwargs):
		self._container = None
		self._isReqAbs = None

		self._indexMapping = None

		super(UniqueMoleculesPartition, self).__init__(*args, **kwargs)


	def initialize(self, moleculeNames, isReqAbs = False):
		self._container = wholecell.utils.bulk_objects_container.BulkObjectsContainer(moleculeNames)
		self._isReqAbs = isReqAbs

		self._indexMapping = np.array(
			self._state._container._getIndexes(moleculeNames)
			)


	def setRequest(self, target):
		self._process.requestBulkMolecules()

		target[self._indexMapping] = self._container._counts


	def counts(self, names = None):
		return self._container.counts(names)


	def countsIs(self, values, names = None):
		self._container.countsIs(values, names)


	def countsInc(self, values, names = None): 
		self._container.countsInc(values, names)


	def countsDec(self, values, names = None):
		self._container.countsDec(values, names)


	def countsView(self, names = None):
		return self._container.countsView(names)


	def count(self, name):
		return self._container.count(name)


	def countIs(self, value, name):
		self._container.countIs(value, name)


	def countInc(self, value, name):
		self._container.countInc(value, name)


	def countDec(self, value, name):
		self._container.countDec(value, name)


	def countView(self, name):
		return self._container.countView(name)
