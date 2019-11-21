"""
test_bulk_molecules.py

@author: Jerry Morrison
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 4/23/2018
"""

from __future__ import absolute_import
from __future__ import division

import unittest

import numpy as np
import numpy.testing as npt

from wholecell.states.internal_state import InternalState
from wholecell.processes.process import Process
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.states.bulk_molecules import BulkMoleculeView, BulkMoleculesView


OBJECT_NAMES = ('ATP', 'glucose', 'glycine')
OBJECT_COUNTS = [100, 20, 10]


class Test_BulkMolecules(unittest.TestCase):

	def setUp(self):
		self.state = InternalState()
		self.state._countsAllocatedFinal = np.array([OBJECT_COUNTS] * 3)
		self.process = Process()
		self.query = set(OBJECT_NAMES)

		self.container = BulkObjectsContainer(OBJECT_NAMES)
		self.container.countsIs(OBJECT_COUNTS)
		self.state.container = self.container

		self.moleculesView = BulkMoleculesView(self.state, self.process,
			self.query)
		self.moleculeView = BulkMoleculeView(self.state, self.process,
			OBJECT_NAMES[0])

	def test_counts(self):
		moleculesView = self.moleculesView
		npt.assert_equal(moleculesView.counts(), [[OBJECT_COUNTS]] * 3)

		incCounts = [10, 20.5, 40.5]
		newCounts = [110, 40, 50]
		moleculesView.countsInc(incCounts)
		npt.assert_equal(moleculesView.counts(), [[newCounts]] * 3)

		decCounts = [1.5, 2, 3.5]
		newCounts = [109, 38, 47]
		moleculesView.countsDec(decCounts)
		npt.assert_equal(moleculesView.counts(), [[newCounts]] * 3)

		moleculeView = self.moleculeView
		npt.assert_equal(moleculeView.count(), [newCounts])

		newCounts = [209, 138, 147]
		moleculeView.countInc(100.9)
		npt.assert_equal(moleculeView.count(), [newCounts])

		newCounts = [206, 135, 144]
		moleculeView.countDec(3.99)
		npt.assert_equal(moleculeView.count(), [newCounts])
