"""
ATPhydrolyzedUsageListener

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/4/15
"""

from __future__ import absolute_import, division, print_function

import wholecell.listeners.listener

class ATPhydrolyzedUsageListener(wholecell.listeners.listener.Listener):
	""" ATPhydrolyzedUsageListener """

	_name = 'ATPhydrolyzedUsageListener'

	def __init__(self, *args, **kwargs):
		super(ATPhydrolyzedUsageListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ATPhydrolyzedUsageListener, self).initialize(sim, sim_data)

		self.atpsHydrolyzed = 0.

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes( # TODO: reconsider attribute names
			atpsHydrolyzed = self.countUnits,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			atpsHydrolyzed = self.atpsHydrolyzed,
			)
