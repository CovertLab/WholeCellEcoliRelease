'''
loadSimulation

Loads a simulation.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/31/14
'''

import os
import cPickle

import wholecell.loggers.shell
import wholecell.sim.simulation as wcSim

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

if __name__ == '__main__':
	TIME_STEP = 50

	sim = wcSim.Simulation.loadSimulation(
		cPickle.load(open(KB_PATH, 'rb')),
		os.path.join('out', 'to_load', 'state.hdf'),
		TIME_STEP
		)

	sim.setOptions({'lengthSec':100})

	sim.loggerAdd(wholecell.loggers.shell.Shell())

	sim.run()
