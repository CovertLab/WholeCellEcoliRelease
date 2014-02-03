'''
loadSimulation

Loads a simulation.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/31/14
'''

import os
import cPickle

import tables

import wholecell.sim.logger.Shell
import wholecell.sim.Simulation
import wholecell.util.Fitter

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

def loadSimulation(simDir, outDir, timeStep):
	kb = cPickle.load(open(KB_PATH, 'rb'))

	sim = wholecell.sim.Simulation.Simulation(kb)

	sim.setOptions({'lengthSec':100})

	# wholecell.util.Fitter.Fitter.FitSimulation(sim, kb)

	# sim.calcInitialConditions()

	with tables.openFile(os.path.join(simDir, 'state.hdf')) as h5file:
		sim.pytablesLoad(h5file, timeStep)
		for state in sim.states.values():
			state.pytablesLoad(h5file, timeStep)

		# hack to get the right time
		time = sim.states['Time']
		time.value = timeStep * sim.timeStepSec

	return sim


if __name__ == '__main__':
	TIME_STEP = 50

	sim = loadSimulation(
		os.path.join('out', 'to_load'),
		os.path.join('out', 'loaded'),
		TIME_STEP
		)

	sim.run([wholecell.sim.logger.Shell.Shell()], TIME_STEP)
