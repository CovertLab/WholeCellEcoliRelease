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

def loadSimulation(simDir, outDir, timePoint):
	sim = cPickle.load(
		open(os.path.join(simDir, 'simulation.cPickle'), 'rb'),
		)

	with tables.openFile(os.path.join(simDir, 'state.hdf')) as h5file:
		for state in sim.states:
			state.pytablesLoad(h5file, timePoint)

		time = sim.getState('Time')

		time.value = timePoint

	return sim

if __name__ == '__main__':
	loadSimulation(
		os.path.join('out', 'to_load'),
		os.path.join('out', 'loaded'),
		50
		)