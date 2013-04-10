#!/usr/bin/env python

"""
runSimulation
Runs and logs whole-cell simulation

Example:
>>> from runSimulation import *
>>> runSimulation()

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/10/2013
"""

import numpy
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

def runSimulation(simOpts = None, diskOpts = None, kbOpts = None):
	## Import classes
	import wholecell.sim.logger.Disk
	import wholecell.sim.logger.Shell
	import wholecell.sim.Simulation

	simOpts = simOpts or {}
	diskOpts = diskOpts or {}
	kbOpts = kbOpts or {}

	## Instantiate knowledge base
	kb = wholecell.kb.KnowledgeBase.KnowledgeBase(**kbOpts)

	## Instantiate loggers
	diskLogger = wholecell.sim.logger.Disk.Disk(**diskOpts)
	shellLogger = wholecell.sim.logger.Shell.Shell()

	loggers = [
		diskLogger,
		shellLogger
	]

	## Run simulation
	sim = wholecell.sim.Simulation.Simulation(kb)
	sim.setOptions(simOpts)
	sim.run(loggers)

	## Plot Result
	# Load data
	time = numpy.transpose(wholecell.sim.logger.Disk.Disk.load(diskLogger.outDir, "Time", "value") / 3600, (0, 2, 1)).reshape(-1)
	mass = numpy.sum(wholecell.sim.logger.Disk.Disk.load(diskLogger.outDir, "Mass", "cell"), axis = 1).reshape(-1)

	plt.plot(time, mass)
	plt.xlabel("Time (h)")
	plt.ylabel("Mass (fg)")
	plt.plot(time[[0, -1]], mass[0] * numpy.ones(2), color = 0.5 * numpy.ones(3), linestyle=":")
	plt.plot(time[[0, -1]], mass[-1] * numpy.ones(2), color = 0.5 * numpy.ones(3), linestyle=":")
	plt.xlim(time[[0, -1]])
	plt.ylim(numpy.array([numpy.min(mass), numpy.max(mass)]) + 0.05 * numpy.ptp(mass) * numpy.array([-1, 1]))
	plt.show()
	plt.close()