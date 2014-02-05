"""
SimulationParser.py

Functions for parsing and running simulations from a JSON file.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/4/2014
"""

import json
import sys
import os
import cPickle

import wholecell.sim.logger.Disk
import wholecell.sim.logger.Shell
import wholecell.sim.Simulation
import wholecell.kb.KnowledgeBase
import wholecell.util.Fitter

# The default values for the parsed JSON.  Also serves as a template for 
# writing JSON files.
DEFAULT_JSON = '''
{
	"fitSimulation":true,
	"useCachedKB":true,
	"cacheKB":true,
	"useShellLogger":true,
	"useDiskLogger":false,
	"diskLoggerPath":null,
	"simOptions":{}
}
'''

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')


# TODO: add support for running many simulations, possibly in parallel, based 
# on one JSON file, one cached simulation, and an integer for the number of 
# seeds
def runSimulationFromJson(filePath):
	options = json.loads(DEFAULT_JSON)
	options.update(json.load(open(filePath, 'r')))

	if not options['useCachedKB'] or not os.path.exists(KB_PATH):
		kb = wholecell.kb.KnowledgeBase.KnowledgeBase(
			dataFileDir = "data/parsed", seqFileName = "data/raw/sequence.txt"
			)

		if options['cacheKB']:
			cPickle.dump(kb, open(KB_PATH, "wb"),
				protocol = cPickle.HIGHEST_PROTOCOL)

	else:
		kb = cPickle.load(open(KB_PATH, "rb"))

	sim = wholecell.sim.Simulation.Simulation()
	sim.initialize(kb)
	sim.setOptions(options['simOptions'])

	if options['fitSimulation']:
		wholecell.util.Fitter.Fitter.FitSimulation(sim, kb)

	# Instantiate loggers
	if options['useShellLogger']:
		sim.loggerAdd(wholecell.sim.logger.Shell.Shell())

	if options['useDiskLogger']:
		if outDir is None:
			raise Exception('No output directory provided.')

		sim.loggerAdd(
			wholecell.sim.logger.Disk.Disk(outDir = options['diskLoggerPath'])
			)

	# Run simulation
	sim.run()


if __name__ == '__main__':
	runSimulationFromJson(sys.argv[1])


