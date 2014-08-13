#!/usr/bin/env python

"""
justSimulation.py

Runs a simulation.  Called from justSimulation.sh

"""

from __future__ import division

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
import os
import json

def main():

	simOpts = getSimOptsFromEnvVars()

	simOpts["kbLocation"] = os.path.join("fixtures", "kb", "KnowledgeBase_Fit.cPickle")

	sim = EcoliSimulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	main()
