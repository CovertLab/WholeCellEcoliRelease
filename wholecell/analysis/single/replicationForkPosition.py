#!/usr/bin/env python
"""
Plot amino acid counts

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "ReplicationForkPosition.hdf"))
	fork0position = h.root.ReplicationForkPosition.col('fork0position')
	fork1position = h.root.ReplicationForkPosition.col('fork1position')
	h.close()

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))
	table = h.root.Mass
	time = np.array([x["time"] for x in table.iterrows()])
	h.close()

	plt.figure(figsize = (8.5, 11))

	plt.plot(time / 60., fork0position, linewidth = 2)
	plt.plot(time / 60., fork1position, linewidth = 2)
	plt.xlabel("Time (min)")
	plt.ylabel("Replication fork position (nt)")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])