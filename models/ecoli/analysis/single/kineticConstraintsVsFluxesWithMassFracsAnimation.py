#!/usr/bin/env python
"""
Plot reaction max rate over course of the simulation.

@date: Created 7/02/2015
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import matplotlib.animation as animation

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

COLORS = [
	[colorValue/255. for colorValue in color]
	for color in COLORS_256
	]

DISABLED = False

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	if DISABLED:
		print "Currently disabled because it is slow."
		return

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	
	enzymeKineticsArray = enzymeKineticsdata.readColumn("reactionConstraints")
	overconstraintMultiples = enzymeKineticsdata.readColumn("overconstraintMultiples")

	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()


	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))

	reactionRatesUnconstrained = fbaData.readColumn('reactionFluxes')

	rateEstimatesArray = enzymeKineticsArray.copy()

	unconstrainedFluxes = fbaData.readColumn('reactionFluxes')
	fluxNames = fbaData.readAttribute('reactionIDs')
	simulationSteps = fbaData.readColumn('simulationStep')
	fbaData.close()

	mass = TableReader(os.path.join(simOutDir, "Mass"))

	protein = mass.readColumn("proteinMass")
	tRna = mass.readColumn("tRnaMass")
	rRna = mass.readColumn("rRnaMass")
	mRna = mass.readColumn("mRnaMass")
	dna = mass.readColumn("dnaMass")
	smallMolecules = mass.readColumn("smallMoleculeMass")

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	masses = np.vstack([
		protein/protein[0],
		rRna/rRna[0],
		tRna/tRna[0],
		mRna/mRna[0],
		dna/dna[0],
		smallMolecules/smallMolecules[0],
		]).T

	massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Mol.s"]

	# Only look at fluxes which had a noninfinite estimate at at least one point
	rateEstimates =  rateEstimatesArray[1:]
	rateEstimates += 1.0
	rateEstimates[np.where(rateEstimates == np.inf)] = 0

	fluxesWithEstimates = unconstrainedFluxes[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	fluxNamesEstimates = np.array([fluxNames[x] for x in np.where(np.sum(rateEstimates, axis=0) > 0)[0]])

	rateEstimates = rateEstimates[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	rateEstimates -= 1.0
	rateEstimates = np.vstack((np.zeros(rateEstimates.shape[1]),rateEstimates))

	relativeRates = fluxesWithEstimates / rateEstimates

	# Only plot rates which are overconstrained
	relativeRates[np.where(relativeRates < 1)] = 0
	relativeRates[np.isnan(relativeRates)] = 0

	num_x = relativeRates.shape[1]

	# First set up the figure, the axis, and the plot element we want to animate
	fig, axarr = plt.subplots(3, figsize=(8.5,11))

	massCompAx = axarr[0]
	massCompAx.set_color_cycle(COLORS)
	lines = []
	for idx, mass in enumerate(np.transpose(masses)):
		line = massCompAx.plot(t / 60., mass, linewidth = 2, label=massLabels[idx])
		lines.append(line[0])
	massCompAx.set_xlim(0, 1.05*(t[-1] / 60.))
	plt.xlabel("Time (min)")
	plt.ylabel("Mass (normalized by t = 0 min)")
	plt.title("Biomass components")

	massCompAx.legend(loc = "best", framealpha=.5, fontsize=8)

	ratesAx =  axarr[1]
	plt.title("Log Reaction Fluxes and Kinetic Rates")

	rects1_fluxes = ratesAx.bar(range(1, num_x+1), np.log10(fluxesWithEstimates[0,:] + 1), .7, color="b", alpha=.7, label="Fluxes")
	rects1_rates = ratesAx.bar(range(1, num_x+1), np.log10(rateEstimates[0,:] + 1), .5, color="r", alpha=.5, label="Kinetic Rates")

	y_max = np.max((np.log10(rateEstimates + 1).max(), np.log10(rateEstimates + 1).max()))

	ratesAx.set_xlim(0, num_x)
	ratesAx.set_ylim(0, y_max)
	
	ratesAx.set_title("Log Reaction Fluxes and Kinetic Rates")
	ratesAx.set_ylabel("Log10 Reaction Rate ({counts_units}/{volume_units}.{time_units})".format(counts_units=COUNTS_UNITS.strUnit(), volume_units=VOLUME_UNITS.strUnit(), time_units=TIME_UNITS.strUnit()))

	ratesAx.legend(["Fluxes", "Kinetic Rates"],framealpha=.5, fontsize=12)

	barax = axarr[2]

	barax.set_xlim(0, num_x)
	barax.set_ylim(0, 1.1*relativeRates.max())
	
	plt.title("Kinetic Rates Divided By Reaction Fluxes")
	plt.xlabel("Reaction")
	plt.ylabel("Fold Difference")

	rects2 = barax.bar(range(1, num_x+1), relativeRates[0,:],  align='center', color="#800000", alpha=.7, label="Fold Difference")

	# Trim the flux names to fit on the plot
	fluxNamesEstimatesTrimmed = [x[:15] for x in fluxNamesEstimates]
	plt.xticks(range(1, num_x+1), fluxNamesEstimatesTrimmed, rotation='vertical', fontsize=7)

	# initialization function: plot the background of each frame
	def init():
		for idx, patch in enumerate(rects2.patches):
			patch.set_height(relativeRates[0,idx])

		for idx, patch in enumerate(rects1_fluxes.patches):
			patch.set_height(np.log10(fluxesWithEstimates[0,idx] + 1))
		for idx, patch in enumerate(rects1_rates.patches):
			patch.set_height(np.log10(rateEstimates[0,idx] + 1))
		return rects2, rects1_fluxes, rects1_rates

	# animation function.  This is called sequentially
	def animate(i):
		if i % 10 == 0:
			print "Frame %i" % (i)

		# Animate lines
		for idx, line in enumerate(lines):
			x = t[:i] / 60.
			y = masses[:i,idx]
			line.set_data(x,y)

		# Fluxes and limits
		for idx, patch in enumerate(rects1_fluxes.patches):
			patch.set_height(np.log10(fluxesWithEstimates[i,idx] + 1))
		for idx, patch in enumerate(rects1_rates.patches):
			patch.set_height(np.log10(rateEstimates[i,idx] + 1))
		
		# Ratios
		for idx, patch in enumerate(rects2.patches):
			patch.set_height(relativeRates[i, idx])

		return lines, rects1_fluxes, rects1_rates, rects2

	# call the animator.  blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(fig, animate, init_func=init,
		frames=100, interval=20, blit=True)

	# anim = animation.FuncAnimation(fig, animate, init_func=init,
	# 	frames=relativeRates.shape[0], interval=20, blit=True)

	# save the animation as an mp4.  This requires ffmpeg or mencoder to be
	# installed.  The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5.  You may need to adjust this for
	# your system: for more information, see
	# http://matplotlib.sourceforge.net/api/animation_api.html
	anim.save(os.path.join(plotOutDir, plotOutFileName) + '.mp4', fps=50, extra_args=['-vcodec', 'libx264'])

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__
	
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
