#!/usr/bin/env python
"""
Runs tests outside the model to determine which kinetic limits can allow survival.

@date: Created 8/18/2015
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
from random import random

from wholecell.utils import units
import cPickle
import ast
import itertools

from wholecell.utils.modular_fba import FluxBalanceAnalysis

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

NUMERICAL_ZERO = 1e-8

DISABLED = True

COLORS_LARGE = ["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]

COLORS_SMALL = ["#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#000000", "#007FFF",
		"#236B8E", "#70DB93", "#B5A642", "#5F9F9F", "#B87333", "#2F4F2F", "#9932CD", "#871F78", "#855E42",
		"#545454", "#8E2323", "#238E23", "#CD7F32", "#527F76",
		"#9F9F5F", "#8E236B", "#2F2F4F", "#CFB53B", "#FF7F00", "#DB70DB",
		"#5959AB", "#8C1717", "#238E68", "#6B4226", "#8E6B23", "#00FF7F",
		"#38B0DE", "#DB9370", "#5C4033", "#4F2F4F", "#CC3299", "#99CC32"]


def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)


	if DISABLED:
		print "Currently disabled because it's slow."
	else:

		# Control which analyses to run
		findSingleConstraints = False
		findDoubleConstraints = False
		findAllButOneConstraint = False
		findRandomRxnZero = False
		findMaxConstraintsGreedy = False
		findMinMaxConstraintsGreedy = False
		findTimeCourseError = True

		# Check that only one plotting script is active.
		if (findRandomRxnZero + findMaxConstraintsGreedy + findMinMaxConstraintsGreedy + findTimeCourseError) > 1:
			print "Can only run one of %s, %s, %s, %s at a time." % ("findRandomRxnZero" + "findMaxConstraintsGreedy" + "findMinMaxConstraintsGreedy" + "findTimeCourseError")
			assert False

		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

		# Read constraints and metabolite concentrations from the listener
		allConstraints = enzymeKineticsdata.readColumn("allConstraintsLimits")
		constraintIDs = enzymeKineticsdata.readAttribute("constraintIDs")
		metaboliteCountsInitRaw = enzymeKineticsdata.readColumn("metaboliteCountsInit")
		modelMetaboliteCountsFinalRaw = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
		enzymeCountsInitRaw = enzymeKineticsdata.readColumn("enzymeCountsInit")
		countsToMolarArray = enzymeKineticsdata.readColumn("countsToMolar")
		reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
		constraintToReactionDict = enzymeKineticsdata.readAttribute("constraintToReactionDict")

		# Read units from listener
		counts_units = enzymeKineticsdata.readColumn("counts_units")[1]
		volume_units = enzymeKineticsdata.readColumn("volume_units")[1]
		mass_units = enzymeKineticsdata.readColumn("mass_units")[1]
		
		# Read time info from the listener
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		timeStepSec = enzymeKineticsdata.readColumn("timeStep")[1] - enzymeKineticsdata.readColumn("timeStep")[0]

		enzymeKineticsdata.close()

		# This script asssumes certain units - throw an error if they are not the same as in the metabolism.py process
		assert (counts_units[2:] == "[mmol]")
		assert (volume_units[2:] == "[L]")
		assert (mass_units[2:] == "[g]")

		COUNTS_UNITS = units.mmol
		VOLUME_UNITS = units.L
		MASS_UNITS = units.g


		# Load data from KB
		kb = cPickle.load(open(kbFile, "rb"))

		# Load constants
		nAvogadro = kb.constants.nAvogadro.asNumber(1 / COUNTS_UNITS)
		cellDensity = kb.constants.cellDensity.asNumber(MASS_UNITS/VOLUME_UNITS)

		metabolitePoolIDs = kb.process.metabolism.metabolitePoolIDs
		targetConcentrations = kb.process.metabolism.metabolitePoolConcentrations.asNumber(COUNTS_UNITS/VOLUME_UNITS)

		# Load enzyme kinetic rate information
		reactionRateInfo = kb.process.metabolism.reactionRateInfo
		enzymesWithKineticInfo = kb.process.metabolism.enzymesWithKineticInfo["enzymes"]
		constraintIDs = kb.process.metabolism.constraintIDs
		constraintToReactionDict = kb.process.metabolism.constraintToReactionDict

		extIDs = kb.process.metabolism.externalExchangeMolecules
		extMoleculeMasses = kb.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)

		moleculeMasses = dict(zip(
			extIDs,
			kb.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)
			))

		initWaterMass = kb.mass.avgCellWaterMassInit
		initDryMass = kb.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		energyCostPerWetMass = kb.constants.darkATP * initDryMass / initCellMass

		# Load the biomass function
		biomassObjective = kb.process.metabolism.biomassFunction

		# Set up FBA solver
		fba = FluxBalanceAnalysis(
			kb.process.metabolism.reactionStoich.copy(),
			kb.process.metabolism.externalExchangeMolecules,
			biomassObjective,
			objectiveType = "standard",
			# internalExchangedMolecules = metabolitePoolIDs,
			reversibleReactions = kb.process.metabolism.reversibleReactions,
			moleculeMasses = moleculeMasses,
			)
		
		# Find external molecules levels
		externalMoleculeIDs = fba.externalMoleculeIDs()

		coefficient = initDryMass / initCellMass * kb.constants.cellDensity * (timeStepSec * units.s)

		externalMoleculeLevels = kb.process.metabolism.exchangeConstraints(
			externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS
			)

		# Set FBA external molecule levels
		fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		# targetConcentrations in same order as fba output (targetConcentrations is ordered differently)
		metaboliteNames = fba.outputMoleculeIDs()

		# Transpose the constraints and metabolites matrices
		allConstraints = np.transpose(allConstraints)
		metaboliteCountsInitRaw = np.transpose(metaboliteCountsInitRaw)
		modelMetaboliteCountsFinalRaw = np.transpose(modelMetaboliteCountsFinalRaw)		

		# Consider only a single timepoint of external molecules levels
		externalMoleculeLevels = externalMoleculeLevels[1]

		# External molecules constraint
		fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		# Build map from metabolites to colors, reusing as few as possible
		colorMap = dict(itertools.izip(metaboliteNames, itertools.cycle(COLORS_SMALL)))
		# Build map from metabolites to linestyles
		linestyles = ['-', '--', '-.',':']
		linestyles_subset = linestyles[:3]
		n_repeats = len(metaboliteNames)//len(linestyles_subset)
		linestylesLong = []
		for linestyle in linestyles_subset:
			linestylesLong.extend([linestyle]*n_repeats)
		linestyleMap = dict(itertools.izip(metaboliteNames, itertools.cycle(linestylesLong)))


		# FBA analysis on finished timepoints
		# timepoints = [100, 1000, 2000, allConstraints.shape[1]-1]
		timepoints = [2000]
		# timepoints = np.arange(allConstraints.shape[1])

		for timepoint in timepoints:
			# Constraints at one arbitrary time point, fairly far into the simulation
			constraints = allConstraints[:,timepoint]

			# Metabolite concentrations
			metaboliteCountsInit = metaboliteCountsInitRaw[:,timepoint]

			# Relationship between molecule counts and concentrations at this time point
			countsToMolar = countsToMolarArray[timepoint]

			metaboliteConcentrations = metaboliteCountsInit * countsToMolar

			fba.internalMoleculeLevelsIs(
				metaboliteConcentrations
				)

			deltaMetabolites = fba.outputMoleculeLevelsChange() / countsToMolar

			metaboliteCountsFinal = metaboliteCountsInit + deltaMetabolites



if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__
	
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])