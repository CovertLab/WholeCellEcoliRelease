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
from wholecell.analysis.plotting_tools import COLORS_LARGE, COLORS_SMALL

NUMERICAL_ZERO = 1e-8

DISABLED = True

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
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
		timeStepSec = enzymeKineticsdata.readColumn("timeStepSec")[1] - enzymeKineticsdata.readColumn("timeStepSec")[0]

		enzymeKineticsdata.close()

		# This script asssumes certain units - throw an error if they are not the same as in the metabolism.py process
		assert (counts_units[2:] == "[mmol]")
		assert (volume_units[2:] == "[L]")
		assert (mass_units[2:] == "[g]")

		COUNTS_UNITS = units.mmol
		VOLUME_UNITS = units.L
		MASS_UNITS = units.g


		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Load constants
		nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / COUNTS_UNITS)
		cellDensity = sim_data.constants.cellDensity.asNumber(MASS_UNITS/VOLUME_UNITS)

		metabolitePoolIDs = sim_data.process.metabolism.metabolitePoolIDs
		targetConcentrations = sim_data.process.metabolism.metabolitePoolConcentrations.asNumber(COUNTS_UNITS/VOLUME_UNITS)

		# Load enzyme kinetic rate information
		reactionRateInfo = sim_data.process.metabolism.reactionRateInfo
		enzymesWithKineticInfo = sim_data.process.metabolism.enzymesWithKineticInfo["enzymes"]
		constraintIDs = sim_data.process.metabolism.constraintIDs
		constraintToReactionDict = sim_data.process.metabolism.constraintToReactionDict

		extIDs = sim_data.process.metabolism.externalExchangeMolecules
		extMoleculeMasses = sim_data.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)

		moleculeMasses = dict(zip(
			extIDs,
			sim_data.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)
			))

		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		# Load the biomass function
		biomassObjective = sim_data.process.metabolism.biomassFunction

		# Set up FBA solver
		fba = FluxBalanceAnalysis(
			sim_data.process.metabolism.reactionStoich,
			sim_data.process.metabolism.externalExchangeMolecules,
			biomassObjective,
			objectiveType = "standard",
			# internalExchangedMolecules = metabolitePoolIDs,
			reversibleReactions = sim_data.process.metabolism.reversibleReactions,
			moleculeMasses = moleculeMasses,
			)
		
		# Find external molecules levels
		externalMoleculeIDs = fba.externalMoleculeIDs()

		coefficient = initDryMass / initCellMass * sim_data.constants.cellDensity * (timeStepSec * units.s)

		externalMoleculeLevels = sim_data.process.metabolism.exchangeConstraints(
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