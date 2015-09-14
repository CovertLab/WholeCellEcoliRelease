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

		objective = dict(zip(
			metabolitePoolIDs,
			targetConcentrations
			))

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

		# Set up FBA solver
		fba = FluxBalanceAnalysis(
			kb.process.metabolism.reactionStoich.copy(), # TODO: copy in class
			kb.process.metabolism.externalExchangeMolecules,
			objective,
			objectiveType = "pools",
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
		desiredConcentrations = [objective[name] for name in metaboliteNames]

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

			oneConstraintOutputFilename = plotOutDir + '/' + plotOutFileName + '-single_constraints.tsv'
			twoConstraintOutputFilename = plotOutDir + '/' + plotOutFileName + '-double_constraints.tsv'
			allbutOneConstraintOutputFilename = plotOutDir + '/' + plotOutFileName + '-single_dropout_constraints.tsv'
			reactionIDsPositiveNegativeErrorsFilename = plotOutDir + '/' + plotOutFileName + '-reactionIDsPositiveNegativeErrors.tsv'
			greedyConstraintSelectionFilename = plotOutDir + '/' + plotOutFileName + '-greedyConstraintSelection'

			if findSingleConstraints:
				determineOneConstraintErrors(oneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			if findDoubleConstraints:
				determineTwoConstraintErrors(twoConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find errors for all constraints and set of all constraints except one
			if findAllButOneConstraint:
				determineAllButOneConstraints(allbutOneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find errors for a random assortment of N reactions set to zero flux (From ALL reactions, not just kinetic reactions)
			if findRandomRxnZero:
				plotNumConstraintsVersusError(reactionIDsPositiveNegativeErrorsFilename, fba, 35, 20, .5, plotOutDir, plotOutFileName, reactionIDs, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint=timepoint)

			## Greedy subset selection
			if findMaxConstraintsGreedy:
				# Constraint subset with which to start
				alreadyConstrained = np.zeros(len(constraints))
				# Select additional constraints
				selectMaxConstraintsSubset(greedyConstraintSelectionFilename, plotOutDir, plotOutFileName, fba, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint=timepoint)

			if findMinMaxConstraintsGreedy:
				# Constraint subset with which to start
				alreadyConstrained = np.zeros(len(constraints))
				# Select additional constraints
				selectMinMaxConstraintsSubset(greedyConstraintSelectionFilename, plotOutDir, plotOutFileName, metadata, fba, .9, 1.1, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint=timepoint)

		timepoints = np.arange(allConstraints.shape[1])[1:]

		if findTimeCourseError:

			# Constraint subset with which to start
			alreadyConstrained = np.zeros(len(constraints))
			# alreadyConstrained = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
			alreadyConstrainedLength = alreadyConstrained.shape[0]

			np.random.seed(3)
			alreadyConstrained = (np.random.random(alreadyConstrainedLength) < .1)

			# alreadyConstrained.fill(1)
			# alreadyConstrained[(alreadyConstrainedLength/4):(alreadyConstrainedLength/3)] = 1
			

			# alreadyConstrained[63] = 1

			# alreadyConstrained[39:44] = 0
			# # alreadyConstrained[:(alreadyConstrainedLength/4)] = 1
			# alreadyConstrained[(alreadyConstrainedLength/2):] = 1

			plotSingleConditionTimeCourse(plotOutDir, plotOutFileName, metadata, fba, .9, 1.1, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInitRaw, desiredConcentrations, metaboliteNames, countsToMolar, timepoints, colorMap, linestyleMap)



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


def evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar):
		
	metaboliteConcentrations = metaboliteCountsInit * countsToMolar

	fba.internalMoleculeLevelsIs(
			metaboliteConcentrations
			)

	deltaMetabolites = fba.outputMoleculeLevelsChange() / countsToMolar

	metaboliteCountsFinal = metaboliteCountsInit + deltaMetabolites

	individualErrors = (metaboliteCountsFinal * countsToMolar) - desiredConcentrations

	totalError = np.sum(np.absolute(individualErrors))

	return individualErrors, totalError


def determineOneConstraintErrors(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "ConstraintID", "ReactionID", "Rate Limit", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndividualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

		for index, constraint in enumerate(constraintIDs):

			# Set the constraint
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)

			# Evaluate the error with this constraint
			individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			individualRelativeErrors = individualErrors - baseIndividualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in individualRelativeErrors]

			output_file.write("\n" + "\t".join(row))

			# Reset the rate limit
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)



def determineTwoConstraintErrors(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "ConstraintID_1", "ConstraintID_2", "ReactionID_1", "ReactionID_2", "Rate Limit 1", "Rate Limit 2", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndividualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

		for index1, constraint1 in enumerate(constraintIDs):
			# Set the first constraint
			fba.maxReactionFluxIs(constraintToReactionDict[constraint1], constraints[index1], raiseForReversible = False)

			for index2, constraint2 in enumerate(constraintIDs):

				# Don't constrain the same reaction twice
				if constraintToReactionDict[constraint1] == constraintToReactionDict[constraint2]:
					continue

				# Don't compute duplicate combined constraints
				if index2<index1:
					continue

				# Set the second constraint
				fba.maxReactionFluxIs(constraintToReactionDict[constraint2], constraints[index2], raiseForReversible = False)

				# Evaluate the error with these constraints
				individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

				# Find the error relative to the unconstrained solution
				individualRelativeErrors = individualErrors - baseIndividualErrors
				totalRelativeError = totalError - baseTotalError

				row = [str(timepoint), constraint1, constraint2, str(constraintToReactionDict[constraint1]), str(constraintToReactionDict[constraint2]), str(constraints[index1]), str(constraints[index2]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in individualRelativeErrors]

				output_file.write("\n" + "\t".join(row))

				# Reset the rate limit
				fba.maxReactionFluxIs(constraintToReactionDict[constraint2], np.inf, raiseForReversible = False)
			# Reset the rate limit
			fba.maxReactionFluxIs(constraintToReactionDict[constraint1], np.inf, raiseForReversible = False)


def determineAllButOneConstraints(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "Removed ConstraintID", "ReactionID", " Removed Rate Limit", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndividualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

		# Apply all constraints
		for index, constraint in enumerate(constraintIDs):
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)


		# Evaluate the error with all constraints
		individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

		# Find the error relative to the unconstrained solution
		individualRelativeErrors = individualErrors - baseIndividualErrors
		totalRelativeError = totalError - baseTotalError

		row = [str(timepoint), "All constraints", "All Reactions", "No constraints removed", str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in individualRelativeErrors]

		output_file.write("\n" + "\t".join(row))


		# Remove each in turn and evaluate
		for index, constraint in enumerate(constraintIDs):
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)

			# Evaluate the error with these constraints
			individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			individualRelativeErrors = individualErrors - baseIndividualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in individualRelativeErrors]

			output_file.write("\n" + "\t".join(row))

			# Put this constraint back in place
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)


def determineErrorNRandomFluxesPercentConstraint(fba, numZeroFluxes, constraintPercentage, reactionIDs, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):
	"""
		determineErrorNRandomFluxesPercentConstraint

		Sets the flux of numZeroFluxes randomly choosen reactions to their current value multiplied by constraintPercentage.

		Negative constraintPercentage is treated as zero constraint. Negative fluxes are set to zero.
	"""

	# Choose a random vector of numZeroFluxes from reactionIDs
	randomFluxes = []
	for x in xrange(1,int(numZeroFluxes)):
		randomFluxes.append(np.random.choice(range(len(reactionIDs))))

	# Set the fluxes of these numZeroFluxes reactionIDs to constraintPercentage of their current value
	for reaction in randomFluxes:
		if constraintPercentage > NUMERICAL_ZERO:
			currentFlux = fba.reactionFlux(reactionIDs[reaction])
			if currentFlux < 0:
				fba.maxReactionFluxIs(reactionIDs[reaction], 0, raiseForReversible = False)
				# The flux is negative, so set the minimum instead
				# fba.minReactionFluxIs(reactionIDs[reaction], constraintPercentage*currentFlux, raiseForReversible = False)
			else:
				# Set the max flux as a percentage of the current flux
				fba.maxReactionFluxIs(reactionIDs[reaction], constraintPercentage*currentFlux, raiseForReversible = False)
		else:
			# If the constraint is essentially zero, just set the constraint to zero
			fba.maxReactionFluxIs(reactionIDs[reaction], 0, raiseForReversible = False)

	# Evaluate the error with these constraints
	individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

	# Remove these constraints again
	for reaction in randomFluxes:
		fba.maxReactionFluxIs(reactionIDs[reaction], np.inf, raiseForReversible = False)

	# Record which reactions were constrained
	constrainedReactions = []
	for reaction in randomFluxes:
		constrainedReactions.append(reactionIDs[reaction])

	return individualErrors, totalError, constrainedReactions



def plotNumConstraintsVersusError(outputFileName, fba, samplesPerPoint, numPointsToSample, constraintPercentage, plotOutDir, plotOutFileName, reactionIDs, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	constraintNum = []
	totalRelativeErrors = []

	# Sample from 0 to full length, inclusive, total of numPointsToSample points
	pointsToSample = np.append(np.arange(0,len(reactionIDs),len(reactionIDs)/(numPointsToSample-1)),len(reactionIDs))

	# Find the unconstrained deviation from the target metabolite concentrations
	baseIndividualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)


	lowErrorReactions = set()
	lowErrorCombinations = set()

	for point in pointsToSample:
		for sample in range(samplesPerPoint):

			individualErrors, totalError, constrainedReactions = determineErrorNRandomFluxesPercentConstraint(fba, point, constraintPercentage, reactionIDs, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find the error relative to the unconstrained solution
			individualRelativeErrors = individualErrors - baseIndividualErrors
			totalRelativeError = totalError - baseTotalError

			# Record this point
			constraintNum.append(int(point))
			totalRelativeErrors.append(totalRelativeError)

			# Record any constraints which produce zero error, as well as all constituents of those combinations
			if totalRelativeError <= .1:
				lowErrorCombinations.add(tuple([totalRelativeError, totalError, tuple(constrainedReactions)]))
				for constrainedReaction in constrainedReactions:
					lowErrorReactions.add(constrainedReaction)

			# For zero and all constraints, this number won't change (sampling 0 from N or N from N only gives one answer)
			if point == 0 or point == len(reactionIDs):
				# Record this value for all the values planned for sampling
				for i in range(samplesPerPoint-1):
					constraintNum.append(int(point))
					totalRelativeErrors.append(totalRelativeError)
				# Don't sample this point any more
				break

	lowErrorReactionsArray = [x for x in lowErrorReactions]

	with open(outputFileName, 'w') as output_file:
		output_file.write("Number of Constraints" + "\t" + "Low Error Combinations" + "\t" + "Total Relative Error" + "\t" + "Total Error" + "\t" + "Zero Error Reactions" + "\n")
		for index, entry in enumerate(lowErrorCombinations):
			# Unpack the error and the tuple of combinations
			entry = list(entry)
			totalRelativeError = entry[0]
			totalError = entry[1]
			combination = entry[2]

			# Record all low error combinations in a spreadsheet ALSO record all reactions which appear in a low error combination
			if index < len(lowErrorReactionsArray):
				output_file.write(str(len(combination)) + "\t" + str(combination) + "\t" + str(totalRelativeError) + "\t" + str(totalError) + "\t" + lowErrorReactionsArray[index])
			else:
				output_file.write(str(len(combination)) + "\t" + str(combination) + "\t" + str(totalRelativeError) + "\t" + str(totalError))
			output_file.write('\n')


	# Find the mean and standard deviation of results at each point
	sampleVector = []
	sampleMeans = []
	sampleMeansAbs = []
	sampleStds = []
	constraintNumAdjusted = []
	for index, point in enumerate(totalRelativeErrors):
		sampleVector.append(totalRelativeErrors[index])
		if len(sampleVector) == samplesPerPoint:
			sampleMeans.append(np.average(sampleVector))
			sampleMeansAbs.append(np.average(np.abs(sampleVector)))
			sampleStds.append(np.std(sampleVector))
			sampleVector = []
			constraintNumAdjusted.append(constraintNum[index])


	zoomedTotalRelativeErrors = []
	zoomedConstraintNum = []

	inRangeOfFullConstraint = (np.array(totalRelativeErrors) < 1.1*totalRelativeErrors[-1])
	nearZeros = (np.array(totalRelativeErrors) < .1)
	for index, relativeError in enumerate(totalRelativeErrors):
		if nearZeros[index]:
			zoomedTotalRelativeErrors.append(relativeError)
			zoomedConstraintNum.append(constraintNum[index])

	# Plot the results
	plt.figure(figsize = (8.5, 11))

	plt.subplot(3,1,1)
	plt.title("Random Constraints (Relative Flux = %s), Timepoint= %d Seconds" % (timepoint, constraintPercentage))
	plt.scatter(constraintNum, totalRelativeErrors, c='b')
	plt.ylabel("jFBA Error Relative to No Constraint")

	plt.subplot(3,1,2)
	plt.scatter(zoomedConstraintNum, zoomedTotalRelativeErrors, c='b')
	plt.ylabel("jFBA Error Relative to No Constraint")

	plt.subplot(3,1,3)
	plt.bar(constraintNumAdjusted, sampleMeansAbs, yerr=sampleStds)
	plt.ylabel("Mean Error (abs value)")

	# plt.subplot(3,1,3)
	# plt.bar(constraintNumAdjusted, sampleStds)
	# plt.ylabel("Std Dev")


	plt.xlabel("Number of Reactions Constrained")


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")



def determineNextConstraint(fba, alreadyConstrained, metaboliteCountsInit, minConstraintMultiple, maxConstraintMultiple, constraints, constraintIDs, constraintToReactionDict, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):
	"""
		determineNextConstraint

		Starting from an fba object, selects the next constraint which minimizes the jFBA error.

		Inputs: alreadyConstrained - a boolean vector of constraints already added

		Returns: the alreadyConstrained vector with one more constraint added, and the name of the constraint to be added.
	"""
		
	# Find the starting jFBA error
	baseIndividualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

	# Keeps track of the smallest value so far
	currentMin = np.inf

	# For every constraint not yet used, evaluate the error
	for index, constraint in enumerate(constraintIDs):
		# Skip constraints already applied
		if alreadyConstrained[index]:
			continue

		# Set the constraint
		fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*maxConstraintMultiple, raiseForReversible = False)
		fba.minReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*minConstraintMultiple, raiseForReversible = False)

		# Evaluate the error with this constraint
		individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

		# Find the error relative to the unconstrained solution
		totalRelativeError = totalError - baseTotalError

		if totalError < currentMin:
			currentMin = totalError
			currentBestConstraintIndex = index

		# Reset the rate limit
		fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)

	alreadyConstrained[currentBestConstraintIndex] = True
	additionalError = currentMin - baseTotalError

	return alreadyConstrained, currentBestConstraintIndex, additionalError


def selectMaxConstraintsSubset(outputFileName, plotOutDir, plotOutFileName, fba, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):
	"""
		selectConstraintsSubset

		Selects constraints to greedily minimize jFBA error.

		Removes any flux constraints currently present, and then applies all those marked 'True' in alreadyConstrained.
	"""
	
	# Constrain those reactions marked in alreadyConstrained, remove any others
	for index, constraint in enumerate(constraintIDs):
		if alreadyConstrained[index]:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)
		else:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)

	# Find the current deviation from the target metabolite concentrations
	notUsed, currentTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

	# Number of constraints already present at start
	startingNumConstraints = np.sum(alreadyConstrained)

	constraintNum = []
	totalRelativeErrors = []
	additionalErrors = []

	# Prepare the output file
	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "Added ConstraintID", "Max Flux Value", "ReactionID", "Total Error", "Additional Error", "Constraints Boolean Array", "Reactions Boolean Array"]
		output_file.write("\t".join(headers))

		# Greedily add constraints one by one
		for point in range(np.amin(len(constraints))):
			alreadyConstrained, newConstraintIndex, additionalError = determineNextConstraint(fba, alreadyConstrained, metaboliteCountsInit, 0, 1, constraints, constraintIDs, constraintToReactionDict, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none')

			newConstraint = constraintIDs[newConstraintIndex]

			# Apply this new constraint
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[newConstraintIndex], raiseForReversible = False)

			# Add the new error to the existing error
			currentTotalError += additionalError

			# Record this point and constraint
			constraintNum.append(int(point))
			totalRelativeErrors.append(currentTotalError)
			additionalErrors.append(additionalError)

			# Write this point to file
			row = [str(timepoint), str(newConstraint), str(constraints[newConstraintIndex]), str(constraintToReactionDict[constraint]), str(currentTotalError), str(additionalError),str([int(x) for x in alreadyConstrained])]
			output_file.write("\n" + "\t".join(row))


	# Plot the results
	plt.figure(figsize = (8.5, 11))

	plt.subplot(2,1,1)
	plt.title("Greedy Constraint Selection, timepoint = %s" % timepoint)
	plt.scatter(constraintNum + startingNumConstraints, totalRelativeErrors, c='b')
	plt.ylabel("jFBA Error")
	plt.xlabel("Number of Reactions Constrained")

	plt.subplot(2,1,2)
	plt.scatter(constraintNum, additionalErrors, c='b')
	plt.ylabel("Additional jFBA Error per New Constraint")
	plt.xlabel("Number of Constraints Added")


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


def selectMinMaxConstraintsSubset(outputFileName, plotOutDir, plotOutFileName, metadata, fba, minConstraintMultiple, maxConstraintMultiple, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInit, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):
	"""
		selectMinMaxConstraintsSubset

		Selects constraints in turn to greedily minimize jFBA error. Applies both a min and a max constraint,
		above and below the predicted value by maxConstraintMultiple and minConstraintMultiple respectively.

		Removes any flux constraints currently present, and then applies all those marked 'True' in alreadyConstrained.
	"""
	
	# Constrain those reactions marked in alreadyConstrained, remove any others
	for index, constraint in enumerate(constraintIDs):
		if alreadyConstrained[index]:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*maxConstraintMultiple, raiseForReversible = False)
			fba.minReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*minConstraintMultiple, raiseForReversible = False)			
		else:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)
			fba.minReactionFluxIs(constraintToReactionDict[constraint], 0, raiseForReversible = False)			

	# Find the current deviation from the target metabolite concentrations
	notUsed, currentTotalError = evaluatePoint(fba, metaboliteCountsInit, desiredConcentrations, countsToMolar)

	# Number of constraints already present at start
	startingNumConstraints = np.sum(alreadyConstrained)

	constraintNum = []
	totalRelativeErrors = []
	additionalErrors = []

	# Prepare the output file
	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "Added ConstraintID", "Flux Value (within %s to %sx of given value)" % (str(minConstraintMultiple), str(maxConstraintMultiple)), "ReactionID", "Total Error", "Additional Error", "Constraints Boolean Array", "Reactions Boolean Array"]
		output_file.write("\t".join(headers))

		# Greedily add constraints one by one
		for point in range(np.amin(len(constraints))):
			alreadyConstrained, newConstraintIndex, additionalError = determineNextConstraint(fba, alreadyConstrained, metaboliteCountsInit, minConstraintMultiple, maxConstraintMultiple, constraints, constraintIDs, constraintToReactionDict, desiredConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none')

			newConstraint = constraintIDs[newConstraintIndex]

			# Apply this new constraint
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[newConstraintIndex]*maxConstraintMultiple, raiseForReversible = False)
			fba.minReactionFluxIs(constraintToReactionDict[constraint], constraints[newConstraintIndex]*minConstraintMultiple, raiseForReversible = False)

			# Add the new error to the existing error
			currentTotalError += additionalError

			# Record this point and constraint
			constraintNum.append(int(point))
			totalRelativeErrors.append(currentTotalError)
			additionalErrors.append(additionalError)

			# Write this point to file
			row = [str(timepoint), str(newConstraint), str(constraints[newConstraintIndex]), str(constraintToReactionDict[constraint]), str(currentTotalError), str(additionalError),str([int(x) for x in alreadyConstrained])]
			output_file.write("\n" + "\t".join(row))


	# Plot the results
	plt.figure(figsize = (8.5, 11))

	plt.subplot(2,1,1)
	plt.title("Greedy Constraint Selection (Min = %s, Max=%s), timepoint = %s" % (str(minConstraintMultiple), str(maxConstraintMultiple), timepoint) )
	plt.scatter(constraintNum + startingNumConstraints, totalRelativeErrors, c='b')
	plt.ylabel("jFBA Error")
	plt.xlabel("Number of Reactions Constrained")

	plt.subplot(2,1,2)
	plt.scatter(constraintNum, additionalErrors, c='b')
	plt.ylabel("Additional jFBA Error per New Constraint")
	plt.xlabel("Number of Constraints Added")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


def plotSingleConditionTimeCourse(plotOutDir, plotOutFileName, metadata, fba, minConstraintMultiple, maxConstraintMultiple, constraints, constraintIDs, constraintToReactionDict, alreadyConstrained, metaboliteCountsInitRaw, desiredConcentrations, metaboliteNames, countsToMolar, timepoints, colorMap, linestyleMap):

	# Constrain those reactions marked in alreadyConstrained, remove any others
	for index, constraint in enumerate(constraintIDs):
		if alreadyConstrained[index]:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*maxConstraintMultiple, raiseForReversible = False)
			fba.minReactionFluxIs(constraintToReactionDict[constraint], constraints[index]*minConstraintMultiple, raiseForReversible = False)			
		else:
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)
			fba.minReactionFluxIs(constraintToReactionDict[constraint], 0, raiseForReversible = False)			


	errorsTimeCourse = np.zeros((timepoints.shape[0],metaboliteCountsInitRaw.shape[0] + 1))
	# Loop through all timepoints, determine the jFBA error for each
	for index, timepoint in enumerate(timepoints):

		# Find the deviation from the target metabolite concentrations
		individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInitRaw[:,timepoint], desiredConcentrations, countsToMolar)

		errorsTimeCourse[index][0] = totalError

		errorsTimeCourse[index][1:] = individualErrors

	# Plot the results
	plt.figure(figsize = (10, 15))

	ax0 = plt.subplot(4,1,1)
	plt.title("jFBA Error Time Course, %i constraints" % (alreadyConstrained.sum()))
	plt.scatter(timepoints, np.transpose(errorsTimeCourse)[0], c='b')
	plt.ylabel("jFBA Total Error")
	plt.xlabel("Time (timeStepSec)")

	ax1 = plt.subplot(4,1,2)
	plt.scatter(timepoints, np.transpose(errorsTimeCourse)[0], c='b')
	plt.ylabel("jFBA Total Error")
	plt.xlabel("Time (timeStepSec)")
	ax1.set_ylim([0,1000])

	lineLabelsTop = []
	plt.subplot(4,1,3)
	plt.hold(True)
	for index, timeCourse in enumerate(np.transpose(errorsTimeCourse)[1:]):
		if np.amax(np.abs(timeCourse)) > .01:
			plt.plot(timeCourse, label=metaboliteNames[index], color=colorMap[metaboliteNames[index]], linestyle=linestyleMap[metaboliteNames[index]])
			lineLabelsTop.append(metaboliteNames[index])
	plt.hold(False)
	plt.ylabel("jFBA individual Errors")
	plt.xlabel("Time (timeStepSec)")
	plt.legend(lineLabelsTop, bbox_to_anchor=(1.05, 3.66), loc=2, borderaxespad=0.)	

	ax = plt.subplot(4,1,4)
	plt.hold(True)
	for index, timeCourse in enumerate(np.transpose(errorsTimeCourse)[1:]):
		if np.amax(np.abs(timeCourse)) > .01:
			plt.plot(timeCourse, label=metaboliteNames[index], color=colorMap[metaboliteNames[index]], linestyle=linestyleMap[metaboliteNames[index]])
	plt.hold(False)
	plt.ylabel("jFBA individual Errors")
	plt.xlabel("Time (timeStepSec)")
	ax.set_ylim([-.1,.02])

	plt.subplots_adjust(left=.12, right=.65, top=0.9, bottom=0.1)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")