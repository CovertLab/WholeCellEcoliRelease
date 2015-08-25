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

from wholecell.utils.modular_fba import FluxBalanceAnalysis

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

DISABLED = True

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


		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

		# Read constraints and metabolite concentrations from the listener
		allConstraints = enzymeKineticsdata.readColumn("allConstraintsLimits")
		constraintIDs = enzymeKineticsdata.readAttribute("constraintIDs")
		metaboliteCountsInitRaw = enzymeKineticsdata.readColumn("metaboliteCountsInit")
		countsToMolarArray = enzymeKineticsdata.readColumn("countsToMolar")
		reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
		constraintToReactionDict = enzymeKineticsdata.readAttribute("constraintToReactionDict")

		# Data needed to reconstruct the FBA object in analysis
		reactionStoich = enzymeKineticsdata.readAttribute("reactionStoich")
		externalExchangeMolecules = enzymeKineticsdata.readAttribute("externalExchangeMolecules")
		objective = enzymeKineticsdata.readAttribute("objective")
		reversibleReactions = enzymeKineticsdata.readAttribute("reversibleReactions")
		moleculeMasses = enzymeKineticsdata.readAttribute("moleculeMasses")
		metabolitePoolIDs = enzymeKineticsdata.readAttribute("metabolitePoolIDs")
		targetConcentrations = enzymeKineticsdata.readAttribute("targetConcentrations")
		
		# Read time info from the listener
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		enzymeKineticsdata.close()

		# Transpose the constraints and metabolites matrices
		allConstraints = np.transpose(allConstraints)
		metaboliteCountsInitRaw = np.transpose(metaboliteCountsInitRaw)


		# Set up FBA solver
		fba = FluxBalanceAnalysis(
			reactionStoich, # TODO: copy in class
			externalExchangeMolecules,
			objective,
			objectiveType = "pools",
			reversibleReactions = reversibleReactions,
			moleculeMasses = moleculeMasses,
			)

		# FBA analysis on finished timepoints
		# timepoints = [100, 1000, 2000, allConstraints.shape[1]-1]
		timepoints = [2000]

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
			reactionIDsPositiveNegativeErrors = plotOutDir + '/' + plotOutFileName + '-reactionIDsPositiveNegativeErrors.tsv'

			if findSingleConstraints:
				determineOneConstraintErrors(oneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			if findDoubleConstraints:
				determineTwoConstraintErrors(twoConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find errors for all constraints and set of all constraints except one
			if findAllButOneConstraint:
				determineAllButOneConstraints(allbutOneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find errors for a random assortment of N reactions set to zero flux (From ALL reactions, not just kinetic reactions)
			if findRandomRxnZero:
				plotNumConstraintsVersusError(reactionIDsPositiveNegativeErrors, fba, 30, 15, plotOutDir, plotOutFileName, reactionIDs, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint=timepoint)



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




def evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar):
		
	metaboliteConcentrations = metaboliteCountsInit * countsToMolar

	fba.internalMoleculeLevelsIs(
			metaboliteConcentrations
			)

	deltaMetabolites = fba.outputMoleculeLevelsChange() / countsToMolar

	metaboliteCountsFinal = metaboliteCountsInit + deltaMetabolites

	individualErrors = (metaboliteCountsFinal * countsToMolar) - targetConcentrations

	totalError = np.sum(np.absolute(individualErrors))

	return individualErrors, totalError


def determineOneConstraintErrors(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "ConstraintID", "ReactionID", "Rate Limit", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndiviualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

		for index, constraint in enumerate(constraintIDs):

			# Set the constraint
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)

			# Evaluate the error with this constraint
			individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			indiviualRelativeErrors = individualErrors - baseIndiviualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in indiviualRelativeErrors]

			output_file.write("\n" + "\t".join(row))

			# Reset the rate limit
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)



def determineTwoConstraintErrors(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "ConstraintID_1", "ConstraintID_2", "ReactionID_1", "ReactionID_2", "Rate Limit 1", "Rate Limit 2", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndiviualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

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
				individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

				# Find the error relative to the unconstrained solution
				indiviualRelativeErrors = individualErrors - baseIndiviualErrors
				totalRelativeError = totalError - baseTotalError

				row = [str(timepoint), constraint1, constraint2, str(constraintToReactionDict[constraint1]), str(constraintToReactionDict[constraint2]), str(constraints[index1]), str(constraints[index2]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in indiviualRelativeErrors]

				output_file.write("\n" + "\t".join(row))

				# Reset the rate limit
				fba.maxReactionFluxIs(constraintToReactionDict[constraint2], np.inf, raiseForReversible = False)
			# Reset the rate limit
			fba.maxReactionFluxIs(constraintToReactionDict[constraint1], np.inf, raiseForReversible = False)


def determineAllButOneConstraints(outputFileName, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	with open(outputFileName, 'w') as output_file:

		headers = ["Time Point (seconds)", "Removed ConstraintID", "ReactionID", " Removed Rate Limit", "Total Error", "Total Relative Error"] + [x + " Error" for x in metabolitePoolIDs] + [x + " Relative Error" for x in metabolitePoolIDs]
		output_file.write("\t".join(headers))

		# Find the unconstrained deviation from the target metabolite concentrations
		baseIndiviualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

		# Apply all constraints
		for index, constraint in enumerate(constraintIDs):
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)


		# Evaluate the error with all constraints
		individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

		# Find the error relative to the unconstrained solution
		indiviualRelativeErrors = individualErrors - baseIndiviualErrors
		totalRelativeError = totalError - baseTotalError

		row = [str(timepoint), "All constraints", "All Reactions", "No constraints removed", str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in indiviualRelativeErrors]

		output_file.write("\n" + "\t".join(row))


		# Remove each in turn and evaluate
		for index, constraint in enumerate(constraintIDs):
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)

			# Evaluate the error with these constraints
			individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			indiviualRelativeErrors = individualErrors - baseIndiviualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in individualErrors] + [str(x) for x in indiviualRelativeErrors]

			output_file.write("\n" + "\t".join(row))

			# Put this constraint back in place
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)


def determineErrorNRandomFluxesZero(fba, numZeroFluxes, constraint, reactionIDs, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):
	"""
		determineErrorNRandomFluxesZero

		Sets the flux of numZeroFluxes randomly choosen reactions to constraintPercentage
	"""

	# Choose a random vector of numZeroFluxes from reactionIDs
	randomFluxes = []
	for x in xrange(1,int(numZeroFluxes)):
		randomFluxes.append(np.random.choice(range(len(reactionIDs))))

	# Set the fluxes of these numZeroFluxes reactionIDs to zero
	for reaction in randomFluxes:
		# fba.maxReactionFluxIs(reactionIDs[reaction], constaintPercentage*fba.reactionFlux(reactionIDs[reaction]), raiseForReversible = False)
		fba.maxReactionFluxIs(reactionIDs[reaction], constraint, raiseForReversible = False)

	# Evaluate the error with these constraints
	individualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

	# Remove these constraints again
	for reaction in randomFluxes:
		fba.maxReactionFluxIs(reactionIDs[reaction], np.inf, raiseForReversible = False)

	# Record which reactions were constrained
	constrainedReactions = []
	for reaction in randomFluxes:
		constrainedReactions.append(reactionIDs[reaction])

	return individualErrors, totalError, constrainedReactions



def plotNumConstraintsVersusError(outputFileName, fba, samplesPerPoint, numPointsToSample, plotOutDir, plotOutFileName, reactionIDs, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint='none'):

	constraintNum = []
	totalRelativeErrors = []

	pointsToSample = np.arange(1,len(reactionIDs),len(reactionIDs)/numPointsToSample)

	# Find the unconstrained deviation from the target metabolite concentrations
	baseIndiviualErrors, baseTotalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

	positiveErrorReactions = set()
	negativeErrorReactions = set()

	for point in pointsToSample:
		for sample in range(samplesPerPoint):

			individualErrors, totalError, constrainedReactions = determineErrorNRandomFluxesZero(fba, point, 0, reactionIDs, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find the error relative to the unconstrained solution
			indiviualRelativeErrors = individualErrors - baseIndiviualErrors
			totalRelativeError = totalError - baseTotalError

			# Record this point
			constraintNum.append(int(point))
			totalRelativeErrors.append(totalRelativeError)

			if totalRelativeError > 0:
				for constrainedReaction in constrainedReactions:
					positiveErrorReactions.add(constrainedReaction)

			if totalRelativeError < 0:
				for constrainedReaction in constrainedReactions:
					negativeErrorReactions.add(constrainedReaction)

	import ipdb; ipdb.set_trace()

	with open(outputFileName, 'w') as output_file:
		output_file.write("Positive Error" + "\t")
		output_file.write('\t'.join([str(x) for x in positiveErrorReactions]))
		output_file.write('\n')
		output_file.write("Negative Error" + "\t")
		output_file.write('\t'.join([str(x) for x in negativeErrorReactions]))
		output_file.write('\n')
		output_file.write("Intersect" + "\t")
		output_file.write('\t'.join([str(x) for x in positiveErrorReactions.intersection(negativeErrorReactions)]))
		output_file.write('\n')
		output_file.write("Unique Positives" + "\t")
		output_file.write('\t'.join([str(x) for x in positiveErrorReactions.difference(positiveErrorReactions.intersection(negativeErrorReactions))]))
		output_file.write('\n')
		output_file.write("Unique Negatives" + "\t")
		output_file.write('\t'.join([str(x) for x in negativeErrorReactions.difference(positiveErrorReactions.intersection(negativeErrorReactions))]))

	# Find the mean and standard deviation of results at each point
	sampleVector = []
	sampleMeans = []
	sampleMeansAbs = []
	sampleStds = []
	sampleStdsAbs = []
	constraintNumAdjusted = []
	for index, point in enumerate(totalRelativeErrors):
		sampleVector.append(totalRelativeErrors[index])
		if len(sampleVector) == samplesPerPoint:
			sampleMeans.append(np.average(sampleVector))
			sampleMeansAbs.append(np.average(np.abs(sampleVector)))
			sampleStds.append(np.std(sampleVector))
			sampleStdsAbs.append(np.std(np.abs(sampleVector)))
			sampleVector = []
			constraintNumAdjusted.append(constraintNum[index])

	# Plot the results
	plt.figure(figsize = (8.5, 11))

	plt.subplot(3,1,1)
	plt.title("Random Constraints (Zero Flux), Timepoint= %d Seconds" % timepoint )
	plt.scatter(constraintNum, totalRelativeErrors, c='b')
	plt.ylabel("jFBA Error Relative to No Constraint)")

	plt.subplot(3,1,2)
	plt.bar(constraintNumAdjusted, sampleMeansAbs, yerr=sampleStdsAbs)
	plt.ylabel("Mean Absolute Value")

	plt.subplot(3,1,3)
	plt.bar(constraintNumAdjusted, sampleStds)
	plt.ylabel("Std Dev (non abs)")


	plt.xlabel("Number of Reactions Constrained")


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

