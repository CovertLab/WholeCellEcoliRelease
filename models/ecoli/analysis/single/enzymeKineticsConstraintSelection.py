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

DISABLED = False

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

			if findSingleConstraints:
				determineOneConstraintErrors(oneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			if findDoubleConstraints:
				determineTwoConstraintErrors(twoConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)

			# Find errors for all constraints and set of all constraints except one
			if findAllButOneConstraint:
				determineAllButOneConstraints(allbutOneConstraintOutputFilename, fba, constraints, constraintIDs, constraintToReactionDict, metaboliteCountsInit, targetConcentrations, metabolitePoolIDs, countsToMolar, timepoint)


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

	indiviualErrors = (metaboliteCountsFinal * countsToMolar) - targetConcentrations

	totalError = np.sum(np.absolute(indiviualErrors))

	return indiviualErrors, totalError


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
			indiviualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			indiviualRelativeErrors = indiviualErrors - baseIndiviualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in indiviualErrors] + [str(x) for x in indiviualRelativeErrors]

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
				indiviualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

				# Find the error relative to the unconstrained solution
				indiviualRelativeErrors = indiviualErrors - baseIndiviualErrors
				totalRelativeError = totalError - baseTotalError

				row = [str(timepoint), constraint1, constraint2, str(constraintToReactionDict[constraint1]), str(constraintToReactionDict[constraint2]), str(constraints[index1]), str(constraints[index2]), str(totalError), str(totalRelativeError)] + [str(x) for x in indiviualErrors] + [str(x) for x in indiviualRelativeErrors]

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
		indiviualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

		# Find the error relative to the unconstrained solution
		indiviualRelativeErrors = indiviualErrors - baseIndiviualErrors
		totalRelativeError = totalError - baseTotalError

		row = [str(timepoint), "All constraints", "All Reactions", "No constraints removed", str(totalError), str(totalRelativeError)] + [str(x) for x in indiviualErrors] + [str(x) for x in indiviualRelativeErrors]

		output_file.write("\n" + "\t".join(row))


		# Remove each in turn and evaluate
		for index, constraint in enumerate(constraintIDs):
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], np.inf, raiseForReversible = False)

			# Evaluate the error with these constraints
			indiviualErrors, totalError = evaluatePoint(fba, metaboliteCountsInit, targetConcentrations, countsToMolar)

			# Find the error relative to the unconstrained solution
			indiviualRelativeErrors = indiviualErrors - baseIndiviualErrors
			totalRelativeError = totalError - baseTotalError

			row = [str(timepoint), constraint, str(constraintToReactionDict[constraint]), str(constraints[index]), str(totalError), str(totalRelativeError)] + [str(x) for x in indiviualErrors] + [str(x) for x in indiviualRelativeErrors]

			output_file.write("\n" + "\t".join(row))

			# Put this constraint back in place
			fba.maxReactionFluxIs(constraintToReactionDict[constraint], constraints[index], raiseForReversible = False)