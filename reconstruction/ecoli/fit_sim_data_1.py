"""
The fitter, aka parameter calculator.
"""

from __future__ import absolute_import
from __future__ import division

import numpy as np
import os
import scipy.optimize
import cPickle

import wholecell
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.utils import filepath
from wholecell.utils import units
from wholecell.utils.fitting import normalize, massesAndCountsToAddForHomeostaticTargets

from cvxpy import Variable, Problem, Minimize, norm

from multiprocessing import Pool

# Tweaks
RNA_POLY_MRNA_DEG_RATE_PER_S = np.log(2) / 30. # half-life of 30 seconds
FRACTION_INCREASE_RIBOSOMAL_PROTEINS = 0.0  # reduce stochasticity from protein expression

NUMERICAL_ZERO = 1e-10

# TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)

FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100
N_SEEDS = 10

BASAL_EXPRESSION_CONDITION = "M9 Glucose minus AAs"

VERBOSE = 1

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

def fitSimData_1(raw_data, cpus=1, debug=False):
	'''
	Fits parameters necessary for the simulation based on the knowledge base

	Inputs:
		raw_data (KnowledgeBaseEcoli) - knowledge base consisting of the
			necessary raw data
		cpus (int) - number of processes to use (if >1 uses multiprocessing)
		debug (bool) - if True, fit only one arbitrarily-chosen transcription
			factor in order to speed up a debug cycle (should not be used for
			an actual simulation)
	'''

	sim_data = SimulationDataEcoli()
	sim_data.initialize(
		raw_data = raw_data,
		basal_expression_condition = BASAL_EXPRESSION_CONDITION,
		)

	# Limit the number of conditions that are being fit so that execution time decreases
	if debug:
		print("Warning: running fitter in debug mode - not all conditions will be fit")
		key = sim_data.tfToActiveInactiveConds.keys()[0]
		sim_data.tfToActiveInactiveConds = {key: sim_data.tfToActiveInactiveConds[key]}

	# Increase RNA poly mRNA deg rates
	setRnaPolymeraseCodingRnaDegradationRates(sim_data)

	# Make adjustments for metabolic enzymes
	setTranslationEfficiencies(sim_data)
	setRNAExpression(sim_data)
	setRNADegRates(sim_data)
	setProteinDegRates(sim_data)

	# Set C-period
	setCPeriod(sim_data)

	cellSpecs = buildBasalCellSpecifications(sim_data)

	# Modify other properties

	# Re-compute Km's
	if sim_data.constants.EndoRNaseCooperation:
		sim_data.process.transcription.rnaData["KmEndoRNase"] = setKmCooperativeEndoRNonLinearRNAdecay(sim_data, cellSpecs["basal"]["bulkContainer"])

	## Calculate and set maintenance values

	# ----- Growth associated maintenance -----
	fitMaintenanceCosts(sim_data, cellSpecs["basal"]["bulkContainer"])

	if cpus > 1:
		print "Start parallel processing with %i processes" % (cpus)
		pool = Pool(processes = cpus)
		results = [pool.apply_async(buildTfConditionCellSpecifications, (sim_data, tf)) for tf in sorted(sim_data.tfToActiveInactiveConds)]
		pool.close()
		pool.join()
		for result in results:
			assert(result.successful())
			cellSpecs.update(result.get())
		print "End parallel processing"
	else:
		for tf in sorted(sim_data.tfToActiveInactiveConds):
			cellSpecs.update(buildTfConditionCellSpecifications(sim_data, tf))

	for conditionKey in cellSpecs:
		if conditionKey == "basal":
			continue

		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

	buildCombinedConditionCellSpecifications(sim_data, cellSpecs)

	sim_data.process.transcription.rnaSynthProbFraction = {}
	sim_data.process.transcription.rnapFractionActiveDict = {}
	sim_data.process.transcription.rnaSynthProbRProtein = {}
	sim_data.process.transcription.rnaSynthProbRnaPolymerase = {}
	sim_data.process.transcription.rnaPolymeraseElongationRateDict = {}
	sim_data.expectedDryMassIncreaseDict = {}
	sim_data.process.translation.ribosomeElongationRateDict = {}
	sim_data.process.translation.ribosomeFractionActiveDict = {}

	# Fit kinetic parameters
	# findKineticCoeffs(sim_data, cellSpecs["basal"]["bulkContainer"])

	if cpus > 1:
		print "Start parallel processing with %i processes" % (cpus)
		pool = Pool(processes = cpus)
		results = [pool.apply_async(fitCondition, (sim_data, cellSpecs[condition], condition)) for condition in sorted(cellSpecs)]
		pool.close()
		pool.join()
		for result in results:
			assert(result.successful())
			cellSpecs.update(result.get())
		print "End parallel processing"
	else:
		for condition in sorted(cellSpecs):
			cellSpecs.update(fitCondition(sim_data, cellSpecs[condition], condition))

	for condition_label in sorted(cellSpecs):
		condition = sim_data.conditions[condition_label]
		if condition["nutrients"] not in sim_data.translationSupplyRate.keys():
			sim_data.translationSupplyRate[condition["nutrients"]] = cellSpecs[condition_label]["translation_aa_supply"]

	rVector = fitPromoterBoundProbability(sim_data, cellSpecs)

	for condition_label in sorted(cellSpecs):

		condition = sim_data.conditions[condition_label]

		if VERBOSE > 0:
			print "Updating mass in condition {}".format(condition_label)
		spec = cellSpecs[condition_label]

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			condition["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[condition_label]))

		avgCellDryMassInit, fitAvgSolublePoolMass = rescaleMassForSolubleMetabolites(
			sim_data, spec["bulkContainer"], concDict, sim_data.conditionToDoublingTime[condition_label]
			)

		if VERBOSE > 0:
			print str(spec["avgCellDryMassInit"]) + " to "  + str(avgCellDryMassInit)

		spec["avgCellDryMassInit"] = avgCellDryMassInit
		spec["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass

		mRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isMRna"]].sum()
		tRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isTRna"]].sum()
		rRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRRna"]].sum()

		if (condition["nutrients"] not in sim_data.process.transcription.rnaSynthProbFraction
				and len(condition["perturbations"]) == 0):
			sim_data.process.transcription.rnaSynthProbFraction[condition["nutrients"]] = {
				"mRna": mRnaSynthProb,
				"tRna": tRnaSynthProb,
				"rRna": rRnaSynthProb,
				}

		if (condition["nutrients"] not in sim_data.process.transcription.rnaSynthProbRProtein
				and len(condition["perturbations"]) == 0):
			sim_data.process.transcription.rnaSynthProbRProtein[condition["nutrients"]] = (
				sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRProtein"]]
				)

		if (condition["nutrients"] not in sim_data.process.transcription.rnaSynthProbRnaPolymerase
				and len(condition["perturbations"]) == 0):
			sim_data.process.transcription.rnaSynthProbRnaPolymerase[condition["nutrients"]] = (
				sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRnap"]]
				)

		if (condition["nutrients"] not in sim_data.process.transcription.rnapFractionActiveDict
				and len(condition["perturbations"]) == 0):
			sim_data.process.transcription.rnapFractionActiveDict[condition["nutrients"]]	= (
				sim_data.growthRateParameters.getFractionActiveRnap(spec["doubling_time"])
				)

		if (condition["nutrients"] not in sim_data.process.transcription.rnaPolymeraseElongationRateDict
				and len(condition["perturbations"]) == 0):
			sim_data.process.transcription.rnaPolymeraseElongationRateDict[condition["nutrients"]] = (
				sim_data.growthRateParameters.getRnapElongationRate(spec["doubling_time"])
				)

		if (condition["nutrients"] not in sim_data.expectedDryMassIncreaseDict
				and len(condition["perturbations"]) == 0):
			sim_data.expectedDryMassIncreaseDict[condition["nutrients"]] = spec["avgCellDryMassInit"]

		if (condition["nutrients"] not in sim_data.process.translation.ribosomeElongationRateDict
				and len(condition["perturbations"]) == 0):
			sim_data.process.translation.ribosomeElongationRateDict[condition["nutrients"]] = (
				sim_data.growthRateParameters.getRibosomeElongationRate(spec["doubling_time"])
				)

		if (condition["nutrients"] not in sim_data.process.translation.ribosomeFractionActiveDict
				and len(condition["perturbations"]) == 0):
			sim_data.process.translation.ribosomeFractionActiveDict[condition["nutrients"]] = (
				sim_data.growthRateParameters.getFractionActiveRibosome(spec["doubling_time"])
				)

	calculateRnapRecruitment(sim_data, cellSpecs, rVector)

	return sim_data


def buildBasalCellSpecifications(sim_data):
	cellSpecs = {}
	cellSpecs["basal"] = {
		"concDict": sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					"minimal"
				),
		"expression": sim_data.process.transcription.rnaExpression["basal"].copy(),
		"doubling_time": sim_data.conditionToDoublingTime["basal"],
		"translation_km": np.zeros(len(sim_data.moleculeGroups.aaIDs))
	}

	expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, _ = expressionConverge(
		sim_data,
		cellSpecs["basal"]["expression"],
		cellSpecs["basal"]["concDict"],
		cellSpecs["basal"]["doubling_time"],
		)

	cellSpecs["basal"]["expression"] = expression
	cellSpecs["basal"]["synthProb"] = synthProb
	cellSpecs["basal"]["avgCellDryMassInit"] = avgCellDryMassInit
	cellSpecs["basal"]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
	cellSpecs["basal"]["bulkContainer"] = bulkContainer

	sim_data.mass.avgCellDryMassInit = avgCellDryMassInit
	sim_data.mass.avgCellDryMass = sim_data.mass.avgCellDryMassInit * sim_data.mass.avgCellToInitialCellConvFactor
	sim_data.mass.avgCellWaterMassInit = sim_data.mass.avgCellDryMassInit / sim_data.mass.cellDryMassFraction * sim_data.mass.cellWaterMassFraction
	sim_data.mass.fitAvgSolubleTargetMolMass = fitAvgSolubleTargetMolMass


	sim_data.process.transcription.rnaExpression["basal"][:] = cellSpecs["basal"]["expression"]
	sim_data.process.transcription.rnaSynthProb["basal"][:] = cellSpecs["basal"]["synthProb"]

	translation_aa_supply = calculateTranslationSupply(
									sim_data,
									cellSpecs["basal"]["doubling_time"],
									cellSpecs["basal"]["bulkContainer"],
									cellSpecs["basal"]["avgCellDryMassInit"],
									)

	return cellSpecs

def calculateTranslationSupply(sim_data, doubling_time, bulkContainer, avgCellDryMassInit):
	aaCounts = sim_data.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomerData["id"])
	nAvogadro = sim_data.constants.nAvogadro

	molAAPerGDCW = (
			units.sum(
				aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * nAvogadro)) *
				(1 / avgCellDryMassInit)
			)
		)

	translation_aa_supply = molAAPerGDCW * np.log(2) / doubling_time

	return translation_aa_supply

def buildTfConditionCellSpecifications(sim_data, tf):
	cellSpecs = {}
	for choice in ["__active", "__inactive"]:
		conditionKey = tf + choice
		conditionValue = sim_data.conditions[conditionKey]

		fcData = {}
		if choice == "__active" and conditionValue != sim_data.conditions["basal"]:
			fcData = sim_data.tfToFC[tf]
		if choice == "__inactive" and conditionValue != sim_data.conditions["basal"]:
			fcDataTmp = sim_data.tfToFC[tf].copy()
			for key, value in fcDataTmp.iteritems():
				fcData[key] = 1. / value
		expression = expressionFromConditionAndFoldChange(
			sim_data.process.transcription.rnaData["id"],
			sim_data.process.transcription.rnaExpression["basal"],
			conditionValue["perturbations"],
			fcData,
		)

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[conditionKey]))

		cellSpecs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.conditionToDoublingTime.get(
				conditionKey,
				sim_data.conditionToDoublingTime["basal"]
			)
		}

		expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cellSpecs[conditionKey]["expression"],
			cellSpecs[conditionKey]["concDict"],
			cellSpecs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rnaData["KmEndoRNase"],
			updateConcDict = True,
			)

		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

	return cellSpecs

		# Uncomment when concDict is actually calculated for non-base [AA]
		# if len(conditionValue["perturbations"]) == 0:
		# 	nutrientLabel = conditionValue["nutrients"]
		# 	sim_data.process.metabolism.nutrientsToInternalConc[nutrientLabel] = concDict

def buildCombinedConditionCellSpecifications(sim_data, cellSpecs):
	fcData = {}
	for conditionKey in sim_data.conditionActiveTfs:
		if conditionKey == "basal":
			continue

		conditionValue = sim_data.conditions[conditionKey]
		for tf in sim_data.conditionActiveTfs[conditionKey]:
			for gene in sim_data.tfToFC[tf]:
				if gene in fcData:
					# TODO: multiply if multiple genes affected?
					# fcData[gene] *= sim_data.tfToFC[tf][gene]
					raise Exception("Derek check/implement this: multiple genes regulated")
				else:
					fcData[gene] = sim_data.tfToFC[tf][gene]

		expression = expressionFromConditionAndFoldChange(
			sim_data.process.transcription.rnaData["id"],
			sim_data.process.transcription.rnaExpression["basal"],
			conditionValue["perturbations"],
			fcData,
		)

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[conditionKey]))

		cellSpecs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.conditionToDoublingTime.get(
				conditionKey,
				sim_data.conditionToDoublingTime["basal"]
			)
		}

		expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cellSpecs[conditionKey]["expression"],
			cellSpecs[conditionKey]["concDict"],
			cellSpecs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rnaData["KmEndoRNase"],
			updateConcDict = True,
			)

		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

		# Uncomment when concDict is actually calculated for non-base [AA]
		# if len(conditionValue["perturbations"]) == 0:
		# 	nutrientLabel = conditionValue["nutrients"]
		# 	sim_data.process.metabolism.nutrientsToInternalConc[nutrientLabel] = concDict


def expressionConverge(sim_data, expression, concDict, doubling_time, Km = None, updateConcDict = False):
	# Fit synthesis probabilities for RNA
	if VERBOSE > 0:
		print "Fitting RNA synthesis probabilities."
	for iteration in xrange(MAX_FITTING_ITERATIONS):
		if VERBOSE > 1: print 'Iteration: {}'.format(iteration)

		initialExpression = expression.copy()

		expression = setInitialRnaExpression(sim_data, expression, doubling_time)

		bulkContainer = createBulkContainer(sim_data, expression, doubling_time)

		avgCellDryMassInit, fitAvgSolubleTargetMolMass = rescaleMassForSolubleMetabolites(sim_data, bulkContainer, concDict, doubling_time)

		setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time)

		setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		# Normalize expression and write out changes

		expression, synthProb = fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		if updateConcDict:
			concDict = concDict.copy() # Calculate non-base condition [AA]

		finalExpression = expression

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - finalExpression)))
		if VERBOSE > 1: print 'degree of fit: {}'.format(degreeOfFit)

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	return expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict

def fitCondition(sim_data, spec, condition):
	if VERBOSE > 0:
		print "Fitting condition {}".format(condition)

	bulkAverageContainer, bulkDeviationContainer, proteinMonomerAverageContainer, proteinMonomerDeviationContainer = calculateBulkDistributions(
		sim_data,
		spec["expression"],
		spec["concDict"],
		spec["avgCellDryMassInit"],
		spec["doubling_time"],
		)
	spec["bulkAverageContainer"] = bulkAverageContainer
	spec["bulkDeviationContainer"] = bulkDeviationContainer
	spec["proteinMonomerAverageContainer"] = proteinMonomerAverageContainer
	spec["proteinMonomerDeviationContainer"] = proteinMonomerDeviationContainer

	spec["translation_aa_supply"] = calculateTranslationSupply(
											sim_data,
											spec["doubling_time"],
											spec["proteinMonomerAverageContainer"],
											spec["avgCellDryMassInit"],
											)

	return {condition: spec}

# Sub-fitting functions

def setRnaPolymeraseCodingRnaDegradationRates(sim_data):
	# Increase RNA poly mRNA deg rates
	# TODO: set this based on transcription unit structure
	# i.e. same synthesis prob. but different deg rates

	rnaPolySubunits = sim_data.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"]
	subunitIndexes = np.array([np.where(sim_data.process.translation.monomerData["id"] == id_)[0].item() for id_ in rnaPolySubunits]) # there has to be a better way...
	mRNA_indexes = sim_data.relation.rnaIndexToMonomerMapping[subunitIndexes]
	sim_data.process.transcription.rnaData.struct_array["degRate"][mRNA_indexes] = RNA_POLY_MRNA_DEG_RATE_PER_S

def setTranslationEfficiencies(sim_data):
	adjustments = {
		"ADCLY-MONOMER[c]": 5,
		"EG12438-MONOMER[c]": 5,
		"EG12298-MONOMER[p]": 5, # for anaerobic condition
		"ACETYL-COA-ACETYLTRANSFER-MONOMER[c]": 5, # for anaerobic condition
		}

	for protein in adjustments:
		idx = np.where(sim_data.process.translation.monomerData["id"] == protein)[0]
		sim_data.process.translation.translationEfficienciesByMonomer[idx] *= adjustments[protein]

def setRNAExpression(sim_data):
	adjustments = {
		"EG11493_RNA[c]": 10,
		"EG12438_RNA[c]": 10,
		"EG10139_RNA[c]": 10,
		"EG12298_RNA[c]": 10, # for anaerobic condition
		"EG11672_RNA[c]": 10, # for anaerobic condition
		}

	for rna in adjustments:
		idx = np.where(sim_data.process.transcription.rnaData["id"] == rna)[0]
		sim_data.process.transcription.rnaExpression["basal"][idx] *= adjustments[rna]

	sim_data.process.transcription.rnaExpression["basal"] /= sim_data.process.transcription.rnaExpression["basal"].sum()

def setRNADegRates(sim_data):
	adjustments = {
		"EG11493_RNA[c]": 2,
		"EG10139_RNA[c]": 2,
		}

	for rna in adjustments:
		idx = np.where(sim_data.process.transcription.rnaData["id"] == rna)[0]
		sim_data.process.transcription.rnaData.struct_array["degRate"][idx] *= adjustments[rna]

def setProteinDegRates(sim_data):
	adjustments = {
		"EG12298-MONOMER[p]": 0.1, # for anaerobic condition
		}

	for protein in adjustments:
		idx = np.where(sim_data.process.translation.monomerData["id"] == protein)[0]
		sim_data.process.translation.monomerData.struct_array["degRate"][idx] *= adjustments[protein]

def setCPeriod(sim_data):
	sim_data.growthRateParameters.c_period = sim_data.process.replication.genome_length * units.nt / sim_data.growthRateParameters.dnaPolymeraseElongationRate / 2

def rescaleMassForSolubleMetabolites(sim_data, bulkMolCntr, concDict, doubling_time):
	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)

	mass = (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	targetMoleculeIds = sorted(concDict)
	targetMoleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in targetMoleculeIds])

	massesToAdd, countsToAdd = massesAndCountsToAddForHomeostaticTargets(
		mass,
		targetMoleculeIds,
		targetMoleculeConcentrations,
		sim_data.getter.getMass(targetMoleculeIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		targetMoleculeIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculetargetMoleculesDryMass = units.hstack((massesToAdd[:targetMoleculeIds.index('WATER[c]')], massesToAdd[targetMoleculeIds.index('WATER[c]') + 1:]))
	newAvgCellDryMassInit = units.sum(mass) + units.sum(smallMoleculetargetMoleculesDryMass)
	fitAvgSolubleTargetMolMass = units.sum(units.hstack((massesToAdd[:targetMoleculeIds.index('WATER[c]')], massesToAdd[targetMoleculeIds.index('WATER[c]') + 1:]))) * sim_data.mass.avgCellToInitialCellConvFactor

	return newAvgCellDryMassInit, fitAvgSolubleTargetMolMass

def setInitialRnaExpression(sim_data, expression, doubling_time):
	# Set expression for all of the noncoding RNAs

	# Load from KB

	## IDs
	ids_rnas = sim_data.process.transcription.rnaData["id"]
	ids_rRNA23S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna23S"]]
	ids_rRNA16S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna16S"]]
	ids_rRNA5S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna5S"]]
	ids_tRNA = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isTRna"]]
	ids_mRNA = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isMRna"]]

	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)

	## Mass fractions
	totalMass_rRNA23S = avgCellFractionMass["rRna23SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA16S = avgCellFractionMass["rRna16SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA5S = avgCellFractionMass["rRna5SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_tRNA = avgCellFractionMass["tRnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_mRNA = avgCellFractionMass["mRnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	## Molecular weights
	individualMasses_RNA = sim_data.getter.getMass(ids_rnas) / sim_data.constants.nAvogadro
	individualMasses_rRNA23S = sim_data.getter.getMass(ids_rRNA23S) / sim_data.constants.nAvogadro
	individualMasses_rRNA16S = sim_data.getter.getMass(ids_rRNA16S) / sim_data.constants.nAvogadro
	individualMasses_rRNA5S = sim_data.getter.getMass(ids_rRNA5S) / sim_data.constants.nAvogadro
	individualMasses_tRNA = sim_data.process.transcription.rnaData["mw"][sim_data.process.transcription.rnaData["isTRna"]] / sim_data.constants.nAvogadro
	individualMasses_mRNA = sim_data.process.transcription.rnaData["mw"][sim_data.process.transcription.rnaData["isMRna"]] / sim_data.constants.nAvogadro

	## Molecule expression distributions
	distribution_rRNA23S = np.array([1.] + [0.] * (ids_rRNA23S.size-1)) # all expression from first rRNA operon
	distribution_rRNA16S = np.array([1.] + [0.] * (ids_rRNA16S.size-1)) # all expression from first rRNA operon
	distribution_rRNA5S = np.array([1.] + [0.] * (ids_rRNA5S.size-1)) # all expression from first rRNA operon
	distribution_tRNA = normalize(sim_data.mass.getTrnaDistribution()['molar_ratio_to_16SrRNA'])
	distribution_mRNA = normalize(expression[sim_data.process.transcription.rnaData['isMRna']])

	# Construct bulk container

	rnaExpressionContainer = BulkObjectsContainer(ids_rnas, dtype = np.float64)

	## Assign rRNA counts based on mass

	totalCount_rRNA23S = totalCountFromMassesAndRatios(
		totalMass_rRNA23S,
		individualMasses_rRNA23S,
		distribution_rRNA23S
		)

	totalCount_rRNA23S.normalize()
	totalCount_rRNA23S.checkNoUnit()

	totalCount_rRNA16S = totalCountFromMassesAndRatios(
		totalMass_rRNA16S,
		individualMasses_rRNA16S,
		distribution_rRNA16S
		)

	totalCount_rRNA16S.normalize()
	totalCount_rRNA16S.checkNoUnit()

	totalCount_rRNA5S = totalCountFromMassesAndRatios(
		totalMass_rRNA5S,
		individualMasses_rRNA5S,
		distribution_rRNA5S
		)

	totalCount_rRNA5S.normalize()
	totalCount_rRNA5S.checkNoUnit()

	totalCount_rRNA_average = sum([totalCount_rRNA23S, totalCount_rRNA16S, totalCount_rRNA5S]) / 3

	counts_rRNA23S = totalCount_rRNA_average * distribution_rRNA23S
	counts_rRNA16S = totalCount_rRNA_average * distribution_rRNA16S
	counts_rRNA5S = totalCount_rRNA_average * distribution_rRNA5S

	rnaExpressionContainer.countsIs(counts_rRNA23S, ids_rRNA23S)
	rnaExpressionContainer.countsIs(counts_rRNA16S, ids_rRNA16S)
	rnaExpressionContainer.countsIs(counts_rRNA5S, ids_rRNA5S)

	## Assign tRNA counts based on mass and relative abundances (see Dong 1996)

	totalCount_tRNA = totalCountFromMassesAndRatios(
		totalMass_tRNA,
		individualMasses_tRNA,
		distribution_tRNA
		)

	totalCount_tRNA.normalize()
	totalCount_tRNA.checkNoUnit()

	counts_tRNA = totalCount_tRNA * distribution_tRNA

	rnaExpressionContainer.countsIs(counts_tRNA, ids_tRNA)

	## Assign mRNA counts based on mass and relative abundances (microarrays)

	totalCount_mRNA = totalCountFromMassesAndRatios(
		totalMass_mRNA,
		individualMasses_mRNA,
		distribution_mRNA
		)

	totalCount_mRNA.normalize()
	totalCount_mRNA.checkNoUnit()

	counts_mRNA = totalCount_mRNA * distribution_mRNA

	rnaExpressionContainer.countsIs(counts_mRNA, ids_mRNA)

	expression = normalize(rnaExpressionContainer.counts())

	return expression
	# Note that now rnaData["synthProb"] does not match "expression"

def totalCountIdDistributionProtein(sim_data, expression, doubling_time):
	ids_protein = sim_data.process.translation.monomerData["id"]
	totalMass_protein = sim_data.mass.getFractionMass(doubling_time)["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_protein = sim_data.process.translation.monomerData["mw"] / sim_data.constants.nAvogadro
	distribution_transcriptsByProtein = normalize(expression[sim_data.relation.rnaIndexToMonomerMapping])
	translation_efficienciesByProtein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	degradationRates = sim_data.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doubling_time, degradationRates)

	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcriptsByProtein,
		translation_efficienciesByProtein,
		netLossRate_protein
		)

	totalCount_protein = totalCountFromMassesAndRatios(
		totalMass_protein,
		individualMasses_protein,
		distribution_protein
		)

	totalCount_protein.normalize()
	totalCount_protein.checkNoUnit()

	return totalCount_protein, ids_protein, distribution_protein

def totalCountIdDistributionRNA(sim_data, expression, doubling_time):
	ids_rnas = sim_data.process.transcription.rnaData["id"]
	totalMass_RNA = sim_data.mass.getFractionMass(doubling_time)["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_RNA = sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro

	distribution_RNA = normalize(expression)

	totalCount_RNA = totalCountFromMassesAndRatios(
		totalMass_RNA,
		individualMasses_RNA,
		distribution_RNA
		)
	totalCount_RNA.normalize()
	totalCount_RNA.checkNoUnit()

	return totalCount_RNA, ids_rnas, distribution_RNA

def createBulkContainer(sim_data, expression, doubling_time):

	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)
	ids_molecules = sim_data.internal_state.bulkMolecules.bulkData["id"]

	## Construct bulk container

	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	## Assign RNA counts based on mass and expression distribution

	counts_RNA = totalCount_RNA * distribution_RNA

	bulkContainer.countsIs(counts_RNA, ids_rnas)

	## Assign protein counts based on mass and mRNA counts

	counts_protein = totalCount_protein * distribution_protein

	bulkContainer.countsIs(counts_protein, ids_protein)

	return bulkContainer


def setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time):
	'''
	setRibosomeCountsConstrainedByPhysiology

	Methodology: Set counts of ribosomal subunits based on three constraints.
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal subunit counts based on expression
	'''
	ribosome30SSubunits = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s30_fullComplex[0])['subunitIds']
	ribosome50SSubunits = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s50_fullComplex[0])['subunitIds']
	ribosome30SStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s30_fullComplex[0])['subunitStoich']
	ribosome50SStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.s50_fullComplex[0])['subunitStoich']

	# -- CONSTRAINT 1: Expected protien distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(sim_data.process.translation.monomerData['aaCounts'], axis = 1)
	proteinDegradationRates =  sim_data.process.translation.monomerData["degRate"]
	proteinCounts =  bulkContainer.counts(sim_data.process.translation.monomerData["id"])

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(
		doubling_time,
		proteinDegradationRates
		)

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
	proteinLengths, sim_data.growthRateParameters.getRibosomeElongationRate(doubling_time), netLossRate_protein, proteinCounts)
	nRibosomesNeeded.normalize() # FIXES NO UNIT BUG
	nRibosomesNeeded.checkNoUnit()
	nRibosomesNeeded = nRibosomesNeeded.asNumber()

	# Minimum number of ribosomes needed
	constraint1_ribosome30SCounts = (
		nRibosomesNeeded * ribosome30SStoich
		) * (1 + FRACTION_INCREASE_RIBOSOMAL_PROTEINS)

	constraint1_ribosome50SCounts = (
		nRibosomesNeeded * ribosome50SStoich
		) * (1 + FRACTION_INCREASE_RIBOSOMAL_PROTEINS)


	# -- CONSTRAINT 2: Measured rRNA mass fraction -- #
	## Calculate exact number of 30S and 50S subunits based on measured mass fractions of
	## 16S, 23S, and 5S rRNA.
	rRna23SCounts = bulkContainer.counts(sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna23S"]])
	rRna16SCounts = bulkContainer.counts(sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna16S"]])
	rRna5SCounts = bulkContainer.counts(sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna5S"]])

	## 16S rRNA is in the 30S subunit
	massFracPredicted_30SCount = rRna16SCounts.sum()
	## 23S and 5S rRNA are in the 50S subunit
	massFracPredicted_50SCount = min(rRna23SCounts.sum(), rRna5SCounts.sum())

	constraint2_ribosome30SCounts = massFracPredicted_30SCount * ribosome30SStoich
	constraint2_ribosome50SCounts = massFracPredicted_50SCount * ribosome50SStoich



	# -- CONSTRAINT 3: Expected ribosomal subunit counts based distribution
	## Calculate fundamental ribosomal subunit count distribution based on RNA expression data
	## Already calculated and stored in bulkContainer
	ribosome30SCounts = bulkContainer.counts(ribosome30SSubunits)
	ribosome50SCounts = bulkContainer.counts(ribosome50SSubunits)

	# -- SET RIBOSOME FUNDAMENTAL SUBUNIT COUNTS TO MAXIMUM CONSTRAINT -- #
	constraint_names = np.array(["Insufficient to double protein counts", "Too small for mass fraction", "Current level OK"])
	nRibosomesNeeded = nRibosomesNeeded * (1 + FRACTION_INCREASE_RIBOSOMAL_PROTEINS)
	rib30lims = np.array([nRibosomesNeeded, massFracPredicted_30SCount, (ribosome30SCounts / ribosome30SStoich).min()])
	rib50lims = np.array([nRibosomesNeeded, massFracPredicted_50SCount, (ribosome50SCounts / ribosome50SStoich).min()])
	if VERBOSE > 1: print '30S limit: {}'.format(constraint_names[np.where(rib30lims.max() == rib30lims)[0]][-1])
	if VERBOSE > 1: print '30S actual count: {}'.format((ribosome30SCounts / ribosome30SStoich).min())
	if VERBOSE > 1: print '30S count set to: {}'.format(rib30lims[np.where(rib30lims.max() == rib30lims)[0]][-1])
	if VERBOSE > 1: print '50S limit: {}'.format(constraint_names[np.where(rib50lims.max() == rib50lims)[0]][-1])
	if VERBOSE > 1: print '50S actual count: {}'.format((ribosome50SCounts / ribosome50SStoich).min())
	if VERBOSE > 1: print '50S count set to: {}'.format(rib50lims[np.where(rib50lims.max() == rib50lims)[0]][-1])

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome30SCounts, constraint1_ribosome30SCounts), constraint2_ribosome30SCounts),
		ribosome30SSubunits
		)

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome50SCounts, constraint1_ribosome50SCounts), constraint2_ribosome50SCounts),
		ribosome50SSubunits
		)

	# Fix rRNA counts
	bulkContainer.countsIs(rRna23SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna23S"]])
	bulkContainer.countsIs(rRna16SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna16S"]])
	bulkContainer.countsIs(rRna5SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna5S"]])


def setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km = None):
	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(sim_data.process.transcription.rnaData['countsACGU'], axis = 1)

	rnaLossRate = None

	if Km is None:
		rnaLossRate = netLossRateFromDilutionAndDegradationRNALinear(
			doubling_time,
			sim_data.process.transcription.rnaData["degRate"],
			bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		)
	else:
		# Get constants to compute countsToMolar factor
		cellDensity = sim_data.constants.cellDensity
		cellVolume = avgCellDryMassInit / cellDensity / 0.3
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

		rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
			doubling_time,
			(1 / countsToMolar) * totalEndoRnaseCapacity,
			Km,
			rnaConc,
			countsToMolar,
			)

	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths, sim_data.growthRateParameters.getRnapElongationRate(doubling_time), rnaLossRate)

	nActiveRnapNeeded = units.convertNoUnitToNumber(nActiveRnapNeeded)
	nRnapsNeeded = nActiveRnapNeeded / sim_data.growthRateParameters.getFractionActiveRnap(doubling_time)

	rnapIds = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])['subunitIds']
	rnapStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])['subunitStoich']

	minRnapSubunitCounts = (
		nRnapsNeeded * rnapStoich # Subunit stoichiometry
		) * (1 + sim_data.growthRateParameters.getFractionIncreaseRnapProteins(doubling_time))

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(rnapIds)

	## -- SET RNAP COUNTS TO MAXIMIM CONSTRAINTS -- #
	constraint_names = np.array(["Current level OK", "Insufficient to double RNA distribution"])
	rnapLims = np.array([(rnapCounts / rnapStoich).min(), (minRnapSubunitCounts / rnapStoich).min()])
	if VERBOSE > 1: print 'rnap limit: {}'.format(constraint_names[np.where(rnapLims.max() == rnapLims)[0]][0])
	if VERBOSE > 1: print 'rnap actual count: {}'.format((rnapCounts / rnapStoich).min())
	if VERBOSE > 1: print 'rnap counts set to: {}'.format(rnapLims[np.where(rnapLims.max() == rnapLims)[0]][0])

	bulkContainer.countsIs(minRnapSubunitCounts, rnapIds)

def fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km = None):

	view_RNA = bulkContainer.countsView(sim_data.process.transcription.rnaData["id"])
	counts_protein = bulkContainer.counts(sim_data.process.translation.monomerData["id"])

	translation_efficienciesByProtein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)
	totalMass_RNA = avgCellFractionMass["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	degradationRates_protein = sim_data.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doubling_time, degradationRates_protein)

	### Modify sim_dataFit to reflect our bulk container ###

	## RNA and monomer expression ##
	rnaExpressionContainer = BulkObjectsContainer(list(sim_data.process.transcription.rnaData["id"]), dtype = np.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(view_RNA.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		sim_data.process.translation.monomerData["rnaId"][sim_data.relation.monomerIndexToRnaMapping] == sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids" # TODO: move to KB tests

	mRnaExpressionView = rnaExpressionContainer.countsView(sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * mRNADistributionFromProtein(
			normalize(counts_protein), translation_efficienciesByProtein, netLossRate_protein
			)[sim_data.relation.monomerIndexToRnaMapping]
		)

	expression = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = totalCountFromMassesAndRatios(
		totalMass_RNA,
		sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro,
		expression
		)

	nRnas.normalize()
	nRnas.checkNoUnit()

	view_RNA.countsIs(nRnas * expression)

	rnaLossRate = None
	if Km is None:
		rnaLossRate = netLossRateFromDilutionAndDegradationRNALinear(
			doubling_time,
			sim_data.process.transcription.rnaData["degRate"],
			bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		)
	else:
		# Get constants to compute countsToMolar factor
		cellDensity = sim_data.constants.cellDensity
		cellVolume = avgCellDryMassInit / cellDensity / 0.3
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

		rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
			doubling_time,
			(1 / countsToMolar) * totalEndoRnaseCapacity,
			Km,
			countsToMolar * view_RNA.counts(),
			countsToMolar,
		)

	synthProb = normalize(rnaLossRate.asNumber(1 / units.min))

	return expression, synthProb


def fitMaintenanceCosts(sim_data, bulkContainer):
	aaCounts = sim_data.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomerData["id"])
	nAvogadro = sim_data.constants.nAvogadro
	avgCellDryMassInit = sim_data.mass.avgCellDryMassInit
	gtpPerTranslation = sim_data.constants.gtpPerTranslation

	# GTPs used for translation (recycled, not incorporated into biomass)
	aaMmolPerGDCW = (
			units.sum(
				aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * nAvogadro)) *
				(1 / avgCellDryMassInit)
			)
		)

	aasUsedOverCellCycle = units.sum(aaMmolPerGDCW)
	gtpUsedOverCellCycleMmolPerGDCW = gtpPerTranslation * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		sim_data.constants.growthAssociatedMaintenance -
		gtpUsedOverCellCycleMmolPerGDCW
		)

	additionalGtpPerTranslation = darkATP / aasUsedOverCellCycle
	additionalGtpPerTranslation.normalize()
	additionalGtpPerTranslation.checkNoUnit()
	additionalGtpPerTranslation = additionalGtpPerTranslation.asNumber()

	# Assign the growth associated "dark energy" to translation
	# TODO: Distribute it amongst growth-related processes
	# sim_data.constants.gtpPerTranslation += additionalGtpPerTranslation

	sim_data.constants.darkATP = darkATP

def calculateBulkDistributions(sim_data, expression, concDict, avgCellDryMassInit, doubling_time):

	# Ids
	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)
	ids_complex = sim_data.process.complexation.moleculeNames
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_twoComponentSystem = sim_data.process.two_component_system.moleculeNames
	ids_metabolites = sorted(concDict)
	conc_metabolites = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in ids_metabolites])
	allMoleculesIDs = sorted(
		set(ids_rnas) | set(ids_protein) | set(ids_complex) | set(ids_equilibrium) | set(ids_twoComponentSystem) | set(ids_metabolites)
		)

	# Data for complexation
	complexationStoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64, order = "F")

	complexationPrebuiltMatrices = mccBuildMatrices(
		complexationStoichMatrix
		)

	# Data for equilibrium binding
	equilibriumDerivatives = sim_data.process.equilibrium.derivatives
	equilibriumDerivativesJacobian = sim_data.process.equilibrium.derivativesJacobian

	# Data for metabolites
	cellDensity = sim_data.constants.cellDensity
	cellVolume = avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction

	# Construct bulk container

	# We want to know something about the distribution of the copy numbers of
	# macromolecules in the cell.  While RNA and protein expression can be
	# approximated using well-described statistical distributions, we need
	# absolute copy numbers to form complexes.  To get a distribution, we must
	# instantiate many cells, form complexes, and finally compute the
	# statistics we will use in the fitting operations.

	bulkContainer = BulkObjectsContainer(sim_data.internal_state.bulkMolecules.bulkData['id'])
	rnaView = bulkContainer.countsView(ids_rnas)
	proteinView = bulkContainer.countsView(ids_protein)
	complexationMoleculesView = bulkContainer.countsView(ids_complex)
	equilibriumMoleculesView = bulkContainer.countsView(ids_equilibrium)
	twoComponentSystemMoleculesView = bulkContainer.countsView(ids_twoComponentSystem)
	metabolitesView = bulkContainer.countsView(ids_metabolites)
	allMoleculesView = bulkContainer.countsView(allMoleculesIDs)

	allMoleculeCounts = np.empty((N_SEEDS, allMoleculesView.counts().size), np.int64)
	proteinMonomerCounts = np.empty((N_SEEDS, proteinView.counts().size), np.int64)

	if VERBOSE > 1:
		print "Bulk distribution seed:"

	for seed in xrange(N_SEEDS):
		if VERBOSE > 1:
			print seed
		randomState = np.random.RandomState(seed)

		allMoleculesView.countsIs(0)

		# rnaView.countsIs(randomState.multinomial(
		# 	totalCount_RNA,
		# 	distribution_RNA
		# 	))

		# proteinView.countsIs(randomState.multinomial(
		# 	totalCount_protein,
		# 	distribution_protein
		# 	))

		rnaView.countsIs(totalCount_RNA * distribution_RNA)

		proteinView.countsIs(totalCount_protein * distribution_protein)
		proteinMonomerCounts[seed, :] = proteinView.counts()
		complexationMoleculeCounts = complexationMoleculesView.counts()

		updatedCompMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			complexationMoleculeCounts,
			seed,
			complexationStoichMatrix,
			*complexationPrebuiltMatrices
			)

		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		metDiffs = np.inf * np.ones_like(metabolitesView.counts())
		nIters = 0

		while(np.linalg.norm(metDiffs, np.inf) > 1):
			# Metabolite concentrations were measured as steady-state values (not initial values)
			# So we run this until we get to steady state
			metCounts = conc_metabolites * cellVolume * sim_data.constants.nAvogadro
			metCounts.normalize()
			metCounts.checkNoUnit()

			metabolitesView.countsIs(
				metCounts.asNumber().round()
				)

			rxnFluxes, _ = sim_data.process.equilibrium.fluxesAndMoleculesToSS(
				equilibriumMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.nAvogadro.asNumber(1 / units.mol),
				)

			equilibriumMoleculesView.countsInc(
				np.dot(sim_data.process.equilibrium.stoichMatrix().astype(np.int64), rxnFluxes)
				)

			assert np.all(equilibriumMoleculesView.counts() >= 0)

			_, moleculeCountChanges = sim_data.process.two_component_system.moleculesToSS(
				twoComponentSystemMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.nAvogadro.asNumber(1 / units.mmol),
				1e6,
				)
			twoComponentSystemMoleculesView.countsInc(moleculeCountChanges)

			metDiffs = metabolitesView.counts() - metCounts.asNumber().round()

			nIters += 1
			if nIters > 100:
				raise Exception, "Equilibrium reactions are not converging!"


		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	bulkAverageContainer = BulkObjectsContainer(sim_data.internal_state.bulkMolecules.bulkData['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(sim_data.internal_state.bulkMolecules.bulkData['id'], np.float64)
	proteinMonomerAverageContainer = BulkObjectsContainer(sim_data.process.translation.monomerData["id"], np.float64)
	proteinMonomerDeviationContainer = BulkObjectsContainer(sim_data.process.translation.monomerData["id"], np.float64)

	bulkAverageContainer.countsIs(allMoleculeCounts.mean(0), allMoleculesIDs)
	bulkDeviationContainer.countsIs(allMoleculeCounts.std(0), allMoleculesIDs)
	proteinMonomerAverageContainer.countsIs(proteinMonomerCounts.mean(0), sim_data.process.translation.monomerData["id"])
	proteinMonomerDeviationContainer.countsIs(proteinMonomerCounts.std(0), sim_data.process.translation.monomerData["id"])

	return bulkAverageContainer, bulkDeviationContainer, proteinMonomerAverageContainer, proteinMonomerDeviationContainer


# Math functions

def totalCountFromMassesAndRatios(totalMass, individualMasses, distribution):
	"""
	Total mass = dot(mass, count)

	Fraction of i:
	f = count / Total counts

	Substituting:
	Total mass = dot(mass, f * Total counts)
	Total mass = Total counts * dot(mass, f)

	Total counts = Total mass / dot(mass, f)
	"""
	assert np.allclose(np.sum(distribution), 1)
	return 1 / units.dot(individualMasses, distribution) * totalMass


def proteinDistributionFrommRNA(distribution_mRNA, translation_efficiencies, netLossRate):
	"""
	dP_i / dt = k * M_i * e_i - P_i * Loss_i

	At steady state:
	P_i = k * M_i * e_i / Loss_i

	Fraction of mRNA for ith gene is defined as:
	f_i = M_i / M_total

	Substituting in:
	P_i = k * f_i * e_i * M_total / Loss_i

	Normalizing P_i by summing over all i cancels out k and M_total
	assuming constant translation rate.
	"""

	assert np.allclose(np.sum(distribution_mRNA), 1)
	assert np.allclose(np.sum(translation_efficiencies), 1)
	distributionUnnormed = 1 / netLossRate * distribution_mRNA * translation_efficiencies
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()
	return distributionNormed.asNumber()


def mRNADistributionFromProtein(distribution_protein, translation_efficiencies, netLossRate):
	"""
	dP_i / dt = k * M_i * e_i - P_i * Loss_i

	At steady state:
	M_i = Loss_i * P_i / (k * e_i)

	Fraction of protein for ith gene is defined as:
	f_i = P_i / P_total

	Substituting in:
	M_i = Loss_i * f_i * P_total / (k * e_i)

	Normalizing M_i by summing over all i cancles out k and P_total
	assuming a constant translation rate.

	"""
	assert np.allclose(np.sum(distribution_protein), 1)
	distributionUnnormed = netLossRate * distribution_protein / translation_efficiencies
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()
	return distributionNormed.asNumber()


def calculateMinPolymerizingEnzymeByProductDistribution(productLengths, elongationRate, netLossRate, productCounts):
	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRate
			* netLossRate
			* productCounts
		)
	return nPolymerizingEnzymeNeeded

def calculateMinPolymerizingEnzymeByProductDistributionRNA(productLengths, elongationRate, netLossRate):
	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRate
			* netLossRate
		)
	return nPolymerizingEnzymeNeeded


def netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates):
	return np.log(2) / doublingTime + degradationRates


def netLossRateFromDilutionAndDegradationRNA(doublingTime, totalEndoRnaseCountsCapacity, Km, rnaConc, countsToMolar):
	fracSaturated = rnaConc / Km / (1 + units.sum(rnaConc / Km))
	rnaCounts = (1 / countsToMolar) * rnaConc
	return (np.log(2) / doublingTime) * rnaCounts + (totalEndoRnaseCountsCapacity * fracSaturated)


def netLossRateFromDilutionAndDegradationRNALinear(doublingTime, degradationRates, rnaCounts):
	return (np.log(2) / doublingTime + degradationRates) * rnaCounts


def expressionFromConditionAndFoldChange(rnaIds, basalExpression, condPerturbations, tfFCs):
	expression = basalExpression.copy()

	rnaIdxs = []
	fcs = []

	for key in sorted(condPerturbations):
		value = condPerturbations[key]
		rnaIdxs.append(np.where(rnaIds == key)[0][0])
		fcs.append(value)

	for key in sorted(tfFCs):
		rnaIdxs.append(np.where(rnaIds == key + "[c]")[0][0])
		fcs.append(tfFCs[key])

	fcs = [fc for (rnaIdx, fc) in sorted(zip(rnaIdxs, fcs), key = lambda pair: pair[0])]
	rnaIdxs = [rnaIdx for (rnaIdx, fc) in sorted(zip(rnaIdxs, fcs), key = lambda pair: pair[0])]

	rnaIdxsBool = np.zeros(len(rnaIds), dtype = np.bool)
	rnaIdxsBool[rnaIdxs] = 1
	fcs = np.array(fcs)
	scaleTheRestBy = (1. - (expression[rnaIdxs] * fcs).sum()) / (1. - (expression[rnaIdxs]).sum())
	expression[rnaIdxsBool] *= fcs
	expression[~rnaIdxsBool] *= scaleTheRestBy

	return expression


def fitTfPromoterKd(sim_data, cellSpecs):
	sim_data.process.equilibrium.ratesRevOrig = sim_data.process.equilibrium.ratesRev.copy()
	cellDensity = sim_data.constants.cellDensity
	rnaIdList = sim_data.process.transcription.rnaData["id"].tolist()

	def alphaGtZero(kdLog10, activeSignalConc, inactiveSignalConc, signalCoeff, activeKSynth, inactiveKSynth, negativeSignal):
		kd = 10**kdLog10
		pPromBoundActive = sim_data.process.transcription_regulation.pPromoterBoundSKd(activeSignalConc, kd, signalCoeff)
		pPromBoundInactive = sim_data.process.transcription_regulation.pPromoterBoundSKd(inactiveSignalConc, kd, signalCoeff)
		if negativeSignal:
			pPromBoundActive = 1. - pPromBoundActive
			pPromBoundInactive = 1. - pPromBoundInactive

		# To have alpha > 0, the following expression must be non-negative
		return -1. * (activeKSynth * pPromBoundInactive - inactiveKSynth * pPromBoundActive)

	def alphaPlusDeltaRGtZero(kdLog10, activeSignalConc, inactiveSignalConc, signalCoeff, activeKSynth, inactiveKSynth, negativeSignal):
		kd = 10**kdLog10
		pPromBoundActive = sim_data.process.transcription_regulation.pPromoterBoundSKd(activeSignalConc, kd, signalCoeff)
		pPromBoundInactive = sim_data.process.transcription_regulation.pPromoterBoundSKd(inactiveSignalConc, kd, signalCoeff)
		if negativeSignal:
			pPromBoundActive = 1. - pPromBoundActive
			pPromBoundInactive = 1. - pPromBoundInactive

		# To have alpha + \delta r > 0, the following expression must be non-negative
		return -1. * (activeKSynth * pPromBoundInactive - inactiveKSynth * pPromBoundActive - (activeKSynth - inactiveKSynth))

	def l1Distance(x, x_init):
		return np.abs(x - x_init)


	for tf in sorted(sim_data.tfToActiveInactiveConds):
		activeKey = tf + "__active"
		inactiveKey = tf + "__inactive"

		boundId = sim_data.process.transcription_regulation.activeToBound[tf]
		negativeSignal = False
		if tf != boundId:
			negativeSignal = True
		kd = sim_data.process.equilibrium.getRevRate(boundId + "[c]") / sim_data.process.equilibrium.getFwdRate(boundId + "[c]")
		tfTargets = sorted(sim_data.tfToFC[tf])
		tfTargetsIdxs = [rnaIdList.index(x + "[c]") for x in tfTargets]

		metabolite = sim_data.process.equilibrium.getMetabolite(boundId + "[c]")
		metaboliteCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(boundId + "[c]")

		activeCellVolume = cellSpecs[activeKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		activeCountsToMolar = 1 / (sim_data.constants.nAvogadro * activeCellVolume)
		activeSignalConc = (activeCountsToMolar * cellSpecs[activeKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol / units.L)
		activeSynthProb = sim_data.process.transcription.rnaSynthProb[activeKey]
		activeSynthProbTargets = activeSynthProb[tfTargetsIdxs]

		inactiveCellVolume = cellSpecs[inactiveKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		inactiveCountsToMolar = 1 / (sim_data.constants.nAvogadro * inactiveCellVolume)
		inactiveSignalConc = (inactiveCountsToMolar * cellSpecs[inactiveKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol / units.L)
		inactiveSynthProb = sim_data.process.transcription.rnaSynthProb[inactiveKey]
		inactiveSynthProbTargets = inactiveSynthProb[tfTargetsIdxs]

		kdLog10Init = np.log10(kd)
		constraints = [
			{"type": "ineq", "fun": lambda logKd, power: (logKd / power) + 12, "args": (metaboliteCoeff,)},
			{"type": "ineq", "fun": lambda logKd, power: -(logKd / power), "args": (metaboliteCoeff,)},
		]
		for activeKSynth, inactiveKSynth in zip(activeSynthProbTargets, inactiveSynthProbTargets):
			args = (activeSignalConc, inactiveSignalConc, metaboliteCoeff, activeKSynth, inactiveKSynth, negativeSignal)
			constraints.append({"type": "ineq", "fun": alphaGtZero, "args": args})
			constraints.append({"type": "ineq", "fun": alphaPlusDeltaRGtZero, "args": args})

		ret = scipy.optimize.minimize(l1Distance, kdLog10Init, args = kdLog10Init, method = "COBYLA", constraints = constraints, options = {"catol": 1e-10})
		if ret.status == 1:
			kdNew = 10**ret.x
			sim_data.process.equilibrium.setRevRate(boundId + "[c]", kdNew * sim_data.process.equilibrium.getFwdRate(boundId + "[c]"))
		else:
			raise Exception, "Can't get positive RNA Polymerase recruitment rate for %s" % tf

def fitPromoterBoundProbability(sim_data, cellSpecs):

	def fromArray(p, pPromoterBound, pPromoterBoundIdxs):
		for condition in sorted(pPromoterBoundIdxs):
			for tf in sorted(pPromoterBoundIdxs[condition]):
				pPromoterBound[condition][tf] = p[pPromoterBoundIdxs[condition][tf]]

	def updateSynthProb(sim_data, kInfo, k):
		for D, value in zip(kInfo, k):
			sim_data.process.transcription.rnaSynthProb[D["condition"]][D["idx"]] = value

		for condition in sim_data.process.transcription.rnaSynthProb:
			assert np.all(sim_data.process.transcription.rnaSynthProb[condition] >= 0)
			sim_data.process.transcription.rnaSynthProb[condition] /= sim_data.process.transcription.rnaSynthProb[condition].sum()

	def makeG(sim_data, pPromoterBound):
		gI, gJ, gV, k, rowNames, colNames, kInfo = [], [], [], [], [], [], []
		for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
			rnaIdNoLoc = rnaId[:-3]

			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			conditions = ["basal"]
			tfsWithData = []
			for tf in tfs:
				if tf not in sorted(sim_data.tfToActiveInactiveConds):
					continue
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)
			for condition in conditions:
				if len(tfsWithData) > 0 and condition == "basal":
					continue
				rowName = rnaIdNoLoc + "__" + condition
				rowNames.append(rowName)
				for tf in tfsWithData:
					colName = rnaIdNoLoc + "__" + tf
					if colName not in colNames:
						colNames.append(colName)
					gI.append(rowNames.index(rowName))
					gJ.append(colNames.index(colName))
					gV.append(pPromoterBound[condition][tf])
				colName = rnaIdNoLoc + "__alpha"
				if colName not in colNames:
					colNames.append(colName)
				gI.append(rowNames.index(rowName))
				gJ.append(colNames.index(colName))
				gV.append(1.)
				k.append(sim_data.process.transcription.rnaSynthProb[condition][idx])
				kInfo.append({"condition": condition, "idx": idx})

		k = np.array(k)
		gI, gJ, gV = np.array(gI), np.array(gJ), np.array(gV)
		G = np.zeros((len(rowNames), len(colNames)), np.float64)
		G[gI, gJ] = gV

		return G, k, rowNames, colNames, kInfo

	def makeZ(sim_data, colNames):
		combinationIdxToColIdxs = {
			0: [0], 1: [0, 1], 2: [0, 2], 3: [0, 1, 2],
			4: [0, 3], 5: [0, 1, 3], 6: [0, 2, 3], 7: [0, 1, 2, 3],
			8: [0, 4], 9: [0, 1, 4], 10: [0, 2, 4], 11: [0, 1, 2, 4],
			12: [0, 3, 4], 13: [0, 1, 3, 4], 14: [0, 2, 3, 4], 15: [0, 1, 2, 3, 4],
			}
		zI, zJ, zV, rowNames = [], [], [], []
		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			tfsWithData = []
			colIdxs = [colNames.index(rnaIdNoLoc + "__alpha")]
			for tf in tfs:
				if tf not in sim_data.tfToActiveInactiveConds:
					continue
				tfsWithData.append(tf)
				colIdxs.append(colNames.index(rnaIdNoLoc + "__" + tf))
			nTfs = len(tfsWithData)
			for combinationIdx in xrange(2**nTfs):
				rowName = rnaIdNoLoc + "__%d" % combinationIdx
				rowNames.append(rowName)
				for colIdx in combinationIdxToColIdxs[combinationIdx]:
					zI.append(rowNames.index(rowName))
					zJ.append(colIdxs[colIdx])
					zV.append(1)

		zI, zJ, zV = np.array(zI), np.array(zJ), np.array(zV)

		Z = np.zeros((zI.max() + 1, zJ.max() + 1), np.float64)
		Z[zI, zJ] = zV

		return Z

	def makeT(sim_data, colNames):
		tI, tJ, tV, rowNamesT = [], [], [], []
		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			tfsWithData = []
			for tf in tfs:
				if tf not in sim_data.tfToActiveInactiveConds:
					continue
				tfsWithData.append(tf)
			for tf in tfsWithData:
				rowName = rnaIdNoLoc + "__" + tf
				rowNamesT.append(rowName)
				colName = rnaIdNoLoc + "__" + tf
				tI.append(rowNamesT.index(rowName))
				tJ.append(colNames.index(colName))
				tV.append(sim_data.tfToDirection[tf][rnaIdNoLoc])
			rowName = rnaIdNoLoc + "__alpha"
			rowNamesT.append(rowName)
			colName = rnaIdNoLoc + "__alpha"
			tI.append(rowNamesT.index(rowName))
			tJ.append(colNames.index(colName))
			tV.append(0)

		tI, tJ, tV = np.array(tI), np.array(tJ), np.array(tV)
		T = np.zeros((tI.max() + 1, tJ.max() + 1), np.float64)
		T[tI, tJ] = tV

		return T

	def makeH(sim_data, colNames, pPromoterBound, r, fixedTFs, cellSpecs):
		rDict = dict([(colName, value) for colName, value in zip(colNames, r)])

		pPromoterBoundIdxs = dict([(condition, {}) for condition in pPromoterBound])
		hI, hJ, hV, rowNames, colNamesH, pInitI, pInitV = [], [], [], [], [], [], []
		for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
			rnaIdNoLoc = rnaId[:-3]
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			conditions = ["basal"]
			tfsWithData = []
			for tf in tfs:
				if tf not in sorted(sim_data.tfToActiveInactiveConds):
					continue
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)
			for condition in conditions:
				if len(tfsWithData) > 0 and condition == "basal":
					continue
				rowName = rnaIdNoLoc + "__" + condition
				rowNames.append(rowName)
				for tf in tfsWithData:
					colName = tf + "__" + condition
					if colName not in colNamesH:
						colNamesH.append(colName)
					hI.append(rowNames.index(rowName))
					hJ.append(colNamesH.index(colName))

					# Handle the case of the TF being knocked out (admittedly not the cleanest solution)
					if cellSpecs[condition]["bulkAverageContainer"].count(tf + "[c]") == 0:
						hV.append(0)
					else:
						hV.append(rDict[rnaIdNoLoc + "__" + tf])
					pInitI.append(colNamesH.index(colName))
					pInitV.append(pPromoterBound[condition][tf])
					pPromoterBoundIdxs[condition][tf] = colNamesH.index(colName)
				colName = rnaIdNoLoc + "__alpha"
				if colName not in colNamesH:
					colNamesH.append(colName)
				hI.append(rowNames.index(rowName))
				hJ.append(colNamesH.index(colName))
				hV.append(rDict[colName])
				pInitI.append(colNamesH.index(colName))
				pInitV.append(1.)

		pInit = np.zeros(len(set(pInitI)))
		pInit[pInitI] = pInitV
		hI, hJ, hV = np.array(hI), np.array(hJ), np.array(hV)
		Hshape = (hI.max() + 1, hJ.max() + 1)
		H = np.zeros(Hshape, np.float64)
		H[hI, hJ] = hV
		pAlphaIdxs = np.array([colNamesH.index(colName) for colName in colNamesH if colName.endswith("__alpha")])
		pNotAlphaIdxs = np.array([colNamesH.index(colName) for colName in colNamesH if not colName.endswith("__alpha")])
		fixedTFIdxs = []
		for idx, colName in enumerate(colNamesH):
			secondElem = colName.split("__")[1]
			if secondElem in fixedTFs:
				fixedTFIdxs.append(idx)
		fixedTFIdxs = np.array(fixedTFIdxs, dtype = np.int)

		return H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, colNamesH

	def makePdiff(sim_data, colNamesH, pPromoterBound):
		PdiffI, PdiffJ, PdiffV = [], [], []
		for rowIdx, tf in enumerate(sorted(sim_data.tfToActiveInactiveConds)):
			condition = tf + "__active"
			colName = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(colNamesH.index(colName))
			PdiffV.append(1)

			condition = tf + "__inactive"
			colName = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(colNamesH.index(colName))
			PdiffV.append(-1)

		PdiffI, PdiffJ, PdiffV = np.array(PdiffI), np.array(PdiffJ), np.array(PdiffV)
		Pdiffshape = (PdiffI.max() + 1, len(colNamesH))
		Pdiff = np.zeros(Pdiffshape, np.float64)
		Pdiff[PdiffI, PdiffJ] = PdiffV

		return Pdiff

	pPromoterBound = calculatePromoterBoundProbability(sim_data, cellSpecs)
	pInit0 = None
	lastNorm = np.inf

	fixedTFs = []
	for tf in sim_data.tfToActiveInactiveConds:
		if sim_data.process.transcription_regulation.tfToTfType[tf] == "2CS":
			fixedTFs.append(tf)
		if sim_data.process.transcription_regulation.tfToTfType[tf] == "1CS" and sim_data.tfToActiveInactiveConds[tf]["active nutrients"] == sim_data.tfToActiveInactiveConds[tf]["inactive nutrients"]:
			fixedTFs.append(tf)

	SCALING = 1e1
	NORM = 1
	for _ in xrange(100):
		G, k, rowNamesG, colNamesG, kInfo = makeG(sim_data, pPromoterBound)
		Z = makeZ(sim_data, colNamesG)
		T = makeT(sim_data, colNamesG)


		R = Variable(G.shape[1])
		prob = Problem(Minimize(norm(G * (SCALING * R) - (SCALING * k), NORM)), [0 <= Z * (SCALING * R), Z * (SCALING * R) <= SCALING * 1, T * (SCALING * R) >= 0])
		prob.solve(solver = "GLPK")
		if prob.status != "optimal":
			raise Exception, "Solver could not find optimal value"
		r = np.array(R.value).reshape(-1)

		print np.linalg.norm(np.dot(G, r) - k, NORM)

		H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, colNamesH = makeH(sim_data, colNamesG, pPromoterBound, r, fixedTFs, cellSpecs)
		Pdiff = makePdiff(sim_data, colNamesH, pPromoterBound)
		if _ == 0:
			pInit0 = pInit.copy()

		print np.linalg.norm(np.dot(H, pInit) - k, NORM)

		P = Variable(H.shape[1])
		D = np.zeros(H.shape[1])
		D[pAlphaIdxs] = 1
		D[fixedTFIdxs] = 1
		Drhs = pInit0.copy()
		Drhs[D != 1] = 0
		prob = Problem(Minimize(norm(H * (SCALING * P) - (SCALING * k), NORM) + 1e-3 * norm(P - pInit0, NORM)), [0 <= (SCALING * P), (SCALING * P) <= SCALING * 1, np.diag(D) * (SCALING * P) == (SCALING * Drhs), Pdiff * (SCALING * P) >= SCALING * 0.1])
		prob.solve(solver = "GLPK")
		if prob.status != "optimal":
			raise Exception, "Solver could not find optimal value"
		pF = np.array(P.value).reshape(-1)
		fromArray(pF, pPromoterBound, pPromoterBoundIdxs)

		print np.linalg.norm(np.dot(H, pF) - k, NORM)

		if np.abs(np.linalg.norm(np.dot(H, pF) - k, NORM) - lastNorm) < 1e-9:
			break
		else:
			lastNorm = np.linalg.norm(np.dot(H, pF) - k, NORM)
	sim_data.pPromoterBound = pPromoterBound
	updateSynthProb(sim_data, kInfo, np.dot(H, pF))

	cellDensity = sim_data.constants.cellDensity
	rnaIdList = sim_data.process.transcription.rnaData["id"].tolist()
	for tf in sorted(sim_data.tfToActiveInactiveConds):
		if sim_data.process.transcription_regulation.tfToTfType[tf] != "1CS":
			continue
		if len(sim_data.tfToActiveInactiveConds[tf]["active genotype perturbations"]) > 0 or len(sim_data.tfToActiveInactiveConds[tf]["inactive genotype perturbations"]) > 0:
			print "Not updating 1CS parameters for %s" % tf
			continue
		activeKey = tf + "__active"
		inactiveKey = tf + "__inactive"

		boundId = sim_data.process.transcription_regulation.activeToBound[tf]
		negativeSignal = False
		if tf != boundId:
			negativeSignal = True
		kd = sim_data.process.equilibrium.getRevRate(boundId + "[c]") / sim_data.process.equilibrium.getFwdRate(boundId + "[c]")
		tfTargets = sorted(sim_data.tfToFC[tf])
		tfTargetsIdxs = [rnaIdList.index(x + "[c]") for x in tfTargets]

		metabolite = sim_data.process.equilibrium.getMetabolite(boundId + "[c]")
		metaboliteCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(boundId + "[c]")

		activeCellVolume = cellSpecs[activeKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		activeCountsToMolar = 1 / (sim_data.constants.nAvogadro * activeCellVolume)
		activeSignalConc = (activeCountsToMolar * cellSpecs[activeKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol / units.L)
		activeSynthProb = sim_data.process.transcription.rnaSynthProb[activeKey]
		activeSynthProbTargets = activeSynthProb[tfTargetsIdxs]

		inactiveCellVolume = cellSpecs[inactiveKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		inactiveCountsToMolar = 1 / (sim_data.constants.nAvogadro * inactiveCellVolume)
		inactiveSignalConc = (inactiveCountsToMolar * cellSpecs[inactiveKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol / units.L)
		inactiveSynthProb = sim_data.process.transcription.rnaSynthProb[inactiveKey]
		inactiveSynthProbTargets = inactiveSynthProb[tfTargetsIdxs]

		kdNew = None
		oldVal = sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite]
		if negativeSignal:
			if 1 - pPromoterBound[activeKey][tf] < 1e-9:
				kdNew = kd
			else:
				kdNew = ((activeSignalConc)**metaboliteCoeff) * pPromoterBound[activeKey][tf] / (1 - pPromoterBound[activeKey][tf])
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (units.mol / units.L) * (kdNew * (1 - pPromoterBound[inactiveKey][tf]) / pPromoterBound[inactiveKey][tf])**(1. / metaboliteCoeff)
		else:
			if pPromoterBound[inactiveKey][tf] < 1e-9:
				kdNew = kd
			else:
				kdNew = ((inactiveSignalConc)**metaboliteCoeff) * (1 - pPromoterBound[inactiveKey][tf]) / pPromoterBound[inactiveKey][tf]
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (units.mol / units.L) * (kdNew * pPromoterBound[activeKey][tf] / (1 - pPromoterBound[activeKey][tf]))**(1. / metaboliteCoeff)
		print metabolite, oldVal, sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite]

		sim_data.process.equilibrium.setRevRate(boundId + "[c]", kdNew * sim_data.process.equilibrium.getFwdRate(boundId + "[c]"))
	return r


def calculatePromoterBoundProbability(sim_data, cellSpecs):
	D = {}
	cellDensity = sim_data.constants.cellDensity
	for conditionKey in sorted(cellSpecs):
		D[conditionKey] = {}

		cellVolume = cellSpecs[conditionKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		for tf in sorted(sim_data.tfToActiveInactiveConds):
			tfType = sim_data.process.transcription_regulation.tfToTfType[tf]

			if tfType == "0CS":
				tfCount = cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")
				if tfCount > 0:
					D[conditionKey][tf] = 1.
				else:
					D[conditionKey][tf] = 0.

			elif tfType == "1CS":
				boundId = sim_data.process.transcription_regulation.activeToBound[tf]
				kd = sim_data.process.equilibrium.getRevRate(boundId + "[c]") / sim_data.process.equilibrium.getFwdRate(boundId + "[c]")
				signal = sim_data.process.equilibrium.getMetabolite(boundId + "[c]")
				signalCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(boundId + "[c]")
				signalConc = (countsToMolar * cellSpecs[conditionKey]["bulkAverageContainer"].count(signal)).asNumber(units.mol / units.L)
				tfConc = (countsToMolar * cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")).asNumber(units.mol / units.L)
				if tf == boundId:
					if tfConc > 0:
						D[conditionKey][tf] = sim_data.process.transcription_regulation.pPromoterBoundSKd(signalConc, kd, signalCoeff)
					else:
						D[conditionKey][tf] = 0.
				else:
					if tfConc > 0:
						D[conditionKey][tf] = 1. - sim_data.process.transcription_regulation.pPromoterBoundSKd(signalConc, kd, signalCoeff)
					else:
						D[conditionKey][tf] = 0.

			elif tfType == "2CS":
				activeTfConc = (countsToMolar * cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")).asNumber(units.mol / units.L)
				inactiveTf = sim_data.process.two_component_system.activeToInactiveTF[tf + "[c]"]
				inactiveTfConc = (countsToMolar * cellSpecs[conditionKey]["bulkAverageContainer"].count(inactiveTf)).asNumber(units.mol / units.L)

				if activeTfConc == 0 and inactiveTfConc == 0:
					D[conditionKey][tf] = 0.
				else:
					D[conditionKey][tf] = activeTfConc / (activeTfConc + inactiveTfConc)

	return D

def cosine_similarity(samples):
	"""
	Finds the cosine similarity between samples.

	samples is a matrix of size (n_samples, sample_size)

	The output is a matrix of size (n_samples, n_samples).

	The cosine similarity is the normalized dot product between two
	vectors.  The name originates from the fact that the normalized dot
	product between two vectors is equal to the cosine of the angle
	formed by the two vectors.
	"""

	magnitudes = np.sqrt(np.sum(np.square(samples), 1))

	normed = samples / magnitudes[:, None]

	return normed.dot(normed.T)

def calculateRnapRecruitment(sim_data, cellSpecs, rVector):
	gI = []
	gJ = []
	gV = []
	k = []
	rowNames = []
	colNames = []
	for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
		rnaIdNoLoc = rnaId[:-3]

		tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
		conditions = ["basal"]
		tfsWithData = []
		for tf in tfs:
			if tf not in sorted(sim_data.tfToActiveInactiveConds):
				continue
			conditions.append(tf + "__active")
			conditions.append(tf + "__inactive")
			tfsWithData.append(tf)
		for condition in conditions:
			if len(tfsWithData) > 0 and condition == "basal":
				continue
			rowName = rnaIdNoLoc + "__" + condition
			rowNames.append(rowName)
			for tf in tfsWithData:
				colName = rnaIdNoLoc + "__" + tf
				if colName not in colNames:
					colNames.append(colName)
				gI.append(rowNames.index(rowName))
				gJ.append(colNames.index(colName))
				gV.append(sim_data.pPromoterBound[condition][tf])
			colName = rnaIdNoLoc + "__alpha"
			if colName not in colNames:
				colNames.append(colName)
			gI.append(rowNames.index(rowName))
			gJ.append(colNames.index(colName))
			gV.append(1.)
			k.append(sim_data.process.transcription.rnaSynthProb[condition][idx])

	gI = np.array(gI)
	gJ = np.array(gJ)
	gV = np.array(gV)
	k = np.array(k)

	shape = (gI.max() + 1, gJ.max() + 1)
	G = np.zeros(shape, np.float64)
	G[gI, gJ] = gV
	S = cosine_similarity(G)
	dupIdxs = np.where((np.tril(S, -1) > 1 - 1e-3).sum(axis = 1))[0]
	uniqueIdxs = [x for x in xrange(G.shape[0]) if x not in dupIdxs]
	G = G[uniqueIdxs]
	k = k[uniqueIdxs]
	rowNames = [rowNames[x] for x in uniqueIdxs]
	r = np.linalg.solve(G, k)
	r = rVector

	# TODO: Delete this once things actually work
	# This is like scaffolding
	sim_data.tfCondToAvgRnapRecruitment = {}
	sim_data.tfCondToAvgRnapRecruitment2 = {}
	for tf in sorted(sim_data.tfToActiveInactiveConds):
		activeCondition = tf + "__active"
		inactiveCondition = tf + "__inactive"
		sim_data.tfCondToAvgRnapRecruitment[activeCondition] = np.ones(G.shape[1])
		sim_data.tfCondToAvgRnapRecruitment[inactiveCondition] = np.ones(G.shape[1])
		sim_data.tfCondToAvgRnapRecruitment2[activeCondition] = np.ones(G.shape[1])
		sim_data.tfCondToAvgRnapRecruitment2[inactiveCondition] = np.ones(G.shape[1])

	for tf in sorted(sim_data.tfToActiveInactiveConds):
		activeCondition = tf + "__active"
		inactiveCondition = tf + "__inactive"

		for idx, colName in enumerate(colNames):
			if not colName.endswith("__alpha"):
				_, thisTf = colName.split("__")
				sim_data.tfCondToAvgRnapRecruitment[activeCondition][idx] = sim_data.pPromoterBound[activeCondition][thisTf]
				sim_data.tfCondToAvgRnapRecruitment[inactiveCondition][idx] = sim_data.pPromoterBound[inactiveCondition][thisTf]
				if tf == thisTf:
					sim_data.tfCondToAvgRnapRecruitment2[activeCondition][idx] = sim_data.pPromoterBound[activeCondition][thisTf]
					sim_data.tfCondToAvgRnapRecruitment2[inactiveCondition][idx] = sim_data.pPromoterBound[inactiveCondition][thisTf]

	# TODO: End delete


	hI = []
	hJ = []
	hV = []
	rowNames = []
	stateMasses = []
	for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
		rnaIdNoLoc = rnaId[:-3]

		tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
		tfsWithData = []
		for tf in tfs:
			if tf not in sorted(sim_data.tfToActiveInactiveConds):
				continue
			tfsWithData.append({"id": tf, "mass_g/mol": sim_data.getter.getMass([tf]).asNumber(units.g / units.mol)})
		rowName = rnaIdNoLoc + "__" + condition
		rowNames.append(rowName)
		for tf in tfsWithData:
			colName = rnaIdNoLoc + "__" + tf["id"]
			hI.append(rowNames.index(rowName))
			hJ.append(colNames.index(colName))
			hV.append(r[colNames.index(colName)])
			stateMasses.append([0.] * 6 + [tf["mass_g/mol"]] + [0.] * 4)
		colName = rnaIdNoLoc + "__alpha"
		hI.append(rowNames.index(rowName))
		hJ.append(colNames.index(colName))
		hV.append(r[colNames.index(colName)])
		stateMasses.append([0.] * 11)

	stateMasses = units.g / units.mol * np.array(stateMasses)
	hI = np.array(hI)
	hJ = np.array(hJ)
	hV = np.array(hV)
	shape = (hI.max() + 1, hJ.max() + 1)
	H = np.zeros(shape, np.float64)
	H[hI, hJ] = hV

	# Deals with numerical tolerance issue of having negative alpha values
	colIdxs = [colNames.index(colName) for colName in colNames if colName.endswith("__alpha")]
	nRows = H[:, colIdxs].shape[0]
	H[range(nRows), colIdxs] -= H[range(nRows), colIdxs].min()
	hV = H[hI, hJ]

	sim_data.internal_state.bulkMolecules.addToBulkState(colNames, stateMasses)
	sim_data.moleculeGroups.bulkMoleculesSetTo1Division = [x for x in colNames if x.endswith("__alpha")]
	sim_data.moleculeGroups.bulkMoleculesBinomialDivision += [x for x in colNames if not x.endswith("__alpha")]
	sim_data.process.transcription_regulation.recruitmentData = {
		"hI": hI,
		"hJ": hJ,
		"hV": hV,
		"shape": shape,
	}
	sim_data.process.transcription_regulation.recruitmentColNames = colNames


def setKmCooperativeEndoRNonLinearRNAdecay(sim_data, bulkContainer):
	cellDensity = sim_data.constants.cellDensity
	cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction
	countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

	degradationRates = sim_data.process.transcription.rnaData["degRate"]
	endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
	kcatEndoRNase = sim_data.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	isRna = np.zeros(len(isMRna))

	endoRnaseRnaIds = sim_data.moleculeGroups.endoRnase_RnaIDs
	isEndoRnase = np.array([x in endoRnaseRnaIds for x in sim_data.process.transcription.rnaData["id"]])

	rnaCounts = bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
	endoCounts = bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)

	rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
	Kmcounts = (( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc).asNumber()
	sim_data.process.rna_decay.KmFirstOrderDecay = Kmcounts

	# Residuals can be written as follows: Res = f(Km) = 0, then Km = g(Km)
	# Compute derivative g(Km) in counts:
	KmQuadratic = 1 / np.power((1 / countsToMolar * Kmcounts).asNumber(), 2)
	denominator = np.power(np.sum(rnaCounts / (1 / countsToMolar * Kmcounts).asNumber()),2)
	numerator = (1 / countsToMolar * totalEndoRnaseCapacity).asNumber() * (denominator - (rnaCounts / (1 / countsToMolar * Kmcounts).asNumber()))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print "Max derivative (counts) = %f" % max(gDerivative)

	# Compute derivative g(Km) in concentrations:
	KmQuadratic = 1 / np.power(Kmcounts, 2)
	denominator = np.power(np.sum(rnaConc.asNumber() / Kmcounts),2)
	numerator = (totalEndoRnaseCapacity).asNumber() * (denominator - (rnaConc.asNumber() / Kmcounts))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print "Max derivative (concentration) = %f" % max(gDerivative)


	# Sensitivity analysis: alpha (regularization term)
	Alphas = []
	if sim_data.constants.SensitivityAnalysisAlpha:
		Alphas = [0.0001, 0.001, 0.01, 0.1, 1, 10]

	for alpha in Alphas:

		if VERBOSE: print 'Alpha = %f' % alpha

		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.kmLossFunction(
				(totalEndoRnaseCapacity).asNumber(units.mol / units.L / units.s),
				(countsToMolar * rnaCounts).asNumber(units.mol / units.L),
				degradationRates.asNumber(1 / units.s),
				isEndoRnase,
				alpha
			)
		KmCooperativeModel = scipy.optimize.fsolve(LossFunction, Kmcounts, fprime = LossFunctionP)
		sim_data.process.rna_decay.SensitivityAnalysisAlphaResidual[alpha] = np.sum(np.abs(R_aux(KmCooperativeModel)))
		sim_data.process.rna_decay.SensitivityAnalysisAlphaRegulariNeg[alpha] = np.sum(np.abs(Rneg(KmCooperativeModel)))

	alpha = 0.5

	# Sensitivity analysis: kcatEndoRNase
	kcatEndo = []
	if sim_data.constants.SensitivityAnalysisKcatEndo:
		kcatEndo = [0.0001, 0.001, 0.01, 0.1, 1, 10]

	for kcat in kcatEndo:

		if VERBOSE: print 'Kcat = %f' % kcat

		totalEndoRNcap = units.sum(endoRNaseConc * kcat)
		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.kmLossFunction(
				(totalEndoRNcap).asNumber(units.mol / units.L),
				(countsToMolar * rnaCounts).asNumber(units.mol / units.L),
				degradationRates.asNumber(1 / units.s),
				isEndoRnase,
				alpha
			)
		KmcountsIni = (( totalEndoRNcap / degradationRates.asNumber() ) - rnaConc).asNumber()
		KmCooperativeModel = scipy.optimize.fsolve(LossFunction, KmcountsIni, fprime = LossFunctionP)
		sim_data.process.rna_decay.SensitivityAnalysisKcat[kcat] = KmCooperativeModel
		sim_data.process.rna_decay.SensitivityAnalysisKcat_ResIni[kcat] = np.sum(np.abs(R_aux(Kmcounts)))
		sim_data.process.rna_decay.SensitivityAnalysisKcat_ResOpt[kcat] = np.sum(np.abs(R_aux(KmCooperativeModel)))


	# Loss function, and derivative
	LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.kmLossFunction(
				(totalEndoRnaseCapacity).asNumber(units.mol / units.L / units.s),
				(countsToMolar * rnaCounts).asNumber(units.mol / units.L),
				degradationRates.asNumber(1 / units.s),
				isEndoRnase,
				alpha
			)

	needToUpdate = False
	fixturesDir = filepath.makedirs(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"fixtures",
			"endo_km"
			)

	if os.path.exists(os.path.join(fixturesDir, "km.cPickle")):
		KmcountsCached = cPickle.load(open(os.path.join(fixturesDir, "km.cPickle"), "rb"))
		if np.sum(np.abs(R_aux(KmcountsCached))) > 1e-15:
			needToUpdate = True
	else:
		needToUpdate = True


	if needToUpdate:
		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		degradationRates = sim_data.process.transcription.rnaData["degRate"]
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)
		Kmcounts = (( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc).asNumber()

		if VERBOSE: print "Running non-linear optimization"
		KmCooperativeModel = scipy.optimize.fsolve(LossFunction, Kmcounts, fprime = LossFunctionP)
		cPickle.dump(KmCooperativeModel, open(os.path.join(fixturesDir, "km.cPickle"), "w"))
	else:
		if VERBOSE: print "Not running non-linear optimization--using cached result"
		KmCooperativeModel = KmcountsCached

	if VERBOSE > 1:
		print "Loss function (Km inital) = %f" % np.sum(np.abs(LossFunction(Kmcounts)))
		print "Loss function (optimized Km) = %f" % np.sum(np.abs(LossFunction(KmCooperativeModel)))

		print "Negative km ratio = %f" % np.sum(np.abs(Rneg(KmCooperativeModel)))

		print "Residuals (Km initial) = %f" % np.sum(np.abs(R(Kmcounts)))
		print "Residuals optimized = %f" % np.sum(np.abs(R(KmCooperativeModel)))

		print "EndoR residuals (Km initial) = %f" % np.sum(np.abs(isEndoRnase * R(Kmcounts)))
		print "EndoR residuals optimized = %f" % np.sum(np.abs(isEndoRnase * R(KmCooperativeModel)))

		print "Residuals (scaled by Kdeg * RNAcounts) Km initial = %f" % np.sum(np.abs(R_aux(Kmcounts)))
		print "Residuals (scaled by Kdeg * RNAcounts) optimized = %f" % np.sum(np.abs(R_aux(KmCooperativeModel)))


	# Evaluate Jacobian around solutions (Kmcounts and KmCooperativeModel)
	JacobDiag = np.diag(Jacob(KmCooperativeModel))
	Jacob_auxDiag = np.diag(Jacob_aux(KmCooperativeModel))

	# Compute convergence of non-linear optimization: g'(Km)
	Gkm = np.abs(1. - JacobDiag)
	Gkm_aux = np.abs(1. - Jacob_auxDiag)
	sim_data.process.rna_decay.KmConvergence = Gkm_aux

	# Convergence is guaranteed if g'(Km) <= K < 1
	if VERBOSE: print "Convergence (Jacobian) = %.0f%% (<K> = %.5f)" % (len(Gkm[Gkm < 1.]) / float(len(Gkm)) * 100., np.mean(Gkm))
	if VERBOSE: print "Convergence (Jacobian_aux) = %.0f%% (<K> = %.5f)" % (len(Gkm_aux[Gkm_aux < 1.]) / float(len(Gkm_aux)) * 100., np.mean(Gkm_aux[Gkm_aux < 1.]))

	# Save statistics KM optimization
	sim_data.process.rna_decay.StatsFit['LossKm'] = np.sum(np.abs(LossFunction(Kmcounts)))
	sim_data.process.rna_decay.StatsFit['LossKmOpt'] = np.sum(np.abs(LossFunction(KmCooperativeModel)))

	sim_data.process.rna_decay.StatsFit['RnegKmOpt'] = np.sum(np.abs(Rneg(KmCooperativeModel)))

	sim_data.process.rna_decay.StatsFit['ResKm'] = np.sum(np.abs(R(Kmcounts)))
	sim_data.process.rna_decay.StatsFit['ResKmOpt'] = np.sum(np.abs(R(KmCooperativeModel)))

	sim_data.process.rna_decay.StatsFit['ResEndoRNKm'] = np.sum(np.abs(isEndoRnase * R(Kmcounts)))
	sim_data.process.rna_decay.StatsFit['ResEndoRNKmOpt'] = np.sum(np.abs(isEndoRnase * R(KmCooperativeModel)))

	sim_data.process.rna_decay.StatsFit['ResScaledKm'] = np.sum(np.abs(R_aux(Kmcounts)))
	sim_data.process.rna_decay.StatsFit['ResScaledKmOpt'] = np.sum(np.abs(R_aux(KmCooperativeModel)))

	return units.mol / units.L * KmCooperativeModel
