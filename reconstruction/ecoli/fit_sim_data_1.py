#!/usr/bin/env python

from __future__ import division

import numpy as np
import os
import scipy.optimize
import time
import cPickle

import wholecell
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.utils import units
from wholecell.utils.fitting import normalize, massesAndCountsToAddForPools
from wholecell.utils.modular_fba import FluxBalanceAnalysis

import cvxpy

# Hacks
RNA_POLY_MRNA_DEG_RATE_PER_S = np.log(2) / 30. # half-life of 30 seconds
FRACTION_INCREASE_RIBOSOMAL_PROTEINS = 0.2  # reduce stochasticity from protein expression
FRACTION_INCREASE_RNAP_PROTEINS = 0.05

NUMERICAL_ZERO = 1e-10

# TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)

FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100
N_SEEDS = 10

BASAL_EXPRESSION_CONDITION = "M9 Glucose minus AAs"

VERBOSE = False

from models.ecoli.processes.metabolism import SECRETION_PENALTY_COEFF

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

def fitSimData_1(raw_data):

	sim_data = SimulationDataEcoli()
	sim_data.initialize(
		raw_data = raw_data,
		basal_expression_condition = BASAL_EXPRESSION_CONDITION,
		)

	# Increase RNA poly mRNA deg rates
	setRnaPolymeraseCodingRnaDegradationRates(sim_data)

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
	
	buildTfConditionCellSpecifications(sim_data, cellSpecs)
	buildCombinedConditionCellSpecifications(sim_data, cellSpecs)

	sim_data.conditionToAvgCellDryMassInit = dict([(cellSpec, cellSpecs[cellSpec]["avgCellDryMassInit"]) for cellSpec in cellSpecs])

	# Fit kinetic parameters
	findKineticCoeffs(sim_data, cellSpecs["basal"]["bulkContainer"])

	for condition in sorted(cellSpecs):
		spec = cellSpecs[condition]
		bulkAverageContainer, bulkDeviationContainer = calculateBulkDistributions(
			sim_data,
			spec["expression"],
			spec["concDict"],
			spec["avgCellDryMassInit"],
			spec["doubling_time"],
			)
		spec["bulkAverageContainer"] = bulkAverageContainer
		spec["bulkDeviationContainer"] = bulkDeviationContainer

	fitTfPromoterKd(sim_data, cellSpecs)

	sim_data.pPromoterBound = calculatePromoterBoundProbability(sim_data, cellSpecs)

	calculateRnapRecruitment(sim_data, cellSpecs)

	return sim_data


def buildBasalCellSpecifications(sim_data):
	cellSpecs = {}
	cellSpecs["basal"] = {
		"concDict": sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					"minimal"
				),
		"expression": sim_data.process.transcription.rnaExpression["basal"].copy(),
		"doubling_time": sim_data.conditionToDoublingTime["basal"],
	}

	expression, synthProb, avgCellDryMassInit, fitAvgSolublePoolMass, bulkContainer, _ = expressionConverge(
		sim_data,
		cellSpecs["basal"]["expression"],
		cellSpecs["basal"]["concDict"],
		cellSpecs["basal"]["doubling_time"],
		)

	cellSpecs["basal"]["expression"] = expression
	cellSpecs["basal"]["synthProb"] = synthProb
	cellSpecs["basal"]["avgCellDryMassInit"] = avgCellDryMassInit
	cellSpecs["basal"]["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass
	cellSpecs["basal"]["bulkContainer"] = bulkContainer

	sim_data.mass.avgCellDryMassInit = avgCellDryMassInit
	sim_data.mass.avgCellDryMass = sim_data.mass.avgCellDryMassInit * sim_data.mass.avgCellToInitialCellConvFactor
	sim_data.mass.avgCellWaterMassInit = sim_data.mass.avgCellDryMassInit / sim_data.mass.cellDryMassFraction * sim_data.mass.cellWaterMassFraction
	sim_data.mass.fitAvgSolublePoolMass = fitAvgSolublePoolMass


	sim_data.process.transcription.rnaExpression["basal"][:] = cellSpecs["basal"]["expression"]
	sim_data.process.transcription.rnaSynthProb["basal"][:] = cellSpecs["basal"]["synthProb"]

	return cellSpecs


def buildTfConditionCellSpecifications(sim_data, cellSpecs):
	for tf in sorted(sim_data.tfToActiveInactiveConds):
		for choice in ["__active", "__inactive"]:
			conditionKey = tf + choice
			conditionValue = sim_data.conditions[conditionKey]

			fcData = {}
			if choice == "__active":
				fcData = sim_data.tfToFC[tf]
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

			expression, synthProb, avgCellDryMassInit, fitAvgSolublePoolMass, bulkContainer, concDict = expressionConverge(
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
			cellSpecs[conditionKey]["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass
			cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

			sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
			sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

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

		expression, synthProb, avgCellDryMassInit, fitAvgSolublePoolMass, bulkContainer, concDict = expressionConverge(
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
		cellSpecs[conditionKey]["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

		# Uncomment when concDict is actually calculated for non-base [AA]
		# if len(conditionValue["perturbations"]) == 0:
		# 	nutrientLabel = conditionValue["nutrients"]
		# 	sim_data.process.metabolism.nutrientsToInternalConc[nutrientLabel] = concDict


def expressionConverge(sim_data, expression, concDict, doubling_time, Km = None, updateConcDict = False):
	# Fit synthesis probabilities for RNA
	for iteration in xrange(MAX_FITTING_ITERATIONS):
		if VERBOSE: print 'Iteration: {}'.format(iteration)

		initialExpression = expression.copy()

		expression = setInitialRnaExpression(sim_data, expression, doubling_time)

		bulkContainer = createBulkContainer(sim_data, expression, doubling_time)

		avgCellDryMassInit, fitAvgSolublePoolMass = rescaleMassForSolubleMetabolites(sim_data, bulkContainer, concDict, doubling_time)

		setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time)

		setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		# Normalize expression and write out changes

		expression, synthProb = fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		if updateConcDict:
			concDict = concDict.copy() # Calculate non-base condition [AA]

		finalExpression = expression

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - finalExpression)))
		if VERBOSE: print 'degree of fit: {}'.format(degreeOfFit)

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	return expression, synthProb, avgCellDryMassInit, fitAvgSolublePoolMass, bulkContainer, concDict


# Sub-fitting functions

def setRnaPolymeraseCodingRnaDegradationRates(sim_data):
	# Increase RNA poly mRNA deg rates
	# TODO: set this based on transcription unit structure
	# i.e. same synthesis prob. but different deg rates

	rnaPolySubunits = sim_data.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"]
	subunitIndexes = np.array([np.where(sim_data.process.translation.monomerData["id"] == id_)[0].item() for id_ in rnaPolySubunits]) # there has to be a better way...
	mRNA_indexes = sim_data.relation.rnaIndexToMonomerMapping[subunitIndexes]
	sim_data.process.transcription.rnaData.struct_array["degRate"][mRNA_indexes] = RNA_POLY_MRNA_DEG_RATE_PER_S

def setCPeriod(sim_data):
	sim_data.growthRateParameters.c_period = sim_data.process.replication.genome_length * units.nt / sim_data.growthRateParameters.dnaPolymeraseElongationRate / 2

def rescaleMassForSolubleMetabolites(sim_data, bulkMolCntr, concDict, doubling_time):
	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)

	mass = (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = sorted(concDict)
	poolConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in poolIds])

	massesToAdd, countsToAdd = massesAndCountsToAddForPools(
		mass,
		poolIds,
		poolConcentrations,
		sim_data.getter.getMass(poolIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculePoolsDryMass = units.hstack((massesToAdd[:poolIds.index('WATER[c]')], massesToAdd[poolIds.index('WATER[c]') + 1:]))
	newAvgCellDryMassInit = units.sum(mass) + units.sum(smallMoleculePoolsDryMass)
	fitAvgSolublePoolMass = units.sum(units.hstack((massesToAdd[:poolIds.index('WATER[c]')], massesToAdd[poolIds.index('WATER[c]') + 1:]))) * sim_data.mass.avgCellToInitialCellConvFactor

	return newAvgCellDryMassInit, fitAvgSolublePoolMass

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
	distribution_rRNA23S = np.array([1.] + [0.] * (ids_rRNA23S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA16S = np.array([1.] + [0.] * (ids_rRNA16S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA5S = np.array([1.] + [0.] * (ids_rRNA5S.size-1)) # currently only expressing first rRNA operon
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
	ids_molecules = sim_data.state.bulkMolecules.bulkData["id"]

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
	if VERBOSE: print '30S limit: {}'.format(constraint_names[np.where(rib30lims.max() == rib30lims)[0]][-1])
	if VERBOSE: print '30S actual count: {}'.format((ribosome30SCounts / ribosome30SStoich).min())
	if VERBOSE: print '30S count set to: {}'.format(rib30lims[np.where(rib30lims.max() == rib30lims)[0]][-1])
	if VERBOSE: print '50S limit: {}'.format(constraint_names[np.where(rib50lims.max() == rib50lims)[0]][-1])
	if VERBOSE: print '50S actual count: {}'.format((ribosome50SCounts / ribosome50SStoich).min())
	if VERBOSE: print '50S count set to: {}'.format(rib50lims[np.where(rib50lims.max() == rib50lims)[0]][-1])

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
		) * (1 + FRACTION_INCREASE_RNAP_PROTEINS)

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(rnapIds)

	## -- SET RNAP COUNTS TO MAXIMIM CONSTRAINTS -- #
	constraint_names = np.array(["Current level OK", "Insufficient to double RNA distribution"])
	rnapLims = np.array([(rnapCounts / rnapStoich).min(), (minRnapSubunitCounts / rnapStoich).min()])
	if VERBOSE: print 'rnap limit: {}'.format(constraint_names[np.where(rnapLims.max() == rnapLims)[0]][0])
	if VERBOSE: print 'rnap actual count: {}'.format((rnapCounts / rnapStoich).min())
	if VERBOSE: print 'rnap counts set to: {}'.format(rnapLims[np.where(rnapLims.max() == rnapLims)[0]][0])

	bulkContainer.countsIs(
		np.fmax(rnapCounts, minRnapSubunitCounts),
		rnapIds
		)


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
	ids_metabolites = sorted(concDict)
	conc_metabolites = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in ids_metabolites])
	allMoleculesIDs = sorted(
		set(ids_rnas) | set(ids_protein) | set(ids_complex) | set(ids_equilibrium) | set(ids_metabolites)
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

	bulkContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'])
	rnaView = bulkContainer.countsView(ids_rnas)
	proteinView = bulkContainer.countsView(ids_protein)
	complexationMoleculesView = bulkContainer.countsView(ids_complex)
	equilibriumMoleculesView = bulkContainer.countsView(ids_equilibrium)
	metabolitesView = bulkContainer.countsView(ids_metabolites)
	allMoleculesView = bulkContainer.countsView(allMoleculesIDs)

	allMoleculeCounts = np.empty((N_SEEDS, allMoleculesView.counts().size), np.int64)


	for seed in xrange(N_SEEDS):
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
			metDiffs = metabolitesView.counts() - metCounts.asNumber().round()

			nIters += 1
			if nIters > 100:
				raise Exception, "Equilibrium reactions are not converging!"

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	bulkAverageContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'], np.float64)

	bulkAverageContainer.countsIs(allMoleculeCounts.mean(0), allMoleculesIDs)
	bulkDeviationContainer.countsIs(allMoleculeCounts.std(0), allMoleculesIDs)

	return bulkAverageContainer, bulkDeviationContainer


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
	fracSaturated = rnaConc / (Km + rnaConc)
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

	def alphaGtZero(kdLog10, activeSignalConc, inactiveSignalConc, signalCoeff, activeKSynth, inactiveKSynth):
		kd = 10**kdLog10
		pPromBoundActive = sim_data.process.transcription_regulation.pPromoterBoundSKd(activeSignalConc, kd, signalCoeff)
		pPromBoundInactive = sim_data.process.transcription_regulation.pPromoterBoundSKd(inactiveSignalConc, kd, signalCoeff)

		# To have alpha > 0, the following expression must be non-negative
		return -1. * (activeKSynth * pPromBoundInactive - inactiveKSynth * pPromBoundActive)

	def alphaPlusDeltaRGtZero(kdLog10, activeSignalConc, inactiveSignalConc, signalCoeff, activeKSynth, inactiveKSynth):
		kd = 10**kdLog10
		pPromBoundActive = sim_data.process.transcription_regulation.pPromoterBoundSKd(activeSignalConc, kd, signalCoeff)
		pPromBoundInactive = sim_data.process.transcription_regulation.pPromoterBoundSKd(inactiveSignalConc, kd, signalCoeff)

		# To have alpha + \delta r > 0, the following expression must be non-negative
		return -1. * (activeKSynth * pPromBoundInactive - inactiveKSynth * pPromBoundActive - (activeKSynth - inactiveKSynth))

	def l1Distance(x, x_init):
		return np.abs(x - x_init)


	for tf in sorted(sim_data.tfToActiveInactiveConds):
		activeKey = tf + "__active"
		inactiveKey = tf + "__inactive"

		kd = sim_data.process.equilibrium.getRevRate(tf + "[c]") / sim_data.process.equilibrium.getFwdRate(tf + "[c]")
		tfTargets = sorted(sim_data.tfToFC[tf])
		tfTargetsIdxs = [rnaIdList.index(x + "[c]") for x in tfTargets]

		metabolite = sim_data.process.equilibrium.getMetabolite(tf + "[c]")
		metaboliteCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(tf + "[c]")

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
			args = (activeSignalConc, inactiveSignalConc, metaboliteCoeff, activeKSynth, inactiveKSynth)
			constraints.append({"type": "ineq", "fun": alphaGtZero, "args": args})
			constraints.append({"type": "ineq", "fun": alphaPlusDeltaRGtZero, "args": args})

		ret = scipy.optimize.minimize(l1Distance, kdLog10Init, args = kdLog10Init, method = "COBYLA", constraints = constraints, options = {"catol": 1e-9})
		if ret.status == 1:
			kdNew = 10**ret.x
			sim_data.process.equilibrium.setRevRate(tf + "[c]", kdNew * sim_data.process.equilibrium.getFwdRate(tf + "[c]"))
		else:
			raise Exception, "Can't get positive RNA Polymerase recruitment rate for %s" % tf


def calculatePromoterBoundProbability(sim_data, cellSpecs):
	D = {}
	cellDensity = sim_data.constants.cellDensity
	for conditionKey in sorted(cellSpecs):
		D[conditionKey] = {}

		cellVolume = cellSpecs[conditionKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		for tf in sorted(sim_data.tfToActiveInactiveConds):
			kd = sim_data.process.equilibrium.getRevRate(tf + "[c]") / sim_data.process.equilibrium.getFwdRate(tf + "[c]")
			signal = sim_data.process.equilibrium.getMetabolite(tf + "[c]")
			signalCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(tf + "[c]")
			signalConc = (countsToMolar * cellSpecs[conditionKey]["bulkAverageContainer"].count(signal)).asNumber(units.mol / units.L)

			D[conditionKey][tf] = sim_data.process.transcription_regulation.pPromoterBoundSKd(signalConc, kd, signalCoeff)
	return D


def calculateRnapRecruitment(sim_data, cellSpecs):
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
	r = np.linalg.solve(G, k)

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

	sim_data.state.bulkMolecules.addToBulkState(colNames, stateMasses)
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

	# Loss function, and derivative 
	LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux = sim_data.process.rna_decay.kmLossFunction(
				(totalEndoRnaseCapacity).asNumber(units.mol / units.L / units.s),
				(countsToMolar * rnaCounts).asNumber(units.mol / units.L),
				degradationRates.asNumber(1 / units.s),
				isEndoRnase
			)

	needToUpdate = False
	fixturesDir = os.path.join(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"fixtures",
			"endo_km"
			)

	if not os.path.exists(fixturesDir):
		needToUpdate = True
		os.makedirs(fixturesDir)

	if os.path.exists(os.path.join(fixturesDir, "km.cPickle")):
		Kmcounts = cPickle.load(open(os.path.join(fixturesDir, "km.cPickle"), "rb"))
		if np.sum(np.abs(R_aux(Kmcounts))) > 1e-15:
			needToUpdate = True
	else:
		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		degradationRates = sim_data.process.transcription.rnaData["degRate"]
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)
		Kmcounts = (( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc).asNumber()
		needToUpdate = True

	if needToUpdate:
		if VERBOSE: print "Running non-linear optimization"
		KmCooperativeModel = scipy.optimize.fsolve(LossFunction, Kmcounts, fprime = LossFunctionP)
		cPickle.dump(KmCooperativeModel, open(os.path.join(fixturesDir, "km.cPickle"), "w"))
	else:
		if VERBOSE: print "Not running non-linear optimization--using cached result"
		KmCooperativeModel = Kmcounts

	if VERBOSE:
		print "Loss function (Km inital) = %f" % np.sum(np.abs(LossFunction(Kmcounts)))
		print "Loss function (optimized Km) = %f" % np.sum(np.abs(LossFunction(KmCooperativeModel)))

		print "Negative km ratio = %f" % np.sum(np.abs(Rneg(KmCooperativeModel)))

		print "Residuals optimized = %f" % np.sum(np.abs(R(KmCooperativeModel)))
		print "Residuals (Km initial) = %f" % np.sum(np.abs(R(Kmcounts)))

		print "EndoR residuals optimized = %f" % np.sum(np.abs(isEndoRnase * R(Kmcounts)))
		print "EndoR residuals optimized = %f" % np.sum(np.abs(isEndoRnase * R(KmCooperativeModel)))

		print "Residuals (scaled by RNAcounts) Km initial = %f" % np.sum(np.abs(R_aux(Kmcounts)))
		print "Residuals (scaled by RNAcounts) optimized = %f" % np.sum(np.abs(R_aux(KmCooperativeModel)))

	return units.mol / units.L * KmCooperativeModel


def findKineticCoeffs(sim_data, bulkContainer):

	# Estimate the volume of the initial cell
	initDryMass = sim_data.mass.avgCellDryMassInit
	dryMassPerCellMass = (1. - sim_data.mass.cellWaterMassFraction)
	totalCellMassInit = initDryMass / dryMassPerCellMass
	cellVolumeInit = totalCellMassInit / sim_data.constants.cellDensity

	energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / totalCellMassInit

	reactionStoich = sim_data.process.metabolism.reactionStoich
	externalExchangeMolecules = sorted(sim_data.nutrientData["externalExchangeMolecules"]["minimal"])
	extMoleculeMasses = sim_data.getter.getMass(externalExchangeMolecules)

	moleculeMasses = dict(zip(
		externalExchangeMolecules,
		sim_data.getter.getMass(externalExchangeMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)
		))

	# Make previously observed external molecule fluxes into coefficients of a biomass reaction
	jfbaBiomassReactionStoich = {molID:-coeff.asNumber(COUNTS_UNITS/VOLUME_UNITS/TIME_UNITS) for molID, coeff in sim_data.process.metabolism.previousBiomassMeans.iteritems() if np.abs(coeff.asNumber(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS)) > NUMERICAL_ZERO}
	jfbaBiomassReactionStoich["biomass"] = 1

	n_iterations = 1
	overallFluxes = None
	modification_probability = .2

	for iteration in xrange(n_iterations):

		reactionStoichWithBiomass = reactionStoich.copy()

		jfbaBiomassReactionStoichMod = {}
		for key, value in jfbaBiomassReactionStoich.iteritems():
			changeAmount = 1.
			if iteration > 0:
				changeBoolean = np.random.random() < modification_probability
				changeAmount = np.random.random()*2 if changeBoolean else 1.
			jfbaBiomassReactionStoichMod[key] = value * changeAmount
		reactionStoichWithBiomass["JFBA-BIOMASS-RXN"] = jfbaBiomassReactionStoichMod

		fbaObject = FluxBalanceAnalysis(
			reactionStoich = reactionStoichWithBiomass,
			externalExchangedMolecules = externalExchangeMolecules,
			objective = {"biomass": 1},
			objectiveType = "standard",
			moleculeMasses=moleculeMasses,
			secretionPenaltyCoeff = SECRETION_PENALTY_COEFF,
			solver = "glpk",
			maintenanceCostGAM = energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			maintenanceReaction = {
				"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1,
				}
			)

		# Implicit one second time-step
		coefficient = dryMassPerCellMass * sim_data.constants.cellDensity * (1. * units.s)

		externalMoleculeLevels, newObjective = sim_data.process.metabolism.exchangeConstraints(
			fbaObject.externalMoleculeIDs(),
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			sim_data.nutrientsTimeSeriesLabel,
			1., # time only matters in changing environments, so any > 0 is fine
			preview=True
			)

		fbaObject.externalMoleculeLevelsIs(externalMoleculeLevels)

		fbaObject.maxReactionFluxIs(fbaObject._reactionID_NGAM, (sim_data.constants.nonGrowthAssociatedMaintenance * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
		fbaObject.minReactionFluxIs(fbaObject._reactionID_NGAM, (sim_data.constants.nonGrowthAssociatedMaintenance * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		arrayModel = fbaObject.getArrayBasedModel()

		S_matrix = arrayModel["S_matrix"]
		reactionNames = arrayModel["Reactions"]
		upperBoundsDict = arrayModel["Upper bounds"]
		lowerBoundsDict = arrayModel["Lower bounds"]
		upperBounds = np.array([upperBoundsDict[x] for x in reactionNames])
		lowerBounds = np.array([lowerBoundsDict[x] for x in reactionNames])

		x = cvxpy.Variable(S_matrix.shape[1])
		c = np.zeros(len(reactionNames))
		c[reactionNames.index("JFBA-BIOMASS-RXN")] = 1

		# Construct the problem.
		objective = cvxpy.Maximize(c*x)
		constraints = [lowerBounds <= x, x <= upperBounds, S_matrix*x == 0]
		prob = cvxpy.Problem(objective, constraints)

		expectedFlux = prob.solve(solver="GLPK")

		predictedFluxes = np.array(x.value)

		if overallFluxes is None:
			overallFluxes = predictedFluxes
			originalPredictedFluxes = predictedFluxes.copy()
		else:
			overallFluxes += predictedFluxes

	overallFluxes /= n_iterations

	fluxNames = np.array(reactionNames).reshape(-1, 1)
	predictedFluxesDict = dict(zip([fluxName[0] for fluxName in fluxNames], [(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * flux[0] for flux in overallFluxes]))

	# Use rabinowitz metabolite concentrations as the estimate
	metaboliteConcentrationsDict = sim_data.process.metabolism.concDict

	# Use fit estimates for protein counts
	proteinIds = sim_data.process.translation.monomerData["id"]
	proteinCounts =  bulkContainer.counts(proteinIds)

	# Complex monomers
	seed = 1 # Any seed is fine
	complexationMoleculeNames = sim_data.process.complexation.moleculeNames
	complexationMolecules = bulkContainer.counts(complexationMoleculeNames).astype(np.int64)
	complexationStoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64, order = "F")
	complexationPrebuiltMatrices = mccBuildMatrices(complexationStoichMatrix)
	updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
		complexationMolecules,
		seed,
		complexationStoichMatrix,
		*complexationPrebuiltMatrices
		)

	# Add updates to concentrations based on complexation
	proteinCountsDict = {}
	for idx, proteinId in enumerate(proteinIds):
		if proteinId in complexationMolecules:
			proteinCountsDict[proteinId] = updatedMoleculeCounts[complexationMolecules.index(proteinId)]
		else:
			proteinCountsDict[proteinId] = proteinCounts[idx]

	proteinConcentrations = (1. / sim_data.constants.nAvogadro / cellVolumeInit) * proteinCounts

	proteinConcDict = {}
	for proteinId, count in proteinCountsDict.iteritems():
		proteinConcDict[proteinId] = ((1. / sim_data.constants.nAvogadro / cellVolumeInit) * count)

	# Build an enzyme-reaction association matrix
	reactionEnzymes  = sim_data.process.metabolism.reactionEnzymes
	reactionIDs = list(fbaObject.reactionIDs())
	reactionIDs.remove("JFBA-BIOMASS-RXN")
	enzymeReactionMatrix = sim_data.process.metabolism.enzymeReactionMatrix(reactionIDs, proteinIds, reactionEnzymes)

	spontaneousIndices = np.where(np.sum(enzymeReactionMatrix, axis=1) == 0)

	# Estimate kinetic rates based on the estimated protein concentrations and the known kcats
	kineticRates = enzymeReactionMatrix.dot(proteinConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))
	kineticRates[spontaneousIndices] = np.inf
	kineticRatesDict = dict(zip(list(reactionIDs), kineticRates))

	# Determine which kinetic constraints overconstrain the predicted metabolic fluxes
	overconstraintRatio = {}
	for rateName, rate in kineticRatesDict.iteritems():
		if rateName in predictedFluxesDict:
			flux = predictedFluxesDict[rateName].asNumber(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS)
			ratio = 0 if flux == 0 else flux / rate
			if ratio >= 1:
				print ratio
			overconstraintRatio[rateName] = ratio

	# Adjust any kcats which are predicted to overconstrain until they no longer do
	coefficientsSet = set()
	for constraintID, reactionInfo in sim_data.process.metabolism.reactionRateInfo.iteritems():
		reactionID = reactionInfo["reactionID"]
		if reactionID in overconstraintRatio:
			if overconstraintRatio[reactionID] > 1:
				coefficientsSet.add((reactionInfo["reactionID"], overconstraintRatio[reactionID]))
				reactionInfo["constraintMultiple"] = np.ceil(overconstraintRatio[reactionID])

	if VERBOSE: print coefficientsSet

	sim_data.process.metabolism.predictedFluxesDict = predictedFluxesDict
