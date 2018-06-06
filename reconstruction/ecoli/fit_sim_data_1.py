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

	for condition in sorted(cellSpecs):
		if sim_data.conditions[condition]["nutrients"] not in sim_data.translationSupplyRate.keys():
			sim_data.translationSupplyRate[sim_data.conditions[condition]["nutrients"]] = cellSpecs[condition]["translation_aa_supply"]

	rVector = fitPromoterBoundProbability(sim_data, cellSpecs)

	for condition in sorted(cellSpecs):
		if VERBOSE > 0:
			print "Updating mass in condition {}".format(condition)
		spec = cellSpecs[condition]

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			sim_data.conditions[condition]["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[condition]))

		avgCellDryMassInit, fitAvgSolublePoolMass = rescaleMassForSolubleMetabolites(sim_data, spec["bulkContainer"], concDict, sim_data.conditionToDoublingTime[condition])

		if VERBOSE > 0:
			print str(spec["avgCellDryMassInit"]) + " to "  + str(avgCellDryMassInit)

		spec["avgCellDryMassInit"] = avgCellDryMassInit
		spec["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass

		mRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition][sim_data.process.transcription.rnaData["isMRna"]].sum()
		tRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition][sim_data.process.transcription.rnaData["isTRna"]].sum()
		rRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition][sim_data.process.transcription.rnaData["isRRna"]].sum()

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.transcription.rnaSynthProbFraction and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.transcription.rnaSynthProbFraction[sim_data.conditions[condition]["nutrients"]] = {
				"mRna": mRnaSynthProb,
				"tRna": tRnaSynthProb,
				"rRna": rRnaSynthProb,
				}

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.transcription.rnaSynthProbRProtein and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.transcription.rnaSynthProbRProtein[sim_data.conditions[condition]["nutrients"]] = sim_data.process.transcription.rnaSynthProb[condition][sim_data.process.transcription.rnaData["isRProtein"]]

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.transcription.rnaSynthProbRnaPolymerase and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.transcription.rnaSynthProbRnaPolymerase[sim_data.conditions[condition]["nutrients"]] = sim_data.process.transcription.rnaSynthProb[condition][sim_data.process.transcription.rnaData["isRnap"]]

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.transcription.rnapFractionActiveDict and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.transcription.rnapFractionActiveDict[sim_data.conditions[condition]["nutrients"]] = sim_data.growthRateParameters.getFractionActiveRnap(spec["doubling_time"])

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.transcription.rnaPolymeraseElongationRateDict and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.transcription.rnaPolymeraseElongationRateDict[sim_data.conditions[condition]["nutrients"]] = sim_data.growthRateParameters.getRnapElongationRate(spec["doubling_time"])

		if sim_data.conditions[condition]["nutrients"] not in sim_data.expectedDryMassIncreaseDict and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.expectedDryMassIncreaseDict[sim_data.conditions[condition]["nutrients"]] = spec["avgCellDryMassInit"]

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.translation.ribosomeElongationRateDict and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.translation.ribosomeElongationRateDict[sim_data.conditions[condition]["nutrients"]] = sim_data.growthRateParameters.getRibosomeElongationRate(spec["doubling_time"])

		if sim_data.conditions[condition]["nutrients"] not in sim_data.process.translation.ribosomeFractionActiveDict and len(sim_data.conditions[condition]["perturbations"]) == 0:
			sim_data.process.translation.ribosomeFractionActiveDict[sim_data.conditions[condition]["nutrients"]] = sim_data.growthRateParameters.getFractionActiveRibosome(spec["doubling_time"])

	calculateRnapRecruitment(sim_data, cellSpecs, rVector)

	return sim_data


def buildBasalCellSpecifications(sim_data):
	"""
	Creates cell specifications for the basal condition by fitting expression.
	Relies on expressionConverge() to set the expression and update masses.

	Requires
	--------
	- Metabolite concentrations based on 'minimal' nutrients
	- 'basal' RNA expression
	- 'basal' doubling time

	Modifies
	--------
	- Average mass values of the cell
	- RNA expression and synthesis probabilities

	Returns
	--------
	- dict {'basal': dict} with the following keys in the dict from key 'basal':
		'concDict' {metabolite_name (str): concentration (float with units)} -
			dictionary of concentrations for each metabolite with a concentration
		'expression' (1D numpy array of floats) - expression for each RNA,
			total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (1D numpy array of floats) - synthesis probability for
			each RNA, total normalized to 1
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - TODO - what is this?
		'bulkContainer' (BulkObjectsContainer object) - expected counts for
			bulk molecules based on expression

	Notes
	-----
	- TODO - bad form to return values and set sim_data values within the function -
	should this be changed?
	"""

	# Create dictionary for basal condition
	cellSpecs = {}
	cellSpecs["basal"] = {
		"concDict": sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients("minimal"),
		"expression": sim_data.process.transcription.rnaExpression["basal"].copy(),
		"doubling_time": sim_data.conditionToDoublingTime["basal"],
		}

	# Determine expression and synthesis probabilities
	expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, _ = expressionConverge(
		sim_data,
		cellSpecs["basal"]["expression"],
		cellSpecs["basal"]["concDict"],
		cellSpecs["basal"]["doubling_time"],
		)

	# Store calculated values
	cellSpecs["basal"]["expression"] = expression
	cellSpecs["basal"]["synthProb"] = synthProb
	cellSpecs["basal"]["avgCellDryMassInit"] = avgCellDryMassInit
	cellSpecs["basal"]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
	cellSpecs["basal"]["bulkContainer"] = bulkContainer

	# Modify sim_data mass
	sim_data.mass.avgCellDryMassInit = avgCellDryMassInit
	sim_data.mass.avgCellDryMass = sim_data.mass.avgCellDryMassInit * sim_data.mass.avgCellToInitialCellConvFactor
	sim_data.mass.avgCellWaterMassInit = sim_data.mass.avgCellDryMassInit / sim_data.mass.cellDryMassFraction * sim_data.mass.cellWaterMassFraction
	sim_data.mass.fitAvgSolubleTargetMolMass = fitAvgSolubleTargetMolMass

	# Modify sim_data expression
	sim_data.process.transcription.rnaExpression["basal"][:] = cellSpecs["basal"]["expression"]
	sim_data.process.transcription.rnaSynthProb["basal"][:] = cellSpecs["basal"]["synthProb"]

	return cellSpecs

def calculateTranslationSupply(sim_data, doubling_time, bulkContainer, avgCellDryMassInit):
	"""
	Returns the supply rates of all amino acids to translation given the desired
	doubling time. This creates a limit on the polypeptide elongation process,
	and thus on growth. The amino acid supply rate is found by calculating the
	concentration of amino acids per gram dry cell weight and multiplying by the
	loss to dilution given doubling time.

	Requires
	--------
	- doubling_time: measured doubling times given the condition, in units of minutes.
	- bulkContainer: a container that tracks the counts of all bulk molecules
	- avgCellDryMassInit: The average initial cell dry mass, in units of fg.

	Notes
	-----
	- The supply of amino acids should not be based on a desired doubling time,
	but should come from a more mechanistic basis. This would allow simulations
	of environmental shifts in which the doubling time is unknown.
	"""

	aaCounts = sim_data.process.translation.monomerData["aaCounts"] # the counts of each amino acid required for each protein
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomerData["id"]) # the counts of all proteins
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

	# Calculate required amino acid supply to translation to counter dilution
	translation_aa_supply = molAAPerGDCW * np.log(2) / doubling_time

	return translation_aa_supply

def buildTfConditionCellSpecifications(sim_data, tf):
	"""
	Creates cell specifications for a given transcription factor by
	fitting expression. Will set for the active and inactive TF condition.
	Relies on expressionConverge() to set the expression and masses.
	Uses fold change data relative to the 'basal' condition to determine
	expression for a given TF.

	Inputs
	------
	- tf (str) - label for the transcription factor to fit (eg. 'CPLX-125')

	Requires
	--------
	- Metabolite concentrations based on nutrients for the TF
	- Adjusted 'basal' RNA expression
	- Doubling time for the TF
	- Fold changes in expression for each gene given the TF

	Returns
	--------
	- dict {tf + '__active'/'__inactive': dict} with the following keys in each dict:
		'concDict' {metabolite_name (str): concentration (float with units)} -
			dictionary of concentrations for each metabolite with a concentration
		'expression' (1D numpy array of floats) - expression for each RNA,
			total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (1D numpy array of floats) - synthesis probability for
			each RNA, total normalized to 1
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - TODO - what is this?
		'bulkContainer' (BulkObjectsContainer object) - expected counts for
			bulk molecules based on expression
	"""

	cellSpecs = {}
	for choice in ["__active", "__inactive"]:
		conditionKey = tf + choice
		conditionValue = sim_data.conditions[conditionKey]

		# Get expression for the condition based on fold changes over 'basal'
		# condition if the condition is not the same as 'basal'
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

		# Get metabolite concentrations for the condition
		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[conditionKey]))

		# Create dictionary for the condition
		cellSpecs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.conditionToDoublingTime.get(
				conditionKey,
				sim_data.conditionToDoublingTime["basal"]
				)
			}

		# Determine expression and synthesis probabilities
		expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cellSpecs[conditionKey]["expression"],
			cellSpecs[conditionKey]["concDict"],
			cellSpecs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rnaData["KmEndoRNase"],
			updateConcDict = True,
			)

		# Store calculated values
		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

	return cellSpecs

	# TODO - what is this and what is the calculation that should be done? remove?
		# Uncomment when concDict is actually calculated for non-base [AA]
		# if len(conditionValue["perturbations"]) == 0:
		# 	nutrientLabel = conditionValue["nutrients"]
		# 	sim_data.process.metabolism.nutrientsToInternalConc[nutrientLabel] = concDict

def buildCombinedConditionCellSpecifications(sim_data, cellSpecs):
	"""
	Creates cell specifications for sets of transcription factors being active.
	These sets include conditions like 'with_aa' or 'no_oxygen' where multiple
	transcription factors will be active at the same time.

	Inputs
	------
	- cellSpecs {condition: dict} - information about each individual
	transcription factor condition

	Requires
	--------
	- Metabolite concentrations based on nutrients for the condition
	- Adjusted 'basal' RNA expression
	- Doubling time for the combined condition
	- Fold changes in expression for each gene given the TF

	Modifies
	--------
	- cellSpecs dictionary for each combined condition
	- RNA expression and synthesis probabilities for each combined condition

	Notes
	-----
	- TODO - include TFs in the no_oxygen condition
	- TODO - determine how to handle fold changes if multiple TFs change the
	same gene (currently, an exception is raised)
	"""

	fcData = {}
	for conditionKey in sim_data.conditionActiveTfs:
		# Skip adjustments if 'basal' condition
		if conditionKey == "basal":
			continue

		# Get expression from fold changes for each TF in the given condition
		conditionValue = sim_data.conditions[conditionKey]
		for tf in sim_data.conditionActiveTfs[conditionKey]:
			for gene in sim_data.tfToFC[tf]:
				if gene in fcData:
					# fcData[gene] *= sim_data.tfToFC[tf][gene]
					raise Exception("Check this implementation: multiple genes regulated")
				else:
					fcData[gene] = sim_data.tfToFC[tf][gene]
		expression = expressionFromConditionAndFoldChange(
			sim_data.process.transcription.rnaData["id"],
			sim_data.process.transcription.rnaExpression["basal"],
			conditionValue["perturbations"],
			fcData,
			)

		# Get metabolite concentrations for the condition
		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[conditionKey]))

		# Create dictionary for the condition
		cellSpecs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.conditionToDoublingTime.get(
				conditionKey,
				sim_data.conditionToDoublingTime["basal"]
				)
			}

		# Determine expression and synthesis probabilities
		expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cellSpecs[conditionKey]["expression"],
			cellSpecs[conditionKey]["concDict"],
			cellSpecs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rnaData["KmEndoRNase"],
			updateConcDict = True,
			)

		# Modify cellSpecs for calculated values
		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

		# Modify sim_data expression
		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

		# TODO - same as TODO comment on block in function above
		# Uncomment when concDict is actually calculated for non-base [AA]
		# if len(conditionValue["perturbations"]) == 0:
		# 	nutrientLabel = conditionValue["nutrients"]
		# 	sim_data.process.metabolism.nutrientsToInternalConc[nutrientLabel] = concDict


def expressionConverge(sim_data, expression, concDict, doubling_time, Km = None, updateConcDict = False):
	"""
	Iteratively fits synthesis probabilities for RNA. Calculates initial
	expression based on gene expression data and makes adjustments to match
	physiological constraints for ribosome and RNAP counts. Relies on
	fitExpression() to converge

	Inputs
	--------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- concDict {metabolite (str): concentration (float with units)} - dictionary
	for concentrations of each metabolite with location tag
	- doubling_time (float with units) - doubling time
	- Km (array of floats with concentration units) - Km for each RNA associated
	with RNases
	- updateConcDict - TODO - remove?

	Requires
	--------
	- MAX_FITTING_ITERATIONS (int) - number of iterations to adjust expression
	before an exception is raised
	- FITNESS_THRESHOLD (float) - acceptable change from one iteration to break
	the fitting loop

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1
	- synthProb (array of floats) - synthesis probability for each RNA which
	accounts for expression and degradation rate, normalized to 1
	- avgCellDryMassInit (float with units) - expected initial dry cell mass
	- fitAvgSolubleTargetMolMass (float with units) - TODO - what is this?
	- bulkContainer (BulkObjectsContainer object) - expected counts for
	bulk molecules based on expression
	- concDict - TODO - remove?

	Notes
	-----
	- TODO - remove updateConcDict? - doesn't get updated at all
	"""

	if VERBOSE > 0:
		print("Fitting RNA synthesis probabilities.")

	for iteration in xrange(MAX_FITTING_ITERATIONS):
		if VERBOSE > 1:
			print('Iteration: {}'.format(iteration))

		initialExpression = expression.copy()
		expression = setInitialRnaExpression(sim_data, expression, doubling_time)
		bulkContainer = createBulkContainer(sim_data, expression, doubling_time)
		avgCellDryMassInit, fitAvgSolubleTargetMolMass = rescaleMassForSolubleMetabolites(sim_data, bulkContainer, concDict, doubling_time)

		setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time)
		setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		# Normalize expression and write out changes
		expression, synthProb = fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		# TODO - remove?
		if updateConcDict:
			concDict = concDict.copy() # Calculate non-base condition [AA]

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - expression)))
		if VERBOSE > 1:
			print('degree of fit: {}'.format(degreeOfFit))

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	return expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict

def fitCondition(sim_data, spec, condition):
	"""
	Takes a given condition and returns the predicted bulk average, bulk deviation,
	protein monomer average, protein monomer deviation, and amino acid supply to
	translation. This relies on calculateBulkDistributions and calculateTranslationSupply.

	Requires
	--------
	- condition: a string specifying the condition.
	- spec: a dictionary with cell specifications for the given condition. This
	function uses the specs "expression", "concDict", "avgCellDryMassInit", and
	"doubling_time"

	Returns
	--------
	- The following specs are updated:
		- bulkAverageContainer: The mean of the bulk counts.
		- bulkDeviationContainer: The standard deviation of the bulk counts.
		- proteinMonomerAverageContainer: The mean of the protein monomer counts.
		- proteinMonomerDeviationContainer: The standard deviation of the protein
		monomer	counts.
		- translation_aa_supply: the supply rates of all amino acids to translation.
	"""

	if VERBOSE > 0:
		print "Fitting condition {}".format(condition)

	# Find bulk and protein distributions
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

	# Find the supply rates of amino acids to translation given doubling time
	spec["translation_aa_supply"] = calculateTranslationSupply(
											sim_data,
											spec["doubling_time"],
											spec["proteinMonomerAverageContainer"],
											spec["avgCellDryMassInit"],
											)

	return {condition: spec}

# Sub-fitting functions

def setRnaPolymeraseCodingRnaDegradationRates(sim_data):
	"""
	Increase the degradation rates for the RNA polymerase mRNAs.  This is done to increase the
	rate of	mRNA synthesis and overall reduce the stochasticity in RNA polymerase subunit
	expression, which would otherwise constrain transcription.

	Requires
	--------
	- RNA_POLY_MRNA_DEG_RATE_PER_S: The new first-order degradation rate, in units of per second.

	Modifies
	--------
	- Degradation rates of RNA polymerase subunit mRNAs.

	Notes
	-----
	- Incorporating transcription unit structure would facilitate co-expression of the subunits
		but might not address the fundamental stochasticity issue.
	"""

	rnaPolySubunits = sim_data.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"] # APORNAP-CPLX[c] is the RNA polymerase complex
	subunitIndexes = np.array([np.where(sim_data.process.translation.monomerData["id"] == id_)[0].item() for id_ in rnaPolySubunits]) # there has to be a better way...
	mRNA_indexes = sim_data.relation.rnaIndexToMonomerMapping[subunitIndexes]

	# Modifies
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

	totalCount_rRNA16S = totalCountFromMassesAndRatios(
		totalMass_rRNA16S,
		individualMasses_rRNA16S,
		distribution_rRNA16S
		)

	totalCount_rRNA5S = totalCountFromMassesAndRatios(
		totalMass_rRNA5S,
		individualMasses_rRNA5S,
		distribution_rRNA5S
		)

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

	counts_tRNA = totalCount_tRNA * distribution_tRNA

	rnaExpressionContainer.countsIs(counts_tRNA, ids_tRNA)

	## Assign mRNA counts based on mass and relative abundances (microarrays)

	totalCount_mRNA = totalCountFromMassesAndRatios(
		totalMass_mRNA,
		individualMasses_mRNA,
		distribution_mRNA
		)

	counts_mRNA = totalCount_mRNA * distribution_mRNA

	rnaExpressionContainer.countsIs(counts_mRNA, ids_mRNA)

	expression = normalize(rnaExpressionContainer.counts())

	return expression
	# Note that now rnaData["synthProb"] does not match "expression"

def totalCountIdDistributionProtein(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of proteins from the relative expression of RNA,
	individual protein mass, and total protein mass. Relies on the math functions
	netLossRateFromDilutionAndDegradationProtein, proteinDistributionFrommRNA,
	totalCountFromMassesAndRatios.

	Requires
	--------
	- expression: relative frequency distribution of RNA expression.
	- doubling_time: measured doubling times given the condition, in units of minutes.

	Returns
	--------
	- Counts and IDs of all proteins.
	"""

	ids_protein = sim_data.process.translation.monomerData["id"]
	totalMass_protein = sim_data.mass.getFractionMass(doubling_time)["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_protein = sim_data.process.translation.monomerData["mw"] / sim_data.constants.nAvogadro
	distribution_transcriptsByProtein = normalize(expression[sim_data.relation.rnaIndexToMonomerMapping])
	translation_efficienciesByProtein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	degradationRates = sim_data.process.translation.monomerData["degRate"]

	# Find the net protein loss
	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doubling_time, degradationRates)

	# Find the protein distribution
	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcriptsByProtein,
		translation_efficienciesByProtein,
		netLossRate_protein
		)

	# Find total protein counts
	totalCount_protein = totalCountFromMassesAndRatios(
		totalMass_protein,
		individualMasses_protein,
		distribution_protein
		)

	return totalCount_protein, ids_protein, distribution_protein

def totalCountIdDistributionRNA(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of RNA from their relative expression, individual
	mass, and total RNA mass. Relies on the math function totalCountFromMassesAndRatios.

	Requires
	--------
	- expression: relative frequency distribution of RNA expression.
	- doubling_time: measured doubling times given the condition, in units of minutes.

	Returns
	--------
	- Counts and IDs of all RNAs.
	"""

	ids_rnas = sim_data.process.transcription.rnaData["id"]
	totalMass_RNA = sim_data.mass.getFractionMass(doubling_time)["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_RNA = sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro

	distribution_RNA = normalize(expression)

	totalCount_RNA = totalCountFromMassesAndRatios(
		totalMass_RNA,
		individualMasses_RNA,
		distribution_RNA
		)

	return totalCount_RNA, ids_rnas, distribution_RNA

def createBulkContainer(sim_data, expression, doubling_time):
	"""
	Creates a container that tracks the counts of all bulk molecules. Relies on
	totalCountIdDistributionRNA and totalCountIdDistributionProtein to set the
	counts and IDs of all RNAs and proteins.

	Requires
	--------
	- expression: the relative frequency distribution of RNA expression.
	- doubling_time: the measured doubling times given the condition, in units
	of minutes.

	Returns
	-------
	- bulkContainer: a wrapper around a NumPy array that tracks the counts of
	bulk molecules.
	"""

	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)
	ids_molecules = sim_data.state.bulkMolecules.bulkData["id"]

	# Construct bulk container
	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	# Assign RNA counts based on mass and expression distribution
	counts_RNA = totalCount_RNA * distribution_RNA
	bulkContainer.countsIs(counts_RNA, ids_rnas)

	# Assign protein counts based on mass and mRNA counts
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
	"""
	Fits the growth-associated maintenance (GAM) cost associated with metabolism.

	The energetic costs associated with growth have been estimated utilizing flux-balance analysis
	and are used with FBA to obtain accurate growth predictions.  In the whole-cell model, some of
	these costs are explicitly associated with the energetic costs of translation, a biomass
	assembly process.  Consequently we must estimate the amount of energy utilized by translation
	per unit of biomass (i.e. dry mass) produced, and subtract that quantity from reported GAM to
	acquire the modified GAM that we use in the metabolic submodel.

	Requires
	--------
	- amino acid counts associated with protein monomers
	- average initial dry mass
	- energetc (GTP) cost of translation (per amino acid polymerized)
	- observed growth-associated maintenance (GAM)
		In dimensions of ATP or ATP equivalents consumed per biomass

	Modifies
	--------
	- the "dark" ATP, i.e. the modified GAM

	Notes
	-----
	As more non-metabolic submodels account for energetic costs, this function should be extended
	to subtract those costs off the observed GAM.

	There also exists, in contrast, non-growth-associated-maintenance (NGAM), which is relative to
	total biomass rather than the biomass accumulation rate.  As the name would imply, this
	accounts for the energetic costs of maintaining the existing biomass.  It is also accounted for
	in the metabolic submodel.

	TODO (John): Rewrite as a true function.

	"""
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
	"""
	Finds a distribution of copy numbers for macromolecules. While RNA and protein
	expression can be approximated using well-described statistical	distributions,
	complexes require absolute copy numbers. To get these distributions, this
	function instantiates many cells with a reduced set of molecules, forms complexes,
	and iterates through equilibrium and two-component system processes until
	metabolite counts reach a steady-state. It then computes the resulting
	statistical distributions.

	Requires
	--------
	- N_SEEDS: The number of instantiated cells.
	- expression: Relative frequency distribution of RNA expression.
	- concDict: A dictionary of a small set of metabolites and their concentrations,
	in units of mol/L.
	- avgCellDryMassInit: The average initial cell dry mass, in units of fg.
	- doubling_time: Measured doubling times given the condition, in units of minutes.

	Returns
	--------
	- bulkAverageContainer: The mean of the bulk counts.
	- bulkDeviationContainer: The standard deviation of the bulk counts.
	- proteinMonomerAverageContainer: The mean of the protein monomer counts.
	- proteinMonomerDeviationContainer: The standard deviation of the protein monomer
	counts.

	TODO (ERAN): How does mccFormComplexesWithPrebuiltMatrices work?
	"""

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
	bulkContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'])
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

	# Instantiate cells to find average copy numbers of macromolecules
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

		# Form complexes
		updatedCompMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			complexationMoleculeCounts,
			seed,
			complexationStoichMatrix,
			*complexationPrebuiltMatrices
			)

		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		metDiffs = np.inf * np.ones_like(metabolitesView.counts())
		nIters = 0

		# Iterate processes until metabolites converge to a steady-state
		while(np.linalg.norm(metDiffs, np.inf) > 1):
			metCounts = conc_metabolites * cellVolume * sim_data.constants.nAvogadro
			metCounts.normalize()
			metCounts.checkNoUnit()
			metabolitesView.countsIs(
				metCounts.asNumber().round()
				)

			# Find reaction fluxes from equilibrium process
			rxnFluxes, _ = sim_data.process.equilibrium.fluxesAndMoleculesToSS(
				equilibriumMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.nAvogadro.asNumber(1 / units.mol),
				)
			equilibriumMoleculesView.countsInc(
				np.dot(sim_data.process.equilibrium.stoichMatrix().astype(np.int64), rxnFluxes)
				)
			assert np.all(equilibriumMoleculesView.counts() >= 0)

			# Find changes from two component system
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


	# Update counts in bulk objects container
	bulkAverageContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(sim_data.state.bulkMolecules.bulkData['id'], np.float64)
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
	Function to determine the expected total counts for a group of molecules
	in order to achieve a total mass with a given distribution of individual
	molecules.

	Math:
		Total mass = dot(mass, count)

		Fraction of i:
		f = count / Total counts

		Substituting:
		Total mass = dot(mass, f * Total counts)
		Total mass = Total counts * dot(mass, f)

		Total counts = Total mass / dot(mass, f)

	Requires
	--------
	- totalMass (float with mass units): total mass of the group of molecules
	- individualMasses (array of floats with mass units): mass for individual
	molecules in the group
	- distribution (array of floats): distribution of individual molecules,
	normalized to 1

	Returns
	--------
	- float of the total counts (does not need to be a whole number)

	"""

	assert np.allclose(np.sum(distribution), 1)
	return (1 / units.dot(individualMasses, distribution) * totalMass).asNumber()


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
	"""
	Adjusts expression of RNA based on fold changes from basal for a given condition.

	Inputs
	------
	- rnaIds (array of str) - name of each RNA with location tag
	- basalExpression (array of floats) - expression for each RNA in the basal
	condition, normalized to 1
	- condPerturbations {RNA ID with location tag (str): fold change (float)} -
	dictionary of fold changes for RNAs based on the given condition
	- tfFCs {RNA ID without location tag (str): fold change (float)} -
	dictionary of fold changes for RNAs based on transcription factors in the
	given condition

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1

	Notes
	-----
	- TODO (Travis) - Might not properly handle if an RNA is adjusted from both a
	perturbation and a transcription factor, currently RNA self regulation is not
	included in tfFCs
	"""

	expression = basalExpression.copy()

	# Gather RNA indices and fold changes for each RNA that will be adjusted
	rnaIdxs = []
	fcs = []
	for key in sorted(condPerturbations):
		value = condPerturbations[key]
		rnaIdxs.append(np.where(rnaIds == key)[0][0])
		fcs.append(value)
	for key in sorted(tfFCs):
		rnaIdxs.append(np.where(rnaIds == key + "[c]")[0][0])
		fcs.append(tfFCs[key])

	# Sort fold changes and indices for the bool array indexing to work properly
	fcs = [fc for (rnaIdx, fc) in sorted(zip(rnaIdxs, fcs), key = lambda pair: pair[0])]
	rnaIdxs = [rnaIdx for (rnaIdx, fc) in sorted(zip(rnaIdxs, fcs), key = lambda pair: pair[0])]

	# Adjust expression based on fold change and normalize
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
	'''
	Calculates the probabilities (P) that each transcription factor will bind
	to its target RNA. This function initially calculates these probabilities
	from the bulk average counts of the TFs and ligands calculated from
	previous steps. Then, values of parameters alpha and r in the equation
	below are fit such that the computed RNA synthesis probabilities converge
	to the measured RNA synthesis probabilities.

	v_{synth, j} = \alpha_j + \sum_{i} P_{T,i}*r_{ij}
	(See supplementary materials on transcription regulation for details)

	Due to constraints applied in the optimization, both v and P need to
	be shifted from their initial values. Using the fit values of P, the set
	concentrations of ligand metabolites and the kd's of the ligand-TF binding
	reactions are also updated.

	Requires
	--------
	- Bulk average counts of transcription factors and associated ligands
	for each condition (in cellSpecs)

	Modifies
	--------
	Probabilities of TFs binding to their promoters
	RNA synthesis probabilities
	Set concentrations of metabolites that act as ligands in 1CS
	kd's of equilibrium reactions in 1CS

	Returns
	--------
	- r: Fit parameters on how the recruitment of a TF affects the expression
	of a gene. High (positive) values of r imply that the TF binding increases
	the probability that the gene is expressed.
	'''
	def makeG(sim_data, pPromoterBound):
		"""
		Construct matrix G that contains probabilities of pPromoterBound as
		elements. Each row of the matrix is named "[RNA]__[condition]", where
		there are two conditions [active/inactive] for each TF that regulates
		the expression of the given RNA. For RNAs that are not regulated by
		any TFs, a single row named "[RNA]__basal" represents the RNA. Each
		column is named "[RNA]__[TF]", for each TF that regulates the
		expression of the given RNA. Each element is set to the value in
		pPromoterBound that corresponds to the condition given by the row,
		and the TF given by the column. For each RNA, there is an additional
		column named "[RNA]__alpha", and all elements in this column that
		corresponds to the rows for the RNA are set to 1.

		Requires
		--------
		- pPromoterBound: Probabilities that a given TF is bound to its
		promoter in a given condition, calculated from bulk average
		concentrations of the TF and its associated ligands.

		Returns
		--------
		- G: Matrix of values in pPromoterBound, rearranged based on each RNA
		- rowNames: List of row names of G as strings
		- colNames: List of column names of G as strings
		- k: List of RNA synthesis probabilities per each RNA and condition
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		"""
		gI, gJ, gV, k, rowNames, colNames, kInfo = [], [], [], [], [], [], []

		for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
			rnaIdNoLoc = rnaId[:-3]  # Strip off compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			# TODO: cache this list of TFs for each RNA
			for tf in tfs:
				if tf not in sorted(sim_data.tfToActiveInactiveConds):
					continue

				# Add conditions for selected TFs
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)

			for condition in conditions:
				# Skip basal conditions, unless the RNA is not regulated by any TFs
				if len(tfsWithData) > 0 and condition == "basal":
					continue

				# Add row for each condition specific to each RNA
				rowName = rnaIdNoLoc + "__" + condition
				rowNames.append(rowName)

				for tf in tfsWithData:
					# Add column for each TF that regulates each RNA
					colName = rnaIdNoLoc + "__" + tf

					if colName not in colNames:
						colNames.append(colName)

					gI.append(rowNames.index(rowName))
					gJ.append(colNames.index(colName))
					gV.append(pPromoterBound[condition][tf])  # Probability that TF is bound in given condition

				# Add alpha column for each RNA
				colName = rnaIdNoLoc + "__alpha"

				if colName not in colNames:
					colNames.append(colName)

				gI.append(rowNames.index(rowName))
				gJ.append(colNames.index(colName))
				gV.append(1.)  # TODO: Why is this set to 1?

				# Also gather RNA synthesis probabilities for each RNA per condition
				k.append(sim_data.process.transcription.rnaSynthProb[condition][idx])
				kInfo.append({"condition": condition, "idx": idx})

		k = np.array(k)
		gI, gJ, gV = np.array(gI), np.array(gJ), np.array(gV)
		G = np.zeros((len(rowNames), len(colNames)), np.float64)
		G[gI, gJ] = gV

		return G, rowNames, colNames, k, kInfo

	def makeZ(sim_data, colNames):
		"""
		Construct matrix Z that connects all possible TF combinations with
		each TF. Each row of the matrix corresponds to an RNA-(TF combination)
		pair, and each column corresponds to an RNA-TF pair, with an additional
		RNA-alpha column for each RNA (identical to matrix G). Matrix values
		are set to one if the TF specified by the column is "active" in the
		combination specified by the row or if the column is an RNA-alpha
		column, and zero otherwise.

		Requires
		--------
		- colNames: List of column names from matrix G.

		Returns
		--------
		- Z: Matrix of zeros and ones, specifying which TFs in the columns
		correspond to combinations in the rows.
		"""
		# TODO: refactor this as a function, not a fixed dictionary
		combinationIdxToColIdxs = {
			0: [0], 1: [0, 1], 2: [0, 2], 3: [0, 1, 2],
			4: [0, 3], 5: [0, 1, 3], 6: [0, 2, 3], 7: [0, 1, 2, 3],
			8: [0, 4], 9: [0, 1, 4], 10: [0, 2, 4], 11: [0, 1, 2, 4],
			12: [0, 3, 4], 13: [0, 1, 3, 4], 14: [0, 2, 3, 4], 15: [0, 1, 2, 3, 4],
			}

		zI, zJ, zV, rowNames = [], [], [], []

		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]  # Strip off compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			tfsWithData = []

			# Get column index of the RNA's alpha column
			colIdxs = [colNames.index(rnaIdNoLoc + "__alpha")]

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sim_data.tfToActiveInactiveConds:
					continue

				tfsWithData.append(tf)

				# Get column index of the RNA-TF pair
				colIdxs.append(colNames.index(rnaIdNoLoc + "__" + tf))

			nTfs = len(tfsWithData)

			# For all possible combinations of TFs
			for combinationIdx in xrange(2**nTfs):
				# Add a row for each combination
				rowName = rnaIdNoLoc + "__%d" % combinationIdx
				rowNames.append(rowName)

				# Set matrix value to one if the TF specified by the column is
				# present in the combination of TFs specified by the row
				for colIdx in combinationIdxToColIdxs[combinationIdx]:
					zI.append(rowNames.index(rowName))
					zJ.append(colIdxs[colIdx])
					zV.append(1)

		# Build matrix Z
		zI, zJ, zV = np.array(zI), np.array(zJ), np.array(zV)
		Z = np.zeros((zI.max() + 1, zJ.max() + 1), np.float64)
		Z[zI, zJ] = zV

		return Z

	def makeT(sim_data, colNames):
		"""
		Construct matrix T that specifies the direction of regulation for each
		RNA-TF pair.

		Requires
		--------
		- colNames: List of column names from matrix G.

		Returns
		--------
		- T: Diagonal matrix. Diagonal value is +1 if the direction of
		regulation by the TF-RNA pair specified by the row is positive, -1 if
		this is negative, and 0 if the row is an RNA_alpha row.
		"""
		tI, tJ, tV, rowNamesT = [], [], [], []

		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]  # Strip off compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sim_data.tfToActiveInactiveConds:
					continue

				tfsWithData.append(tf)

			for tf in tfsWithData:
				# Add row for TF and find column for TF in colNames
				rowName = rnaIdNoLoc + "__" + tf
				rowNamesT.append(rowName)
				colName = rnaIdNoLoc + "__" + tf

				# Set matrix value to regulation direction (+1 or -1)
				tI.append(rowNamesT.index(rowName))
				tJ.append(colNames.index(colName))
				tV.append(sim_data.tfToDirection[tf][rnaIdNoLoc])

			# Add RNA_alpha rows and columns, and set matrix value to zero
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
		"""
		Construct matrix H that contains values of vector r as elements.
		Each row of the matrix is named "[RNA]__[condition]", where
		there are two conditions [active/inactive] for each TF that regulates
		the expression of the given RNA. For RNAs that are not regulated by
		any TFs, a single row named "[RNA]__basal" represents the RNA. Each
		column is named "[TF]__[condition]", for each TF that regulates the
		expression of the given RNA, and each condition given for the RNA.
		Each element is set to the optimized value in r that corresponds to
		the RNA given by the row, and the TF given by the column. For each RNA,
		there is an additional column named "[RNA]__alpha", and all elements in
		this column are set to the value of parameter alpha for the RNA
		optimized in r.

		Requires
		--------
		- colNames: List of column names from matrix G.
		- pPromoterBound: Probabilities that a given TF is bound to its
		promoter in a given condition, calculated from bulk average
		concentrations of the TF and its associated ligands.
		- r: Optimized values of \Delta r (the effect of TF on RNAP
		recruitment) and alpha (basal recruitment of RNAP)
		- fixedTFs: List of TFs whose activities do not change with the
		nutrient conditions

		Returns
		--------
		- H: Matrix of values in optimized r, rearranged for each RNA and
		condition
		- pInit: Vector of values in pPromoterBound, rearranged to the ordering
		of the columns of H
		- pAlphaIdxs: Indexes of columns that correspond to alpha's in H and pInit
		- pNotAlphaIdxs: Indexes of columns that correspond to r's in H and pInit
		- fixedTFIdxs: Indexes of columns that correspond to fixed TFs in H and pInit
		- pPromoterBoundIdxs: Dictionary of indexes to pInit.
		- colNamesH: List of column names of H as strings
		"""
		rDict = dict([(colName, value) for colName, value in zip(colNames, r)])

		pPromoterBoundIdxs = dict([(condition, {}) for condition in pPromoterBound])
		hI, hJ, hV, rowNames, colNamesH, pInitI, pInitV = [], [], [], [], [], [], []

		for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
			rnaIdNoLoc = rnaId[:-3]  # Strip off compartment ID from RNA ID

			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sorted(sim_data.tfToActiveInactiveConds):
					continue

				# Add conditions for selected TFs
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)

			for condition in conditions:
				# Skip basal conditions, unless the RNA is not regulated by any TFs
				if len(tfsWithData) > 0 and condition == "basal":
					continue

				# Add row for each condition specific to each RNA
				rowName = rnaIdNoLoc + "__" + condition
				rowNames.append(rowName)

				for tf in tfsWithData:
					# Add column for each TF and condition
					colName = tf + "__" + condition

					if colName not in colNamesH:
						colNamesH.append(colName)

					hI.append(rowNames.index(rowName))
					hJ.append(colNamesH.index(colName))

					# Handle the case of the TF being knocked out (admittedly not the cleanest solution)
					if cellSpecs[condition]["bulkAverageContainer"].count(tf + "[c]") == 0:
						hV.append(0)  # TF is knocked out in the given condition
					else:
						hV.append(rDict[rnaIdNoLoc + "__" + tf])  # Optimized r value for TF-RNA pair

					# Rearrange values in pPromoterBound in the same order
					# given by the columns of H
					pInitI.append(colNamesH.index(colName))
					pInitV.append(pPromoterBound[condition][tf])
					pPromoterBoundIdxs[condition][tf] = colNamesH.index(colName)

				# Add alpha column for each RNA
				colName = rnaIdNoLoc + "__alpha"

				if colName not in colNamesH:
					colNamesH.append(colName)

				# Add optimized value of alpha in r to H
				hI.append(rowNames.index(rowName))
				hJ.append(colNamesH.index(colName))
				hV.append(rDict[colName])

				# Set corresponding value in pInit to one
				pInitI.append(colNamesH.index(colName))
				pInitV.append(1.)

		# Build vector pInit and matrix H
		pInit = np.zeros(len(set(pInitI)))
		pInit[pInitI] = pInitV

		hI, hJ, hV = np.array(hI), np.array(hJ), np.array(hV)
		Hshape = (hI.max() + 1, hJ.max() + 1)
		H = np.zeros(Hshape, np.float64)
		H[hI, hJ] = hV

		# Get indexes of alpha and non-alpha columns in pInit and H
		pAlphaIdxs = np.array([colNamesH.index(colName) for colName in colNamesH if colName.endswith("__alpha")])
		pNotAlphaIdxs = np.array([colNamesH.index(colName) for colName in colNamesH if not colName.endswith("__alpha")])

		# Get indexes of columns that correspond to fixed TFs
		fixedTFIdxs = []
		for idx, colName in enumerate(colNamesH):
			secondElem = colName.split("__")[1]

			if secondElem in fixedTFs:
				fixedTFIdxs.append(idx)

		fixedTFIdxs = np.array(fixedTFIdxs, dtype=np.int)

		return H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, colNamesH

	def makePdiff(sim_data, colNamesH):
		"""
		Construct matrix Pdiff that specifies the indexes of corresponding
		TFs and conditions.

		Requires
		--------
		- colNamesH: List of column names from matrix H.

		Returns
		--------
		- Pdiff: Matrix with [TF] as rows and [TF]_[condition] as columns.
		Matrix value is set to 1 when the TF of the column matches with the TF
		of the row, and the condition is TF__active. Matrix value is set to -1
		when the TF of the column matches with the TF of the row, and the
		condition is TF__inactive.
		"""
		PdiffI, PdiffJ, PdiffV = [], [], []

		for rowIdx, tf in enumerate(sorted(sim_data.tfToActiveInactiveConds)):
			# For each TF, find condition [TF]__[TF]__active and set element to 1
			condition = tf + "__active"
			colName = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(colNamesH.index(colName))
			PdiffV.append(1)

			# Find condition [TF]__[TF]__inactive and set element to -1
			condition = tf + "__inactive"
			colName = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(colNamesH.index(colName))
			PdiffV.append(-1)

		# Build matrix Pdiff
		PdiffI, PdiffJ, PdiffV = np.array(PdiffI), np.array(PdiffJ), np.array(PdiffV)
		Pdiffshape = (PdiffI.max() + 1, len(colNamesH))
		Pdiff = np.zeros(Pdiffshape, np.float64)
		Pdiff[PdiffI, PdiffJ] = PdiffV

		return Pdiff

	def fromArray(p, pPromoterBound, pPromoterBoundIdxs):
		'''
		Updates values in pPromoterBound with fitted probabilities.

		Requires
		--------
		- p: Vector of probabilities optimized in the current step.
		- pPromoterBoundIdxs: Dictionary of indexes to p

		Modifies
		--------
		Values in pPromoterBound - probabilities that each transcription factor
		is bound to its promoter for each growth condition
		'''
		for condition in sorted(pPromoterBoundIdxs):
			for tf in sorted(pPromoterBoundIdxs[condition]):
				pPromoterBound[condition][tf] = p[pPromoterBoundIdxs[condition][tf]]

	def updateSynthProb(sim_data, kInfo, k):
		'''
		Updates RNA synthesis probabilities with fitted values of P and R.

		Requires
		--------
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		- k: RNA synthesis probabilities computed from fitted P and R.

		Modifies
		--------
		- RNA synthesis probabilities
		'''
		# Update sim_data values with fitted values
		for D, value in zip(kInfo, k):
			sim_data.process.transcription.rnaSynthProb[D["condition"]][D["idx"]] = value

		# Normalize values such that probabilities for each condition sum to one
		for condition in sim_data.process.transcription.rnaSynthProb:
			assert np.all(sim_data.process.transcription.rnaSynthProb[condition] >= 0)
			sim_data.process.transcription.rnaSynthProb[condition] /= sim_data.process.transcription.rnaSynthProb[condition].sum()

	# Initialize pPromoterBound using mean TF and ligand concentrations
	pPromoterBound = calculatePromoterBoundProbability(sim_data, cellSpecs)
	pInit0 = None
	lastNorm = np.inf

	fixedTFs = []
	for tf in sim_data.tfToActiveInactiveConds:
		if sim_data.process.transcription_regulation.tfToTfType[tf] == "2CS":
			fixedTFs.append(tf)
		if sim_data.process.transcription_regulation.tfToTfType[tf] == "1CS" and sim_data.tfToActiveInactiveConds[tf]["active nutrients"] == sim_data.tfToActiveInactiveConds[tf]["inactive nutrients"]:
			fixedTFs.append(tf)

	SCALING = 10  # Multiplied to all matrices for numerical stability
	NORM_TYPE = 1  # Matrix 1-norm
	MOMENTUM_COEFF = 1e-3  # How strongly should the probabilities adhere to original values?
	# TODO: A better name for this constant?

	# Repeat for a fixed maximum number of iterations
	for i in xrange(100):
		# Construct matrices used in optimizing R
		# TODO: Make these matrix names more meaningful
		# TODO: separate the routine that gets k and kInfo
		G, rowNamesG, colNamesG, k, kInfo = makeG(sim_data, pPromoterBound)
		Z = makeZ(sim_data, colNamesG)
		T = makeT(sim_data, colNamesG)

		# Optimize R such that the RNA polymerase initiation probabilities
		# computed from current values of P in matrix G are close to fitted
		# values.
		R = Variable(G.shape[1])  # Vector of r's and alpha's

		# Objective: minimize difference between k (fitted RNAP initiation
		# probabilities) and G*R (computed initiation probabilities)
		objective_r = Minimize(norm(G*(SCALING*R) - SCALING*k, NORM_TYPE))

		# Constraints
		# 1) 0 <= Z*R <= 1 : Assuming P = 1 for all TFs, all possible
		# combinations of TFs should yield a valid RNAP initialization
		# probability value between zero and one.
		# 2) T*R >= 0 : Values of r for positive regulation should be positive,
		# and values of r for negative regulation should be negative.
		constraint_r = [0 <= Z*(SCALING*R), Z*(SCALING*R) <= SCALING*1, T*(SCALING*R) >= 0]

		# Solve optimization problem
		prob_r = Problem(objective_r, constraint_r)
		prob_r.solve(solver = "GLPK")

		if prob_r.status != "optimal":
			raise Exception, "Solver could not find optimal value"

		# Get optimal value of R
		r = np.array(R.value).reshape(-1)

		# Use optimal value of R to construct matrix H and vector Pdiff
		# TODO: separate the routine that gets pInit
		H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, colNamesH = makeH(sim_data, colNamesG, pPromoterBound, r, fixedTFs, cellSpecs)
		Pdiff = makePdiff(sim_data, colNamesH)

		# On first iteration, save the value of the initial p
		if i == 0:
			pInit0 = pInit.copy()

		# Optimize P such that the RNA polymerase initiation probabilities
		# computed from current values of R in matrix H are close to fitted
		# values.
		P = Variable(H.shape[1])

		# Boolean vector marking columns of H that correspond to alpha's and fixed TFs
		D = np.zeros(H.shape[1])
		D[pAlphaIdxs] = 1
		D[fixedTFIdxs] = 1

		# Initial p masked by D
		Drhs = pInit0.copy()
		Drhs[D != 1] = 0

		# Objective: minimize difference between k (fitted RNAP initiation
		# probabilities) and H*P (computed initiation probabilities) while
		# also minimizing deviation of P from the original value calculated
		# from mean TF and ligand concentrations
		objective_p = Minimize(norm(H*(SCALING*P) - SCALING*k, NORM_TYPE) + MOMENTUM_COEFF*norm(P - pInit0, NORM_TYPE))

		# Constraints
		# 1) 0 <= P <= 1 : All DNA-bound probabilities should be between zero
		# and one.
		# 2) D*P == Drhs : Values of P that correspond to alpha's and fixed TFs
		# should not change.
		# 3) Pdiff*P >= 0.1 : There must be at least a difference of 0.1
		# between binding probabilities of a TF in conditions TF__active and
		# TF__inactive
		# TODO: 0.1 should also be a parameter
		constraint_p = [0 <= SCALING*P, SCALING*P <= SCALING*1, np.diag(D)*(SCALING*P) == SCALING*Drhs, Pdiff*(SCALING*P) >= SCALING*0.1]

		# Solve optimization problem
		prob_p = Problem(objective_p, constraint_p)
		prob_p.solve(solver = "GLPK")

		if prob_p.status != "optimal":
			raise Exception, "Solver could not find optimal value"

		# Get optimal value of P
		p = np.array(P.value).reshape(-1)
		fromArray(p, pPromoterBound, pPromoterBoundIdxs)  # Update pPromoterBound with fitted p

		# Break from loop if parameters have converged
		if np.abs(np.linalg.norm(np.dot(H, p) - k, NORM_TYPE) - lastNorm) < 1e-9:
			break
		else:
			lastNorm = np.linalg.norm(np.dot(H, p) - k, NORM_TYPE)

	# Update sim_data with fitted bound probabilities and RNAP initiation
	# probabilities computed from these bound probabilities
	sim_data.pPromoterBound = pPromoterBound
	updateSynthProb(sim_data, kInfo, np.dot(H, p))

	# TODO: This should be a separate function - fitLigandConcentration
	cellDensity = sim_data.constants.cellDensity

	for tf in sorted(sim_data.tfToActiveInactiveConds):
		# Skip TFs that are not 1CS or are linked to genotypic perturbations
		if sim_data.process.transcription_regulation.tfToTfType[tf] != "1CS":
			continue
		if len(sim_data.tfToActiveInactiveConds[tf]["active genotype perturbations"]) > 0 or len(sim_data.tfToActiveInactiveConds[tf]["inactive genotype perturbations"]) > 0:
			continue

		activeKey = tf + "__active"
		inactiveKey = tf + "__inactive"

		# Determine if metabolite-bound form of the TF is the active form
		boundId = sim_data.process.transcription_regulation.activeToBound[tf]
		negativeSignal = (tf != boundId)  # True if unbound form is the active TF

		# Calculate kd of bound TF
		fwdRate = sim_data.process.equilibrium.getFwdRate(boundId + "[c]")
		revRate = sim_data.process.equilibrium.getRevRate(boundId + "[c]")
		kd = revRate/fwdRate

		# Get the metabolite that binds to the TF and its stoich coefficient
		metabolite = sim_data.process.equilibrium.getMetabolite(boundId + "[c]")
		metaboliteCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(boundId + "[c]")

		# Calculate the concentrations of the metabolite under conditions where
		# TF is active and inactive
		activeCellVolume = cellSpecs[activeKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		activeCountsToMolar = 1 / (sim_data.constants.nAvogadro * activeCellVolume)
		activeSignalConc = (activeCountsToMolar * cellSpecs[activeKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol/units.L)
		inactiveCellVolume = cellSpecs[inactiveKey]["avgCellDryMassInit"] / cellDensity / sim_data.mass.cellDryMassFraction
		inactiveCountsToMolar = 1 / (sim_data.constants.nAvogadro * inactiveCellVolume)
		inactiveSignalConc = (inactiveCountsToMolar * cellSpecs[inactiveKey]["bulkAverageContainer"].count(metabolite)).asNumber(units.mol/units.L)

		# Update kd with fitted values of P and the bulk average concentrations
		# of the metabolite, and use this fitted kd to recalculate the set
		# amounts of the metabolite in metabolism
		if negativeSignal:
			P = pPromoterBound[activeKey][tf]
			if 1 - P < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = (activeSignalConc**metaboliteCoeff) * P/(1 - P)

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (kdNew*(1 - P)/P)**(1./metaboliteCoeff)*(units.mol/units.L)

		else:
			P = pPromoterBound[inactiveKey][tf]
			if P < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = (inactiveSignalConc**metaboliteCoeff) * (1 - P)/P

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (kdNew*P/(1 - P))**(1./metaboliteCoeff)*(units.mol/units.L)

		# Fit reverse rate in line with fitted kd
		sim_data.process.equilibrium.setRevRate(boundId + "[c]", kdNew*fwdRate)

	return r


def calculatePromoterBoundProbability(sim_data, cellSpecs):
	"""
	Calculate the probability that a transcription factor is bound to its
	associated promoter for all simulated growth conditions. The bulk
	average concentrations calculated for TFs and their ligands are used to
	compute the probabilities based on the type (0CS, 1CS, 2CS) of the TF.

	Requires
	--------
	- Bulk average counts of transcription factors and associated ligands
	for each condition (in cellSpecs)

	Returns
	--------
	- pPromoterBound: Probability that a transcription factor is bound to
	its promoter, per growth condition and TF. Each probability is indexed by
	pPromoterBound[condition][TF].
	"""
	pPromoterBound = {}  # Initialize return value
	cellDensity = sim_data.constants.cellDensity

	for conditionKey in sorted(cellSpecs):
		pPromoterBound[conditionKey] = {}

		cellVolume = cellSpecs[conditionKey]["avgCellDryMassInit"]/cellDensity/sim_data.mass.cellDryMassFraction
		countsToMolar = 1/(sim_data.constants.nAvogadro*cellVolume)

		for tf in sorted(sim_data.tfToActiveInactiveConds):
			tfType = sim_data.process.transcription_regulation.tfToTfType[tf]

			if tfType == "0CS":
				tfCount = cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")

				if tfCount > 0:
					pPromoterBound[conditionKey][tf] = 1.  # If TF exists, the promoter is always bound to the TF
				else:
					pPromoterBound[conditionKey][tf] = 0.

			elif tfType == "1CS":
				boundId = sim_data.process.transcription_regulation.activeToBound[tf]  # ID of TF bound to ligand
				kd = sim_data.process.equilibrium.getRevRate(boundId + "[c]")/sim_data.process.equilibrium.getFwdRate(boundId + "[c]")

				signal = sim_data.process.equilibrium.getMetabolite(boundId + "[c]")  # ID of ligand that binds to TF
				signalCoeff = sim_data.process.equilibrium.getMetaboliteCoeff(boundId + "[c]")  # Stoichiometric coefficient of ligand

				# Get bulk average concentrations of ligand and TF
				signalConc = (countsToMolar*cellSpecs[conditionKey]["bulkAverageContainer"].count(signal)).asNumber(units.mol/units.L)
				tfConc = (countsToMolar*cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")).asNumber(units.mol/units.L)

				# If TF is active in its bound state
				if tf == boundId:
					if tfConc > 0:
						pPromoterBound[conditionKey][tf] = sim_data.process.transcription_regulation.pPromoterBoundSKd(signalConc, kd, signalCoeff)
					else:
						pPromoterBound[conditionKey][tf] = 0.

				# If TF is active in its unbound state
				else:
					if tfConc > 0:
						pPromoterBound[conditionKey][tf] = 1. - sim_data.process.transcription_regulation.pPromoterBoundSKd(signalConc, kd, signalCoeff)
					else:
						pPromoterBound[conditionKey][tf] = 0.

			elif tfType == "2CS":
				# Get bulk average concentrations of active and inactive TF
				activeTfConc = (countsToMolar*cellSpecs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")).asNumber(units.mol/units.L)
				inactiveTf = sim_data.process.two_component_system.activeToInactiveTF[tf + "[c]"]
				inactiveTfConc = (countsToMolar*cellSpecs[conditionKey]["bulkAverageContainer"].count(inactiveTf)).asNumber(units.mol/units.L)

				if activeTfConc == 0 and inactiveTfConc == 0:
					pPromoterBound[conditionKey][tf] = 0.
				else:
					pPromoterBound[conditionKey][tf] = activeTfConc/(activeTfConc + inactiveTfConc)

	return pPromoterBound

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
	"""
	Fits the affinities (Michaelis-Menten constants) for RNAs binding to endoRNAses.

	EndoRNAses perform the first step of RNA decay by cleaving a whole RNA somewhere inside its
	extent.  This results in RNA fragments, which are then digested into monomers by exoRNAses. To
	model endoRNAse activity, we need to determine an affinity (Michaelis-Menten constant) for each
	RNA that is consistent with experimentally observed half-lives.  The Michaelis-Menten constants
	must be determined simultaneously, as the RNAs must compete for the active site of the
	endoRNAse.  (See the RnaDegradation Process class for more information about the dynamical
	model.)  The parameters are estimated using a root solver (scipy.optimize.fsolve).  (See the
	sim_data.process.rna_decay.kmLossFunction method for more information about the optimization
	problem.)

	Requires
	--------
	- cell density, dry mass fraction, and average initial dry mass
		Used to calculate the cell volume, which in turn is used to calculate concentrations.
	- observed RNA degradation rates (half-lives)
	- enoRNAse counts
	- endoRNAse catalytic rate constants
	- RNA counts
	- boolean options that enable sensitivity analyses (see Notes below)

	Modifies
	--------
	- Michaelis-Menten constants for first-order decay
		TODO (John): Determine the purpose of these values - legacy?
	- Several optimization-related values
		Sensitivity analyses (optional, see Notes below)
		Terminal values for optimization-related functions

	Returns
	-------
	- enoRNAse Km values, in units of M

	Notes
	-----
	If certain options are set, a sensitivity analysis will be performed using a range of
	metaparameters. TODO (John): Determine default behavior.

	Outputs will be cached and utilized instead of running the optimization if possible.
	TODO (John): Determine if caching is functional, and consider removing functionality.

	The function that generates the optimization functions is defined under sim_data but has no
	dependency on sim_data, and therefore could be moved here or elsewhere. (TODO)

	TODO (John): Refactor as a pure function.

	TODO (John): Why is this function called 'cooperative'?  It seems to instead assume and model
		competitive binding.

	TODO (John): Determine what part (if any) of the 'linear' parameter fitting should be retained.
	"""

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
