"""
The parca, aka parameter calculator.

TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)
TODO: functionalize so that values are not both set and returned from some methods
"""

from __future__ import absolute_import, division, print_function

import os
import multiprocessing as mp
import sys
import traceback
from typing import Callable, List

from arrow import StochasticSystem
from cvxpy import Variable, Problem, Minimize, norm
import numpy as np
import scipy.optimize
import six
from six.moves import cPickle, range, zip

from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils import filepath, parallelization, units
from wholecell.utils.fitting import normalize, masses_and_counts_for_homeostatic_target


# Tweaks
# Adjustments to get protein expression for certain enzymes required for metabolism
TRANSLATION_EFFICIENCIES_ADJUSTMENTS = {
	"ADCLY-MONOMER[c]": 5,  # pabC, aminodeoxychorismate lyase
	"EG12438-MONOMER[c]": 5,  # menH, 2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthetase
	"EG12298-MONOMER[p]": 5,  # yibQ, Predicted polysaccharide deacetylase; This RNA is fit for the anaerobic condition viability
	"ACETYL-COA-ACETYLTRANSFER-MONOMER[c]": 5,  # atoB; This RNA is fit for the anaerobic condition viability
	}
RNA_EXPRESSION_ADJUSTMENTS = {
	"EG11493_RNA[c]": 10,  # pabC, aminodeoxychorismate lyase
	"EG12438_RNA[c]": 10,  # menH, 2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthetase
	"EG12298_RNA[c]": 10,  # yibQ, Predicted polysaccharide deacetylase; This RNA is fit for the anaerobic condition viability
	"EG11672_RNA[c]": 10,  # atoB, acetyl-CoA acetyltransferase; This RNA is fit for the anaerobic condition viability
	"EG10238_RNA[c]": 10,  # dnaE, DNA polymerase III subunit alpha; This RNA is fit for the sims to produce enough DNAPs for timely replication
	"EG11673_RNA[c]": 10,  # folB, dihydroneopterin aldolase; needed for growth (METHYLENE-THF) in acetate condition
	"EG10808_RNA[c]": 4,  # pyrE, orotate phosphoribosyltransferase; Needed for UTP synthesis, transcriptional regulation by UTP is not included in the model
	}
RNA_DEG_RATES_ADJUSTMENTS = {
	"EG11493_RNA[c]": 2,  # pabC, aminodeoxychorismate lyase
	"EG10709_RNA[c]": 2,  # pheS, phenylalanine synthetase subunit; for tRNA charging in anaerobic condition
	"EG10710_RNA[c]": 2,  # pheT, phenylalanine synthetase subunit; for tRNA charging in anaerobic condition
	}
PROTEIN_DEG_RATES_ADJUSTMENTS = {
	"ADENYLATECYC-MONOMER[c]": 2. / 600,  # CyaA, adenylate cyclase; convert from 2 min to 10 hr half life to get expression in acetate condition (required for cAMP)
	"SPOT-MONOMER[c]": 2. / 600,  # SpoT, ppGpp phosphatase; convert from 2 min to 10 hr half life to better match expected protein counts
	"EG12298-MONOMER[p]": 0.1, # yibQ, Predicted polysaccharide deacetylase; This protein is fit for the anaerobic condition
	}

# Fitting parameters
FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100
N_SEEDS = 10

# Parameters used in fitPromoterBoundProbability()
PROMOTER_PDIFF_THRESHOLD = 0.1  # Minimum difference between binding probabilities of a TF in conditions where TF is active and inactive
PROMOTER_REG_COEFF = 1e-3  # Optimization weight on how much probability should stay close to original values
PROMOTER_SCALING = 10  # Multiplied to all matrices for numerical stability
PROMOTER_NORM_TYPE = 1  # Matrix 1-norm
PROMOTER_MAX_ITERATIONS = 100
PROMOTER_CONVERGENCE_THRESHOLD = 1e-9
ECOS_0_TOLERANCE = 1e-12  # Tolerance to adjust solver output to 0

BASAL_EXPRESSION_CONDITION = "M9 Glucose minus AAs"

VERBOSE = 1

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s


def fitSimData_1(
		raw_data,
		cpus=1,
		debug=False,
		variable_elongation_transcription=False,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False):
	"""
	Fits parameters necessary for the simulation based on the knowledge base

	Inputs:
		raw_data (KnowledgeBaseEcoli) - knowledge base consisting of the
			necessary raw data
		cpus (int) - number of processes to use (if >1 uses multiprocessing)
		debug (bool) - if True, fit only one arbitrarily-chosen transcription
			factor in order to speed up a debug cycle (should not be used for
			an actual simulation)
		variable_elongation_transcription (bool) - enable variable elongation
			for transcription
		variable_elongation_translation (bool) - enable variable elongation for
			translation
		disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
			is not fit to protein synthesis demands
		disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
			expression is not fit to protein synthesis demands
	"""

	sim_data = SimulationDataEcoli()
	sim_data.initialize(
		raw_data = raw_data,
		basal_expression_condition = BASAL_EXPRESSION_CONDITION,
		)

	# Limit the number of conditions that are being fit so that execution time decreases
	if debug:
		print("Warning: Running the Parca in debug mode - not all conditions will be fit")
		key = list(sim_data.tfToActiveInactiveConds.keys())[0]
		sim_data.tfToActiveInactiveConds = {key: sim_data.tfToActiveInactiveConds[key]}

	# Make adjustments for metabolic enzymes
	setTranslationEfficiencies(sim_data)
	setRNAExpression(sim_data)
	setRNADegRates(sim_data)
	setProteinDegRates(sim_data)

	# Set C-period
	setCPeriod(sim_data)

	cellSpecs = buildBasalCellSpecifications(
		sim_data,
		disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting
		)

	# Set expression based on ppGpp regulation from basal expression
	sim_data.process.transcription.set_ppgpp_expression(sim_data)
	# TODO (Travis): use ppGpp expression in condition fitting below

	# Modify other properties

	# Re-compute Km's
	if sim_data.constants.EndoRNaseCooperation:
		sim_data.process.transcription.rnaData["KmEndoRNase"] = setKmCooperativeEndoRNonLinearRNAdecay(sim_data, cellSpecs["basal"]["bulkContainer"])

	## Calculate and set maintenance values

	# ----- Growth associated maintenance -----
	fitMaintenanceCosts(sim_data, cellSpecs["basal"]["bulkContainer"])

	# NOTE: multiprocessing `fork` seems to work here even on macOS, so override the
	# cpus() safety check for now. Be careful calling native libraries that
	# use threads and other resources which don't play well with `fork`.
	# See Issue #392.
	cpus = parallelization.cpus(cpus, advice='mac override')

	# Apply updates to cellSpecs from buildTfConditionCellSpecifications for each TF condition
	conditions = list(sorted(sim_data.tfToActiveInactiveConds))
	args = [(sim_data, tf, disable_ribosome_capacity_fitting, disable_rnapoly_capacity_fitting)
		for tf in conditions]
	apply_updates(buildTfConditionCellSpecifications, args, conditions, cellSpecs, cpus)

	for conditionKey in cellSpecs:
		if conditionKey == "basal":
			continue

		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

	buildCombinedConditionCellSpecifications(
		sim_data,
		cellSpecs,
		variable_elongation_transcription,
		variable_elongation_translation,
		disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting)

	sim_data.process.transcription.rnaSynthProbFraction = {}
	sim_data.process.transcription.rnapFractionActiveDict = {}
	sim_data.process.transcription.rnaSynthProbRProtein = {}
	sim_data.process.transcription.rnaSynthProbRnaPolymerase = {}
	sim_data.process.transcription.rnaPolymeraseElongationRateDict = {}
	sim_data.expectedDryMassIncreaseDict = {}
	sim_data.process.translation.ribosomeElongationRateDict = {}
	sim_data.process.translation.ribosomeFractionActiveDict = {}

	# Apply updates from fitCondition to cellSpecs for each fit condition
	conditions = list(sorted(cellSpecs))
	args = [(sim_data, cellSpecs[condition], condition)
		for condition in conditions]
	apply_updates(fitCondition, args, conditions, cellSpecs, cpus)

	for condition_label in sorted(cellSpecs):
		nutrients = sim_data.conditions[condition_label]["nutrients"]
		if nutrients not in sim_data.translationSupplyRate:
			sim_data.translationSupplyRate[nutrients] = cellSpecs[condition_label]["translation_aa_supply"]

	if VERBOSE > 0:
		print('Fitting promoter binding')
	# noinspection PyTypeChecker
	rVector = fitPromoterBoundProbability(sim_data, cellSpecs)
	# noinspection PyTypeChecker
	fitLigandConcentrations(sim_data, cellSpecs)

	# Adjust ppGpp regulated expression after conditions have been fit for physiological constraints
	sim_data.process.transcription.adjust_polymerizing_ppgpp_expression(sim_data)

	for condition_label in sorted(cellSpecs):
		condition = sim_data.conditions[condition_label]
		nutrients = condition["nutrients"]

		if VERBOSE > 0:
			print("Updating mass in condition {}".format(condition_label))
		spec = cellSpecs[condition_label]

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(nutrients)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[condition_label]))

		avgCellDryMassInit, fitAvgSolublePoolMass = rescaleMassForSolubleMetabolites(
			sim_data, spec["bulkContainer"], concDict, sim_data.conditionToDoublingTime[condition_label]
			)

		if VERBOSE > 0:
			print('{} to {}'.format(spec["avgCellDryMassInit"], avgCellDryMassInit))

		spec["avgCellDryMassInit"] = avgCellDryMassInit
		spec["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass

		mRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isMRna"]].sum()
		tRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isTRna"]].sum()
		rRnaSynthProb = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRRna"]].sum()

		if len(condition["perturbations"]) == 0:
			if nutrients not in sim_data.process.transcription.rnaSynthProbFraction:
				sim_data.process.transcription.rnaSynthProbFraction[nutrients] = {
					"mRna": mRnaSynthProb,
					"tRna": tRnaSynthProb,
					"rRna": rRnaSynthProb,
					}

			if nutrients not in sim_data.process.transcription.rnaSynthProbRProtein:
				prob = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRProtein"]]
				sim_data.process.transcription.rnaSynthProbRProtein[nutrients] = prob

			if nutrients not in sim_data.process.transcription.rnaSynthProbRnaPolymerase:
				prob = sim_data.process.transcription.rnaSynthProb[condition_label][sim_data.process.transcription.rnaData["isRnap"]]
				sim_data.process.transcription.rnaSynthProbRnaPolymerase[nutrients] = prob

			if nutrients not in sim_data.process.transcription.rnapFractionActiveDict:
				frac = sim_data.growthRateParameters.getFractionActiveRnap(spec["doubling_time"])
				sim_data.process.transcription.rnapFractionActiveDict[nutrients] = frac

			if nutrients not in sim_data.process.transcription.rnaPolymeraseElongationRateDict:
				rate = sim_data.growthRateParameters.getRnapElongationRate(spec["doubling_time"])
				sim_data.process.transcription.rnaPolymeraseElongationRateDict[nutrients] = rate

			if nutrients not in sim_data.expectedDryMassIncreaseDict:
				sim_data.expectedDryMassIncreaseDict[nutrients] = spec["avgCellDryMassInit"]

			if nutrients not in sim_data.process.translation.ribosomeElongationRateDict:
				rate = sim_data.growthRateParameters.getRibosomeElongationRate(spec["doubling_time"])
				sim_data.process.translation.ribosomeElongationRateDict[nutrients] = rate

			if nutrients not in sim_data.process.translation.ribosomeFractionActiveDict:
				frac = sim_data.growthRateParameters.getFractionActiveRibosome(spec["doubling_time"])
				sim_data.process.translation.ribosomeFractionActiveDict[nutrients] = frac

	calculateRnapRecruitment(sim_data, rVector)
	sim_data.process.metabolism.set_supply_constants(sim_data)

	return sim_data

def apply_updates(func, args, labels, dest, cpus):
	# type: (Callable[..., dict], List[tuple], List[str], dict, int) -> None
	"""
	Use multiprocessing (if cpus > 1) to apply args to a function to get
	dictionary updates for a destination dictionary.

	Args:
		func: function to call with args
		args: list of args to apply to func
		labels: label for each set of args for exception information
		dest: destination dictionary that will be updated with results
			from each function call
		cpus: number of cpus to use
	"""

	if cpus > 1:
		print("Starting {} Parca processes".format(cpus))

		# Apply args to func
		pool = mp.Pool(processes=cpus)
		results = {
			label: pool.apply_async(func, a)
			for label, a in zip(labels, args)
			}
		pool.close()
		pool.join()

		# Check results from function calls and update dest
		failed = []
		for label, result in results.items():
			if result.successful():
				dest.update(result.get())
			else:
				# noinspection PyBroadException
				try:
					result.get()
				except Exception as e:
					traceback.print_exc()
					failed.append(label)

		# Cleanup
		if failed:
			raise RuntimeError('Error(s) raised for {} while using multiple processes'
				.format(', '.join(failed)))
		pool = None
		print("End parallel processing")
	else:
		for a in args:
			dest.update(func(*a))

def buildBasalCellSpecifications(
		sim_data,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False
		):
	"""
	Creates cell specifications for the basal condition by fitting expression.
	Relies on expressionConverge() to set the expression and update masses.

	Inputs
	------
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

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
		'expression' (array of floats) - expression for each RNA, total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (array of floats) - synthesis probability for each RNA,
			total normalized to 1
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - the adjusted dry mass
			of the soluble fraction of a cell
		'bulkContainer' (BulkObjectsContainer object) - expected counts for
			bulk molecules based on expression

	Notes
	-----
	- TODO - sets sim_data attributes and returns values - change to only return values
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
		disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting
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

@parallelization.full_traceback
def buildTfConditionCellSpecifications(
		sim_data,
		tf,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False
		):
	"""
	Creates cell specifications for a given transcription factor by
	fitting expression. Will set for the active and inactive TF condition.
	Relies on expressionConverge() to set the expression and masses.
	Uses fold change data relative to the 'basal' condition to determine
	expression for a given TF.

	Inputs
	------
	- tf (str) - label for the transcription factor to fit (eg. 'CPLX-125')
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

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
		'expression' (array of floats) - expression for each RNA, total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (array of floats) - synthesis probability for each RNA,
			total normalized to 1
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - the adjusted dry mass
			of the soluble fraction of a cell
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
			for key, value in six.viewitems(fcDataTmp):
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
			disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
			disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting
			)

		# Store calculated values
		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

	return cellSpecs

def buildCombinedConditionCellSpecifications(
		sim_data,
		cellSpecs,
		variable_elongation_transcription=False,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False):
	"""
	Creates cell specifications for sets of transcription factors being active.
	These sets include conditions like 'with_aa' or 'no_oxygen' where multiple
	transcription factors will be active at the same time.

	Inputs
	------
	- cellSpecs {condition (str): dict} - information about each individual
	transcription factor condition
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

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
	- TODO - determine how to handle fold changes when multiple TFs change the
	same gene because multiplying both fold changes together might not be
	appropriate
	"""

	for conditionKey in sim_data.conditionActiveTfs:
		# Skip adjustments if 'basal' condition
		if conditionKey == "basal":
			continue

		# Get expression from fold changes for each TF in the given condition
		fcData = {}
		conditionValue = sim_data.conditions[conditionKey]
		for tf in sim_data.conditionActiveTfs[conditionKey]:
			for gene, fc in sim_data.tfToFC[tf].items():
				fcData[gene] = fcData.get(gene, 1) * fc

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
			variable_elongation_transcription = variable_elongation_transcription,
			variable_elongation_translation = variable_elongation_translation,
			disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
			disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting)

		# Modify cellSpecs for calculated values
		cellSpecs[conditionKey]["expression"] = expression
		cellSpecs[conditionKey]["synthProb"] = synthProb
		cellSpecs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cellSpecs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cellSpecs[conditionKey]["bulkContainer"] = bulkContainer

		# Modify sim_data expression
		sim_data.process.transcription.rnaExpression[conditionKey] = cellSpecs[conditionKey]["expression"]
		sim_data.process.transcription.rnaSynthProb[conditionKey] = cellSpecs[conditionKey]["synthProb"]

def expressionConverge(
		sim_data,
		expression,
		concDict,
		doubling_time,
		Km=None,
		variable_elongation_transcription=False,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False):
	"""
	Iteratively fits synthesis probabilities for RNA. Calculates initial
	expression based on gene expression data and makes adjustments to match
	physiological constraints for ribosome and RNAP counts. Relies on
	fitExpression() to converge

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- concDict {metabolite (str): concentration (float with units of mol/volume)} -
	dictionary for concentrations of each metabolite with location tag
	- doubling_time (float with units of time) - doubling time
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

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
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- fitAvgSolubleTargetMolMass (float with units of mass) - the adjusted dry mass
	of the soluble fraction of a cell
	- bulkContainer (BulkObjectsContainer object) - expected counts for
	bulk molecules based on expression
	"""

	if VERBOSE > 0:
		print("Fitting RNA synthesis probabilities.")

	for iteration in range(MAX_FITTING_ITERATIONS):
		if VERBOSE > 1:
			print('Iteration: {}'.format(iteration))

		initialExpression = expression.copy()
		expression = setInitialRnaExpression(sim_data, expression, doubling_time)

		bulkContainer = createBulkContainer(sim_data, expression, doubling_time)
		avgCellDryMassInit, fitAvgSolubleTargetMolMass = rescaleMassForSolubleMetabolites(sim_data, bulkContainer, concDict, doubling_time)

		if not disable_rnapoly_capacity_fitting:
			setRNAPCountsConstrainedByPhysiology(
				sim_data,
				bulkContainer,
				doubling_time,
				avgCellDryMassInit,
				variable_elongation_transcription,
				Km)

		if not disable_ribosome_capacity_fitting:
			setRibosomeCountsConstrainedByPhysiology(
				sim_data,
				bulkContainer,
				doubling_time,
				variable_elongation_translation)

		# Normalize expression and write out changes
		expression, synthProb = fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - expression)))
		if VERBOSE > 1:
			print('degree of fit: {}'.format(degreeOfFit))

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	return expression, synthProb, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict

@parallelization.full_traceback
def fitCondition(sim_data, spec, condition):
	"""
	Takes a given condition and returns the predicted bulk average, bulk deviation,
	protein monomer average, protein monomer deviation, and amino acid supply to
	translation. This relies on calculateBulkDistributions and calculateTranslationSupply.

	Inputs
	------
	- condition (str) - condition to fit (eg 'CPLX0-7705__active')
	- spec {property (str): property values} - cell specifications for the given condition.
	This function uses the specs "expression", "concDict", "avgCellDryMassInit",
	and "doubling_time"

	Returns
	--------
	- A dictionary {condition (str): spec (dict)} with the updated spec dictionary
	with the following values updated:
		- bulkAverageContainer (BulkObjectsContainer object) - the mean of the bulk counts
		- bulkDeviationContainer (BulkObjectsContainer object) - the standard deviation of the bulk counts
		- proteinMonomerAverageContainer (BulkObjectsContainer object) - the mean of the protein monomer counts
		- proteinMonomerDeviationContainer (BulkObjectsContainer object) - the standard deviation of the protein
		monomer counts
		- translation_aa_supply (array with units of mol/(mass.time)) - the supply rates
		for each amino acid to translation
	"""

	if VERBOSE > 0:
		print("Fitting condition {}".format(condition))

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

def calculateTranslationSupply(sim_data, doubling_time, bulkContainer, avgCellDryMassInit):
	"""
	Returns the supply rates of all amino acids to translation given the desired
	doubling time. This creates a limit on the polypeptide elongation process,
	and thus on growth. The amino acid supply rate is found by calculating the
	concentration of amino acids per gram dry cell weight and multiplying by the
	loss to dilution given doubling time.

	Inputs
	------
	- doubling_time (float with units of time) - measured doubling times given the condition
	- bulkContainer (BulkObjectsContainer object) - a container that tracks the counts of all bulk molecules
	- avgCellDryMassInit (float with units of mass) - the average initial cell dry mass

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

# Sub-fitting functions

def setTranslationEfficiencies(sim_data):
	"""
	This function's goal is to set translation efficiencies for a subset of metabolic proteins.
	It first gathers the index of the proteins it wants to modify, then changes the monomer
	translation efficiencies based on the adjustment that is specified.
	These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each protein that needs to be modified, it takes in a hard coded adjustment factor.

	Modifies
	--------
	- This function modifies, for a subset of proteins, their translational efficiencies in sim_data.
	It takes their current efficiency and multiplies them by the factor specified in adjustments.
	"""

	for protein in TRANSLATION_EFFICIENCIES_ADJUSTMENTS:
		idx = np.where(sim_data.process.translation.monomerData["id"] == protein)[0]
		sim_data.process.translation.translationEfficienciesByMonomer[idx] *= TRANSLATION_EFFICIENCIES_ADJUSTMENTS[protein]

def setRNAExpression(sim_data):
	"""
	This function's goal is to set expression levels for a subset of RNAs.
	It first gathers the index of the RNA's it wants to modify, then changes
	the expression levels of those RNAs, within sim_data, based on the
	specified adjustment factor.
	These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each RNA that needs to be modified, it takes in a hard coded
	adjustment factor.

	Modifies
	--------
	- This function modifies the basal RNA expression levels set in sim_data,
	for the chosen RNAs. It takes their current basal expression and multiplies
	them by the factor specified in adjustments.
	- After updating the basal expression levels for the given genes, the
	function normalizes all the basal expression levels.
	"""

	for rna in RNA_EXPRESSION_ADJUSTMENTS:
		idx = np.where(sim_data.process.transcription.rnaData["id"] == rna)[0]
		sim_data.process.transcription.rnaExpression["basal"][idx] *= RNA_EXPRESSION_ADJUSTMENTS[rna]

	sim_data.process.transcription.rnaExpression["basal"] /= sim_data.process.transcription.rnaExpression["basal"].sum()

def setRNADegRates(sim_data):
	"""
	This function's goal is to set the degradation rates for a subset of metabolic RNA's.
	It first gathers the index of the RNA's it wants to modify, then changes the degradation
	rates of those RNAs. These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each RNA that needs to be modified, it takes in a hard coded adjustment factor

	Modifies
	--------
	- This function modifies the RNA degradation rates for the chosen RNAs in sim_data.
	It takes their current degradation rate and multiplies them by the factor specified in adjustments.
	"""

	for rna in RNA_DEG_RATES_ADJUSTMENTS:
		idx = np.where(sim_data.process.transcription.rnaData["id"] == rna)[0]
		sim_data.process.transcription.rnaData.struct_array["degRate"][idx] *= RNA_DEG_RATES_ADJUSTMENTS[rna]

def setProteinDegRates(sim_data):
	"""
	This function's goal is to set the degradation rates for a subset of proteins.
	It first gathers the index of the proteins it wants to modify, then changes the degradation
	rates of those proteins. These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each protein that needs to be modified it take in a hard coded adjustment factor.

	Modifies
	--------
	- This function modifies the protein degradation rates for the chosen proteins in sim_data.
	It takes their current degradation rate and multiplies them by the factor specified in adjustments.
	"""

	for protein in PROTEIN_DEG_RATES_ADJUSTMENTS:
		idx = np.where(sim_data.process.translation.monomerData["id"] == protein)[0]
		sim_data.process.translation.monomerData.struct_array["degRate"][idx] *= PROTEIN_DEG_RATES_ADJUSTMENTS[protein]

def setCPeriod(sim_data):
	"""
	The C period is the time the cell takes to replicate its chromosome. This function calculates the C period
	based on knowledge of the length of the genome (in nucleotides) and the elongation rate.
	Dividing the genome length by the elongation rate alone will give the time to replicate that many nucleotides,
	this value is further divided by two since replication can take place in two directions.

	Requires
	--------
	- Genome length (nt) and the DNA polymerase elongation rate (nt/s).

	Modifies
	--------
	- This function modifies sim_data to contain the c_period.
	"""

	sim_data.growthRateParameters.c_period = sim_data.process.replication.genome_length * units.nt / sim_data.growthRateParameters.dnaPolymeraseElongationRate / 2
	sim_data.process.replication._c_period = sim_data.growthRateParameters.c_period.asNumber(units.min)

def rescaleMassForSolubleMetabolites(sim_data, bulkMolCntr, concDict, doubling_time):
	"""
	Adjust the cell's mass to accomodate target small molecule concentrations.

	Inputs
	------
	- bulkMolCntr (BulkObjectsContainer object) - a container that tracks the counts of all bulk molecules
	- concDict (dict) - a dictionary of metabolite ID (string) : concentration (unit'd number, dimensions of concentration) pairs
	- doubling_time (float with units of time) - measured doubling times given the condition

	Requires
	--------
	- Cell mass fraction data at a given doubling time.
	- Average cell density.
	- The conversion factor for transforming from the size of an average cell to the size of a cell
	  immediately following division.
	- Avogadro's number.
	- Concentrations of small molecules (including both dry mass components and water).

	Modifies
	--------
	- Adds small molecule counts to bulkMolCntr.

	Returns
	-------
	- newAvgCellDryMassInit, the adjusted dry mass of a cell immediately following division.
	- fitAvgSolubleTargetMolMass, the adjusted dry mass of the soluble fraction of a cell
	"""

	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)

	non_small_molecule_initial_cell_mass = (
		avgCellFractionMass["proteinMass"]
		+ avgCellFractionMass["rnaMass"]
		+ avgCellFractionMass["dnaMass"]
		) / sim_data.mass.avgCellToInitialCellConvFactor

	molar_units = units.mol / units.L

	targetMoleculeIds = sorted(concDict)
	targetMoleculeConcentrations = molar_units * np.array([
		concDict[key].asNumber(molar_units) for key in targetMoleculeIds
		]) # Have to strip and replace units to obtain the proper array data type

	assert np.all(targetMoleculeConcentrations.asNumber(molar_units) > 0), 'Homeostatic dFBA objective requires non-zero (positive) concentrations'

	molecular_weights = sim_data.getter.getMass(targetMoleculeIds)

	massesToAdd, countsToAdd = masses_and_counts_for_homeostatic_target(
		non_small_molecule_initial_cell_mass,
		targetMoleculeConcentrations,
		molecular_weights,
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		targetMoleculeIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculetargetMoleculesDryMass = units.hstack((
		massesToAdd[:targetMoleculeIds.index('WATER[c]')],
		massesToAdd[targetMoleculeIds.index('WATER[c]') + 1:]
		)) # remove water since it's not part of the dry mass

	newAvgCellDryMassInit = non_small_molecule_initial_cell_mass + units.sum(smallMoleculetargetMoleculesDryMass)
	fitAvgSolubleTargetMolMass = units.sum(smallMoleculetargetMoleculesDryMass) * sim_data.mass.avgCellToInitialCellConvFactor

	return newAvgCellDryMassInit, fitAvgSolubleTargetMolMass

def setInitialRnaExpression(sim_data, expression, doubling_time):
	"""
	Creates a container that with the initial count and ID of each RNA, calculated based on the mass fraction,
	molecular weight, and expression distribution of each RNA. For rRNA the counts are set based on mass, while for
	tRNA and mRNA the counts are set based on mass and relative abundance. Relies on the math function
	totalCountFromMassesAndRatios.

	Requires
	--------
	- Needs information from the knowledge base about the mass fraction, molecular weight, and distribution of each
	RNA species.

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- doubling_time (float with units of time) - doubling time for condition

	Returns
	--------
	- expression (array of floats) - contains the adjusted RNA expression,
	normalized to 1

	Notes
	-----
	- Now rnaData["synthProb"] does not match "expression"

	"""

	# Load from sim_data
	n_avogadro = sim_data.constants.nAvogadro
	rna_data = sim_data.process.transcription.rnaData
	get_average_copy_number = sim_data.process.replication.get_average_copy_number
	rna_mw = rna_data['mw']
	rna_coord = rna_data['replicationCoordinate']

	## Mask arrays for rRNAs
	is_rRNA23S = rna_data["isRRna23S"]
	is_rRNA16S = rna_data["isRRna16S"]
	is_rRNA5S = rna_data["isRRna5S"]
	is_tRNA = rna_data["isTRna"]
	is_mRNA = rna_data["isMRna"]

	## IDs
	ids_rnas = rna_data["id"]
	ids_rRNA23S = ids_rnas[is_rRNA23S]
	ids_rRNA16S = ids_rnas[is_rRNA16S]
	ids_rRNA5S = ids_rnas[is_rRNA5S]
	ids_mRNA = ids_rnas[is_mRNA]

	## Mass fractions
	initial_rna_mass = (sim_data.mass.getFractionMass(doubling_time)['rnaMass']
		/ sim_data.mass.avgCellToInitialCellConvFactor)
	ppgpp = sim_data.growthRateParameters.getppGppConc(doubling_time)
	rna_fractions = sim_data.process.transcription.get_rna_fractions(ppgpp)
	total_mass_rRNA23S = initial_rna_mass * rna_fractions['23S']
	total_mass_rRNA16S = initial_rna_mass * rna_fractions['16S']
	total_mass_rRNA5S = initial_rna_mass * rna_fractions['5S']
	total_mass_tRNA = initial_rna_mass * rna_fractions['trna']
	total_mass_mRNA = initial_rna_mass * rna_fractions['mrna']

	## Molecular weights
	individual_masses_rRNA23S = rna_mw[is_rRNA23S] / n_avogadro
	individual_masses_rRNA16S = rna_mw[is_rRNA16S] / n_avogadro
	individual_masses_rRNA5S = rna_mw[is_rRNA5S] / n_avogadro
	individual_masses_tRNA = rna_mw[is_tRNA] / n_avogadro
	individual_masses_mRNA = rna_mw[is_mRNA] / n_avogadro

	# Molecule expression distributions
	tau = doubling_time.asNumber(units.min)

	## Get replication coordinates of rRNA genes
	coord_rRNA23S = rna_coord[is_rRNA23S]
	coord_rRNA16S = rna_coord[is_rRNA16S]
	coord_rRNA5S = rna_coord[is_rRNA5S]

	## Get average copy numbers for all rRNA genes
	n_avg_copy_rRNA23S = get_average_copy_number(tau, coord_rRNA23S)
	n_avg_copy_rRNA16S = get_average_copy_number(tau, coord_rRNA16S)
	n_avg_copy_rRNA5S = get_average_copy_number(tau, coord_rRNA5S)

	# For rRNAs it is assumed that all operons have the same per-copy
	# transcription probabilities. Since the positions of the operons and thus
	# the average copy numbers of each rRNA gene are different, the
	# distribution is given as the normalized ratio of the average copy numbers
	# of each rRNA gene.
	distribution_rRNA23S = normalize(n_avg_copy_rRNA23S)
	distribution_rRNA16S = normalize(n_avg_copy_rRNA16S)
	distribution_rRNA5S = normalize(n_avg_copy_rRNA5S)

	trna_distribution = sim_data.mass.getTrnaDistribution(doubling_time)
	ids_tRNA = trna_distribution['id']
	distribution_tRNA = normalize(trna_distribution['molar_ratio_to_16SrRNA'])
	distribution_mRNA = normalize(expression[is_mRNA])

	# Construct bulk container
	rna_expression_container = BulkObjectsContainer(ids_rnas, dtype=np.float64)

	## Assign rRNA counts based on mass
	total_count_rRNA23S = totalCountFromMassesAndRatios(
		total_mass_rRNA23S,
		individual_masses_rRNA23S,
		distribution_rRNA23S
		)
	total_count_rRNA16S = totalCountFromMassesAndRatios(
		total_mass_rRNA16S,
		individual_masses_rRNA16S,
		distribution_rRNA16S
		)
	total_count_rRNA5S = totalCountFromMassesAndRatios(
		total_mass_rRNA5S,
		individual_masses_rRNA5S,
		distribution_rRNA5S
		)

	# Mass weighted average rRNA count to set rRNA subunit counts equal to each
	# other but keep the same expected total rRNA mass
	mass_weighting_rRNA23S = individual_masses_rRNA23S * distribution_rRNA23S
	mass_weighting_rRNA16S = individual_masses_rRNA16S * distribution_rRNA16S
	mass_weighting_rRNA5S = individual_masses_rRNA5S * distribution_rRNA5S
	total_count_rRNA_average = (
		units.sum(mass_weighting_rRNA23S * total_count_rRNA23S)
		+ units.sum(mass_weighting_rRNA16S * total_count_rRNA16S)
		+ units.sum(mass_weighting_rRNA5S * total_count_rRNA5S)
		) / (
		units.sum(mass_weighting_rRNA23S)
		+ units.sum(mass_weighting_rRNA16S)
		+ units.sum(mass_weighting_rRNA5S)
		)

	counts_rRNA23S = total_count_rRNA_average * distribution_rRNA23S
	counts_rRNA16S = total_count_rRNA_average * distribution_rRNA16S
	counts_rRNA5S = total_count_rRNA_average * distribution_rRNA5S

	rna_expression_container.countsIs(counts_rRNA23S, ids_rRNA23S)
	rna_expression_container.countsIs(counts_rRNA16S, ids_rRNA16S)
	rna_expression_container.countsIs(counts_rRNA5S, ids_rRNA5S)

	## Assign tRNA counts based on mass and relative abundances (see Dong 1996)
	total_count_tRNA = totalCountFromMassesAndRatios(
		total_mass_tRNA,
		individual_masses_tRNA,
		distribution_tRNA
		)
	counts_tRNA = total_count_tRNA * distribution_tRNA
	rna_expression_container.countsIs(counts_tRNA, ids_tRNA)

	## Assign mRNA counts based on mass and relative abundances (microarrays)
	total_count_mRNA = totalCountFromMassesAndRatios(
		total_mass_mRNA,
		individual_masses_mRNA,
		distribution_mRNA
		)
	counts_mRNA = total_count_mRNA * distribution_mRNA
	rna_expression_container.countsIs(counts_mRNA, ids_mRNA)

	expression = normalize(rna_expression_container.counts())

	return expression

def totalCountIdDistributionProtein(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of proteins from the relative expression of RNA,
	individual protein mass, and total protein mass. Relies on the math functions
	netLossRateFromDilutionAndDegradationProtein, proteinDistributionFrommRNA,
	totalCountFromMassesAndRatios.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA expression
	- doubling_time (float with units of time) - measured doubling time given the condition

	Returns
	--------
	- total_count_protein (float) - total number of proteins
	- ids_protein (array of str) - name of each protein with location tag
	- distribution_protein (array of floats) - distribution for each protein,
	normalized to 1
	"""

	ids_protein = sim_data.process.translation.monomerData["id"]
	total_mass_protein = sim_data.mass.getFractionMass(doubling_time)["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individual_masses_protein = sim_data.process.translation.monomerData["mw"] / sim_data.constants.nAvogadro
	distribution_transcripts_by_protein = normalize(expression[sim_data.relation.rnaIndexToMonomerMapping])
	translation_efficiencies_by_protein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	degradationRates = sim_data.process.translation.monomerData["degRate"]

	# Find the net protein loss
	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doubling_time, degradationRates)

	# Find the protein distribution
	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcripts_by_protein,
		translation_efficiencies_by_protein,
		netLossRate_protein
		)

	# Find total protein counts
	total_count_protein = totalCountFromMassesAndRatios(
		total_mass_protein,
		individual_masses_protein,
		distribution_protein
		)

	return total_count_protein, ids_protein, distribution_protein

def totalCountIdDistributionRNA(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of RNA from their relative expression, individual
	mass, and total RNA mass. Relies on the math function totalCountFromMassesAndRatios.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA expression
	- doubling_time (float with units of time) - measured doubling time given the condition

	Returns
	--------
	- total_count_RNA (float) - total number of RNAs
	- ids_rnas (array of str) - name of each RNA with location tag
	- distribution_RNA (array of floats) - distribution for each RNA,
	normalized to 1
	"""

	ids_rnas = sim_data.process.transcription.rnaData["id"]
	total_mass_RNA = sim_data.mass.getFractionMass(doubling_time)["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individual_masses_RNA = sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro

	distribution_RNA = normalize(expression)

	total_count_RNA = totalCountFromMassesAndRatios(
		total_mass_RNA,
		individual_masses_RNA,
		distribution_RNA
		)

	return total_count_RNA, ids_rnas, distribution_RNA

def createBulkContainer(sim_data, expression, doubling_time):
	"""
	Creates a container that tracks the counts of all bulk molecules. Relies on
	totalCountIdDistributionRNA and totalCountIdDistributionProtein to set the
	counts and IDs of all RNAs and proteins.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA expression
	- doubling_time (float with units of time) - measured doubling time given the condition

	Returns
	-------
	- bulkContainer (BulkObjectsContainer object) - a wrapper around a NumPy
	array that tracks the counts of bulk molecules
	"""

	total_count_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)

	total_count_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)
	ids_molecules = sim_data.internal_state.bulkMolecules.bulkData["id"]

	# Construct bulk container
	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	# Assign RNA counts based on mass and expression distribution
	counts_RNA = total_count_RNA * distribution_RNA
	bulkContainer.countsIs(counts_RNA, ids_rnas)

	# Assign protein counts based on mass and mRNA counts
	counts_protein = total_count_protein * distribution_protein
	bulkContainer.countsIs(counts_protein, ids_protein)

	return bulkContainer

def setRibosomeCountsConstrainedByPhysiology(
		sim_data,
		bulkContainer,
		doubling_time,
		variable_elongation_translation):
	"""
	Set counts of ribosomal subunits based on three constraints:
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal subunit counts based on RNA expression data

	Inputs
	------
	bulkContainer (BulkObjectsContainer object) - counts of bulk molecules
	doubling_time (float with units of time) - doubling time given the condition
	variable_elongation_translation (bool) - whether there is variable elongation for translation

	Modifies
	--------
	- counts of ribosomal protein subunits in bulkContainer
	"""

	active_fraction = sim_data.growthRateParameters.getFractionActiveRibosome(doubling_time)

	# Get IDs and stoichiometry of ribosome subunits
	ribosome30SSubunits = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)['subunitIds']
	ribosome50SSubunits = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)['subunitIds']
	ribosome30SStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)['subunitStoich']
	ribosome50SStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)['subunitStoich']

	# -- CONSTRAINT 1: Expected protein distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(sim_data.process.translation.monomerData['aaCounts'], axis = 1)
	proteinDegradationRates = sim_data.process.translation.monomerData["degRate"]
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomerData["id"])

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(
		doubling_time,
		proteinDegradationRates,
		)

	elongation_rates = sim_data.process.translation.make_elongation_rates(
		None,
		sim_data.growthRateParameters.getRibosomeElongationRate(doubling_time).asNumber(units.aa / units.s),
		1,
		variable_elongation_translation)

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
		proteinLengths,
		elongation_rates,
		netLossRate_protein,
		proteinCounts).asNumber(units.aa / units.s) / active_fraction

	# Minimum number of ribosomes needed
	constraint1_ribosome30SCounts = (
		nRibosomesNeeded * ribosome30SStoich
		)

	constraint1_ribosome50SCounts = (
		nRibosomesNeeded * ribosome50SStoich
		)


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


	# -- CONSTRAINT 3: Expected ribosomal subunit counts based expression
	## Calculate fundamental ribosomal subunit count distribution based on RNA expression data
	## Already calculated and stored in bulkContainer
	ribosome30SCounts = bulkContainer.counts(ribosome30SSubunits)
	ribosome50SCounts = bulkContainer.counts(ribosome50SSubunits)

	# -- SET RIBOSOME FUNDAMENTAL SUBUNIT COUNTS TO MAXIMUM CONSTRAINT -- #
	constraint_names = np.array(["Insufficient to double protein counts", "Too small for mass fraction", "Current level OK"])
	rib30lims = np.array([nRibosomesNeeded, massFracPredicted_30SCount, (ribosome30SCounts / ribosome30SStoich).min()])
	rib50lims = np.array([nRibosomesNeeded, massFracPredicted_50SCount, (ribosome50SCounts / ribosome50SStoich).min()])
	if VERBOSE > 1:
		print('30S limit: {}'.format(constraint_names[np.where(rib30lims.max() == rib30lims)[0]][-1]))
		print('30S actual count: {}'.format((ribosome30SCounts / ribosome30SStoich).min()))
		print('30S count set to: {}'.format(rib30lims[np.where(rib30lims.max() == rib30lims)[0]][-1]))
		print('50S limit: {}'.format(constraint_names[np.where(rib50lims.max() == rib50lims)[0]][-1]))
		print('50S actual count: {}'.format((ribosome50SCounts / ribosome50SStoich).min()))
		print('50S count set to: {}'.format(rib50lims[np.where(rib50lims.max() == rib50lims)[0]][-1]))

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome30SCounts, constraint1_ribosome30SCounts), constraint2_ribosome30SCounts),
		ribosome30SSubunits
		)

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome50SCounts, constraint1_ribosome50SCounts), constraint2_ribosome50SCounts),
		ribosome50SSubunits
		)

	# Return rRNA counts to value in sim_data
	bulkContainer.countsIs(rRna23SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna23S"]])
	bulkContainer.countsIs(rRna16SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna16S"]])
	bulkContainer.countsIs(rRna5SCounts, sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna5S"]])

def setRNAPCountsConstrainedByPhysiology(
		sim_data,
		bulkContainer,
		doubling_time,
		avgCellDryMassInit,
		variable_elongation_transcription,
		Km=None):
	"""
	Set counts of RNA polymerase based on two constraints:
	(1) Number of RNAP subunits required to maintain steady state of mRNAs
	(2) Expected RNAP subunit counts based on (mRNA) distribution recorded in
		bulkContainer

	Inputs
	------
	- bulkContainer (BulkObjectsContainer object) - counts of bulk molecules
	- doubling_time (float with units of time) - doubling time given the condition
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases

	Modifies
	--------
	- bulkContainer (BulkObjectsContainer object) - the counts of RNA polymerase
	subunits are set according to Constraint 1

	Notes
	-----
	- Constraint 2 is not being used -- see final line of this function.
	"""

	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(sim_data.process.transcription.rnaData['countsACGU'], axis = 1)

	rnaLossRate = None

	if Km is None:
		# RNA loss rate is in units of counts/time, and computed by summing the
		# contributions of degradation and dilution.
		rnaLossRate = netLossRateFromDilutionAndDegradationRNALinear(
			doubling_time,
			sim_data.process.transcription.rnaData["degRate"],
			bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		)
	else:
		# Get constants to compute countsToMolar factor
		cellDensity = sim_data.constants.cellDensity
		cellVolume = avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		# Gompute input arguments for netLossRateFromDilutionAndDegradationRNA()
		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

		# RNA loss rate is in units of counts/time, and computed by accounting
		# for the competitive inhibition of RNase by other RNA targets.
		rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
			doubling_time,
			(1 / countsToMolar) * totalEndoRnaseCapacity,
			Km,
			rnaConc,
			countsToMolar,
			)

	# Compute number of RNA polymerases required to maintain steady state of mRNA
	elongation_rates = sim_data.process.transcription.make_elongation_rates(
		None,
		sim_data.growthRateParameters.getRnapElongationRate(doubling_time).asNumber(units.nt / units.s),
		1,
		variable_elongation_transcription)

	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths,
		elongation_rates,
		rnaLossRate).asNumber(units.nt / units.s)

	nRnapsNeeded = nActiveRnapNeeded / sim_data.growthRateParameters.getFractionActiveRnap(doubling_time)

	# Convert nRnapsNeeded to the number of RNA polymerase subunits required
	rnapIds = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)['subunitIds']
	rnapStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)['subunitStoich']
	minRnapSubunitCounts = nRnapsNeeded * rnapStoich

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(rnapIds)

	## -- SET RNAP COUNTS TO MAXIMUM CONSTRAINTS -- #
	constraint_names = np.array(["Current level OK", "Insufficient to double RNA distribution"])
	rnapLims = np.array([(rnapCounts / rnapStoich).min(), (minRnapSubunitCounts / rnapStoich).min()])
	if VERBOSE > 1:
		print('rnap limit: {}'.format(constraint_names[np.where(rnapLims.max() == rnapLims)[0]][0]))
		print('rnap actual count: {}'.format((rnapCounts / rnapStoich).min()))
		print('rnap counts set to: {}'.format(rnapLims[np.where(rnapLims.max() == rnapLims)[0]][0]))

	if np.any(minRnapSubunitCounts < 0):
		raise ValueError('RNAP protein counts must be positive.')

	bulkContainer.countsIs(minRnapSubunitCounts, rnapIds)

def fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km=None):
	"""
	Determines expression and synthesis probabilities for RNA molecules to fit
	protein levels and RNA degradation rates. Assumes a steady state analysis
	where the RNA synthesis probability will be the same as the degradation rate.
	If no Km is given, then RNA degradation is assumed to be linear otherwise
	degradation is calculated based on saturation with RNases.

	Inputs
	------
	- bulkContainer (BulkObjectsContainer object) - expected counts for
	bulk molecules based on expression
	- doubling_time (float with units of time) - doubling time
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases

	Modifies
	--------
	- bulkContainer counts of RNA and proteins

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1
	- synthProb (array of floats) - synthesis probability for each RNA which
	accounts for expression and degradation rate, normalized to 1

	Notes
	-----
	- TODO - sets bulkContainer counts and returns values - change to only return values
	"""

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
		dryMassFraction = sim_data.mass.cellDryMassFraction
		cellVolume = avgCellDryMassInit / cellDensity / dryMassFraction
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

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
	- energetic (GTP) cost of translation (per amino acid polymerized)
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
	atp_per_charge = 2  # ATP -> AMP is explicitly used in charging reactions so can remove from GAM

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
	explicit_mmol_maintenance_per_gdcw = (atp_per_charge + gtpPerTranslation) * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		sim_data.constants.growthAssociatedMaintenance -
		explicit_mmol_maintenance_per_gdcw
		)

	# We do not want to create energy with growth by having a negative darkATP
	# value. GAM measurements have some error so it's possible explicit
	# accounting could be more accurate or the GAM value used is too low which
	# would lead to a negative value. Easy fix is setting darkATP = 0 if this
	# error is raised.
	if darkATP.asNumber() < 0:
		raise ValueError('GAM has been adjusted too low. Explicit energy accounting should not exceed GAM.'
			' Consider setting darkATP to 0 if energy corrections are accurate.')

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
	- N_SEEDS (int) - the number of instantiated cells

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- concDict {metabolite (str): concentration (float with units of mol/volume)} -
	dictionary for concentrations of each metabolite with location tag
	- avgCellDryMassInit (float with units of mass) - initial dry cell mass
	- doubling_time (float with units of time) - doubling time for condition

	Returns
	--------
	- bulkAverageContainer (BulkObjectsContainer object) - the mean of the bulk counts
	- bulkDeviationContainer (BulkObjectsContainer object) - the standard deviation of the bulk counts
	- proteinMonomerAverageContainer (BulkObjectsContainer object) - the mean of the protein monomer counts
	- proteinMonomerDeviationContainer (BulkObjectsContainer object) - the standard deviation of the protein
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
	# Data for equilibrium binding
	# equilibriumDerivatives = sim_data.process.equilibrium.derivatives
	# equilibriumDerivativesJacobian = sim_data.process.equilibrium.derivativesJacobian

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
		print("Bulk distribution seed:")

	# Instantiate cells to find average copy numbers of macromolecules
	for seed in range(N_SEEDS):
		if VERBOSE > 1:
			print('seed = {}'.format(seed))

		allMoleculesView.countsIs(0)

		# randomState = np.random.RandomState(seed)

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
		time_step = 2**31 # don't stop until all complexes are formed.
		complexation_rates = sim_data.process.complexation.rates
		system = StochasticSystem(complexationStoichMatrix.T, random_seed=seed)
		complexation_result = system.evolve(
			time_step, complexationMoleculeCounts, complexation_rates)

		updatedCompMoleculeCounts = complexation_result['outcome']
		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		metDiffs = np.inf * np.ones_like(metabolitesView.counts())
		nIters = 0

		# Iterate processes until metabolites converge to a steady-state
		while np.linalg.norm(metDiffs, np.inf) > 1:
			random_state = np.random.RandomState(seed)
			metCounts = conc_metabolites * cellVolume * sim_data.constants.nAvogadro
			metCounts.normalize()
			metCounts.checkNoUnit()
			metabolitesView.countsIs(
				metCounts.asNumber().round()
				)

			# Find reaction fluxes from equilibrium process
			# Do not use jit to avoid compiling time (especially when running
			# in parallel since sim_data needs to be pickled and reconstructed
			# each time)
			rxnFluxes, _ = sim_data.process.equilibrium.fluxesAndMoleculesToSS(
				equilibriumMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.nAvogadro.asNumber(1 / units.mol),
				random_state, jit=False,
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
				raise Exception("Equilibrium reactions are not converging!")

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	# Update counts in bulk objects container
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
	- counts (float): total counts (does not need to be a whole number)
	"""

	assert np.allclose(np.sum(distribution), 1)
	counts = 1 / units.dot(individualMasses, distribution) * totalMass
	return units.strip_empty_units(counts)

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

	Inputs
	------
	- distribution_mRNA (array of floats) - distribution for each mRNA,
	normalized to 1
	- translation_efficiencies (array of floats) - translational efficiency for each mRNA,
	normalized to 1
	- netLossRate (array of floats with units of 1/time) - rate of loss for each protein

	Returns
	-------
	- array of floats for the distribution of each protein, normalized to 1
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

	Normalizing M_i by summing over all i cancels out k and P_total
	assuming a constant translation rate.

	Inputs
	------
	- distribution_protein (array of floats) - distribution for each protein,
	normalized to 1
	- translation_efficiencies (array of floats) - translational efficiency for each mRNA,
	normalized to 1
	- netLossRate (array of floats with units of 1/time) - rate of loss for each protein

	Returns
	-------
	- array of floats for the distribution of each mRNA, normalized to 1
	"""

	assert np.allclose(np.sum(distribution_protein), 1)
	distributionUnnormed = netLossRate * distribution_protein / translation_efficiencies
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()

	return distributionNormed.asNumber()

def calculateMinPolymerizingEnzymeByProductDistribution(
		productLengths, elongationRates, netLossRate, productCounts):
	"""
	Compute the number of ribosomes required to maintain steady state.

	dP/dt = production rate - loss rate
	dP/dt = e_r * (1/L) * R - (k_loss * P)

	At steady state: dP/dt = 0
	R = sum over i ((L_i / e_r) * k_loss_i * P_i)

	Multiplying both sides by volume gives an equation in terms of counts.

	P = protein concentration
	e_r = polypeptide elongation rate per ribosome
	L = protein length
	R = ribosome concentration
	k_loss = net protein loss rate
	i = ith protein

	Inputs
	------
	- productLengths (array of ints with units of amino_acids) - L, protein lengths
	- elongationRates (array of ints with units of amino_acid/time) e_r, polypeptide elongation rate
	- netLossRate (array of floats with units of 1/time) - k_loss, protein loss rate
	- productCounts (array of floats) - P, protein counts

	Returns
	--------
	- float with dimensionless units for the number of ribosomes required to
	maintain steady state
	"""

	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRates
		* netLossRate
		* productCounts)
	return nPolymerizingEnzymeNeeded

def calculateMinPolymerizingEnzymeByProductDistributionRNA(productLengths, elongationRates, netLossRate):
	"""
	Compute the number of RNA polymerases required to maintain steady state of mRNA.

	dR/dt = production rate - loss rate
	dR/dt = e_r * (1/L) * RNAp - k_loss

	At steady state: dR/dt = 0
	RNAp = sum over i ((L_i / e_r) * k_loss_i)

	Multiplying both sides by volume gives an equation in terms of counts.

	R = mRNA transcript concentration
	e_r = transcript elongation rate per RNAp
	L = transcript length
	RNAp = RNAp concentration
	k_loss = net transcript loss rate (unit: concentration / time)
	i = ith transcript

	Inputs
	------
	- productLengths (array of ints with units of nucleotides) - L, transcript lengths
	- elongationRates (array of ints with units of nucleotide/time) - e_r, transcript elongation rate
	- netLossRate (array of floats with units of 1/time) - k_loss, transcript loss rate

	Returns
	-------
	- float with dimensionless units for the number of RNA polymerases required to
	maintain steady state
	"""

	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRates
		* netLossRate)
	return nPolymerizingEnzymeNeeded

def netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates):
	"""
	Compute total loss rate (summed contributions of degradation and dilution).

	Inputs
	------
	- doublingTime (float with units of time) - doubling time of the cell
	- degradationRates (array of floats with units of 1/time) - protein degradation rate

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each protein
	"""

	return np.log(2) / doublingTime + degradationRates

def netLossRateFromDilutionAndDegradationRNA(doublingTime, totalEndoRnaseCountsCapacity, Km, rnaConc, countsToMolar):
	"""
	Compute total loss rate (summed impact of degradation and dilution).
	Returns the loss rate in units of (counts/time) in preparation for use in
	the steady state analysis in fitExpression() and
	setRNAPCountsConstrainedByPhysiology()
	(see calculateMinPolymerizingEnzymeByProductDistributionRNA()).

	Derived from steady state analysis of Michaelis-Menten enzyme kinetics with
	competitive inhibition: for a given RNA, all other RNAs compete for RNase.

	V_i = k_cat * [ES_i]
	v_i = k_cat * [E]0 * ([S_i]/Km_i) / (1 + sum over j genes([S_j] / Km_j))

	Inputs
	------
	- doublingTime (float with units of time) - doubling time of the cell
	- totalEndoRnaseCountsCapacity (float with units of 1/time) total kinetic
	capacity of all RNases in the cell
	- Km (array of floats with units of mol/volume) - Michaelis-Menten constant
	for each RNA
	- rnaConc (array of floats with units of mol/volume) - concentration for each RNA
	- countsToMolar (float with units of mol/volume) - conversion between counts and molar

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each RNA
	"""

	fracSaturated = rnaConc / Km / (1 + units.sum(rnaConc / Km))
	rnaCounts = (1 / countsToMolar) * rnaConc

	return (np.log(2) / doublingTime) * rnaCounts + (totalEndoRnaseCountsCapacity * fracSaturated)

def netLossRateFromDilutionAndDegradationRNALinear(doublingTime, degradationRates, rnaCounts):
	"""
	Compute total loss rate (summed contributions of degradation and dilution).
	Returns the loss rate in units of (counts/time) in preparation for use in
	the steady state analysis in fitExpression() and
	setRNAPCountsConstrainedByPhysiology()
	(see calculateMinPolymerizingEnzymeByProductDistributionRNA()).

	Requires
	--------
	- doublingTime (float with units of time) - doubling time of the cell
	- degradationRates (array of floats with units of 1/time) - degradation rate
	for each RNA
	- rnaCounts (array of floats) - counts for each RNA

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each RNA
	"""

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

def fitPromoterBoundProbability(sim_data, cellSpecs):
	r"""
	Calculates the probabilities (P) that each transcription factor will bind
	to its target RNA. This function initially calculates these probabilities
	from the bulk average counts of the TFs and ligands calculated from
	previous steps. Then, values of parameters alpha and r in the equation
	below are fit such that the computed RNA synthesis probabilities converge
	to the measured RNA synthesis probabilities.

	v_{synth, j} = \alpha_j + \sum_{i} P_{T,i}*r_{ij}

	Due to constraints applied in the optimization, both v and P need to
	be shifted from their initial values.

	Requires
	--------
	- Bulk average counts of transcription factors and associated ligands
	for each condition (in cellSpecs)

	Inputs
	------
	- cellSpecs {condition (str): dict} - information about each condition

	Modifies
	--------
	- Probabilities of TFs binding to their promoters
	- RNA synthesis probabilities

	Returns
	--------
	- r: Fit parameters on how the recruitment of a TF affects the expression
	of a gene. High (positive) values of r indicate that the TF binding
	increases the probability that the gene is expressed.

	Notes
	--------
	See supplementary materials on transcription regulation for details on
	the parameters being fit.
	"""

	def build_vector_k(sim_data, cellSpecs):
		"""
		Construct vector k that contains existing fit transcription
		probabilities of RNAs in each relevant condition, normalized by the
		average copy number of the gene encoding the RNA while the cell grows
		in that condition.

		Returns
		--------
		- k: List of RNA synthesis probabilities for each RNA and condition,
		normalized by gene copy number.
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		"""

		k, kInfo = [], []

		for idx, (rnaId, rnaCoordinate) in enumerate(
				zip(sim_data.process.transcription.rnaData["id"],
				sim_data.process.transcription.rnaData["replicationCoordinate"])):
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			# TODO (Gwanggyu): cache this list of conditions for each RNA
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

				# Get specific doubling time for this condition
				tau = cellSpecs[condition]["doubling_time"].asNumber(units.min)

				# Calculate average copy number of gene for this condition
				n_avg_copy = sim_data.process.replication.get_average_copy_number(tau, rnaCoordinate)

				# Compute synthesis probability per gene copy
				prob_per_copy = sim_data.process.transcription.rnaSynthProb[condition][idx] / n_avg_copy

				# Gather RNA synthesis probabilities for each RNA per condition
				k.append(prob_per_copy)
				kInfo.append({"condition": condition, "idx": idx})

		k = np.array(k)

		return k, kInfo

	def build_matrix_G(sim_data, pPromoterBound):
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

		Inputs
		------
		- pPromoterBound: Probabilities that a given TF is bound to its
		promoter in a given condition, calculated from bulk average
		concentrations of the TF and its associated ligands.

		Returns
		--------
		- G: Matrix of values in pPromoterBound, rearranged based on each RNA
		- rowNames: List of row names of G as strings
		- colNames: List of column names of G as strings
		"""

		gI, gJ, gV, rowNames, colNames = [], [], [], [], []

		for idx, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
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
					# Add column for each TF that regulates each RNA
					colName = rnaIdNoLoc + "__" + tf

					# TODO (Gwanggyu): Are these checks necessary?
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
				gV.append(1.)

		gI, gJ, gV = np.array(gI), np.array(gJ), np.array(gV)
		G = np.zeros((len(rowNames), len(colNames)), np.float64)
		G[gI, gJ] = gV

		return G, rowNames, colNames

	def build_matrix_Z(sim_data, colNames):
		"""
		Construct matrix Z that connects all possible TF combinations with
		each TF. Each row of the matrix corresponds to an RNA-(TF combination)
		pair, and each column corresponds to an RNA-TF pair, with an additional
		RNA-alpha column for each RNA (identical to matrix G). Matrix values
		are set to one if the TF specified by the column is "active" in the
		combination specified by the row or if the column is an RNA-alpha
		column, and zero otherwise.

		Inputs
		------
		- colNames: List of column names from matrix G.

		Returns
		--------
		- Z: Matrix of zeros and ones, specifying which TFs in the columns
		correspond to combinations in the rows.
		"""

		combinationIdxToColIdxs = {
			0: [0], 1: [0, 1], 2: [0, 2], 3: [0, 1, 2],
			4: [0, 3], 5: [0, 1, 3], 6: [0, 2, 3], 7: [0, 1, 2, 3],
			8: [0, 4], 9: [0, 1, 4], 10: [0, 2, 4], 11: [0, 1, 2, 4],
			12: [0, 3, 4], 13: [0, 1, 3, 4], 14: [0, 2, 3, 4], 15: [0, 1, 2, 3, 4],
			}

		zI, zJ, zV, rowNames = [], [], [], []

		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

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
			for combinationIdx in range(2**nTfs):
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

	def build_matrix_T(sim_data, colNames):
		"""
		Construct matrix T that specifies the direction of regulation for each
		RNA-TF pair.

		Inputs
		------
		- colNames: List of column names from matrix G.

		Returns
		--------
		- T: Diagonal matrix. Diagonal value is +1 if the direction of
		regulation by the TF-RNA pair specified by the row is positive, -1 if
		this is negative, and 0 if the row is an RNA_alpha row.
		"""

		tI, tJ, tV, rowNamesT = [], [], [], []

		for rnaId in sim_data.process.transcription.rnaData["id"]:
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

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

	def build_matrix_H(sim_data, colNames, pPromoterBound, r, fixedTFs, cellSpecs):
		r"""
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

		Inputs
		------
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
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

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

	def build_matrix_pdiff(sim_data, colNamesH):
		"""
		Construct matrix Pdiff that specifies the indexes of corresponding
		TFs and conditions.

		Inputs
		------
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
		"""
		Updates values in pPromoterBound with fit probabilities.

		Inputs
		------
		- p: Vector of probabilities optimized in the current step.
		- pPromoterBoundIdxs: Dictionary of indexes to p

		Modifies
		--------
		Values in pPromoterBound - probabilities that each transcription factor
		is bound to its promoter for each growth condition
		"""

		for condition in sorted(pPromoterBoundIdxs):
			for tf in sorted(pPromoterBoundIdxs[condition]):
				pPromoterBound[condition][tf] = p[pPromoterBoundIdxs[condition][tf]]

	def updateSynthProb(sim_data, cellSpecs, kInfo, k):
		"""
		Updates RNA synthesis probabilities with fit values of P and R. The
		expected average copy number of genes for the condition are multiplied
		to the per-copy probabilities.

		Inputs
		------
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		- k: RNA synthesis probabilities computed from fit P and R.

		Modifies
		--------
		- RNA synthesis probabilities

		Notes
		--------
		These values are used to calculate pre-set probabilities for
		pools of mRNA, tRNA, and rRNA, and the genes that encode RNAP and
		ribosome subunits.
		"""

		# Get replication coordinates of each RNA
		replicationCoordinate = sim_data.process.transcription.rnaData["replicationCoordinate"]

		# Update sim_data values with fit values
		for D, k_value in zip(kInfo, k):
			condition = D["condition"]
			rna_idx = D["idx"]

			# Get coordinate of RNA
			rnaCoordinate = replicationCoordinate[rna_idx]

			# Get specific doubling time for this condition
			tau = cellSpecs[condition]["doubling_time"].asNumber(units.min)

			# Calculate average copy number of gene for this condition
			n_avg_copy = sim_data.process.replication.get_average_copy_number(tau, rnaCoordinate)

			sim_data.process.transcription.rnaSynthProb[condition][rna_idx] = k_value*n_avg_copy

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
		if (sim_data.process.transcription_regulation.tfToTfType[tf] == "1CS"
				and sim_data.tfToActiveInactiveConds[tf]["active nutrients"] ==
				sim_data.tfToActiveInactiveConds[tf]["inactive nutrients"]):
			fixedTFs.append(tf)

	# Build vector of existing fit transcription probabilities
	k, kInfo = build_vector_k(sim_data, cellSpecs)

	# Repeat for a fixed maximum number of iterations
	for i in range(PROMOTER_MAX_ITERATIONS):
		# Build matrices used in optimizing R
		G, rowNamesG, colNamesG = build_matrix_G(sim_data, pPromoterBound)
		Z = build_matrix_Z(sim_data, colNamesG)
		T = build_matrix_T(sim_data, colNamesG)

		# Optimize R such that transcription initiation probabilities computed
		# from existing values of P in matrix G are close to fit values.
		R = Variable(G.shape[1])  # Vector of r's and alpha's

		# Objective: minimize difference between k and G @ R
		objective_r = Minimize(
			norm(G @ (PROMOTER_SCALING * R) - PROMOTER_SCALING * k, PROMOTER_NORM_TYPE))

		# Optimization constraints
		# 1) 0 <= Z @ R <= 1 : Assuming P = 1 for all TFs, all possible
		# combinations of TFs should yield a valid transcription probability
		# value between zero and one.
		# 2) T @ R >= 0 : Values of r for positive regulation should be positive,
		# and values of r for negative regulation should be negative.
		constraint_r = [
			0 <= Z @ (PROMOTER_SCALING * R), Z @ (PROMOTER_SCALING * R) <= PROMOTER_SCALING,
			T @ (PROMOTER_SCALING * R) >= 0]

		# Solve optimization problem
		prob_r = Problem(objective_r, constraint_r)
		prob_r.solve(solver='ECOS')

		if prob_r.status != "optimal":
			raise Exception("Solver could not find optimal value")

		# Get optimal value of R
		r = np.array(R.value).reshape(-1)
		r[np.abs(r) < ECOS_0_TOLERANCE] = 0  # Adjust to 0 for small values from solver tolerance

		# Use optimal value of R to construct matrix H and vector Pdiff
		H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, colNamesH = build_matrix_H(
			sim_data, colNamesG, pPromoterBound, r, fixedTFs, cellSpecs)
		pdiff = build_matrix_pdiff(sim_data, colNamesH)

		# On first iteration, save the value of the initial p
		if i == 0:
			pInit0 = pInit.copy()

		# Optimize P such that the transcription probabilities computed from
		# current values of R in matrix H are close to fit values.
		P = Variable(H.shape[1])

		# Construct a boolean vector that marks column indexes of H
		# corresponding to alpha's and fixed TFs
		D = np.zeros(H.shape[1])
		D[pAlphaIdxs] = 1
		D[fixedTFIdxs] = 1

		# Mask initial p with boolean vector constructed above
		Drhs = pInit0.copy()
		Drhs[D != 1] = 0

		# Objective: minimize difference between k (fit RNAP initiation
		# probabilities) and H @ P (computed initiation probabilities) while
		# also minimizing deviation of P from the original value calculated
		# from mean TF and ligand concentrations
		objective_p = Minimize(
			norm(H @ (PROMOTER_SCALING * P) - PROMOTER_SCALING * k, PROMOTER_NORM_TYPE)
			+ PROMOTER_REG_COEFF * norm(P - pInit0, PROMOTER_NORM_TYPE))

		# Constraints
		# 1) 0 <= P <= 1 : All DNA-bound probabilities should be between zero
		# and one.
		# 2) D @ P == Drhs : Values of P that correspond to alpha's and fixed TFs
		# should not change.
		# 3) pdiff @ P >= 0.1 : There must be at least a difference of 0.1
		# between binding probabilities of a TF in conditions TF__active and
		# TF__inactive
		constraint_p = [
			0 <= PROMOTER_SCALING * P, PROMOTER_SCALING * P <= PROMOTER_SCALING,
			np.diag(D) @ (PROMOTER_SCALING * P) == PROMOTER_SCALING * Drhs,
			pdiff @ (PROMOTER_SCALING * P) >= PROMOTER_SCALING * PROMOTER_PDIFF_THRESHOLD]

		# Solve optimization problem
		prob_p = Problem(objective_p, constraint_p)
		prob_p.solve(solver='ECOS')

		if prob_p.status != "optimal":
			raise Exception("Solver could not find optimal value")

		# Get optimal value of P
		p = np.array(P.value).reshape(-1)
		# Adjust for solver tolerance over bounds to get proper probabilities
		p[p < 0] = 0
		p[p > 1] = 1

		# Update pPromoterBound with fit p
		fromArray(p, pPromoterBound, pPromoterBoundIdxs)

		# Break from loop if parameters have converged
		if np.abs(np.linalg.norm(np.dot(H, p) - k, PROMOTER_NORM_TYPE) - lastNorm) < PROMOTER_CONVERGENCE_THRESHOLD:
			break
		else:
			lastNorm = np.linalg.norm(np.dot(H, p) - k, PROMOTER_NORM_TYPE)

	# Update sim_data with fit bound probabilities and RNAP initiation
	# probabilities computed from these bound probabilities
	sim_data.pPromoterBound = pPromoterBound
	updateSynthProb(sim_data, cellSpecs, kInfo, np.dot(H, p))

	return r


def fitLigandConcentrations(sim_data, cellSpecs):
	"""
	Using the fit values of pPromoterBound, updates the set concentrations of
	ligand metabolites and the kd's of the ligand-TF binding reactions.

	Requires
	--------
	- Fitted pPromoterBound: probabilities that a TF will bind to its promoter,
	fit by function fitPromoterBoundProbability().

	Inputs
	------
	- cellSpecs {condition (str): dict} - information about each condition

	Modifies
	--------
	- Set concentrations of metabolites that are ligands in 1CS
	- kd's of equilibrium reactions in 1CS
	"""
	cellDensity = sim_data.constants.cellDensity
	pPromoterBound = sim_data.pPromoterBound

	for tf in sorted(sim_data.tfToActiveInactiveConds):
		# Skip TFs that are not 1CS or are linked to genotypic perturbations
		if sim_data.process.transcription_regulation.tfToTfType[tf] != "1CS":
			continue
		if (len(sim_data.tfToActiveInactiveConds[tf]["active genotype perturbations"]) > 0
				or len(sim_data.tfToActiveInactiveConds[tf]["inactive genotype perturbations"]) > 0):
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
		activeCellVolume = (cellSpecs[activeKey]["avgCellDryMassInit"] /
			cellDensity / sim_data.mass.cellDryMassFraction)
		activeCountsToMolar = 1/(sim_data.constants.nAvogadro*activeCellVolume)
		activeSignalConc = ((activeCountsToMolar * cellSpecs[activeKey]["bulkAverageContainer"].count(metabolite))
			.asNumber(units.mol/units.L))
		inactiveCellVolume = (cellSpecs[inactiveKey]["avgCellDryMassInit"] /
			cellDensity / sim_data.mass.cellDryMassFraction)
		inactiveCountsToMolar = 1/(sim_data.constants.nAvogadro * inactiveCellVolume)
		inactiveSignalConc = ((inactiveCountsToMolar * cellSpecs[inactiveKey]["bulkAverageContainer"].count(metabolite))
			.asNumber(units.mol/units.L))

		# Update kd with fitted values of P and the bulk average concentrations
		# of the metabolite, and use this fitted kd to recalculate the set
		# amounts of the metabolite in metabolism
		p_active = pPromoterBound[activeKey][tf]
		p_inactive = pPromoterBound[inactiveKey][tf]

		if negativeSignal:
			if p_inactive == 0:
				raise ValueError('Inf ligand concentration from p_inactive = 0.'
					' Check results from fitPromoterBoundProbability and Kd values.')
			if 1 - p_active < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = (activeSignalConc**metaboliteCoeff) * p_active/(1 - p_active)

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (
				(kdNew*(1 - p_inactive)/p_inactive)**(1./metaboliteCoeff)*(units.mol/units.L))

		else:
			if p_active == 1:
				raise ValueError('Inf ligand concentration from p_active = 1.'
					' Check results from fitPromoterBoundProbability and Kd values.')
			if p_inactive < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = (inactiveSignalConc**metaboliteCoeff) * (1 - p_inactive)/p_inactive

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentrationUpdates.moleculeSetAmounts[metabolite] = (
				(kdNew*p_active/(1 - p_active))**(1./metaboliteCoeff)*(units.mol/units.L))

		# Fit reverse rate in line with fitted kd
		sim_data.process.equilibrium.setRevRate(boundId + "[c]", kdNew*fwdRate)


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


def calculateRnapRecruitment(sim_data, r):
	"""
	Constructs the basal_prob vector and delta_prob matrix from values of r.
	The basal_prob vector holds the basal transcription probabilities of each
	transcription unit. The delta_prob matrix holds the differences in
	transcription probabilities when transcription factors bind to the
	promoters of each transcription unit. Both values are stored in sim_data.

	Requires
	--------
	- r: Fit parameters on how the recruitment of a TF affects the expression
	of a gene. High (positive) values of r indicate that the TF binding
	increases the probability that the gene is expressed.

	Modifies
	--------
	- Rescales values in basal_prob such that all values are positive
	- Adds basal_prob and delta_prob arrays to sim_data
	"""

	colNames = []

	# Get list of transcription units and TF IDs
	all_TUs = sim_data.process.transcription.rnaData["id"]
	all_tfs = sim_data.process.transcription_regulation.tf_ids

	# Initialize basal_prob vector and delta_prob sparse matrix
	basal_prob = np.zeros(len(all_TUs))
	deltaI, deltaJ, deltaV = [], [], []

	for rna_idx, rnaId in enumerate(all_TUs):
		rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

		tfs = sim_data.process.transcription_regulation.targetTf.get(rnaIdNoLoc, [])
		tfsWithData = []

		# Take only those TFs with active/inactive conditions data
		for tf in tfs:
			if tf not in sorted(sim_data.tfToActiveInactiveConds):
				continue

			tfsWithData.append(
				{"id": tf, "mass_g/mol": sim_data.getter.getMass([tf]).asNumber(units.g/units.mol)}
				)

		# Add one column for each TF that regulates the RNA
		for tf in tfsWithData:
			colName = rnaIdNoLoc + "__" + tf["id"]
			if colName not in colNames:
				colNames.append(colName)

			# Set element in delta to value in r that corresponds to the
			# transcription unit of the row, and the TF of the column
			deltaI.append(rna_idx)
			deltaJ.append(all_tfs.index(tf["id"]))
			deltaV.append(r[colNames.index(colName)])

		# Add alpha column for each RNA
		colName = rnaIdNoLoc + "__alpha"
		if colName not in colNames:
			colNames.append(colName)

		# Set element in basal_prob to the transcription unit's value for alpha
		basal_prob[rna_idx] = r[colNames.index(colName)]

	# Convert to arrays
	deltaI, deltaJ, deltaV = np.array(deltaI), np.array(deltaJ), np.array(deltaV)
	delta_shape = (len(all_TUs), len(all_tfs))

	# Adjust any negative basal probabilities to 0
	basal_prob[basal_prob < 0] = 0

	# Add basal_prob vector and delta_prob matrix to sim_data
	sim_data.process.transcription_regulation.basal_prob = basal_prob
	sim_data.process.transcription_regulation.delta_prob = {
		"deltaI": deltaI,
		"deltaJ": deltaJ,
		"deltaV": deltaV,
		"shape": delta_shape,
		}


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

	# isMRna = sim_data.process.transcription.rnaData["isMRna"]
	# isRna = np.zeros(len(isMRna))

	endoRnaseRnaIds = sim_data.moleculeGroups.endoRnase_RnaIDs
	isEndoRnase = np.array([x in endoRnaseRnaIds for x in sim_data.process.transcription.rnaData["id"]])

	rnaCounts = bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
	# endoCounts = bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)

	rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
	Kmcounts = (( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc).asNumber()
	sim_data.process.rna_decay.KmFirstOrderDecay = Kmcounts

	# Residuals can be written as follows: Res = f(Km) = 0, then Km = g(Km)
	# Compute derivative g(Km) in counts:
	KmQuadratic = 1 / np.power((1 / countsToMolar * Kmcounts).asNumber(), 2)
	denominator = np.power(np.sum(rnaCounts / (1 / countsToMolar * Kmcounts).asNumber()),2)
	numerator = (1 / countsToMolar * totalEndoRnaseCapacity).asNumber() * (denominator - (rnaCounts / (1 / countsToMolar * Kmcounts).asNumber()))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print("Max derivative (counts) = %f" % max(gDerivative))

	# Compute derivative g(Km) in concentrations:
	KmQuadratic = 1 / np.power(Kmcounts, 2)
	denominator = np.power(np.sum(rnaConc.asNumber() / Kmcounts),2)
	numerator = totalEndoRnaseCapacity.asNumber() * (denominator - (rnaConc.asNumber() / Kmcounts))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print("Max derivative (concentration) = %f" % max(gDerivative))


	# Sensitivity analysis: alpha (regularization term)
	Alphas = []
	if sim_data.constants.SensitivityAnalysisAlpha:
		Alphas = [0.0001, 0.001, 0.01, 0.1, 1, 10]

	for alpha in Alphas:

		if VERBOSE: print('Alpha = %f' % alpha)

		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.kmLossFunction(
				totalEndoRnaseCapacity.asNumber(units.mol / units.L / units.s),
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

		if VERBOSE: print('Kcat = %f' % kcat)

		totalEndoRNcap = units.sum(endoRNaseConc * kcat)
		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.kmLossFunction(
				totalEndoRNcap.asNumber(units.mol / units.L),
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
				totalEndoRnaseCapacity.asNumber(units.mol / units.L / units.s),
				(countsToMolar * rnaCounts).asNumber(units.mol / units.L),
				degradationRates.asNumber(1 / units.s),
				isEndoRnase,
				alpha
			)

	needToUpdate = False
	fixturesDir = filepath.makedirs(filepath.ROOT_PATH, "fixtures", "endo_km")
	# Numpy 'U' fields make these files incompatible with older code, so change
	# the filename. No need to make files compatible between Python 2 & 3; we'd
	# have to set the same protocol version and set Python 3-only args like
	# encoding='latin1'.
	km_filepath = os.path.join(fixturesDir, 'km{}.cPickle'.format(sys.version_info[0]))

	if os.path.exists(km_filepath):
		with open(km_filepath, "rb") as f:
			KmcountsCached = cPickle.load(f)

		# KmcountsCached fits a set of Km values to give the expected degradation rates.
		# It takes 1.5 - 3 minutes to recompute.
		# R_aux calculates the difference of the degradation rate based on these
		# Km values and the expected rate so this sum seems like a reliable test of
		# whether the cache fits current input data.
		if Kmcounts.shape != KmcountsCached.shape or np.sum(np.abs(R_aux(KmcountsCached))) > 1e-15:
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

		if VERBOSE: print("Running non-linear optimization")
		KmCooperativeModel = scipy.optimize.fsolve(LossFunction, Kmcounts, fprime = LossFunctionP)

		with open(km_filepath, "wb") as f:
			cPickle.dump(KmCooperativeModel, f, protocol=cPickle.HIGHEST_PROTOCOL)
	else:
		if VERBOSE:
			print("Not running non-linear optimization--using cached result {}".format(km_filepath))
		KmCooperativeModel = KmcountsCached

	if VERBOSE > 1:
		print("Loss function (Km inital) = %f" % np.sum(np.abs(LossFunction(Kmcounts))))
		print("Loss function (optimized Km) = %f" % np.sum(np.abs(LossFunction(KmCooperativeModel))))

		print("Negative km ratio = %f" % np.sum(np.abs(Rneg(KmCooperativeModel))))

		print("Residuals (Km initial) = %f" % np.sum(np.abs(R(Kmcounts))))
		print("Residuals optimized = %f" % np.sum(np.abs(R(KmCooperativeModel))))

		print("EndoR residuals (Km initial) = %f" % np.sum(np.abs(isEndoRnase * R(Kmcounts))))
		print("EndoR residuals optimized = %f" % np.sum(np.abs(isEndoRnase * R(KmCooperativeModel))))

		print("Residuals (scaled by Kdeg * RNAcounts) Km initial = %f" % np.sum(np.abs(R_aux(Kmcounts))))
		print("Residuals (scaled by Kdeg * RNAcounts) optimized = %f" % np.sum(np.abs(R_aux(KmCooperativeModel))))


	# Evaluate Jacobian around solutions (Kmcounts and KmCooperativeModel)
	JacobDiag = np.diag(Jacob(KmCooperativeModel))
	Jacob_auxDiag = np.diag(Jacob_aux(KmCooperativeModel))

	# Compute convergence of non-linear optimization: g'(Km)
	Gkm = np.abs(1. - JacobDiag)
	Gkm_aux = np.abs(1. - Jacob_auxDiag)
	sim_data.process.rna_decay.KmConvergence = Gkm_aux

	# Convergence is guaranteed if g'(Km) <= K < 1
	if VERBOSE:
		print("Convergence (Jacobian) = %.0f%% (<K> = %.5f)" % (len(Gkm[Gkm < 1.]) / float(len(Gkm)) * 100., np.mean(Gkm)))
		print("Convergence (Jacobian_aux) = %.0f%% (<K> = %.5f)" % (len(Gkm_aux[Gkm_aux < 1.]) / float(len(Gkm_aux)) * 100., np.mean(Gkm_aux[Gkm_aux < 1.])))

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
