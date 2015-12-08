#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.utils import units
from wholecell.utils.fitting import normalize

# Hacks
RNA_POLY_MRNA_DEG_RATE_PER_S = np.log(2) / 30. # half-life of 30 seconds
FRACTION_INCREASE_RIBOSOMAL_PROTEINS = 0  # reduce stochasticity from protein expression

# TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)

FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100
N_SEEDS = 20

DOUBLING_TIME = 60. * units.min
MEDIA_CONDITIONS = "M9 Glucose minus AAs"
TIME_STEP_SEC = None # If this is None the time step will be fit for the simulation in fitTimeStep

VERBOSE = False

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g

def fitSimData_1(raw_data):
	# Initialize simulation data with growth rate
	sim_data = SimulationDataEcoli()
	sim_data.initialize(doubling_time = DOUBLING_TIME, raw_data = raw_data, time_step_sec = TIME_STEP_SEC, media_conditions = MEDIA_CONDITIONS)

	# Increase RNA poly mRNA deg rates
	setRnaPolymeraseCodingRnaDegradationRates(sim_data)

	# Set C-period
	setCPeriod(sim_data)

	unfitExpression = sim_data.process.transcription.rnaData["expression"].copy()

	# Fit synthesis probabilities for RNA
	for iteration in xrange(MAX_FITTING_ITERATIONS):
		if VERBOSE: print 'Iteration: {}'.format(iteration)

		initialExpression = sim_data.process.transcription.rnaData["expression"].copy()

		setInitialRnaExpression(sim_data)

		bulkContainer = createBulkContainer(sim_data)

		rescaleMassForSoluableMetabolites(sim_data, bulkContainer)

		setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer)

		setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer)

		# Normalize expression and write out changes

		fitExpression(sim_data, bulkContainer)

		finalExpression = sim_data.process.transcription.rnaData["expression"]

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - finalExpression)))
		if VERBOSE: print 'degree of fit: {}'.format(degreeOfFit)

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	# Modify other properties

	fitRNAPolyTransitionRates(sim_data, bulkContainer)

	## Calculate and set maintenance values

	# ----- Growth associated maintenance -----

	fitMaintenanceCosts(sim_data, bulkContainer)

	fitTimeStep(sim_data, bulkContainer)


	calculateBulkDistributions(sim_data)

	return sim_data

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

def rescaleMassForSoluableMetabolites(sim_data, bulkMolCntr):
	avgCellSubMass = sim_data.mass.avgCellSubMass

	mass = (avgCellSubMass["proteinMass"] + avgCellSubMass["rnaMass"] + avgCellSubMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(sim_data.process.metabolism.metabolitePoolIDs) if sim_data.process.metabolism.metabolitePoolConcentrations.asNumber()[idx] > 0]
	poolConcentrations = (units.mol / units.L) * np.array([x for x in sim_data.process.metabolism.metabolitePoolConcentrations.asNumber() if x > 0])

	cellDensity = sim_data.constants.cellDensity
	mws = sim_data.getter.getMass(poolIds)
	concentrations = poolConcentrations.copy()

	diag = (cellDensity / (mws * concentrations) - 1).asNumber()
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass.asNumber(units.g) * np.ones(diag.size)

	massesToAdd = units.g * np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * sim_data.constants.nAvogadro

	V = (mass + units.sum(massesToAdd)) / cellDensity

	assert np.allclose(
		(countsToAdd / sim_data.constants.nAvogadro / V).asNumber(units.mol / units.L),
		(poolConcentrations).asNumber(units.mol / units.L)
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculePoolsDryMass = units.hstack((massesToAdd[:poolIds.index('WATER[c]')], massesToAdd[poolIds.index('WATER[c]') + 1:]))
	newAvgCellDryMassInit = units.sum(mass) + units.sum(smallMoleculePoolsDryMass)

	sim_data.mass.avgCellDryMassInit = newAvgCellDryMassInit
	sim_data.mass.avgCellDryMass = sim_data.mass.avgCellDryMassInit * sim_data.mass.avgCellToInitialCellConvFactor
	sim_data.mass.avgCellWaterMassInit = sim_data.mass.avgCellDryMass / 0.3 * 0.7

def setInitialRnaExpression(sim_data):
	# Set expression for all of the noncoding RNAs

	# Load from KB

	## IDs
	ids_rnas = sim_data.process.transcription.rnaData["id"]
	ids_rRNA23S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna23S"]]
	ids_rRNA16S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna16S"]]
	ids_rRNA5S = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isRRna5S"]]
	ids_tRNA = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isTRna"]]
	ids_mRNA = sim_data.process.transcription.rnaData["id"][sim_data.process.transcription.rnaData["isMRna"]]

	avgCellSubMass = sim_data.mass.avgCellSubMass

	## Mass fractions
	totalMass_rRNA23S = avgCellSubMass["rRna23SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA16S = avgCellSubMass["rRna16SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA5S = avgCellSubMass["rRna5SMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_tRNA = avgCellSubMass["tRnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	totalMass_mRNA = avgCellSubMass["mRnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

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
	distribution_mRNA = normalize(sim_data.process.transcription.rnaData["expression"][sim_data.process.transcription.rnaData['isMRna']])

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

	sim_data.process.transcription.rnaData["expression"] = normalize(rnaExpressionContainer.counts())
	# Note that now rnaData["synthProb"] does not match rnaData["expression"]

def totalCountIdDistributionProtein(sim_data):
	ids_protein = sim_data.process.translation.monomerData["id"]
	totalMass_protein = sim_data.mass.avgCellSubMass["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_protein = sim_data.process.translation.monomerData["mw"] / sim_data.constants.nAvogadro
	distribution_transcriptsByProtein = normalize(sim_data.process.transcription.rnaData["expression"][sim_data.relation.rnaIndexToMonomerMapping])
	translation_efficienciesByProtein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	degradationRates = sim_data.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(sim_data.doubling_time, degradationRates)

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

def totalCountIdDistributionRNA(sim_data):
	ids_rnas = sim_data.process.transcription.rnaData["id"]
	totalMass_RNA = sim_data.mass.avgCellSubMass["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	individualMasses_RNA = sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro

	distribution_RNA = normalize(sim_data.process.transcription.rnaData["expression"])

	totalCount_RNA = totalCountFromMassesAndRatios(
		totalMass_RNA,
		individualMasses_RNA,
		distribution_RNA
		)
	totalCount_RNA.normalize()
	totalCount_RNA.checkNoUnit()

	return totalCount_RNA, ids_rnas, distribution_RNA

def createBulkContainer(sim_data):

	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data)
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


def setRibosomeCountsConstrainedByPhysiology(sim_data, bulkContainer):
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
		sim_data.doubling_time,
		proteinDegradationRates
		)

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
	proteinLengths, sim_data.growthRateParameters.ribosomeElongationRate, netLossRate_protein, proteinCounts)
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
	if VERBOSE: print '30S limit: {}'.format(constraint_names[np.where(rib30lims.max() == rib30lims)[0]][0])
	if VERBOSE: print '30S actual count: {}'.format((ribosome30SCounts / ribosome30SStoich).min())
	if VERBOSE: print '30S count set to: {}'.format(rib30lims[np.where(rib30lims.max() == rib30lims)[0]][0])
	if VERBOSE: print '50S limit: {}'.format(constraint_names[np.where(rib50lims.max() == rib50lims)[0]][0])
	if VERBOSE: print '50S actual count: {}'.format((ribosome50SCounts / ribosome50SStoich).min())
	if VERBOSE: print '50S count set to: {}'.format(rib50lims[np.where(rib50lims.max() == rib50lims)[0]][0])

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


def setRNAPCountsConstrainedByPhysiology(sim_data, bulkContainer):
	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(sim_data.process.transcription.rnaData['countsACGU'], axis = 1)

	# Get constants to compute countsToMolar factor
	cellDensity = sim_data.constants.cellDensity
	cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / 0.3
	countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

	# Compute Km's
	rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rnaData['id'])
	degradationRates = sim_data.process.transcription.rnaData["degRate"]
	endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
	kcatEndoRNase = sim_data.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)
	Km = ( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc
	
	# Set Km's
	sim_data.process.transcription.rnaData["KmEndoRNase"] = Km

	rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
		sim_data.doubling_time,
		(1 / countsToMolar) * totalEndoRnaseCapacity,
		Km, 
		rnaConc,
		countsToMolar,
		)
	
	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths, sim_data.growthRateParameters.rnaPolymeraseElongationRate, rnaLossRate)

	nActiveRnapNeeded = units.convertNoUnitToNumber(nActiveRnapNeeded)
	nRnapsNeeded = nActiveRnapNeeded / sim_data.growthRateParameters.fractionActiveRnap

	rnapIds = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])['subunitIds']
	rnapStoich = sim_data.process.complexation.getMonomers(sim_data.moleculeGroups.rnapFull[0])['subunitStoich']

	minRnapSubunitCounts = (
		nRnapsNeeded * rnapStoich # Subunit stoichiometry
		)

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(rnapIds)

	## -- SET RNAP COUNTS TO MAXIMIM CONSTRAINTS -- #
	constraint_names = np.array(["Current level OK", "Insufficient to double RNA distribution"])
	rnapLims = np.array([(rnapCounts / rnapStoich).min(), (minRnapSubunitCounts / rnapStoich).min()])
	if VERBOSE: print 'rnap limit: {}'.format(constraint_names[np.where(rnapLims.max() == rnapLims)[0]][0])
	if VERBOSE: print 'rnap actual count: {}'.format((rnapCounts / rnapStoich).min())
	if VERBOSE: print 'rnap counts set to: {}'.format(rnapLims[np.where(rnapLims.max() == rnapLims)[0]][0])

	bulkContainer.countsIs(np.fmax(rnapCounts, minRnapSubunitCounts), rnapIds)



def fitExpression(sim_data, bulkContainer):

	view_RNA = bulkContainer.countsView(sim_data.process.transcription.rnaData["id"])
	counts_protein = bulkContainer.counts(sim_data.process.translation.monomerData["id"])

	translation_efficienciesByProtein = normalize(sim_data.process.translation.translationEfficienciesByMonomer)

	avgCellSubMass = sim_data.mass.avgCellSubMass
	totalMass_RNA = avgCellSubMass["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	doublingTime = sim_data.doubling_time
	degradationRates_protein = sim_data.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates_protein)

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

	sim_data.process.transcription.rnaData["expression"] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = totalCountFromMassesAndRatios(
		totalMass_RNA,
		sim_data.process.transcription.rnaData["mw"] / sim_data.constants.nAvogadro,
		sim_data.process.transcription.rnaData["expression"]
		)

	nRnas.normalize()
	nRnas.checkNoUnit()

	view_RNA.countsIs(nRnas * sim_data.process.transcription.rnaData["expression"])

	# Get constants to compute countsToMolar factor
	cellDensity = sim_data.constants.cellDensity
	cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / 0.3
	countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

	# Compute total endornase maximum capacity
	endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRnaseIds)
	kcatEndoRNase = sim_data.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

	rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
		sim_data.doubling_time,
		(1 / countsToMolar) * totalEndoRnaseCapacity,
		sim_data.process.transcription.rnaData["KmEndoRNase"],
		countsToMolar * view_RNA.counts(),
		countsToMolar,
		)

	synthProb = normalize(rnaLossRate.asNumber(1 / units.min))

	sim_data.process.transcription.rnaData["synthProb"][:] = synthProb


def fitRNAPolyTransitionRates(sim_data, bulkContainer):
	## Transcription activation rate

	synthProb = sim_data.process.transcription.rnaData["synthProb"]

	rnaLengths = sim_data.process.transcription.rnaData["length"]

	elngRate = sim_data.growthRateParameters.rnaPolymeraseElongationRate

	# In our simplified model of RNA polymerase state transition, RNAp can be
	# active (transcribing) or inactive (free-floating).  To solve for the
	# rate of activation, we need to calculate the average rate of termination,
	# which is a function of the average transcript length and the
	# transcription rate.

	averageTranscriptLength = units.dot(synthProb, rnaLengths)

	expectedTerminationRate = elngRate / averageTranscriptLength

	sim_data.transcriptionActivationRate = expectedTerminationRate * sim_data.growthRateParameters.fractionActiveRnap / (1 - sim_data.growthRateParameters.fractionActiveRnap)

	sim_data.fracActiveRnap = sim_data.growthRateParameters.fractionActiveRnap




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
	sim_data.constants.gtpPerTranslation += additionalGtpPerTranslation

	sim_data.constants.darkATP = darkATP

def fitTimeStep(sim_data, bulkContainer):
	'''
	Assumes that major limitor of growth will be translation associated
	resources, specifically AAs or GTP.

	Basic idea is that the rate of usage scales at the same rate as the size of the
	pool of resources.

	[Polymerized resource] * ln(2) * dt / doubling_time < [Pool of resource]

	'''
	aaCounts = sim_data.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomerData["id"])
	aasInProteins = units.sum(aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)), axis = 0)

	# USE IF AA LIMITING - When metabolism is implementing GAM
	# aaPools = units.aa * bulkContainer.counts(sim_data.moleculeGroups.aaIDs)
	# avgCellDryMassInit = sim_data.mass.avgCellDryMassInit
	# cellDensity = sim_data.constants.cellDensity
	# cellVolume = avgCellDryMassInit / cellDensity

	# aaPoolConcentration = aaPools / cellVolume
	# aaPolymerizedConcentration = aasInProteins / cellVolume

	# time_step = (aaPoolConcentration / aaPolymerizedConcentration) * sim_data.doubling_time / np.log(2)

	# USE IF GTP LIMITING - When GAM is incorperated into GTP/aa polymerized
	gtpPool = bulkContainer.counts(['GTP[c]'])
	gtpPolymerizedPool = (aasInProteins.asNumber(units.aa) * sim_data.constants.gtpPerTranslation).sum()
	timeStep = ((gtpPool / gtpPolymerizedPool) * sim_data.doubling_time.asNumber(units.s) / np.log(2))[0]
	timeStep = np.floor(timeStep * 100) / 100.0 # Round down to 2nd decimal

	if sim_data.timeStepSec != None:
		raise Exception("timeStepSec was set to a specific value!")
	else:
		sim_data.timeStepSec = timeStep * 0.7

def calculateBulkDistributions(sim_data):

	# Ids
	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data)
	ids_complex = sim_data.process.complexation.moleculeNames
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_metabolites = [x for idx, x in enumerate(sim_data.process.metabolism.metabolitePoolIDs) if sim_data.process.metabolism.metabolitePoolConcentrations.asNumber()[idx] > 0]
	conc_metabolites = (units.mol / units.L) * np.array([x for x in sim_data.process.metabolism.metabolitePoolConcentrations.asNumber() if x > 0])
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
	cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity

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

		rnaView.countsIs(randomState.multinomial(
			totalCount_RNA,
			distribution_RNA
			))

		proteinView.countsIs(randomState.multinomial(
			totalCount_protein,
			distribution_protein
			))

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

