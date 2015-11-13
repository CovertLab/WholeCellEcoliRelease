#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
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

def fitKb_1(kb):
	# Initialize simulation data with growth rate
	raw_data = KnowledgeBaseEcoli()
	kb.initialize(doubling_time = DOUBLING_TIME, raw_data = raw_data, time_step_sec = TIME_STEP_SEC, media_conditions = MEDIA_CONDITIONS)

	# Increase RNA poly mRNA deg rates
	setRnaPolymeraseCodingRnaDegradationRates(kb)

	# Set C-period
	setCPeriod(kb)

	unfitExpression = kb.process.transcription.rnaData["expression"].copy()

	# Fit synthesis probabilities for RNA
	for iteration in xrange(MAX_FITTING_ITERATIONS):
		if VERBOSE: print 'Iteration: {}'.format(iteration)

		initialExpression = kb.process.transcription.rnaData["expression"].copy()

		setInitialRnaExpression(kb)

		bulkContainer = createBulkContainer(kb)

		rescaleMassForSoluableMetabolites(kb, bulkContainer)

		setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer)

		setRNAPCountsConstrainedByPhysiology(kb, bulkContainer)

		# Normalize expression and write out changes

		fitExpression(kb, bulkContainer)

		finalExpression = kb.process.transcription.rnaData["expression"]

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - finalExpression)))
		if VERBOSE: print 'degree of fit: {}'.format(degreeOfFit)

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	# Modify other properties

	fitRNAPolyTransitionRates(kb, bulkContainer)

	## Calculate and set maintenance values

	# ----- Growth associated maintenance -----

	fitMaintenanceCosts(kb, bulkContainer)

	fitTimeStep(kb, bulkContainer)


	calculateBulkDistributions(kb)

# Sub-fitting functions

def setRnaPolymeraseCodingRnaDegradationRates(kb):
	# Increase RNA poly mRNA deg rates
	# TODO: set this based on transcription unit structure
	# i.e. same synthesis prob. but different deg rates

	rnaPolySubunits = kb.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"]
	subunitIndexes = np.array([np.where(kb.process.translation.monomerData["id"] == id_)[0].item() for id_ in rnaPolySubunits]) # there has to be a better way...
	mRNA_indexes = kb.relation.rnaIndexToMonomerMapping[subunitIndexes]
	kb.process.transcription.rnaData.struct_array["degRate"][mRNA_indexes] = RNA_POLY_MRNA_DEG_RATE_PER_S

def setCPeriod(kb):
	kb.growthRateParameters.c_period = kb.process.replication.genome_length * units.nt / kb.growthRateParameters.dnaPolymeraseElongationRate / 2

def rescaleMassForSoluableMetabolites(kb, bulkMolCntr):
	avgCellSubMass = kb.mass.avgCellSubMass

	mass = (avgCellSubMass["proteinMass"] + avgCellSubMass["rnaMass"] + avgCellSubMass["dnaMass"]) / kb.mass.avgCellToInitialCellConvFactor

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(kb.process.metabolism.metabolitePoolIDs) if kb.process.metabolism.metabolitePoolConcentrations.asNumber()[idx] > 0]
	poolConcentrations = (units.mol / units.L) * np.array([x for x in kb.process.metabolism.metabolitePoolConcentrations.asNumber() if x > 0])

	cellDensity = kb.constants.cellDensity
	mws = kb.getter.getMass(poolIds)
	concentrations = poolConcentrations.copy()

	diag = (cellDensity / (mws * concentrations) - 1).asNumber()
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass.asNumber(units.g) * np.ones(diag.size)

	massesToAdd = units.g * np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * kb.constants.nAvogadro

	V = (mass + units.sum(massesToAdd)) / cellDensity

	assert np.allclose(
		(countsToAdd / kb.constants.nAvogadro / V).asNumber(units.mol / units.L),
		(poolConcentrations).asNumber(units.mol / units.L)
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculePoolsDryMass = units.hstack((massesToAdd[:poolIds.index('WATER[c]')], massesToAdd[poolIds.index('WATER[c]') + 1:]))
	newAvgCellDryMassInit = units.sum(mass) + units.sum(smallMoleculePoolsDryMass)

	kb.mass.avgCellDryMassInit = newAvgCellDryMassInit
	kb.mass.avgCellDryMass = kb.mass.avgCellDryMassInit * kb.mass.avgCellToInitialCellConvFactor

def setInitialRnaExpression(kb):
	# Set expression for all of the noncoding RNAs

	# Load from KB

	## IDs
	ids_rnas = kb.process.transcription.rnaData["id"]
	ids_rRNA23S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]]
	ids_rRNA16S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]]
	ids_rRNA5S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]]
	ids_tRNA = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isTRna"]]
	ids_mRNA = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]]

	avgCellSubMass = kb.mass.avgCellSubMass

	## Mass fractions
	totalMass_rRNA23S = avgCellSubMass["rRna23SMass"] / kb.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA16S = avgCellSubMass["rRna16SMass"] / kb.mass.avgCellToInitialCellConvFactor
	totalMass_rRNA5S = avgCellSubMass["rRna5SMass"] / kb.mass.avgCellToInitialCellConvFactor
	totalMass_tRNA = avgCellSubMass["tRnaMass"] / kb.mass.avgCellToInitialCellConvFactor
	totalMass_mRNA = avgCellSubMass["mRnaMass"] / kb.mass.avgCellToInitialCellConvFactor

	## Molecular weights
	individualMasses_RNA = kb.getter.getMass(ids_rnas) / kb.constants.nAvogadro
	individualMasses_rRNA23S = kb.getter.getMass(ids_rRNA23S) / kb.constants.nAvogadro
	individualMasses_rRNA16S = kb.getter.getMass(ids_rRNA16S) / kb.constants.nAvogadro
	individualMasses_rRNA5S = kb.getter.getMass(ids_rRNA5S) / kb.constants.nAvogadro
	individualMasses_tRNA = kb.process.transcription.rnaData["mw"][kb.process.transcription.rnaData["isTRna"]] / kb.constants.nAvogadro
	individualMasses_mRNA = kb.process.transcription.rnaData["mw"][kb.process.transcription.rnaData["isMRna"]] / kb.constants.nAvogadro

	## Molecule expression distributions
	distribution_rRNA23S = np.array([1.] + [0.] * (ids_rRNA23S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA16S = np.array([1.] + [0.] * (ids_rRNA16S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA5S = np.array([1.] + [0.] * (ids_rRNA5S.size-1)) # currently only expressing first rRNA operon
	distribution_tRNA = normalize(kb.mass.getTrnaDistribution()['molar_ratio_to_16SrRNA'])
	distribution_mRNA = normalize(kb.process.transcription.rnaData["expression"][kb.process.transcription.rnaData['isMRna']])

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

	kb.process.transcription.rnaData["expression"] = normalize(rnaExpressionContainer.counts())
	# Note that now rnaData["synthProb"] does not match rnaData["expression"]

def totalCountIdDistributionProtein(kb):
	ids_protein = kb.process.translation.monomerData["id"]
	totalMass_protein = kb.mass.avgCellSubMass["proteinMass"] / kb.mass.avgCellToInitialCellConvFactor
	individualMasses_protein = kb.process.translation.monomerData["mw"] / kb.constants.nAvogadro
	distribution_transcriptsByProtein = normalize(kb.process.transcription.rnaData["expression"][kb.relation.rnaIndexToMonomerMapping])

	degradationRates = kb.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(kb.doubling_time, degradationRates)

	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcriptsByProtein,
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

def totalCountIdDistributionRNA(kb):
	ids_rnas = kb.process.transcription.rnaData["id"]
	totalMass_RNA = kb.mass.avgCellSubMass["rnaMass"] / kb.mass.avgCellToInitialCellConvFactor
	individualMasses_RNA = kb.process.transcription.rnaData["mw"] / kb.constants.nAvogadro

	distribution_RNA = normalize(kb.process.transcription.rnaData["expression"])

	totalCount_RNA = totalCountFromMassesAndRatios(
		totalMass_RNA,
		individualMasses_RNA,
		distribution_RNA
		)
	totalCount_RNA.normalize()
	totalCount_RNA.checkNoUnit()

	return totalCount_RNA, ids_rnas, distribution_RNA

def createBulkContainer(kb):

	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(kb)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(kb)
	ids_molecules = kb.state.bulkMolecules.bulkData["id"]

	## Construct bulk container

	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	## Assign RNA counts based on mass and expression distribution

	counts_RNA = totalCount_RNA * distribution_RNA

	bulkContainer.countsIs(counts_RNA, ids_rnas)

	## Assign protein counts based on mass and mRNA counts

	counts_protein = totalCount_protein * distribution_protein

	bulkContainer.countsIs(counts_protein, ids_protein)

	return bulkContainer


def setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer):
	'''
	setRibosomeCountsConstrainedByPhysiology

	Methodology: Set counts of ribosomal subunits based on three constraints.
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal subunit counts based on expression
	'''
	ribosome30SSubunits = kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds']
	ribosome50SSubunits = kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds']
	ribosome30SStoich = kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitStoich']
	ribosome50SStoich = kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitStoich']

	# -- CONSTRAINT 1: Expected protien distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(kb.process.translation.monomerData['aaCounts'], axis = 1)
	proteinDegradationRates =  kb.process.translation.monomerData["degRate"]
	proteinCounts =  bulkContainer.counts(kb.process.translation.monomerData["id"])

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(
		kb.doubling_time,
		proteinDegradationRates
		)

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
	proteinLengths, kb.growthRateParameters.ribosomeElongationRate, netLossRate_protein, proteinCounts)
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
	rRna23SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]])
	rRna16SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]])
	rRna5SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]])

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
	bulkContainer.countsIs(rRna23SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]])
	bulkContainer.countsIs(rRna16SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]])
	bulkContainer.countsIs(rRna5SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]])


def setRNAPCountsConstrainedByPhysiology(kb, bulkContainer):
	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(kb.process.transcription.rnaData['countsACGU'], axis = 1)

	# Get constants to compute countsToMolar factor
	cellDensity = kb.constants.cellDensity
	cellVolume = kb.mass.avgCellDryMassInit / cellDensity
	countsToMolar = 1 / (kb.constants.nAvogadro * cellVolume)

	# Compute Km's
	rnaConc = countsToMolar * bulkContainer.counts(kb.process.transcription.rnaData['id'])
	degradationRates = kb.process.transcription.rnaData["degRate"]
	endoRNaseConc = countsToMolar * bulkContainer.counts(kb.process.rna_decay.endoRnaseIds)
	kcatEndoRNase = kb.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)
	Km = ( 1 / degradationRates * totalEndoRnaseCapacity ) - rnaConc
	
	# Set Km's
	kb.process.transcription.rnaData["KmEndoRNase"][:] = Km

	rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
		kb.doubling_time,
		(1 / countsToMolar) * totalEndoRnaseCapacity,
		Km, 
		rnaConc,
		countsToMolar,
		)
	
	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths, kb.growthRateParameters.rnaPolymeraseElongationRate, rnaLossRate)

	nActiveRnapNeeded = units.convertNoUnitToNumber(nActiveRnapNeeded)
	nRnapsNeeded = nActiveRnapNeeded / kb.growthRateParameters.fractionActiveRnap

	rnapIds = kb.process.complexation.getMonomers(kb.moleculeGroups.rnapFull[0])['subunitIds']
	rnapStoich = kb.process.complexation.getMonomers(kb.moleculeGroups.rnapFull[0])['subunitStoich']

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


def fitExpression(kb, bulkContainer):

	view_RNA = bulkContainer.countsView(kb.process.transcription.rnaData["id"])
	counts_protein = bulkContainer.counts(kb.process.translation.monomerData["id"])

	avgCellSubMass = kb.mass.avgCellSubMass
	totalMass_RNA = avgCellSubMass["rnaMass"] / kb.mass.avgCellToInitialCellConvFactor

	doublingTime = kb.doubling_time
	degradationRates_protein = kb.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates_protein)

	### Modify kbFit to reflect our bulk container ###

	## RNA and monomer expression ##
	rnaExpressionContainer = BulkObjectsContainer(list(kb.process.transcription.rnaData["id"]), dtype = np.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(view_RNA.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		kb.process.translation.monomerData["rnaId"][kb.relation.monomerIndexToRnaMapping] == kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids" # TODO: move to KB tests

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * mRNADistributionFromProtein(
			normalize(counts_protein), netLossRate_protein
			)[kb.relation.monomerIndexToRnaMapping]
		)

	kb.process.transcription.rnaData["expression"] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = totalCountFromMassesAndRatios(
		totalMass_RNA,
		kb.process.transcription.rnaData["mw"] / kb.constants.nAvogadro,
		kb.process.transcription.rnaData["expression"]
		)

	nRnas.normalize()
	nRnas.checkNoUnit()

	view_RNA.countsIs(nRnas * kb.process.transcription.rnaData["expression"])

	# Get constants to compute countsToMolar factor
	cellDensity = kb.constants.cellDensity
	cellVolume = kb.mass.avgCellDryMassInit / cellDensity
	countsToMolar = 1 / (kb.constants.nAvogadro * cellVolume)

	# Compute total endornase maximum capacity
	endoRNaseConc = countsToMolar * bulkContainer.counts(kb.process.rna_decay.endoRnaseIds)
	kcatEndoRNase = kb.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

	rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
		kb.doubling_time,
		(1 / countsToMolar) * totalEndoRnaseCapacity,
		kb.process.transcription.rnaData["KmEndoRNase"], 
		countsToMolar * view_RNA.counts(),
		countsToMolar,
		)

	synthProb = normalize(rnaLossRate.asNumber(1 / units.min))

	kb.process.transcription.rnaData["synthProb"][:] = synthProb


def fitRNAPolyTransitionRates(kb, bulkContainer):
	## Transcription activation rate

	synthProb = kb.process.transcription.rnaData["synthProb"]

	rnaLengths = kb.process.transcription.rnaData["length"]

	elngRate = kb.growthRateParameters.rnaPolymeraseElongationRate

	# In our simplified model of RNA polymerase state transition, RNAp can be
	# active (transcribing) or inactive (free-floating).  To solve for the
	# rate of activation, we need to calculate the average rate of termination,
	# which is a function of the average transcript length and the
	# transcription rate.

	averageTranscriptLength = units.dot(synthProb, rnaLengths)

	expectedTerminationRate = elngRate / averageTranscriptLength

	kb.transcriptionActivationRate = expectedTerminationRate * kb.growthRateParameters.fractionActiveRnap / (1 - kb.growthRateParameters.fractionActiveRnap)

	kb.fracActiveRnap = kb.growthRateParameters.fractionActiveRnap


def fitMaintenanceCosts(kb, bulkContainer):
	aaCounts = kb.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(kb.process.translation.monomerData["id"])
	nAvogadro = kb.constants.nAvogadro
	avgCellDryMassInit = kb.mass.avgCellDryMassInit
	gtpPerTranslation = kb.constants.gtpPerTranslation

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
		kb.constants.growthAssociatedMaintenance -
		gtpUsedOverCellCycleMmolPerGDCW
		)

	additionalGtpPerTranslation = darkATP / aasUsedOverCellCycle
	additionalGtpPerTranslation.normalize()
	additionalGtpPerTranslation.checkNoUnit()
	additionalGtpPerTranslation = additionalGtpPerTranslation.asNumber()

	# Assign the growth associated "dark energy" to translation
	# TODO: Distribute it amongst growth-related processes
	kb.constants.gtpPerTranslation += additionalGtpPerTranslation

	kb.constants.darkATP = darkATP

def fitTimeStep(kb, bulkContainer):
	'''
	Assumes that major limitor of growth will be translation associated
	resources, specifically AAs or GTP.

	Basic idea is that the rate of usage scales at the same rate as the size of the
	pool of resources.

	[Polymerized resource] * ln(2) * dt / doubling_time < [Pool of resource]

	'''
	aaCounts = kb.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(kb.process.translation.monomerData["id"])
	aasInProteins = units.sum(aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)), axis = 0)

	# USE IF AA LIMITING - When metabolism is implementing GAM
	# aaPools = units.aa * bulkContainer.counts(kb.moleculeGroups.aaIDs)
	# avgCellDryMassInit = kb.mass.avgCellDryMassInit
	# cellDensity = kb.constants.cellDensity
	# cellVolume = avgCellDryMassInit / cellDensity

	# aaPoolConcentration = aaPools / cellVolume
	# aaPolymerizedConcentration = aasInProteins / cellVolume

	# time_step = (aaPoolConcentration / aaPolymerizedConcentration) * kb.doubling_time / np.log(2)

	# USE IF GTP LIMITING - When GAM is incorperated into GTP/aa polymerized
	gtpPool = bulkContainer.counts(['GTP[c]'])
	gtpPolymerizedPool = (aasInProteins.asNumber(units.aa) * kb.constants.gtpPerTranslation).sum()
	timeStep = ((gtpPool / gtpPolymerizedPool) * kb.doubling_time.asNumber(units.s) / np.log(2))[0]
	timeStep = np.floor(timeStep * 100) / 100.0 # Round down to 2nd decimal

	if kb.timeStepSec != None:
		raise Exception("timeStepSec was set to a specific value!")
	else:
		kb.timeStepSec = timeStep * 0.7

def calculateBulkDistributions(kb):

	# Ids
	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(kb)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(kb)
	ids_complex = kb.process.complexation.moleculeNames
	ids_equilibrium = kb.process.equilibrium.moleculeNames
	allMoleculesIDs = sorted(
		set(ids_rnas) | set(ids_protein) | set(ids_complex) | set(ids_equilibrium)
		)

	# Data for complexation

	complexationStoichMatrix = kb.process.complexation.stoichMatrix().astype(np.int64, order = "F")

	complexationPrebuiltMatrices = mccBuildMatrices(
		complexationStoichMatrix
		)

	# Data for equilibrium binding
	equilibriumDerivatives = kb.process.equilibrium.derivatives
	equilibriumDerivativesJacobian = kb.process.equilibrium.derivativesJacobian

	# Construct bulk container

	# We want to know something about the distribution of the copy numbers of
	# macromolecules in the cell.  While RNA and protein expression can be
	# approximated using well-described statistical distributions, we need
	# absolute copy numbers to form complexes.  To get a distribution, we must
	# instantiate many cells, form complexes, and finally compute the
	# statistics we will use in the fitting operations.

	bulkContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'])
	rnaView = bulkContainer.countsView(ids_rnas)
	proteinView = bulkContainer.countsView(ids_protein)
	complexationMoleculesView = bulkContainer.countsView(ids_complex)
	equilibriumMoleculesView = bulkContainer.countsView(ids_equilibrium)
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

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	bulkAverageContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'], np.float64)

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


def proteinDistributionFrommRNA(distribution_mRNA, netLossRate):
	"""
	dP_i / dt = k * M_i - P_i * Loss_i

	At steady state:
	P_i = k * M_i / Loss_i

	Fraction of mRNA for ith gene is defined as:
	f_i = M_i / M_total

	Substituting in:
	P_i = k * f_i * M_total / Loss_i

	Normalizing P_i by summing over all i cancels out k and M_total
	assuming constant translation rate.
	"""

	assert np.allclose(np.sum(distribution_mRNA), 1)
	distributionUnnormed = 1 / netLossRate * distribution_mRNA
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()
	return distributionNormed.asNumber()


def mRNADistributionFromProtein(distribution_protein, netLossRate):
	"""
	dP_i / dt = k * M_i - P_i * Loss_i

	At steady state:
	M_i = Loss_i * P_i / k

	Fraction of protein for ith gene is defined as:
	f_i = P_i / P_total

	Substituting in:
	M_i = Loss_i * f_i * P_total / k

	Normalizing M_i by summing over all i cancles out k and P_total
	assuming a constant translation rate.

	"""
	assert np.allclose(np.sum(distribution_protein), 1)
	distributionUnnormed = netLossRate * distribution_protein
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
	return ((np.log(2) / doublingTime) + totalEndoRnaseCountsCapacity * fracSaturated) * rnaCounts

