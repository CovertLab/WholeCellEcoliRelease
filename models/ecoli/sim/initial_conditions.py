
"""

TODO:
- document math
- raise/warn if physiological metabolite concentration targets appear to be smaller than what
  is needed at this time step size

"""

from __future__ import division

from itertools import izip
import scipy.sparse

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, calcProteinCounts, massesAndCountsToAddForHomeostaticTargets
from wholecell.utils.polymerize import buildSequences, computeMassIncrease
from wholecell.utils import units
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices


from wholecell.io.tablereader import TableReader

def calcInitialConditions(sim, sim_data):
	assert sim._inheritedStatePath == None
	randomState = sim.randomState

	massCoeff = 1.0
	if sim._massDistribution:
		massCoeff = randomState.normal(loc = 1.0, scale = 0.1)

	bulkMolCntr = sim.states['BulkMolecules'].container
	uniqueMolCntr = sim.states["UniqueMolecules"].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff)
	initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

def initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff):

	## Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set RNA counts from expression
	initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set DNA
	initializeDNA(bulkMolCntr, sim_data, randomState)

	## Set other biomass components
	initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set constitutive expression
	initializeConstitutiveExpression(bulkMolCntr, sim_data, randomState)

	## Form complexes
	initializeComplexation(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState):

	initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data)

	# Activate rna polys, with fraction based on environmental conditions
	initializeRNApolymerase(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Activate ribosomes, with fraction based on environmental conditions
	initializeRibosomes(bulkMolCntr, uniqueMolCntr, sim_data, randomState)


def initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff):

	monomersView = bulkMolCntr.countsView(sim_data.process.translation.monomerData["id"])
	monomerMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	# TODO: unify this logic with the fitter so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	monomerExpression = normalize(
		sim_data.process.transcription.rnaExpression[sim_data.condition][sim_data.relation.rnaIndexToMonomerMapping] *
		sim_data.process.translation.translationEfficienciesByMonomer /
		(np.log(2) / sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		sim_data.process.translation.monomerData["mw"].asNumber(units.g/units.mol),
		monomerExpression,
		sim_data.constants.nAvogadro.asNumber(1/units.mol)
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)

def initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff):

	rnaView = bulkMolCntr.countsView(sim_data.process.transcription.rnaData["id"])
	rnaMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	rnaExpression = normalize(sim_data.process.transcription.rnaExpression[sim_data.condition])

	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		sim_data.process.transcription.rnaData["mw"].asNumber(units.g/units.mol),
		rnaExpression,
		sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		)

	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

def initializeDNA(bulkMolCntr, sim_data, randomState):

	chromosomeView = bulkMolCntr.countsView(sim_data.moleculeGroups.fullChromosome)
	chromosomeView.countsIs([1])

# TODO: remove checks for zero concentrations (change to assertion)
# TODO: move any rescaling logic to KB/fitting
def initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff):
	avgCellFractionMass = sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])

	mass = massCoeff * (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
		sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel][0][1]
		)
	concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition]))
	moleculeIds = sorted(concDict)
	moleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in moleculeIds])

	massesToAdd, countsToAdd = massesAndCountsToAddForHomeostaticTargets(
		mass,
		moleculeIds,
		moleculeConcentrations,
		sim_data.getter.getMass(moleculeIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		moleculeIds
		)

def initializeConstitutiveExpression(bulkMolCntr, sim_data, randomState):
	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	alphaNames = [x for x in recruitmentColNames if x.endswith("__alpha")]
	alphaView = bulkMolCntr.countsView(alphaNames)
	alphaView.countsIs(1)

def initializeComplexation(bulkMolCntr, sim_data, randomState):
	moleculeNames = sim_data.process.complexation.moleculeNames
	moleculeView = bulkMolCntr.countsView(moleculeNames)

	# save rnase counts for after complexation, need to be monomers
	rnases = sim_data.process.rna_decay.endoRnaseIds
	rnaseCounts = bulkMolCntr.countsView(rnases).counts()
	bulkMolCntr.countsIs(0, rnases)

	stoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64, order = "F")
	prebuiltMatrices = mccBuildMatrices(stoichMatrix)

	# form complexes until no new complexes form (some complexes are complexes of complexes)
	while True:
		moleculeCounts = moleculeView.counts()
		updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			randomState.randint(1000),
			stoichMatrix,
			*prebuiltMatrices
			)

		bulkMolCntr.countsIs(
			updatedMoleculeCounts,
			moleculeNames,
			)

		if not np.any(moleculeCounts - updatedMoleculeCounts):
			break

	bulkMolCntr.countsIs(rnaseCounts, rnases)


def initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	initializeReplication

	Purpose: Create the appropriate number of replication forks given the cell growth rate.
	"""

	## Determine the number and location of replication forks at the start of the cell cycle
	# Find growth rate constants
	C = sim_data.growthRateParameters.c_period.asUnit(units.min)
	D = sim_data.growthRateParameters.d_period.asUnit(units.min)
	tau = sim_data.conditionToDoublingTime[sim_data.condition].asUnit(units.min)
	genome_length = sim_data.process.replication.genome_length
	replication_length = np.ceil(.5*genome_length) * units.nt

	# Generate arrays specifying appropriate replication conditions
	sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

	# Determine the number of OriC's currently present in the cell
	numOric = determineNumOriC(C, D, tau)
	oriC = uniqueMolCntr.objectsNew('originOfReplication', numOric)

	# Check that sequenceIdx, sequenceLength, replicationRound, and
	# chromosomeIndex are equal length, numOric should be half the
	# size (4 DNAP/fork, only 2 oriC/fork)
	assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 4*(numOric - 1))

	# Return if no replication is occuring at all
	if(len(sequenceIdx) == 0):
		return

	## Update polymerases mass to account for already completed DNA
	# Determine the sequences of already-replicated DNA
	sequences = sim_data.process.replication.replication_sequences
	sequenceElongations = np.array(sequenceLength, dtype=np.int64)
	massIncreaseDna = computeMassIncrease(
			np.tile(sequences,(len(sequenceIdx) // 4,1)),
			sequenceElongations,
			sim_data.process.replication.replicationMonomerWeights.asNumber(units.fg)
			)

	# Update the attributes of replicating DNA polymerases
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', len(sequenceIdx))
	dnaPoly.attrIs(
		sequenceIdx = np.array(sequenceIdx),
		sequenceLength = np.array(sequenceLength),
		replicationRound = np.array(replicationRound),
		chromosomeIndex = np.array(chromosomeIndex),
		massDiff_DNA = massIncreaseDna,
		)


def initializeRNApolymerase(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	"""
	Purpose: Activates RNA polymerases as unique molecules, and distributes them along length of genes,
	decreases counts of unactivated RNA polymerases (APORNAP-CPLX[c]).

	Normalizes RNA poly placement per length of completed RNA, with synthesis probability based on each environmental condition
	"""

	# Load parameters
	nAvogadro = sim_data.constants.nAvogadro
	rnaLengths = sim_data.process.transcription.rnaData['length'].asNumber()
	currentNutrients = sim_data.conditions[sim_data.condition]['nutrients']
	fracActiveRnap = sim_data.process.transcription.rnapFractionActiveDict[currentNutrients]
	inactiveRnaPolyCounts = bulkMolCntr.countsView(['APORNAP-CPLX[c]']).counts()[0]
	rnaSequences = sim_data.process.transcription.transcriptionSequences
	ntWeights = sim_data.process.transcription.transcriptionMonomerWeights
	endWeight = sim_data.process.transcription.transcriptionEndWeight

	# Number of rnaPoly to activate
	rnaPolyToActivate = np.int64(fracActiveRnap * inactiveRnaPolyCounts)

	# Parameters for rnaSynthProb
	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	recruitmentView = bulkMolCntr.counts(recruitmentColNames)
	recruitmentData = sim_data.process.transcription_regulation.recruitmentData
	recruitmentMatrix = scipy.sparse.csr_matrix(
			(recruitmentData['hV'], (recruitmentData['hI'], recruitmentData['hJ'])),
			shape = recruitmentData['shape']
		)

	# Synthesis probabilities for different categories of genes
	rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
	rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
	rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

	# Determine changes from genetic perturbations
	genetic_perturbations = {}
	perturbations = getattr(sim_data, 'genetic_perturbations', {})
	if len(perturbations) > 0:
		rnaIdxs, synthProbs = zip(*[(int(np.where(sim_data.process.transcription.rnaData['id'] == rnaId)[0]), synthProb) for rnaId, synthProb in sim_data.genetic_perturbations.iteritems()])
		fixedSynthProbs = [synthProb for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
		fixedRnaIdxs = [rnaIdx for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
		genetic_perturbations = {'fixedRnaIdxs': fixedRnaIdxs, 'fixedSynthProbs': fixedSynthProbs}

	# If initiationShuffleIdxs does not exist, set value to None
	shuffleIdxs = getattr(sim_data.process.transcription, 'initiationShuffleIdxs', None)

	# ID Groups
	isRRna = sim_data.process.transcription.rnaData['isRRna']
	isMRna = sim_data.process.transcription.rnaData['isMRna']
	isTRna = sim_data.process.transcription.rnaData['isTRna']
	isRProtein = sim_data.process.transcription.rnaData['isRProtein']
	isRnap = sim_data.process.transcription.rnaData['isRnap']
	isRegulated = np.array([1 if x[:-3] in sim_data.process.transcription_regulation.targetTf or x in perturbations else 0 for x in sim_data.process.transcription.rnaData["id"]], dtype = np.bool)
	setIdxs = isRRna | isTRna | isRProtein | isRnap | isRegulated

	# Calculate synthesis probabilities based on transcription regulation
	rnaSynthProb = recruitmentMatrix.dot(recruitmentView)
	if len(genetic_perturbations) > 0:
		rnaSynthProb[genetic_perturbations['fixedRnaIdxs']] = genetic_perturbations['fixedSynthProbs']
	regProbs = rnaSynthProb[isRegulated]

	# Adjust probabilities to not be negative
	rnaSynthProb[rnaSynthProb < 0] = 0.0
	rnaSynthProb /= rnaSynthProb.sum()
	if np.any(rnaSynthProb < 0):
		raise Exception("Have negative RNA synthesis probabilities")

	# Adjust synthesis probabilities depending on environment
	synthProbFractions = rnaSynthProbFractions[currentNutrients]
	rnaSynthProb[isMRna] *= synthProbFractions['mRna'] / rnaSynthProb[isMRna].sum()
	rnaSynthProb[isTRna] *= synthProbFractions['tRna'] / rnaSynthProb[isTRna].sum()
	rnaSynthProb[isRRna] *= synthProbFractions['rRna'] / rnaSynthProb[isRRna].sum()
	rnaSynthProb[isRegulated] = regProbs
	rnaSynthProb[isRProtein] = rnaSynthProbRProtein[currentNutrients]
	rnaSynthProb[isRnap] = rnaSynthProbRnaPolymerase[currentNutrients]
	rnaSynthProb[rnaSynthProb < 0] = 0 # to avoid precision issue
	scaleTheRestBy = (1. - rnaSynthProb[setIdxs].sum()) / rnaSynthProb[~setIdxs].sum()
	rnaSynthProb[~setIdxs] *= scaleTheRestBy

	# Shuffle initiation rates if we're running the variant that calls this
	if shuffleIdxs is not None:
		rnaSynthProb = rnaSynthProb[shuffleIdxs]

	# normalize to length of rna
	synthProbLengthAdjusted = rnaSynthProb * rnaLengths
	synthProbNormalized = synthProbLengthAdjusted / synthProbLengthAdjusted.sum()

	# Sample a multinomial distribution of synthesis probabilities to determine what RNA are initialized
	nNewRnas = randomState.multinomial(rnaPolyToActivate, synthProbNormalized)

	# RNA Indices
	rnaIndices = np.empty(rnaPolyToActivate, np.int64)
	startIndex = 0
	nonzeroCount = (nNewRnas > 0)
	for rnaIndex, counts in izip(np.arange(nNewRnas.size)[nonzeroCount], nNewRnas[nonzeroCount]):
		rnaIndices[startIndex:startIndex+counts] = rnaIndex
		startIndex += counts

	# TODO (Eran) -- make sure there aren't any rnapolys at same location on same gene
	updatedLengths = np.array(randomState.rand(rnaPolyToActivate) * rnaLengths[rnaIndices], dtype=np.int)

	# update mass
	sequences = rnaSequences[rnaIndices]
	massIncreaseRna = computeMassIncrease(sequences, updatedLengths, ntWeights)
	massIncreaseRna[updatedLengths != 0] += endWeight  # add endWeight to all new Rna

	#update molecules. Attributes include which rnas are being transcribed, and the position (length)
	activeRnaPolys = uniqueMolCntr.objectsNew('activeRnaPoly', rnaPolyToActivate)
	activeRnaPolys.attrIs(
		rnaIndex = rnaIndices,
		transcriptLength = updatedLengths,
		massDiff_mRNA = massIncreaseRna,
		)
	bulkMolCntr.countsIs(inactiveRnaPolyCounts - rnaPolyToActivate, ['APORNAP-CPLX[c]'])


def initializeRibosomes(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	"""
	Purpose: Activates ribosomes as unique molecules, and distributes them along length of RNA,
	decreases counts of unactivated ribosomal subunits (ribosome30S and ribosome50S).

	Normalizes ribosomes placement per length of protein
	"""

	# Load parameters
	nAvogadro = sim_data.constants.nAvogadro
	currentNutrients = sim_data.conditions[sim_data.condition]['nutrients']
	fracActiveRibosome = sim_data.process.translation.ribosomeFractionActiveDict[currentNutrients]
	mrnaIds = sim_data.process.translation.monomerData['rnaId']
	proteinLengths = sim_data.process.translation.monomerData['length'].asNumber()
	proteinSequences = sim_data.process.translation.translationSequences
	translationEfficiencies = normalize(sim_data.process.translation.translationEfficienciesByMonomer)
	mRnas = bulkMolCntr.countsView(mrnaIds)
	aaWeightsIncorporated = sim_data.process.translation.translationMonomerWeights
	endWeight = sim_data.process.translation.translationEndWeight

	#find number of ribosomes to activate
	ribosome30S = bulkMolCntr.countsView(sim_data.moleculeGroups.s30_fullComplex).counts()[0]
	ribosome50S = bulkMolCntr.countsView(sim_data.moleculeGroups.s50_fullComplex).counts()[0]
	inactiveRibosomeCount = np.minimum(ribosome30S, ribosome50S)
	ribosomeToActivate = np.int64(fracActiveRibosome * inactiveRibosomeCount)

	# protein synthesis probabilities
	proteinInitProb = normalize(mRnas.counts() * translationEfficiencies)

	# normalize to protein length
	probLengthAdjusted = proteinInitProb * proteinLengths
	probNormalized = probLengthAdjusted / probLengthAdjusted.sum()

	# Sample a multinomial distribution of synthesis probabilities to determine what RNA are initialized
	nNewProteins = randomState.multinomial(ribosomeToActivate, probNormalized)

	# protein Indices
	proteinIndices = np.empty(ribosomeToActivate, np.int64)
	startIndex = 0
	nonzeroCount = (nNewProteins > 0)
	for proteinIndex, counts in izip(np.arange(nNewProteins.size)[nonzeroCount], nNewProteins[nonzeroCount]):
		proteinIndices[startIndex:startIndex+counts] = proteinIndex
		startIndex += counts

	# TODO (Eran) -- make sure there aren't any peptides at same location on same rna
	updatedLengths = np.array(randomState.rand(ribosomeToActivate) * proteinLengths[proteinIndices], dtype=np.int)

	# update mass
	sequences = proteinSequences[proteinIndices]
	massIncreaseProtein = computeMassIncrease(sequences, updatedLengths, aaWeightsIncorporated)
	massIncreaseProtein[updatedLengths != 0] += endWeight  # add endWeight to all new Rna

	# Create active 70S ribosomes and assign their protein Indices calculated above
	activeRibosomes = uniqueMolCntr.objectsNew('activeRibosome', ribosomeToActivate)
	activeRibosomes.attrIs(
		proteinIndex = proteinIndices,
		peptideLength = updatedLengths,
		massDiff_protein = massIncreaseProtein,
		)

	# decrease free 30S and 50S ribosomal subunit counts
	bulkMolCntr.countsIs(ribosome30S - ribosomeToActivate, sim_data.moleculeGroups.s30_fullComplex)
	bulkMolCntr.countsIs(ribosome50S - ribosomeToActivate, sim_data.moleculeGroups.s50_fullComplex)


def setDaughterInitialConditions(sim, sim_data):
	assert sim._inheritedStatePath != None

	import cPickle
	isDead = cPickle.load(open(os.path.join(sim._inheritedStatePath, "IsDead.cPickle"), "rb"))
	sim._isDead = isDead

	elngRate = cPickle.load(open(os.path.join(sim._inheritedStatePath, "ElngRate.cPickle"), "rb"))
	elng_rate_factor = cPickle.load(open(os.path.join(sim._inheritedStatePath, "elng_rate_factor.cPickle"), "rb"))
	if sim._growthRateNoise:
		sim.processes["PolypeptideElongation"].setElngRate = elngRate
		sim.processes["PolypeptideElongation"].elngRateFactor = elng_rate_factor

	bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkMolecules"))
	sim.states["BulkMolecules"].tableLoad(bulk_table_reader, 0)

	unique_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "UniqueMolecules"))
	sim.states["UniqueMolecules"].tableLoad(unique_table_reader, 0)

	time_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "Time"))
	initialTime = TableReader(os.path.join(sim._inheritedStatePath, "Time")).readAttribute("initialTime")
	sim._initialTime = initialTime


def determineChromosomeState(C, D, tau, replication_length):
	"""
	determineChromosomeState

	Purpose: calculates the number and position of replicating DNA polymerases
			 at the beginning of the cell cycle.

	Inputs: C  - the C period of the cell, the length of time between
			replication initiation and replication completion.
			D  - the D period of the cell, the length of time between
			completing replication of the chromosome and division of the cell.
			tau - the doubling time of the cell
			replication_length - the amount of DNA to be replicated per fork,
			usually half of the genome, in base-pairs

	Outputs: a tuple of vectors for input into the dnaPoly.attrIs() function
			of the initializeReplication() function in initial_conditions.py.
			These  vectors are, in order:
			sequenceIdx - an index for each of the four types/directions of DNA
						replication - leading and lagging strand of the forward
						and the reverse fork = 4 total. This vector is always
						simply [0,1,2,3] repeated once for each replication
						event. Ie for three active replication events (6 forks,
						12 polymerases) sequenceIdx = [0,1,2,3,0,1,2,3,0,1,2,3]
			sequenceLength - the position in the genome that each polymerase
						referenced in sequenceIdx has reached, in base-pairs.
						This is handled such that even though in reality some
						polymerases are replicating in different directions all
						inputs here are as though each starts at 0 and goes up
						the the number of base-pairs to be replicated.
			replicationRound - a unique integer stating in which replication
						generation the polymerase referenced by sequenceIdx.
						Each time all origins of replication in the cell fire,
						a new replication generation has started. This array
						is integer-valued, and counts from 0 (the oldest
						generation) up to n (the most recent initiation) event.
			chromosomeIndex - indicator variable for which daughter cell
						should inherit which polymerase at division. This
						array is only relevant to draw distinctions within a
						generaation of replicationRound, now between rounds.
						Within each generation in replicationRound (run of
						numbers with the same value), half should have
						chromosomeIndex = 0, and half chromosomeIndex
						= 1. This should be contiguous halves, ie
						[0,0,0,0,1,1,1,1], NOT interspersed like
						[0,1,0,1,0,1,0,1]. The half-and-half rule is excepted
						for any replication generation with only 4 polymerases
						(the oldest replication generation should be the only
						one	with fewer than 8 polymerases). In this case
						chromosomeIndex doesn't matter/is effectively NaN,
						but is set to all 0's to prevent conceptually dividing
						a single chromosome between two daughter cells.

	Notes: if NO polymerases are active at the start of the cell cycle,
			equivalent to the C + D periods being shorter than the doubling
			rate tau, then this function returns empty lists. dnaPoly.attrIs()
			should not be run in this case, as no DNA replication will be
			underway.
	"""

	## Error check inputs

	# Check that all inputs have units
	assert (units.hasUnit(C)), 'C must have units'
	assert (units.hasUnit(D)), 'D must have units'
	assert (units.hasUnit(tau)), 'tau must have units'
	assert (units.hasUnit(replication_length)) and (units.getUnit(replication_length).strUnit() == units.nt.strUnit()), 'replication_length must have units of units.nt.'

	# All inputs must be positive numbers
	assert (C.asNumber(units.min) >= 0), "C value can't be negative."
	assert (D.asNumber(units.min) >= 0), "D value can't be negative."
	assert (tau.asNumber(units.min) >= 0), "tau value can't be negative."
	assert (replication_length.asNumber(units.nt) >= 0), "replication_length value can't be negative."

	# Require that D is shorter than tau - time between completing replication
	# and dividing must be shorter than the time between two divisions.
	assert (D.asNumber(units.min) < tau.asNumber(units.min)), 'The D period must be shorter than the doubling time tau.'


	# Number active replication generations (can be many initiations per gen.)
	limit = np.floor((C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min))

	# Initialize arrays to be returned
	sequenceIdx = []
	sequenceLength = []
	replicationRound = []
	chromosomeIndex = []

	# Loop through the generations of replication which are active (if limit = 0
	# skips loop entirely --> no active replication generations)
	n = 1;
	while n <= limit:
		# Determine at what base each strand of a given replication event should start
		# Replication forks should be at base (1 - (n*tau - D)/(C))(basepairs to be replicated)

		ratio = (1 - ((n*tau - D)/(C)))
		ratio = units.convertNoUnitToNumber(ratio)
		fork_location = np.floor(ratio*(replication_length.asNumber(units.nt)))

		# Add 2^(n-1) replication events (two forks, four strands per inintiaion event)
		num_events = 2 ** (n-1)
		# sequenceIdx refers to the type of elongation - ie forward and
		# reverse, lagging and leading strands.
		sequenceIdx += [0,1,2,3]*num_events
		# sequenceLength refers to how far along in the replication process
		# this event already is - at what basepair currently. All four are
		# assumed to be equally far along.
		sequenceLength += [fork_location]*4*num_events
		# replicationRound is defined the same way as n - each time a
		# replication round starts, all origins in the cell fire, and are then
		# of the same generation and round number
		replicationRound += [(n-1)]*4*num_events
		# WITHIN each round, chromosomeIndex uniquely identifies individual
		# origin initaion points. Loop through each intiation event in this
		# generation (2 forks, 4 polymerases each), assign it an increaing,
		# unique number, starting at zero.
		chromosomeIndex += [0]*2*num_events + [0]*2*num_events

		n += 1

	# The first replication generation should not be divided, so set all values
	# to 0 (effectively NaN, the first four values are not used)
	if len(chromosomeIndex):
		chromosomeIndex[:4] = [0,0,0,0]

	return (sequenceIdx, sequenceLength, replicationRound, chromosomeIndex)



def determineNumOriC(C, D, tau):
	"""
	determineNumOriC

	Purpose: calculates the number of OriC's in a cell upon initiation, determined by the replication state of the chromosome.

	Inputs: C  - the C period of the cell, the length of time between
			replication initiation and replication completion.
			D  - the D period of the cell, the length of time between
			completing replication of the chromosome and division of the cell.
			tau - the doubling time of the cell
			replication_length - the amount of DNA to be replicated per fork,
			usually half of the genome, in base-pairs

	Outputs: the number of OriC's in the cell at initiation.
	"""

	# Number active replication generations (can be many initiations per gen.)
	total_active_initiations = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())

	return int(2 ** (total_active_initiations))
