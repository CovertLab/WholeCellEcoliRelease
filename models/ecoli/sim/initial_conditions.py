"""
TODO: document math
TODO: raise/warn if physiological metabolite concentration targets appear to be smaller than what
  is needed at this time step size
"""

from __future__ import absolute_import, division, print_function

from itertools import izip
import scipy.sparse

import numpy as np

from wholecell.utils.fitting import normalize, countsFromMassAndExpression, masses_and_counts_for_homeostatic_target
from wholecell.utils.polymerize import computeMassIncrease
from wholecell.utils import units
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.sim.divide_cell import load_inherited_state


def calcInitialConditions(sim, sim_data):
	assert sim._inheritedStatePath is None
	randomState = sim.randomState

	massCoeff = 1.0
	if sim._massDistribution:
		massCoeff = randomState.normal(loc = 1.0, scale = 0.1)

	bulkMolCntr = sim.internal_states['BulkMolecules'].container
	uniqueMolCntr = sim.internal_states["UniqueMolecules"].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff)
	initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Must be called after unique and bulk molecules are initialized to get
	# concentrations for ribosomes, tRNA, synthetases etc from cell volume
	if sim._trna_charging:
		initialize_trna_charging(sim_data, sim.internal_states, sim.processes['PolypeptideElongation'].calculate_trna_charging)

def initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff):

	# Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff)

	# Set RNA counts from expression
	initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff)

	# Set DNA
	initializeDNA(bulkMolCntr, sim_data, randomState)

	# Set other biomass components
	initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff)

	# Form complexes
	initializeComplexation(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState):

	initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data)

	# Activate rna polys, with fraction based on environmental conditions
	initializeRNApolymerase(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Activate ribosomes, with fraction based on environmental conditions
	initializeRibosomes(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

def initialize_trna_charging(sim_data, states, calc_charging):
	'''
	Initializes charged tRNA from uncharged tRNA and amino acids

	Inputs:
		sim_data (SimulationDataEcoli object)
		states (dict with internal_state objects as values) - internal states of sim
		calc_charging (function) - function to calculate charging of tRNA

	Notes:
		Does not adjust for mass of amino acids on charged tRNA (~0.01% of cell mass)
	'''

	# Calculate cell volume for concentrations
	mass = 0
	for state in states.values():
		state.calculatePreEvolveStateMass()
		mass += np.sum(state._masses)
	cell_volume = units.fg * mass / sim_data.constants.cellDensity
	counts_to_molar = 1 / (sim_data.constants.nAvogadro * cell_volume)

	# Get molecule views and concentrations
	transcription = sim_data.process.transcription
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna
	bulk_molecules = states['BulkMolecules'].container
	synthetases = bulk_molecules.countsView(transcription.synthetase_names)
	uncharged_trna = bulk_molecules.countsView(transcription.rnaData['id'][transcription.rnaData['isTRna']])
	charged_trna = bulk_molecules.countsView(transcription.charged_trna_names)
	aas = bulk_molecules.countsView(sim_data.moleculeGroups.aaIDs)
	ribosome_counts = states['UniqueMolecules'].container.counts(['activeRibosome'])

	synthetase_conc = counts_to_molar * np.dot(aa_from_synthetase, synthetases.counts())
	uncharged_trna_conc = counts_to_molar * np.dot(aa_from_trna, uncharged_trna.counts())
	charged_trna_conc = counts_to_molar * np.dot(aa_from_trna, charged_trna.counts())
	aa_conc = counts_to_molar * aas.counts()
	ribosome_conc = counts_to_molar * ribosome_counts

	# Estimate fraction of amino acids from sequences, excluding first index for padding of -1
	_, aas_in_sequences = np.unique(sim_data.process.translation.translationSequences, return_counts=True)
	f = aas_in_sequences[1:] / np.sum(aas_in_sequences[1:])

	# Estimate initial charging state
	fraction_charged, _ = calc_charging(synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f)

	# Update counts of tRNA to match charging
	total_trna_counts = uncharged_trna.counts() + charged_trna.counts()
	charged_trna_counts = np.round(total_trna_counts * np.dot(fraction_charged, aa_from_trna))
	uncharged_trna_counts = total_trna_counts - charged_trna_counts
	charged_trna.countsIs(charged_trna_counts)
	uncharged_trna.countsIs(uncharged_trna_counts)

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

	chromosomeView = bulkMolCntr.countsView([sim_data.moleculeIds.fullChromosome])
	chromosomeView.countsIs([1])

# TODO: remove checks for zero concentrations (change to assertion)
# TODO: move any rescaling logic to KB/fitting
def initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff):
	avgCellFractionMass = sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])

	mass = massCoeff * (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
		sim_data.external_state.environment.nutrients_time_series[sim_data.external_state.environment.nutrients_time_series_label][0][1]
		)
	concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition]))
	moleculeIds = sorted(concDict)
	moleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in moleculeIds])

	massesToAdd, countsToAdd = masses_and_counts_for_homeostatic_target(
		mass,
		moleculeConcentrations,
		sim_data.getter.getMass(moleculeIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		moleculeIds
		)

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
		updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
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
	Initializes replication by creating an appropriate number of replication
	forks given the cell growth rate. This also initializes the gene dosage
	bulk counts using the initial locations of the forks.
	"""

	# Determine the number and location of replication forks at the start of
	# the cell cycle
	# Get growth rate constants
	C = sim_data.growthRateParameters.c_period.asUnit(units.min)
	D = sim_data.growthRateParameters.d_period.asUnit(units.min)
	tau = sim_data.conditionToDoublingTime[sim_data.condition].asUnit(units.min)

	# Calculate replication length
	genome_length = sim_data.process.replication.genome_length
	replication_length = np.ceil(0.5*genome_length) * units.nt

	# Generate arrays specifying appropriate initial replication conditions
	n_oric, chromosomeIndexOriC = determineOriCState(C, D, tau)
	sequenceIdx, sequenceLength, replicationRound, chromosomeIndexPolymerase = determineChromosomeState(
		C, D, tau, replication_length)
	n_dnap = sequenceIdx.size

	oriC = uniqueMolCntr.objectsNew('originOfReplication', n_oric)
	oriC.attrIs(chromosomeIndex = chromosomeIndexOriC)

	if n_dnap != 0:
		# Update mass to account for DNA strands that have already been elongated
		# Determine the sequences of already-replicated DNA
		sequences = sim_data.process.replication.replication_sequences
		sequenceElongations = sequenceLength.astype(np.int64)
		massIncreaseDna = computeMassIncrease(
				np.tile(sequences, (n_dnap//4, 1)),
				sequenceElongations,
				sim_data.process.replication.replicationMonomerWeights.asNumber(units.fg)
				)

		# Add replicating DNA polymerases as unique molecules and set attributes
		dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', n_dnap)
		dnaPoly.attrIs(
			sequenceIdx = sequenceIdx,
			sequenceLength = sequenceLength,
			replicationRound = replicationRound,
			chromosomeIndex = chromosomeIndexPolymerase,
			massDiff_DNA = massIncreaseDna,
			)

	# Initialize gene dosage
	geneCopyNumberColNames = sim_data.process.transcription_regulation.geneCopyNumberColNames
	geneCopyNumberView = bulkMolCntr.countsView(geneCopyNumberColNames)
	replicationCoordinate = sim_data.process.transcription.rnaData["replicationCoordinate"]

	# Set all copy numbers to one initially
	initialGeneCopyNumber = np.ones(len(geneCopyNumberColNames))

	# Get coordinates of forks in both directions
	forward_fork_coordinates = sequenceLength[sequenceIdx == 0]
	reverse_fork_coordinates = np.negative(sequenceLength[sequenceIdx == 1])

	assert len(forward_fork_coordinates) == len(reverse_fork_coordinates)

	# Increment copy number by one for any gene that lies between two forks
	for (forward, reverse) in izip(forward_fork_coordinates,
			reverse_fork_coordinates):
		initialGeneCopyNumber[
			np.logical_and(replicationCoordinate < forward,
				replicationCoordinate > reverse)
			] += 1

	geneCopyNumberView.countsIs(initialGeneCopyNumber)


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
	ribosome30S = bulkMolCntr.countsView([sim_data.moleculeIds.s30_fullComplex]).counts()[0]
	ribosome50S = bulkMolCntr.countsView([sim_data.moleculeIds.s50_fullComplex]).counts()[0]
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
	bulkMolCntr.countsIs(ribosome30S - ribosomeToActivate, [sim_data.moleculeIds.s30_fullComplex])
	bulkMolCntr.countsIs(ribosome50S - ribosomeToActivate, [sim_data.moleculeIds.s50_fullComplex])


def setDaughterInitialConditions(sim, sim_data):
	inherited_state_path = sim._inheritedStatePath
	assert inherited_state_path is not None

	inherited_state = load_inherited_state(inherited_state_path)

	sim._isDead = inherited_state['is_dead']

	elngRate = inherited_state['elng_rate']
	elng_rate_factor = inherited_state['elng_rate_factor']
	if sim._growthRateNoise:
		sim.processes["PolypeptideElongation"].setElngRate = elngRate
		sim.processes["PolypeptideElongation"].elngRateFactor = elng_rate_factor

	sim.internal_states["BulkMolecules"].loadSnapshot(inherited_state['bulk_molecules'])
	sim.internal_states["UniqueMolecules"].loadSnapshot(inherited_state['unique_molecules'])

	sim._initialTime = inherited_state['initial_time']


def determineChromosomeState(C, D, tau, replication_length):
	"""
	Calculates the number and position of replicating DNA polymerases at the
	beginning of the cell cycle.

	Inputs
	--------
	- C: the C period of the cell, the length of time between replication
	initiation and replication completion.
	- D: the D period of the cell, the length of time between completing
	replication of the chromosome and division of the cell.
	- tau: the doubling time of the cell
	- replication_length: the amount of DNA to be replicated per fork, usually
	half of the genome, in base-pairs

	Returns
	--------
	- n_oric: the number of OriC's in the cell at initiation.
	- sequenceIdx: an index for each of the four types/directions of DNA
	replication - leading and lagging strand of the forward and the reverse
	fork = 4 total. This vector is always simply [0,1,2,3] repeated once for
	each replication event. i.e. for three active replication events (6 forks,
	12 polymerases) sequenceIdx = [0,1,2,3,0,1,2,3,0,1,2,3]
	- sequenceLength: the position in the genome that each polymerase
	referenced in sequenceIdx has reached, in base-pairs. This is handled such
	that even though in reality some polymerases are replicating in different
	directions,	all values here are calculated as though each starts at 0 and
	goes up to the total number of base-pairs to be replicated.
	- replicationRound: an integer stating in which replication round the
	polymerase referenced by sequenceIdx has been initiated. Each time all
	origins of replication in the cell fire, a new replication round has
	started. This array is integer-valued, and counts from 0 (the oldest round)
	up to n (the most recent round).
	- chromosomeIndex: indicator variable for which chromosome the polymerases
	are associated with and therefore which daughter cell should inherit each
	polymerase. Since there is only one chromosome initially, all indexes are
	set to zero.

	Notes
	--------
	If NO polymerases are active at the start of the cell cycle, equivalent to
	the C + D periods being shorter than the doubling rate tau, then this
	function returns empty lists. dnaPoly.attrIs() should not be run in this
	case, as no DNA replication will be underway.
	"""

	# All inputs must be positive numbers
	assert C.asNumber(units.min) >= 0, "C value can't be negative."
	assert D.asNumber(units.min) >= 0, "D value can't be negative."
	assert tau.asNumber(units.min) >= 0, "tau value can't be negative."
	assert replication_length.asNumber(units.nt) >= 0, "replication_length value can't be negative."

	# Require that D is shorter than tau - time between completing DNA
	# replication and cell division must be shorter than the time between two
	# cell divisions.

	assert D.asNumber(units.min) < tau.asNumber(units.min), "The D period must be shorter than the doubling time tau."

	# Calculate the number of active replication rounds
	n_round = int(np.floor(
		(C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min)))

	# Initialize arrays to be returned
	sequenceIdx = []
	sequenceLength = []
	replicationRound = []
	chromosomeIndex = []

	# Loop through active replication rounds, starting from the oldest round.
	# If n_round = 0 skip loop entirely - no active replication round.
	for n in xrange(n_round):
		# Determine at which location (base) of the chromosome the replication
		# forks should be initialized to
		rel_location = 1 - (((n + 1)*tau - D)/C)
		rel_location = units.convertNoUnitToNumber(rel_location)
		fork_location = np.floor(rel_location*(
			replication_length.asNumber(units.nt)))

		# Add 2^n initiation events per round. A single initiation event
		# generates two replication forks and four elongating strands
		# (polymerases).
		n_event = 2**n

		# sequenceIdx refers to the type of the strands that the polymerase is
		# elongating - i.e. forward and reverse, lagging and leading strands.
		sequenceIdx += [0, 1, 2, 3] * n_event

		# sequenceLength refers to how far along the chromosome the polymerases
		# have elongated to. All four are assumed to have elongated up to the
		# replication fork.
		sequenceLength += [fork_location] * (4*n_event)

		# replicationRound is the index of the replication round that the
		# polymerases belong to - for each replication round, all origins in
		# the cell are fired simultaneously, and the initiated polymerases
		# share the same round index
		replicationRound += [n] * (4*n_event)

		# chromosomeIndex indicates which daughter cell will inherit the
		# polymerase. Since there is only one initial chromosome, all
		# polymerases are initially given index zero.
		chromosomeIndex += [0] * (4*n_event)

	# Convert to numpy arrays
	sequenceIdx = np.array(sequenceIdx, dtype=np.int8)
	sequenceLength = np.array(sequenceLength, dtype=np.int64)
	replicationRound = np.array(replicationRound, dtype=np.int64)
	chromosomeIndex = np.array(chromosomeIndex, dtype=np.int64)

	return sequenceIdx, sequenceLength, replicationRound, chromosomeIndex


def determineOriCState(C, D, tau):
	"""
	Calculates the number of OriC's in a cell upon initiation and the indexes
	of chromosomes that the OriC's belong to, determined by the replication
	state of the chromosome.

	Inputs
	--------
	- C: the C period of the cell, the length of time between replication
	initiation and replication completion.
	- D: the D period of the cell, the length of time between completing
	replication of the chromosome and division of the cell.
	- tau: the doubling time of the cell

	Returns
	--------
	- n_oric: the number of OriC's in the cell at initiation.
	- chromosomeIndex: indicator variable for which chromosome the oriC's are
	associated with and therefore which daughter cell should inherit each oriC.
	Since there is only one chromosome initially, all indexes are set to zero.
	"""

	# Number active replication generations (can be many initiations per gen.)
	n_round = int(np.floor(
		(C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min)))
	n_oric = 2**n_round
	chromosomeIndex = np.zeros(n_oric, dtype=np.int)

	return n_oric, chromosomeIndex
