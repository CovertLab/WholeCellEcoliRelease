"""
TODO: document math
TODO: raise/warn if physiological metabolite concentration targets appear to be smaller than what
  is needed at this time step size
"""

from __future__ import absolute_import, division, print_function

from itertools import izip
import scipy.sparse

import numpy as np
from arrow import StochasticSystem

from wholecell.utils.fitting import normalize, countsFromMassAndExpression, masses_and_counts_for_homeostatic_target
from wholecell.utils.polymerize import computeMassIncrease
from wholecell.utils import units
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.sim.divide_cell import load_inherited_state

RAND_MAX = 2**31

def calcInitialConditions(sim, sim_data):
	'''Calculate the initial conditions for a new cell without inherited state
	from a parent cell.

	Params:
		sim: The simulation object, which must not have an _inheritedStatePath,
			otherwise the simulation should use setDaughterInitialConditions()
		sim_data: to initialize state from
	'''
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

	# Set other biomass components
	initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff)

	# Form complexes
	initializeComplexation(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	# Initialize counts of full chromosomes
	initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data)

	# Initialize unique molecules relevant to replication
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

	stoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64)

	time_step = 1
	seed = randomState.randint(RAND_MAX)
	complexation_rates = sim_data.process.complexation.rates
	system = StochasticSystem(stoichMatrix.T, complexation_rates, random_seed=seed)

	# form complexes until no new complexes form (some complexes are complexes of complexes)
	while True:
		moleculeCounts = moleculeView.counts()
		complexation_result = system.evolve(time_step, moleculeCounts)

		updatedMoleculeCounts = complexation_result['outcome']
		complexationEvents = complexation_result['occurrences']

		bulkMolCntr.countsIs(
			updatedMoleculeCounts,
			moleculeNames)

		if not np.any(moleculeCounts - updatedMoleculeCounts):
			break

	bulkMolCntr.countsIs(rnaseCounts, rnases)


def initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	Initializes the counts of full chromosomes to one. The division_time of
	this initial chromosome is set to be zero for consistency.
	"""
	full_chromosome = uniqueMolCntr.objectsNew("fullChromosome", 1)
	full_chromosome.attrIs(
		division_time = 0.0,
		has_induced_division = True,
		domain_index = 0,
		)


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

	# Calculate length of replichore
	genome_length = sim_data.process.replication.genome_length
	replichore_length = np.ceil(0.5*genome_length) * units.nt

	# Generate arrays specifying appropriate initial replication conditions
	oric_state, replisome_state, domain_state = determine_chromosome_state(
		C, D, tau, replichore_length, sim_data.process.replication.no_child_place_holder
		)

	n_oric = oric_state["domain_index"].size
	n_replisome = replisome_state["domain_index"].size
	n_domain = domain_state["domain_index"].size

	# Add OriC molecules with the proposed attributes
	oriC = uniqueMolCntr.objectsNew('originOfReplication', n_oric)
	oriC.attrIs(
		domain_index = oric_state["domain_index"],
		)

	# Add chromosome domain molecules with the proposed attributes
	chromosome_domains = uniqueMolCntr.objectsNew("chromosome_domain", n_domain)
	chromosome_domains.attrIs(
		domain_index = domain_state["domain_index"],
		child_domains = domain_state["child_domains"],
		)

	if n_replisome != 0:
		# Update mass to account for DNA strands that have already been
		# elongated.
		sequences = sim_data.process.replication.replication_sequences
		fork_coordinates = replisome_state["coordinates"]
		sequence_elongations = np.abs(np.repeat(fork_coordinates, 2))

		mass_increase_dna = computeMassIncrease(
				np.tile(sequences, (n_replisome//2, 1)),
				sequence_elongations,
				sim_data.process.replication.replicationMonomerWeights.asNumber(units.fg)
				)

		# Add active replisomes as unique molecules and set attributes
		active_replisomes = uniqueMolCntr.objectsNew(
			'active_replisome', n_replisome
			)
		active_replisomes.attrIs(
			coordinates = replisome_state["coordinates"],
			right_replichore = replisome_state["right_replichore"],
			domain_index = replisome_state["domain_index"],
			massDiff_DNA = mass_increase_dna[0::2] + mass_increase_dna[1::2],
		)

		# Remove replisome subunits from bulk molecules
		bulkMolCntr.countsDec(3*n_replisome, sim_data.moleculeGroups.replisome_trimer_subunits)
		bulkMolCntr.countsDec(n_replisome, sim_data.moleculeGroups.replisome_monomer_subunits)

	# Initialize gene dosage
	geneCopyNumberColNames = sim_data.process.transcription_regulation.geneCopyNumberColNames
	geneCopyNumberView = bulkMolCntr.countsView(geneCopyNumberColNames)
	replicationCoordinate = sim_data.process.transcription.rnaData["replicationCoordinate"]

	# Set all copy numbers to one initially
	initialGeneCopyNumber = np.ones(len(geneCopyNumberColNames))

	# Get coordinates of forks in both directions
	forward_fork_coordinates = replisome_state["coordinates"][
		replisome_state["right_replichore"]
	]
	reverse_fork_coordinates = replisome_state["coordinates"][
		np.logical_not(replisome_state["right_replichore"])
	]

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
	setIdxs = isRRna | isTRna | isRProtein | isRnap

	# Calculate synthesis probabilities based on transcription regulation
	rnaSynthProb = recruitmentMatrix.dot(recruitmentView)
	if len(genetic_perturbations) > 0:
		rnaSynthProb[genetic_perturbations['fixedRnaIdxs']] = genetic_perturbations['fixedSynthProbs']

	# Adjust probabilities to not be negative
	rnaSynthProb[rnaSynthProb < 0] = 0.0
	rnaSynthProb /= rnaSynthProb.sum()
	if np.any(rnaSynthProb < 0):
		raise Exception("Have negative RNA synthesis probabilities")

	# Adjust synthesis probabilities depending on environment
	synthProbFractions = rnaSynthProbFractions[currentNutrients]

	# Allocate synthesis probabilities based on type of RNA
	rnaSynthProb[isMRna] *= synthProbFractions['mRna'] / rnaSynthProb[isMRna].sum()
	rnaSynthProb[isTRna] *= synthProbFractions['tRna'] / rnaSynthProb[isTRna].sum()
	rnaSynthProb[isRRna] *= synthProbFractions['rRna'] / rnaSynthProb[isRRna].sum()

	# Set fixed synthesis probabilities for RProteins and RNAPs
	rnaSynthProb[isRProtein] = rnaSynthProbRProtein[currentNutrients]
	rnaSynthProb[isRnap] = rnaSynthProbRnaPolymerase[currentNutrients]

	assert rnaSynthProb[setIdxs].sum() < 1.0

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
	'''Calculate the initial conditions for a new cell from state inherited
	from a parent cell, stored in a file.

	Params:
		sim: The simulation object, which must have an _inheritedStatePath to
			name the file directory to load the state from, otherwise the
			simulation should use calcInitialConditions()
		sim_data: Unused argument needed to conform to an informal interface
	'''
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


def determine_chromosome_state(C, D, tau, replichore_length, place_holder):
	"""
	Calculates the attributes of oriC's, replisomes, and chromosome domains on
	the chromosomes at the beginning of the cell cycle.

	Inputs
	--------
	- C: the C period of the cell, the length of time between replication
	initiation and replication termination.
	- D: the D period of the cell, the length of time between completing
	replication of the chromosome and division of the cell.
	- tau: the doubling time of the cell
	- replichore_length: the amount of DNA to be replicated per fork, usually
	half of the genome, in base-pairs
	- place_holder: placeholder value for chromosome domains without child
	domains

	Returns
	--------
	oric_state: dictionary of attributes for the oriC molecules with the
	following keys.
	- domain_index: a vector of integers indicating which chromosome domain the
	oriC sequence belongs to.

	replisome_state: dictionary of attributes for the replisome molecules
	with the following keys.
	- coordinates: a vector of integers that indicates where the replisomes
	are located on the chromosome relative to the origin, in base pairs.
	- right_replichore: a vector of boolean values that indicates whether the
	replisome is on the right replichore (True) or the left replichore (False).
	- domain_index: a vector of integers indicating which chromosome domain the
	replisomes belong to. The index of the "mother" domain of the replication
	fork is assigned to the replisome.

	domain_state: dictionary of attributes for the chromosome domains with the
	following keys.
	- domain_index: the indexes of the domains.
	- child_domains: the (n_domain X 2) array of the domain indexes of the two
	children domains that are connected on the oriC side with the given domain.
	"""

	# All inputs must be positive numbers
	assert C.asNumber(units.min) >= 0, "C value can't be negative."
	assert D.asNumber(units.min) >= 0, "D value can't be negative."
	assert tau.asNumber(units.min) >= 0, "tau value can't be negative."
	assert replichore_length.asNumber(units.nt) >= 0, "replichore_length value can't be negative."

	# Require that D is shorter than tau - time between completing DNA
	# replication and cell division must be shorter than the time between two
	# cell divisions.
	assert D.asNumber(units.min) < tau.asNumber(units.min), "The D period must be shorter than the doubling time tau."

	# Calculate the number of active replication rounds
	n_round = int(np.floor(
		(C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min)))

	# Initialize arrays for replisomes
	n_replisomes = 2*(2**n_round - 1)
	coordinates = np.zeros(n_replisomes, dtype=np.int64)
	right_replichore_replisome = np.zeros(n_replisomes, dtype=np.bool)
	domain_index_replisome = np.zeros(n_replisomes, dtype=np.int32)

	# Initialize child domain array for chromosome domains
	n_domains = 2**(n_round + 1) - 1
	child_domains = np.full((n_domains, 2), place_holder, dtype=np.int32)

	# Set domain_index attribute of oriC's and chromosome domains
	domain_index_oric = np.arange(2**n_round - 1, 2**(n_round + 1) - 1, dtype=np.int32)
	domain_index_domains = np.arange(0, n_domains, dtype=np.int32)

	def n_events_before_this_round(round_idx):
		"""
		Calculates the number of replication events that happen before the
		replication round index given as an argument. Since 2**i events happen
		at each round i = 0, 1, ..., the sum of the number of events before
		round j is 2**j - 1.
		"""
		return 2**round_idx - 1

	# Loop through active replication rounds, starting from the oldest round.
	# If n_round = 0 skip loop entirely - no active replication round.
	for round_idx in np.arange(n_round):
		# Determine at which location (base) of the chromosome the replication
		# forks should be initialized to
		rel_location = 1.0 - (((round_idx + 1.0)*tau - D)/C)
		rel_location = units.convertNoUnitToNumber(rel_location)
		fork_location = np.floor(rel_location*(
			replichore_length.asNumber(units.nt)))

		# Add 2^n initiation events per round. A single initiation event
		# generates two replication forks.
		n_events_this_round = 2**round_idx

		# Set attributes of replisomes for this replication round
		coordinates[
			2*n_events_before_this_round(round_idx):
			2*n_events_before_this_round(round_idx + 1)
			] = np.tile(np.array([fork_location, -fork_location]), n_events_this_round)

		right_replichore_replisome[
			2*n_events_before_this_round(round_idx):
			2*n_events_before_this_round(round_idx + 1)
			] = np.tile(np.array([True, False]), n_events_this_round)

		for i, domain_index in enumerate(
				np.arange(n_events_before_this_round(round_idx),
					n_events_before_this_round(round_idx+1))):
			domain_index_replisome[
				2*n_events_before_this_round(round_idx) + 2*i:
				2*n_events_before_this_round(round_idx) + 2*(i+1)
				] = np.repeat(domain_index, 2)

		# Set attributes of chromosome domains for this replication round
		for i, domain_index in enumerate(
				np.arange(n_events_before_this_round(round_idx + 1),
					n_events_before_this_round(round_idx + 2), 2)):
			child_domains[
				n_events_before_this_round(round_idx) + i, :
				] = np.array([domain_index, domain_index + 1])

	# Convert to numpy arrays and wrap into dictionaries
	oric_state = {
		"domain_index": domain_index_oric,
		}

	replisome_state = {
		"coordinates": coordinates,
		"right_replichore": right_replichore_replisome,
		"domain_index": domain_index_replisome,
		}

	domain_state = {
		"child_domains": child_domains,
		"domain_index": domain_index_domains,
		}

	return oric_state, replisome_state, domain_state
