#!/usr/bin/env python

"""
PolypeptideElongation

Translation elongation sub-model.

TODO:
- see the initiation process for more TODOs

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

MICROMOLAR_UNITS = units.umol / units.L


class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideElongation, self).initialize(sim, sim_data)

		constants = sim_data.constants
		translation = sim_data.process.translation
		transcription = sim_data.process.transcription

		# Load parameters
		self.nAvogadro = constants.nAvogadro
		self.cellDensity = constants.cellDensity
		self.aaNames = sim_data.moleculeGroups.aaIDs
		proteinIds = translation.monomerData['id']
		self.proteinLengths = translation.monomerData["length"].asNumber()
		self.proteinSequences = translation.translationSequences
		self.aaWeightsIncorporated = translation.translationMonomerWeights
		self.endWeight = translation.translationEndWeight
		self.gtpPerElongation = constants.gtpPerTranslation

		self.maxRibosomeElongationRate = float(constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		self.ribosomeElongationRateDict = translation.ribosomeElongationRateDict

		self.translation_aa_supply = sim_data.translationSupplyRate

		# Used for figure in publication
		self.trpAIndex = np.where(proteinIds == "TRYPSYN-APROTEIN[c]")[0][0]

		# Create view onto activly elongating 70S ribosomes
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# Create views onto 30S and 70S ribosomal subunits for termination
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Create view onto all proteins
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		# Create views onto all polymerization reaction small molecules
		self.aas = self.bulkMoleculesView(self.aaNames)
		self.water = self.bulkMoleculeView('WATER[c]')
		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.proton = self.bulkMoleculeView("PROTON[c]")

		# Set for timestep calculation
		self.gtpUsed = 0
		self.gtpAvailable = 0

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Simulation options
		self.translationSupply = sim._translationSupply
		self.use_trna_charging = sim._trna_charging

		self.elngRateFactor = 1.

		# Data structures for charging
		self.aa_from_synthetase = transcription.aa_from_synthetase
		self.aa_from_trna = transcription.aa_from_trna
		self.charging_stoich_matrix = transcription.charging_stoich_matrix()

		# Names of molecules associated with tRNA charging
		self.uncharged_trna_names = transcription.rnaData['id'][transcription.rnaData['isTRna']]
		self.charged_trna_names = transcription.charged_trna_names
		self.charging_molecule_names = transcription.charging_molecules
		self.synthetase_names = transcription.synthetase_names

		# Create views for tRNA charging molecules
		self.uncharged_trna = self.bulkMoleculesView(self.uncharged_trna_names)
		self.charged_trna = self.bulkMoleculesView(self.charged_trna_names)
		self.charging_molecules = self.bulkMoleculesView(self.charging_molecule_names)
		self.synthetases = self.bulkMoleculesView(self.synthetase_names)

		# ppGpp parameters for tRNA charging and ribosome elongation
		self.kS = constants.synthetase_charging_rate.asNumber(1 / units.s)
		self.KMtf = constants.Km_synthetase_uncharged_trna.asNumber(MICROMOLAR_UNITS)
		self.KMaa = constants.Km_synthetase_amino_acid.asNumber(MICROMOLAR_UNITS)
		self.krta = constants.Kdissociation_charged_trna_ribosome.asNumber(MICROMOLAR_UNITS)
		self.krtf = constants.Kdissociation_uncharged_trna_ribosome.asNumber(MICROMOLAR_UNITS)

	def calculateRequest(self):
		# Set ribosome elongation rate based on simulation medium environment and elongation rate factor
		# which is used to create single-cell variability in growth rate
		# The maximum number of amino acids that can be elongated in a single timestep is set to 22 intentionally as the minimum number of padding values
		# on the protein sequence matrix is set to 22. If timesteps longer than 1.0s are used, this feature will lead to errors in the effective ribosome
		# elongation rate.

		current_nutrients = self._external_states['Environment'].nutrients

		if self.translationSupply:
			self.ribosomeElongationRate = np.min([self.maxRibosomeElongationRate, int(stochasticRound(self.randomState,
				self.maxRibosomeElongationRate * self.timeStepSec()))]) # Will be set to maxRibosomeElongationRate if timeStepSec > 1.0s
		else:
			self.ribosomeElongationRate = np.min([22, int(stochasticRound(self.randomState,
				self.elngRateFactor * self.ribosomeElongationRateDict[current_nutrients].asNumber(units.aa / units.s) * self.timeStepSec()))])

		# Request all active ribosomes
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		if len(activeRibosomes) == 0:
			return

		# Build sequences to request appropriate amount of amino acids to
		# polymerize for next timestep
		proteinIndexes, peptideLengths = activeRibosomes.attrs(
					'proteinIndex', 'peptideLength'
					)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.ribosomeElongationRate,
			)

		sequenceHasAA = (sequences != polymerize.PAD_VALUE)
		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		if self.use_trna_charging:
			# Conversion from counts to molarity
			cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
			cell_volume = cell_mass / self.cellDensity
			counts_to_molar = 1 / (self.nAvogadro * cell_volume)

			# Get counts and convert synthetase and tRNA to a per AA basis
			synthetase_counts = np.dot(self.aa_from_synthetase, self.synthetases.total())
			aa_counts = self.aas.total()
			uncharged_trna_counts = np.dot(self.aa_from_trna, self.uncharged_trna.total())
			charged_trna_counts = np.dot(self.aa_from_trna, self.charged_trna.total())
			ribosome_counts = len(self.activeRibosomes.allMolecules())

			# Get concentration
			f = aasInSequences / aasInSequences.sum()
			synthetase_conc = counts_to_molar * synthetase_counts
			aa_conc = counts_to_molar * aa_counts
			uncharged_trna_conc = counts_to_molar * uncharged_trna_counts
			charged_trna_conc = counts_to_molar * charged_trna_counts
			ribosome_conc = counts_to_molar * ribosome_counts

			# Calculate steady state tRNA levels and resulting elongation rate
			fraction_charged, v_rib = self.calculate_trna_charging(
				synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f, self.timeStepSec()
				)

			aa_counts_for_translation = v_rib * f * self._sim.timeStepSec() / counts_to_molar.asNumber(MICROMOLAR_UNITS)

			total_trna = self.charged_trna.total() + self.uncharged_trna.total()
			final_charged_trna = np.dot(fraction_charged, self.aa_from_trna * total_trna)

			charged_trna_request = self.charged_trna.total() - final_charged_trna
			charged_trna_request[charged_trna_request < 0] = 0
			uncharged_trna_request = final_charged_trna - self.charged_trna.total()
			uncharged_trna_request[uncharged_trna_request < 0] = 0

			self.aa_counts_for_translation = np.array(aa_counts_for_translation)

			fraction_trna_per_aa = total_trna / np.dot(np.dot(self.aa_from_trna, total_trna), self.aa_from_trna)
			total_charging_reactions = np.dot(aa_counts_for_translation, self.aa_from_trna) * fraction_trna_per_aa + uncharged_trna_request

			# Only request molecules that will be consumed in the charging reactions
			requested_molecules = -np.dot(self.charging_stoich_matrix, total_charging_reactions)
			requested_molecules[requested_molecules < 0] = 0
			self.charging_molecules.requestIs(requested_molecules)

			# Request charged tRNA that will become uncharged
			self.charged_trna.requestIs(charged_trna_request)

			# Request water for transfer of AA from tRNA for initial polypeptide.
			# This is severe overestimate assuming the worst case that every
			# elongation is initializing a polypeptide. This excess of water
			# shouldn't matter though.
			self.water.requestIs(aa_counts_for_translation.sum())
		else:
			if self.translationSupply:
				translationSupplyRate = self.translation_aa_supply[current_nutrients] * self.elngRateFactor

				self.writeToListener("RibosomeData", "translationSupply", translationSupplyRate.asNumber())

				dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

				molAasRequested = translationSupplyRate * dryMass * self.timeStepSec() * units.s

				aa_counts_for_translation = units.convertNoUnitToNumber(molAasRequested * self.nAvogadro)

				aa_counts_for_translation = np.fmin(aa_counts_for_translation, aasInSequences) # Check if this is required. It is a better request but there may be fewer elongations.
			else:
				aa_counts_for_translation = aasInSequences

			self.aas.requestIs(aa_counts_for_translation)

			# Not modeling charging so set fraction charged to 0 for all tRNA
			fraction_charged = np.zeros(len(self.aaNames))

		self.writeToListener("GrowthLimits", "fraction_trna_charged", np.dot(fraction_charged, self.aa_from_trna))
		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total())
		self.writeToListener("GrowthLimits", "aaRequestSize", aa_counts_for_translation)

		# Request GTP for polymerization based on sequences
		gtpsHydrolyzed = np.int64(np.ceil(self.gtpPerElongation * aa_counts_for_translation.sum()))

		self.writeToListener("GrowthLimits", "gtpPoolSize", self.gtp.total()[0])
		self.writeToListener("GrowthLimits", "gtpRequestSize", gtpsHydrolyzed)

		# GTP hydrolysis is carried out in Metabolism process for growth associated maintenence
		# THis is set here for metabolism to use
		self.gtpRequest = gtpsHydrolyzed

	def evolveState(self):
		# Write allocation data to listener
		self.writeToListener("GrowthLimits", "gtpAllocated", self.gtp.count())
		self.writeToListener("GrowthLimits", "aaAllocated", self.aas.counts())

		# Get number of active ribosomes
		activeRibosomes = self.activeRibosomes.molecules()

		self.writeToListener("GrowthLimits", "activeRibosomeAllocated", len(activeRibosomes))

		if len(activeRibosomes) == 0:
			return

		# Build amino acids sequences for each ribosome to polymerize
		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiff_protein'
			)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.ribosomeElongationRate,
			)

		if sequences.size == 0:
			return

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != polymerize.PAD_VALUE)])
		total_aa_counts = self.aas.counts()
		if self.use_trna_charging:
			aa_counts_for_translation = np.fmin(total_aa_counts, self.aa_counts_for_translation)
		else:
			aa_counts_for_translation = total_aa_counts

		# Using polymerization algorithm elongate each ribosome up to the limits
		# of amino acids, sequence, and GTP
		result = polymerize(
			sequences,
			aa_counts_for_translation,
			10000000, # Set to a large number, the limit is now taken care of in metabolism
			self.randomState
			)

		sequenceElongations = result.sequenceElongation
		aas_used = result.monomerUsages
		nElongations = result.nReactions

		# Update masses of ribosomes attached to polymerizing polypeptides
		massIncreaseProtein = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.aaWeightsIncorporated
			)

		updatedLengths = peptideLengths + sequenceElongations

		didInitialize = (
			(sequenceElongations > 0) &
			(peptideLengths == 0)
			)

		updatedMass = massDiffProtein + massIncreaseProtein

		updatedMass[didInitialize] += self.endWeight

		# Write current average elongation to listener
		currElongRate = (sequenceElongations.sum() / len(activeRibosomes)) / self.timeStepSec()
		self.writeToListener("RibosomeData", "effectiveElongationRate", currElongRate)

		# Update active ribosomes, terminating if neccessary
		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiff_protein = updatedMass
			)

		# Ribosomes that reach the end of their sequences are terminated and
		# dissociated into 30S and 50S subunits. The polypeptide that they are polymerizing
		# is converted into a protein in BulkMolecules
		terminalLengths = self.proteinLengths[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])
		self.bulkMonomers.countsInc(terminatedProteins)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

		if self.use_trna_charging:
			# Get tRNA counts
			uncharged_trna = self.uncharged_trna.counts()
			charged_trna = self.charged_trna.counts()
			total_trna = uncharged_trna + charged_trna

			# Adjust molecules for number of charging reactions that occurred
			## Net charged is tRNA that can be charged minus allocated charged tRNA for uncharging
			aa_for_charging = total_aa_counts - aas_used
			n_aa_charged = np.fmin(aa_for_charging, np.dot(self.aa_from_trna, uncharged_trna))
			n_trna_charged = self.distribution_from_aa(n_aa_charged, uncharged_trna, True)
			net_charged = n_trna_charged - charged_trna

			## Reactions that are charged and elongated in same time step
			charged_and_elongated = self.distribution_from_aa(aas_used, total_trna)
			total_charging_reactions = charged_and_elongated + net_charged
			self.charging_molecules.countsInc(np.dot(self.charging_stoich_matrix, total_charging_reactions))

			## Account for uncharging of tRNA during elongation
			self.charged_trna.countsDec(charged_and_elongated)
			self.uncharged_trna.countsInc(charged_and_elongated)

			# Update proton counts to reflect polymerization reactions and transfer of AA from tRNA
			# Peptide bond formation releases a water but transferring AA from tRNA consumes a OH-
			# Net production of H+ for each elongation, consume extra water for each initialization
			# since a peptide bond doesn't form
			self.proton.countInc(nElongations)
			self.water.countDec(nInitialized)
		else:
			# Update counts of amino acids and water to reflect polymerization reactions
			self.aas.countsDec(aas_used)
			self.water.countInc(nElongations - nInitialized)
			net_charged = np.zeros(len(self.uncharged_trna_names))

		# Write data to listeners
		self.writeToListener("GrowthLimits", "net_charged", net_charged)
		self.writeToListener("GrowthLimits", "aasUsed", aas_used)
		self.writeToListener("GrowthLimits", "gtpUsed", self.gtpUsed)

		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aa_counts_for_translation)

		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())
		self.writeToListener("RibosomeData", "actualElongationHist", np.histogram(sequenceElongations, bins = np.arange(0,23))[0])
		self.writeToListener("RibosomeData", "elongationsNonTerminatingHist", np.histogram(sequenceElongations[~didTerminate], bins=np.arange(0,23))[0])

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptideLengths)[didTerminate].sum())
		self.writeToListener("RibosomeData", "numTrpATerminated", terminatedProteins[self.trpAIndex])

		self.writeToListener("RibosomeData", "processElongationRate", self.ribosomeElongationRate / self.timeStepSec())

	def calculate_trna_charging(self, synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f, time_limit=None):
		'''
		Calculates the steady state value of tRNA based on charging and incorporation through polypeptide elongation.
		The fraction of charged/uncharged is also used to determine how quickly the ribosome is elongating.

		Inputs:
			synthetase_conc (array of floats with concentration units) - concentration of synthetases associated
				with each amino acid
			uncharged_trna_conc (array of floats with concentration units) - concentration of uncharged tRNA associated
				with each amino acid
			charged_trna_conc (array of floats with concentration units) - concentration of charged tRNA associated
				with each amino acid
			aa_conc (array of floats with concentration units) - concentration of each amino acid
			ribosome_conc (float with concentration units) - concentration of active ribosomes
			f (array of floats) - fraction of each amino acid to be incorporated to total amino acids incorporated
			time_limit (float) - if not None, time limit to reach steady state, will warn if reached

		Returns:
			fraction_charged (array of floats) - fraction of total tRNA that is charged for each tRNA species
			v_rib (float) - ribosomal elongation rate in units of uM/s
		'''

		synthetase_conc = synthetase_conc.asNumber(MICROMOLAR_UNITS)
		uncharged_trna_conc = uncharged_trna_conc.asNumber(MICROMOLAR_UNITS)
		charged_trna_conc = charged_trna_conc.asNumber(MICROMOLAR_UNITS)
		aa_conc = aa_conc.asNumber(MICROMOLAR_UNITS)
		ribosome_conc = ribosome_conc.asNumber(MICROMOLAR_UNITS)

		# Solve to steady state with short time steps
		t = 0
		dt = 0.001
		diff = 1
		while diff > 1e-3:
			v_charging = (self.kS * synthetase_conc * uncharged_trna_conc * aa_conc
				/ (self.KMaa * self.KMtf *
				(1 + uncharged_trna_conc/self.KMtf + aa_conc/self.KMaa + uncharged_trna_conc*aa_conc/self.KMtf/self.KMaa))
				)
			numerator_ribosome = 1 + np.sum(f * (self.krta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * self.krta / self.krtf))
			v_rib = self.maxRibosomeElongationRate * ribosome_conc / numerator_ribosome

			# Handle case when f is 0 and charged_trna_conc is 0
			if not np.isfinite(v_rib):
				v_rib = 0

			delta_conc = (v_charging - v_rib*f) * dt
			uncharged_trna_conc -= delta_conc
			charged_trna_conc += delta_conc
			diff = np.linalg.norm(delta_conc)

			t += dt
			if time_limit is not None and t >= time_limit:
				print('Warning: time limit reached for tRNA charging, norm is {:0.4f}'.format(diff))
				break

		fraction_charged = charged_trna_conc / (uncharged_trna_conc + charged_trna_conc)
		return fraction_charged, v_rib

	def distribution_from_aa(self, n_aa, n_trna, limited=False):
		'''
		Distributes counts of amino acids to tRNAs that are associated with each amino acid.
		Uses self.aa_from_trna mapping to distribute from amino acids to tRNA based on the
		fraction that each tRNA species makes up for all tRNA species that code for the
		same amino acid.

		Inputs:
			n_aa (array of ints) - counts of each amino acid to distribute to each tRNA
			n_trna (array of ints) - counts of each tRNA to determine the distribution
			limited (bool) - optional, if True, limits the amino acids distributed to
				each tRNA to the number of tRNA that are available (n_trna)

		Returns:
			array of ints - distributed counts for each tRNA
		'''

		# Determine the fraction each tRNA species makes up out of all tRNA of the
		# associated amino acid
		f_trna = n_trna / np.dot(np.dot(self.aa_from_trna, n_trna), self.aa_from_trna)
		f_trna[~np.isfinite(f_trna)] = 0

		trna_counts = np.zeros(f_trna.shape, np.int64)
		for count, row in izip(n_aa, self.aa_from_trna):
			idx = (row == 1)
			frac = f_trna[idx]

			counts = np.floor(frac * count)
			diff = int(count - counts.sum())

			# Add additional counts to get up to counts to distribute
			# Prevent adding over the number of tRNA available if limited
			if diff > 0:
				if limited:
					for _ in range(diff):
						frac[(n_trna[idx] - counts) == 0] = 0
						frac /= frac.sum()  # normalize for multinomial distribution
						adjustment = self.randomState.multinomial(1, frac)
						counts += adjustment
				else:
					adjustment = self.randomState.multinomial(diff, frac)
					counts += adjustment

			trna_counts[idx] = counts

		return trna_counts

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		"""
		Assumes GTP is the readout for failed translation with respect to the timestep.
		"""

		# Until more padding values are added to the protein sequence matrix, limit the maximum timestep length to 1 second
		# Since the current upper limit on a.a's elongated by ribosomes during a single timestep is set to 22, timesteps
		# longer than 1.0s do not lead to errors, but does slow down the ribosome elongation rate of the resulting simulation.
		# Must be modified if timesteps longer than 1.0s are desired.
		if inputTimeStep > 1.0:
			return False

		activeRibosomes = float(self.activeRibosomes.total()[0])
		self.gtpAvailable = float(self.gtp.total()[0])

		# Without an estimate on ribosome counts, require a short timestep until estimates available
		if activeRibosomes == 0:
			if inputTimeStep <= .2:
				return True
			else:
				return False

		dt = inputTimeStep * timeStepSafetyFraction
		gtpExpectedUsage = activeRibosomes * self.ribosomeElongationRate * self.gtpPerElongation * dt

		if gtpExpectedUsage < self.gtpAvailable:
			return True
		else:
			return False

	def wasTimeStepShortEnough(self):
		"""
		If translation used more than 90 percent of gtp, timeStep was too short.
		"""

		# If gtpAvailable is 0 and the timeStep is short, use the gtp produced this timeStep as the estimate
		if self.gtpAvailable == 0 and self.timeStepSec() <= .2:
			self.gtpAvailable = self.gtp.total()[0]

		if (self.gtpAvailable * .9) < self.gtpUsed:
			return False
		else:
			return True
