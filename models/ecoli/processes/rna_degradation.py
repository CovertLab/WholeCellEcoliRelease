"""
Submodel for RNA degradation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/26/2015 - Updated 22/3/2017

Mathematical formulation:

dr/dt = Kb - Kd * r
or,

dr/dt = Kb - kcatEndoRNase * EndoRNase * r / (Km + r)
or,

dr/dt = Kb - kcatEndoRNase * EndoRNase * r/Km / (1 + Sum(r/Km))

	where	r = RNA counts
			Kb = RNA production given a RNAP synthesis rate 
			tau = doubling time
			kcatEndoRNase = enzymatic activity for EndoRNases
			kd = RNA degradation rates 
			Km = Michaelis-Menten constants fitted to recapitulate first-order
			RNA decay:
				kd * r = kcatEndoRNase * EndoRNase * r / (Km + r),
				    non-cooperative EndoRNases
				kd * r = kcatEndoRNase * EndoRNase * r/Km / (1 + sum(r/Km)),
				    cooperation

This sub-model encodes molecular simulation of RNA degradation as two main
steps guided by RNases, "endonucleolytic cleavage" and "exonucleolytic
digestion":

1. Compute total counts of RNA to be degraded (D) and total capacity for
endo-cleavage (C) at each time point
2. Evaluate C and D. If C > D, then define a fraction of active endoRNases 
3. Dissect RNA degraded into different species (mRNA, tRNA, and rRNA) by
accounting endoRNases specificity
4. Update RNA fragments (assumption: fragments are represented as a pool of
nucleotides) created because of endonucleolytic cleavage
5. Compute total capacity of exoRNases and determine fraction of nucleotides
that can be diggested
6. Update pool of metabolites (H and H2O) created because of exonucleolytic
digestion
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaDegradation, self).initialize(sim, sim_data)

		rnaIds = sim_data.process.transcription.rnaData['id']

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		# Load RNase kinetic data
		endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
		exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
		self.KcatExoRNase = sim_data.constants.KcatExoRNase
		self.KcatEndoRNases = sim_data.process.rna_decay.kcats

		# Load RNase specificity
		self.TargetEndoRNasesFullMRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullMRNA
		self.TargetEndoRNasesFullTRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullTRNA
		self.TargetEndoRNasesFullRRNA_indexes = sim_data.process.rna_decay.TargetEndoRNasesFullRRNA

		self.mrna_index = sim_data.process.rna_decay.mrna_index
		self.trna_index = sim_data.process.rna_decay.trna_index
		self.rrna_index = sim_data.process.rna_decay.rrna_index
		self.rtrna_index = sim_data.process.rna_decay.rtrna_index

		# Load information about charged tRNA
		self.charged_trna_names = sim_data.process.transcription.charged_trna_names
		self.charged_trna = self.bulkMoleculesView(self.charged_trna_names)

		# Load first-order RNA degradation rates (estimated by mRNA half-life data)
		self.rnaDegRates = sim_data.process.transcription.rnaData['degRate']

		shuffleIdxs = None
		if hasattr(sim_data.process.transcription, "rnaDegRateShuffleIdxs") and sim_data.process.transcription.rnaDegRateShuffleIdxs != None:
			shuffleIdxs = sim_data.process.transcription.rnaDegRateShuffleIdxs
			self.rnaDegRates = self.rnaDegRates[shuffleIdxs]

		self.isMRna = sim_data.process.transcription.rnaData["isMRna"]
		self.isRRna = sim_data.process.transcription.rnaData["isRRna"]
		self.isTRna = sim_data.process.transcription.rnaData["isTRna"]

		self.rna_lengths = sim_data.process.transcription.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		endCleavageMetaboliteIds = [id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs]
		endCleavageMetaboliteIds.extend(["WATER[c]", "PPI[c]", "PROTON[c]"])
		nmpIdxs = range(4)
		h2oIdx = endCleavageMetaboliteIds.index("WATER[c]")
		ppiIdx = endCleavageMetaboliteIds.index("PPI[c]")
		hIdx = endCleavageMetaboliteIds.index("PROTON[c]")
		self.endoDegradationSMatrix = np.zeros((len(endCleavageMetaboliteIds), len(rnaIds)), np.int64)
		self.endoDegradationSMatrix[nmpIdxs, :] = units.transpose(sim_data.process.transcription.rnaData['countsACGU']).asNumber()
		self.endoDegradationSMatrix[h2oIdx, :] = 0
		self.endoDegradationSMatrix[ppiIdx, :] = 1
		self.endoDegradationSMatrix[hIdx, :] = 0

		# Build Views
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.h2o = self.bulkMoleculesView(["WATER[c]"])
		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.h = self.bulkMoleculesView(["PROTON[c]"])

		self.fragmentMetabolites = self.bulkMoleculesView(endCleavageMetaboliteIds)
		self.fragmentBases = self.bulkMoleculesView([id_ + "[c]" for id_ in sim_data.moleculeGroups.fragmentNT_IDs])

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.rrfaIdx = sim_data.process.transcription.rnaData["id"].tolist().index("RRFA-RRNA[c]")
		self.rrlaIdx = sim_data.process.transcription.rnaData["id"].tolist().index("RRLA-RRNA[c]")
		self.rrsaIdx = sim_data.process.transcription.rnaData["id"].tolist().index("RRSA-RRNA[c]")


		self.endoRnases = self.bulkMoleculesView(endoRnaseIds)
		self.exoRnases = self.bulkMoleculesView(exoRnaseIds)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		# Load Michaelis-Menten constants fitted to recapitulate first-order RNA decay model
		self.Km = sim_data.process.transcription.rnaData["KmEndoRNase"]

		self.EndoRNaseCoop = sim_data.constants.EndoRNaseCooperation
		self.EndoRNaseFunc = sim_data.constants.EndoRNaseFunction


	def calculateRequest(self):

		# Compute factor that convert counts into concentration, and vice versa
		cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
		cell_volume = cell_mass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cell_volume)

		# Get total counts of RNAs including rRNAs and charged tRNAs
		rna_counts = self.rnas.total().copy()
		rna_counts[self.rrsaIdx] += self.ribosome30S.total()
		rna_counts[[self.rrlaIdx, self.rrfaIdx]] += self.ribosome50S.total()
		rna_counts[[self.rrlaIdx, self.rrfaIdx, self.rrsaIdx]] += self.activeRibosomes.total()
		rna_counts[self.isTRna] += self.charged_trna.total()

		# Get counts of endoRNases
		endornase_counts = self.endoRnases.total().copy()
		total_kcat_endornase = units.dot(self.KcatEndoRNases, endornase_counts)

		# Calculate the fraction of active endoRNases for each RNA based on
		# Michaelis-Menten kinetics
		if self.EndoRNaseCoop:
			frac_endornase_saturated = (
				countsToMolar * rna_counts / self.Km / (
				1 + units.sum((countsToMolar * rna_counts) / self.Km))
			).asNumber()
		else:
			frac_endornase_saturated = (
				countsToMolar * rna_counts / (self.Km + (
				countsToMolar * rna_counts))
			).asNumber()

		# Calculate difference in degradation rates from first-order decay
		# and the number of EndoRNases per one molecule of RNA
		total_endornase_counts = np.sum(endornase_counts)
		diff_relative_first_order_decay = units.sum(
			units.abs(self.rnaDegRates * rna_counts -
				total_kcat_endornase * frac_endornase_saturated)
			)
		endornase_per_rna = total_endornase_counts / np.sum(rna_counts)

		self.writeToListener("RnaDegradationListener",
			"FractionActiveEndoRNases",
			np.sum(frac_endornase_saturated)
			)
		self.writeToListener("RnaDegradationListener",
			"DiffRelativeFirstOrderDecay",
			diff_relative_first_order_decay.asNumber()
			)
		self.writeToListener(
			"RnaDegradationListener",
			"FractEndoRRnaCounts",
			endornase_per_rna)

		if self.EndoRNaseFunc:
			# Dissect RNAse specificity into mRNA, tRNA, and rRNA
			mrna_specificity = np.dot(frac_endornase_saturated, self.isMRna)
			trna_specificity = np.dot(frac_endornase_saturated, self.isTRna)
			rrna_specificity = np.dot(frac_endornase_saturated, self.isRRna)
	
			n_total_mrnas_to_degrade = self._calculate_total_n_to_degrade(
				mrna_specificity,
				total_kcat_endornase
				)
			n_total_trnas_to_degrade = self._calculate_total_n_to_degrade(
				trna_specificity,
				total_kcat_endornase
				)
			n_total_rrnas_to_degrade = self._calculate_total_n_to_degrade(
				rrna_specificity,
				total_kcat_endornase
				)
	
			# Compute RNAse specificity
			rna_specificity = frac_endornase_saturated / np.sum(frac_endornase_saturated)
	
			# Boolean variable that tracks existence of each RNA
			rna_exists = rna_counts.astype(np.bool)

			# Compute degradation probabilities of each RNA: for mRNAs, this
			# is based on the specificity of each mRNA. For tRNAs and rRNAs,
			# this is distributed evenly.
			mrna_deg_probs = 1. / np.dot(rna_specificity, self.isMRna * rna_exists) * rna_specificity * self.isMRna * rna_exists
			trna_deg_probs = 1. / np.sum(self.isTRna * rna_exists) * self.isTRna * rna_exists
			rrna_deg_probs = 1. / np.sum(self.isRRna * rna_exists) * self.isRRna * rna_exists

			# Mask RNA counts into each class of RNAs
			mrna_counts = rna_counts * self.isMRna
			trna_counts = rna_counts * self.isTRna
			rrna_counts = rna_counts * self.isRRna
	
			# Determine number of individual RNAs to be degraded for each class
			# of RNA.
			n_mrnas_to_degrade = self._get_rnas_to_degrade(
				n_total_mrnas_to_degrade, mrna_deg_probs, mrna_counts
				)

			n_trnas_to_degrade = self._get_rnas_to_degrade(
				n_total_trnas_to_degrade, trna_deg_probs, trna_counts
				)

			n_rrnas_to_degrade = self._get_rnas_to_degrade(
				n_total_rrnas_to_degrade, rrna_deg_probs, rrna_counts
				)
	
			n_rnas_to_degrade = n_mrnas_to_degrade + n_trnas_to_degrade + n_rrnas_to_degrade

		# First order decay with non-functional EndoRNase activity 
		# Determine mRNAs to be degraded by sampling a Poisson distribution
		# (Kdeg * RNA)
		else:
			n_rnas_to_degrade = np.fmin(
				self.randomState.poisson(
					(self.rnaDegRates * rna_counts).asNumber()
					),
				self.rnas.total()
				)

		self.rnas.requestIs(n_rnas_to_degrade)
		self.endoRnases.requestAll()
		self.exoRnases.requestAll()
		self.fragmentBases.requestAll()

		# Calculate the amount of water required for total RNA hydrolysis by
		# endo and exonucleases. Assuming complete hydrolysis for now. Note
		# that one additional water molecule is needed for each RNA to
		# hydrolyze the 5' diphosphate.
		waterForNewRnas = np.dot(n_rnas_to_degrade, self.rna_lengths)
		waterForLeftOverFragments = self.fragmentBases.total().sum()
		self.h2o.requestIs(waterForNewRnas + waterForLeftOverFragments)
		

	def evolveState(self):
		n_degraded_rna = self.rnas.counts()

		self.writeToListener(
			"RnaDegradationListener", "countRnaDegraded", n_degraded_rna
			)
		self.writeToListener(
			"RnaDegradationListener", "nucleotidesFromDegradation",
			np.dot(n_degraded_rna, self.rna_lengths)
			)

		# Calculate endolytic cleavage events

		# Modeling assumption: Once a RNA is cleaved by an endonuclease its
		# resulting nucleotides are lumped together as "polymerized fragments".
		# These fragments can carry over from previous timesteps. We are also
		# assuming that during endonucleolytic cleavage the 5'terminal
		# phosphate is removed. This is modeled as all of the fragments being
		# one long linear chain of "fragment bases".

		# Example:
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-OH
		# ==>
		# Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + PPi
		# Note: Lack of -OH on 3' end of chain

		metabolitesEndoCleavage = np.dot(
			self.endoDegradationSMatrix, n_degraded_rna)
		self.rnas.countsIs(0)
		self.fragmentMetabolites.countsInc(metabolitesEndoCleavage)

		# Check if exonucleolytic digestion can happen 
		if self.fragmentBases.counts().sum() == 0:
			return

		# Calculate exolytic cleavage events

		# Modeling assumption: We model fragments as one long fragment chain of
		# polymerized nucleotides. We are also assuming that there is no
		# sequence specificity or bias towards which nucleotides are
		# hydrolyzed.

		# Example:
		# Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + 3 H2O
		# ==>
		# 3 NMP + 3 H(+)
		# Note: Lack of -OH on 3' end of chain

		n_exornases = self.exoRnases.counts()
		n_fragment_bases = self.fragmentBases.counts()
		n_fragment_bases_sum = n_fragment_bases.sum()

		exornase_capacity = n_exornases.sum() * self.KcatExoRNase * (
				units.s * self.timeStepSec()
			)

		if exornase_capacity >= n_fragment_bases_sum:
			self.nmps.countsInc(n_fragment_bases)
			self.h2o.countsDec(n_fragment_bases_sum)
			self.h.countsInc(n_fragment_bases_sum)
			self.fragmentBases.countsIs(0)

			total_fragment_bases_digested = n_fragment_bases_sum

		else:
			fragmentSpecificity = n_fragment_bases / n_fragment_bases_sum
			possibleBasesToDigest = self.randomState.multinomial(
				exornase_capacity, fragmentSpecificity)
			n_fragment_bases_digested = n_fragment_bases - np.fmax(
				n_fragment_bases - possibleBasesToDigest, 0)

			total_fragment_bases_digested = n_fragment_bases_digested.sum()

			self.nmps.countsInc(n_fragment_bases_digested)
			self.h2o.countsDec(total_fragment_bases_digested)
			self.h.countsInc(total_fragment_bases_digested)
			self.fragmentBases.countsDec(n_fragment_bases_digested)

		self.writeToListener("RnaDegradationListener",
			"fragmentBasesDigested", total_fragment_bases_digested)


	def _calculate_total_n_to_degrade(self, specificity, total_kcat_endornase):
		"""
		Calculate the total number of RNAs to degrade for a specific class of
		RNAs, based on the specificity of endoRNases on that specific class and
		the total kcat value of the endoRNases.

		Args:
			specificity: Sum of fraction of active endoRNases for all RNAs
			in a given class
			total_kcat_endornase: The summed kcat of all existing endoRNases
		Returns:
			Total number of RNAs to degrade for the given class of RNAs
		"""
		return np.round(
			(specificity * total_kcat_endornase
			 * (units.s * self.timeStepSec())).asNumber()
			)


	def _get_rnas_to_degrade(self, n_total_rnas_to_degrade, rna_deg_probs,
			rna_counts):
		"""
		Distributes the total count of RNAs to degrade for each class of RNAs
		into individual RNAs, based on the given degradation probabilities
		of individual RNAs. The upper bound is set by the current count of the
		specific RNA.

		Args:
			n_total_rnas_to_degrade: Total number of RNAs to degrade for the
			given class of RNAs (integer, scalar)
			rna_deg_probs: Degradation probabilities of each RNA (vector of
			equal length to the total number of different RNAs)
			rna_counts: Current counts of each RNA molecule (vector of equal
			length to the total number of different RNAs)
		Returns:
			Vector of equal length to rna_counts, specifying the number of
			molecules to degrade for each RNA
		"""
		n_rnas_to_degrade = np.zeros_like(rna_counts)

		if rna_counts.sum() != 0:
			while n_rnas_to_degrade.sum() < n_total_rnas_to_degrade:
				n_rnas_to_degrade += np.fmin(
					self.randomState.multinomial(
						n_total_rnas_to_degrade - n_rnas_to_degrade.sum(),
						rna_deg_probs
						),
					rna_counts
					)

		return n_rnas_to_degrade