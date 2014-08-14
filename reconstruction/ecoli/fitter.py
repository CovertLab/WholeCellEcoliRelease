#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

TODO:
- document the math
- compute and use activation rates for RNA poly, ribosomes
- fit metabolism enzyme expression
- replace fake metabolite pools with measured metabolite pools

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

from __future__ import division

import numpy as np
import os
import copy
import collections

import wholecell.states.bulk_molecules
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from reconstruction.units.unit_registration import UREG
from reconstruction.units.unit_registration import Q_
import pint
pint._DEFAULT_REGISTRY = UREG

# Constants (should be moved to KB)
RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA
GROWTH_ASSOCIATED_MAINTENANCE = 59.81 # mmol/gDCW (from Feist)
NON_GROWTH_ASSOCIATED_MAINTENANCE = 8.39 # mmol/gDCW/hr (from Feist)

# Correction factors
EXCESS_RNAP_CAPACITY = 1.
EXCESS_FREE_DNTP_CAPACITY = 1.3
# If RNA-poly capacity exactly matches the amount needed to double RNAs over a 
# cell cycle, the simulation will be unable to double RNAs since a small number
# of RNA-polymerases must be turned over following termination.  It may be 
# possible to choose an excess capacity coefficient rationally based on 
# diffusive limitations, i.e., one that does not depend on simulation 
# particulars, but this has yet to be explored.

# Fitter logic 
# TODO: confirm this with Derek
# TODO: split off these subroutines in the main fitter function
# 1) Assign expected quantities based on dry mass composition, expression, and sequences
# 2) Ensure that there is enough RNAP/ribosome capacity for (1), and adjust if needed
# 3) Update the metabolism FBA objective based on expression

def fitAtLevel(fitLevel, kb, simOutDir):
	# TODO: Obviously make this more sophisticated
	if fitLevel == 1:
		fitKb(kb)

	if fitLevel == 2:
		print simOutDir


def fitKb(kb):

	# Construct bulk container

	bulkContainer = wholecell.states.bulk_molecules.bulkObjectsContainer(kb, dtype = np.dtype("float64"))

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	mRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	miscRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMiscRna"]])
	rRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna"]])
	tRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isTRna"]])

	rRna23SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna23S"]])
	rRna16SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna16S"]])
	rRna5SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna5S"]])

	adjustDryCompositionBasedOnChromosomeSeq(bulkContainer, kb)
	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60]

	### RNA Mass fraction ###
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit.to('DCW_g') * rnaMassFraction
	setRNACounts(
		kb, rnaMass, mRnaView,
		rRna23SView, rRna16SView, rRna5SView, tRnaView
		)


	### Protein Mass fraction ###

	monomersView = bulkContainer.countsView(kb.monomerData["id"])

	monomerMassFraction = float(dryComposition60min["proteinMassFraction"])
	monomerMass = kb.avgCellDryMassInit * monomerMassFraction
	setMonomerCounts(kb, monomerMass, monomersView)

	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapView = bulkContainer.countsView(kb.rnapIds)

	## Number of ribosomes needed ##
	monomerLengths = np.sum(kb.monomerData['aaCounts'], axis = 1)
	nRibosomesNeeded = np.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]
			) * monomersView.counts()
		).to('dimensionless').magnitude
	
	if np.sum(rRna23SView.counts()) < nRibosomesNeeded:
		raise NotImplementedError, "Cannot handle having too few ribosomes"

	## Number of RNA Polymerases ##
	rnaLengths = np.sum(kb.rnaData['countsACGU'], axis = 1)

	nRnapsNeeded = np.sum(
		rnaLengths / kb.rnaPolymeraseElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
			) * rnaView.counts()
		).to('dimensionless').magnitude * EXCESS_RNAP_CAPACITY

	minRnapCounts = (
		nRnapsNeeded * np.array([2, 1, 1, 1]) # Subunit stoichiometry
		)

	rnapView.countsIs(
		np.fmax(rnapView.counts(), minRnapCounts)
		)

	
	### Modify kbFit to reflect our bulk container ###

	## Fraction of active Ribosomes ##
	kb.fracActiveRibosomes = Q_(float(nRibosomesNeeded) / np.sum(rRna23SView.counts()), "dimensionless")

	## RNA and monomer expression ##
	rnaExpressionContainer = wholecell.containers.bulk_objects_container.BulkObjectsContainer(list(kb.rnaData["id"]), dtype = np.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(rnaView.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		kb.monomerData["rnaId"][kb.monomerIndexToRnaMapping] == kb.rnaData["id"][kb.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids"

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * normalize(monomersView.counts()[kb.monomerIndexToRnaMapping])
		)

	kb.rnaExpression['expression'] = Q_(rnaExpressionContainer.counts(),'dimensionless')

	# Set number of RNAs based on expression we just set
	nRnas = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude,
		kb.rnaData["mw"].to('g/mol').magnitude,
		kb.rnaExpression['expression'].to('dimensionless').magnitude,
		kb.nAvogadro.to('1/mol').magnitude
		)

	rnaView.countsIs(nRnas * kb.rnaExpression['expression'].to('dimensionless').magnitude)

	## Synthesis probabilities ##
	synthProb = normalize(
			(
			Q_(1, 'second') * (
				np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
				) * rnaView.counts()
			).to('dimensionless').magnitude
		)

	kb.rnaData["synthProb"][:] = synthProb


	## Calculate and set maintenance values

	# ----- Non growth associated maintenance -----
	kb.NGAM = Q_(NON_GROWTH_ASSOCIATED_MAINTENANCE, "mmol/DCW_g/hr")

	# ----- Growth associated maintenance -----

	# GTPs used for translation (recycled, not incorporated into biomass)

	aaMmolPerGDCW = (
			np.sum(
				kb.monomerData["aaCounts"] *
				np.tile(monomersView.counts().reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / kb.nAvogadro.to('amino_acid/mmol')) *
				(1 / kb.avgCellDryMassInit)
			)
		).to('mmol/DCW_g')

	aasUsedOverCellCycle = aaMmolPerGDCW.magnitude.sum()
	gtpUsedOverCellCycleMmolPerGDCW = kb.gtpPerTranslation * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		GROWTH_ASSOCIATED_MAINTENANCE -
		gtpUsedOverCellCycleMmolPerGDCW
		)

	# Assign the growth associated "dark energy" to translation
	# TODO: Distribute it amongst growth-related processes
	kb.gtpPerTranslation += darkATP / aaMmolPerGDCW.magnitude.sum()

def normalize(array):
	return np.array(array).astype("float") / np.linalg.norm(array, 1)

def countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro):
	assert np.allclose(np.sum(relativeExpression), 1)
	assert type(mass) != Q_
	assert type(mws) != Q_
	assert type(relativeExpression) != Q_
	assert type(nAvogadro) != Q_
	return mass / np.dot(mws / nAvogadro, relativeExpression)

def setRNACounts(kb, rnaMass, mRnaView, rRna23SView, rRna16SView, rRna5SView, tRnaView):

	## 23S rRNA Mass Fractions ##

	# Assume all 23S rRNAs are expressed equally
	rRna23SExpression = normalize(np.ones(rRna23SView.counts().size))

	nRRna23Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * RRNA23S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna23S"]].to('g/mol').magnitude,
		rRna23SExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	## 16S rRNA Mass Fractions ##

	# Assume all 16S rRNAs are expressed equally
	rRna16SExpression = normalize(np.ones(rRna16SView.counts().size))

	nRRna16Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * RRNA16S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna16S"]].to('g/mol').magnitude,
		rRna16SExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	## 5S rRNA Mass Fractions ##

	# Assume all 5S rRNAs are expressed equally
	rRna5SExpression = normalize(np.ones(rRna5SView.counts().size))

	nRRna5Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * RRNA5S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna5S"]].to('g/mol').magnitude,
		rRna5SExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	# ## Correct numbers of 23S, 16S, 5S rRNAs so that they are all equal
	# # TODO: Maybe don't need to do this at some point (i.e., when the model is more sophisticated)
	nRRna23Ss = nRRna16Ss = nRRna5Ss = np.mean((nRRna23Ss, nRRna16Ss, nRRna5Ss))

	# TODO: Remove this hack once complexation is working
	rRna23SExpression[:] = 0.
	rRna23SExpression[0] = 1.

	rRna16SExpression[:] = 0.
	rRna16SExpression[0] = 1.

	rRna5SExpression[:] = 0.
	rRna5SExpression[0] = 1.

	rRna23SView.countsIs((nRRna23Ss * rRna23SExpression))
	rRna16SView.countsIs((nRRna16Ss * rRna16SExpression))
	rRna5SView.countsIs((nRRna5Ss * rRna5SExpression))

	## tRNA Mass Fractions ##

	# Assume all tRNAs are expressed equally (TODO: Change this based on monomer expression!)
	tRnaExpression = normalize(np.ones(tRnaView.counts().size))

	nTRnas = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * TRNA_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isTRna"]].to('g/mol').magnitude,
		tRnaExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	tRnaView.countsIs((nTRnas * tRnaExpression))

	## mRNA Mass Fractions ##

	mRnaExpression = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])

	nMRnas = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * MRNA_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isMRna"]].to('g/mol').magnitude,
		mRnaExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	mRnaView.countsIs((nMRnas * mRnaExpression))

def setMonomerCounts(kb, monomerMass, monomersView):

	monomerExpression = normalize(
		kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping].magnitude /
		(np.log(2) / kb.cellCycleLen.to("s").magnitude + kb.monomerData["degRate"].to("1/s").magnitude)
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.to("DCW_g").magnitude,
		kb.monomerData["mw"].to('g/mol').magnitude,
		monomerExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	monomersView.countsIs((nMonomers * monomerExpression))

def calcChromosomeMass(numA, numC, numG, numT, kb):
	weights = collections.OrderedDict({
		"A": (
			float(kb.getMass(["DAMP[n]"]).magnitude)
			),
		"C": (
			float(kb.getMass(["DCMP[n]"]).magnitude)
			),
		"G": (
			float(kb.getMass(["DGMP[n]"]).magnitude)
			),
		"T": (
			float(kb.getMass(["DTMP[n]"]).magnitude)
			),
		})

	seqLen = numA + numC + numG + numT

	return (
		weights["A"] * numA +
		weights["C"] * numC +
		weights["G"] * numG +
		weights["T"] * numT# -
		# TODO: Ask Nick about this line (below)
		#seqLen * 17.01 # Note: no factor of 2 is needed because the num variables account for double-strandedness
		)


def calcNumDntpsDnmps(kb, tau_d):
	if tau_d != 60:
		raise NotImplementedError, "This function currently only works for the special case of 60 min doubling time."

	nPolymerases = 4
	k_elng = kb.dnaPolymeraseElongationRate.to("nucleotide / s").magnitude

	seqLen = len(kb.genomeSeq)
	t_C = seqLen / 2. / k_elng # Length of C period (approximate)
	tau_d = kb.cellCycleLen.to("s").magnitude # Doubling time
	N_p = 2 * seqLen	# Number of polymerized dNMPs (DNA is double-stranded, thus the factor of 2)
	dt = kb.timeStep.to("s").magnitude

	return np.fmax(
		2 * N_p / np.exp((np.log(2) / tau_d) * t_C),
		(nPolymerases * k_elng * dt) / (np.exp((np.log(2) / tau_d) * dt) - 1)
		)


def adjustDryCompositionBasedOnChromosomeSeq(bulkContainer, kb):

	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60]
	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit * dnaMassFraction
	chromMass = calcChromosomeMass(
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A"),
		kb) / kb.nAvogadro.magnitude

	dnaMassCalc = chromMass

	fracDifference = (dnaMass.magnitude - dnaMassCalc) / kb.avgCellDryMassInit.magnitude
	# if fracDifference < 0:
	# 	raise NotImplementedError, "Have to add DNA mass. Make sure you want to do this."
	idx60Min = np.where(kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60)
	dNtpCompositionIdx = 3 # TODO: Get this from code somehow
	nElems = 9 # TODO: Get this from code somehow
	nonDNtpsIdxs = [x for x in range(1, nElems + 1) if x != dNtpCompositionIdx]
	amountToAdd = fracDifference / len(nonDNtpsIdxs)
	kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, nonDNtpsIdxs] += amountToAdd
	kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, dNtpCompositionIdx] = dnaMassCalc / kb.avgCellDryMassInit.magnitude
	assert np.allclose(1, kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, 1:].sum()), "Composition fractions must sum to 1!"


def setMetaboliteCountsFromBiomassFraction(kb, biomassContainer, fractionMetaboliteIds, fractionOfDryMass, fractionComposition):
	massFractionView = biomassContainer.countsView(
		list(fractionMetaboliteIds)
		)

	fractionOfDryMass = float(fractionOfDryMass)
	mass = kb.avgCellDryMassInit * fractionOfDryMass

	mws = kb.getMass(fractionMetaboliteIds)

	fractionMmolPerGDCW = (
			(
			mass * fractionComposition
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	massFractionView.countsInc(
		fractionMmolPerGDCW.magnitude
		)

if __name__ == "__main__":
	import wholecell.utils.constants

	kb = cPickle.load(
		open(os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			), "rb")
		)
	
	fitKb(kb)
