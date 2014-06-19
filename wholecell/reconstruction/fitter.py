#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

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

from units.unit_registration import UREG
from units.unit_registration import Q_
import pint
pint._DEFAULT_REGISTRY = UREG

# Constants (should be moved to KB)
RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA

# Correction factors
EXCESS_RNAP_CAPACITY = 2
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

	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60]

	adjustCompositionBasedOnChromosomeSeq(bulkContainer, kb)
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

	### DNA Mass fraction ###
	dNtpIds = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	dNmpIds = ["DAMP[c]", "DCMP[c]", "DGMP[c]", "DTMP[c]"]
	dNtpsView = bulkContainer.countsView(dNtpIds)

	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit * dnaMassFraction

	dNtpRelativeAmounts = normalize(np.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		]))

	dNtpIdxs = [np.where(kb.bulkMolecules["moleculeId"] == idx)[0][0] for idx in dNtpIds]
	dNmpIdxs = [np.where(kb.bulkMolecules["moleculeId"] == idx)[0][0] for idx in dNmpIds]


	dNtpMws = kb.bulkMolecules["mass"][dNtpIdxs].sum(axis = 1)
	dNmpMws = kb.bulkMolecules["mass"][dNmpIdxs].sum(axis = 1)

	nDNtps = countsFromMassAndExpression(
		dnaMass.to('DCW_g').magnitude,
		dNmpMws.to('g/mol').magnitude - 17.01,
		dNtpRelativeAmounts,
		kb.nAvogadro.to('1/mol').magnitude
		)

	dNtpsView.countsIs(nDNtps * dNtpRelativeAmounts)

	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]
	rnapView = bulkContainer.countsView(rnapIds)

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
	kb.parameters["fracActiveRibosomes"] = Q_(float(nRibosomesNeeded) / np.sum(rRna23SView.counts()), "dimensionless")
	kb.__dict__.update(kb.parameters)

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


	## Full WT Biomass function ##

	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)

	# Amino acid fraction
	aminoAcidView = biomassContainer.countsView(kb.aaIDs)

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
	
	aminoAcidView.countsIs(
		aaMmolPerGDCW.magnitude
		)

	# RNA fraction
	ntpView = biomassContainer.countsView(
		["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
		)

	ntpMmolPerGDCW = (
			np.sum(
				kb.rnaData["countsACGU"] *
				np.tile(rnaView.counts().reshape(-1, 1), (1, 4)),
				axis = 0
			) * (
				(1 / kb.nAvogadro.to('nucleotide / mmol')) *
				(1 / kb.avgCellDryMassInit)
			)
		).to('mmol/DCW_g')

	ntpView.countsIs(
		ntpMmolPerGDCW.magnitude
		)

	# DNA fraction

	dNtpView = biomassContainer.countsView(		# TODO: Better name so as not to confuse with bulkContainer view
		["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
		)

	dNtpMmolPerGDCW = (
			Q_(dNtpsView.counts(),'nucleotide') * (
			(1 / kb.nAvogadro.to('nucleotide/mmol')) *
			(1 / kb.avgCellDryMassInit)
			)
		).to('mmol/DCW_g')

	dNtpView.countsIs(
		dNtpMmolPerGDCW.magnitude
		)
	
	# Glycogen fraction

	glycogenView = biomassContainer.countsView(
		list(kb.cellGlycogenFractionData["metaboliteId"])
		)

	glycogenMassFraction = float(dryComposition60min["glycogenMassFraction"])
	glycogenMass = kb.avgCellDryMassInit * glycogenMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellGlycogenFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB

	glycogenMmolPerGDCW = (
			(
				glycogenMass * kb.cellGlycogenFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
				1 / kb.avgCellDryMassInit
			)
		).to('mmol/DCW_g')
	
	glycogenView.countsIs(
		glycogenMmolPerGDCW.magnitude
		)

	# Murein fraction

	mureinView = biomassContainer.countsView(
		list(kb.cellMureinFractionData["metaboliteId"])
		)

	mureinMassFraction = float(dryComposition60min["mureinMassFraction"])
	mureinMass = kb.avgCellDryMassInit * mureinMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellMureinFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB

	mureinMmolPerGDCW = (
			(
			mureinMass * kb.cellMureinFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	mureinView.countsIs(
		mureinMmolPerGDCW.magnitude
		)

	# LPS fraction

	lpsView = biomassContainer.countsView(
		list(kb.cellLPSFractionData["metaboliteId"])
		)

	lpsMassFraction = float(dryComposition60min["lpsMassFraction"])
	lpsMass = kb.avgCellDryMassInit * lpsMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLPSFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB

	lpsMmolPerGDCW = (
			(
			lpsMass * kb.cellLPSFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	lpsView.countsIs(
		lpsMmolPerGDCW.magnitude
		)

	# Lipid fraction

	lipidView = biomassContainer.countsView(
		list(kb.cellLipidFractionData["metaboliteId"])
		)

	lipidMassFraction = float(dryComposition60min["lipidMassFraction"])
	lipidMass = kb.avgCellDryMassInit * lipidMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLipidFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB

	lipidMmolPerGDCW = (
			(
			lipidMass * kb.cellLipidFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	lipidView.countsIs(
		lipidMmolPerGDCW.magnitude
		)

	# Inorganic ion fraction

	inorganicIonView = biomassContainer.countsView(
		list(kb.cellInorganicIonFractionData["metaboliteId"])
		)

	inorganicIonMassFraction = float(dryComposition60min["inorganicIonMassFraction"])
	inorganicIonMass = kb.avgCellDryMassInit * inorganicIonMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellInorganicIonFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB
	
	inorganicIonMmolPerGDCW = (
			(
			inorganicIonMass * kb.cellInorganicIonFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	inorganicIonView.countsIs(
		inorganicIonMmolPerGDCW.magnitude
		)

	# Soluble pool fraction

	solublePoolView = biomassContainer.countsView(
		list(kb.cellSolublePoolFractionData["metaboliteId"])
		)

	solublePoolMassFraction = float(dryComposition60min["solublePoolMassFraction"])
	solublePoolMass = kb.avgCellDryMassInit * solublePoolMassFraction

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellSolublePoolFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1) # TOKB
	
	solublePoolMmolPerGDCW = (
			(
			solublePoolMass * kb.cellSolublePoolFractionData["massFraction"]
			) / mws.to('DCW_g/mmol') * (
			1 / kb.avgCellDryMassInit)
		).to('mmol/DCW_g')

	solublePoolView.countsIs(
		solublePoolMmolPerGDCW.magnitude
		)
	
	# Validate

	bulkMoleculesIdxs = np.array([
		np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in biomassContainer._objectNames
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].sum(1).magnitude # TOKB
	# TODO

	# aaIdxs = np.array([
	# 	np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in aaIDs
	# 	])
	# aaMws = kb.bulkMolecules["mass"][aaIdxs].magnitude

	# import ipdb; ipdb.set_trace()

	# TODO: Unhack this
	kb.wildtypeBiomass.struct_array["biomassFlux"] = biomassContainer.counts()

	import ipdb; ipdb.set_trace()
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
		# Handles reverse complement
		"A": (
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DAMP[c]"]["mass"].sum().magnitude) +
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DTMP[c]"]["mass"].sum().magnitude)
			),
		"C": (
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DCMP[c]"]["mass"].sum().magnitude) +
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DGMP[c]"]["mass"].sum().magnitude)
			),
		"G": (
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DGMP[c]"]["mass"].sum().magnitude) +
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DCMP[c]"]["mass"].sum().magnitude)
			),
		"T": (
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DTMP[c]"]["mass"].sum().magnitude) +
			float(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == "DAMP[c]"]["mass"].sum().magnitude)
			),
		})

	seqLen = numA + numC + numG + numT

	return (
		weights["A"] * numA +
		weights["C"] * numC +
		weights["G"] * numG +
		weights["T"] * numT -
		2 * seqLen * 17.01 # The "2" is because DNA is double stranded (need to account for both)
		)


def adjustCompositionBasedOnChromosomeSeq(bulkContainer, kb):

	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60]
	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit * dnaMassFraction
	dnaMassCalc = calcChromosomeMass(
		kb.genomeSeq.count("A"),
		kb.genomeSeq.count("C"),
		kb.genomeSeq.count("G"),
		kb.genomeSeq.count("T"),
		kb) / kb.nAvogadro.magnitude
	calcNumDntps(kb, 60)
	fracDifference = (dnaMass.magnitude - dnaMassCalc) / kb.avgCellDryMassInit.magnitude
	if fracDifference < 0:
		raise NotImplementedError, "Have to add DNA mass. Make sure you want to do this."
	idx60Min = np.where(kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60)
	dNtpCompositionIdx = 3 # TODO: Get this from code somehow
	nElems = 9 # TODO: Get this from code somehow
	nonDNtpsIdxs = [x for x in range(1, nElems + 1) if x != dNtpCompositionIdx]
	amountToAdd = fracDifference / len(nonDNtpsIdxs)
	kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, nonDNtpsIdxs] += amountToAdd
	kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, dNtpCompositionIdx] = dnaMassCalc / kb.avgCellDryMassInit.magnitude
	assert np.allclose(1, kb.cellDryMassComposition.struct_array.view((np.float, 10))[idx60Min, 1:].sum()), "Composition fractions must sum to 1!"

if __name__ == "__main__":
	import wholecell.utils.constants

	kb = cPickle.load(
		open(os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			), "rb")
		)
	
	fitKb(kb)
