#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

from __future__ import division

import numpy
import os
import copy

import wholecell.states.bulk_molecules
import wholecell.utils.rand_stream
from units.unit_registration import Q_

def fitKb(kb):

	# Construct bulk container

	bulkContainer = wholecell.states.bulk_molecules.bulkObjectsContainer(kb, dtype = numpy.dtype("float64"))

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	mRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	miscRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMiscRna"]])
	rRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna"]])
	tRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isTRna"]])

	rRna23SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna23S"]])
	rRna16SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna16S"]])
	rRna5SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna5S"]])

	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].magnitude == 60]

	### RNA Mass Fractions ###
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit.magnitude * rnaMassFraction

	## 23S rRNA Mass Fractions ##
	rRna23SMassFraction = 0.525 # TOKB # This is the fraction of RNA that is 23S rRNA

	# Assume all 23S rRNAs are expressed equally
	rRna23SExpression = normalize(numpy.ones(rRna23SView.counts().size))

	nRRna23Ss = countsFromMassAndExpression(
		rnaMass * rRna23SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna23S"]],
		rRna23SExpression,
		kb.nAvogadro.magnitude
		)

	## 16S rRNA Mass Fractions ##
	rRna16SMassFraction = 0.271 # TOKB # This is the fraction of RNA that is 16S rRNA

	# Assume all 16S rRNAs are expressed equally
	rRna16SExpression = normalize(numpy.ones(rRna16SView.counts().size))

	nRRna16Ss = countsFromMassAndExpression(
		rnaMass * rRna16SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna16S"]],
		rRna16SExpression,
		kb.nAvogadro.magnitude
		)

	## 5S rRNA Mass Fractions ##
	rRna5SMassFraction = 0.017 # TOKB # This is the fraction of RNA that is 5S rRNA

	# Assume all 5S rRNAs are expressed equally
	rRna5SExpression = normalize(numpy.ones(rRna5SView.counts().size))

	nRRna5Ss = countsFromMassAndExpression(
		rnaMass * rRna5SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna5S"]],
		rRna5SExpression,
		kb.nAvogadro.magnitude
		)

	# ## Correct numbers of 23S, 16S, 5S rRNAs so that they are all equal
	# # TODO: Maybe don't need to do this at some point (i.e., when the model is more sophisticated)
	nRRna23Ss = nRRna16Ss = nRRna5Ss = numpy.mean((nRRna23Ss, nRRna16Ss, nRRna5Ss))

	rRna23SView.countsIs((nRRna23Ss * rRna23SExpression))
	rRna16SView.countsIs((nRRna16Ss * rRna16SExpression))
	rRna5SView.countsIs((nRRna5Ss * rRna5SExpression))

	## tRNA Mass Fractions ##
	tRnaMassFraction = 0.146 # TOKB # This is the fraction of RNA that is tRNA

	# Assume all tRNAs are expressed equally (TODO: Change this based on monomer expression!)
	tRnaExpression = normalize(numpy.ones(tRnaView.counts().size))

	nTRnas = countsFromMassAndExpression(
		rnaMass * tRnaMassFraction,
		kb.rnaData["mw"][kb.rnaData["isTRna"]],
		tRnaExpression,
		kb.nAvogadro.magnitude
		)

	tRnaView.countsIs((nTRnas * tRnaExpression))

	## mRNA Mass Fractions ##
	mRnaMassFraction = 0.041 # TOKB # This is the fraction of RNA that is mRNA

	mRnaExpression = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])

	nMRnas = countsFromMassAndExpression(
		rnaMass * mRnaMassFraction,
		kb.rnaData["mw"][kb.rnaData["isMRna"]],
		mRnaExpression,
		kb.nAvogadro.magnitude
		)

	mRnaView.countsIs((nMRnas * mRnaExpression))


	### Protein Mass fraction ###

	monomersView = bulkContainer.countsView(kb.monomerData["id"])

	monomerMassFraction = float(dryComposition60min["proteinMassFraction"])
	monomerMass = kb.avgCellDryMassInit.magnitude * monomerMassFraction

	monomerExpression = normalize(kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping])

	nMonomers = countsFromMassAndExpression(
		monomerMass,
		kb.monomerData["mw"],
		monomerExpression,
		kb.nAvogadro.magnitude
		)

	monomersView.countsIs((nMonomers * monomerExpression))


	### DNA Mass fraction ###
	# TODO: Don't just keep dNTPs in the soluble pool
	# TODO: Compute based on chromosome sequence
	dNtpsView = bulkContainer.countsView(["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"])

	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit.magnitude * dnaMassFraction

	dNtpRelativeAmounts = normalize(numpy.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		]))

	dNtpMws = numpy.array([
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "DATP[c]"].magnitude[0],
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "DCTP[c]"].magnitude[0],
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "DGTP[c]"].magnitude[0],
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "DTTP[c]"].magnitude[0]
		])

	nDNtps = countsFromMassAndExpression(
		dnaMass,
		dNtpMws,
		dNtpRelativeAmounts,
		kb.nAvogadro.magnitude
		)

	dNtpsView.countsIs(nDNtps * dNtpRelativeAmounts)



	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]
	rnapView = bulkContainer.countsView(rnapIds)

	## Number of ribosomes needed ##
	monomerLengths = numpy.sum(kb.monomerData['aaCounts'], axis = 1)
	nRibosomesNeeded = numpy.sum(
		monomerLengths / kb.ribosomeElongationRate.magnitude * (
			numpy.log(2) / kb.cellCycleLen.magnitude
			) * monomersView.counts()
		)

	if numpy.sum(rRna23SView.counts()) < nRibosomesNeeded:
		raise NotImplementedError, "Cannot handle having too few ribosomes"

	## Number of RNA Polymerases ##
	rnaLengths = numpy.sum(kb.rnaData['countsAUCG'], axis = 1)
	nRnapsNeeded = numpy.sum(
		rnaLengths / kb.rnaPolymeraseElongationRate.magnitude * (
			numpy.log(2) / kb.cellCycleLen.magnitude + kb.rnaData["degRate"]
			) * rnaView.counts()
		)

	minRnapCounts = (
		nRnapsNeeded * numpy.array([2, 1, 1, 1]) # Subunit stoichiometry
		)

	rnapView.countsIs(
		numpy.fmax(rnapView.counts(), minRnapCounts)
		)


	### Modify kbFit to reflect our bulk container ###

	## Fraction of active Ribosomes ##
	kb.parameters["fracActiveRibosomes"] = Q_(float(nRibosomesNeeded) / numpy.sum(rRna23SView.counts()), "dimensionless")
	# kb.parameters["fracActiveRibosomes"] = Q_(0.7, "dimensionless")
	kb.__dict__.update(kb.parameters)

	## RNA and monomer expression ##
	rnaExpressionContainer = wholecell.containers.bulk_objects_container.BulkObjectsContainer(list(kb.rnaData["id"]), dtype = numpy.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(rnaView.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert numpy.all(
		kb.monomerData["rnaId"][kb.monomerIndexToRnaMapping] == kb.rnaData["id"][kb.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids"

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	mRnaExpressionFrac = numpy.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * normalize(monomersView.counts()[kb.monomerIndexToRnaMapping])
		)

	# TODO: Remove this hack! Keep track of units.
	kb.rnaExpression.struct_array['expression'] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = countsFromMassAndExpression(
		rnaMass,
		kb.rnaData["mw"],
		kb.rnaExpression['expression'].magnitude,
		kb.nAvogadro.magnitude
		)

	rnaView.countsIs(nRnas * kb.rnaExpression['expression'].magnitude)


	## Synthesis probabilities ##

	synthProb = normalize(
		rnaLengths / kb.rnaPolymeraseElongationRate.magnitude * (
			numpy.log(2) / kb.cellCycleLen.magnitude + kb.rnaData["degRate"]
			) * rnaView.counts()
		)

	kb.rnaData["synthProb"][:] = synthProb


	## Full WT Biomass function ##

	biomassContainer = wholecell.containers.bulk_objects_container.BulkObjectsContainer(list(kb.wildtypeBiomass["metaboliteId"]), dtype = numpy.dtype("float64"))

	# Amino acid fraction
	oneToThreeMapping = dict((
		("A", "ALA-L[c]"), ("R", "ARG-L[c]"), ("N", "ASN-L[c]"), ("D", "ASP-L[c]"),
		("C", "CYS-L[c]"), ("E", "GLU-L[c]"), ("Q", "GLN-L[c]"), ("G", "GLY[c]"),
		("H", "HIS-L[c]"), ("I", "ILE-L[c]"), ("L", "LEU-L[c]"), ("K", "LYS-L[c]"),
		("M", "MET-L[c]"), ("F", "PHE-L[c]"), ("P", "PRO-L[c]"), ("S", "SER-L[c]"),
		("T", "THR-L[c]"), ("W", "TRP-L[c]"), ("Y", "TYR-L[c]"), ("U", "SEC-L[c]"),
		("V", "VAL-L[c]")
	)) # TOKB

	aminoAcidView = biomassContainer.countsView(
		[oneToThreeMapping[x] for x in kb._aaWeights.iterkeys() if x != "U"] # Ignore selenocysteine (TODO: Include it)
		)

	aaMmolPerGDCW = numpy.sum(
		kb.monomerData["aaCounts"] *
		numpy.tile(monomersView.counts().reshape(-1, 1), (1, 21)),
		axis = 0
		) * (
		(1 / kb.nAvogadro.magnitude) *
		(1000 / kb.avgCellDryMassInit.magnitude)
		)

	aminoAcidView.countsIs(
		aaMmolPerGDCW[range(19) + [20]]
		)

	# RNA fraction
	ntpView = biomassContainer.countsView(
		["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"]
		)

	ntpPerGDCW = numpy.sum(
		kb.rnaData["countsAUCG"] *
		numpy.tile(rnaView.counts().reshape(-1, 1), (1, 4)),
		axis = 0
		) * (
		(1 / kb.nAvogadro.magnitude) *
		(1000 / kb.avgCellDryMassInit.magnitude)
		)

	ntpView.countsIs(
		ntpPerGDCW
		)

	# DNA fraction

	dNtpView = biomassContainer.countsView(		# TODO: Better name so as not to confuse with bulkContainer view
		["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
		)

	dNtpPerGDCW = dNtpsView.counts() * (
		(1 / kb.nAvogadro.magnitude) *
		(1000 / kb.avgCellDryMassInit.magnitude)
		)

	dNtpView.countsIs(
		dNtpPerGDCW
		)

	# Glycogen fraction

	glycogenView = biomassContainer.countsView(
		list(kb.cellGlycogenFractionData["metaboliteId"])
		)

	glycogenMassFraction = float(dryComposition60min["glycogenMassFraction"])
	glycogenMass = kb.avgCellDryMassInit.magnitude * glycogenMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellGlycogenFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	glycogenPerGDCW = (
		glycogenMass * kb.cellGlycogenFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	glycogenView.countsIs(
		glycogenPerGDCW
		)

	# Murein fraction

	mureinView = biomassContainer.countsView(
		list(kb.cellMureinFractionData["metaboliteId"])
		)

	mureinMassFraction = float(dryComposition60min["mureinMassFraction"])
	mureinMass = kb.avgCellDryMassInit.magnitude * mureinMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellMureinFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	mureinPerGDCW = (
		mureinMass * kb.cellMureinFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	mureinView.countsIs(
		mureinPerGDCW
		)

	# LPS fraction

	lpsView = biomassContainer.countsView(
		list(kb.cellLPSFractionData["metaboliteId"])
		)

	lpsMassFraction = float(dryComposition60min["lpsMassFraction"])
	lpsMass = kb.avgCellDryMassInit.magnitude * lpsMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLPSFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	lpsPerGDCW = (
		lpsMass * kb.cellLPSFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	lpsView.countsIs(
		lpsPerGDCW
		)

	# Lipid fraction

	lipidView = biomassContainer.countsView(
		list(kb.cellLipidFractionData["metaboliteId"])
		)

	lipidMassFraction = float(dryComposition60min["lipidMassFraction"])
	lipidMass = kb.avgCellDryMassInit.magnitude * lipidMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLipidFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	lipidPerGDCW = (
		lipidMass * kb.cellLipidFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	lipidView.countsIs(
		lipidPerGDCW
		)

	# Inorganic ion fraction

	inorganicIonView = biomassContainer.countsView(
		list(kb.cellInorganicIonFractionData["metaboliteId"])
		)

	inorganicIonMassFraction = float(dryComposition60min["inorganicIonMassFraction"])
	inorganicIonMass = kb.avgCellDryMassInit.magnitude * inorganicIonMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellInorganicIonFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	inorganicIonPerGDCW = (
		inorganicIonMass * kb.cellInorganicIonFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	inorganicIonView.countsIs(
		inorganicIonPerGDCW
		)

	# Soluble pool fraction

	solublePoolView = biomassContainer.countsView(
		list(kb.cellSolublePoolFractionData["metaboliteId"])
		)

	solublePoolMassFraction = float(dryComposition60min["solublePoolMassFraction"])
	solublePoolMass = kb.avgCellDryMassInit.magnitude * solublePoolMassFraction

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellSolublePoolFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB

	solublePoolPerGDCW = (
		solublePoolMass * kb.cellSolublePoolFractionData["massFraction"]
		) / mws * (
		1000 / kb.avgCellDryMassInit.magnitude)

	solublePoolView.countsIs(
		solublePoolPerGDCW
		)

	# Validate

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in biomassContainer._objectNames
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB
	# TODO


	# TODO: Unhack this
	kb.wildtypeBiomass.struct_array["biomassFlux"] = biomassContainer.counts()

def normalize(array):
	return numpy.array(array).astype("float") / numpy.linalg.norm(array, 1)

def countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro):
	assert numpy.allclose(numpy.sum(relativeExpression), 1)
	return mass / numpy.dot(mws / nAvogadro, relativeExpression)

if __name__ == "__main__":
	import wholecell.utils.config
	import wholecell.utils.knowledgebase_fixture_manager

	kbDir = wholecell.utils.config.SIM_FIXTURE_DIR
	kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(
				os.path.join(kbDir, 'KnowledgeBase.cPickle'))
	kbFit = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(
				os.path.join(kbDir, 'KnowledgeBase.cPickle'))
	fitKb(kb)
