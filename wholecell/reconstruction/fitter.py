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
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
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

	dryComposition60min = kb.cellDryMassComposition[kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60]

	### RNA Mass Fractions ###
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit.to('DCW_g') * rnaMassFraction

	## 23S rRNA Mass Fractions ##
	rRna23SMassFraction = 0.525 # TOKB # This is the fraction of RNA that is 23S rRNA

	# Assume all 23S rRNAs are expressed equally
	rRna23SExpression = normalize(numpy.ones(rRna23SView.counts().size))

	nRRna23Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * rRna23SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna23S"]].to('g/mol').magnitude,
		rRna23SExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)
	
	## 16S rRNA Mass Fractions ##
	rRna16SMassFraction = 0.271 # TOKB # This is the fraction of RNA that is 16S rRNA

	# Assume all 16S rRNAs are expressed equally
	rRna16SExpression = normalize(numpy.ones(rRna16SView.counts().size))

	nRRna16Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * rRna16SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna16S"]].to('g/mol').magnitude,
		rRna16SExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	## 5S rRNA Mass Fractions ##
	rRna5SMassFraction = 0.017 # TOKB # This is the fraction of RNA that is 5S rRNA

	# Assume all 5S rRNAs are expressed equally
	rRna5SExpression = normalize(numpy.ones(rRna5SView.counts().size))

	nRRna5Ss = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * rRna5SMassFraction,
		kb.rnaData["mw"][kb.rnaData["isRRna5S"]].to('g/mol').magnitude,
		rRna5SExpression,
		kb.nAvogadro.to('1/mol').magnitude
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
		rnaMass.to('DCW_g').magnitude * tRnaMassFraction,
		kb.rnaData["mw"][kb.rnaData["isTRna"]].to('g/mol').magnitude,
		tRnaExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	tRnaView.countsIs((nTRnas * tRnaExpression))

	## mRNA Mass Fractions ##
	mRnaMassFraction = 0.041 # TOKB # This is the fraction of RNA that is mRNA

	mRnaExpression = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])

	nMRnas = countsFromMassAndExpression(
		rnaMass.to('DCW_g').magnitude * mRnaMassFraction,
		kb.rnaData["mw"][kb.rnaData["isMRna"]].to('g/mol').magnitude,
		mRnaExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	mRnaView.countsIs((nMRnas * mRnaExpression))


	### Protein Mass fraction ###

	monomersView = bulkContainer.countsView(kb.monomerData["id"])

	monomerMassFraction = float(dryComposition60min["proteinMassFraction"])
	monomerMass = kb.avgCellDryMassInit * monomerMassFraction

	monomerExpression = normalize(kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping])

	nMonomers = countsFromMassAndExpression(
		monomerMass.to("DCW_g").magnitude,
		kb.monomerData["mw"].to('g/mol').magnitude,
		monomerExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	monomersView.countsIs((nMonomers * monomerExpression))

	
	### DNA Mass fraction ###
	# TODO: Don't just keep dNTPs in the soluble pool
	# TODO: Compute based on chromosome sequence
	dNtpIds = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	dNtpsView = bulkContainer.countsView(dNtpIds)

	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit * dnaMassFraction

	dNtpRelativeAmounts = normalize(numpy.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		]))

	dNtpIdxs = [numpy.where(kb.bulkMolecules["moleculeId"] == idx)[0][0] for idx in dNtpIds]

	dNtpMws = kb.bulkMolecules["mass"][dNtpIdxs]

	nDNtps = countsFromMassAndExpression(
		dnaMass.to('DCW_g').magnitude,
		dNtpMws.to('g/mol').magnitude,
		dNtpRelativeAmounts,
		kb.nAvogadro.to('1/mol').magnitude
		)

	dNtpsView.countsIs(nDNtps * dNtpRelativeAmounts)

	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]
	rnapView = bulkContainer.countsView(rnapIds)

	## Number of ribosomes needed ##
	monomerLengths = numpy.sum(kb.monomerData['aaCounts'], axis = 1)
	nRibosomesNeeded = numpy.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			numpy.log(2) / kb.cellCycleLen
			) * monomersView.counts()
		).to('dimensionless').magnitude
	
	if numpy.sum(rRna23SView.counts()) < nRibosomesNeeded:
		raise NotImplementedError, "Cannot handle having too few ribosomes"

	## Number of RNA Polymerases ##
	rnaLengths = numpy.sum(kb.rnaData['countsAUCG'], axis = 1)

	nRnapsNeeded = numpy.sum(
		rnaLengths / kb.rnaPolymeraseElongationRate * (
			numpy.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
			) * rnaView.counts()
		).to('dimensionless').magnitude

	# nRnapsNeeded = 600	

	minRnapCounts = (
		nRnapsNeeded * numpy.array([2, 1, 1, 1]) # Subunit stoichiometry
		)

	rnapView.countsIs(
		numpy.fmax(rnapView.counts(), minRnapCounts)
		)

	
	### Modify kbFit to reflect our bulk container ###

	## Fraction of active Ribosomes ##
	kb.parameters["fracActiveRibosomes"] = Q_(float(nRibosomesNeeded) / numpy.sum(rRna23SView.counts()), "dimensionless")
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
				numpy.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
				) * rnaView.counts()
			).to('dimensionless').magnitude
		)

	kb.rnaData["synthProb"][:] = synthProb


	## Full WT Biomass function ##

	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = numpy.dtype("float64")
		)

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
		) # TODO: Don't use private variable from KB. Use AMINO_ACID_1_TO_3_ORDERED order.

	aaMmolPerGDCW = (
			numpy.sum(
				kb.monomerData["aaCounts"] *
				numpy.tile(monomersView.counts().reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / kb.nAvogadro.to('amino_acid/mmol')) *
				(1 / kb.avgCellDryMassInit)
			)
		).to('mmol/DCW_g')
	
	# TODO: Skipping selenocystine (U) here. Re add this!
	aminoAcidView.countsIs(
		aaMmolPerGDCW[range(19) + [20]].magnitude
		)

	# RNA fraction
	ntpView = biomassContainer.countsView(
		["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"]
		)

	ntpMmolPerGDCW = (
			numpy.sum(
				kb.rnaData["countsAUCG"] *
				numpy.tile(rnaView.counts().reshape(-1, 1), (1, 4)),
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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellGlycogenFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB

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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellMureinFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB

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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLPSFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB

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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellLipidFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB

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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellInorganicIonFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB
	
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

	bulkMoleculesIdxs = numpy.array([
		numpy.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.cellSolublePoolFractionData["metaboliteId"]
		])
	mws = kb.bulkMolecules["mass"][bulkMoleculesIdxs] # TOKB
	
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
	assert type(mass) != Q_
	assert type(mws) != Q_
	assert type(relativeExpression) != Q_
	assert type(nAvogadro) != Q_
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
