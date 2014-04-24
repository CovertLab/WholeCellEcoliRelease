
import numpy as np
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.reconstruction.fitter import countsFromMassAndExpression
from wholecell.reconstruction.fitter import normalize


def calcInitialConditions(sim, kb):
	randStream = sim.randStream

	bulk = sim.states['BulkMolecules']

	initializeBulk(bulk.container, kb, randStream)


def initializeBulk(bulkContainer, kb, randStream):

	## Set protein counts from expression
	initializeProteinMonomers(bulkContainer, kb, randStream)

	## Set RNA counts from expression
	initializeRNA(bulkContainer, kb, randStream)

	## Set DNA
	initializeDNA(bulkContainer, kb, randStream)

	## Set other biomass components
	initializeBulkComponents(bulkContainer, kb, randStream)

	## Set pools
	initializePools(bulkContainer, kb, randStream)

	## Set water
	initializeBulkWater(kb, bulkContainer, randStream)


def initializeProteinMonomers(bulkContainer, kb, randStream):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

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

	monomersView.countsIs(
		randStream.mnrnd(nMonomers, monomerExpression)
		)

	# monomersView.countsIs(nMonomers * monomerExpression)


def initializeRNA(bulkContainer, kb, randStream):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit.magnitude * rnaMassFraction

	rnaExpression = normalize(kb.rnaExpression['expression'])

	nRnas = countsFromMassAndExpression(
		rnaMass,
		kb.rnaData["mw"],
		rnaExpression,
		kb.nAvogadro.magnitude
		)

	rnaView.countsIs(
		randStream.mnrnd(nRnas, rnaExpression)
		)

	# rnaView.countsIs(nRnas * rnaExpression)


def initializeDNA(bulkContainer, kb, randStream):

	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	dNTPs = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	dnmpIds = ["DAMP[n]", "DCMP[n]", "DGMP[n]", "DTMP[n]"]

	dnmpsView = bulkContainer.countsView(dnmpIds)
	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit.magnitude * dnaMassFraction

	dnaExpression = normalize(np.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		], dtype = np.float64))

	mws = np.array([
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == x][0].magnitude for x in dnmpIds]
		) # This is a hack. Without a real chromosome, though, it's all a hack

	nDntps = countsFromMassAndExpression(
		dnaMass,
		mws,
		dnaExpression,
		kb.nAvogadro.magnitude
		)

	dnmpsView.countsIs(
		randStream.mnrnd(nDntps, dnaExpression)
		)


def initializeBulkComponents(bulkContainer, kb, randStream):
	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	notPRDMetabolites = (
		list(kb.cellGlycogenFractionData["metaboliteId"]) +
		list(kb.cellMureinFractionData["metaboliteId"]) +
		list(kb.cellLPSFractionData["metaboliteId"]) +
		list(kb.cellLipidFractionData["metaboliteId"]) +
		list(kb.cellInorganicIonFractionData["metaboliteId"]) +
		list(kb.cellSolublePoolFractionData["metaboliteId"])
		)

	notPRDBulkView = bulkContainer.countsView(notPRDMetabolites)

	notPRDBiomassView = biomassContainer.countsView(notPRDMetabolites)

	notPRDBulkView.countsIs((
		kb.avgCellDryMassInit.to("DCW_gram").magnitude *
		notPRDBiomassView.counts() *
		kb.nAvogadro.to("1 / millimole").magnitude
		))


def initializePools(bulkContainer, kb, randStream):
	# Note: This is adding dry biomass, so the cell will appear heavier

	from wholecell.reconstruction.knowledge_base_ecoli import AMINO_ACID_1_TO_3_ORDERED

	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
	dntpIds = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	aaIds = [x for x in AMINO_ACID_1_TO_3_ORDERED.values() if x != "SEC-L[c]"]

	ntpsBiomassView = biomassContainer.countsView(ntpIds)
	dntpsBiomassView = biomassContainer.countsView(dntpIds)
	aasBiomassView = biomassContainer.countsView(aaIds)

	ppiBulkView = bulkContainer.countView("PPI[c]")
	ntpsBulkView = bulkContainer.countsView(ntpIds)
	dntpsBulkView = bulkContainer.countsView(dntpIds)
	aasBulkView = bulkContainer.countsView(aaIds)

	dt = kb.timeStep.to("second").magnitude
	tau_d = kb.cellCycleLen.to("second").magnitude

	## NTPs
	ntpsFromLastStep = (
		ntpsBiomassView.counts() *
		(1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	ntpsBulkView.countsIs(ntpsFromLastStep)

	## dNTPs
	dntpsFromLastStep = (
		dntpsBiomassView.counts() * (1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	dntpsBulkView.countsIs(dntpsFromLastStep)

	## Amino Acids
	aasFromLastStep = (
		aasBiomassView.counts() * (1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	aasBulkView.countsIs(aasFromLastStep)

	## PPI
	# Assumption:
	# NTPs and dNTPs from two steps ago got completely used up in the last step
	ntpsFromTwoStepsAgo = (
		ntpsBiomassView.counts() *
		(np.exp(-np.log(2) / tau_d * dt) - np.exp(-np.log(2) / tau_d * 2 * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)

	dntpsFromTwoStepsAgo = (
		dntpsBiomassView.counts() *
		(np.exp(-np.log(2) / tau_d * dt) - np.exp(-np.log(2) / tau_d * 2 * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)

	ppiFromNtps = np.round(np.sum(
		ntpsFromLastStep
		))

	ppiFromDntps = np.round(np.sum(
		dntpsFromLastStep
		))

	ppiBulkView.countIs(ppiFromNtps + ppiFromDntps)


def initializeBulkWater(kb, bulkContainer, randStream):
	h2oView = bulkContainer.countView('H2O[c]')

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	mwH2O = kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "H2O[c]"].magnitude
	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	h2oView.countIs(
		(avgCellWaterMassInit) / mwH2O * nAvogadro
		)