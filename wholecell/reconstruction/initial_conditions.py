
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

	## Set dNTPs
	initializeDNA(bulkContainer, kb, randStream)

	## Set other biomass components
	initializeBulkComponents(bulkContainer, kb, randStream)

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
	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	dNTPs = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]

	dNTPBulkView = bulkContainer.countsView(dNTPs)

	dNTPBiomassView = biomassContainer.countsView(dNTPs)

	dNTPBulkView.countsIs((
		kb.avgCellDryMassInit.to("DCW_gram").magnitude *
		dNTPBiomassView.counts() *
		kb.nAvogadro.to("1 / millimole").magnitude
		))


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


def initializeBulkWater(kb, bulkContainer, randStream):
	h2oView = bulkContainer.countView('H2O[c]')

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	mwH2O = kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "H2O[c]"].magnitude
	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	h2oView.countIs(
		(avgCellWaterMassInit) / mwH2O * nAvogadro
		)