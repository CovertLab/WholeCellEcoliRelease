
import numpy as np

ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
aaIds = ["ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]",
	"GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
	"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]",
	"TRP-L[c]", "TYR-L[c]", "VAL-L[c]"]


def calcInitialConditions(sim, kb):
	randStream = sim.randStream

	bulk = sim.states['BulkMolecules']

	initializeBulk(bulk.container, kb, randStream)


def initializeBulk(bulkContainer, kb, randStream):

	## Set metabolite counts from Feist biomass
	initializeBulkBiomass(kb, bulkContainer, randStream)

	## Set water
	initializeBulkWater(kb, bulkContainer, randStream)

	## Set RNA counts from expression levels
	initializeBulkRNA(kb, bulkContainer, randStream)
	initializeBulkNTPs(kb, bulkContainer, randStream)

	## Set protein counts from expression levels
	initializeBulkMonomers(kb, bulkContainer, randStream)
	initializeBulkAAs(kb, bulkContainer, randStream)


def initializeBulkWater(kb, bulkContainer, randStream):
	h2oView = bulkContainer.countView('H2O[c]')

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	mwH2O = kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "H2O[c]"].magnitude
	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	h2oView.countIs(
		(avgCellWaterMassInit + randStream.normal(0, 1e-15)) / mwH2O * nAvogadro
		)


def initializeBulkBiomass(kb, bulkContainer, randStream):
	biomassMetabolites = kb.coreBiomass['metaboliteId']

	feistCoreView = bulkContainer.countsView(biomassMetabolites)

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	initialDryMass = kb.avgCellDryMassInit.to('g').magnitude + randStream.normal(0.0, 1e-15)
	biomassFlux = kb.wildtypeBiomass['biomassFlux'].to('mol/(DCW_g)').magnitude


	feistCoreView.countsIs(
		np.round(
			np.fmax(biomassFlux, 0) * nAvogadro * initialDryMass
			)
		)


def initializeBulkRNA(kb, bulkContainer, randStream):
	rnaIds = kb.bulkMolecules['moleculeId'][kb.bulkMolecules['isRnaMonomer']]

	ntpsView = bulkContainer.countsView(ntpIds)
	rnaView = bulkContainer.countsView(rnaIds)

	fracInitFreeNTPs = kb.fracInitFreeNTPs.to('dimensionless').magnitude
	rnaLength = np.sum(kb.rnaData['countsAUCG'], axis = 1)
	rnaExpression = kb.rnaExpression['expression'].to('dimensionless').magnitude
	rnaExpression /= np.sum(rnaExpression)

	ntpsToPolym = np.round(
		(1 - fracInitFreeNTPs) * np.sum(ntpsView.counts())
		)
	try:
		rnaCnts = randStream.mnrnd(
			np.round(ntpsToPolym / (np.dot(rnaExpression, rnaLength))),
			rnaExpression
			)
	except:
		import ipdb; ipdb.set_trace()

	rnaView.countsIs(rnaCnts)


def initializeBulkNTPs(kb, bulkContainer, randStream):
	ntpsView = bulkContainer.countsView(ntpIds)

	fracInitFreeNTPs = kb.fracInitFreeNTPs.to('dimensionless').magnitude

	ntpsView.countsIs(
		np.round(
			fracInitFreeNTPs * ntpsView.counts()
			)
		)


def initializeBulkMonomers(kb, bulkContainer, randStream):
	## Monomers are not complexes and not modified
	monomers = kb.bulkMolecules[kb.bulkMolecules['isProteinMonomer'] & ~kb.bulkMolecules['isModifiedForm']]
	rnapIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

	aasView = bulkContainer.countsView(aaIds)
	monomersView = bulkContainer.countsView(monomers['moleculeId'])
	rnapView = bulkContainer.countsView(rnapIds)

	fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude

	monomerExpression = kb.rnaExpression[kb.rnaIndexToMonomerMapping].magnitude
	monomerExpression /= np.sum(monomerExpression)

	monomerLength = np.sum(kb.proteinMonomerAACounts, axis = 1)

	aasToPolym = np.round(
		(1 - fracInitFreeAAs) * np.sum(aasView.counts())
		)

	monCnts = randStream.mnrnd(
		np.round(aasToPolym / (np.dot(monomerExpression, monomerLength))),
		monomerExpression
		)

	monomersView.countsIs(monCnts)
	import ipdb; ipdb.set_trace()
	kb.rnaExpression[kb.monomerIndexToRnaMapping[[np.where(x == kb.monomerData["id"])[0][0] for x in rnapIds]]]


def initializeBulkAAs(kb, bulkContainer, randStream):
	aasView = bulkContainer.countsView(aaIds)

	fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude

	aasView.countsIs(
		np.round(
			fracInitFreeAAs * aasView.counts()
			)
		)