
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
	# Data from KB

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	initialDryMass = kb.avgCellDryMassInit.to('g').magnitude
	mwH2O = kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "H2O[c]"].magnitude

	biomassFlux = kb.coreBiomass['biomassFlux'].to('mol/(DCW_g*hr)').magnitude
	biomassMetabolites = kb.coreBiomass['metaboliteId']

	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	# Set bulk molecules

	initialDryMass = initialDryMass + randStream.normal(0.0, 1e-15)

	feistCoreView = bulkContainer.countsView(biomassMetabolites)
	h2oView = bulkContainer.countView('H2O[c]')
	ntpsView = bulkContainer.countsView(ntpIds)

	## Set metabolite counts from Feist biomass
	feistCoreView.countsIs(
		np.round(
			np.fmax(biomassFlux, 0) * nAvogadro * initialDryMass
			)
		)

	## Set water
	h2oView.countIs( # TODO: get mass of water directly from the KB
		(avgCellWaterMassInit + randStream.normal(0, 1e-15)) / mwH2O * nAvogadro
		) # TOKB

	## Set RNA counts from expression levels
	initializeBulkRNA(kb, bulkContainer, randStream)
	initializeBulkNTPs(kb, bulkContainer, randStream)

	## Set protein counts from expression levels
	initializeBulkMonomers(kb, bulkContainer, randStream)
	initializeBulkAAs(kb, bulkContainer, randStream)


def initializeBulkRNA(kb, bulkContainer, randStream):
	rnaIds = kb.bulkMolecules['moleculeId'][kb.bulkMolecules['isRna']]

	ntpsView = bulkContainer.countsView(ntpIds)
	rnaView = bulkContainer.countsView(rnaIds)

	fracInitFreeNTPs = kb.fracInitFreeNTPs.to('dimensionless').magnitude
	rnaLength = np.sum(kb.rnaNTCounts, axis = 1)
	rnaExpression = kb.rnaExpression.to('dimensionless').magnitude
	rnaExpression /= np.sum(rnaExpression)

	ntpsToPolym = np.round(
		(1 - fracInitFreeNTPs) * np.sum(ntpsView.counts())
		)

	rnaCnts = randStream.mnrnd(
		np.round(ntpsToPolym / (np.dot(rnaExpression, rnaLength))),
		rnaExpression
		)

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

	aasView = bulkContainer.countsView(aaIds)
	monomersView = bulkContainer.countsView(monomers['moleculeId'])

	fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude
	monomerExpression = np.zeros(len(kb._proteinMonomerData), dtype = float)
	monomerExpression = [
		kb._rnaData['expression'][np.where(rnaId == kb._rnaData['id'])[0][0]]
		for rnaId in kb._proteinMonomerData['rnaId']
		] # TODO: I'm almost certain this data is already parsed in the KB, get it from there

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

def initializeBulkAAs(kb, bulkContainer, randStream):
	aasView = bulkContainer.countsView(aaIds)

	fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude

	aasView.countsIs(
		np.round(
			fracInitFreeAAs * aasView.counts()
			)
		)