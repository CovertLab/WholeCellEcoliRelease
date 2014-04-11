
import numpy as np

def calcInitialConditions(sim, kb):
	# Data from KB

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	initialDryMass = kb.avgCellDryMassInit.to('g').magnitude

	rnaLength = np.sum(kb.rnaNTCounts, axis = 1)
	rnaExpression = kb.rnaExpression.to('dimensionless').magnitude
	rnaExpression /= np.sum(rnaExpression)

	fracInitFreeNTPs = kb.fracInitFreeNTPs.to('dimensionless').magnitude
	fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude
	biomass = kb.coreBiomass

	# Monomers are not complexes and not modified
	monomers = kb.bulkMolecules[kb.bulkMolecules['isProteinMonomer'] & np.logical_not(kb.bulkMolecules['isModifiedForm'])]
	monomerLength = np.sum(kb.proteinMonomerAACounts, axis = 1)

	monomerExpression = np.zeros(len(kb._proteinMonomerData), dtype = float)
	monomerExpression = [
		kb._rnaData['expression'][np.where(rnaId == kb._rnaData['id'])[0][0]]
		for rnaId in kb._proteinMonomerData['rnaId']
		]

	monomerExpression /= np.sum(monomerExpression)

	rnaIds = kb.bulkMolecules['moleculeId'][kb.bulkMolecules['isRna']]

	# Set bulk molecules

	bulk = sim.states['BulkMolecules']

	initialDryMass = initialDryMass + bulk.randStream.normal(0.0, 1e-15)

	feistCoreView = bulk.container.countsView(biomass['metaboliteId'])
	h2oView = bulk.container.countView('H2O[c]')
	ntpsView = bulk.container.countsView(IDS['ntps'])
	rnaView = bulk.container.countsView(rnaIds)
	aasView = bulk.container.countsView(IDS['aas'])
	monomersView = bulk.container.countsView(monomers['moleculeId'])

	# Set metabolite counts from Feist biomass
	feistCoreView.countsIs(
		np.round(
			np.fmax(biomass['biomassFlux'].to('mol/(DCWg*hr)').magnitude,0) * nAvogadro * initialDryMass
			)
		)

	# Set water
	h2oView.countIs(
		(6.7e-13 / 1.36 + bulk.randStream.normal(0, 1e-15)) / bulk._moleculeMass[bulk._moleculeIDs == 'H2O[c]'] * nAvogadro
		) # TOKB

	# Set RNA counts from expression levels
	ntpsToPolym = np.round(
		(1 - fracInitFreeNTPs) * np.sum(ntpsView.counts())
		)

	rnaCnts = bulk.randStream.mnrnd(
		np.round(ntpsToPolym / (np.dot(rnaExpression, rnaLength))),
		rnaExpression
		)

	ntpsView.countsIs(
		np.round(
			fracInitFreeNTPs * ntpsView.counts()
			)
		)

	rnaView.countsIs(rnaCnts)

	# Set protein counts from expression levels
	aasToPolym = np.round(
		(1 - fracInitFreeAAs) * np.sum(aasView.counts())
		)

	monCnts = bulk.randStream.mnrnd(
		np.round(aasToPolym / (np.dot(monomerExpression, monomerLength))),
		monomerExpression
		)

	aasView.countsIs(
		np.round(
			fracInitFreeAAs * aasView.counts()
			)
		)

	monomersView.countsIs(monCnts)

IDS = {
	'ntps':["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"],
	'ndps':["ADP[c]", "CDP[c]", "GDP[c]", "UDP[c]"],
	'nmps':["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"],
	'aas':[
		"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
		"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]"
		],
	'h2o':"H2O[c]",
	'h':"H[c]",
	'ppi':"PPI[c]",
	'adp':"ADP[c]",
	'pi':"PI[c]",
	'tRnas':[
		"gltV-tRNA", "gltT-tRNA", "gltW-tRNA", "gltU-tRNA", "glnU-tRNA", "glnW-tRNA", "glnX-tRNA", "glnV-tRNA", "serT-tRNA", "serW-tRNA", "selC-tRNA",
		"serU-tRNA", "serV-tRNA", "serX-tRNA", "RNA0-302", "lysV-tRNA", "RNA0-303", "RNA0-301", "lysW-tRNA", "lysT-tRNA", "RNA0-306", "metY-tRNA",
		"metW-tRNA", "metZ-tRNA", "metU-tRNA", "metT-tRNA", "thrW-tRNA", "thrV-tRNA", "thrU-tRNA", "thrT-tRNA", "trpT-tRNA", "pheV-tRNA",
		"pheU-tRNA", "glyV-tRNA", "glyY-tRNA", "glyU-tRNA", "glyT-tRNA", "glyX-tRNA", "glyW-tRNA", "proL-tRNA", "proK-tRNA", "proM-tRNA",
		"RNA0-300", "valU-tRNA", "valV-tRNA", "valX-tRNA", "valY-tRNA", "valT-tRNA", "valW-tRNA", "hisR-tRNA", "ileX-tRNA", "RNA0-305",
		"ileV-tRNA", "ileT-tRNA", "ileU-tRNA", "tyrV-tRNA", "tyrU-tRNA", "tyrT-tRNA", "alaX-tRNA", "alaW-tRNA", "alaT-tRNA", "alaV-tRNA",
		"alaU-tRNA", "argY-tRNA", "argZ-tRNA", "argX-tRNA", "argU-tRNA", "argV-tRNA", "argQ-tRNA", "argW-tRNA", "aspV-tRNA", "aspU-tRNA",
		"aspT-tRNA", "RNA0-304", "asnV-tRNA", "asnU-tRNA", "asnT-tRNA", "leuU-tRNA", "leuQ-tRNA", "leuX-tRNA", "leuV-tRNA", "leuT-tRNA",
		"leuZ-tRNA", "leuW-tRNA", "leuP-tRNA", "cysT-tRNA"
		],
	'rRnas':[
		"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]", "RRLD-RRNA[c]", "RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]",
		"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]", "RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]",
		"RRFA-RRNA[c]", "RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]", "RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
		],
	'rRna23Ss':[
		"RRLA-RRNA[c]", "RRLB-RRNA[c]", "RRLC-RRNA[c]", "RRLD-RRNA[c]", "RRLE-RRNA[c]", "RRLG-RRNA[c]", "RRLH-RRNA[c]",
		],
	'rRna16ss':[
		"RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]", "RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]",
		],
	'rRna5Ss':[
		"RRFA-RRNA[c]", "RRFB-RRNA[c]", "RRFC-RRNA[c]", "RRFD-RRNA[c]", "RRFE-RRNA[c]", "RRFF-RRNA[c]", "RRFG-RRNA[c]", "RRFH-RRNA[c]"
		],
	'FeistCore':[
		"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLN-L[c]", "GLU-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
		"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
		"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]", "CTP[c]", "GTP[c]", "UTP[c]", "ATP[c]", "MUREIN5PX4P[p]", "KDO2LIPID4[o]",
		"PE160[c]", "PE161[c]", "K[c]", "NH4[c]", "MG2[c]", "CA2[c]", "FE2[c]", "FE3[c]", "CU2[c]", "MN2[c]",
		"MOBD[c]", "COBALT2[c]", "ZN2[c]", "CL[c]", "SO4[c]", "PI[c]", "COA[c]", "NAD[c]", "NADP[c]", "FAD[c]",
		"THF[c]", "MLTHF[c]", "10FTHF[c]", "THMPP[c]", "PYDX5P[c]", "PHEME[c]", "SHEME[c]", "UDCPDP[c]", "AMET[c]", "2OHPH[c]",
		"RIBFLV[c]"
		],
	} # TOKB
