#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

# TODO: break this into multiple files, it's becoming unbearably long

import re

import numpy as np
import tables

import wholecell.states.state
import wholecell.states.partition
import wholecell.utils.bulk_objects_container


class BulkMolecules(wholecell.states.state.State):
	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'BulkMolecules',
			'name':'Bulk Molecules',
			'dynamics':[],
			'units':{}
			}

		self.time = None
		self.partitionClass = BulkMoleculesPartition

		self._container = None

		self._moleculeMass = None

		self._moleculeIDs = None
		self._compartmentIDs = None

		self._nCompartments = None

		self._countsRequested = None
		self._countsAllocated = None
		self._countsReturned = None
		self._countsUnallocated = None

		self._typeIdxs = None
		self._typeLocalizations = None

		self._rnaLens = None
		self._rnaExp = None

		self._monLens = None
		self._monExp = None

		super(BulkMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(BulkMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# HACK

		# kb does not yet have feist core vals!
		if hasattr(kb, 'feistCoreVals'):
			self.feistCoreVals = kb.feistCoreVals

		else:
			self.feistCoreVals = FEIST_CORE_VALS.copy()

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!

		self._moleculeIDs = []

		self._moleculeIDs += [molecule['id'] for molecule in kb.metabolites]
		self._moleculeIDs += [molecule['id'] for molecule in kb.rnas]
		self._moleculeIDs += [molecule['id'] for molecule in kb.proteins]

		self._compartmentIDs = [compartment['id'] for compartment in COMPARTMENTS]
		self._nCompartments = len(self._compartmentIDs)

		self._moleculeMass = np.array(sum(
			(
				[molecule['mw7.2'] for molecule in kb.metabolites],
				[molecule['mw'] for molecule in kb.rnas],
				[molecule['mw'] for molecule in kb.proteins]
			),
			[]), np.float64)

		self._moleculeMass[np.where(self._moleculeMass < 0)] = 0

		# Information needed for calcInitialConditions

		self._typeIdxs = {}

		self._typeIdxs.update({
			'metabolites':np.arange(len(kb.metabolites)),
			'rnas':np.arange(len(kb.rnas))+len(kb.metabolites),
			'proteins':np.arange(len(kb.proteins))+len(kb.metabolites)+len(kb.rnas),
			})

		self._typeIdxs['water'] = self._moleculeIDs.index('H2O')

		self._typeLocalizations = {}

		self._typeLocalizations.update({
			'proteins':[mol['location'] for mol in kb.proteins],
			})

		# Values needed for calcInitialConditions
		self._rnaLens = np.array([np.sum(rna["ntCount"]) for rna in kb.rnas])
		self._rnaExp = np.array([x["expression"] for x in kb.rnas])
		self._rnaExp /= np.sum(self._rnaExp)

		monomers = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]
		self._monLens = np.array([np.sum(monomer["aaCount"]) for monomer in monomers])
		rnaIdToExp = dict([(x["id"], x["expression"]) for x in kb.rnas if x["monomerId"] != None])
		self._monExp = np.array([rnaIdToExp[monomer["rnaId"]] for monomer in monomers])
		self._monExp /= np.sum(self._monExp)

		self._typeIdxs['monomers'] = np.array([self._moleculeIDs.index(x["id"]) for x in monomers])
		self._typeLocalizations['monomers'] = [monomer['location'] for monomer in monomers]


	def allocate(self):
		super(BulkMolecules, self).allocate() # Allocates partitions

		# Create the container for molecule counts
		self._container = wholecell.utils.bulk_objects_container.BulkObjectsContainer([
			'{}[{}]'.format(moleculeID, compartmentID)
			for moleculeID in self._moleculeIDs
			for compartmentID in self._compartmentIDs
			])

		nMolecules = self._container._counts.size
		nPartitions = len(self.partitions)
		dtype = self._container._counts.dtype

		# Arrays for tracking values related to partitioning
		self._countsRequested = np.zeros((nMolecules, nPartitions), dtype)
		self._countsAllocated = np.zeros((nMolecules, nPartitions), dtype)
		self._countsReturned = np.zeros((nMolecules, nPartitions), dtype)
		self._countsUnallocated = np.zeros(nMolecules, dtype)

	
	def calcInitialConditions(self):
		from wholecell.utils.constants import Constants

		initialDryMass = INITIAL_DRY_MASS + self.randStream.normal(0.0, 1e-15)

		feistCoreView = self._container.countsView(IDS['FeistCore'])
		h2oView = self._container.countView('H2O[c]')
		ntpsView = self._container.countsView(IDS['ntps'])
		matureRnaView = self._container.countsView(
			[self._moleculeIDs[i] + '[c]' for i in self._typeIdxs['rnas']])
		aasView = self._container.countsView(IDS['aas'])
		monomersView = self._container.countsView([
			self._moleculeIDs[index] + '[{}]'.format(self._typeLocalizations['monomers'][i])
			for i, index in enumerate(self._typeIdxs['monomers'])
			])

		# Set metabolite counts from Feist core
		feistCoreView.countsIs(
			np.round(
				self.feistCoreVals * 1e-3 * Constants.nAvogadro * initialDryMass
				)
			)

		# Set water
		h2oView.countIs(
			(6.7e-13 / 1.36 + self.randStream.normal(0, 1e-15)) / self._moleculeMass[self._typeIdxs['water']] * Constants.nAvogadro
			) # TOKB

		# Set RNA counts from expression levels
		ntpsToPolym = np.round(
			(1 - FRAC_INIT_FREE_NTPS) * np.sum(ntpsView.counts())
			)

		rnaCnts = self.randStream.mnrnd(
			np.round(ntpsToPolym / (np.dot(self._rnaExp, self._rnaLens))),
			self._rnaExp
			)

		ntpsView.countsIs(
			np.round(
				FRAC_INIT_FREE_NTPS * ntpsView.counts()
				)
			)

		matureRnaView.countsIs(rnaCnts)

		# Set protein counts from expression levels
		aasToPolym = np.round(
			(1 - FRAC_INIT_FREE_AAS) * np.sum(aasView.counts())
			)

		monCnts = self.randStream.mnrnd(
			np.round(aasToPolym / (np.dot(self._monExp, self._monLens))),
			self._monExp
			)

		aasView.countsIs(
			np.round(
				FRAC_INIT_FREE_AAS * aasView.counts()
				)
			)

		monomersView.countsIs(monCnts)


	def partition(self):
		# TODO
		pass


	def merge(self):
		# TODO
		pass

	
	def mass(self, typeKey = None):
		if typeKey is None:
			indexes = np.s_[:]

		else:
			indexes = self._typeIdxs[typeKey]

		return np.dot(
			self._moleculeMass[indexes],
			self._container._counts.view().reshape((-1, self._nCompartments))[indexes, :]
			)


	def pytablesCreate(self, h5file, expectedRows):
		# TODO
		pass


	def pytablesAppend(self, h5file):
		# TODO
		pass


	def pytablesLoad(self, h5file, timePoint):
		# TODO
		pass


class BulkMoleculesPartition(wholecell.states.partition.Partition):
	def __init__(self, *args, **kwargs):
		self._container = None
		self._isReqAbs = None

		self._indexMapping = None

		super(UniqueMoleculesPartition, self).__init__(*args, **kwargs)


	def initialize(self, moleculeNames, isReqAbs = False):
		self._container = wholecell.utils.bulk_objects_container.BulkObjectsContainer(
			moleculeNames)
		self._isReqAbs = isReqAbs

		self._indexMapping = np.array(
			self._state._container._getIndexes(moleculeNames)
			)


	def request(self, target):
		self._process.requestBulkMolecules()

		target[self._indexMapping] = self._container._counts # direct reference to avoid a copy operation


	def allocationIs(self, source):
		self._container.countsIs(source[self._indexMapping])


	def returned(self, target):
		target[self._indexMapping] = self._container._counts


# Constants (should be taken from the KB)

FEIST_CORE_VALS = np.array([ # TODO: This needs to go in the KB
	0.513689, 0.295792, 0.241055, 0.241055, 0.091580, 0.263160, 0.263160, 0.612638, 0.094738, 0.290529,
	0.450531, 0.343161, 0.153686, 0.185265, 0.221055, 0.215792, 0.253687, 0.056843, 0.137896, 0.423162,
	0.026166, 0.027017, 0.027017, 0.026166, 0.133508, 0.215096, 0.144104, 0.174831, 0.013894, 0.019456,
	0.063814, 0.075214, 0.177645, 0.011843, 0.007895, 0.004737, 0.007106, 0.007106, 0.003158, 0.003158,
	0.003158, 0.003158, 0.003158, 0.004737, 0.003948, 0.003948, 0.000576, 0.001831, 0.000447, 0.000223,
	0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000055, 0.000223, 0.000223,
	0.000223		# mmol/gDCW (supp info 3, "biomass_core", column G)
	]) # TOKB

INITIAL_DRY_MASS = 2.8e-13 / 1.36 # TOKB

COMPARTMENTS = [ # TODO: move to KB
	{"id": "c", "name": "Cytosol"},
	{"id": "e", "name": "Extracellular space"},
	{"id": "i", "name": "Inner membrane"},
	{"id": "j", "name": "Projection"},
	{"id": "l", "name": "Pilus"},
	{"id": "m", "name": "Membrane"},
	{"id": "n", "name": "Nucleoid"},
	{"id": "o", "name": "Outer membrane"},
	{"id": "p", "name": "Periplasm"},
	{"id": "w", "name": "Cell wall"}
	]

FRAC_INIT_FREE_NTPS = 0.0015
FRAC_INIT_FREE_AAS = 0.001

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
		]
	} # TOKB
