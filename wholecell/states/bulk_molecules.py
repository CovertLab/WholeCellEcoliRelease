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

from __future__ import division

import re

import numpy as np
import tables

import wholecell.states.state
import wholecell.views.view
from wholecell.containers.bulk_objects_container import BulkObjectsContainer


class BulkMolecules(wholecell.states.state.State):
	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'BulkMolecules',
			'name':'Bulk Molecules',
			'dynamics':[],
			'units':{}
			}

		self.time = None

		self.container = None

		self._moleculeMass = None #

		self._moleculeIDs = None #
		self._rnaIds = None
		self._monomers = None
		self._compartmentIDs = None #

		self._nCompartments = None #

		self._isRequestAbsolute = None

		self._countsRequested = None
		self._countsAllocatedInitial = None
		self._countsAllocatedFinal = None
		self._countsUnallocated = None

		self._typeIdxs = None

		self._rnaLength = None #
		self._rnaExpression = None #

		self._monomerLength = None #
		self._monomerExpression = None #

		# Move '#' to KB

		super(BulkMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(BulkMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# Load constants
		self.nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
		self.initialDryMass = kb.avgCellDryMassInit.to('g').magnitude
		self.fracInitFreeNTPs = kb.fracInitFreeNTPs.to('dimensionless').magnitude
		self.fracInitFreeAAs = kb.fracInitFreeAAs.to('dimensionless').magnitude
		self.biomass = kb.coreBiomass

		self._moleculeIDs = kb.bulkMolecules['moleculeId']
		self._rnaIds = kb.bulkMolecules['moleculeId'][kb.bulkMolecules['isRna']]
		self._compartmentIDs = kb.compartments['compartmentAbbreviation']
		self._nCompartments = kb.nCompartments

		self._moleculeMass = kb.bulkMolecules['mass'].to('g/mol').magnitude

		self._typeIdxs = {'metabolites'	:	kb.bulkMolecules['isMetabolite'],
							'rnas'		:	kb.bulkMolecules['isRna'],
							'proteins'	:	kb.bulkMolecules['isProteinMonomer'],
							'water'		:	kb.bulkMolecules['isWater']}

		# Create the container for molecule counts
		self.container = BulkObjectsContainer(self._moleculeIDs)

		# Values needed for calcInitialConditions
		self._rnaLength = np.sum(kb.rnaNTCounts, axis = 1)
		self._rnaExpression = kb.rnaExpression.to('dimensionless').magnitude
		self._rnaExpression /= np.sum(self._rnaExpression)

		# Monomers are not complexes and not modified
		self._monomers = kb.bulkMolecules[kb.bulkMolecules['isProteinMonomer'] & np.logical_not(kb.bulkMolecules['isModifiedForm'])]
		self._monomerLength = np.sum(kb.proteinMonomerAACounts, axis = 1)

		self._monomerExpression = np.zeros(len(kb._proteinMonomerData), dtype = float)
		self._monomerExpression = [
			kb._rnaData['expression'][np.where(rnaId == kb._rnaData['id'])[0][0]]
			for rnaId in kb._proteinMonomerData['rnaId']
			]

		self._monomerExpression /= np.sum(self._monomerExpression)
		
		# TODO: restore this behavior or replace it with something bettter

		self._isRequestAbsolute = np.zeros(self._nProcesses, np.bool)
		try:
			self._isRequestAbsolute[sim.processes['RnaDegradation']._processIndex] = True

		except KeyError:
			pass

	def allocate(self):
		super(BulkMolecules, self).allocate() # Allocates partitions

		nMolecules = self.container._counts.size
		dtype = self.container._counts.dtype

		# Arrays for tracking values related to partitioning
		self._countsRequested = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedInitial = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsAllocatedFinal = np.zeros((nMolecules, self._nProcesses), dtype)
		self._countsUnallocated = np.zeros(nMolecules, dtype)

	
	def calcInitialConditions(self):
		initialDryMass = self.initialDryMass + self.randStream.normal(0.0, 1e-15)

		feistCoreView = self.container.countsView(self.biomass['metaboliteId'])
		h2oView = self.container.countView('H2O[c]')
		ntpsView = self.container.countsView(IDS['ntps'])
		rnaView = self.container.countsView(self._rnaIds)
		aasView = self.container.countsView(IDS['aas'])
		monomersView = self.container.countsView(self._monomers['moleculeId'])

		# Set metabolite counts from Feist biomass
		feistCoreView.countsIs(
			np.round(
				np.fmax(self.biomass['biomassFlux'].to('mol/(DCWg*hr)').magnitude,0) * self.nAvogadro * initialDryMass
				)
			)

		# Set water
		h2oView.countIs(
			(6.7e-13 / 1.36 + self.randStream.normal(0, 1e-15)) / self._moleculeMass[self._moleculeIDs == 'H2O[c]'] * self.nAvogadro
			) # TOKB

		# Set RNA counts from expression levels
		ntpsToPolym = np.round(
			(1 - self.fracInitFreeNTPs) * np.sum(ntpsView.counts())
			)

		rnaCnts = self.randStream.mnrnd(
			np.round(ntpsToPolym / (np.dot(self._rnaExpression, self._rnaLength))),
			self._rnaExpression
			)

		ntpsView.countsIs(
			np.round(
				self.fracInitFreeNTPs * ntpsView.counts()
				)
			)

		rnaView.countsIs(rnaCnts)

		# Set protein counts from expression levels
		aasToPolym = np.round(
			(1 - self.fracInitFreeAAs) * np.sum(aasView.counts())
			)

		monCnts = self.randStream.mnrnd(
			np.round(aasToPolym / (np.dot(self._monomerExpression, self._monomerLength))),
			self._monomerExpression
			)

		aasView.countsIs(
			np.round(
				self.fracInitFreeAAs * aasView.counts()
				)
			)

		monomersView.countsIs(monCnts)

	def updateQueries(self):
		for view in self._views:
			view._totalIs(self.container._counts[view.containerIndexes])


	def partition(self):
		if self._nProcesses:
			# Calculate and store requests
			self._countsRequested[:] = 0

			for view in self._views:
				self._countsRequested[view.containerIndexes, view._processIndex] += view._request()

			calculatePartition(self._isRequestAbsolute, self._countsRequested, self.container._counts, self._countsAllocatedInitial)
			
			# Record unpartitioned counts for later merging
			self._countsUnallocated = self.container._counts - np.sum(self._countsAllocatedInitial, axis = -1)

			self._countsAllocatedFinal[:] = self._countsAllocatedInitial

		else:
			self._countsUnallocated = self.container._counts


	def merge(self):
		self.container.countsIs(
			self._countsUnallocated + self._countsAllocatedFinal.sum(axis = -1)
			)

	def mass(self, typeKey = None):
		if typeKey is None:
			indexes = np.s_[:]

		else:
			indexes = self._typeIdxs[typeKey]

		return np.dot(
			self._moleculeMass[indexes],
			self.container._counts[indexes]
			)


	def pytablesCreate(self, h5file, expectedRows):
		countsShape = self.container._counts.shape
		partitionsShape = self._countsRequested.shape

		# Columns
		d = {
			"time": tables.Int64Col(),
			"counts":tables.UInt64Col(countsShape),
			"countsRequested":tables.UInt64Col(partitionsShape),
			"countsAllocatedInitial":tables.UInt64Col(partitionsShape),
			"countsAllocatedFinal":tables.UInt64Col(partitionsShape),
			"countsUnallocated":tables.UInt64Col(countsShape),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self.meta["id"],
			d,
			title = self.meta["name"],
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)
	
		groupNames = h5file.create_group(h5file.root,
			'names', 'Molecule, compartment, and process names')

		h5file.create_array(groupNames, 'moleculeIDs', [str(s) for s in self._moleculeIDs]) # pytables doesn't support unicode
		h5file.create_array(groupNames, 'compartmentIDs', [str(s) for s in self._compartmentIDs])

		groupIdxs = h5file.create_group(h5file.root,
			'indexes', 'Indexes for various groups of molecules')

		for type_, indexes in self._typeIdxs.viewitems():
			h5file.create_array(groupIdxs, type_, indexes)


	def pytablesAppend(self, h5file):
		simTime = self.time.value

		t = h5file.get_node("/", self.meta["id"])
		entry = t.row

		entry["time"] = simTime
		entry['counts'] = self.container._counts
		entry['countsRequested'] = self._countsRequested
		entry['countsAllocatedInitial'] = self._countsAllocatedInitial
		entry['countsAllocatedFinal'] = self._countsAllocatedFinal
		entry['countsUnallocated'] = self._countsUnallocated
		
		entry.append()

		t.flush()


	def pytablesLoad(self, h5file, timePoint):
		entry = h5file.get_node('/', self.meta['id'])[timePoint]

		self.container.countsIs(entry['counts'])
		
		if self._nProcesses:
			self._countsRequested[:] = entry['countsRequested']
			self._countsAllocatedInitial[:] = entry['countsAllocatedInitial']
			self._countsAllocatedFinal[:] = entry['countsAllocatedFinal']
			self._countsUnallocated[:] = entry['countsUnallocated']


def calculatePartition(isRequestAbsolute, countsBulkRequested, countsBulk, countsBulkPartitioned):
	requestsAbsolute = np.sum(countsBulkRequested[..., isRequestAbsolute], axis = -1)
	requestsRelative = np.sum(countsBulkRequested[..., ~isRequestAbsolute], axis = -1)

	# TODO: Remove the warnings filter or move it elsewhere
	# there may also be a way to avoid these warnings by only evaluating 
	# division "sparsely", which should be faster anyway - JM
	oldSettings = np.seterr(invalid = 'ignore', divide = 'ignore') # Ignore divides-by-zero errors

	scaleAbsolute = np.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
		np.minimum(1, # Restrict requests to at most 100% (absolute requests can do strange things)
			np.minimum(countsBulk, requestsAbsolute) / requestsAbsolute) # Divide requests amongst partitions proportionally
		)

	scaleRelative = np.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
		np.maximum(0, countsBulk - requestsAbsolute) / requestsRelative # Divide remaining requests amongst partitions proportionally
		)

	scaleRelative[requestsRelative == 0] = 0 # nan handling?

	np.seterr(**oldSettings) # Restore error handling to the previous state

	# Compute allocations and assign counts to the partitions
	for iPartition in range(countsBulkPartitioned.shape[-1]):
		scale = scaleAbsolute if isRequestAbsolute[iPartition] else scaleRelative
		allocation = np.floor(countsBulkRequested[..., iPartition] * scale)
		countsBulkPartitioned[..., iPartition] = allocation

class BulkMoleculesViewBase(wholecell.views.view.View):
	_stateID = 'BulkMolecules'

	def _counts(self):
		return self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex].copy()


	def _countsIs(self, values):
		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] = values


	def _countsInc(self, values):
		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] += values


	def _countsDec(self, values):
		self._state._countsAllocatedFinal[self.containerIndexes, self._processIndex] -= values


class BulkMoleculesView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculesView, self).__init__(*args, **kwargs)

		# State references
		self.containerIndexes = self._state.container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def counts(self):
		return self._counts()


	def countsIs(self, values):
		self._countsIs(values)


	def countsInc(self, values):
		self._countsInc(values)


	def countsDec(self, values):
		self._countsDec(values)


class BulkMoleculeView(BulkMoleculesViewBase):
	def __init__(self, *args, **kwargs):
		super(BulkMoleculeView, self).__init__(*args, **kwargs)

		# State references
		self.containerIndexes = self._state.container._namesToIndexes((self._query,))


	def count(self):
		return self._counts()


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)

FEIST_CORE_VALS = np.array([ # TODO: This needs to go in the KB
	0.513689, 0.295792, 0.241055, 0.241055, 0.091580, 0.263160, 0.263160, 0.612638, 0.094738, 0.290529,
	0.450531, 0.343161, 0.153686, 0.185265, 0.221055, 0.215792, 0.253687, 0.056843, 0.137896, 0.423162,
	0.026166, 0.027017, 0.027017, 0.026166, 0.133508, 0.215096, 0.144104, 0.174831, 0.013894, 0.019456,
	0.063814, 0.075214, 0.177645, 0.011843, 0.007895, 0.004737, 0.007106, 0.007106, 0.003158, 0.003158,
	0.003158, 0.003158, 0.003158, 0.004737, 0.003948, 0.003948, 0.000576, 0.001831, 0.000447, 0.000223,
	0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000055, 0.000223, 0.000223,
	0.000223		# mmol/gDCW (supp info 3, "biomass_core", column G)
	]) # TOKB

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
