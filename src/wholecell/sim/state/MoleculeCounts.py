#!/usr/bin/env python

"""
MoleculeCounts.py

State which represents for a class of molecules the bulk copy numbers as an 
array and unique instances in a series of lists sharing a common index.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import re

import numpy
import wholecell.sim.state.State as wcState
import wholecell.sim.state.Partition as wcPartition

DEFAULT_FORM = ':mature' # TODO: reconcile "forms" concept

# TODO: get from KB
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
		"RRLA-RRNA:mature[c]", "RRLB-RRNA:mature[c]", "RRLC-RRNA:mature[c]", "RRLD-RRNA:mature[c]", "RRLE-RRNA:mature[c]", "RRLG-RRNA:mature[c]", "RRLH-RRNA:mature[c]",
		"RRSA-RRNA:mature[c]", "RRSB-RRNA:mature[c]", "RRSC-RRNA:mature[c]", "RRSD-RRNA:mature[c]", "RRSE-RRNA:mature[c]", "RRSG-RRNA:mature[c]", "RRSH-RRNA:mature[c]",
		"RRFA-RRNA:mature[c]", "RRFB-RRNA:mature[c]", "RRFC-RRNA:mature[c]", "RRFD-RRNA:mature[c]", "RRFE-RRNA:mature[c]", "RRFF-RRNA:mature[c]", "RRFG-RRNA:mature[c]", "RRFH-RRNA:mature[c]"
		],
	'rRna23Ss':[
		"RRLA-RRNA:mature[c]", "RRLB-RRNA:mature[c]", "RRLC-RRNA:mature[c]", "RRLD-RRNA:mature[c]", "RRLE-RRNA:mature[c]", "RRLG-RRNA:mature[c]", "RRLH-RRNA:mature[c]",
		],
	'rRna16ss':[
		"RRSA-RRNA:mature[c]", "RRSB-RRNA:mature[c]", "RRSC-RRNA:mature[c]", "RRSD-RRNA:mature[c]", "RRSE-RRNA:mature[c]", "RRSG-RRNA:mature[c]", "RRSH-RRNA:mature[c]",
		],
	'rRna5Ss':[
		"RRFA-RRNA:mature[c]", "RRFB-RRNA:mature[c]", "RRFC-RRNA:mature[c]", "RRFD-RRNA:mature[c]", "RRFE-RRNA:mature[c]", "RRFF-RRNA:mature[c]", "RRFG-RRNA:mature[c]", "RRFH-RRNA:mature[c]"
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
	}

# TODO: make most of these classes _private
# TODO: make a base class for bulk counts?
class MoleculeCountsBase(object):
	'''
	MoleculeCountsBase

	Base object for the MoleculeCounts partition and state.  Manages indexing 
	of molecules to support unique instances and bulk quantities.
	'''

	_nMolIDs = None
	_nCompartments = None

	_countsBulk = None
	# _countsUnique = None

	_molecules = None

	_molIDIndex = None
	_compartmentIndex = None

	_molIDs = None
	_compartments = None

	_uniqueDict = None


	def molecule(self, wid, comp):
		wid = wid.rstrip(DEFAULT_FORM)

		if (wid, comp) not in self._molecules:
			self._molecules[wid, comp] = _Molecule(self, self._molIDIndex[wid], self._compartmentIndex[comp], wid)

		return self._molecules[wid, comp]


	def countsBulk(self, ids = None):
		return self.countsBulkViewNew(ids).countsBulk()


	def countsBulkIs(self, counts, ids = None):
		return self.countsBulkViewNew(ids).countsBulkIs(counts)


	def countsBulkInc(self, counts, ids = None):
		return self.countsBulkViewNew(ids).countsBulkInc(counts)


	def _getIndices(self, ids):
		# TODO: better support for form values
		# TODO: use internally in place of self._molIDs?
		# TODO: better support for compartments

		molecules = []
		compartments = []

		for id_ in ids:
			match = re.match("^(?P<molecule>[^:\[\]]+)(?P<form>:[^:\[\]]+)*(?P<compartment>\[[^:\[\]]+\])*$", id_)

			if match is None:
				raise Exception('Invalid ID: {}'.format(id_))

			if match.group('form') is not None and match.group('form') != DEFAULT_FORM:
				#raise NotImplementedError()

				molecules.append(match.group('molecule') + match.group('form'))

			else:
				molecules.append(match.group("molecule"))

			if match.group("compartment") is None:
				compartments.append(self._compartments[0])

			else:
				compartments.append(match.group("compartment")[1])

		try:
			molIdxs = numpy.array([self._molIDIndex[m] for m in molecules])

		except ValueError:
			raise Exception('Invalid molecule: {}'.format(m))

		try:
			compIdxs = numpy.array([self._compartmentIndex[c] for c in compartments])

		except ValueError:
			raise Exception('Invalid compartment: {}'.format(c))

		idxs = numpy.ravel_multi_index(
			numpy.array([molIdxs, compIdxs]),
			(len(self._molIDs), len(self._compartments))
			)

		return idxs, molIdxs, compIdxs


	def _getIndex(self, id_):
		return [values[0] for values in self._getIndices((id_,))]


	def countsBulkViewNew(self, ids = None):
		if ids is None:
			return CountsBulkView(self)

		else:
			return CountsBulkView(self, self._getIndices(ids)[1:])


class CountsBulkView(object):
	'''
	CountsBulkView

	A "view" into a MoleculeCountsBase subclass's bulk counts.  This allows for 
	easy caching of access to and mutation of bulk quantities of molecules, 
	which is helpful for certain calculations.
	'''

	_parent = None
	_indices = None

	def __init__(self, parent, indices = None):
		self._parent = parent

		self._indices = indices

	def countsBulk(self):
		if self._indices is None:
			return self._parent._countsBulk

		else:
			return self._parent._countsBulk[self._indices]

	def countsBulkIs(self, counts):
		if self._indices is None:
			if type(counts) == numpy.ndarray and counts.ndim == 1:
				counts = counts[:, numpy.newaxis] # fixes broadcasting from 1D arrays

			self._parent._countsBulk[:] = counts

			# TODO: determine if there is a better solution for this
			# TODO: consider filing a report for NumPy's broadcasting rules

		else:
			self._parent._countsBulk[self._indices] = counts


	def countsBulkInc(self, counts):
		if self._indices is None:
			self._parent._countsBulk += counts

		else:
			self._parent._countsBulk[self._indices] += counts


class MoleculeCounts(wcState.State, MoleculeCountsBase):
	'''
	MoleculeCounts

	State for MoleculeCounts.  Adds support for partitioning and 
	initialization.
	'''

	compartments = [ # TODO: move to KB
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

	_typeIdxs = None

	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "MoleculeCounts",
			"name": "Molecules Container",
			"dynamics": ["_countsBulk"],#, "_countsUnique", "_uniqueDict"],
			"units": {
				"_countsBulk"   : "molecules",
				# "_countsUnique" : "molecules",
				# "_uniqueDict"	: "molecules"
				}
			}

		self._countsUnique = None # Counts of molecules with unique properties

		# Molecule mass
		self._molMass = None    # Mass of a single bulk molecule
		self._massSingle = None # Mass of a single bulk molecule, for each compartment
		self._dmass = None      # Deviation from typical mass due to some unique property

		# Object references
		self._molecules = {}    # Molecule objects
		self._uniqueDict = None # Record of unique attributes TODO: verify function, rename?

		# Partitioning
		self._countsBulkRequested = None
		self._countsBulkPartitioned = None
		self._countsBulkUnpartitioned = None

		self.partitionClass = MoleculeCountsPartition

		# Reference attributes
		self._typeIdxs = {} 

		super(MoleculeCounts, self).__init__(*args, **kwargs)


	def initialize_old(self, kb):
		self._nMolIDs = len(kb.molecules)
		self._nCompartments = len(kb.compartments)

		# TODO: make this section legible, rewrite as a list comprehension
		self._uniqueDict = []

		for mol in kb.molecules:
			if mol["uniqueAttrs"] is not None:
				self._uniqueDict.append([dict(dict(zip(mol["uniqueAttrs"] + ["objects"], [[] for x in xrange(len(mol["uniqueAttrs"]) + 1)]))) for x in kb.compartments])
			
			else:
				self._uniqueDict.append([{} for x in kb.compartments])

		self._molIDs = [molecule["id"] for molecule in kb.molecules]
		self._compartments = [compartment["id"] for compartment in kb.compartments]

		self._molMass = numpy.array([molecule["mass"] for molecule in kb.molecules], float)

		self._molIDIndex = {wid:i for i, wid in enumerate(self._molIDs)}
		self._compartmentIndex = {c:i for i, c in enumerate(self._compartments)}

		# MoleculeCounts expects the knowledge base to pass metabolites, 
		# proteins, RNAs; not just "molecules"

		# Also expects "forms" (mature vs. nascent) for proteins and RNAs
		# Type assignments (metabolite/rna/protein)

		# References to simulation states (complexation, trans/trans)

		# Track lumped indices i.e. ntps

		# Core biomass function

		# Protein and RNA primary structure data

		# Complexes...?

		# Prefered localization

		# Media and biomass concentrations


	def initialize(self, sim, kb):
		super(MoleculeCounts, self).initialize(sim, kb)

		# Because I'm not sure how we want to handle forms, and the names won't
		# hash properly without them, I'm combining IDs and form values

		# Molecules
		self._molIDs = []

		self._molIDs += [x['id'] for x in kb.metabolites]
		self._molIDs += [x['id'] + ':nascent' for x in kb.rnas]
		self._molIDs += [x['id'] for x in kb.rnas]
		self._molIDs += [x['id'] + ':nascent' for x in kb.proteins]
		self._molIDs += [x['id'] for x in kb.proteins]

		self._molIDIndex = {wid:i for i, wid in enumerate(self._molIDs)}

		self._nMolIDs = len(self._molIDs)

		# Compartments
		self._compartments = [x['id'] for x in self.compartments]
		# self._compartments = [x['id'] for x in kb.compartments]

		self._compartmentIndex = {c:i for i, c in enumerate(self._compartments)}

		self._nCompartments = len(self._compartments)

		# Masses
		molMass = []

		molMass += [x['mw7.2'] for x in kb.metabolites]
		molMass += [x['mw'] for x in kb.rnas]*2
		molMass += [x['mw'] for x in kb.proteins]*2

		self._molMass = numpy.array(molMass, float)

		self._typeIdxs.update({
			'metabolites':numpy.arange(len(kb.metabolites)),
			'rnas':numpy.arange(2*len(kb.rnas))+len(kb.metabolites),
			'proteins':numpy.arange(2*len(kb.proteins))+len(kb.metabolites)+2*len(kb.rnas)
			})

		self._typeIdxs['water'] = self._molIDs.index('H2O')

		# Unique instances
		self._uniqueDict = []

		# TODO: add unique attributes to KB and test
		# TODO: figure out what this block is really doing/refactor
		# for mol in kb.molecules:
		# 	if mol["uniqueAttrs"] is not None:
		# 		self._uniqueDict.append([dict(dict(zip(mol["uniqueAttrs"] + ["objects"], [[] for x in xrange(len(mol["uniqueAttrs"]) + 1)]))) for x in kb.compartments])

		# 	else:
		# 		self._uniqueDict.append([{} for x in kb.compartments])

		for mol in self._molIDs:
			self._uniqueDict.append([{} for x in self._compartments])


	def calcInitialConditions(self):
		self._countsBulk[:] = 0


	def allocate(self):
		super(MoleculeCounts, self).allocate() # Allocates partitions

		self._countsBulk = numpy.zeros((self._nMolIDs, self._nCompartments), float)
		self._massSingle = numpy.tile(self._molMass, [self._nCompartments, 1]).transpose() # Repeat for each compartment

		self._countsUnique = numpy.zeros_like(self._countsBulk)
		self._dmass = numpy.zeros_like(self._countsBulk)

		self._countsBulkRequested = numpy.zeros_like(self._countsBulk)
		self._countsBulkPartitioned = numpy.zeros((self._nMolIDs, self._nCompartments, len(self.partitions)))
		self._countsBulkUnpartitioned = numpy.zeros_like(self._countsBulk)


	# Partitioning

	def addPartition(self, process, reqMols, reqFunc, isReqAbs = False):
		# TODO: warning for adding a partition after allocation
		partition = super(MoleculeCounts, self).addPartition(process)
		
		# TODO: move this to partition constructor?
		partition.reqFunc = reqFunc
		partition.isReqAbs = isReqAbs

		mapping, iMolecule = self._getIndices(reqMols)[:2]

		if len(set(mapping)) < len(mapping):
			raise Exception('Partition request cannot contain duplicate IDs')

		partition.mapping = mapping

		partition._molIDs = [self._molIDs[i] for i in iMolecule]
		partition._molIDIndex = {wid:i for i, wid in enumerate(partition._molIDs)}

		# TODO: determine how compartments should be handled here...
		partition._compartments = ["merged"] # "merged"
		partition._compartmentIndex = {"merged":0} # "merged"

		partition._nMolIDs = len(partition._molIDs)
		partition._nCompartments = len(partition._compartments)

		return partition


	def prepartition(self):
		# Clear out the existing partitions in preparation for the requests
		for partition in self.partitions:
			partition.countsBulkIs(0)


	def partition(self):
		if self.partitions:
			# TODO: partitioning of unique instances (for both specific and nonspecific requests)
			requestsShape = (
				self._countsBulk.shape[0],
				self._countsBulk.shape[1],
				len(self.partitions)
				)

			requests = numpy.zeros(requestsShape)

			# Calculate and store requests
			for iPartition, partition in enumerate(self.partitions):
				# Call request function and record requests
				requests[..., iPartition].flat[partition.mapping] = numpy.maximum(0, partition.request().flatten())

			isRequestAbsolute = numpy.array([x.isReqAbs for x in self.partitions], bool)
			requestsAbsolute = numpy.sum(requests[..., isRequestAbsolute], axis = 2)
			requestsRelative = numpy.sum(requests[..., ~isRequestAbsolute], axis = 2)

			self._countsBulkRequested = numpy.sum(requests, axis = 2)

			# TODO: Remove the warnings filter or move it elsewhere
			# there may also be a way to avoid these warnings by only evaluating 
			# division "sparsely", which should be faster anyway - JM
			oldSettings = numpy.seterr(invalid = 'ignore', divide = 'ignore') # Ignore divides-by-zero errors

			scaleAbsolute = numpy.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
				numpy.minimum(1, # Restrict requests to at most 100% (absolute requests can do strange things)
					numpy.minimum(self._countsBulk, requestsAbsolute) / requestsAbsolute) # Divide requests amongst partitions proportionally
				)

			scaleRelative = numpy.fmax(0, # Restrict requests to at least 0% (fmax replaces nan's)
				numpy.maximum(0, self._countsBulk - requestsAbsolute) / requestsRelative # Divide remaining requests amongst partitions proportionally
				)

			scaleRelative[requestsRelative == 0] = 0 # nan handling?

			numpy.seterr(**oldSettings) # Restore error handling to the previous state

			# Compute allocations and assign counts to the partitions
			for iPartition, partition in enumerate(self.partitions):
				scale = scaleAbsolute if partition.isReqAbs else scaleRelative

				allocation = numpy.floor(requests[..., iPartition] * scale)

				self._countsBulkPartitioned[..., iPartition] = allocation
				partition.countsBulkIs(
					allocation.flat[partition.mapping]
					)
			
			# Record unpartitioned counts for later merging
			self._countsBulkUnpartitioned = self._countsBulk - numpy.sum(self._countsBulkPartitioned, axis = 2)

		else:
			self._countsBulkUnpartitioned = self._countsBulk


	def merge(self):
		self._countsBulk = self._countsBulkUnpartitioned

		for partition in self.partitions:
			self._countsBulk.flat[partition.mapping] += partition.countsBulk().flatten()


	def massAll(self, typeKey = None):
		if typeKey is None:
			return numpy.dot(self._molMass, self._countsBulk)

		else:
			return numpy.dot(
				self._molMass[self._typeIdxs[typeKey]],
				self._countsBulk[self._typeIdxs[typeKey], :]
				)

		# else:
		# 	return numpy.dot(self._molMass[])
		# TODO: unique counts & dMass


class MoleculeCountsPartition(wcPartition.Partition, MoleculeCountsBase):
	'''
	MoleculeCountsPartition

	Partition for MoleculeCounts.  Acts mostly as a container class, with some
	methods for indexing and creating views.  Partitions act as containers for 
	requests prior to partitioning, as well as containers for the counts that
	are ultimately partitioned to the state.
	'''

	mapping = None

	def __init__(self, *args, **kwargs):
		super(MoleculeCountsPartition, self).__init__(*args, **kwargs)
		
		# hack; uniques instances are unimplemented
		self._molecules = {}
		self._uniqueDict = self.state()._uniqueDict


	def allocate(self):
		self._countsBulk = numpy.zeros((self._nMolIDs, self._nCompartments), float)


	def request(self):
		self.reqFunc(self)

		return self._countsBulk


def _uniqueInit(self, uniqueIdx):
	# Default initialization method for _Molecule objects
	self._uniqueIdx = uniqueIdx


def _makeGetter(attr):
	# Used to construct attribute getters for new _Molecules
	def attrGetter(self):
		molCont, uniqueIdx = self._container, self._uniqueIdx
		molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

		return molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx]

	return attrGetter


def _makeSetter(attr):
	# Used to construct attribute setters for new _Molecule objects
	def attrSetter(self, newVal):
		molCont, uniqueIdx = self._container, self._uniqueIdx
		molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

		molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx] = newVal

	return attrSetter


class _Molecule(object):
	uniqueClassRegistry = {}
	def __init__(self, container, rowIdx, colIdx, wid):
		self._container = container # Parent MoleculeCounts object
		self._rowIdx = rowIdx
		self._colIdx = colIdx
		self._wid = wid

		if len(self._container._uniqueDict[self._rowIdx][self._colIdx]) == 0:
			# Molecule has no attributes
			pass

		elif self._wid in self.uniqueClassRegistry:
			# Molecule has been registered as a unique instance; use that class definition
			self._MoleculeUnique = self.uniqueClassRegistry[self._wid]
			self._MoleculeUnique._container = self._container
			self._MoleculeUnique._molRowIdx = self._rowIdx
			self._MoleculeUnique._molColIdx = self._colIdx

			for attr in self._container._uniqueDict[self._rowIdx][self._colIdx]:
				if not hasattr(self._MoleculeUnique, attr):
					setattr(self._MoleculeUnique, attr, _makeGetter(attr))

				if not hasattr(self._MoleculeUnique, attr + "Is"):
					setattr(self._MoleculeUnique, attr + "Is", _makeGetter(attr + "Is"))

		else:
			# Molecule has unique attributes, but isn't registered
			uniqueClassDefDict = {}
			uniqueClassDefDict["_container"] = self._container
			uniqueClassDefDict["_molRowIdx"] = self._rowIdx
			uniqueClassDefDict["_molColIdx"] = self._colIdx

			uniqueClassDefDict["__init__"] = _uniqueInit

			for attr in self._container._uniqueDict[self._rowIdx][self._colIdx]:
				uniqueClassDefDict[attr] = _makeGetter(attr)
				uniqueClassDefDict[attr + "Is"] = _makeSetter(attr)

			self._MoleculeUnique = type("MoleculeUnique", (), uniqueClassDefDict)

	# Interface methods

	def countBulk(self):
		# Returns bulk count of molecule as a float
		return self._container._countsBulk[self._rowIdx, self._colIdx]

	def countBulkIs(self, newVal):
		# Sets bulk count of molecule
		self._container._countsBulk[self._rowIdx, self._colIdx] = newVal

	def countBulkInc(self, incVal):
		# Increments counts bulk by incVal
		self._container._countsBulk[self._rowIdx, self._colIdx] += incVal

	def countBulkDec(self, decVal):
		# Decrements counts bulk by decVal
		self.countBulkInc(-1 * decVal)

	def countUnique(self):
		# Returns unique count of molecule as a float
		return self._container._countsUnique[self._rowIdx, self._colIdx]
	
	def dMassIs(self, newVal):
		# Sets the value for dMass
		self._container._dmass[self._rowIdx, self._colIdx] = newVal

	def dMassInc(self, incVal):
		# Increments dmass by incVal
		self._container._dmass[self._rowIdx, self._colIdx] += incVal

	def dMassDec(self, decVal):
		# Decrements count of dmass by decVal
		self.dMassInc(-1 * decVal)

	def massSingle(self):
		# Returns mass of single object
		return self._container._massSingle[self._rowIdx, self._colIdx]

	def massAll(self):
		# Returns mass of all objects bulk and unique
		return (self.countBulk() + self.countUnique()) * self.massSingle() + self._container._dmass[self._rowIdx, self._colIdx]

	def uniqueNew(self, attrs = None):
		# Creates new unique object with attributes defined by attrs
		# attrs should be in format: {"attr1" : value1, "attr2" : value2, ...}
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]
		
		if not len(uniqueDict):
			raise uniqueException('Attempting to create unique from object with no unique attributes!\n')

		if attrs is not None and len(set(attrs).difference(set(uniqueDict.keys()))): # TODO: change to (set(...) - uniqueDict.viewkeys())
			raise uniqueException('A specified attribute is not included in knoweldge base for this unique object!\n')

		for attr in uniqueDict:
			if attrs is not None and attr in attrs:
				uniqueDict[attr].append(attrs[attr])
			else:
				uniqueDict[attr].append(None)

		uniqueIdx = len(uniqueDict["objects"]) - 1
		uniqueDict["objects"][uniqueIdx] = self._MoleculeUnique(uniqueIdx)
		self._container._countsUnique[self._rowIdx, self._colIdx] += 1

		return uniqueDict["objects"][uniqueIdx]

	def uniquesWithAttrs(self, attrs = None):
		# Returns list of objects with attributes specified in attrs
		# attrs should be in format: {"attr1" : value1, "attr2" : value2, ...}
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]

		if attrs is not None and len(set(attrs).difference(set(uniqueDict.keys()))): # TODO: change to (set(...) - uniqueDict.viewkeys())
			raise uniqueException('A specified attribute is not included in knoweldge base for this unique object!\n')

		if attrs is None or len(attrs) == 0 or (hasattr(attrs, "lower") and attrs.lower() == "all"):
			return uniqueDict["objects"][:]

		L = []
		for i in xrange(len(uniqueDict["objects"])):
			addThis = True
			for attr in attrs:
				if uniqueDict[attr][i] != attrs[attr]:
					addThis = False
			if addThis:
				L.append(uniqueDict["objects"][i])
		return L

	def uniqueDel(self, uniqueObj):
		# Deletes unique object uniqueObj and decrements unique count
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]
		uniqueIdx = uniqueObj._uniqueIdx

		if id(uniqueObj) != id(uniqueDict["objects"][uniqueIdx]):
			raise uniqueException('Unique object to delete does not match row in unique table!\n')

		for i in xrange(uniqueIdx + 1, len(uniqueDict["objects"])):
			uniqueDict["objects"][i]._uniqueIdx -= 1
		for attr in uniqueDict:
			del uniqueDict[attr][uniqueIdx]
		self._container._countsUnique[self._rowIdx, self._colIdx] -= 1			


class uniqueException(Exception):
	'''
	uniqueException
	'''


class MoleculeUniqueMeta(type):
	def __new__(cls, name, bases, attrs):
		attrs.update({"_container": None, "_molRowIdx": None, "_molColIdx": None})

		if '__init__' not in attrs:
			attrs['__init__'] = _uniqueInit

		newClass =  super(MoleculeUniqueMeta, cls).__new__(cls, name, bases, attrs)
		_Molecule.uniqueClassRegistry[attrs["registrationId"]] = newClass
		return newClass
