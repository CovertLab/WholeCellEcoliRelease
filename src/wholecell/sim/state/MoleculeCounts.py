#!/usr/bin/env python

"""
MoleculeCounts

State which represents the copy numbers of a class of molecules as an array

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import numpy
import re

import wholecell.sim.state.State

class MoleculeCounts(wholecell.sim.state.State.State):
	""" MoleculeCounts """

	compartments = [
		{"id": "c", "name": "Cytosol"},
		{"id": "e", "name": "Extracellular space"},
		{"id": "m", "name": "Membrane"}
	]
	cIdx = {"c": 0, "e": 1, "m": 2}

	# Form values
	formVals = {"nascent": 1, "mature": 0}
	typeVals = {"metabolite": 0, "rna": 1, "protein": 2}

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "MoleculeCounts",
			"name": "Molecule Counts",
			"dynamics": ["counts"],
			"units": {"counts": "molecules"}
		}

		# References to processes
		self.complexation = None

		# Constants
		self.ids = None				# Molecule ids
		self.names = None			# Molecule names
		self.forms = None			# Molecule forms (e.g., nascent, mature, etc.)
		self.types = None			# Molecule type (e.g., metabolite, RNA, protein)
		self.mws = None				# Molecular weights (Da)
		self.localizations = None	# Preferred intracellular locations
		self.idx = {}				# Indices over molecules

		self.chamberVolume = 4.3667e-12	# L

		self.metMediaConc = None		# Metabolite media concentration (mM)
		self.metBiomassConc = None		# Metabolite biomass concentration (molecules/cell)
		self.fracInitFreeNMPs = 0.03
		self.fracInitFreeAAs = 0.001

		self.rnaLens = None			# RNA lengths
		self.rnaExp = None			# mature RNA expression

		self.monLens = None			# Protein monomer lengths
		self.monExp = None			# Mature protein monomer expression

		# Dynamical properties
		self.counts = None	# Molecule counts (molecules x compartments)

		# -- Partitioning --
		# Used by parent
		self.partitionedCounts = None		# Molecules partitioned at each time step (molecules x compartments x partitions)
		self.unpartitionedCounts = None		# Molecules not partitioned at each time step (molecules x compartments)

		# Used by children
		self.fullCounts = None				# Full count in parent
		self.mapping = None					# Index mapping between parent, partition
		self.reqFunc = None					# Request function handle
		self.isReqAbs = None				# Requesting absolute copy number (True/False)

		super(MoleculeCounts, self).__init__(*args, **kwargs)

	# Calculate constants
	def initialize(self, sim, kb):
		super(MoleculeCounts, self).initialize(sim, kb)

		self.complexation = sim.getProcess("Complexation")

		# Molecule identities:
		self.ids = \
			[x["id"] for x in kb.metabolites] + \
			[x["id"] for x in kb.rnas] + \
			[x["id"] for x in kb.rnas] + \
			[x["id"] for x in kb.proteins] + \
			[x["id"] for x in kb.proteins]
		
		self.forms = numpy.array(
			[self.formVals["mature"]] * len(kb.metabolites) + \
			[self.formVals["nascent"]] * len(kb.rnas) + \
			[self.formVals["mature"]] * len(kb.rnas) + \
			[self.formVals["nascent"]] * len(kb.proteins) + \
			[self.formVals["mature"]] * len(kb.proteins)
		)
		
		self.types = numpy.array(
			[self.typeVals["metabolite"]] * len(kb.metabolites) + \
			[self.typeVals["rna"]] * len(kb.rnas) + \
			[self.typeVals["rna"]] * len(kb.rnas) + \
			[self.typeVals["protein"]] * len(kb.proteins) + \
			[self.typeVals["protein"]] * len(kb.proteins) 
		)
		
		self.names = \
			[x["name"] for x in kb.metabolites] + \
			[x["name"] for x in kb.rnas] + \
			[x["name"] for x in kb.rnas] + \
			[x["name"] for x in kb.proteins] + \
			[x["name"] for x in kb.proteins]
		
		self.mws = numpy.array(
			[x["mw"] for x in kb.metabolites] + \
			[x["mw"] for x in kb.rnas] + \
			[x["mw"] for x in kb.rnas] + \
			[x["mw"] for x in kb.proteins] + \
			[x["mw"] for x in kb.proteins]
		)

		self.idx["ntps"] = self.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[1]
		self.idx["ndps"] = self.getIndex(["ADP[c]", "CDP[c]", "GDP[c]", "UDP[c]"])[1]
		self.idx["nmps"] = self.getIndex(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])[1]
		self.idx["aas"] = self.getIndex([
			"ALA[c]", "ARG[c]", "ASN[c]", "ASP[c]", "CYS[c]", "GLU[c]", "GLN[c]", "GLY[c]", "HIS[c]", "ILE[c]",  "LEU[c]",
			"LYS[c]", "MET[c]", "PHE[c]", "PRO[c]", "SER[c]", "THR[c]", "TRP[c]", "TYR[c]", "VAL[c]"
			])[1]
		self.idx["h2o"] = self.getIndex("H2O[c]")[1]

		# Localizations
		metLocs = numpy.zeros(len(kb.metabolites))
		metLocs[numpy.array([x["hydrophobic"] for x in kb.metabolites])] = self.cIdx["m"]
		metLocs[numpy.array([not x["hydrophobic"] for x in kb.metabolites])] = self.cIdx["c"]
		protLocs = numpy.array(map(lambda x, lookupTable = self.cIdx: lookupTable[x], [x["compartment"] for x in kb.proteins]))
		self.localizations = numpy.concatenate((
			metLocs,
			numpy.array([self.cIdx["c"]] * len(kb.rnas)),
			numpy.array([self.cIdx["c"]] * len(kb.rnas)),
			numpy.array([self.cIdx["c"]] * len(kb.proteins)),
			protLocs
		))

		# Composition
		self.rnaLens = map(lambda rna: numpy.sum(rna["ntCount"]), kb.rnas)
		self.rnaExp = numpy.array([x["exp"] for x in kb.rnas])
		self.rnaExp /= numpy.sum(self.rnaExp)
		self.idx["nascentRna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["nascent"]: 
													tup[0] == typeVal and tup[1] == formVal, 
												zip(self.types, self.forms)))[0]
		self.idx["matureRna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["mature"]:
													tup[0] == typeVal and tup[1] == formVal,
												zip(self.types, self.forms)))[0]
		self.idx["nascentMrna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["nascent"], validIds = [x["id"] for x in kb.rnas if x["type"] == "mRNA"]:
													tup[0] == typeVal and tup[1] == formVal and tup[2] in validIds,
												zip(self.types, self.forms, self.ids)))[0]
		self.idx["matureMrna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["mature"], validIds = [x["id"] for x in kb.rnas if x["type"] == "mRNA"]:
													tup[0] == typeVal and tup[1] == formVal and tup[2] in validIds,
												zip(self.types, self.forms, self.ids)))[0]

		mons = [x for x in kb.proteins if x["monomer"] == True]
		self.monLens = map(lambda mon: numpy.sum(mon["aaCount"]), mons)
		self.monExp = numpy.array([x["exp"] for x in kb.rnas if x["type"] == "mRNA"])
		self.monExp /= numpy.sum(self.monExp)
		self.idx["matureMonomers"] = self.getIndex([x["id"] + ":mature[" + x["compartment"] + "]" for x in mons])[1]

		cpxs = [x for x in kb.proteins if x["monomer"] == False]
		self.idx["matureComplexes"] = self.getIndex([x["id"] + ":mature[" + x["compartment"] + "]" for x in cpxs])[1]

		self.metMediaConc = numpy.array([x["mediaConc"] for x in kb.metabolites])
		self.metBiomassConc = numpy.array([x["biomassConc"] for x in kb.metabolites])

	# Allocate memory
	def allocate(self):
		super(MoleculeCounts, self).allocate()

		if self.parentState == None:
			self.counts = numpy.zeros((len(self.ids), len(self.compartments)))
			self.partitionedCounts = numpy.zeros((len(self.ids), len(self.compartments), len(self.partitions)))
			self.unpartitionedCounts = numpy.zeros((len(self.ids), len(self.compartments)))
		else:
			self.counts = numpy.zeros(len(self.ids))
			self.fullCounts = numpy.zeros((len(self.ids), len(self.compartments)))
	
	# Calculate initial conditions
	def calcInitialConditions(self):
		from wholecell.util.Constants import Constants

		self.counts[:] = 0

		# Media metabolites
		self.counts[self.types == self.typeVals["metabolite"], self.cIdx["e"]] = numpy.round(self.metMediaConc * self.chamberVolume * Constants.nAvogadro * 1e-3)

		# Biomass metabolites
		metIdx = numpy.where(self.types == self.typeVals["metabolite"])[0]
		self.counts[metIdx, self.localizations[metIdx].astype('int')] = numpy.round(self.metBiomassConc)

		# RNA
		rnaCnts = self.randStream.mnrnd(numpy.round((1 - self.fracInitFreeNMPs) * numpy.sum(self.counts[self.idx["nmps"], self.cIdx["c"]]) / (numpy.dot(self.rnaExp, self.rnaLens))), self.rnaExp)
		self.counts[self.idx["nmps"], self.cIdx["c"]] = numpy.round(self.fracInitFreeNMPs * self.counts[self.idx["nmps"], self.cIdx["c"]])
		self.counts[self.idx["matureRna"], self.localizations[self.idx["matureRna"]].astype('int')] = rnaCnts

		# Protein Monomers
		monCnts = self.randStream.mnrnd(numpy.round((1 - self.fracInitFreeAAs) * numpy.sum(self.counts[self.idx["aas"], self.cIdx["c"]]) / (numpy.dot(self.monExp, self.monLens))), self.monExp)
		self.counts[self.idx["aas"], self.cIdx["c"]] = numpy.round(self.fracInitFreeAAs * self.counts[self.idx["aas"], self.cIdx["c"]])
		self.counts[self.idx["matureMonomers"], self.localizations[self.idx["matureMonomers"]].astype('int')] = monCnts

		# Macromolecular complexation
		c = self.complexation

		c.subunit.counts = self.counts[numpy.unravel_index(c.subunit.mapping, self.counts.shape)]
		c.complex.counts = self.counts[numpy.unravel_index(c.complex.mapping, self.counts.shape)]

		c.subunit.counts, c.complex.counts = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 100)
		c.subunit.counts, c.complex.counts = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 10)
		c.subunit.counts, c.complex.counts = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 1)

		self.counts[numpy.unravel_index(c.subunit.mapping, self.counts.shape)] = c.subunit.counts
		self.counts[numpy.unravel_index(c.complex.mapping, self.counts.shape)] = c.complex.counts

	# -- Partitioning into substates --

	def addPartition(self, process, reqMols, reqFunc, isReqAbs = False):
		# Super class
		partition = super(MoleculeCounts, self).addPartition(process)

		# Clear inherited properties only valid on parent
		partition.compartments = {"id": "merged__", "name": "merged"}			# TODO: Make this a list and fix allocate()
		partition.partitionedCounts = None
		partition.unpartitionedCounts = None
		partition.idx = {}

		# Set child properties for mapping to parent
		iMolFormComp, iMolForm = self.getIndex(reqMols)[0:2]
		if len(set(iMolFormComp)) < len(iMolFormComp):
			raise Exception, "Partition request cannot contain duplicate ids"

		partition.ids = [self.ids[i] for i in iMolForm]
		partition.names = [self.ids[i] for i in iMolForm]
		partition.forms = [self.forms[i] for i in iMolForm]
		partition.types = [self.types[i] for i in iMolForm]
		partition.mws = numpy.array([self.mws[i] for i in iMolForm])

		partition.mapping = iMolFormComp
		partition.reqFunc = reqFunc
		partition.isReqAbs = isReqAbs

		return partition

	# Prepare to partition state among processes
	def prepartition(self):
		for partition in self.partitions:
			partition.fullCounts = self.counts[numpy.unravel_index(partition.mapping, self.counts.shape)]

	# Partition state among processes
	def partition(self):
		# Calculate requests
		touchs = numpy.zeros(self.counts.shape + (len(self.partitions),))
		reqs = numpy.zeros(self.counts.shape + (len(self.partitions),))

		for iPartition in xrange(len(self.partitions)):
			partition = self.partitions[iPartition]

			if not partition.isReqAbs:
				touch = numpy.zeros(self.counts.shape)
				touch[numpy.unravel_index(partition.mapping, touch.shape)] = 1
				touchs[:, :, iPartition] = touch

			req = numpy.zeros(self.counts.shape)
			req[numpy.unravel_index(partition.mapping, req.shape)] = numpy.maximum(0, partition.reqFunc())	# TODO: Fix this line depending on reqFunc's return statement
			reqs[:, :, iPartition] = req

		tmp = numpy.array([x.isReqAbs for x in self.partitions])
		absReqs = numpy.sum(reqs[:, :, tmp], axis = 2)
		relReqs = numpy.sum(reqs[:, :, numpy.logical_not(tmp)], axis = 2)

		absScale = numpy.fmax(0, numpy.minimum(numpy.minimum(self.counts, absReqs) / absReqs, 1))
		relScale = numpy.fmax(0, numpy.maximum(0, self.counts - absReqs) / relReqs)
		relScale[relReqs == 0] = 0

		unReqs = numpy.fmax(0, self.counts - absReqs) / numpy.sum(touchs, axis = 2) * (relReqs == 0)
		unReqs[numpy.sum(touchs, axis = 2) == 0] = 0

		for iPartition in xrange(len(self.partitions)):
			partition = self.partitions[iPartition]

			if partition.isReqAbs:
				scale = absScale
			else:
				scale = relScale

			alloc = numpy.floor(reqs[:, :, iPartition] * scale + unReqs * touchs[:, :, iPartition])
			self.partitionedCounts[:, :, iPartition] = alloc
			partition.counts = alloc[numpy.unravel_index(partition.mapping, alloc.shape)]

		# TODO: Allocate unpartitioned molecules
		self.unpartitionedCounts = self.counts - numpy.sum(self.partitionedCounts, axis = 2)

	# Merge sub-states partitioned to processes
	def merge(self):
		self.counts = self.unpartitionedCounts
		for partition in self.partitions:
			cnt = numpy.zeros(self.counts.shape)
			cnt[numpy.unravel_index(partition.mapping, cnt.shape)] = partition.counts
			self.counts += cnt


	# Get index of molecule by id (id, form, compartment)
	def getIndex(self, ids):
		if self.parentState == None:
			return self.getIndex_parent(ids)
		else:
			mappingList = list(self.mapping)
			try:
				idxs = numpy.array([mappingList.index(x) for x in self.parentState.getIndex(ids)[0]])
			except ValueError, e:
				raise Exception, "Invalid index:\n%s" % x
			compIdxs = numpy.ones(idxs.shape)
			return idxs, idxs, compIdxs

	def getIndex_parent(self, ids):
		if type(ids) == str:
			ids = [ids]

		idForms = []
		comps = []
		for thisId in ids:
			match = re.match("^(?P<molecule>[^:\[\]]+)(?P<form>:[^:\[\]]+)*(?P<compartment>\[[^:\[\]]+\])*$", thisId)
			if match == None:
				raise Exception, "Invalid id: %s" % thisId
			
			if match.group("form") == None:
				idForm = [match.group("molecule"), 0]
			else:
				idForm = [match.group("molecule"), self.formVals[match.group("form")[1:]]]
			idForms.append(idForm)

			if match.group("compartment") == None:
				comps.append(self.compartments[0]["id"])
			else:
				comps.append(match.group("compartment")[1:-1])

		compIds =[x["id"] for x in self.compartments] 
		try:
			compIdxs = numpy.array([compIds.index(x) for x in comps])
		except ValueError, e:
			raise Exception, "Invalid compartment: \n%s" % x

		idFormStr = [x[0] + ":" + str(x[1]) for x in idForms]
		allIds = [x[0] + ":" + str(x[1]) for x in zip(self.ids, self.forms)]
		try:
			idFormIdxs = [allIds.index(x) for x in idFormStr]
		except Exception, e:
			raise Exception, "Invalid id/form: \n%s" % x

		idxs = numpy.ravel_multi_index(numpy.array([idFormIdxs, compIdxs]), (len(self.ids), len(self.compartments)))
		#idxs = numpy.reshape(idxs, ids.shape)

		return idxs, idFormIdxs, compIdxs