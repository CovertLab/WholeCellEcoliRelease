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
		{"id": "i", "name": "Inner membrane"},
		{"id": "j", "name": "Projection"},
		{"id": "l", "name": "Pilus"},
		{"id": "m", "name": "Membrane"},
		{"id": "n", "name": "Nucleoid"},
		{"id": "o", "name": "Outer membrane"},
		{"id": "p", "name": "Periplasm"},
		{"id": "w", "name": "Cell wall"}
	]
	cIdx = {"c": 0, "e": 1, "i": 2, "j": 3, "l": 4, "m": 5, "n": 6, "o": 7, "p": 8, "w": 9}

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
			[x["mw7.2"] for x in kb.metabolites] + \
			[x["mw"] for x in kb.rnas] + \
			[x["mw"] for x in kb.rnas] + \
			[x["mw"] for x in kb.proteins] + \
			[x["mw"] for x in kb.proteins]
		)

		self.idx["ntps"] = self.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[1]
		self.idx["ndps"] = self.getIndex(["ADP[c]", "CDP[c]", "GDP[c]", "UDP[c]"])[1]
		self.idx["nmps"] = self.getIndex(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])[1]
		self.idx["aas"] = self.getIndex([
			"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",  "LEU-L[c]",
			"LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]"
			])[1]
		self.idx["h2o"] = self.getIndex("H2O[c]")[1]
		self.idx["tRnas"] = self.getIndex([
			"gltV-tRNA", "gltT-tRNA", "gltW-tRNA", "gltU-tRNA", "glnU-tRNA", "glnW-tRNA", "glnX-tRNA", "glnV-tRNA", "serT-tRNA", "serW-tRNA", "selC-tRNA",
			"serU-tRNA", "serV-tRNA", "serX-tRNA", "RNA0-302", "lysV-tRNA", "RNA0-303", "RNA0-301", "lysW-tRNA", "lysT-tRNA", "RNA0-306", "metY-tRNA",
			"metW-tRNA", "metZ-tRNA", "metU-tRNA", "metT-tRNA", "thrW-tRNA", "thrV-tRNA", "thrU-tRNA", "thrT-tRNA", "trpT-tRNA", "pheV-tRNA",
			"pheU-tRNA", "glyV-tRNA", "glyY-tRNA", "glyU-tRNA", "glyT-tRNA", "glyX-tRNA", "glyW-tRNA", "proL-tRNA", "proK-tRNA", "proM-tRNA",
			"RNA0-300", "valU-tRNA", "valV-tRNA", "valX-tRNA", "valY-tRNA", "valT-tRNA", "valW-tRNA", "hisR-tRNA", "ileX-tRNA", "RNA0-305",
			"ileV-tRNA", "ileT-tRNA", "ileU-tRNA", "tyrV-tRNA", "tyrU-tRNA", "tyrT-tRNA", "alaX-tRNA", "alaW-tRNA", "alaT-tRNA", "alaV-tRNA",
			"alaU-tRNA", "argY-tRNA", "argZ-tRNA", "argX-tRNA", "argU-tRNA", "argV-tRNA", "argQ-tRNA", "argW-tRNA", "aspV-tRNA", "aspU-tRNA",
			"aspT-tRNA", "RNA0-304", "asnV-tRNA", "asnU-tRNA", "asnT-tRNA", "leuU-tRNA", "leuQ-tRNA", "leuX-tRNA", "leuV-tRNA", "leuT-tRNA",
			"leuZ-tRNA", "leuW-tRNA", "leuP-tRNA", "cysT-tRNA"
			])[1]

		# Localizations
		metLocs = numpy.array([self.cIdx[x["biomassLoc"]] if x["biomassConc"] > 0 else self.cIdx["c"] for x in kb.metabolites])
		# metLocs = -1 * numpy.ones(len(kb.metabolites))
		# metLocs[numpy.array([x["hydrophobic"] for x in kb.metabolites])] = self.cIdx["m"]
		# metLocs[numpy.array([not x["hydrophobic"] for x in kb.metabolites])] = self.cIdx["c"]
		protLocs = numpy.array(map(lambda x, lookupTable = self.cIdx: lookupTable[x], [x["location"] for x in kb.proteins]))
		# protLocs = numpy.array(map(lambda x, lookupTable = self.cIdx: lookupTable[x], [x["compartment"] for x in kb.proteins]))
		self.localizations = numpy.concatenate((
			metLocs,
			numpy.array([self.cIdx["c"]] * len(kb.rnas)),
			numpy.array([self.cIdx["c"]] * len(kb.rnas)),
			numpy.array([self.cIdx["c"]] * len(kb.proteins)),
			protLocs
		))

		# Composition
		self.rnaLens = map(lambda rna: numpy.sum(rna["ntCount"]), kb.rnas)
		self.rnaExp = numpy.array([x["expression"] for x in kb.rnas])
		self.rnaExp /= numpy.sum(self.rnaExp)
		self.idx["nascentRna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["nascent"]: 
													tup[0] == typeVal and tup[1] == formVal, 
												zip(self.types, self.forms)))[0]
		self.idx["matureRna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["mature"]:
													tup[0] == typeVal and tup[1] == formVal,
												zip(self.types, self.forms)))[0]
		self.idx["nascentMrna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["nascent"], validIds = [x["id"] for x in kb.rnas if x["monomerId"] != None]:
													tup[0] == typeVal and tup[1] == formVal and tup[2] in validIds,
												zip(self.types, self.forms, self.ids)))[0]
		self.idx["matureMrna"] = numpy.where(map(lambda tup, typeVal = self.typeVals["rna"], formVal = self.formVals["mature"], validIds = [x["id"] for x in kb.rnas if x["monomerId"] != None]:
													tup[0] == typeVal and tup[1] == formVal and tup[2] in validIds,
												zip(self.types, self.forms, self.ids)))[0]

		mons = [x for x in kb.proteins if x["monomer"] == True and x["modifiedForm"] == False]
		self.monLens = map(lambda mon: numpy.sum(mon["aaCount"]), mons)
		self.monExp = numpy.array([x["expression"] for x in kb.rnas if x["monomerId"] != None])
		self.monExp /= numpy.sum(self.monExp)
		self.idx["matureMonomers"] = self.getIndex([x["id"] + ":mature[" + x["location"] + "]" for x in mons])[1]

		cpxs = [x for x in kb.proteins if x["monomer"] == False]
		self.idx["matureComplexes"] = self.getIndex([x["id"] + ":mature[" + x["location"] + "]" for x in cpxs])[1]

		self.metMediaConc = numpy.array([x["mediaConc"] for x in kb.metabolites])
		self.metBiomassConc = numpy.array([numpy.maximum(x["biomassConc"], 0) if x["id"] != "ATP" else x["biomassConc"] * 0.003 for x in kb.metabolites])

	# Allocate memory
	def allocate(self):
		super(MoleculeCounts, self).allocate()

		if self.parentState == None:
			self.counts = numpy.zeros((len(self.ids), len(self.compartments)))
			self.partitionedCounts = numpy.zeros((len(self.ids), len(self.compartments), len(self.partitions)))
			self.unpartitionedCounts = numpy.zeros((len(self.ids), len(self.compartments)))
			self.requestedCounts = numpy.zeros((len(self.ids), len(self.compartments)))
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
		# TODO: Fix this initialization
		metIdx = numpy.where(self.types == self.typeVals["metabolite"])[0]
		self.counts[metIdx, self.localizations[metIdx].astype('int')] = numpy.round(self.metBiomassConc)

		# RNA
		self.counts[self.getIndex(["RRLA-RRNA", "RRLB-RRNA", "RRLC-RRNA", "RRLD-RRNA", "RRLE-RRNA", "RRLG-RRNA", "RRLH-RRNA"])[1], self.cIdx["c"]] = numpy.round(18700 / 1.36 / 7)
		self.counts[self.getIndex(["RRSA-RRNA", "RRSB-RRNA", "RRSC-RRNA", "RRSD-RRNA", "RRSE-RRNA", "RRSG-RRNA", "RRSH-RRNA"])[1], self.cIdx["c"]] = numpy.round(18700 / 1.36 / 7)
		self.counts[self.getIndex(["RRFB-RRNA", "RRFC-RRNA", "RRFD-RRNA", "RRFE-RRNA", "RRFF-RRNA", "RRFG-RRNA", "RRFH-RRNA"])[1], self.cIdx["c"]] = numpy.round(18700 / 1.36 / 7)
		self.counts[self.idx["tRnas"], self.cIdx["c"]] = numpy.round(205000 / 1.36 / len(self.idx["tRnas"]))
		self.counts[self.idx["matureMrna"], self.cIdx["c"]] = self.randStream.mnrnd(numpy.round(1380 / 1.36), self.rnaExp)[self.idx["matureMrna"] - self.idx["matureRna"][0]]
		# rnaCnts = self.randStream.mnrnd(numpy.round((1 - self.fracInitFreeNMPs) * numpy.sum(self.counts[self.idx["nmps"], self.cIdx["c"]]) / (numpy.dot(self.rnaExp, self.rnaLens))), self.rnaExp)
		# self.counts[self.idx["nmps"], self.cIdx["c"]] = numpy.round(self.fracInitFreeNMPs * self.counts[self.idx["nmps"], self.cIdx["c"]])
		# self.counts[self.idx["matureRna"], self.localizations[self.idx["matureRna"]].astype('int')] = rnaCnts

		# Protein Monomers
		self.counts[self.idx["matureMonomers"], self.localizations[self.idx["matureMonomers"]]] = self.randStream.mnrnd(numpy.round(2360000 / 1.36), self.monExp)
		import ipdb
		ipdb.set_trace()
		# monCnts = self.randStream.mnrnd(numpy.round((1 - self.fracInitFreeAAs) * numpy.sum(self.counts[self.idx["aas"], self.cIdx["c"]]) / (numpy.dot(self.monExp, self.monLens))), self.monExp)
		# self.counts[self.idx["aas"], self.cIdx["c"]] = numpy.round(self.fracInitFreeAAs * self.counts[self.idx["aas"], self.cIdx["c"]])
		# self.counts[self.idx["matureMonomers"], self.localizations[self.idx["matureMonomers"]].astype('int')] = monCnts

		# Macromolecular complexation
		c = self.complexation

		c.subunit.counts = self.counts[numpy.unravel_index(c.subunit.mapping, self.counts.shape)]
		c.complex.counts = self.counts[numpy.unravel_index(c.complex.mapping, self.counts.shape)]

		c.subunit.counts, c.complex.counts = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 1000)
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
		self.requestedCounts = numpy.sum(reqs, axis = 2)
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

	def pytablesCreate(self, h5file, sim):
		import tables

		# TODO: Look into using enumerated data types in PyTables

		# Columns
		d = {
			"time": tables.Int64Col(),
			"id": tables.StringCol(max([len(x) for x in self.ids])),
			"form": tables.StringCol(max([len(x) for x in self.formVals.keys()])),
			"type": tables.StringCol(max([len(x) for x in self.typeVals.keys()])),
			"name": tables.StringCol(max([len(x) for x in self.names if x != None])),
			"compartment": tables.StringCol(max([len(x) for x in self.cIdx.keys()])),
			"counts": tables.Float64Col(),
			"requested": tables.Float64Col()
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(h5file.root, self.meta["id"], d, title = self.meta["name"], filters = tables.Filters(complevel = 9, complib="zlib"))

		# The following lines make querying much faster, but make simulation run-time considerably slower
		# t.cols.id.create_index()
		# t.cols.compartment.create_index()

		# Store units as metadata
		t.attrs.counts_units = self.meta["units"]["counts"]

	def pytablesAppend(self, h5file, sim):
		import tables

		simTime = sim.getState("Time").value
		t = h5file.get_node("/", self.meta["id"])
		entry = t.row

		for i in xrange(self.counts.shape[0]):
			for j in xrange(self.counts.shape[1]):
				entry["time"] = simTime
				entry["id"] = self.ids[i]
				entry["form"] = [key for key,val in self.formVals.iteritems() if val == self.forms[i]][0] # TODO: maybe create another dictionary to avoid this
				entry["type"] = [key for key,val in self.typeVals.iteritems() if val == self.types[i]][0] # TODO: maybe create another dictionary to avoid this
				# entry["name"] = self.names[i].encode("ascii", "ignore")
				entry["compartment"] = [key for key,val in self.cIdx.iteritems() if val == j][0]
				entry["counts"] = self.counts[i, j]
				entry["requested"] = self.requestedCounts[i, j]
				entry.append()

		t.flush()
