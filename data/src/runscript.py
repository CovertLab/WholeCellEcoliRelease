#!/usr/bin/env python
import os
import json
import csv
import re
import numpy as np
import sets
import ipdb
import generateSequencesForGRAVY as gravy
import urllib
import time
from SOAPpy import WSDL
import xml.dom.minidom
import itertools
import copy

t = time.strftime("%Y-%m-%d_%H_%M_%S", time.localtime())

def getEcocycChildren(frameid, inst = []):
	websvcUrl = "http://websvc.biocyc.org/getxml?ECOLI:%s" % frameid
	dom = xml.dom.minidom.parse(urllib.urlopen(websvcUrl))

	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "RNA"):
		tagname = 'RNA'
	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "Protein"):
		tagname = 'Protein'
	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "Compound"):
		tagname = 'Compound'

	for subclass in dom.getElementsByTagName("subclass"):
		rnaSub = subclass.getElementsByTagName('RNA')
		protSub = subclass.getElementsByTagName('Protein')
		cmpdSub = subclass.getElementsByTagName('Compound')
		if len(rnaSub):
			fid = subclass.getElementsByTagName('RNA')[0].getAttribute('frameid')
		elif len(protSub):
			fid = subclass.getElementsByTagName('Protein')[0].getAttribute('frameid')
		elif len(cmpdSub):
			fid = subclass.getElementsByTagName('Compound')[0].getAttribute('frameid')
		getEcocycChildren(fid, inst)
	for instance in dom.getElementsByTagName("instance"):
		rnaSub = instance.getElementsByTagName('RNA')
		protSub = instance.getElementsByTagName('Protein')
		cmpdSub = instance.getElementsByTagName('Compound')
		if len(rnaSub):
			fid = instance.getElementsByTagName('RNA')[0].getAttribute('frameid')
		elif len(protSub):
			fid = instance.getElementsByTagName('Protein')[0].getAttribute('frameid')
		elif len(cmpdSub):
			fid = instance.getElementsByTagName('Compound')[0].getAttribute('frameid')
		inst.append(fid)

def getEcocycParents(frameid, parents):
	websvcUrl = "http://websvc.biocyc.org/getxml?ECOLI:%s" % frameid
	dom = xml.dom.minidom.parse(urllib.urlopen(websvcUrl))

	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "RNA"):
		tagname = 'RNA'
	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "Protein"):
		tagname = 'Protein'
	if any(x for x in dom.childNodes[0].childNodes if x.nodeType == x.ELEMENT_NODE and x.tagName == "Compound"):
		tagname = 'Compound'

	for parentClass in dom.getElementsByTagName("parent"):
		if len(parentClass.getElementsByTagName(tagname)):
			fid = parentClass.getElementsByTagName(tagname)[0].getAttribute('frameid')
			parents.append(fid)
			# TODO: Check - Right now with this commented we are only going one deep!
			#getEcocycParents(fid, parents)

def getEcocycModFormReactions(frameid):
	websvcUrl = "http://websvc.biocyc.org/getxml?ECOLI:%s" % frameid
	dom = xml.dom.minidom.parse(urllib.urlopen(websvcUrl))
	L = []
	for rxn in dom.getElementsByTagName("appears-in-right-side-of"):
		elemRxn = rxn.getElementsByTagName("Reaction")
		if len(elemRxn) > 0:
			fId = elemRxn[0].getAttribute("frameid")
		else:
			raise Exception, "Don't have a reaction frame id."
		L.append(getEcocycReactionStoich(fId))
	for rxn in dom.getElementsByTagName("appears-in-left-side-of"):
		elemRxn = rxn.getElementsByTagName("Reaction")
		if len(elemRxn) > 0:
			fId = elemRxn[0].getAttribute("frameid")
		else:
			raise Exception, "Don't have a reaction frame id."
		L.append(getEcocycReactionStoich(fId))
	return L

def getEcocycReactionStoich(rxn):
	websvcUrl = "http://websvc.biocyc.org/getxml?ECOLI:%s" % rxn
	dom = xml.dom.minidom.parse(urllib.urlopen(websvcUrl))
	L = []

	direction = dom.getElementsByTagName("reaction-direction")
	if len(direction):
		direction = direction[0]
		rxnDir = direction.firstChild.data
	else:
		rxnDir = 'UNKNOWN'

	enzyme = dom.getElementsByTagName("enzyme")
	enz = []
	for e in enzyme:
		enz.append(e.getElementsByTagName("Protein")[0].getAttribute('frameid'))

	for left in dom.getElementsByTagName("left"):
		elemProt = left.getElementsByTagName("Protein")
		elemRna = left.getElementsByTagName("RNA")
		elemCmpnd = left.getElementsByTagName("Compound")
		isclass = False
		if len(elemProt) > 0:
			fId = elemProt[0].getAttribute("frameid")
			if elemProt[0].getAttribute("class") == 'true':
				isclass = True
		elif len(elemRna) > 0:
			fId = elemRna[0].getAttribute("frameid")
			if elemRna[0].getAttribute("class") == 'true':
				isclass = True
		elif len(elemCmpnd) > 0:
			fId = elemCmpnd[0].getAttribute("frameid")
			if elemCmpnd[0].getAttribute("class") == 'true':
				isclass = True
		else:
			raise Exception, "Don't have a frame id for LHS reactant."
		elemCoeff = left.getElementsByTagName("coefficient")
		if len(elemCoeff) > 0:
			coeff = unicode(-1 * float(elemCoeff[0].childNodes[0].data))
		else:
			coeff = u"-1"
		L.append({'id' : fId, 'coeff' : coeff, 'isclass' : isclass})
	for right in dom.getElementsByTagName("right"):
		elemProt = right.getElementsByTagName("Protein")
		elemRna = right.getElementsByTagName("RNA")
		elemCmpnd = right.getElementsByTagName("Compound")
		isclass = False
		if len(elemProt) > 0:
			fId = elemProt[0].getAttribute("frameid")
			if elemProt[0].getAttribute("class") == 'true':
				isclass = True
		elif len(elemRna) > 0:
			fId = elemRna[0].getAttribute("frameid")
			if elemRna[0].getAttribute("class") == 'true':
				isclass = True
		elif len(elemCmpnd) > 0:
			fId = elemCmpnd[0].getAttribute("frameid")
			if elemCmpnd[0].getAttribute("class") == 'true':
				isclass = True
		elif right.firstChild.data == 'a tRNA':
			# Stupid idot forgot to tag these things in the XML files. Doing a manual catch here.
			fId = 'tRNA-Holder'
			isclass = True
		else:
			raise Exception, "Don't have a frame id for RHS reactant."
		elemCoeff = right.getElementsByTagName("coefficient")
		if len(elemCoeff) > 0:
			coeff = unicode(1 * float(elemCoeff[0].childNodes[0].data))
		else:
			coeff = u"1"
		L.append({'id' : fId, 'coeff' : coeff, 'isclass' : isclass})
	return {'id' : rxn, 'direction' : rxnDir, 'components' : L, 'enzyme' : enz}

def main():
	initalizeLog()
	getEcocyc(fetchNew = False)
	parseIntermediateFiles()

	parseGenes()
	parseLocations()
	parseProteinMonomers()
	parseProteinMonomers_modified()
	parseRna()
	parseRNA_modified()
	parseComplexes()
	parseComplexes_modified()
	parseTranscriptionUnits()
	parseMetabolites()
	parseReactions()
	parseEnzymeKinetics()

def initalizeLog():
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log')):
		os.makedirs(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log'))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'wb') as txtfile:
		txtfile.write(t + '\n')

# Ecocyc flat file creation
def generateEcocycFlatFile(query, outFile):
	# TODO: Use ecocyc web services (http://ecocyc.org/web-services.shtml)
	params = urllib.urlencode({"object": "(\"TABULATED\" \"%s\")" % query})
	h = urllib.urlopen("http://ecocyc.org/query", params)
	s = re.sub("\"", "", h.read())
	f = open(outFile, "w")
	f.write(s)
	f.close()

def getEcocyc(fetchNew = False):
	if not fetchNew:
		return

	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'a')

	# Build protein-protein complexes
	bioVeloQuery = 'sort [(x^frame-id, x^name, components, [z^frame-id : z <- x^modified-form]): x <- ecoli^^protein-complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_protein_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build protein-rna complexes
	bioVeloQuery = 'sort [(x^frame-id, x^name, components, [z^frame-id : z <- x^modified-form]): x <- ecoli^^Protein-RNA-Complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_rna_protein_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build protein-small molecule complexes
	bioVeloQuery = 'sort [(x^frame-id, x^name, components, [z^frame-id : z <- x^modified-form]): x <- ecoli^^Protein-Small-Molecule-Complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_protein_small_molecule_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build genes
	bioVeloQuery = 'sort [(G^FRAME-ID, G^NAME, G^LEFT-END-POSITION, G^RIGHT-END-POSITION, G^TRANSCRIPTION-DIRECTION, [(TMP^FRAME-ID, TMP^NAME): TMP <- P], [c^FRAME-ID : c <- G^component-of, c isa transcription-units]) : G<-ECOLI^^All-Genes, P := [TMP: TMP <- G^PRODUCT], #P > 0]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_genes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	bioVeloQuery = 'sort [(G^NAME, G^FRAME-ID, G^SYNONYMS, G^ACCESSION-1, G^ACCESSION-2) : G<-ECOLI^^All-Genes, P := [TMP: TMP <- G^PRODUCT], #P > 0]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_gene_synonyms.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build locations
	bioVeloQuery = '{([(U^NAME, U^FRAME-ID): U<-Z1^LOCATIONS]): Z1<-ECOLI^^Proteins}'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_locations.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build RNA
	bioVeloQuery = 'sort [(Z1^NAME, Z1^FRAME-ID, Z1^GENE, [G^FRAME-ID : G <- Z1^GENE], Z1^UNMODIFIED-FORM, [UM^FRAME-ID : UM <- Z1^UNMODIFIED-FORM], Z1^MODIFIED-FORM, [UZ^FRAME-ID : UZ <- Z1^MODIFIED-FORM]) :  Z1<-ECOLI^^RNAs]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_rna.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build protein monomers
	bioVeloQuery = 'sort [(Z1^NAME, Z1^FRAME-ID, Z1^UNMODIFIED-FORM, [U^FRAME-ID: U<-Z1^UNMODIFIED-FORM], Z1^MODIFIED-FORM, [M^FRAME-ID: M<-Z1^MODIFIED-FORM], Z1^GENE, [G^FRAME-ID: G<-Z1^GENE], [(U^NAME, U^FRAME-ID): U<-Z1^LOCATIONS]): Z1<-ECOLI^^Proteins]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_proteins.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build transcription units
	bioVeloQuery = 'sort [(t^name,t^FRAME-ID,[c^FRAME-ID : c <- t^components, c isa promoters],[c^FRAME-ID : c <- t^components, c isa terminators],[c^FRAME-ID : c <- t^components, c isa all-genes]) : t <- ecoli^^transcription-units]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_transcriptionUnits.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build promoters
	bioVeloQuery = 'sort [(Z1^FRAME-ID, Z1^NAME, Z1^BINDS-SIGMA-FACTOR, Z1^ABSOLUTE-PLUS-1-POS, [c^FRAME-ID : c <- Z1^component-of, c isa transcription-units]) :  Z1<-ECOLI^^Promoters]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_promoters.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	# Build terminators
	bioVeloQuery = 'sort [(Z1^FRAME-ID, Z1^NAME, Z1^LEFT-END-POSITION, Z1^RIGHT-END-POSITION, [c^FRAME-ID : c <- Z1^component-of, c isa transcription-units]) :  Z1<-ECOLI^^Rho-Independent-Terminators]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rhoIndepTerm.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	bioVeloQuery = 'sort [(Z1^FRAME-ID, Z1^NAME, Z1^LEFT-END-POSITION, Z1^RIGHT-END-POSITION, [c^FRAME-ID : c <- Z1^component-of, c isa transcription-units]) :  Z1<-ECOLI^^Rho-Dependent-Terminators]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rhoDepTerm.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)
	writeOut(bioVeloQuery, logFile)

	logFile.close()

# Intermediate file functions
def parseIntermediateFiles():
	# Load and save gene synonym dictionary
	parseGeneSynonymDictionary()
 	parseGeneProductUnmodifiedForm()
	parseRnaTypes()

def parseGeneSynonymDictionary():
	synDict = {}
	synDictFrameId = {}

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_gene_synonyms.csv')) as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			name = re.sub('<[^<]+?>', '', row[0])
			name = name.replace('-','')
			frameId = row[1]
			synRaw = row[2]
			synRaw = synRaw[1:-2]
			synRaw = synRaw.replace('"', '')
			synRaw = synRaw.replace('-', '')
			synList = synRaw.split(' ')
			if synList == ['']:
				synList = []

			if len(row) > 3:
				synList.append(row[3].lower())
			if len(row) > 4:
				synList.append(row[4].lower())

			synDict[name.lower()] = name.lower()
			synDictFrameId[name.lower()] = frameId
			if synList != ['']:
				for syn in synList:
					if syn != '':
						# Make sure names are not a key and a value
						if synDict.has_key(syn.lower()):
							pass
						else:
							synDict[syn.lower()] = name.lower()
							synDictFrameId[syn.lower()] = frameId

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_name_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDict, indent = 4))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_frameId_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDictFrameId, indent = 4))

def parseGeneProductUnmodifiedForm():
	unmodifiedForm = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[3] != '':
				if row[4] != '':
					unmodifiedForm[row[1]] = True
				else:
					unmodifiedForm[row[1]] = False

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_proteins.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[2] != '':
				unmodifiedForm[row[1]] = True
			else:
				unmodifiedForm[row[1]] = False

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_product_unmodifiedForm.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(unmodifiedForm))

def parseRnaTypes():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	rnaType = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[3] != '':
				if row[1].count('RRNA') > 0:
					rnaType[row[1]] = 'rRNA'
				elif (row[1].count('tRNA') > 0 or row[0].count('tRNA')) and not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'tRNA'
				elif not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'miscRNA'

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'rna_types.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(rnaType))

# Parse genes
def parseGenes():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'a')

	# Load unmodified forms of RNA and proteins
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	# Load RNA types
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'rna_types.json'),'rb') as jsonfile:
		rnaType = json.loads(jsonfile.read())
	
	# Load synonym dictionaries
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_name_synonyms.json'),'rb') as jsonfile:
		synDict = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_frameId_synonyms.json'),'rb') as jsonfile:
		synDictFrameId = json.loads(jsonfile.read())

	# Parse basic information, RNA type, and product
	geneDict = {}
	parameters = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')
		for row in csvreader:
			# Check for products. One gene --> one product (even if gene is duplicated in name)
			allProducts = splitBigBracket(row[5])
			for product in allProducts:
				props = splitSmallBracket(product)
				if not unmodifiedForm[props['frameId']]:
					# If there is no unmodified form then it is the base product of a gene
					# Build new gene
					newGene = gene()

					# Add names
					newGene.productFrameId = props['frameId']
					newGene.name = props['description']
					newGene.frameId = row[0]
					newGene.symbol = row[1]

					# Add direction
					newGene.direction = row[4]

					# Add locations
					if row[2] != '':
						if newGene.direction == '+':
							newGene.coordinate = int(row[2])
							newGene.length = int(row[3]) - int(row[2]) + 1
						else:
							newGene.coordinate = int(row[3])
							newGene.length = int(row[3]) - int(row[2]) + 1

						# Pick new gene name for product if gene name is already used
						# for another valid product
						if geneDict.has_key(newGene.frameId):
							count = 0
							while geneDict.has_key(newGene.frameId + '_' + str(count)):
								count += 1
							geneDict[newGene.frameId + '_' + str(count)] = newGene
							newGene.frameId = newGene.frameId + '_' + str(count)
						else:
							geneDict[newGene.frameId] = newGene

						# Add RNA type
						if rnaType.has_key(newGene.productFrameId):
							newGene.type = rnaType[newGene.productFrameId]
						else:
							newGene.type = 'mRNA'

	# Parse splicing information
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'gene_coordinates.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			frameId = row[0]
			coord = row[1]
			if 'join' in coord:
				spliceCoordsRaw = re.findall("(?P<numbers>[0-9]+\.\.[0-9]+)", coord)
				allSplices = []
				for s in spliceCoordsRaw:
					splitString = s.split('..')
					splice = [int(splitString[0]), int(splitString[1])]
					splice.sort()
					splice = tuple(splice)
					allSplices.append(splice)

				allSplices.sort()
				allSplices = tuple(allSplices)
				geneDict[frameId].splices = allSplices

	# Write random sequence modification information manually
	# 'EG11227' has a selenocystine in place of a stop codon
	geneDict['EG11227'].sequenceSubstitution = (196, '*', 'U')
	# 'EG11858' has a selenocystine in place of a stop codon
	geneDict['EG11858'].sequenceSubstitution = (196, '*', 'U')
	# 'EG10285' has a selenocystine in place of a stop codon
	geneDict['EG10285'].sequenceSubstitution = (140, '*', 'U')
	# 'G8205' has a sequence conflict
	geneDict['G8205'].sequenceSubstitution = (39, '*', 'Y')

	# Parse half life information
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','other_parameters.csv')) as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			parameters[row[0]] = {'value' : row[1], 'units' : row[2]}

	rRNAhl = []
	tRNAhl = []
	miscRNAhl = []
	mRNAhl = []

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Bernstein 2002.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		lineCnt = 0
		startRead = False
		for row in csvreader:
			if not startRead:
				lineCnt += 1
				if lineCnt >= 9:
					startRead = True
			elif row[4] != '':
				read_bnum = row[0].lower()
				read_name = row[1].lower()
				same = False
				if read_name == read_bnum:
					same = True
				read_bnum = 'b' + read_bnum[1:]
				if same:
					read_name = 'b' + read_name[1:]
				read_halfLife = row[4]

				name = [read_bnum, read_name]
				halfLife = float(read_halfLife)*60. # seconds

				if synDictFrameId.has_key(read_bnum.lower()):
					geneName = synDictFrameId[read_bnum.lower()]
				elif synDictFrameId.has_key(read_name.lower()):
					geneName = synDictFrameId[read_name.lower()]
				else:
					s = 'Gene half life not found ' + read_name + ' ' + read_bnum
					writeOut(s, logFile)

				geneDict[geneName].halfLife = halfLife

				# Calculate average half lives of each type of gene
				if geneDict[geneName].type == 'rRNA':
					rRNAhl.append(halfLife)
				elif geneDict[geneName].type == 'tRNA':
					tRNAhl.append(halfLife)
				elif geneDict[geneName].type == 'miscRNA':
					miscRNAhl.append(halfLife)
				else:
					mRNAhl.append(halfLife)

		s = 'Half lives found for ' + str(len(mRNAhl)) + ' mRNAs'
		writeOut(s, logFile)
		s = 'Half lives found for ' + str(len(rRNAhl)) + ' rRNAs'
		writeOut(s, logFile)
		s = 'Half lives found for ' + str(len(tRNAhl)) + ' tRNAs'
		writeOut(s, logFile)
		s = 'Half lives found for ' + str(len(miscRNAhl)) + ' miscRNAs'
		writeOut(s, logFile)

		# Average half lives added for genes that half life was not measured
		mrnaAverage = np.around(np.average(mRNAhl),decimals=2)

		if len(rRNAhl) == 0:
			s = 'No rRNA half lives measured. Using parameter.'
			writeOut(s, logFile)
			rrnaAverage = int(parameters['rRNA half life']['value'])*24*60*60 # seconds
		if len(tRNAhl) == 0:
			s = 'No tRNA half lives measured. Using parameter.'
			writeOut(s, logFile)
			trnaAverage = int(parameters['tRNA half life']['value'])*24*60*60 # seconds
		if len(miscRNAhl) == 0:
			s = 'No miscRNA half lives measured. Using parameter.'
			writeOut(s, logFile)
			miscrnaAverage = mrnaAverage

		for geneId in geneDict.iterkeys():
			geneOfInterest = geneDict[geneId]

			if geneOfInterest.halfLife == None:
				if geneOfInterest.type == 'rRNA':
					geneOfInterest.halfLife = rrnaAverage
				if geneOfInterest.type == 'tRNA':
					geneOfInterest.halfLife = trnaAverage
				if geneOfInterest.type == 'miscRNA':
					geneOfInterest.halfLife = mrnaAverage
				if geneOfInterest.type == 'mRNA':
					geneOfInterest.halfLife = mrnaAverage

	# Parse expression information
	expressionDict = {}

	rRNAexp = []
	tRNAexp = []
	miscRNAexp = []
	mRNAexp = []

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Blattner 2005.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		lineCnt = 0
		startRead = False
		for row in csvreader:
			if not startRead:
				lineCnt += 1
				if lineCnt >= 98:
					startRead = True
			else:
				skip = False
				if row[27] == '':
					skip = True

				if not skip:
					bnum = row[27]
					name = row[1]
					glucoseValues = row[2:7]
					glucoseValues = [float(x) for x in glucoseValues]
					expression = np.average(glucoseValues)

				if synDictFrameId.has_key(name.lower()):
					geneName = synDictFrameId[name.lower()]
				elif synDictFrameId.has_key(bnum.lower()):
					geneName = synDictFrameId[bnum.lower()]
				else:
					s = 'Gene expression not found ' + name + ' ' + bnum
					writeOut(s, logFile)


				expressionDict[geneName] = expression

				# Calculat average expression for each type of RNA
				if geneDict[geneName].type == 'rRNA':
					rRNAexp.append(expression)
				elif geneDict[geneName].type == 'tRNA':
					tRNAexp.append(expression)
				elif geneDict[geneName].type == 'miscRNA':
					miscRNAexp.append(expression)
				else:
					mRNAexp.append(expression)

	s = 'Expression found for ' + str(len(mRNAexp)) + ' mRNAs'
	writeOut(s, logFile)
	s = 'Expression found for ' + str(len(rRNAexp)) + ' rRNAs'
	writeOut(s, logFile)
	s = 'Expression found for ' + str(len(tRNAexp)) + ' tRNAs'
	writeOut(s, logFile)
	s = 'Expression found for ' + str(len(miscRNAexp)) + ' miscRNAs'
	writeOut(s, logFile)

	# Average expression for gene types where expression was not measured
	mrnaAverage = np.around(np.average(mRNAexp),decimals=2)

	if len(rRNAexp) == 0:
		s = 'No rRNA expression measured.'
		writeOut(s, logFile)
		rrnaAverage = mrnaAverage
	else:
		rrnaAverage = np.average(rRNAexp)
	if len(tRNAexp) == 0:
		s = 'No tRNA expression measured.'
		writeOut(s, logFile)
		trnaAverage = mrnaAverage
	else:
		trnaAverage = np.average(tRNAexp)
	if len(miscRNAexp) == 0:
		s = 'No miscRNA expression measured.'
		writeOut(s, logFile)
		miscrnaAverage = mrnaAverage
	else:
		miscrnaAverage = np.average(miscRNAexp)


	for geneId in geneDict.iterkeys():
		if not expressionDict.has_key(geneId):
			if geneDict[geneId].type == 'mRNA':
				expressionDict[geneId] = mrnaAverage
			elif geneDict[geneId].type == 'rRNA':
				expressionDict[geneId] = rrnaAverage
			elif geneDict[geneId].type == 'tRNA':
				expressionDict[geneId] = trnaAverage
			elif geneDict[geneId].type == 'miscRNA':
				expressionDict[geneId] = miscrnaAverage

	total = np.sum(expressionDict.values())

	for key in expressionDict.iterkeys():
		geneDict[key].expression = expressionDict[key] / total

	# Calculate GRAVY and save output
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'protein_monomer_gravy.csv')):
		s = 'Calculating gravy for all genes'
		writeOut(s, logFile)
		gravy.main()
	else:
		s = 'Gravy already exists'
		writeOut(s, logFile)

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['ID', 'Name', 'Symbol', 'Type', 'Coordinate', 'Length', 'Direction', 'Expression', 'Half life (s)', 'Product', 'Splices', '(absolute nt position, old, new)', 'Comments'])

		keys = geneDict.keys()
		keys.sort()

		for key in keys:
			g = geneDict[key]
			csvwriter.writerow([g.frameId, g.name, g.symbol, g.type, g.coordinate, g.length, g.direction, "%0.10f" % g.expression, g.halfLife, g.productFrameId, json.dumps(g.splices), json.dumps(g.sequenceSubstitution), g.comments])


	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'genes.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(geneDict.keys()))

	logFile.close()

# Parse Locations
def parseLocations():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'a')

	locationDict = {}
	# Finds unique set of location frameId's in Ecocyc. Creates a dict so that any locaitons
	# read from another file can be translated into the ones used in the model.
	locationList = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_locations.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row != []:
				locations = splitBigBracket(row[0])

				for location in locations:
					param = splitSmallBracket(location)
					locationList.append(param['description'])

	locationSet = sets.Set(locationList)
	for item in locationSet:
		if item == 'CCO-RIBOSOME':
			locationDict[item] = 'CCO-CYTOSOL'
		elif item == 'CCO-MIT-LUM':
			locationDict[item] = 'CCO-CYTOSOL'
		elif item == 'CCO-MIT-MEM':
			locationDict[item] = 'CCO-PM-BAC-NEG'
		elif item == 'CCO-CYTOSKELETON':
			locationDict[item] = 'CCO-CYTOSOL'
		elif item == 'CCO-ENVELOPE':
			locationDict[item] = 'CCO-OUTER-MEM'
		else:
			locationDict[item] = item

	# Location keys to remove
	toRemove = ['CCO-RIBOSOME', 'CCO-MIT-LUM', 'CCO-MIT-MEM', 'CCO-CYTOSKELETON', 'CCO-ENVELOPE']
	# Assign one letter abbreviations to locations
	abbrevDict = {'CCO-BAC-NUCLEOID' 		: 'n',
					'CCO-CELL-PROJECTION' 	: 'j',
					'CCO-CW-BAC-NEG' 		: 'w',
					'CCO-CYTOSOL'			: 'c',
					'CCO-EXTRACELLULAR'		: 'e',
					'CCO-MEMBRANE'			: 'm',
					'CCO-OUTER-MEM' 		: 'o',
					'CCO-PERI-BAC'			: 'p',
					'CCO-PILUS'				: 'l',
					'CCO-PM-BAC-NEG'		: 'i'}

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'locations_equivalent_names.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(locationDict))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		for item in toRemove:
			locationDict.pop(item)
		keys = locationDict.keys()
		keys.sort()
		csvwriter.writerow(['ID', 'Abbreviation'])
		for key in keys:
			csvwriter.writerow([key, abbrevDict[key]])

# Parse protein monomers
def parseProteinMonomers():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'a')

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_name_synonyms.json'),'rb') as jsonfile:
		synDict = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_frameId_synonyms.json'),'rb') as jsonfile:
		synDictFrameId = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'genes.json'),'rb') as jsonfile:
		geneIdList = json.loads(jsonfile.read())
	# TODO: Might not need to use this last dict depending on the datasource
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'locations_equivalent_names.json'),'rb') as jsonfile:
		locationEquivDict = json.loads(jsonfile.read())

	proteinMonomerDict = {}
	geneToProteinMonomerDict = {}
	unknownGenes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_proteins.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			# Check for unmodified forms. If they exist then skip it.
			unmodifiedForms = True
			if row[2] == '':
				unmodifiedForms = False

			# Check for a known associated gene. If one does not exist skip it.
			knownGene = True
			if row[7] == '':
				knownGene = False
			elif row[7][1:-1] not in geneIdList:
				knownGene = False

			if not knownGene:
				unknownGenes.append(row[1])

			if not unmodifiedForms and knownGene:
				pMono = proteinMonomer()

				pMono.frameId = row[1]
				pMono.name = re.sub('<[^<]+?>', '', row[0])
				pMono.gene = row[7][1:-1]

				modifiedForm = row[5][1:-1].split(' ')
				if modifiedForm == ['']:
					modifiedForm = []
				pMono.modifiedForm = modifiedForm

				proteinMonomerDict[pMono.frameId] = pMono
				geneToProteinMonomerDict[pMono.gene] = pMono.frameId

	# Write out protein monomers with unknown genes
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'protein_monomers_unknown_genes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Frame ID', 'Comments'])
		for u in unknownGenes:
			csvwriter.writerow([u])

	# Add location information
	locationSynDict = {'C'	:	'CCO-CYTOSOL',
						'IM':	'CCO-PM-BAC-NEG',
						'OM':	'CCO-OUTER-MEM',
						'P'	:	'CCO-PERI-BAC',
						'F'	:	'CCO-CELL-PROJECTION',
						'BF':	'CCO-CELL-PROJECTION',
						'E'	:	'CCO-EXTRACELLULAR'}

	# Start with experimentally determined locations in E. coli K-12
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Lopez Campistrous 2005.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()

		for row in csvreader:
			name = row[0].lower()
			bnum = row[1].lower()

			if synDictFrameId.has_key(name):
				geneFrameId = synDictFrameId[name]
			elif synDictFrameId.has_key(bnum):
				geneFrameId = synDictFrameId[bnum]
			else:
				print 'Location parsing Lopez Campistrous 2005: No name found for ' + name + ' ' + bnum

			if geneToProteinMonomerDict.has_key(geneFrameId):
				proteinMonomerFrameId = geneToProteinMonomerDict[geneFrameId]
			else:
				s = 'Location parsing Lopez Campistrous 2005: No name found for gene ' + geneFrameId
				writeOut(s, logFile)

			location = row[2]
			if location != '?':
				location = [locationSynDict[location]]
				proteinMonomerDict[proteinMonomerFrameId].comments += 'Location information from Lopez Campistrous 2005.'
			else:
				location = []
			proteinMonomerDict[proteinMonomerFrameId].location = location

	s = 'Locations found in E. coli K-12 for ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))
	writeOut(s, logFile)

	# Fill in more with locaitons computationally inferred in E. coli B
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Han 2011.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		# Skip header
		csvreader.next()
		for row in csvreader:
			hasData = True
			# Checks to make sure that it either has a gene name or a Blattner number
			if row[0] == '' and row[1] == '':
				hasData = False

			if hasData:
				name = row[0].lower()
				bnumList = row[1].lower().split('-')

				found = False
				geneFrameId = None
				if synDictFrameId.has_key(name):
					geneFrameId = synDictFrameId[name]
					found = True
				else:
					for bnum in bnumList:
						if synDictFrameId.has_key(bnum):
							geneFrameId = synDictFrameId[bnum]
							found = True
							break
				if not found:
					s = 'Location parsing Han 2011: No name found for ' + name + ' ' + bnum
					writeOut(s, logFile)

				if geneToProteinMonomerDict.has_key(geneFrameId) and found:
					proteinMonomerFrameId = geneToProteinMonomerDict[geneFrameId]
				elif found:
					s = 'Location parsing Han 2011: Gene found but no corresponding protein monomer for ' + geneFrameId + ' ' + name + ' ' + str(bnumList)
					writeOut(s, logFile)

				location = row[2]
				if proteinMonomerDict[proteinMonomerFrameId].location == [] and found:
					location = [locationSynDict[location]]
					proteinMonomerDict[proteinMonomerFrameId].comments += 'Location information from Han 2011.'
					proteinMonomerDict[proteinMonomerFrameId].location = location

	s = 'Locations found in E. coli K-12 and B for ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))
	writeOut(s, logFile)

	# Fill in rest from Ecocyc that are possible
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_proteins.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			# Check for unmodified forms. If they exist then skip it.
			unmodifiedForms = True
			if row[2] == '':
				unmodifiedForms = False

			# Check for a known associated gene. If one does not exist skip it.
			knownGene = True
			if row[6] == '':
				knownGene = False

			if not unmodifiedForms and knownGene:
				frameId = row[1]
				rawLocations = splitBigBracket(row[8])
				parsedLocations= []
				for loc in rawLocations:
					if loc != '':
						param = splitSmallBracket(loc)
						if param['description'] != []:
							parsedLocations.append(locationEquivDict[param['description']])
				if len(parsedLocations) == 1:
					if proteinMonomerDict.has_key(frameId):
						if proteinMonomerDict[frameId].location == []:
							proteinMonomerDict[frameId].location = parsedLocations
							proteinMonomerDict[frameId].comments += 'Localization generated from unambiguous Ecocyc data.'

	s = 'Unambiguous locations found in Ecocyc for ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))
	writeOut(s, logFile)

	# Fill in the rest using GRAVY calculation
	gravy = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'protein_monomer_gravy.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			gravy[row[0]] = float(row[1])

	for key in proteinMonomerDict.iterkeys():
		if proteinMonomerDict[key].location == []:
			if gravy[key] < 0:
				proteinMonomerDict[key].location = ['CCO-CYTOSOL']
			else:
				proteinMonomerDict[key].location = ['CCO-MEMBRANE']
			proteinMonomerDict[key].comments += 'Location calculated to be either CCO-CYTOSOL or CCO-MEMBRANE based on GRAVY.'

	s = 'Gravy used to fill in rest ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))
	writeOut(s, logFile)

	# Manually set some locations
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'manual_protein_location.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			proteinMonomerDict[row[0]].location = [row[1]]
			proteinMonomerDict[row[0]].comments += 'Location manually re-set.'

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinMonomerDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Name', 'Gene', 'Location', 'Modified form', 'Comments'])
		for key in keys:
			pm = proteinMonomerDict[key]
			csvwriter.writerow([pm.frameId, pm.name, pm.gene, json.dumps(pm.location), json.dumps(pm.modifiedForm), pm.comments])
	logFile.close()

def parseProteinMonomers_modified():
	# Load location abbreviations
	locationAbbrevDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			locationAbbrevDict[row['ID']] = row['Abbreviation']

	# Load conversion between Ecocyc metabolite frame id's and metabolite id's from Feist. This is used for small-molecule/protein complexes.
	metaboliteEcocycToFeistIdConversion = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'ecocyc_to_feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			if row[1] == '+':
				metaboliteEcocycToFeistIdConversion[row[0]] = row[0].lower()
			else:
				metaboliteEcocycToFeistIdConversion[row[0]] = row[1]

	# Build cache of modified form reactions
	rebuild = False
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_monomer_modification_reactions.json')) or rebuild:
		modFormRxn = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'rb') as csvfile:
			dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
			for row in dictreader:
				if len(json.loads(row['Modified form'])):
					for frameId in json.loads(row['Modified form']):
						rxn = getEcocycModFormReactions(frameId)
						if rxn == []:
							print 'No reaction for ' + frameId
						print 'Loaded ' + frameId + ' formation reaction'
						modFormRxn[frameId] = rxn

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_monomer_modification_reactions.json'),'wb') as jsonfile:
			jsonfile.write(json.dumps(modFormRxn, indent = 4))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_monomer_modification_reactions.json'),'rb') as jsonfile:
		modFormRxn = json.loads(jsonfile.read())

	proteinMonomerDict_modified = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			if len(json.loads(row['Modified form'])):
				for frameId in json.loads(row['Modified form']):
					pm = proteinMonomer()
					pm.frameId = frameId
					pm.unmodifiedForm = row['Frame ID']
					pm.location = json.loads(row['Location'])

					rxn_raw = modFormRxn[pm.frameId]

					if rxn_raw != []:
						production_reaction = []

						for rxn in rxn_raw:
							for rxn_species in rxn['components']:
								if int(float(rxn_species['coeff'])) > 0 and rxn_species['id'] == pm.frameId and rxn['direction'] == 'LEFT-TO-RIGHT':
									production_reaction.append(rxn)
								if int(float(rxn_species['coeff'])) < 0 and rxn_species['id'] == pm.frameId and rxn['direction'] == 'RIGHT-TO-LEFT':
									production_reaction.append(rxn)
								elif rxn['direction'] == 'UNKNOWN' and rxn_species['id'] == pm.frameId:
									production_reaction.append(rxn)

						for rxn in production_reaction:
							fillInReaction(pm, rxn, locationAbbrevDict, metaboliteEcocycToFeistIdConversion)

					proteinMonomerDict_modified[pm.frameId] = pm

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers_modified.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinMonomerDict_modified.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Unmodified Form', 'Location', 'Reaction ID', 'Reaction enzyme', 'Reaction', 'Comments'])
		for key in keys:
			pm = proteinMonomerDict_modified[key]
			csvwriter.writerow([pm.frameId, pm.unmodifiedForm, json.dumps(pm.location), json.dumps(pm.reactionId), json.dumps(pm.reactionEnzymes), json.dumps(pm.reaction), pm.comments])

def fillInReaction(obj, rxn, locationAbbrevDict, metaboliteEcocycToFeistIdConversion):
	obj.reactionId.append(rxn['id'])
	obj.reaction.append('')
	obj.reactionEnzymes.append(rxn['enzyme'])

	reactants = []
	products = []
	for species in rxn['components']:
		if int(float(species['coeff'])) < 0:
			reactants.append(species)
		else:
			products.append(species)

	obj.reaction[-1] += '[' + locationAbbrevDict[obj.location[0]] + ']: '

	for i,r in enumerate(reactants):
		rst = abs(int(float(r['coeff'])))
		rid = str(r['id'])
		if metaboliteEcocycToFeistIdConversion.has_key(rid):
			rid = metaboliteEcocycToFeistIdConversion[rid]

		if abs(rst) > 1:
			obj.reaction[-1] += '(' + str(rst) + ') '
		obj.reaction[-1] += rid
		if i < len(reactants) - 1:
			obj.reaction[-1] += ' + '

	obj.reaction[-1] += ' ==> '

	for i,p in enumerate(products):
		pst = abs(int(float(p['coeff'])))
		pid = str(p['id'])
		if metaboliteEcocycToFeistIdConversion.has_key(pid):
			pid = metaboliteEcocycToFeistIdConversion[pid]

		if pst > 1:
			obj.reaction[-1] += '(' + str(pst) + ') '
		obj.reaction[-1] += pid
		if i < len(products) - 1:
			obj.reaction[-1] += ' + '

# Parse RNA
def parseRna():
	rnaDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			# Check for unmodified forms. If they exist then skip it.
			unmodifiedForms = True
			if row[4] == '':
				unmodifiedForms = False

			# Check for a known associated gene. If one does not exist skip it.
			knownGene = True
			if row[3] == '':
				knownGene = False

			if not unmodifiedForms and knownGene:
				newRna = rna()

				newRna.frameId = row[1]
				newRna.name = re.sub('<[^<]+?>', '', row[0])
				newRna.gene = row[3][1:-1]

				modifiedForm = row[7][1:-1].split(' ')
				if modifiedForm == ['']:
					modifiedForm = []
				newRna.modifiedForm = modifiedForm

				newRna.location = ['CCO-CYTOSOL']
				newRna.comments += 'Assumed in CCO-CYTOSOL'

				rnaDict[newRna.frameId] = newRna

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = rnaDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Name', 'Gene', 'Location', 'Modified form', 'Comments'])
		for key in keys:
			rnaToPrint = rnaDict[key]
			csvwriter.writerow([rnaToPrint.frameId, rnaToPrint.name, rnaToPrint.gene, json.dumps(rnaToPrint.location), json.dumps(rnaToPrint.modifiedForm), rnaToPrint.comments])

def parseRNA_modified():
	# Load location abbreviations
	locationAbbrevDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			locationAbbrevDict[row['ID']] = row['Abbreviation']

	# Load conversion between Ecocyc metabolite frame id's and metabolite id's from Feist. This is used for small-molecule/protein complexes.
	metaboliteEcocycToFeistIdConversion = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'ecocyc_to_feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			if row[1] == '+':
				metaboliteEcocycToFeistIdConversion[row[0]] = row[0].lower()
			else:
				metaboliteEcocycToFeistIdConversion[row[0]] = row[1]

	# Build cache of modified form reactions
	rebuild = True
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_rna_modification_reactions.json')) or rebuild:
		modFormRxn = {}
		parentsDict = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'rb') as csvfile:
			dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
			for row in dictreader:
				if len(json.loads(row['Modified form'])):
					for frameId in json.loads(row['Modified form']):
						unmodified_form = row['Frame ID']
						formation_reactions = getFormationReactions(frameId, unmodified_form)
						print 'Checked for class sub-species for ' + frameId

						if formation_reactions == []:
							print 'No reaction for ' + frameId
						else:
							print 'Loaded ' + frameId + ' formation reaction | ' + str(formation_reactions)
						modFormRxn[frameId] = formation_reactions
						parentsDict[frameId] = parents

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_rna_modification_reactions.json'),'wb') as jsonfile:
			jsonfile.write(json.dumps(modFormRxn, indent = 4))
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_rna_modification_parentclass.json'),'wb') as jsonfile:
			jsonfile.write(json.dumps(parentsDict, indent = 4))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_rna_modification_reactions.json'),'rb') as jsonfile:
		modFormRxn = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_rna_modification_parentclass.json'),'rb') as jsonfile:
		parentsDict = json.loads(jsonfile.read())

	rnaDict_modified = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			if len(json.loads(row['Modified form'])):
				for frameId in json.loads(row['Modified form']):
					RNA = rna()
					RNA.frameId = frameId
					RNA.unmodifiedForm = row['Frame ID']
					RNA.location = json.loads(row['Location'])

					rxn_raw = modFormRxn[RNA.frameId]

					if rxn_raw != []:
						production_reaction = []
						for rxn in rxn_raw:
							for rxn_species in rxn[2]:
								if int(float(rxn_species[1])) > 0 and rxn_species['id'] == RNA.frameId and rxn[1] == 'LEFT-TO-RIGHT':
									production_reaction.append(rxn)
								if int(float(rxn_species[1])) < 0 and rxn_species['id'] == RNA.frameId and rxn[1] == 'RIGHT-TO-LEFT':
									production_reaction.append(rxn)
								elif rxn[1] == 'UNKNOWN' and rxn_species['id'] == RNA.frameId:
									production_reaction.append(rxn)

						for rxn in production_reaction:
							fillInReaction(RNA, rxn, locationAbbrevDict, metaboliteEcocycToFeistIdConversion)

					rnaDict_modified[RNA.frameId] = RNA
	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna_modified.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = rnaDict_modified.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Unmodified Form', 'Location', 'Reaction ID', 'Reaction Enzyme', 'Reaction', 'Comments'])
		for key in keys:
			ribonuc = rnaDict_modified[key]
			csvwriter.writerow([ribonuc.frameId, ribonuc.unmodifiedForm, json.dumps(ribonuc.location), json.dumps(ribonuc.reactionId), json.dumps(ribonuc.reactionEnzymes), json.dumps(ribonuc.reaction), ribonuc.comments])

def getFormationReactions(frameId, unmodified_form):
	# Look for reactions that include the parent classes of frameid in reaction
	parents = []
	formation_reactions_raw = []
	getEcocycParents(str(frameId), parents)
	for p in parents:
		rxn = getEcocycModFormReactions(p)
		formation_reactions_raw.extend(rxn)
	print 'Checked for parents for ' + frameId

	# Look for class species in reaction and fill in with instance species
	formation_reactions = []
	for rxn in formation_reactions_raw:
		# Builds list of:
		# [[{'classid' : class id 1, 'instance id' : instance id 1}, {'classid' : class id 1, 'instance id' : instance id 2},...]
		#  ['classid' : class id 2, 'instance id' : instance id 3,...]]
		# This makes it easy to build a list of all possible unique combinations of instance frame id's so that we can build
		# a complete list of unique instance based reactions
		components_children = buildReactionInstanceFromClassList(rxn, frameId, unmodified_form)

		# Forms all cartesian-products
		# Example: ({'classid': 'VAL-tRNAs', 'instanceid': 'valT-tRNA'}, {'classid': 'Charged-VAL-tRNAs', 'instanceid': 'charged-valT-tRNA'})
		for cart_product in itertools.product(*components_children[:]):
			# Builds new reaction with class species replaced with instance species from the cartesian product
			new_rxn = buildInstanceReaction(cart_product, rxn)
			formation_reactions.append(new_rxn)
	return formation_reactions

def buildReactionInstanceFromClassList(rxn, modified_form, unmodified_form):
	components_children = []
	noUnmod = True
	for class_comp in [x for x in rxn['components'] if x['isclass'] == True]:
		children = []
		getEcocycChildren(class_comp['id'], children)
		if unmodified_form in children:
			children = [unmodified_form]
			noUnmod = False
		if modified_form in children:
			children = [modified_form]
		components_children.append([{'classid' : class_comp['id'], 'instanceid' : x} for x in children])
	if noUnmod:
		raise Exception, 'No unmodified form found!\n'	
	return components_children

def buildInstanceReaction(pairs_to_replace, rxn):
	# New reaction for cartesian-product
	new_rxn = copy.deepcopy(rxn)

	# Looks over all species in reaction
	for species in new_rxn['components']:
		addOriginal = False
		# Pulls out cartesian-product
		for species_to_replace in pairs_to_replace:
			count = 0
			if species['id'] == species_to_replace['classid']:
				# If species is the class-id from the cartesian product then replace with instance from product
				species['id'] = species_to_replace['instanceid']
				species['isclass'] = False
	return new_rxn

# Parse protein complexes
def parseComplexes():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.log'),'a')

	# Load compartment id --> single letter abbreviation data
	compartmentAbbrev = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			compartmentAbbrev[row[0]] = row[1]

	# Load monomer localization
	monomerCompartment = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			monomerCompartment[row[0]] = json.loads(row[3])

			# TODO: Might not need this here come back and delete. Modified forms have been removed.
			modifiedForm = json.loads(row[4])
			if modifiedForm != []:
				for m in modifiedForm:
					monomerCompartment[m] = json.loads(row[3])

	# Load conversion between Ecocyc metabolite frame id's and metabolite id's from Feist. This is used for small-molecule/protein complexes.
	metaboliteEcocycToFeistIdConversion = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'ecocyc_to_feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			if row[1] == '+':
				metaboliteEcocycToFeistIdConversion[row[0]] = row[0].lower()
			else:
				metaboliteEcocycToFeistIdConversion[row[0]] = row[1]

	# Build list of RNAs that could be included in complexes
	rnaList = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			rnaList.append(row[0])
			modifiedForm = json.loads(row[4])
			if modifiedForm != []:
				for mf in modifiedForm:
					rnaList.append(mf)


	# Build one complete list of protein complexes (includes protein-protein, protein-RNA, and protein-small molecule)
	# Add correct stoichiometry in this file (BioVelo query downloads dependencies)
	rebuild = False
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_protein_complexes_correct_stoich.csv')) or rebuild:
		newRows = []
		modifiedForms = []
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				newRows.append(row)
				modifiedForms.extend(row[3][1:-1].split(' '))
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_protein_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				newRows.append(row)
				modifiedForms.extend(row[3][1:-1].split(' '))
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_small_molecule_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				newRows.append(row)
				modifiedForms.extend(row[3][1:-1].split(' '))
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'Ecocyc_protein_complexes_merge.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')
			csvwriter.writerow(['Frame ID', 'Name', 'Stoichiometry', 'Modified form', 'Comments'])
			for row in newRows:
				if row[0] not in modifiedForms:
					csvwriter.writerow(row)

		rebuiltRow = []
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'Ecocyc_protein_complexes_merge.csv'),'rb') as csvfile:
			dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
			for row in dictreader:
				subunits = getEcocycComplexComponents(row['Frame ID'])
				if subunits == []:
					print 'ERROR!!! ' + row['Frame ID']
				print 'fetched subunits for ' + row['Frame ID']
				subunit_string = '('
				for i,s in enumerate(subunits):
					subunit_string += '('
					subunit_string += s[0]
					subunit_string += ', '
					subunit_string += str(s[1])
					if i == len(subunits) - 1:
						subunit_string += ')'
					else:
						subunit_string += ') '
				subunit_string += ')'
				row['Stoichiometry'] = subunit_string
				rebuiltRow.append([row['Frame ID'], row['Name'], row['Stoichiometry'], row['Modified form']])

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_protein_complexes_correct_stoich.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')
			csvwriter.writerow(['Frame ID', 'Name', 'Stoichiometry', 'Modified form', 'Comments'])
			rebuiltRow.sort()
			for row in rebuiltRow:
				csvwriter.writerow(row)

	# Build list of protein-protein complexes
	proteinComplexes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_protein_complexes_correct_stoich.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			proteinComplexes.append(row[0])

	# Build list of modified forms that could be included in complexes
	modifiedForms ={}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_protein_complexes_correct_stoich.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		
		for row in dictreader:
			modform = row['Modified form'][1:-1].split(' ')
			if modform != [""]:
				for m in modform:
					modifiedForms[m] = row['Frame ID']

	# Build list of protein monomers that have no known gene
	proteinMonomerUnknownGene = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'protein_monomers_unknown_genes.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			unkgene = row['Frame ID']
			proteinMonomerUnknownGene.append(unkgene)

	# Parse protein complex information
	proCompDict = {}
	saveRow = {}
	hasComplexSubunit = []
	proteinComplexUnknownGene = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_protein_complexes_correct_stoich.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			comp = proteinComplex()
			comp.frameId = row[0]
			comp.name = re.sub('<[^<]+?>', '', row[1])
			comp.modifiedForm = row[3][1:-1].split(' ')
			if comp.modifiedForm == [""]:
				comp.modifiedForm = []

			foundAllComponents = True
			components = row[2][2:-2].replace('"','').split(') (')
			if row[2] != '':
				for c in components:
					info = c.split(', ')
					frameId = info[0]
					stoich = int(info[1])

					if (frameId in proteinComplexes) or (frameId in modifiedForms.keys()):
						if frameId not in [x[0] for x in hasComplexSubunit]:
							hasComplexSubunit.append(comp.frameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break
					elif monomerCompartment.has_key(frameId):
						location = monomerCompartment[frameId]
						comp.addReactant(frameId, stoich, location)
					elif proCompDict.has_key(frameId):
						location = proCompDict[frameId].composition['product'][frameId]['compartment']
						comp.addReactant(frameId, stoich, location)
					elif metaboliteEcocycToFeistIdConversion.has_key(frameId):
						# TODO: Check location is correct
						location = ['CCO-CYTOSOL']
						comp.addReactant(metaboliteEcocycToFeistIdConversion[frameId], stoich, location)
					elif (frameId in rnaList):
						location = ['CCO-CYTOSOL']
						comp.addReactant(frameId, stoich, location)
					elif frameId in proteinMonomerUnknownGene:
						proteinComplexUnknownGene.append(comp.frameId)
						foundAllComponents = False
						break
					else:
						foundAllComponents = False
						s = 'Did not create a protein-protein complex for ' + comp.frameId
						writeOut(s, logFile)

				if foundAllComponents:
					comp.addProduct(comp.frameId, 1)
					comp.calculateLocation()
					comp.buildStringComposition(compartmentAbbrev)

					proCompDict[comp.frameId] = comp

			else:
				raise Exception, 'No stoichiometry found!\n'
	
	# Deal with protein complexes that have other protein complexes as subunits
	prev = 0
	breakCount = 0
	while len(hasComplexSubunit):
		this = len(hasComplexSubunit)
		#print this
		if prev == this:
			breakCount += 1
		if breakCount > 100:
			ipdb.set_trace()
		prev = this

		row = saveRow[hasComplexSubunit[0]]

		comp = proteinComplex()
		comp.frameId = hasComplexSubunit[0]
		comp.name = re.sub('<[^<]+?>', '', row[1])

		foundAllComponents = True
		components = row[2][2:-2].split(') (')
		if row[2] != '':
			for c in components:
				info = c.split(', ')
				frameId = info[0]
				stoich = int(info[1])

				if monomerCompartment.has_key(frameId):
					location = monomerCompartment[frameId]
					comp.addReactant(frameId, stoich, location)
				elif proCompDict.has_key(frameId):
					location = proCompDict[frameId].composition['product'][frameId]['compartment']
					comp.addReactant(frameId, stoich, location)
				elif metaboliteEcocycToFeistIdConversion.has_key(frameId):
					# TODO: Check location is correct
					location = ['CCO-CYTOSOL']
					comp.addReactant(metaboliteEcocycToFeistIdConversion[frameId], stoich, location)
				elif (frameId in rnaList):
					location = ['CCO-CYTOSOL']
					comp.addReactant(frameId, stoich, location)
				elif frameId in modifiedForms.keys():
					unmodifiedForm_frameId = modifiedForms[frameId]
					location = proCompDict[unmodifiedForm_frameId].composition['product'][unmodifiedForm_frameId]['compartment']
					comp.addReactant(frameId, stoich, location)
				else:
					foundAllComponents = False
					savePC = hasComplexSubunit.pop(0)
					hasComplexSubunit.append(savePC)

			if foundAllComponents:
				comp.addProduct(comp.frameId, 1)
				comp.calculateLocation()
				comp.buildStringComposition(compartmentAbbrev)

				proCompDict[comp.frameId] = comp
				hasComplexSubunit.pop(0)

	# Write complexes with subunits that have no known genes
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'protein_complexes_unknown_genes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Frame ID', 'Comments'])
		for u in proteinComplexUnknownGene:
			csvwriter.writerow([u])

	# Write complexes
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Frame ID', 'Name', 'Location', 'Composition', 'Composition', 'Modified form', 'Formation process', 'Comments'])
		
		keys = proCompDict.keys()
		keys.sort()
		for key in keys:
			c = proCompDict[key]
			csvwriter.writerow([c.frameId, c.name, json.dumps(c.composition['product'][c.frameId]['compartment']), c.compositionString, json.dumps(c.composition), json.dumps(c.modifiedForm), c.formationProcess, c.comments])

	logFile.close()

def parseComplexes_modified():
	# Load location abbreviations
	locationAbbrevDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			locationAbbrevDict[row['ID']] = row['Abbreviation']

	# Load conversion between Ecocyc metabolite frame id's and metabolite id's from Feist. This is used for small-molecule/protein complexes.
	metaboliteEcocycToFeistIdConversion = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'ecocyc_to_feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			if row[1] == '+':
				metaboliteEcocycToFeistIdConversion[row[0]] = row[0].lower()
			else:
				metaboliteEcocycToFeistIdConversion[row[0]] = row[1]

	# Build cache of modified form reactions
	rebuild = False
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_complex_modification_reactions.json')) or rebuild:
		modFormRxn = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'rb') as csvfile:
			dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
			for row in dictreader:
				if len(json.loads(row['Modified form'])):
					for frameId in json.loads(row['Modified form']):
						rxn = getEcocycModFormReactions(frameId)
						if rxn == []:
							print 'No reaction for ' + frameId
						print 'Loaded ' + frameId + ' formation reaction'
						modFormRxn[frameId] = rxn

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_complex_modification_reactions.json'),'wb') as jsonfile:
			jsonfile.write(json.dumps(modFormRxn, indent = 4))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'ecocyc_prot_complex_modification_reactions.json'),'rb') as jsonfile:
		modFormRxn = json.loads(jsonfile.read())


	proteinComplexDict_modifiedForm = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			if len(json.loads(row['Modified form'])):
				for frameId in json.loads(row['Modified form']):
					pc = proteinComplex()
					pc.frameId = frameId
					pc.unmodifiedForm = row['Frame ID']
					pc.location = json.loads(row['Location'])

					rxn_raw = modFormRxn[pc.frameId]

					if rxn_raw != []:
						production_reaction = []

						for rxn in rxn_raw:
							for rxn_species in rxn['components']:
								if int(float(rxn_species['coeff'])) > 0 and rxn_species['id'] == pc.frameId and rxn['direction'] == 'LEFT-TO-RIGHT':
									production_reaction.append(rxn)
								if int(float(rxn_species['coeff'])) < 0 and rxn_species['id'] == pc.frameId and rxn['direction'] == 'RIGHT-TO-LEFT':
									production_reaction.append(rxn)
								elif rxn['direction'] == 'UNKNOWN' and rxn_species['id'] == pc.frameId:
									production_reaction.append(rxn)

						for rxn in production_reaction:
							fillInReaction(pc, rxn, locationAbbrevDict, metaboliteEcocycToFeistIdConversion)


					proteinComplexDict_modifiedForm[pc.frameId] = pc

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes_modified.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinComplexDict_modifiedForm.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Unmodified Form', 'Location', 'Reaction ID', 'Reaction Enzyme', 'Reaction', 'Comments'])
		for key in keys:
			pc = proteinComplexDict_modifiedForm[key]
			csvwriter.writerow([pc.frameId, pc.unmodifiedForm, json.dumps(pc.location), json.dumps(pc.reactionId), json.dumps(pc.reactionEnzymes), json.dumps(pc.reaction), pc.comments])

def getEcocycComplexComponents(cmplx):
	websvcUrl = "http://websvc.biocyc.org/getxml?ECOLI:%s" % cmplx
	dom = xml.dom.minidom.parse(urllib.urlopen(websvcUrl))
	L = []
	for component in dom.getElementsByTagName("component"):
		elemProt = component.getElementsByTagName("Protein")
		elemRna = component.getElementsByTagName("RNA")
		elemCmpnd = component.getElementsByTagName("Compound")
		if len(elemProt) > 0:
			fId = elemProt[0].getAttribute("frameid")
		elif len(elemRna) > 0:
			fId = elemRna[0].getAttribute("frameid")
		elif len(elemCmpnd) > 0:
			fId = elemCmpnd[0].getAttribute("frameid")
		else:
			raise Exception, "Don't have a frame id."
		elemCoeff = component.getElementsByTagName("coefficient")
		if len(elemCoeff) > 0:
			coeff = elemCoeff[0].childNodes[0].data
		else:
			coeff = u"1"
		L.append((fId, coeff))
	return L

# Parse transcription units
def parseTranscriptionUnits():
	# Load necessary gene information
	geneDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			newGene = gene()
			newGene.frameId = row[0]
			newGene.coordinate = int(row[4])
			newGene.length = int(row[5])
			newGene.direction = row[6]
			if newGene.direction == '+':
				newGene.left = newGene.coordinate
				newGene.right = newGene.coordinate + newGene.length
			else:
				newGene.left = newGene.coordinate - newGene.length
				newGene.right = newGene.coordinate
			geneDict[newGene.frameId] = newGene

	# Load promoters
	# Ignore promoters that are not associated with any transcription units
	promoterDict = featureDictionary()
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_promoters.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			if row[4] != '':
				newPro = promoter()
				newPro.name = row[1]
				newPro.frameId = row[0]
				newPro.sigma = parseSigmaFactors(row[2][1:-1])
				if row[3] != '':
					newPro.tss = int(row[3])
				promoterDict[newPro.frameId] = newPro

	# Load terminators
	terminatorDict = featureDictionary()
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rhoIndepTerm.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			newTerm = terminator()
			newTerm.name = row[1]
			newTerm.frameId = row[0]
			newTerm.left = int(row[2])
			newTerm.right = int(row[3])
			newTerm.rho = False
			terminatorDict[newTerm.frameId] = newTerm
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rhoDepTerm.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			newTerm = terminator()
			newTerm.name = row[1]
			newTerm.frameId = row[0]
			newTerm.left = int(row[2])
			newTerm.right = int(row[3])
			newTerm.rho = True
			terminatorDict[newTerm.frameId] = newTerm

	# Load transcription units
	transcriptionUnitDict = featureDictionary()
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_transcriptionUnits.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			# Get name/frameId
			frameId = row[1]
			tuName = row[0]
			# Get genes
			geneIdsList = row[4][1:-1].split(' ')
			geneList = [geneDict[x] for x in geneIdsList if geneDict.has_key(x)]
			# Get promoter
			if row[2] != '':
				pro = promoterDict[row[2][1:-1]]
			else:
				pro = None
			# Get terminator
			if row[3] != '':
				terminatorList = row[3][1:-1].split(' ')
				terminatorList = [terminatorDict[x] for x in terminatorList]
			else:
				terminatorList = []
			# Check if TU name already used
			if transcriptionUnitDict.has_key(frameId):
				raise chromosomeException, 'ID already used!\n'
			# Build transcription unit
			newTU = buildTranscriptionUnit(tuName, frameId, pro, terminatorList, geneList, promoterDict, terminatorDict)
			if newTU != None:
				transcriptionUnitDict[newTU.frameId] = newTU

	# Write promoters
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'promoters.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['frameId', 'Name', 'Absolute +1 Position', 'Sigma', 'Direction', 'Component of'])
		
		keys = promoterDict.keys()
		keys.sort()
		for key in keys:
			p = promoterDict[key]
			csvwriter.writerow([p.frameId, p.name, p.tss, json.dumps(p.sigma), p.direction, json.dumps(p.cmpOf)])

	# Write terminators
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'terminators.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['frameId', 'Name', 'Left', 'Right', 'Rho Dependent', 'Component of'])
		
		keys = terminatorDict.keys()
		keys.sort()
		for key in keys:
			t = terminatorDict[key]
			csvwriter.writerow([t.frameId, t.name, t.left, t.right, t.rho, json.dumps(t.cmpOf)])

	# Write TUs
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'transcriptionUnits.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['frameId', 'Name', 'Left', 'Right', 'Direction', 'Genes', 'Promoter', 'Terminators'])
		
		keys = transcriptionUnitDict.keys()
		keys.sort()
		for key in keys:
			t = transcriptionUnitDict[key]
			csvwriter.writerow([t.frameId, t.name, t.left, t.right, t.direction, json.dumps(t.genes), t.promoter, json.dumps(t.terminators)])

def buildTranscriptionUnit(tuName, tuFrameId, pro, terminatorList, geneList, promoterDict, terminatorDict):
	# Check which components are known from Ecocyc
	hasPromoter = False
	hasSigma = False
	hasTss = False
	if pro != None:
		hasPromoter = True
	if hasPromoter and len(pro.sigma):
		hasSigma = True
	if hasPromoter and pro.tss != None:
		hasTss = True

	hasTerminator = False
	if len(terminatorList):
		hasTerminator = True

	hasGenes = False
	if len(geneList):
		hasGenes = True

	# Create new transcription unit
	newTU = transcriptionUnit()
	newTU.name = tuName
	newTU.frameId = tuFrameId

	# Add genes
	for gene in geneList:
		newTU.genes.append(gene.frameId)

	# Figure out direction, left, and right from genes
	# Assume:
	# - that TU direction is the same as first gene
	# - if TSS is known that TU start is at TSS
	# - if no TSS is known for promoter that it is the first nucleotide of the first gene
	# - if no terminator known then end of TU is end of last gene
	# - if terminator is known then and of TU is end of last terminator
	# - TU with only psuedo-genes in them are not included
	# - If no sigma factor is known then assume it is Sigma D

	if hasGenes:
		newTU.direction = geneList[0].direction

		# Figure out start
		if hasPromoter and hasTss:
			if newTU.direction == '+':
				newTU.left = pro.tss
				newTU.start = newTU.left
			elif newTU.direction == '-':
				newTU.right = pro.tss
				newTU.start = newTU.right
		elif hasPromoter and not hasTss:
			if newTU.direction == '+':
				newTU.left = getMinCoord(geneList)
				newTU.start = newTU.left
				pro.tss = newTU.start
			elif newTU.direction == '-':
				newTU.right = getMaxCoord(geneList)
				newTU.start = newTU.right
				pro.tss = newTU.start
		else:
			if newTU.direction == '+':
				newTU.left = getMinCoord(geneList)
				newTU.start = newTU.left
			elif newTU.direction == '-':
				newTU.right = getMaxCoord(geneList)
				newTU.start = newTU.right
		# Figure out end
		if hasTerminator:
			if newTU.direction == '+':
				newTU.right = max([term.right for term in terminatorList])
				newTU.end = newTU.right
			elif newTU.direction == '-':
				newTU.left = min([term.left for term in terminatorList])
				newTU.end = newTU.left
		else:
			if newTU.direction == '+':
				newTU.right = getMaxCoord(geneList) + 1
				newTU.end = newTU.right
			elif newTU.direction == '-':
				newTU.left = getMinCoord(geneList)- 1
				newTU.end = newTU.left

		# Add promoter
		if hasPromoter:
			# If a promoter exists for this transcription unit use it
			newTU.promoter = pro.frameId
			pro.direction = newTU.direction
			pro.cmpOf.append(newTU.frameId)
			if not hasSigma:
				pro.sigma = ['D']
		else:
			# Otherwise place upstream of first gene by one nucleotide
			newPro = promoter()
			newPro.name = 'p_WC_' + newTU.name
			newPro.frameId = generatePromoterFrame(promoterDict)
			newPro.direction = newTU.direction
			newPro.tss = newTU.start
			newPro.sigma = ['D']
			newPro.cmpOf = [newTU.frameId]
			newTU.promoter = newPro.frameId
			promoterDict[newPro.frameId] = newPro

		# Add terminator
		if hasTerminator:
			for term in terminatorList:
				newTU.terminators.append(term.frameId)
			for term in terminatorList:
				term.cmpOf.append(newTU.frameId)
		else:
			# If we are building terminator here it is only associated
			# with this transcription unit so direction is known.
			newTerm = terminator()
			newTerm.name = 'TERM_WC_' + newTU.name
			newTerm.frameId = generateTerminatorFrame(terminatorDict)
			if newTU.direction == '+':
				newTerm.left = newTU.right + 1
				newTerm.right = newTU.left
			elif newTU.direction == '-':
				newTerm.left = newTU.left - 1
				newTerm.right = newTU.left
			newTerm.cmpOf = [newTU.frameId]
			newTU.terminators.append(newTerm.frameId)
			terminatorDict[newTerm.frameId] = newTerm

		return newTU

# Parse metabolites
def parseMetabolites():
	metDict = {}
	# Parse Feist metabolites
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()

		for row in csvreader:
			newMet = metabolite()
			newMet.frameId = row[0]
			newMet.name = row[1]
			if row[5] != '':
				newMet.neutralFormula = row[5]
			else:
				newMet.neutralFormula = row[2]
			# Properties at pH 7
			newMet.addPHProp(pH = 7,formula = row[8], charge = row[9])
			# Properties at pH 7.2
			newMet.addPHProp(pH = 7.2, formula = row[11], charge = row[12])

			metDict[newMet.frameId] = newMet

	# Parse metabolites in Ecocyc and needed for complexation but not in Feist, and metabolites that Feist need to have added
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'ecocyc_to_feist_metabolites.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

		for row in csvreader:
			if row[3] != '':
				newMet = metabolite()
				newMet.frameId = row[0].lower()
				if row[2] != '':
					newMet.neutralFormula = row[2]
				else:
					newMet.neutralFormula = row[3]
				# Properties at pH 7
				newMet.addPHProp(pH = 7,formula = row[3], charge = row[4])
				# Properties at pH 7.2
				newMet.addPHProp(pH = 7.2, formula = row[5], charge = row[6])

				metDict[newMet.frameId] = newMet

	# Remove metabolites not being modeled as metabolites
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'metabolites_to_remove.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			metDict.pop(row[0])

	# Parse objective function
	# - Reactants
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Feist_objective.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		csvreader.next()
		csvreader.next()
		for row in csvreader:
			metId = row[7].replace('*','')
			objectiveRead = float(row[6]) # mmol/gDSW
			g1 = 10**-3 # mol/mmol
			g2 = 6.02*(10**23) # molecules/mol
			g3 = 2.8*(10**-13) # gDSW/cell
			objective = objectiveRead * g1 * g2 * g3 # molecules / cell
			metDict[metId].biomassConc = objective
			bmc = None
			if row[8] == 'cytoplasm':
				bmc = 'c'
			elif row[8] == 'periplasm':
				bmc = 'p'
			elif row[8] == 'extra-cellular face':
				bmc = 'o'
			elif row[8] == 'cytoplasm / periplasm':
				bmc = 'c'
			metDict[metId].biomassCompartment = bmc
			if metId == 'h2o':
				break

	# - Products
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Feist_objective.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for i in range(70):
			csvreader.next()
		for row in csvreader:
			metId = row[7].replace('*','')
			objectiveRead = float(row[6]) # mmol/gDSW
			g1 = 10**-3 # mol/mmol
			g2 = 6.02*(10**23) # molecules/mol
			g3 = 2.8*(10**-13) # gDSW/cell
			objective = objectiveRead * g1 * g2 * g3 * (-1.) # molecules / cell
			metDict[metId].biomassConc = objective
			bmc = None
			if row[8] == 'cytoplasm':
				bmc = 'c'
			elif row[8] == 'periplasm':
				bmc = 'p'
			elif row[8] == 'extra-cellular face':
				bmc = 'o'
			elif row[8] == 'cytoplasm / periplasm':
				bmc = 'c'
			metDict[metId].biomassCompartment = bmc
			if metId == 'pi':
				break

	# Write a file of all the fake metabolites
	proteinLocations = loadMonomerAndComplexLocations()
	locationAbbrev = loadLocationAbbrev()
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'pseudo_metabolites.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = metDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Name', 'Neutral formula','What is this?','Import/export frame ID'])
		for key in keys:
			m = metDict[key]
			if m.notRealMetabolte == True:
				cofactorName = None
				exchangeFrameId = None

				# Pseudo-metabolite - Need to add exchange reaction with cytoplasm or whichever compartment is required
				if m.frameId.count('ACP'):
					cofactorName = 'ACP'
					exchangeFrameId = ['EG50003-MONOMER']

				if m.frameId == 'alpp':
					cofactorName = 'apolipoprotein'
					exchangeFrameId = ['EG10544-MONOMER']

				if m.frameId == 'glutrna':
					cofactorName = 'Glu-tRNA'
					exchangeFrameId = ['charged-gltT-tRNA','charged-gltU-tRNA','charged-gltV-tRNA','charged-gltW-tRNA']
				if m.frameId == 'trnaglu':
					cofactorName = 'Glu-tRNA'
					exchangeFrameId = ['gltT-tRNA','gltU-tRNA','gltV-tRNA','gltW-tRNA']

				if m.frameId == 'CPD0-2342':
					cofactorName = 'CPD0-2342'
					exchangeFrameId = ['CPD0-2342']

				if m.frameId == 'dsbdrd':
					cofactorName = 'dsbD'
					exchangeFrameId = ['DSBD-MONOMER']
				if m.frameId == 'dsbdox':
					cofactorName = 'dsbD'
					exchangeFrameId = ['DSBDOXI-MONOMER']

				if m.frameId == 'dsbard':
					cofactorName = 'dsbA'
					exchangeFrameId = ['DISULFOXRED-MONOMER']
				if m.frameId == 'dsbaox':
					cofactorName = 'dsbA'
					exchangeFrameId = ['MONOMER0-4152']

				if m.frameId == 'dsbcrd':
					cofactorName = 'dsbC'
					exchangeFrameId = ['DSBC-CPLX']
				if m.frameId == 'dsbcox':
					cofactorName = 'dsbC'
					exchangeFrameId = ['CPLX0-8002']

				if m.frameId == 'dsbgrd':
					cofactorName = 'dsbG'
					exchangeFrameId = ['DSBG-CPLX']
				if m.frameId == 'dsbgox':
					cofactorName = 'dsbG'
					exchangeFrameId = ['CPLX0-8004']

				if m.frameId == 'fldox':
					cofactorName = 'flavodoxin'
					exchangeFrameId = ['OX-FLAVODOXIN1','OX-FLAVODOXIN2']
				if m.frameId == 'fldrd':
					cofactorName = 'flavodixin'
					exchangeFrameId = ['FLAVODOXIN1-MONOMER','FLAVODOXIN2-MONOMER']

				if m.frameId == 'grxox':
					cofactorName = 'glutaredoxin'
					exchangeFrameId = ['GLUTAREDOXIN-MONOMER','OX-GLUTAREDOXIN-B','OX-GLUTAREDOXIN-C','EG12181-MONOMER']
				if m.frameId == 'grxrd':
					cofactorName = 'glutaredoxin'
					exchangeFrameId = ['RED-GLUTAREDOXIN','GRXB-MONOMER','GRXC-MONOMER','EG12181-MONOMER']

				if m.frameId == 'lpp':
					cofactorName = 'lipoprotein'
					exchangeFrameId = ['G7644-MONOMER']

				if m.frameId == 'trdox':
					cofactorName = 'thioredoxin'
					exchangeFrameId = ['OX-THIOREDOXIN-MONOMER','OX-THIOREDOXIN2-MONOMER']
				if m.frameId == 'trdrd':
					cofactorName = 'thioredoxin'
					exchangeFrameId = ['RED-THIOREDOXIN-MONOMER','RED-THIOREDOXIN2-MONOMER']

				if m.frameId == 'flvubrdox':
					cofactorName = 'flavorubredoxin'
					exchangeFrameId = ['CPLX0-2','CPLX0-1']
				if m.frameId == 'flvubrdrd':
					cofactorName = 'flavorubredoxin'
					exchangeFrameId = ['CPLX0-2','CPLX0-1']
				
				if exchangeFrameId != None:
					for ee in exchangeFrameId:
						location = locationAbbrev[proteinLocations[ee][0]]
						m.equivalentEnzyme.append(ee + '[' + location + ']')

				csvwriter.writerow([m.frameId, m.name, m.neutralFormula, cofactorName, json.dumps(exchangeFrameId)])

	# Write output for metabolites
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'metabolites.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = metDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Name', 'Neutral formula', 'pH 7.2 formula', 'pH 7.2 charge', 'pH 7.2 Weight', 'Media Concentration (mM)', 'Biomass concentration (molecules/cell)', 'Biomass location', 'Maximum exchange rate (mmol/gDSW/hr)', 'Fake metabolite', 'Equivalent enzyme frameId', 'Comments'])
		for key in keys:
			m = metDict[key]
			csvwriter.writerow([m.frameId, m.name, m.neutralFormula, m.pHProps[7.2]['formula'], m.pHProps[7.2]['charge'], m.pHProps[7.2]['weight'], m.mediaConc, m.biomassConc, m.biomassCompartment, m.exchangeRate, m.notRealMetabolte, json.dumps(m.equivalentEnzyme), m.comments])

def loadMonomerAndComplexLocations():
	proteinLocations = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			proteinLocations[row[0]] = json.loads(row[3])
			modifiedForm = json.loads(row[4])
			if len(modifiedForm):
				for m in modifiedForm:
					proteinLocations[m] = json.loads(row[3])
	proteinLocations['CPD0-2342'] = ['CCO-CYTOSOL']

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			proteinLocations[row[0]] = json.loads(row[2])
			modifiedForm = json.loads(row[5])
			if len(modifiedForm):
				for m in modifiedForm:
					proteinLocations[m] = json.loads(row[2])

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			proteinLocations[row[0]] = json.loads(row[3])
			modifiedForm = json.loads(row[4])
			if len(modifiedForm):
				for m in modifiedForm:
					proteinLocations[m] = json.loads(row[3])

	return proteinLocations

def loadLocationAbbrev():
	# Load location abbreviations
	locationAbbrev = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			locationAbbrev[row[0]] = row[1]
	return locationAbbrev

# Parse reactions
def parseReactions():
	# Load reactions
	reactDict = {}
	rp = reactionParser()

	# Add Feist reactions
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Feist_reactions.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()

		for row in csvreader:
			reac = buildReaction(rp,row)
			reactDict[reac.frameId] = reac

	# Add reactions not in Fesit
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'reactions_to_add.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()

		for row in csvreader:
			reac = buildReaction(rp,row)
			reactDict[reac.frameId] = reac

	# Remove reactions
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'reactions_to_remove.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		csvreader.next()
		for row in csvreader:
			reactDict.pop(row[0])

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'reactions.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = reactDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'Name', 'Process', 'EC', 'Stoichiometry (pH 7.2)', 'Enzyme', 'Direction','Comments'])
		for key in keys:
			r = reactDict[key]
			csvwriter.writerow([r.frameId, r.name, r.process, r.EC, r.stoich, json.dumps(r.enzyme), r.direction, r.comments])

def buildReaction(rp,row):
	reac = reaction()
	reac.frameId = row[0]
	reac.name = row[1]
	reac.process = 'Metabolism'
	if row[5] != '':
		reac.EC = row[5]
	reac.stoich = row[2]
	reac.direction = row[3]

	# Figure out enzymes
	if row[6] == '':
		# Nothin known
		reac.enzyme = None
	elif re.match("b([0-9])", row[6]) != None:
		bnum = row[6]

		pMFrameId = rp.getPMFrame(bnum)

		reac.enzyme = []
		for e in pMFrameId:
			reac.enzyme.append([e])
	else:
		if rp.manualAnnotationDict.has_key(row[0]):
			enzymes = rp.findEnzymeManualCuration(rp.manualAnnotationDict[row[0]]['annotation'])
			cofactors = []
		else:
			enzymeInfo = rp.findEnzyme(row[6], row)
			enzymes = enzymeInfo['enzymes']
			cofactors = enzymeInfo['cofactors']
		if len(enzymes[0]):
			for i,e in enumerate(enzymes):
				reac.enzyme.append(e)
		else:
			# Catches spontanious reactions
			reac.enzyme = None
		for c in cofactors:
			reac.requiredCofactors.append(c)
		reac.requiredCofactors.sort()
	return reac


# Parse enzyme kinetics
def parseEnzymeKinetics():
	enzKinDict = {}

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'kinetics.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = enzKinDict.keys()
		keys.sort()
		csvwriter.writerow(['Frame ID', 'kcat forward', 'kcat forward units', 'kcat reverse', 'kcat reverse units','Comments'])
		for key in keys:
			e = enzKinDict[key]
			csvwriter.writerow([e.frameId, e.kcat_forward, e.forward_units, e.kcat_reverse, e.reverse_units, e.comments])

# Utility functions
def splitBigBracket(s):
	s = s[2:-2]
	s = s.replace('"','')
	s = s.split(') (')
	return s

def splitSmallBracket(s):
	s = s.split(', ')
	frameId = s[0]
	description = s[1]
	return {'frameId' : frameId, 'description' : description}

def writeOut(s, file):
	print s
	file.write(s + '\n')

class featureDictionary(dict):
	def __init__(self):
		self.count = 0

def generatePromoterFrame(promoterFeatureDict):
	return generateUniqueFrame('PM_WC-', promoterFeatureDict)

def generateTerminatorFrame(terminatorFeatureDict):
	return generateUniqueFrame('TERM_WC-', terminatorFeatureDict)

def generateUniqueFrame(startStr, featureDict):
	frameId = startStr + str(featureDict.count)
	featureDict.count += 1
	return frameId

def parseSigmaFactors(line):
	if line == '':
		return []
	else:
		a = re.findall("sigma (?P<name>[A-Z]+)", line)
		return a

def getMinCoord(geneList):
	return min([gene.left for gene in geneList])

def getMaxCoord(geneList):
	return max([gene.right for gene in geneList])

# Define data type classes
class enzyme:
	def __init__(self):
		self.frameId = None
		self.kcat_forward = None
		self.kcat_reverse = None
		self.forward_units = None
		self.reverse_units = None
		self.comments = ''

class metabolite:
	def __init__(self):
		# Read in properties
		self.frameId = None
		self.name = None
		self.neutralFormula = None
		self.pHProps = {}

		self.mediaConc = None
		self.biomassConc = None
		self.biomassCompartment = None
		self.biomassNewFlux = None
		self.biomassRecycle = None
		self.exchangeRate = None
		self.notRealMetabolte = None
		self.equivalentEnzyme = []
		self.comments = ''

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'elements.json'),'rb') as jsonfile:
			self.elementDict = json.loads(jsonfile.read())

	def addPHProp(self, pH, formula, charge):
		self.pHProps[pH] = {'formula' : formula, 'charge' : int(charge), 'weight' : "%0.10f" % self.calculateWeight(formula)}

	def calculateWeight(self, formula):
		weight = 0.
		element_stoich = re.findall("[A-Z][a-z]*[0-9]*", formula)
		for value in element_stoich:
			m = re.search("(?P<letters>[A-Za-z]*)(?P<numbers>[0-9]*)", value)
			element = m.group('letters')
			stoich = m.group('numbers')
			if stoich == '':
				stoich = 1
			if self.elementDict.has_key(element):
				weight += int(stoich) * self.elementDict[element]['mass']
			else:
				if 'Non-element in formula. Weight calculated for known elements.\n' not in self.comments:
					self.comments += 'Non-element in formula. Weight calculated for known elements.\n'
					self.notRealMetabolte = True
		return weight

class reaction:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.process = None
		self.EC = None
		self.stoich = None
		self.enzyme = []
		self.requiredCofactors = []
		self.direction = None
		self.comments = ''

class reactionParser:
	def __init__(self):
		self.synDictFrameId = self.loadSynDict()
		self.protMonomerFrameId = self.loadProteinMonomerFrameIds()
		self.monomerToComplex = self.loadMonomerToComplex()
		self.fakeMetaboliteFrameIds = self.loadFakeMetaboliteFrameIds()
		self.fakeMetaboliteDict = self.loadFakeMetabolites()
		self.manualAnnotationDict = self.loadManualAnnotation()

	def loadSynDict(self):
		# Load gene frameId synonym dictionary for blatter numbers in Fiest 
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'gene_frameId_synonyms.json'),'rb') as jsonfile:
			synDictFrameId = json.loads(jsonfile.read())
		return synDictFrameId

	def loadProteinMonomerFrameIds(self):
		protMonomerFrameId = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			for row in csvreader:
				if not protMonomerFrameId.has_key(row[2]):
					protMonomerFrameId[row[2]] = [row[0]]
				else:
					protMonomerFrameId[row[2]].append(row[0])
		return protMonomerFrameId

	def loadMonomerToComplex(self):
		monomerOrComplexToComplex = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			for row in csvreader:
				compositionDict = json.loads(row[4])
				subunits = compositionDict['reactant'].keys()

				cmplxFrameId = row[0]
				monomerOrComplexToComplex[cmplxFrameId] = subunits
		
		monomerToComplex = {}
		for cmplxFrameId in monomerOrComplexToComplex.iterkeys():
			monomers = []
			# Edits monomers in place recursivly
			self.iterateTree(cmplxFrameId, monomers, monomerOrComplexToComplex)
			monomers = set(monomers)
			monomers = list(monomers)
			monomers.sort()
			monomers = tuple(monomers)
			monomerToComplex[monomers] = cmplxFrameId

		return monomerToComplex

	def iterateTree(self, cmplxFrameId, monomers, monomerOrComplexToComplex):
		for subunit in monomerOrComplexToComplex[cmplxFrameId]:
			if monomerOrComplexToComplex.has_key(subunit):
				self.iterateTree(subunit, monomers, monomerOrComplexToComplex)
			else:
				monomers.append(subunit)

	def loadFakeMetaboliteFrameIds(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'pseudo_metabolites.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			frameIdList = []
			for row in csvreader:
				frameId = json.loads(row[4])
				if frameId != None:
					frameIdList.extend(frameId)
		frameIdSet = set(frameIdList)
		uniqueFrameIdList = list(frameIdSet)
		return uniqueFrameIdList

	def loadFakeMetabolites(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', 'pseudo_metabolites.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			fakeMetaboliteDict = {}
			for row in csvreader:
				fakeMetaboliteDict[row[0]] = json.loads(row[4])
		return fakeMetaboliteDict

	def loadManualAnnotation(self):
		# Add manual annotation of reactino enzymes from Feist
		manualAnnotateDict = {}
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'reaction_enzyme_association.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			for row in csvreader:
				manualAnnotateDict[row[0]] = {'annotation' : row[1], 'comments' : row[2]}
		return manualAnnotateDict

	def getPMFrame(self, bnum):
		if self.synDictFrameId.has_key(bnum):
			geneFrameId = self.synDictFrameId[bnum]
		else:
			print 'bnum not found for ' + bnum
			return

		if self.protMonomerFrameId.has_key(geneFrameId):
			pMFrameId = self.protMonomerFrameId[geneFrameId]
		else:
			print 'protein monomer not found for ' + geneFrameId
			return

		return pMFrameId

	def findEnzyme(self, line, row = 'No row specified'):
		enzymes = []
		cofactors = []
		enzymesRaw = line.split('or')
		for i in range(len(enzymesRaw)):
			enzymes.append([])

		for i,e in enumerate(enzymesRaw):
			# Check for spontanious reaction
			if e.count('s0001') == 0:
				# Find all bnumbers and their corresponding monomers
				bnums = re.findall("(b[0-9]+)", e)
				monomers = []
				for b in bnums:
					monomers.extend(self.getPMFrame(b))

				# Check to see if any monomers are actually fake metabolites/cofactors
				for j,m in enumerate(monomers):
					if m in self.fakeMetaboliteFrameIds:
						cofac = monomers.pop(j)
						cofactors.append(cofac)

				# Sort monomers and cast to tuple for hash
				monomers.sort()
				monomers = tuple(monomers)
				if self.monomerToComplex.has_key(monomers) and len(monomers) > 1:
					# If this is actually a complex formed from more than on bnumber
					enzymes[i].append(self.monomerToComplex[monomers])
				elif len(monomers) == 1:
					# This is just a monomer in an OR statement
					enzymes[i].append(monomers[0])
				elif len(monomers) == 0 and len(cofactors) > 0:
					pass
				else:
					enzymes[i].append('UNKNOWN')
					print 'No enzyme complex found for subunits: ' + str(monomers)
					print str(row[:3])
					print str(e)
					print '---'
		return {'enzymes' : enzymes, 'cofactors' : cofactors}

	def findEnzymeManualCuration(self, line):
		if line == '':
			# Catches reactions where enzyme was removed
			return [[]]
		enzymes = []
		cofactors = []
		enzymesRaw = line.split('or')
		for i in range(len(enzymesRaw)):
			enzymes.append([])

		for i,e_raw in enumerate(enzymesRaw):
			e_raw = e_raw.replace(' ','')
			e_raw = e_raw.replace('(','')
			e_raw = e_raw.replace(')','')
			e_final = e_raw.split('and')
			for ee in e_final:
				enzymes[i].append(ee)
		return enzymes

class gene:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.symbol = None
		self.type = None
		self.coordinate = None
		self.length = None
		self.direction = None
		self.expression = 0.
		self.halfLife = None
		self.productFrameId = None
		self.splices = ()
		self.sequenceSubstitution = ()
		self.comments = None

class transcriptionUnit:
	def __init__(self):
		self.frameId = None
		self.name = None

		self.genes = []
		self.promoter = None
		self.terminators = []

		self.left = 0
		self.right = 0
		self.direction = ''
		self.start = 0
		self.end = 0

class promoter:
	def __init__(self):
		self.name = None
		self.direction = ''
		self.tss = None
		self.frameId = None
		self.cmpOf = []

		self.sigma = []
		self.affinity = 0
		self.rnas = {}

class terminator:
	def __init__(self):
		self.name = None
		self.left = 0
		self.right = 0
		self.frameId = None
		self.cmpOf = []
		self.rho = None

class proteinMonomer:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = []
		self.gene = None
		self.modifiedForm = None
		self.unmodifiedForm = None
		self.reactionId = []
		self.reaction = []
		self.reactionEnzymes = []
		self.comments = ''

class rna:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = []
		self.gene = None
		self.modifiedForm = None
		self.unmodifiedForm = None
		self.reactionId = []
		self.reaction = []
		self.reactionEnzymes = []
		self.comments = ''

class proteinComplex:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.composition = {'reactant' : {}, 'product' : {}}
		self.compositionString = ''
		self.formationProcess = 'Complexation'
		self.modifiedForm = []
		self.unmodifiedForm = None
		self.reactionId = []
		self.reaction = []
		self.location = []
		self.reactionEnzymes = []
		self.comments = ''

	def addReactant(self, name, stoich, location):
		self.composition['reactant'][name] = {'stoichiometry' : None, 'compartment' : None}
		self.composition['reactant'][name]['stoichiometry'] = stoich
		self.composition['reactant'][name]['compartment'] = location

	def addProduct(self, name, stoich):
		self.composition['product'][name] = {'stoichiometry' : None, 'compartment' : None}
		self.composition['product'][name]['stoichiometry'] = stoich
		self.composition['product'][name]['compartment'] = None

	def buildStringComposition(self, compartmentDict):
		s = ''
		subComp = self.composition['reactant'].keys()
		locationSet = [self.composition['reactant'][key]['compartment'][0] for key in self.composition['reactant'].iterkeys()]
		locationSetProduct = [self.composition['product'][key]['compartment'][0] for key in self.composition['product'].iterkeys()]
		locationSet.extend(locationSetProduct)
		locationSet = sets.Set(locationSet)
		sameLocation = False
		if len(locationSet) == 1:
			sameLocation = True

		if sameLocation:
			s += '[' + compartmentDict[self.composition['product'][self.frameId]['compartment'][0]] + ']: '

		for i in range(len(subComp)):
			c = subComp[i]
			if sameLocation:
				compartment = ''
			else:
				compartment = '[' + compartmentDict[self.composition['reactant'][c]['compartment'][0]] + ']'

			if self.composition['reactant'][c]['stoichiometry'] == 1:
				s += c + compartment + ' '
			else:
				stoich = self.composition['reactant'][c]['stoichiometry']
				s += '(' + str(stoich) + ') ' + c + compartment + ' '
			if i != len(subComp) - 1:
				s += '+ '
			elif sameLocation:
				s += '==> ' + self.frameId
			else:
				s += '==> ' + self.frameId + '[' + compartmentDict[self.composition['product'][self.frameId]['compartment'][0]] + ']'
		self.compositionString = s

	def calculateLocation(self):
		locationPossible = []
		for reactantId in self.composition['reactant'].iterkeys():
			reactantLocation = self.composition['reactant'][reactantId]['compartment']
			for loc in reactantLocation:
				locationPossible.append(loc)

		locationSet = sets.Set(locationPossible)
		locationList = []
		for item in locationSet:
			locationList.append(item)

		if len(locationList) == 1:
			location = locationList[0]
		elif 'CCO-PM-BAC-NEG' in locationList:
			location = 'CCO-PM-BAC-NEG'
		elif 'CCO-OUTER-MEM' in locationList:
			location = 'CCO-OUTER-MEM'
		elif 'CCO-CW-BAC-NEG' in locationList:
			location = 'CCO-CW-BAC-NEG'
		elif 'CCO-MEMBRANE' in locationList:
			location = 'CCO-MEMBRANE'
		elif 'CCO-CELL-PROJECTION' in locationList:
			location = 'CCO-CELL-PROJECTION'
		else:
			print 'NEED LOCATION HIERARCHY FOR ' + self.frameId

		self.composition['product'][self.frameId]['compartment'] = [location]

if __name__ == "__main__":
    main()
