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

t = time.strftime("%Y-%m-%d_%H_%M_%S", time.localtime())

def main():
	initalizeLog()
	getEcocyc(fetchNew = True)
	parseIntermediateFiles()

	parseGenes()
	parseLocations()
	parseProteinMonomers()
	parseRna()
	parseComplexes()

def initalizeLog():
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log')):
		os.makedirs(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log'))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.txt'),'wb') as txtfile:
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
	# Build protein-protein complexes
	bioVeloQuery = '[(x^frame-id, x^name, components): x <- ecoli^^protein-complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_protein_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)

	# Build protein-rna complexes
	bioVeloQuery = '[(x^frame-id, x^name, components): x <- ecoli^^Protein-RNA-Complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_rna_protein_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)

	# Build protein-small molecule complexes
	bioVeloQuery = '[(x^frame-id, x^name, components): x <- ecoli^^Protein-Small-Molecule-Complexes, components := [(c1^frame-id, c2): (c1, c2) <- protein-to-components x]]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_protein_small_molecule_complexes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)

	# Build genes
	bioVeloQuery = '[(G^FRAME-ID, G^NAME, G^LEFT-END-POSITION, G^RIGHT-END-POSITION, G^TRANSCRIPTION-DIRECTION, [(TMP^FRAME-ID, TMP^NAME): TMP <- P]) : G<-ECOLI^^All-Genes, P := [TMP: TMP <- G^PRODUCT], #P > 0]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_genes.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)

	bioVeloQuery = '[(G^NAME, G^FRAME-ID, G^SYNONYMS, G^ACCESSION-1, G^ACCESSION-2) : G<-ECOLI^^All-Genes, P := [TMP: TMP <- G^PRODUCT], #P > 0]'
	outFile = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_gene_synonyms.csv')
	generateEcocycFlatFile(bioVeloQuery, outFile)

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

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_name_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDict, indent = 4))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_frameId_synonyms.json'),'wb') as jsonfile:
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

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(unmodifiedForm))

def parseRnaTypes():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	rnaType = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[3] != '':
				if row[1].count('RRNA') > 0:
					rnaType[row[1]] = 'rRNA'
				elif row[1].count('tRNA') > 0 and not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'tRNA'
				elif not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'miscRNA'

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'rnaTypes.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(rnaType))

# Parse genes
def parseGenes():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.txt'),'a')

	# Load unmodified forms of RNA and proteins
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	# Load RNA types
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'rnaTypes.json'),'rb') as jsonfile:
		rnaType = json.loads(jsonfile.read())

	# Load synonym dictionaries
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_name_synonyms.json'),'rb') as jsonfile:
		synDict = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_frameId_synonyms.json'),'rb') as jsonfile:
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
							while geneDict.has_key(newGene.frameId + str(count)):
								count += 1
							geneDict[newGene.frameId + str(count)] = newGene
						else:
							geneDict[newGene.frameId] = newGene

						# Add RNA type
						if rnaType.has_key(newGene.productFrameId):
							newGene.type = rnaType[newGene.productFrameId]
						else:
							newGene.type = 'mRNA'

	# Parse splicing information
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'geneCoordinates.csv'),'rb') as csvfile:
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
			rrnaAverage = parameters['rRNA half life']
		if len(tRNAhl) == 0:
			s = 'No tRNA half lives measured. Using parameter.'
			writeOut(s, logFile)
			trnaAverage = parameters['tRNA half life']
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
	if not os.path.exists(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'proteinMonomerGravy.csv')):
		s = 'Calculating gravy for all genes'
		writeOut(s, logFile)
		gravy.main()
	else:
		s = 'Gravy already exists'
		writeOut(s, logFile)

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['ID', 'Name', 'Symbol', 'Type', 'Coordinate', 'Length', 'Direction', 'Expression', 'Half life', 'Product', 'Splices', '(absolute nt position, old, new)', 'Comments'])

		keys = geneDict.keys()
		keys.sort()

		for key in keys:
			g = geneDict[key]
			csvwriter.writerow([g.frameId, g.name, g.symbol, g.type, g.coordinate, g.length, g.direction, "%0.10f" % g.expression, g.halfLife, g.productFrameId, json.dumps(g.splices), json.dumps(g.sequenceSubstitution), g.comments])


	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'genes.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(geneDict.keys()))

	logFile.close()

# Parse Locations
def parseLocations():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.txt'),'a')

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

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'locations_equivalent_names.json'),'wb') as jsonfile:
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
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.txt'),'a')

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_name_synonyms.json'),'rb') as jsonfile:
		synDict = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_frameId_synonyms.json'),'rb') as jsonfile:
		synDictFrameId = json.loads(jsonfile.read())
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'genes.json'),'rb') as jsonfile:
		geneIdList = json.loads(jsonfile.read())
	# TODO: Might not need to use this last dict depending on the datasource
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'locations_equivalent_names.json'),'rb') as jsonfile:
		locationEquivDict = json.loads(jsonfile.read())

	proteinMonomerDict = {}
	geneToProteinMonomerDict = {}
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
				proteinMonomerDict[proteinMonomerFrameId].comments += 'Location information from Lopez Campistrous 2005.\n'
			else:
				location = []
			proteinMonomerDict[proteinMonomerFrameId].location = location

	print 'Locations found in E. coli K-12 for ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))

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
					proteinMonomerDict[proteinMonomerFrameId].comments += 'Location information from Han 2011.\n'
					proteinMonomerDict[proteinMonomerFrameId].location = location

	print 'Locations found in E. coli K-12 and B for ' + str(len([1 for pM in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if pM.location != []]))

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
							proteinMonomerDict[frameId].comments += 'Localization generated from unambiguous Ecocyc data.\n'

	# Fill in the rest using GRAVY calculation
	gravy = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'proteinMonomerGravy.csv'),'rb') as csvfile:
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
			proteinMonomerDict[key].comments = 'Location calculated to be either CCO-CYTOSOL or CCO-MEMBRANE based on GRAVY.\n'

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinMonomerDict.keys()
		keys.sort()
		csvwriter.writerow(['ID', 'Name', 'Gene', 'Location', 'Modified form', 'Comments'])
		for key in keys:
			pm = proteinMonomerDict[key]
			csvwriter.writerow([pm.frameId, pm.name, pm.gene, json.dumps(pm.location), json.dumps(pm.modifiedForm), pm.comments])
	logFile.close()

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
				r = rna()

				r.frameId = row[1]
				r.name = re.sub('<[^<]+?>', '', row[0])
				r.gene = row[3][1:-1]

				modifiedForm = row[7][1:-1].split('" "')
				if modifiedForm == ['']:
					modifiedForm = []
				r.modifiedForm = modifiedForm

				r.location = ['CCO-CYTOSOL']
				r.comments += 'Assumed in CCO-CYTOSOL\n'

				rnaDict[r.frameId] = r

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = rnaDict.keys()
		keys.sort()
		csvwriter.writerow(['ID', 'Name', 'Gene', 'Location', 'Modified form', 'Comments'])
		for key in keys:
			rnaToPrint = rnaDict[key]
			csvwriter.writerow([rnaToPrint.frameId, rnaToPrint.name, rnaToPrint.gene, json.dumps(rnaToPrint.location), json.dumps(rnaToPrint.modifiedForm), rnaToPrint.comments])

# Parse protein complexes
def parseComplexes():
	# Open log file
	logFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'log','log_' + t + '.txt'),'a')

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
			modifiedForm = json.loads(row[4])
			if modifiedForm != []:
				for m in modifiedForm:
					monomerCompartment[m] = json.loads(row[3])

	# Build list of protein-protein complexes
	proteinComplexes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			proteinComplexes.append(row[0])

	# Build list of protein-rna complexes
	rnaProteinComplexes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_protein_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			rnaProteinComplexes.append(row[0])

	# Build list of protein-small molecule complexes
	smallMolecProteinComplexes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_small_molecule_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			smallMolecProteinComplexes.append(row[0])

	# Parse protein complex information
	proCompDict = {}
	saveRow = {}
	hasProteinComplexSubunit = []
	hasRnaProteinComplexSubunit = []
	hasSmallMolecProteinComplexSubunit = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			comp = proteinComplex()
			comp.frameId = row[0]
			comp.name = re.sub('<[^<]+?>', '', row[1])

			foundAllComponents = True
			components = row[2][2:-2].replace('"','').split(') (')
			if row[2] != '':
				for c in components:
					info = c.split(', ')
					frameId = info[0]
					stoich = int(info[1])
					if (frameId in proteinComplexes):
						if frameId not in [x[0] for x in hasProteinComplexSubunit]:
							hasProteinComplexSubunit.append(comp.frameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break
					elif (frameId in rnaProteinComplexes):
						if frameId not in [x[0] for x in hasRnaProteinComplexSubunit]:
							hasRnaProteinComplexSubunit.append(comp.crameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break
					elif (frameId in smallMolecProteinComplexes):
						if frameId not in [x[0] for x in hasSmallMolecProteinComplexSubunit]:
							hasSmallMolecProteinComplexSubunit.append(comp.frameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break
					elif monomerCompartment.has_key(frameId):
						location = monomerCompartment[frameId]
						comp.addReactant(frameId, stoich, location)
					elif proCompDict.has_key(frameId):
						location = proCompDict[frameId].composition['product'][frameId]['compartment']
						comp.addReactant(frameId, stoich, location)
					else:
						foundAllComponents = False
						s = 'Did not create a protein-protein complex for ' + comp.frameId
						writeOut(s, logFile)

				if foundAllComponents:
					comp.addProduct(comp.frameId, 1)
					comp.calculateLocation()
					comp.buildStringComposition(compartmentAbbrev)

					proCompDict[comp.frameId] = comp

	# Deal with protein complexes that have other protein complexes as subunits
	# BUT NOT ones with RNA-protein or small molecule - protein subunits.
	pcSubunit = sets.Set(hasProteinComplexSubunit)
	rpcSubunit = sets.Set(hasRnaProteinComplexSubunit)
	smpcSubunit = sets.Set(hasSmallMolecProteinComplexSubunit)

	if len(pcSubunit.intersection(smpcSubunit)) != 0 or len(pcSubunit.intersection(rpcSubunit)) != 0:
		# Make sure that rna-protein complexes and small molecule-protein complexes do not depend on protein-protein complexes
		ipdb.set_trace()
	
	prev = 0
	breakCount = 0
	while len(hasProteinComplexSubunit):
		this = len(hasProteinComplexSubunit)
		#print this
		if prev == this:
			breakCount += 1
		if breakCount > 10000:
			ipdb.set_trace()
		prev = this

		row = saveRow[hasProteinComplexSubunit[0]]

		comp = proteinComplex()
		comp.frameId = hasProteinComplexSubunit[0]
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
				else:
					foundAllComponents = False
					savePC = hasProteinComplexSubunit.pop(0)
					hasProteinComplexSubunit.append(savePC)

			if foundAllComponents:
				comp.addProduct(comp.frameId, 1)
				comp.calculateLocation()
				comp.buildStringComposition(compartmentAbbrev)

				proCompDict[comp.frameId] = comp
				hasProteinComplexSubunit.pop(0)

	# Parse small-molecule-protein complxes
	ecocycToFeistId = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'Ecocyc_to_Feist.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			if row[1] == '+':
				ecocycToFeistId[row[0]] = row[0]
			else:
				ecocycToFeistId[row[0]] = row[1]

	smallMolecProCompDict = {}
	saveRow = {}
	hasProteinComplexSubunit = []
	hasRnaProteinComplexSubunit = []
	hasSmallMolecProteinComplexSubunit = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_small_molecule_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			comp = proteinComplex()
			comp.frameId = row[0]
			comp.name = re.sub('<[^<]+?>', '', row[1])

			foundAllComponents = True
			components = row[2][2:-2].replace('"','').split(') (')
			if row[2] != '':
				for c in components:
					info = c.split(', ')
					frameId = info[0]
					stoich = int(info[1])
					if (frameId in rnaProteinComplexes):
						if frameId not in [x[0] for x in hasRnaProteinComplexSubunit]:
							hasRnaProteinComplexSubunit.append(comp.crameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break
					elif (frameId in smallMolecProteinComplexes):
						if frameId not in [x[0] for x in hasSmallMolecProteinComplexSubunit]:
							hasSmallMolecProteinComplexSubunit.append(comp.frameId)
							saveRow[comp.frameId] = row
						foundAllComponents = False
						break

					elif ecocycToFeistId.has_key(frameId):
						# TODO: Check location is correct
						location = ['CCO-CYTOSOL']
						comp.addReactant(ecocycToFeistId[frameId], stoich, location)
					elif proCompDict.has_key(frameId):
						location = proCompDict[frameId].composition['product'][frameId]['compartment']
						comp.addReactant(frameId, stoich, location)
					elif monomerCompartment.has_key(frameId):
						location = monomerCompartment[frameId]
						comp.addReactant(frameId, stoich, location)
					elif smallMolecProCompDict.has_key(frameId):
						location = smallMolecProCompDict[frameId].composition['product'][frameId]['compartment']
						comp.addReactant(frameId, stoich, location)
					else:
						foundAllComponents = False
						s = 'Did not create a small molecule-protein complex for ' + comp.frameId + ' could not find ' + frameId
						writeOut(s, logFile)

				if foundAllComponents:
					comp.addProduct(comp.frameId, 1)
					comp.calculateLocation()
					comp.buildStringComposition(compartmentAbbrev)

					smallMolecProCompDict[comp.frameId] = comp

	# Deal with small molecule-protein complexes that have other small-molecule protein as subunits
	# BUT NOT ones with RNA-protein
	pcSubunit = sets.Set(hasProteinComplexSubunit)
	rpcSubunit = sets.Set(hasRnaProteinComplexSubunit)
	smpcSubunit = sets.Set(hasSmallMolecProteinComplexSubunit)

	if len(smpcSubunit.intersection(pcSubunit)) != 0 or len(smpcSubunit.intersection(rpcSubunit)) != 0:
		# Make sure that rna-protein complexes and small molecule-protein complexes do not depend on protein-protein complexes
		ipdb.set_trace()

	

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['frameId', 'Name', 'Location', 'Composition', 'Composition', 'Formation process', 'Comments'])
		
		keys = proCompDict.keys()
		keys.sort()
		for key in keys:
			c = proCompDict[key]
			csvwriter.writerow([c.frameId, c.name, json.dumps(c.composition['product'][c.frameId]['compartment']), c.compositionString, json.dumps(c.composition), c.formationProcess])

	logFile.close()
	ipdb.set_trace()


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

# Define data type classes
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

class proteinMonomer:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = []
		self.gene = None
		self.modifiedForm = None
		self.comments = ''

class rna:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = []
		self.gene = None
		self.modifiedForm = None
		self.comments = ''

class proteinComplex:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = []
		self.composition = {'reactant' : {}, 'product' : {}}
		self.compositionString = ''
		self.formationProcess = 'Complexation'

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
