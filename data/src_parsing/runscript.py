#!/usr/bin/env python
import os
import json
import csv
import re
import numpy as np
import sets
import ipdb
import generateSequencesForGRAVY as gravy

def main():
	parseIntermediateFiles()

	parseGenes()
	parseLocations()
	parseProteinMonomers()
	parseRna()
	parseProteinComplexes()

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
			if len(row) > 3:
				synList.append(row[3].lower())
			if len(row) > 4:
				synList.append(row[4].lower())

			synDict[name.lower()] = name.lower()
			synDictFrameId[name.lower()] = frameId
			if synList != ['']:
				for syn in synList:
					# Make sure names are not a key and a value
					if synDict.has_key(syn.lower()):
						pass
					else:
						synDict[syn.lower()] = name.lower()
						synDictFrameId[syn.lower()] = frameId

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_name_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDict))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_frameId_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDictFrameId))

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
						if newGene.direction == 'forward':
							newGene.coordinate = int(row[2])
							newGene.length = int(row[3]) - int(row[2])
						else:
							newGene.coordinate = int(row[3])
							newGene.length = int(row[3]) - int(row[2])

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
					print 'Gene half life not found ' + read_name + ' ' + read_bnum

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

		print 'Half lives found for ' + str(len(mRNAhl)) + ' mRNAs'
		print 'Half lives found for ' + str(len(rRNAhl)) + ' rRNAs'
		print 'Half lives found for ' + str(len(tRNAhl)) + ' tRNAs'
		print 'Half lives found for ' + str(len(miscRNAhl)) + ' miscRNAs'

		# Average half lives added for genes that half life was not measured
		mrnaAverage = np.around(np.average(mRNAhl),decimals=2)

		if len(rRNAhl) == 0:
			print 'No rRNA half lives measured. Using parameter.'
			rrnaAverage = parameters['rRNA half life']
		if len(tRNAhl) == 0:
			print 'No tRNA half lives measured. Using parameter.'
			trnaAverage = parameters['tRNA half life']
		if len(miscRNAhl) == 0:
			print 'No miscRNA half lives measured. Using parameter.'
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
					print 'Gene expression not found ' + name + ' ' + bnum

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

	print 'Expression found for ' + str(len(mRNAexp)) + ' mRNAs'
	print 'Expression found for ' + str(len(rRNAexp)) + ' rRNAs'
	print 'Expression found for ' + str(len(tRNAexp)) + ' tRNAs'
	print 'Expression found for ' + str(len(miscRNAexp)) + ' miscRNAs'

	# Average expression for gene types where expression was not measured
	mrnaAverage = np.around(np.average(mRNAexp),decimals=2)

	if len(rRNAexp) == 0:
		print 'No rRNA expression measured.'
		rrnaAverage = mrnaAverage
	else:
		rrnaAverage = np.average(rRNAexp)
	if len(tRNAexp) == 0:
		print 'No tRNA expression measured.'
		trnaAverage = mrnaAverage
	else:
		trnaAverage = np.average(tRNAexp)
	if len(miscRNAexp) == 0:
		print 'No miscRNA expression measured.'
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
		print 'Calculating gravy for all genes'
		gravy.main()
	else:
		print 'Gravy already exists'

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

# Parse Locations
def parseLocations():
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
				print 'Location parsing Lopez Campistrous 2005: No name found for gene ' + geneFrameId

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
			if row[0] == '' and row[1] == '':
				hasData = False

			if hasData:
				name = row[0].lower()
				bnumList = row[1].lower().split('-')

				found = False
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
					print 'Location parsing Han 2011: No name found for ' + name + ' ' + bnum

				if geneToProteinMonomerDict.has_key(geneFrameId):
					proteinMonomerFrameId = geneToProteinMonomerDict[geneFrameId]
				else:
					print 'Location parsing Han 2011: No name found for gene ' + geneFrameId

				location = row[2]
				if proteinMonomerDict[proteinMonomerFrameId].location == []:
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

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinMonomerDict.keys()
		keys.sort()
		csvwriter.writerow(['ID', 'Name', 'Gene', 'Location', 'Modified form', 'Comments'])
		for key in keys:
			pm = proteinMonomerDict[key]
			csvwriter.writerow([pm.frameId, pm.name, pm.gene, json.dumps(pm.location), json.dumps(pm.modifiedForm), pm.comments])

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
def parseProteinComplexes():
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
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_small_molecule_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			smallMolecProteinComplexes.append(row[0])

	# Parse protein complex information
	proCompDict = {}
	hasComplexSubunit = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_complexes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			comp = proteinComplex()
			comp.frameId = row[0]
			comp.name = re.sub('<[^<]+?>', '', row[1])

			foundAllComponents = False
			components = row[2][2:-2].split(') (')
			if row[2] != '':
				for c in components:
					info = c.split(', ')
					frameId = info[0]
					stoich = int(info[1])
					if (frameId in proteinComplexes) or (frameId in rnaProteinComplexes) or (frameId in smallMolecProteinComplexes):
						hasComplexSubunit.append(frameId)
						break
					elif monomerCompartment.has_key(frameId):
						location = monomerCompartment[frameId]
						comp.addReactant(frameId, stoich, location)
						foundAllComponents = True
					else:
						print 'Did not create a complex for ' + comp.frameId

			if foundAllComponents:
				comp.addProduct(comp.frameId, 1)
				comp.calculateLocation()

				proCompDict[comp.frameId] = comp

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['frameId', 'Name', 'Location', 'Composition', 'Composition', 'Formation process', 'Comments'])
		
		keys = proCompDict.keys()
		keys.sort()
		for key in keys:
			c = proCompDict[key]
			csvwriter.writerow([c.frameId, c.name, json.dumps(c.composition['product'][c.frameId]['compartment']), c.compositionString, json.dumps(c.composition), c.formationProcess])

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

		sameLocation = False
		if len(self.composition['product'][self.frameId]['compartment']) == 1:
			sameLocation = True

		if sameLocation:
			s += '[' + compartmentDict[self.composition['product'][self.frameId]['compartment'][0]] + ']: '

		for i in range(len(subComp)):
			c = subComp[i]
			if self.composition['reactant'][c]['stoichiometry'] == 1:
				s += c + ' '
			else:
				stoich = self.composition['reactant'][c]['stoichiometry']
				s += '(' + str(stoich) + ') ' + c + ' '
			if i != len(subComp) - 1:
				s += '+ '
			else:
				s += '==> ' + self.frameId
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