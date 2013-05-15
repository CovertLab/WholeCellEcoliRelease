#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np

class parse_genes:
	def __init__(self):
		self.synDict = {}
		self.synDictFrameId = {}

		self.geneDict = {}
		self.protLocDict = {}
		self.parameters = {}

		self.loadConstants()
		self.loadSynDict()
		self.parseLocations()
		self.parseGeneInformation()
		self.loadHalfLife()
		self.loadExpression()

		self.writeGeneCSV()

	def loadConstantData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','other_parameters.csv')) as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
				for row in csvreader:
					self.parameters[row[0]] = {'value' : row[1], 'units' : row[2]}


	def loadSynDict(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_gene_synonyms.csv')) as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				name = row[0]
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

				self.synDict[name.lower()] = name.lower()
				self.synDictFrameId[name.lower()] = frameId
				if synList != ['']:
					for syn in synList:
						# Make sure names are not a key and a value
						if self.synDict.has_key(syn.lower()):
							pass
						else:
							self.synDict[syn.lower()] = name.lower()
							self.synDictFrameId[syn.lower()] = frameId

	def parseLocations(self):
		locationList = []
		locationDict = {}

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_locations.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')

			for row in csvreader:
				if row != []:
					locations = self.splitBigBracket(row[0])

					for location in locations:
						param = self.splitSmallBracket(location)
						locationList.append(param['description'])

		locationSet = sets.Set(locationList)
		for item in locationSet:
			if item == 'CCO-RIBOSOME':
				locationDict[item] = 'CCO-CYTOSOL'
			else:
				locationDict[item] = item

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_location.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')

			for row in csvreader:
				protFrameId = row[1]
				locations = self.splitBigBracket(row[2])
				locString = ''
				for loc in locations:
					if loc != '':
						param = self.splitSmallBracket(loc)
						leader = ''
						if locString != '':
							leader = ':'
						locString += leader + locationDict[param['description']]

				self.protLocDict[protFrameId] = locString

	def parseGeneInformation(self):
		unmodifiedForm = {}
		rnaType = {}

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_mod_tree.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')

			for row in csvreader:
				if row[2] != '':
					unmodifiedForm[row[1]] = True
				else:
					unmodifiedForm[row[1]] = False

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_mod_tree.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')

			for row in csvreader:
				if row[3] != '':
					if row[4] != '':
						unmodifiedForm[row[1]] = True
					else:
						unmodifiedForm[row[1]] = False

					if row[1].count('RRNA') > 0:
						rnaType[row[1]] = 'rRNA'
					elif row[1].count('tRNA') > 0 and not unmodifiedForm[row[1]]:
						rnaType[row[1]] = 'tRNA'
					elif not unmodifiedForm[row[1]]:
						rnaType[row[1]] = 'miscRNA'

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')
			for row in csvreader:
				# Check for products. One gene --> one product (even if gene is duplicated in name)
				allProducts = self.splitBigBracket(row[5])
				for product in allProducts:
					props = self.splitSmallBracket(product)
					if not unmodifiedForm[props['frameId']]:
						# If there is no unmodified form then it is the base product of a gene
						# Build new gene
						newGene = gene()
						# Add names
						newGene.productFrameId = props['frameId']
						newGene.name = props['description']
						newGene.frameId = row[0]
						newGene.symbol = row[1]
						# Add locations
						if row[2] != '':
							newGene.coordinate = int(row[2])
							newGene.length = int(row[3]) - int(row[2])
						else:
							newGene.coordinate = None
							newGene.length = None
						# Add direction
						newGene.direction = row[4]
						# Pick new gene name for product if gene name is already used
						# for another valid product
						if self.geneDict.has_key(newGene.frameId):
							count = 0
							while self.geneDict.has_key(newGene.frameId + str(count)):
								count += 1
							self.geneDict[newGene.frameId + str(count)] = newGene
						else:
							self.geneDict[newGene.frameId] = newGene
						# Add RNA type
						if rnaType.has_key(newGene.productFrameId):
							newGene.type = rnaType[newGene.productFrameId]
						else:
							newGene.type = 'mRNA'
						# Add localization
						if self.protLocDict.has_key(newGene.productFrameId):
							newGene.localization = self.protLocDict[newGene.productFrameId]


	def splitBigBracket(self, s):
		s = s[2:-2]
		s = s.replace('"','')
		s = s.split(') (')
		return s


	def splitSmallBracket(self, s):
		s = s.split(', ')
		frameId = s[0]
		description = s[1]
		return {'frameId' : frameId, 'description' : description}

	def loadHalfLife(self):
		rRNAhl = []
		tRNAhl = []
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
					halfLife = float(read_halfLife)*60 # seconds

					if self.synDictFrameId.has_key(name[0].lower()):
						geneName = self.synDictFrameId[name[0].lower()]
					elif self.synDictFrameId.has_key(name[1].lower()):
						geneName = self.synDictFrameId[name[1].lower()]
					else:
						print 'Gene half life not found ' + name[0] + ' ' + name[1]
						break

					if self.geneDict[geneName].type == 'rRNA':
						rRNAhl.append(self.parameters['rRNA half life']*3600*24)
					elif self.geneDict[geneName].type == 'tRNA':
						tRNAhl.append(self.parameters['tRNA half life']*3600*24)
					else:
						mRNAhl.append(halfLife)

					self.geneDict[geneName].halfLife = halfLife

	def loadExpression(self):
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

					if self.synDictFrameId.has_key(name.lower()):
						geneName = self.synDictFrameId[name.lower()]
					elif self.synDictFrameId.has_key(bnum.lower()):
						geneName = self.synDictFrameId[bnum.lower()]
					else:
						print 'Gene expression not found ' + name + ' ' + bnum
						break

					self.geneDict[geneName].expression = expression

	def writeGeneCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['ID', 'Name', 'Symbol', 'Type', 'Coordinate', 'Length', 'Direction', 'Expression', 'Half life', 'Localization', 'Product','Comments'])

			keys = self.geneDict.keys()
			keys.sort()

			for key in keys:
				g = self.geneDict[key]
				csvwriter.writerow([g.frameId, g.name, g.symbol, g.type, g.coordinate, g.length, g.direction, g.expression, g.halfLife, g.localization, g.productFrameId])


class gene:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.symbol = None
		self.type = None
		self.coordinate = None
		self.length = None
		self.direction = None
		self.expression = None
		self.halfLife = None
		self.localization = None
		self.productFrameId = None
