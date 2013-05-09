#!/usr/bin/env python

import csv
import os
import ipdb

class parse_genes:
	def __init__(self):
		self.synDict = {}
		self.synDictFrameId = {}

		self.geneDict = {}

		self.loadSynDict()
		self.parseGeneInformation()
		self.loadHalfLife()

		self.writeGeneCSV()

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


	def parseGeneInformation(self):
		unmodifiedForm = {}

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
				if row[4] != '':
					unmodifiedForm[row[1]] = True
				else:
					unmodifiedForm[row[1]] = False

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')
			for row in csvreader:
				allProducts = self.splitBigBracket(row[5])
				for product in allProducts:
					props = self.splitSmallBracket(product)
					if not unmodifiedForm[props['frameId']]:
						newGene = gene()
						newGene.productFrameId = props['frameId']
						newGene.frameId = row[0]
						newGene.symbol = row[1]
						if row[2] != '':
							newGene.coordinate = int(row[2])
							newGene.length = int(row[3]) - int(row[2])
						else:
							newGene.coordinate = None
							newGene.length = None
						newGene.direction = row[4]
						if self.geneDict.has_key(newGene.frameId):
							count = 0
							while self.geneDict.has_key(newGene.frameId + str(count)):
								count += 1
							self.geneDict[newGene.frameId + str(count)] = newGene
						else:
							self.geneDict[newGene.frameId] = newGene

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
					halfLife = float(read_halfLife)

					if self.synDictFrameId.has_key(name[0].lower()):
						geneName = self.synDictFrameId[name[0].lower()]
					elif self.synDictFrameId.has_key(name[1].lower()):
						geneName = self.synDictFrameId[name[1].lower()]
					else:
						break

					self.geneDict[geneName].halfLife = halfLife


	def writeGeneCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'genes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['ID', 'Name', 'Symbol', 'Type', 'Coordinate', 'Length', 'Direction', 'Expression', 'Half life', 'Localization', 'Product','Coments'])

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
