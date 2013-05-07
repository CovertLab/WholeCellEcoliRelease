#!/usr/bin/env python

import csv
import os

class parse_genes(Object):
	def __init__(self):
		self.synDict = {}
		self.synDictFrameId = {}

		self.loadSynDict()
		self.loadHalfLife()


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
		pass



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

class gene(Object):
	def __init__(self):
		self.frameId = None
		self.name = None
		self.symbol = None
		self.type = None
		self.e = None
		self.length = None
		self.direction = None
		self.expression = None
		self.halfLife = None
		self.localization = None
		