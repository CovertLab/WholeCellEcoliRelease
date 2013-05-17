#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np
import json
import re

class parse_complexes:
	def __init__(self):
		self.complexDict = {}

		self.loadComplexData()
		self.writeComplexesCSV()

	def loadComplexData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			for row in csvreader:
				comp = proteinComplex()
				comp.frameId = row[0]
				comp.name = re.sub('<[^<]+?>', '', row[1])

				components = row[2][2:-2].split(') (')
				if components[0] != '':
					for c in components:
						info = c.split(', ')
						frameId = info[0]
						try:
							stoich = int(info[1])
						except:
							ipdb.set_trace()
						comp.addReactant(frameId, stoich)
				comp.addProduct(comp.frameId, 1)
				self.complexDict[comp.frameId] = comp

	def writeComplexesCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'complexes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['frameId', 'Name', 'Composition', 'Formation process', 'Comments'])
			
			keys = self.complexDict.keys()
			keys.sort()
			for key in keys:
				c = self.complexDict[key]
				csvwriter.writerow([c.frameId, c.name, json.dumps(c.composition), c.formationProcess])



class proteinComplex:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.composition = {'reactant' : {}, 'product' : {}}
		self.formationProcess = 'Complexation'

	def addReactant(self, name, stoich):
		self.composition['reactant'][name] = stoich

	def addProduct(self, name, stoich):
		self.composition['product'][name] = stoich
