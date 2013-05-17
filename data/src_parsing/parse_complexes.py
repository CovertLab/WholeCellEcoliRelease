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

				print row[2][1:-1]






	def writeComplexesCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'complexes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['frameId', 'Name', 'Composition', 'Formation process', 'Comments'])
			
			keys = self.complexDict.keys()
			keys.sort()
			for key in keys:
				m = self.complexDict[key]
				#csvwriter.writerow([m.metId, m.name, m.frameId, m.empiricalFormula, json.dumps(m.compartments), m.weight, m.SMILES, m.hydrophobic, m.mediaConc, m.biomassConc, m.exchangeRate])



class complex:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.composition = None
		self.formationProcess = 'Complexation'
		