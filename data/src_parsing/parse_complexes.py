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

		self.writeComplexesCSV()


	def writeComplexesCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'complexes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['frameId', 'Name', 'Composition', 'Formation process', 'Comments'])
			
			keys = self.complexDict.keys()
			keys.sort()
			for key in keys:
				m = self.metDict[key]
				#csvwriter.writerow([m.metId, m.name, m.frameId, m.empiricalFormula, json.dumps(m.compartments), m.weight, m.SMILES, m.hydrophobic, m.mediaConc, m.biomassConc, m.exchangeRate])



class complex:
	def __init__(self):
