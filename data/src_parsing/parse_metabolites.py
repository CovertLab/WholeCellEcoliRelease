#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np
import json

class parse_metabolites:
	def __init__(self):
		self.metDict = {}
		self.compartmentDict = {}

		self.defineCompartments()
		self.loadOrthData()

		self.writeMetaboliteCSV()


	def defineCompartments(self):
		self.compartmentDict['Cytosol'] = 'CCO-CYTOSOL'
		self.compartmentDict['Periplasm'] = 'CCO-PERI-BAC'
		self.compartmentDict['Extra-organism'] = 'CCO-EXTRACELLULAR'

	def loadOrthData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Orth_metabolites.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			firstRow = True
			for row in csvreader:
				if firstRow:
					firstRow = False
				else:
					metId = row[0][:-3]
					if self.metDict.has_key(metId):
						# Append this badboy
						pass
					else:
						met = metabolite()
						met.metId = row[0][:-3]
						met.name = row[1]
						met.empiricalFormula = row[2]
						met.compartments[self.compartmentDict[row[5]]] = {'charge' : int(row[4]), 'charged form' : row[3]}

					self.metDict[met.metId] = met


	def writeMetaboliteCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'metabolites.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['ID', 'frameId', 'Empirical formula', 'Compartment data', 'Weight', 'SMILES', 'Hydrophobic', 'Media concentration', 'Maximum exchange rate (mmol/gDCW/h)','Comments'])
			
			keys = self.metDict.keys()
			keys.sort()
			for key in keys:
				m = self.metDict[key]
				csvwriter.writerow([m.metId, m.frameId, m.empiricalFormula, json.dumps(m.compartments), m.weight, m.SMILES, m.hydrophobic, m.mediaConc, m.biomassConc, m.exchangeRate])



class metabolite:
	def __init__(self):
		self.frameId = None
		self.metId = None
		self.name = None
		self.empiricalFormula = None
		self.compartments = {}
		self.hydrophobic = None
		self.mediaConc = None
		self.biomassConc = None
		self.exchangeRate = None

		self.SMILES = None
		self.weight = None


