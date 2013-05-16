#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np

class parse_metabolites:
	def __init__(self):
		self.metDict = {}
		self.compartmentDict = {}

		self.defineCompartments()
		self.loadOrthData()

		self.writeMetaboliteCSV()


	def defineCompartments(self):
		self.compartmentDict['Cytosol'] == 'CCO-CYTOSOL'
		self.compartmentDict['Periplasm'] == 'CCO-PERI-BAC'
		self.compartmentDict['Extra-organism'] == 'CCO-EXTRACELLULAR'

	def loadOrthData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Orth_metabolites.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			firstRow = True
			for row in csvreader:
				if firstRow:
					firstRow = False
				else:
					metId = row[0][:-2]
					if self.metDict.has_key(metName):
						# Append this badboy
						pass
					else:
						met = metabolite()
						met.metId = row[0][:-2]
						met.name = row[1]
						met.empricalFormula = row[2]
						self.compartments


	def writeMetaboliteCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'metabolites.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['ID', 'frameId', 'Empirical formula', 'Type', 'SMILES', 'Charge', 'Weight', 'Hydophobic', 'Media concentration', 'Maximum exchange rate (mmol/gDCW/h)','Comments'])




class metabolite:
	def __init__(self):
		self.frameId = None
		self.metId = None
		self.name = None
		self.empiricalFormula = None
		self.SMILES = None
		self.charge = None
		self.weight = None
		self.hydrophobic = None
		self.mediaConc = None
		self.biomassConc = None
		self.exchangeRate = None
