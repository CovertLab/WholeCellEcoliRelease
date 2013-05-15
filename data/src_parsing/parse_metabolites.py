#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np

class parse_metabolites:
	def __init__(self):
		self.metDict = {}




		self.writeMetaboliteCSV()



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
