#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np
import json
import re
import urllib

class parse_metabolites:
	def __init__(self):
		self.metDict = {}
		self.equivDict = {}
		self.compartmentDict = {}

		self.defineCompartments()
		self.loadOrthData()
		self.loadEcocycData()

		self.writeMetaboliteCSV()


	def defineCompartments(self):
		self.compartmentDict['Cytosol'] = 'CCO-CYTOSOL'
		self.compartmentDict['Periplasm'] = 'CCO-PERI-BAC'
		self.compartmentDict['Extra-organism'] = 'CCO-EXTRACELLULAR'

	def kegg2eco(self, query):
		h = urllib.urlopen("http://ecocyc.org/ECOLI/ajax-frame-search?type=COMPOUND&max=2000&object=%s" % query)
		try:
			return json.loads(h.read())["Results"][0]["id"]
		except TypeError, e:
			return None

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
						met = self.metDict[metId]
						met.compartments[self.compartmentDict[row[5]]] = {'charge' : int(row[4]), 'charged form' : row[3]}
					else:
						met = metabolite()
						met.frameId = self.kegg2eco(row[6])
						met.metId = row[0][:-3]
						met.name = row[1]
						met.empiricalFormula = row[2]
						met.compartments[self.compartmentDict[row[5]]] = {'charge' : int(row[4]), 'charged form' : row[3]}

					self.metDict[met.metId] = met
					self.equivDict[met.name.lower()] = met.metId
					self.equivDict[met.frameId] = met.metId

					synL = row[8].split('/ ')
					for syn in synL:
						self.equivDict[syn.lower()] = met.metId

	def loadEcocycData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_metabolites.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			usedMetId = []
			notUsedMetId = []
			for row in csvreader:
				iupacName = re.sub('<[^<]+?>', '', row[0]).lower()
				otherNames = re.sub('<[^<]+?>', '', row[1][1:-1]).split(', ')
				otherNames = [name.lower() for name in otherNames]
				frameId = row[2]

				if self.equivDict.has_key(iupacName):
					if self.equivDict[iupacName] not in usedMetId:
						usedMetId.append(self.equivDict[iupacName])
						metId = self.equivDict[iupacName]
						# self.metDict[metId].frameId = row[2]
						self.metDict[metId].SMILES = row[4]
				elif self.equivDict.has_key(frameId):
					if self.equivDict[frameId] not in usedMetId:
						usedMetId.append(self.equivDict[frameId])
						metId = self.equivDict[frameId]
						self.metDict[metId].SMILES = row[4]
				else:
					notUsedMetId.append(iupacName)
					# for name in otherNames:
					# 	if self.equivDict.has_key(name):
					# 		if self.equivDict[name] not in usedMetId:
					# 			usedMetId.append(self.equivDict[name])

			print len(usedMetId)
			ipdb.set_trace()





	def writeMetaboliteCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'metabolites.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['ID', 'Name', 'frameId', 'Empirical formula', 'Compartment data', 'Weight', 'SMILES', 'Hydrophobic', 'Media concentration', 'Maximum exchange rate (mmol/gDCW/h)','Comments'])
			
			keys = self.metDict.keys()
			keys.sort()
			for key in keys:
				m = self.metDict[key]
				csvwriter.writerow([m.metId, m.name, m.frameId, m.empiricalFormula, json.dumps(m.compartments), m.weight, m.SMILES, m.hydrophobic, m.mediaConc, m.biomassConc, m.exchangeRate])



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


