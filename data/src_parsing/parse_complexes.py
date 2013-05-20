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

		self.geneProdLocalDict = {}
		self.compartmentDict = {}

		self.loadCompartmentData()
		self.loadMonomerData()
		self.loadComplexData()
		self.writeComplexesCSV()

	def loadCompartmentData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'compartments.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			firstRow = True
			for row in csvreader:
				if firstRow:
					firstRow = False
				else:
					compartmentId = row[0]
					compartmentAbbrev = row[1]
					self.compartmentDict[compartmentId] = compartmentAbbrev

	def loadMonomerData(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			firstRow = True

			for row in csvreader:
				if firstRow:
					firstRow = False
				else:
					product = row[10]
					compartment = json.loads(row[9])
					self.geneProdLocalDict[product] = compartment

	def loadComplexData(self):
		rowsToDo = []

		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')

			for row in csvreader:
				comp = proteinComplex()
				comp.frameId = row[0]
				comp.name = re.sub('<[^<]+?>', '', row[1])

				foundAllSubunits = True
				components = row[2][2:-2].split(') (')
				if components[0] != '':
					for c in components:
						info = c.split(', ')
						frameId = info[0]
						stoich = int(info[1])
						if self.geneProdLocalDict.has_key(frameId):
							location = self.geneProdLocalDict[frameId]
							comp.addReactant(frameId, stoich, location)
						elif row[3] == '':
							rowsToDo.append(row)
							foundAllSubunits = False
							break

				if foundAllSubunits:
					comp.addProduct(comp.frameId, 1)
					comp.buildStringComposition()
					if row[3] == '':
						comp.calculateLocation()
					else:
						comp.location = row[3][1:-1].split(' ')

					self.complexDict[comp.frameId] = comp
		#ipdb.set_trace()

	def writeComplexesCSV(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'complexes.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			csvwriter.writerow(['frameId', 'Name', 'Location', 'Composition', 'Composition', 'Formation process', 'Comments'])
			
			keys = self.complexDict.keys()
			keys.sort()
			for key in keys:
				c = self.complexDict[key]
				csvwriter.writerow([c.frameId, c.name, json.dumps(c.location), c.compositionString, json.dumps(c.composition), c.formationProcess])



class proteinComplex:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.location = None
		self.composition = {'reactant' : {}, 'product' : {}}
		self.compositionString = ''
		self.formationProcess = 'Complexation'

	def addReactant(self, name, stoich, location):
		self.composition['reactant'][name] = {'stoichiometry' : None, 'compartment' : None}
		self.composition['reactant'][name]['stoichiometry'] = stoich
		self.composition['reactant'][name]['compartment'] = location

	def addProduct(self, name, stoich):
		self.composition['product'][name] = {'stoichiometry' : None, 'compartment' : None}
		self.composition['product'][name]['stoichiometry'] = stoich
		self.composition['product'][name]['compartment'] = ''

	def buildStringComposition(self):
		s = ''
		subComp = self.composition['reactant'].keys()

		for i in range(len(subComp)):
			c = subComp[i]
			if self.composition['reactant'][c]['stoichiometry'] == 1:
				s += c + ' '
			else:
				stoich = self.composition['reactant'][c]['stoichiometry']
				s += '(' + str(stoich) + ') ' + c + ' '
			if i != len(subComp) - 1:
				s += '+ '
			else:
				s += '==> ' + self.frameId
		self.compositionString = s

	def calculateLocation(self):
		location = []
		for reactantId in self.composition['reactant'].iterkeys():
			location.append(self.composition['reactant'][reactantId]['compartment'])

		locationSet = sets.Set(location)
		locationList = []
		for item in locationSet:
			locationList.append(item)
		self.location = locationList


