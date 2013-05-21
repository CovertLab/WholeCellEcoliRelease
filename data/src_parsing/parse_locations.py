#!/usr/bin/env python

import csv
import os
import ipdb
import sets
import numpy as np
import json
import re
import copy

class parse_locations:
	def __init__(self):
		self.locationDict = {}
		self.parseLocations()
		self.writeLocationsCSV()
		self.writeLocationsDict()

	def parseLocations(self):
		# Finds unique set of location frameId's in Ecocyc. Creates a dict so that any locaitons
		# read from another file can be translated into the ones used in the model.
		locationList = []
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_locations.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t')

			for row in csvreader:
				if row != []:
					locations = self.splitBigBracket(row[0])

					for location in locations:
						param = self.splitSmallBracket(location)
						locationList.append(param['description'])

		locationSet = sets.Set(locationList)
		for item in locationSet:
			if item == 'CCO-RIBOSOME':
				self.locationDict[item] = 'CCO-CYTOSOL'
			elif item == 'CCO-MIT-LUM':
				self.locationDict[item] = 'CCO-CYTOSOL'
			elif item == 'CCO-MIT-MEM':
				self.locationDict[item] = 'CCO-PM-BAC-NEG'
			elif item == 'CCO-CYTOSKELETON':
				self.locationDict[item] = 'CCO-CYTOSOL'
			elif item == 'CCO-ENVELOPE':
				self.locationDict[item] = 'CCO-OUTER-MEM'
			else:
				self.locationDict[item] = item

	def writeLocationsCSV(self):
		toRemove = ['CCO-RIBOSOME', 'CCO-MIT-LUM', 'CCO-MIT-MEM', 'CCO-CYTOSKELETON', 'CCO-ENVELOPE']
		abbrevDict = {'CCO-BAC-NUCLEOID' 		: 'n',
						'CCO-CELL-PROJECTION' 	: 'j',
						'CCO-CW-BAC-NEG' 		: 'w',
						'CCO-CYTOSOL'			: 'c',
						'CCO-EXTRACELLULAR'		: 'e',
						'CCO-MEMBRANE'			: 'm',
						'CCO-OUTER-MEM' 		: 'o',
						'CCO-PERI-BAC'			: 'p',
						'CCO-PILUS'				: 'l',
						'CCO-PM-BAC-NEG'		: 'i'}


		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations.csv'),'wb') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

			toPrintDict = copy.copy(self.locationDict)
			for item in toRemove:
				toPrintDict.pop(item)
			keys = toPrintDict.keys()
			keys.sort()
			csvwriter.writerow(['ID', 'Abbreviation'])
			for key in keys:
				csvwriter.writerow([key, abbrevDict[key]])

	def writeLocationsDict(self):
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'locations_dict.txt'),'wb') as jsonfile:
			jsonfile.write(json.dumps(self.locationDict))

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
