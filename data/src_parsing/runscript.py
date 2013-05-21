#!/usr/bin/env python
import os
import json
import csv
import re


def main():
	parseIntermediateFiles()

	parseGenes()

# Intermediate file functions
def parseIntermediateFiles():
	# Load and save gene synonym dictionary
	parseGeneSynonymDictionary()
	parseProteionMonomerLocations()
	parseGeneProductUnmodifiedForm()
	parseRnaTypes()

def parseGeneSynonymDictionary():
	synDict = {}
	synDictFrameId = {}

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Ecocyc_gene_synonyms.csv')) as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
		for row in csvreader:
			name = re.sub('<[^<]+?>', '', row[0])
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

			synDict[name.lower()] = name.lower()
			synDictFrameId[name.lower()] = frameId
			if synList != ['']:
				for syn in synList:
					# Make sure names are not a key and a value
					if synDict.has_key(syn.lower()):
						pass
					else:
						synDict[syn.lower()] = name.lower()
						synDictFrameId[syn.lower()] = frameId

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_name_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDict))

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_frameId_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDictFrameId))

def parseProteionMonomerLocations():
	proteinLocationDict = {}

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_location.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			protFrameId = row[1]
			locationInformation = splitBigBracket(row[2])
			locationList= []
			for loc in locationInformation:
				if loc != '':
					param = splitSmallBracket(loc)
					if param['description'] != []:
						locationList.append(param['description'])

					proteinLocationDict[protFrameId] = locationList

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'proteinLocations.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(proteinLocationDict))

def parseGeneProductUnmodifiedForm():
	unmodifiedForm = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_mod_tree.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[3] != '':
				if row[4] != '':
					unmodifiedForm[row[1]] = True
				else:
					unmodifiedForm[row[1]] = False

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_mod_tree.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[2] != '':
				unmodifiedForm[row[1]] = True
			else:
				unmodifiedForm[row[1]] = False

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(unmodifiedForm))

def parseRnaTypes():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	rnaType = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_mod_tree.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')

		for row in csvreader:
			if row[3] != '':
				if row[1].count('RRNA') > 0:
					rnaType[row[1]] = 'rRNA'
				elif row[1].count('tRNA') > 0 and not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'tRNA'
				elif not unmodifiedForm[row[1]]:
					rnaType[row[1]] = 'miscRNA'

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'rnaTypes.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(rnaType))


# Parse genes
def parseGenes():
	# Load unmodified forms of RNA and proteins
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_product_unmodifiedForm.json'),'rb') as jsonfile:
		unmodifiedForm = json.loads(jsonfile.read())

	# Load RNA types
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'rnaTypes.json'),'rb') as jsonfile:
		rnaType = json.loads(jsonfile.read())

	# Parse basic information, RNA type, and product
	geneDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')
		for row in csvreader:
			# Check for products. One gene --> one product (even if gene is duplicated in name)
			allProducts = splitBigBracket(row[5])
			for product in allProducts:
				props = splitSmallBracket(product)
				if not unmodifiedForm[props['frameId']]:
					# If there is no unmodified form then it is the base product of a gene
					# Build new gene
					newGene = gene()

					# Add names
					newGene.productFrameId = props['frameId']
					newGene.name = props['description']
					newGene.frameId = row[0]
					newGene.symbol = row[1]

					# Add locations
					if row[2] != '':
						newGene.coordinate = int(row[2])
						newGene.length = int(row[3]) - int(row[2])
					else:
						newGene.coordinate = None
						newGene.length = None

					# Add direction
					newGene.direction = row[4]

					# Pick new gene name for product if gene name is already used
					# for another valid product
					if geneDict.has_key(newGene.frameId):
						count = 0
						while geneDict.has_key(newGene.frameId + str(count)):
							count += 1
						geneDict[newGene.frameId + str(count)] = newGene
					else:
						geneDict[newGene.frameId] = newGene

					# Add RNA type
					if rnaType.has_key(newGene.productFrameId):
						newGene.type = rnaType[newGene.productFrameId]
					else:
						newGene.type = 'mRNA'

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['ID', 'Name', 'Symbol', 'Type', 'Coordinate', 'Length', 'Direction', 'Expression', 'Half life', 'Localization', 'Product','Comments'])

		keys = geneDict.keys()
		keys.sort()

		for key in keys:
			g = geneDict[key]
			csvwriter.writerow([g.frameId, g.name, g.symbol, g.type, g.coordinate, g.length, g.direction, "%0.10f" % g.expression, g.halfLife, json.dumps(g.localization), g.productFrameId, g.comments])



def parseProteinMonomers():
	pass
	# # Add localization
	# if self.protLocDict.has_key(newGene.productFrameId):
	# 	newGene.localization = self.protLocDict[newGene.productFrameId]
	# else:
	# 	newGene.localization = ['CCO-CYTOSOL']
	# 	newGene.comments = 'No localization known. Assume CO-CYTOSOL.\n'


# Utility functions
def splitBigBracket(s):
	s = s[2:-2]
	s = s.replace('"','')
	s = s.split(') (')
	return s

def splitSmallBracket(s):
	s = s.split(', ')
	frameId = s[0]
	description = s[1]
	return {'frameId' : frameId, 'description' : description}

# Define data type classes
class gene:
	def __init__(self):
		self.frameId = None
		self.name = None
		self.symbol = None
		self.type = None
		self.coordinate = None
		self.length = None
		self.direction = None
		self.expression = 0.
		self.halfLife = None
		self.localization = None
		self.productFrameId = None
		self.comments = None

if __name__ == "__main__":
    main()