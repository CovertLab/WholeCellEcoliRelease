#!/usr/bin/env python
import os
import json
import csv
import re


def main():
	# Load and save gene synonym dictionary
	parseGeneSynonymDictionary()




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

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'gene_synonyms.json'),'wb') as jsonfile:
		jsonfile.write(json.dumps(synDict))
		jsonfile.write(json.dumps(synDictFrameId))

if __name__ == "__main__":
    main()