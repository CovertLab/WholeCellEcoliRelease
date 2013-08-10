#!/usr/bin/env python
import os
import csv
import json
import ipdb

def main():
	modforms = []

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes_modified.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			modforms.append(row['Frame ID'])


	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'rna_modified.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			modforms.append(row['Frame ID'])

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinMonomers_modified.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			modforms.append(row['Frame ID'])

	modforms = set(modforms)

	enzymes = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'reactions.csv'),'rb') as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		for row in dictreader:
			enz_nestList = json.loads(row['Enzyme'])
			lookLower(enz_nestList, enzymes)
	enzyems = set(enzymes)

	ipdb.set_trace()


def lookLower(line, enzymes):
	if isinstance(line, list):
		for li in line:
			lookLower(li, enzymes)
	else:
		enzymes.append(line)

if __name__ == "__main__":m
	main()