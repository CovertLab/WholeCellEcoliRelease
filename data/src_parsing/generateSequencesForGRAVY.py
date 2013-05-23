#!/usr/bin/env python

import os
import csv
import ipdb

def main():
	sequence = loadSequence()

	geneDict
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')
		csvreader.next()
		for row in csvreader:
			pass


def loadSequence():
	sequence = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'sequence.txt')) as txtfile:
		firstLine = True
		for row in txtfile:
			if firstLine:
				firstLine = False
			else:
				sequence.append(row.strip())
	return "".join(sequence)


	
if __name__ == "__main__":
    main()