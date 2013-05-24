#!/usr/bin/env python

import os
import csv
import ipdb

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO


def main():
	sequence = loadSequence()

	proteinMonomerDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')
		csvreader.next()
		for row in csvreader:
			if row[3] == 'mRNA':
				pM = proteinMonomer()
				pM.frameId = row[7]
				pM.coordinate = int(row[4])
				pM.length = int(row[5])
				pM.direction = row[6]


def loadSequence():
	seq_record = SeqIO.read(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'sequence.txt'), "fasta", IUPAC.unambiguous_dna)
	return seq_record.seq

class proteinMonomer():
	def __init__(self):
		self.frameId = ''
		self.sequence = ''
		self.coordinate = 0
		self.length = 0
		self.direction = ''
	
if __name__ == "__main__":
    main()