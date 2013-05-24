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
				pM.frameId = row[9]
				if row[4] != '':
					pM.coordinate = int(row[4])
					pM.length = int(row[5])
				pM.direction = row[6]

				if pM.direction == 'forward':
					pM.sequence = sequence[pM.coordinate: pM.coordinate + pM.length].transcribe()
				elif pM.direction == 'reverse':
					pM.sequence = sequence[pM.coordinate - pM.length: pM.coordinate].reverse_complement().transcribe()

				print pM.frameId + '\t\t' + pM.sequence.tostring()
				proteinMonomerDict[pM.frameId] = pM
	ipdb.set_trace()


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