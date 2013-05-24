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
					coordinate = int(row[4])
					length = int(row[5])
					pM.direction = row[6]

				if pM.direction == 'forward':
					pM.left = coordinate
					pM.right = coordinate + length
				elif pM.direction == 'right':
					pM.left = coordinate - length
					pM.right = coordinate

				if pM.direction == 'forward':
					pM.sequence = sequence[pM.left - 1: pM.right].transcribe().translate(table = 11)
				elif pM.direction == 'reverse':
					pM.sequence = sequence[pM.left - 1: pM.right].reverse_complement()#.transcribe()#.translate(table = 11)

				if row[4] != '':
					print pM.frameId + '\t\t' + pM.direction + '\t\t' + pM.sequence
					proteinMonomerDict[pM.frameId] = pM
	ipdb.set_trace()


def loadSequence():
	seq_record = SeqIO.read(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'sequence.txt'), "fasta", IUPAC.ambiguous_dna)
	return seq_record.seq

class proteinMonomer():
	def __init__(self):
		self.frameId = ''
		self.sequence = ''
		self.direction = ''
		self.left = 0
		self.right = 0
	
if __name__ == "__main__":
    main()