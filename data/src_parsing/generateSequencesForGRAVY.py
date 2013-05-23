#!/usr/bin/env python

import os
import csv
import ipdb

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO


def main():
	sequence = loadSequence()

	geneDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_genes.csv'),'rb') as csvfile:
		csvreader = csv.reader(csvfile, delimiter='\t')
		csvreader.next()
		for row in csvreader:
			pass




def loadSequence():
	seq_record = SeqIO.read(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'sequence.txt'), "fasta", IUPAC.unambiguous_dna)
	return seq_record.seq
	
if __name__ == "__main__":
    main()