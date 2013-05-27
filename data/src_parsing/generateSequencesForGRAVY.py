#!/usr/bin/env python

import os
import csv
import ipdb

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO


def main():
	sequence = loadSequence()

	# Generate sequence information
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
				elif pM.direction == 'reverse':
					pM.left = coordinate - length
					pM.right = coordinate

				if pM.direction == 'forward':
					pM.ntSequence = sequence[pM.left - 1: pM.right]
					pM.sequence = sequence[pM.left - 1: pM.right].transcribe().translate(table = 11)[:-1]
				elif pM.direction == 'reverse':
					pM.ntSequence = sequence[pM.left - 1: pM.right].reverse_complement()
					pM.sequence = sequence[pM.left - 1: pM.right].reverse_complement().transcribe().translate(table = 11)[:-1]

				if row[4] != '':
					#print pM.frameId + '\t\t' + pM.direction + '\t\t' + pM.sequence
					proteinMonomerDict[pM.frameId] = pM


	badSeq = [actPm.frameId for actPm in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if actPm.sequence.count('*') > 0]
	badSeq.sort()
	
	for bS in badSeq:
		proteinMonomerDict.pop(bS)


	# toCompareProteinMonomerDict = {}
	# handle = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'protseq.fasta'), "rU")
	# for record in SeqIO.parse(handle, "fasta") :
	# 	frameId = record.id.split('|')[-1]
	# 	toCompareProteinMonomerDict[frameId] = record.seq
	# handle.close()

	# notMatchSeq = []
	# for frameId in proteinMonomerDict.iterkeys():
	# 	if not toCompareProteinMonomerDict.has_key(frameId):
	# 		print 'fasta does not exist for ' + frameId

	# 	elif toCompareProteinMonomerDict[frameId] != proteinMonomerDict[frameId].sequence:
	# 		notMatchSeq.append(frameId)

	# notMatchSeq.sort()
	# print notMatchSeq
	# ipdb.set_trace()




	for protId in proteinMonomerDict.iterkeys():
		gravy = calculateGravy(proteinMonomerDict[protId])
		proteinMonomerDict[protId].gravy = gravy
	ipdb.set_trace()

def loadSequence():
	seq_record = SeqIO.read(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'sequence.txt'), "fasta", IUPAC.ambiguous_dna)
	return seq_record.seq

def calculateGravy(pMObj):
	length = pMObj.right - pMObj.left

	# Author(s): Kyte J., Doolittle R.F.
	# Reference: J. Mol. Biol. 157:105-132(1982).
	hydropathyValue = {	'A' : 1.800,
						'R' : -4.500,
						'N' : -3.500,
						'D' : -3.500,
						'C' : 2.500,
						'Q' : -3.500,
						'E' : -3.500,
						'G' : -0.400,
						'H' : -3.200,
						'I' : 4.500,
						'L' : 3.800,
						'K' : -3.900,
						'M' : 1.900,
						'F' : 2.800,
						'P' : -1.600,
						'S' : -0.800,
						'T' : -0.700,
						'W' : -0.900,
						'Y' : -1.300,
						'V' : 4.200}

	value = 0.
	for i in range(len(pMObj.sequence)):
		try:
			value += hydropathyValue[pMObj.sequence[i]]
		except:
			ipdb.set_trace()
	value = value / len(pMObj.sequence)
	return value

class proteinMonomer():
	def __init__(self):
		self.frameId = ''
		self.sequence = ''
		self.ntSequence = ''
		self.direction = ''
		self.left = 0
		self.right = 0
		self.gravy = 0.
	
if __name__ == "__main__":
    main()