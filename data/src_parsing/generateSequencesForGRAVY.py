#!/usr/bin/env python

import os
import csv
import ipdb
import json

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

				if row[10] != '[]':
					pM.splice = json.loads(row[10])

				if pM.splice == None:
					baseSequence = sequence[pM.left - 1: pM.right]
				else:
					baseSequence = Seq('', IUPAC.ambiguous_dna)
					for splice in pM.splice:
						baseSequence += sequence[splice[0] - 1: splice[1]]

				if pM.direction == 'forward':
					pM.ntSequence = baseSequence
					pM.sequence = baseSequence.transcribe().translate(table = 11)[:-1]
				elif pM.direction == 'reverse':
					pM.ntSequence = baseSequence.reverse_complement()
					pM.sequence = baseSequence.reverse_complement().transcribe().translate(table = 11)[:-1]

				if row[11] != '[]':
					subRaw = json.loads(row[11])
					position = subRaw[0]
					before = subRaw[1]
					after = subRaw[2]

					saveSequence = list(pM.sequence.tostring())
					if saveSequence[position - 1] != before:
						ipdb.set_trace()
					else:
						saveSequence[position - 1] = after
					newSequence = ''.join(saveSequence)
					pM.sequence = Seq(newSequence, IUPAC.protein)

				if row[4] != '':
					#print pM.frameId + '\t\t' + pM.direction + '\t\t' + pM.sequence
					proteinMonomerDict[pM.frameId] = pM


	badSeq = [actPm.frameId for actPm in [proteinMonomerDict[pmId] for pmId in proteinMonomerDict.iterkeys()] if actPm.sequence.count('*') > 0]
	badSeq.sort()
	print badSeq
	
	ipdb.set_trace()

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

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'proteinMonomerGravy.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		keys = proteinMonomerDict.keys()
		keys.sort()
		csvwriter.writerow(['ID', 'GRAVY'])
		for key in keys:
			pm = proteinMonomerDict[key]
			csvwriter.writerow([pm.frameId, "%0.10000f" % pm.gravy])

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
		value += hydropathyValue[pMObj.sequence[i]]
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
		self.splice = None

if __name__ == "__main__":
    main()