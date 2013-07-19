import csv
import os
import json

def buildTurnoverTable():
	enzymeDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed','reactions.csv')) as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		dictreader.next()
		for row in dictreader:
			if row['EC'] != '' and row['Enzyme'] != 'null':
				isozymes = json.loads(row['Enzyme'])
				for iso in isozymes:
					iso = tuple(iso)
					if not enzymeDict.has_key(iso):
						newEnzyme = enzyme()
						enzymeDict[iso] = newEnzyme
					enzymeDict[iso].frameId = iso
					enzymeDict[iso].EC.append(row['EC'])
					enzymeDict[iso].reacID.append(row['Frame ID'])
					enzymeDict[iso].reacStoich.append(row['Stoichiometry (pH 7.2)'])
					enzymeDict[iso].direction.append(row['Direction'])

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'turnover_annotation.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Enzyme Frame ID', 'EC', 'Reaction ID', 'Reaction stoichiometry', 'Direction', 'Turnover Forward (s^-1)', 'Turnover Reverse (s^-1)', 'Comments'])

		keys = enzymeDict.keys()
		keys.sort()

		for key in keys:
			e = enzymeDict[key]
			for i in range(len(e.EC)):
				reverseSet = None
				if e.direction[i] == 'forward only':
					reverseSet = 0
				csvwriter.writerow([json.dumps(e.frameId), e.EC[i], e.reacID[i], e.reacStoich[i], e.direction[i], None, reverseSet])

class enzyme():
	def __init__(self):
		self.frameId = None
		self.EC = []
		self.reacID = []
		self.reacStoich = []
		self.direction = []