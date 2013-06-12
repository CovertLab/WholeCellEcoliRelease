import urllib
import re
import os
import csv

def main():

	csvInFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'rb')
	csvreader = csv.reader(csvInFile, delimiter='\t')
	csvOutFile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'geneCoordinates.csv'),'wb')
	csvwriter = csv.writer(csvOutFile, delimiter='\t', quotechar='"')

	csvreader.next()
	for row in csvreader:
		geneId = row[0]
		print geneId
		geneInfo = ecosplice(geneId)
		# If ecocyc puts up one of it's "Pathway Tools Tip" divs, the regular expressions in ecosplice() don't work
		# One way to check if this has occurred is to see if the geneId is in the returned string
		while geneId in geneInfo:
			geneInfo = ecosplice(geneId)
		csvwriter.writerow([geneId, geneInfo])
		print geneInfo

	csvOutFile.close()
	csvInFile.close()

def ecosplice(query):
	h = urllib.urlopen("http://ecocyc.org/ECOLI/sequence-spliced?type=GENE&object=%s" % query)
	s = h.read()
	m = re.search("(?P<numbers>[0-9]+\.\.[0-9]+) ", s)
	if m != None:
		return m.group("numbers")
	m = re.search("(?P<numbers>[a-z]*\(.*[0-9]+\.\.[0-9]+.*\)) ", s)
	return m.group("numbers")

if __name__ == "__main__":
    main()
