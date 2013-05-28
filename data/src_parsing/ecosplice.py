import urllib
import re
import os
import csv

def main():

	csvInfile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'genes.csv'),'rb')
	csvreader = csv.reader(csvInfile, delimiter='\t')
	csvOutfile = open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'geneCoordinates.csv'),'wb')
	csvwriter = csv.writer(csvOutfile, delimiter='\t', quotechar='"')

	csvreader.next()
	for row in csvreader:
		geneId = row[0]
		print geneId
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