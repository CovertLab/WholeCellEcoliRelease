import urllib
import re

def ecosplice(query):
	h = urllib.urlopen("http://ecocyc.org/ECOLI/sequence-spliced?type=GENE&object=%s" % query)
	s = h.read()
	m = re.search("(?P<numbers>[0-9]+\.\.[0-9]+) ", s)
	if m != None:
		return m.group("numbers")
	m = re.search("(?P<numbers>[a-z]*\(.*[0-9]+\.\.[0-9]+.*\)) ", s)
	return m.group("numbers")
