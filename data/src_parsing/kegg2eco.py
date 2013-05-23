import urllib
import json

def kegg2eco(query):
	h = urllib.urlopen("http://ecocyc.org/ECOLI/ajax-frame-search?type=COMPOUND&max=2000&object=%s" % query)
	try:
		return json.loads(h.read())["Results"][0]["id"]
	except TypeError, e:
		return None
