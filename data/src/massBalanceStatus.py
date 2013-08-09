#!/usr/bin/env python

import re
import urllib2

def massBalanceStatus(rxnId):
	BASE_URL = "http://ecocyc.org/ECOLI/NEW-IMAGE?object="
	addr = "%s%s" % (BASE_URL, rxnId)
	site = urllib2.urlopen(addr)
	html = site.read()
	tmp = re.search("Mass balance status: (?P<status>.*)<", html)
	if tmp == None:
		return "Information not found."
	return tmp.group("status")