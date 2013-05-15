import string
from SOAPpy import WSDL
wsdl = "http://www.brenda-enzymes.org/soap2/brenda.wsdl"
client = WSDL.Proxy(wsdl)
resultString = client.getTurnoverNumber("organism*Escherichia coli")
entries = []
for entry in resultString.split("!"):
	D = {}
	for field in entry.split("#"):
		if field != "":
			m = re.match("(?P<key>.*)\*(?P<value>.*)", field)
			D[m.group("key")] = m.group("value")
	entries.append(D)
