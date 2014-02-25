'''

Prototype for unique instances refactor

'''

import numpy as np

KB = {
	'DNA polymerase':{'chromosomeIndex':'uint32', 'chromosomeBoundAt':'uint32'},
	'RNA polymerase':{'chromosomeIndex':'uint32', 'chromosomeBoundAt':'uint32', 'transcriptsIndex':'uint32', 'transcriptsBoundAt':'uint32'},
	'Ribosome':{'transcriptsIndex':'uint32', 'transcriptsBoundAt':'uint32'},
	'mRNA transcript':{'startsAt':'uint32', 'extent':'uint32'}
	}

ALLOCATION = 1000 # for now, allocate a large number.  eventually only allocate a few and extend as needed

uniqueMolAttrs = {}

for molecule, attributes in KB.viewitems():
	uniqueMolAttrs[molecule] = np.zeros(
		ALLOCATION,
		dtype = [(attribute, dataType) for attribute, dataType in attributes.viewitems()]
		)

import ipdb
ipdb.set_trace()
