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

'''

Tasks:

Framework
-determine whether to split unique instances away from MoleculeCounts

Knowledge base
-types of attributes (simple numbers, object ID references, states)
-adding unique instance attrs to KB (can be fake for now)
-list-of-arrays instantiation from KB

Accessing
-attribute type accessors
-indexing rules (adding, removing, extending allocated space)

Querying (partition setup and requests)
-query format
-passing queries
-evaluating queries

Paritioning
-trivial partitioning (only one process requests)
-deferred partitioning (handling by other states)
-request collision handling

Disk
-table creation
-table appending
-table loading

Misc
-build actual Chromosome State
-Replication Process
-update Transcription/Translation (clean up logic, add comments, change to TUs instead of gene/protein pairs)
-Transcripts State (better name?)

'''
