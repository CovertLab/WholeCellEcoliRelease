import numpy

from KnowledgeBase import KnowledgeBase
kbOpts = {}
kb =  KnowledgeBase(**kbOpts)

from KnowledgeBaseEcoli import KnowledgeBaseEcoli
kbEOpts = {}
kbE =  KnowledgeBaseEcoli(**kbEOpts)

if kb == kbE:
	print 'yes'
else:
	print 'no'
##
if kb.compartments == kbE.compartments:
	print 'compartments okay\n'
else:
	print 'NOT same compartments\n'
##
if kb.compIdToAbbrev == kbE.compIdToAbbrev:
	print 'compIdToAbbrev okay\n'
else:
	print 'NOT same compIdToAbbrev\n'
##
if kb.metabolites == kbE.metabolites:
	print 'metabolites okay\n'
else:
	print 'NOT same metabolites\n'
##
if kb.genomeSeq == kbE.genomeSeq:
	print 'genomeSeq okay\n'
else:
	print 'NOT same genomeSeq\n'
##
flag = 0
x = 0
while x<len(kb.genes):
	if kb.genes[x]  != kbE.genes[x]:
		for i in kb.genes[x]:
			if kb.genes[x][i] != kbE.genes[x][i]:
				if i != 'name':				
					#print i, kb.genes[x][i]
					flag = flag + 1
	x = x + 1
if flag>0:
	print 'NOT same genes\n', flag
else:
	print 'genes okay, although name attribute has some punctuation mismatch\n'
##
flag = 0
if len(kb.rnas)!= len(kbE.rnas):
	print 'LEN NOT MATCH, RNAS\n'

else:
	x = 0
	while x < len(kb.rnas):
		kbid = kb.rnas[x]["id"]
		xi = 0
		while xi < len(kbE.rnas):
			if kbid == kbE.rnas[xi]["id"]:
				break
			xi =xi+1

		for i in kb.rnas[x]:
			if i == 'ntCount':
				if numpy.array_equiv(kb.rnas[x][i],kbE.rnas[xi][i]) == False:
					flag = 1
					break
				else:
					continue

			if i == 'name':
				continue
		
			if i == 'modifiedForms':
				if len(kb.rnas[x][i]) != len(kbE.rnas[xi][i]):
					flag = 2
					break
			
				for j in kb.rnas[x][i]:
					if j not in kbE.rnas[xi][i]:
						flag = 3
						print j	
						break
				continue

			if i == 'halfLife':
				if round(kb.rnas[x][i]) != round(kbE.rnas[xi][i]):
					flag = 4
					break
				else:
					continue

			if kb.rnas[x][i]!= kbE.rnas[xi][i]:
				flag = 5
				break
		
		if flag > 0 :
			print flag, i, x,len(kb.rnas),len(kbE.rnas), kb.rnas[x][i],'\n',kbE.rnas[x][i]
			break		
		x = x +1

if flag > 0 :
	print 'NOT same RNA\n'
else:
	print 'RNA okay, although half life might have floating point error, name have punctuation error, modified ford has different ordered list\n'

##
flag = 0
if len(kb.proteins)!= len(kbE.proteins):
	print 'LEN NOT MATCH, proteins\n',len(kb.proteins), len(kbE.proteins)
	x = 0
	while x < len(kb.proteins):
		y = 0
		while y < len(kbE.proteins):
			if kbE.proteins[y]['id'] == kb.proteins[x]['id']:
				break
			y = y+1
		if y ==  len(kbE.proteins): 
			print kb.proteins[x]['id']
		x = x+1
else:
	x = 0
	while x < len(kb.proteins):

		kbid = kb.proteins[x]["id"]
		xi = 0
		while xi < len(kbE.proteins):
			if kbid == kbE.proteins[xi]["id"]:
				break
			xi =xi+1

		for i in kb.proteins[x]:
			if i == 'ntCount' or i == 'aaCount':
				if numpy.array_equiv(kb.proteins[x][i],kbE.proteins[xi][i]) == False:
					flag = 1
					break
				else:
					continue

			if i == 'modifiedForms':
				if len(kb.proteins[x][i]) != len(kbE.proteins[xi][i]):
					flag = 2
					break
			
				for j in kb.proteins[x][i]:
					if j not in kbE.proteins[xi][i]:
						flag = 3
						print j	
						break
				continue


			if kb.proteins[x][i]!= kbE.proteins[xi][i]:
				flag = 4
				break
		
		if flag > 0:
			print flag, i, '\n', kb.proteins[x],'\n',kbE.proteins[xi]
			break		
		x = x +1

if flag > 0 :
	print 'NOT same proteins\n'
else:
	print 'proteins okay\n'
		
##
if kb.reactions == kbE.reactions:
	print 'reactions okay'
else:
	print 'NOT same reactions'
##
if kb.reactionsExchange == kbE.reactionsExchange:
	print 'reactionsExchange okay'
else:
	print 'NOT same reactionsExchange'


