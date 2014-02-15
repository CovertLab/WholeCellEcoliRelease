import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, Molecule, Location, Comment, ProteinMonomers, Rna, Metabolite, ProteinComplex, ProteinComplexModified, ProteinMonomerModified, RnaModified, RelationStoichiometry, ProteinComplexReactionRelation,ProteinComplexModifiedReaction,ProteinComplexModReactionRelation, ProteinComplexModReactionEnzyme, ProteinMonomerModifiedReaction, ProteinMonomerModReactionEnzyme, ProteinMonomerModReactionRelation, RnaModifiedReaction, RnaModReactionEnzyme, RnaModifiedReactionRelation, MetaboliteReaction, MetaboliteReactionEnzyme, MetaboliteReactionRelation

import re
import rxnParser as rp

GENOME_SIZE = 4639675

def read_csv(filen, columns):
	f = open(filen,'r')
	a = {}
	id_a = 1
	header = re.split('\t',f.readline())

	for line in f:
		t = re.split('\t',line.strip())
		
		a[id_a] = {}
		for i in range(0,columns):			
			if (header[i].strip() =='Comments') and (len(t) < columns):	
				a[id_a][header[i].strip()] = ''
			else:	
				a[id_a][header[i].strip()] = t[i]
		id_a = id_a + 1
	
	f.close()
	return a


def check_database(relation_sch,r):
	for i in relation_sch:
		if (i.reactant_fk_id == r['molecule_id']) and (i.coefficient == r['coeff']) and (i.location_fk_id == r['location_id']):
			return i.id
	return -1000	

def input_proteinComplex_reaction(filename):

	f = open(filename)
	column = 4
    

	#foreign key check

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.abbreviation]= i.id
	
	#all reactant relation
	relation_sch = RelationStoichiometry.objects.all()

	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id
	
	#ProteinComplex
	protein_complex = {}
	all_proteincomplex = ProteinComplex.objects.all()
	for i in all_proteincomplex:
		protein_complex[i.protein_complex_id] = i.id
	######################################################

	#input from file	
	all_reactions = []	
	total_reaction = 0
    
	f.readline()
	for line in f:
		temp = re.split('\t',line)
		frame_id = temp[0].strip()
		reaction_str = temp[column-1]
        
		x = reaction_str.find('>')
		if (len(reaction_str.strip())-1 <=x):
			reaction,direct = rp.parseLeakReaction(reaction_str.strip())
		else:
			reaction,direct = rp.parseReaction(reaction_str.strip())

		total_reaction = total_reaction + len(reaction)

		for i in reaction:            
			i['molecule_id'] = molecules[i['molecule']]
			i['location_id'] = locations[i['location']]
			i['frame_id'] = protein_complex[molecules[frame_id]]
			all_reactions.append(i)

	print total_reaction, len(all_reactions)
	f.close()

	##################################################
	# input for reaction coef

	reaction_coef = []

	for i in all_reactions:
		db = check_database(relation_sch,i)		
				
		if db < 0:
			temp = {'molecule_id':i['molecule_id'],'coeff':i['coeff'],'location_id':i['location_id']}
			if temp not in	reaction_coef:
				reaction_coef.append(temp)

	print len(reaction_coef)	
	#############################################

	########input database
	#for i in reaction_coef:
		#p = RelationStoichiometry(reactant_fk_id=i['molecule_id'],coefficient=i['coeff'],location_fk_id=i['location_id'])
		#print 'reactant_fk_id=',i['molecule_id'],'coefficient=',i['coeff'],'location_fk_id=',i['location_id']
		
		### never input into database as nothing is unique 
		### done(total = 2812)p .  s  a ve()

	##	
	relation_sch = RelationStoichiometry.objects.all()
	for i in all_reactions:
		sch_id = check_database(relation_sch,i)
		if sch_id < 0:
			print 'Error: Stoichiometry relation absent',i
			exit(0)

		p = ProteinComplexReactionRelation(protein_complex_fk_id = i['frame_id'], reactant_relation_id = sch_id)	
		print 'protein_complex_fk_id =', i['frame_id'], 'reactant_relation_id =', sch_id
		### never input into database again as nothing is unique
		############	###done (total = 2896) p . s  ave()


def input_ReactionsModified(filename):
	
	modified_csv = read_csv(filename,10)	#proteincomplex_modified
	
	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.abbreviation]= i.id
			
	#Modified data
	my_mod = {}
	all_my_mod = ProteinComplexModified.objects.all() # protein complex
	
	for i in all_my_mod:
		my_mod[i.protein_complex_mod_id] = i.id #protein complex
		
	##
	reaction_enzyme = []
	reaction_string = []
	reaction_info = []

	for i in modified_csv:
		if modified_csv[i]['Reaction ID'].strip() == '[]':
			continue
		
		all_rxn = modified_csv[i]['Reaction ID']
		all_rxn = all_rxn[2:len(all_rxn)-2]
		all_rxn = all_rxn.replace('""','')
		temp_rxn = re.split(',', all_rxn)

		######### enzyme input   ################################################
		all_rxn_enz = modified_csv[i]['Reaction Enzyme']
		x = all_rxn_enz.replace(',','')
		x = x.replace('[','')
		x = x.replace(']','')

		if len(x.strip()) > 0:
			all_rxn_enz = all_rxn_enz[2:len(all_rxn_enz)-2]
			all_rxn_enz = all_rxn_enz.replace('""','')
			temp_rxn_enz = re.split(',', all_rxn_enz)
		
			j = 0
			while j< len(temp_rxn):
				x = temp_rxn_enz[j][1:len(temp_rxn_enz[j])-1] 
				x = re.split(',',x)
				for k in x:
					k = k.replace('[','')
					k = k.replace(']','')
					if k == '':
						continue
					temp = {}
					temp['reaction_id'] = temp_rxn[j].strip()
					temp['enz'] = molecules[k.strip()]
					reaction_enzyme.append(temp)
				j = j + 1
		############### done input for reaction enzyme ##########################

	
		all_rxn_str = modified_csv[i]['Reaction']
		all_rxn_str = all_rxn_str[2:len(all_rxn_str)-2]
		all_rxn_str = all_rxn_str.replace('""','')
		temp_rxn_str = re.split(',', all_rxn_str)
		
		j = 0
		while j<len(temp_rxn):
			reaction_str = temp_rxn_str[j].strip()			
			x = reaction_str.find('>')
			if (len(reaction_str.strip())-1 <=x):
				reaction,direct = rp.parseLeakReaction(reaction_str.strip())
			else:
				reaction,direct = rp.parseReaction(reaction_str.strip())
			for k in reaction:
				k['molecule_id'] = molecules[k['molecule']]
				k['location_id'] = locations[k['location']]
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['rxn'] = reaction
			reaction_string.append(temp)
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['direction'] = direct
			temp['frame_id'] = my_mod[molecules[modified_csv[i]['Frame ID']]]
			reaction_info.append(temp)			

			j = j + 1
	#################################
	##################################################
	# input for reaction coef

	reaction_coef = []
	relation_sch = RelationStoichiometry.objects.all()

	for i in reaction_string:
		x = i['rxn']
		for j in x:
			db = check_database(relation_sch,j)		
				
			if db < 0:
				temp = {'molecule_id':j['molecule_id'],'coeff':j['coeff'],'location_id':j['location_id']}
				if temp not in	reaction_coef:
					reaction_coef.append(temp)

	print len(reaction_coef)	
	#############################################
	return
	
	########----------------------------input database--------------###########################

	#for i in reaction_coef:
		#p = RelationStoichiometry(reactant_fk_id=i['molecule_id'],coefficient=i['coeff'],location_fk_id=i['location_id'])
		#print 'reactant_fk_id=',i['molecule_id'],'coefficient=',i['coeff'],'location_fk_id=',i['location_id']		
		### never input into database as nothing is unique 
		#######done(total = 2888-2812)p .   s   a ve()

	##	
	for i in reaction_info:
		p = ProteinComplexModifiedReaction(protein_complex_mod_fk_id = i['frame_id'], reaction_id = i['reaction_id'], reaction_direction = str(i['direction']))
		#print 'protein_complex_mod_fk = ',i['frame_id'], 'reaction_id = ',i['reaction_id'], 'reaction_direction =', str(i['direction'])
		####donep . s a v   e()
	
	#again input from datbase for ids
	relation_sch = RelationStoichiometry.objects.all()
	#
	reaction_modified_protein_complex = {}
	all_reaction_modified_protein_complex = ProteinComplexModifiedReaction.objects.all()		
	for i in all_reaction_modified_protein_complex:
		reaction_modified_protein_complex[i.reaction_id] = i.id

	########input db
	for i in reaction_string:
		frame_id = reaction_modified_protein_complex[i['reaction_id']]
		x = i['rxn']

		for j in x:
			sch_id = check_database(relation_sch,j)
			if sch_id < 0:
				print 'Error: Stoichiometry relation absent',j
				exit(0)

			p = ProteinComplexModReactionRelation(complex_mod_reaction_fk_id = frame_id, reactant_relation_id = sch_id)	
			print 'reaction_id =', frame_id, 'reactant_relation_id =', sch_id
			### never input into database again as nothing is unique
			##########donep .  s  a ve()

	#######
	for i in reaction_enzyme:
		frame_id = reaction_modified_protein_complex[i['reaction_id']]
		p = ProteinComplexModReactionEnzyme(complex_mod_reaction_fk_id = frame_id, reaction_enzyme_fk_id = i['enz'])	
		print 'reaction_id =', frame_id, 'enzyme id=',i['enz']
		### never input into database again as nothing is unique
		###############donep .  s a  v e()

def input_monomerReactionsModified(filename):
	
	modified_csv = read_csv(filename,10)	#monomer_modified
	
	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.abbreviation]= i.id
			
	#Modified data
	my_mod = {}
	all_my_mod = ProteinMonomerModified.objects.all() # monomer
	
	for i in all_my_mod:
		my_mod[i.protein_monomer_mod_id] = i.id #monomer
		
	##
	reaction_enzyme = []
	reaction_string = []
	reaction_info = []

	for i in modified_csv:
		if modified_csv[i]['Reaction ID'].strip() == '[]':
			continue
		
		all_rxn = modified_csv[i]['Reaction ID']
		all_rxn = all_rxn[2:len(all_rxn)-2]
		all_rxn = all_rxn.replace('""','')
		temp_rxn = re.split(',', all_rxn)

		######### enzyme input   ################################################
		all_rxn_enz = modified_csv[i]['Reaction enzyme']
		x = all_rxn_enz.replace(',','')
		x = x.replace('[','')
		x = x.replace(']','')

		if len(x.strip()) > 0:
			all_rxn_enz = all_rxn_enz.replace('"','')
			temp_rxn_enz = re.split('],', all_rxn_enz)

			j = 0
			while j< len(temp_rxn):
				x = temp_rxn_enz[j]#[1:len(temp_rxn_enz[j])-1] 
				x = re.split(',',x)
				
				for k in x:
					k = k.strip()
					k = k.replace('[','')
					k = k.replace(']','')
					if k == '':
						continue
					temp = {}
					temp['reaction_id'] = temp_rxn[j].strip()
					temp['enz'] = molecules[k.strip()]
					reaction_enzyme.append(temp)
				j = j + 1
		############### done input for reaction enzyme ##########################

	
		all_rxn_str = modified_csv[i]['Reaction']
		all_rxn_str = all_rxn_str[2:len(all_rxn_str)-2]
		all_rxn_str = all_rxn_str.replace('""','')
		temp_rxn_str = re.split(',', all_rxn_str)
		
		j = 0
		while j<len(temp_rxn):
			reaction_str = temp_rxn_str[j].strip()			
			x = reaction_str.find('>')
			if (len(reaction_str.strip())-1 <=x):
				reaction,direct = rp.parseLeakReaction(reaction_str.strip())
			else:
				reaction,direct = rp.parseReaction(reaction_str.strip())
			for k in reaction:
				if k['molecule'] not in molecules:
					print modified_csv[i]['Frame ID'],k['molecule']
					continue
				k['molecule_id'] = molecules[k['molecule']]
				k['location_id'] = locations[k['location']]
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['rxn'] = reaction
			reaction_string.append(temp)
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['direction'] = direct
			temp['frame_id'] = my_mod[molecules[modified_csv[i]['Frame ID']]]
			if temp['frame_id'] == 'OX-THIOREDOXIN2-MONOMER':
				print temp
			reaction_info.append(temp)			

			j = j + 1
	#################################
	##################################################
	# input for reaction coef

	reaction_coef = []
	relation_sch = RelationStoichiometry.objects.all()
	flag = 0
	flagin = 0
	for i in reaction_string:
		x = i['rxn']
		for j in x:
			db = check_database(relation_sch,j)		
				
			if db < 0:
				temp = {'molecule_id':j['molecule_id'],'coeff':j['coeff'],'location_id':j['location_id']}
				if temp not in	reaction_coef:
					reaction_coef.append(temp)
				else:
					flag = flag + 1
			else:
				flagin = flagin + 1		
	print len(reaction_info), len(reaction_coef), 'common coef',flag,'already in sch table', flagin	

	#############################################

	########----------------------------input database--------------###########################

	#for i in reaction_coef:
		#p = RelationStoichiometry(reactant_fk_id=i['molecule_id'],coefficient=i['coeff'],location_fk_id=i['location_id'])
		#print 'reactant_fk_id=',i['molecule_id'],'coefficient=',i['coeff'],'location_fk_id=',i['location_id']		
		### never input into database as nothing is unique 
		#######done(total = 6172-5907)p .   s   a ve()
	
	
	mytemp = []
	for i in reaction_info:
		p = ProteinMonomerModifiedReaction(protein_monomer_mod_fk_id = i['frame_id'], reaction_id = i['reaction_id'], reaction_direction = str(i['direction']))
		###print 'protein_monomer_mod_fk = ',i['frame_id'], 'reaction_id = ',i['reaction_id'], 'reaction_direction =', str(i['direction'])
		####donep . s a v   e()

	#again input from datbase for ids
	relation_sch = RelationStoichiometry.objects.all()
	#
	reaction_modified_monomer = {}
	all_reaction_modified_monomer = ProteinMonomerModifiedReaction.objects.all()		
	for i in all_reaction_modified_monomer:
		reaction_modified_monomer[i.reaction_id] = i.id
	
	########input db
	for i in reaction_string:
		frame_id = reaction_modified_monomer[i['reaction_id']]
		x = i['rxn']

		for j in x:
			sch_id = check_database(relation_sch,j)
			if sch_id < 0:
				print 'Error: Stoichiometry relation absent',j
				exit(0)

			p = ProteinMonomerModReactionRelation(reaction_fk_id = frame_id, reactant_relation_id = sch_id)	
			#print 'reaction_id =', frame_id, 'reactant_relation_id =', sch_id
			### never input into database again as nothing is unique
			##########donep .  s  a ve()
	print len(reaction_string)	
	#######
	for i in reaction_enzyme:
		frame_id = reaction_modified_monomer[i['reaction_id']]
		p = ProteinMonomerModReactionEnzyme(reaction_fk_id = frame_id, reaction_enzyme_fk_id = i['enz'])	
		print 'reaction_id =', frame_id, 'enzyme id=',i['enz']
		### never input into database again as nothing is unique
		###############donep .  s a  v e()
		
def input_rnaReactionsModified(filename):
	
	modified_csv = read_csv(filename,9)	#rna_modified
	
	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.abbreviation]= i.id
			
	#Modified data
	my_mod = {}
	all_my_mod = RnaModified.objects.all() # rna

	for i in all_my_mod:
		my_mod[i.rna_mod_id] = i.id #rna
	
	##
	reaction_enzyme = []
	reaction_string = []
	reaction_info = []

	for i in modified_csv:
		if modified_csv[i]['Reaction ID'].strip() == '[]':
			continue
		
		all_rxn = modified_csv[i]['Reaction ID']
		all_rxn = all_rxn[2:len(all_rxn)-2]
		all_rxn = all_rxn.replace('""','')
		temp_rxn = re.split(',', all_rxn)

		######### enzyme input   ################################################
		all_rxn_enz = modified_csv[i]['Reaction Enzyme']
		x = all_rxn_enz.replace(',','')
		x = x.replace('[','')
		x = x.replace(']','')

		if len(x.strip()) > 0:
			all_rxn_enz = all_rxn_enz[2:len(all_rxn_enz)-2]
			all_rxn_enz = all_rxn_enz.replace('""','')
			temp_rxn_enz = re.split(',', all_rxn_enz)
		
			j = 0
			while j< len(temp_rxn):
				x = temp_rxn_enz[j] #for rna				
				x = re.split(',',x)
				for k in x:
					k = k.replace('[','')
					k = k.replace(']','')
					if k == '':
						continue
					temp = {}
					temp['reaction_id'] = temp_rxn[j].strip()
					temp['enz'] = molecules[k.strip()]
					reaction_enzyme.append(temp)
				j = j + 1
		############### done input for reaction enzyme ##########################

	
		#--------------------- only for rna (ec input)------------------------#
		ec = modified_csv[i]['EC']
		ec = ec.replace('"','')
		ec = ec.replace('[','')
		ec = ec.replace(']','')
		temp_ec = re.split(',', ec)
		#---------------------only for rna done ------------------------------#
			
		all_rxn_str = modified_csv[i]['Reaction']
		all_rxn_str = all_rxn_str[2:len(all_rxn_str)-2]
		all_rxn_str = all_rxn_str.replace('""','')
		temp_rxn_str = re.split(',', all_rxn_str)
		
		j = 0
		while j<len(temp_rxn):
			reaction_str = temp_rxn_str[j].strip()			
			x = reaction_str.find('>')
			if (len(reaction_str.strip())-1 <=x):
				reaction,direct = rp.parseLeakReaction(reaction_str.strip())
			else:
				reaction,direct = rp.parseReaction(reaction_str.strip())
			for k in reaction:
				k['molecule_id'] = molecules[k['molecule']]
				k['location_id'] = locations[k['location']]
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['rxn'] = reaction
			reaction_string.append(temp)
			
			temp = {}
			temp['reaction_id'] = temp_rxn[j].strip()
			temp['direction'] = direct
			temp['frame_id'] = my_mod[molecules[modified_csv[i]['Frame ID']]]
			
			#only for rna
			if 'null' not in temp_ec[j]:
				temp['ec'] = temp_ec[j].strip()
			else:
				temp['ec'] = ''

			#done only for RNA		
			reaction_info.append(temp)			

			j = j + 1
	#################################
	##################################################
	# input for reaction coef

	reaction_coef = []
	relation_sch = RelationStoichiometry.objects.all()
	flag = 0
	flagin = 0
	for i in reaction_string:
		x = i['rxn']
		for j in x:
			db = check_database(relation_sch,j)		
				
			if db < 0:
				temp = {'molecule_id':j['molecule_id'],'coeff':j['coeff'],'location_id':j['location_id']}
				if temp not in	reaction_coef:
					reaction_coef.append(temp)
				else:
					flag = flag + 1
			else:
				flagin = flagin + 1		
	print len(reaction_info), len(reaction_coef), 'common coef',flag,'already in sch table', flagin	
	#############################################
	
	
	########----------------------------input database--------------###########################

	#for i in reaction_coef:
		#p = RelationStoichiometry(reactant_fk_id=i['molecule_id'],coefficient=i['coeff'],location_fk_id=i['location_id'])
		#print 'reactant_fk_id=',i['molecule_id'],'coefficient=',i['coeff'],'location_fk_id=',i['location_id']		
		### never input into database as nothing is unique 
		#######done(total = 3085-2888=197)psave)
	##	
	for i in reaction_info:
		p = RnaModifiedReaction(rna_mod_fk_id = i['frame_id'], reaction_id = i['reaction_id'], reaction_direction = str(i['direction']), ec = i['ec'])
		#print 'rna_mod_fk = ',i['frame_id'], 'reaction_id = ',i['reaction_id'], 'reaction_direction =', str(i['direction']), 'ec =', i['ec']
		#####donep . s a   ve() 
	
	#again input from datbase for ids
	relation_sch = RelationStoichiometry.objects.all()
	#
	reaction_modified_rna = {}
	all_reaction_modified_rna = RnaModifiedReaction.objects.all()		
	for i in all_reaction_modified_rna:
		reaction_modified_rna[i.reaction_id] = i.id


	########input db
	for i in reaction_string:
		frame_id = reaction_modified_rna[i['reaction_id']]
		x = i['rxn']

		for j in x:
			sch_id = check_database(relation_sch,j)
			if sch_id < 0:
				print 'Error: Stoichiometry relation absent',j
				exit(0)

			p = RnaModifiedReactionRelation(rna_mod_reaction_fk_id = frame_id, reactant_relation_id = sch_id)	
			print 'reaction_id =', frame_id, 'reactant_relation_id =', sch_id
			### never input into database again as nothing is unique
			####donep . s  av e)

	#######
	for i in reaction_enzyme:
		frame_id = reaction_modified_rna[i['reaction_id']]
		p = RnaModReactionEnzyme(reaction_fk_id = frame_id, reaction_enzyme_fk_id = i['enz'])	
		print 'reaction_id =', frame_id, 'enzyme id=',i['enz']
		### never input into database again as nothing is unique
		########donep. sa v()
	
	
def input_metabolite_reactions(filename):

	mt_file = read_csv(filename,9)	

	#foreign key check

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.abbreviation]= i.id
	
	#all reactant relation
	relation_sch = RelationStoichiometry.objects.all()

	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id
		
	
	#input from file	
	all_reactions = []	
	all_enzymes = []
	total_reaction = 0
    
	for line in mt_file:

		#comment
		mt_file[line]['comment_id'] = comments[mt_file[line]['Comments']]

		#reaction
		reaction_str = mt_file[line]['Stoichiometry (pH 7.2)']
        
		x = reaction_str.find('>')
		if (len(reaction_str.strip())-1 <=x):
			reaction,direct = rp.parseLeakReaction(reaction_str.strip())
		else:
			reaction,direct = rp.parseReaction(reaction_str.strip())

		mt_file[line]['direction'] = direct		
		total_reaction = total_reaction + len(reaction)

		for i in reaction:            
			i['molecule_id'] = molecules[i['molecule']]
			i['location_id'] = locations[i['location']]
			i['frame_id'] = mt_file[line]['Frame ID']
			all_reactions.append(i)
		
		#enzyme
		all_enz = mt_file[line]['Enzyme']
		if 'null' in all_enz:
			continue
		all_enz = all_enz.replace('"','')
		all_enz = all_enz.replace('[','')
		all_enz = all_enz.replace(']','')
		temp_enz = re.split(',',all_enz)

		for i in temp_enz:
			temp = {}
			temp['frame_id'] = mt_file[line]['Frame ID']
			temp['enz'] = molecules[i.strip()]
			all_enzymes.append(temp)
		

	print total_reaction, len(all_reactions), len(all_enzymes)

	##################################################
	# input for reaction coef

	reaction_coef = []

	for i in all_reactions:
		db = check_database(relation_sch,i)		
				
		if db < 0:
			temp = {'molecule_id':i['molecule_id'],'coeff':i['coeff'],'location_id':i['location_id']}
			if temp not in	reaction_coef:
				reaction_coef.append(temp)

	print len(reaction_coef)	
	
	###-----------------------------	insert database ----------------------------------#########
	#for i in reaction_coef:
		#p = RelationStoichiometry(reactant_fk_id=i['molecule_id'],coefficient=i['coeff'],location_fk_id=i['location_id'])
		#print 'reactant_fk_id=',i['molecule_id'],'coefficient=',i['coeff'],'location_fk_id=',i['location_id']		
		### never input into database as nothing is unique 
		#######done(total =5907 -3085)p.save()
	##
	
	for i in mt_file:
		p = MetaboliteReaction(frame_id = mt_file[i]['Frame ID'],name = mt_file[i]['Name'],ec = mt_file[i]['EC'],reaction_direction = mt_file[i]['direction'], upper_bound = mt_file[i]['Upper bound (mmol/gDCW-hr)'], lower_bound = mt_file[i]['Lower bound (mmol/gDCW-hr)'], comment_fk_id = mt_file[i]['comment_id'])
		##print 'frame_id = ',mt_file[i]['Frame ID'],'name =', mt_file[i]['Name'],'ec = ',mt_file[i]['EC'],'reaction_direction = ',mt_file[i]['direction'], 'upper_bound = ',mt_file[i]['Upper bound (mmol/gDCW-hr)'],' lower_bound = ',mt_file[i]['Lower bound (mmol/gDCW-hr)'], 'comment_fk_id =', mt_file[i]['comment_id']
		###donepsa ve)

	
	###----------------------------------------input from db ------------------------------###############
	relation_sch = RelationStoichiometry.objects.all()
	#
	reaction_mt = {}
	all_reaction_mt = MetaboliteReaction.objects.all()		
	for i in all_reaction_mt:
		reaction_mt[i.frame_id] = i.id

	for i in all_reactions:
		frame_id = reaction_mt[i['frame_id']]
		
		sch_id = check_database(relation_sch,i)
		if sch_id < 0:
			print 'Error: Stoichiometry relation absent',i
			continue

		p = MetaboliteReactionRelation(metabolite_reaction_fk_id = frame_id, reactant_relation_id = sch_id)	
		#print 'reaction_id =', frame_id, 'reactant_relation_id =', sch_id
		### never input into database again as nothing is unique
		####donep . s  av e)

	#######
	#############DONE above this#############################
	for i in all_enzymes:
		frame_id = reaction_mt[i['frame_id']]
		p = MetaboliteReactionEnzyme(metabolite_reaction_fk_id = frame_id, enzyme_fk_id = i['enz'])	
		print 'reaction_id =', frame_id, 'enzyme id=',i['enz']
		### never input into database again as nothing is unique
		########donep. sa v()
	print len(all_enzymes)
			


######################################################################################################################
### never input into database as nothing is unique input_proteinComplex_reaction(home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes.csv')

### never input into database as nothing is unique #input_ReactionsModified(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes_modified.csv')

### never input into database as nothing is unique #input_rnaReactionsModified(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna_modified.csv')

### never input into database as nothing is unique #####input_metabolite_reactions(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/reactions.csv')

### never input into database as nothing is unique #####input_monomerReactionsModified(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers_modified.csv')


