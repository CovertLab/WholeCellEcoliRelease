import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, Molecule, Location, Comment, ProteinMonomers, Rna, Metabolite, MetaboliteEquivalentEnzyme, ProteinComplex, ProteinComplexModified, ProteinMonomerModified, RnaModified, MetaboliteBiomass

import re
#import rxnParser as rp
import json

GENOME_SIZE = 4639675

def read_csv(filen, columns):
	f = open(filen,'r')
	a = {}
	id_a = 1
	header = re.split('\t',f.readline())

	for line in f:
		t = re.split('\t',line.strip())
		if len(t) <2:
			if line!='':
				print line
			continue
		a[id_a] = {}
		for i in range(0,columns):			
			if (header[i].strip() =='Comments') and (len(t) < columns):	
				a[id_a][header[i].strip()] = ''
			else:	
				#print i, len(t), id_a
				a[id_a][header[i].strip()] = t[i]
		id_a = id_a + 1
		 
	f.close()
	return a


def input_molecules(filename, column):
	all_data = read_csv(filename,column)

	all_molecules = Molecule.objects.all()
	existing_product = []

	for i in all_molecules:
		existing_product.append(i.product)
	
	for i in all_data:
		if all_data[i]['Frame ID'] in existing_product:
			print all_data[i]['Frame ID'], 'already in molecule'
			exit(0)
	
	for i in all_data:
		p = Molecule(product= all_data[i]['Frame ID'])
		print all_data[i]['Frame ID']
		############donep . s a ve()
	print len(all_data)
	

def input_comments(filename,column,header):
	all_data = read_csv(filename,column)
	new_comments = []

	existing_comments = []

	all_comments = Comment.objects.all()

	for i in all_comments:
		existing_comments.append(i.comment_str)
	
	for i in all_data:
		new = all_data[i][header].strip()
		if (new not in existing_comments) and (new not in new_comments):
			new_comments.append(new)
	
	for i in new_comments:
		p = Comment(comment_str= i)
		print i
		###donep .s  ave()

def input_monomers(filename):

	monomers = read_csv(filename,6)	
	## already in public_molecule table; inserted from geneProduct

	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#gene
	genes = {}
	all_genes = Gene.objects.all()
	for i in all_genes:
		genes[i.frame_id]= i.id

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.location_id]= i.id
			
	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id

	#check foreign key
	for i in monomers:
		frame_id = monomers[i]['Frame ID']	
		if frame_id not in molecules:
			print frame_id, 'not in molecules'
			exit(0) 
		monomers[i]['molecule_id'] = molecules[frame_id]
		##			
		gene_id = monomers[i]['Gene']	
		if gene_id not in genes:
			print gene_id, 'not in genes'
			exit(0) 
		monomers[i]['gene_id'] = genes[gene_id]
		##
		location_id = monomers[i]['Location']
		location_id = location_id.replace('"','')
		location_id = location_id.replace(']','')
		location_id = location_id.replace('[','')	

		if location_id not in locations:
			print location_id, 'not in locations'
			exit(0) 
		monomers[i]['location_id'] = locations[location_id]
		##
		if monomers[i]['Modified form'].strip() == '[]':
			monomers[i]['Modified form'] = '0'
		else:
			monomers[i]['Modified form'] = '1'
		##
		comment_id = monomers[i]['Comments']	
		if comment_id not in comments:
			print comment_id, 'not in comments'
			exit(0) 
		monomers[i]['comment_id'] = comments[comment_id]

	#input into database
	for i in monomers:
		p = ProteinMonomers(frame_id_id = monomers[i]['molecule_id'], name=monomers[i]['Name'], gene_fk_id=monomers[i]['gene_id'], location_fk_id=monomers[i]['location_id'], is_modified=monomers[i]['Modified form'], comment_fk_id=monomers[i]['comment_id'])
		
		print 'frame_id_fk =', monomers[i]['molecule_id'], 'name=',monomers[i]['Name'], 'gene_fk=',monomers[i]['gene_id'], 'location_fk=',monomers[i]['location_id'], 'is_modified=',monomers[i]['Modified form'], 'comment_fk=',monomers[i]['comment_id']
		
		###########already in db p .save()
	print len(monomers)

def input_rna(filename):

	rnas = read_csv(filename,6)	
	## already in public_molecule table; inserted from geneProduct

	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	#gene
	genes = {}
	all_genes = Gene.objects.all()
	for i in all_genes:
		genes[i.frame_id]= i.id

	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.location_id]= i.id
			
	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id

	#check foreign key
	for i in rnas:
		frame_id = rnas[i]['Frame ID']	
		if frame_id not in molecules:
			print frame_id, 'not in molecules'
			exit(0) 
		rnas[i]['molecule_id'] = molecules[frame_id]
		##			
		gene_id = rnas[i]['Gene']	
		if gene_id not in genes:
			print gene_id, 'not in genes'
			exit(0) 
		rnas[i]['gene_id'] = genes[gene_id]
		##
		location_id = rnas[i]['Location']
		location_id = location_id.replace('"','')
		location_id = location_id.replace(']','')
		location_id = location_id.replace('[','')	

		if location_id not in locations:
			print location_id, 'not in locations'
			exit(0) 
		rnas[i]['location_id'] = locations[location_id]
		##
		if rnas[i]['Modified form'].strip() == '[]':
			rnas[i]['Modified form'] = '0'
		else:
			rnas[i]['Modified form'] = '1'
		##
		comment_id = rnas[i]['Comments']	
		if comment_id not in comments:
			print comment_id, 'not in comments'
			exit(0) 
		rnas[i]['comment_id'] = comments[comment_id]

	#input into database
	for i in rnas:
		p = Rna(frame_id_id = rnas[i]['molecule_id'], name=rnas[i]['Name'], gene_fk_id=rnas[i]['gene_id'], location_fk_id=rnas[i]['location_id'], is_modified=rnas[i]['Modified form'], comment_fk_id=rnas[i]['comment_id'])
		
		print 'frame_id_fk =', rnas[i]['molecule_id'], 'name=',rnas[i]['Name'], 'gene_fk=',rnas[i]['gene_id'], 'location_fk=',rnas[i]['location_id'], 'is_modified=',rnas[i]['Modified form'], 'comment_fk=',rnas[i]['comment_id']
		
		##already in db p .save()
	print len(rnas)

def input_metabolites(filename):

	metabolites = read_csv(filename,12)	

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
			
	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id	
	
	#check foreign key and input for other 2 tables
	biomass = {}
	biomass_id = 1
	equ_enzymes = {}
	equ_enzymes_id = 1

	for i in metabolites:
		frame_id = metabolites[i]['Frame ID']
		if frame_id not in molecules:
			print frame_id, 'not in molecules'
			exit(0) 
		metabolites[i]['molecule_id'] = molecules[frame_id]
		##		
		if metabolites[i]['Maximum exchange rate (mmol/gDSW/hr)'].strip() == '':
			metabolites[i]['Maximum exchange rate (mmol/gDSW/hr)'] = 0
		##
		if metabolites[i]['Fake metabolite'].strip() == '':
			metabolites[i]['Fake metabolite'] = False
		##
		comment_id = metabolites[i]['Comments']	
		comment_id = comment_id.replace('"','')
		if comment_id not in comments:
			print comment_id, 'not in comments'
			exit(0) 
		metabolites[i]['comment_id'] = comments[comment_id]
		##
		if metabolites[i]['Media Concentration (mM)'].strip() == '':
			metabolites[i]['Media Concentration (mM)'] = 0		
		##
		if '""core"": [], ""wildtype"": []}' in metabolites[i]['Biomass information']:
			metabolites[i]['is_Biomass'] = False
		else:
			metabolites[i]['is_Biomass'] = True
			temp = metabolites[i]['Biomass information'].strip()
			temp = temp.replace('""','"') 
			temp = temp[1:len(temp)-1]

			json_obj = json.loads(temp)

			biomass_text = 'core'
			if json_obj['core'] == json_obj['wildtype']:
				biomass_text = 'both'
			else:
				print 'diff biomass info:', biomass[biomass_id]['metabolite_id'] 
		
			for biomass_temp in json_obj['core']:
				biomass[biomass_id] = {}
				biomass[biomass_id]['metabolite_id'] = metabolites[i]['molecule_id']
 				
				biomass[biomass_id]['concentration'] = biomass_temp['mmol/gDCW']
				if biomass_temp['location'] not in locations:		
					print 'metabolites location not valid'
					exit(0)
 				biomass[biomass_id]['location'] = locations[biomass_temp['location']]
				biomass[biomass_id]['type'] = biomass_text
				biomass_id = biomass_id + 1

			if biomass_text == 'core':
				
				for biomass_temp in json_obj['wildtype']:
					biomass[biomass_id] = {}
					biomass[biomass_id]['metabolite_id'] = metabolites[i]['molecule_id']
 				
					biomass[biomass_id]['concentration'] = biomass_temp['mmol/gDCW']
					if biomass_temp['location'] not in locations:		
						print 'metabolites location not valid'
						exit(0)
 					biomass[biomass_id]['location'] = locations[biomass_temp['location']]
					biomass[biomass_id]['type'] = 'wildtype'
					biomass_id = biomass_id + 1			
				
		##					
		if metabolites[i]['Equivalent enzyme frameId'].strip() == '[]':
			metabolites[i]['has_equ_enzyme'] = False
		else:
			metabolites[i]['has_equ_enzyme'] = True
	
			temp_enzymes = metabolites[i]['Equivalent enzyme frameId'][2:len(metabolites[i]['Equivalent enzyme frameId'])-2]
			temp_enzymes = temp_enzymes.replace('"','')  
			enzymes = re.split(',',temp_enzymes)

			for j in enzymes:
				j = j.strip()
				enz = j[:len(j)-3]
				l = j[len(j)-2]
				if l not in locations:
					print 'equivalent enzymes location, ',l
					exit(0)
				if enz not in molecules:
					print 'equivalent enzymes enzyme, ',enz	
					exit(0)
				equ_enzymes[equ_enzymes_id] = {}
				equ_enzymes[equ_enzymes_id]['metabolites_id'] =  metabolites[i]['molecule_id']
				equ_enzymes[equ_enzymes_id]['location'] =  locations[l]
				equ_enzymes[equ_enzymes_id]['enzyme'] =  molecules[enz]
				equ_enzymes_id = equ_enzymes_id + 1
	
	#### input databases
	
	###########
	#for i in metabolites:
		#p = Metabolite(metabolite_id_id = metabolites[i]['molecule_id'], name = metabolites[i]['Name'], feist_formula = metabolites[i]['Feist formula'], ph_formula = metabolites[i]['pH 7.2 formula'], ph_charge = metabolites[i]['pH 7.2 charge'] , ph_weight =  metabolites[i]['pH 7.2 Weight'], media_concentration = metabolites[i]['Media Concentration (mM)'], has_biomass = metabolites[i]['is_Biomass'], maximum_exchange_rate =metabolites[i]['Maximum exchange rate (mmol/gDSW/hr)'] , fake_metabolite = metabolites[i]['Fake metabolite'],has_equivalent_enzyme = metabolites[i]['has_equ_enzyme'], comment_fk_id = metabolites[i]['comment_id'])

		#print 'metabolite_id_fk = ',metabolites[i]['molecule_id'], 'name =', metabolites[i]['Name'], 'feist_formula =', metabolites[i]['Feist formula'], 'ph_formula =', metabolites[i]['pH 7.2 formula'], 'ph_charge = ',metabolites[i]['pH 7.2 charge'] , 'ph_weight =  ',metabolites[i]['pH 7.2 Weight'], 'media_concentration = ',metabolites[i]['Media Concentration (mM)'], 'has_biomass = ',metabolites[i]['is_Biomass'],' maximum_exchange_rate =',metabolites[i]['Maximum exchange rate (mmol/gDSW/hr)'] , 'fake_metabolite = ',metabolites[i]['Fake metabolite'],'has_equivalent_enzyme =', metabolites[i]['has_equ_enzyme'], 'comment_fk =', metabolites[i]['comment_id']
		##done p .s a ve ( )

	print len(metabolites)
	
	#########

	#### input metabolites id 
	metab = {}
	all_metab = Metabolite.objects.all()
	for i in all_metab:
		metab[i.metabolite_id_id] = i.id


	#for i in equ_enzymes:
		#p = MetaboliteEquivalentEnzyme(metabolite_id_fk_id = metab[equ_enzymes[i]['metabolites_id']], location_fk_id = equ_enzymes[i]['location'], equivalent_enzyme_id_fk_id = equ_enzymes[i]['enzyme'])
		#print 'metabolite_id_fk = ',metab[equ_enzymes[i]['metabolites_id']], 'location_fk = ',equ_enzymes[i]['location'], 'equivalent_enzyme_id_fk = ',equ_enzymes[i]['enzyme']	

		######done p .s a ve()
	
	for i in biomass:
	
		core = False
		wildtype = False
		if biomass[i]['type'] == 'core':
			core = True
			wildtype = False
		if biomass[i]['type'] == 'wildtype':
			core = False
			wildtype = True
		if biomass[i]['type'] == 'both':
			core = True
			wildtype = True

		#p = MetaboliteBiomass(metabolite_id_fk_id = metab[biomass[i]['metabolite_id']],biomass_concentration = biomass[i]['concentration'],biomass_location_fk_id = biomass[i]['location'], is_core = core, is_wildtype = wildtype)	

		#print 'metabolite_id_fk = ',metab[biomass[i]['metabolite_id']],'biomass_concentration =', biomass[i]['concentration'],'biomass_location_fk = ',biomass[i]['location'],core,wildtype

		#####done p .s a ve()


def input_proteinComplexes(filename):
	proteinComplex = read_csv(filename,8)	

	#check foreign key

	#take input from database
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id
	print len(molecules)
	#location 
	locations = {}
	all_locations = Location.objects.all()
	for i in all_locations:
		locations[i.location_id]= i.id
			
	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id	

	#check
	for i in proteinComplex:
		frame_id = proteinComplex[i]['Frame ID']	
		if frame_id not in molecules:
			print frame_id, 'not in molecules'
			exit(0) 
		proteinComplex[i]['molecule_id'] = molecules[frame_id]
		##
		location_id = proteinComplex[i]['Location']
		location_id = location_id.replace('"','')
		location_id = location_id.replace(']','')
		location_id = location_id.replace('[','')	

		if location_id not in locations:
			print location_id, 'not in locations', locations
			exit(0) 
		proteinComplex[i]['location_id'] = locations[location_id]
		##
		if proteinComplex[i]['Modified form'].strip() == '[]':
			proteinComplex[i]['Modified form'] = False
		else:
			proteinComplex[i]['Modified form'] = True
		##
		comment_id = proteinComplex[i]['Comments']	
		if comment_id not in comments:
			print comment_id, 'not in comments'
			exit(0) 
		proteinComplex[i]['comment_id'] = comments[comment_id]
		##
		reaction_str = proteinComplex[i]['Composition']
        
        x = reaction_str.find('>')
        if (len(reaction_str.strip())-1 <=x):
			reaction,direction = rp.parseLeakReaction(reaction_str.strip())
			proteinComplex[i]['direction'] = direction
        else:
			reaction,direction = rp.parseReaction(reaction_str.strip())	
			proteinComplex[i]['direction'] = direction
		
	for i in proteinComplex:
		p = ProteinComplex(protein_complex_id = proteinComplex[i]['molecule_id'], name = proteinComplex[i]['Name'], location_fk_id = proteinComplex[i]['location_id'], modified_form = proteinComplex[i]['Modified form'], comment_fk_id = proteinComplex[i]['comment_id'], reaction_direction = direction)
		print 'protein_complex_id',  proteinComplex[i]['molecule_id'], 'name =', proteinComplex[i]['Name'], 'location_fk_id',  proteinComplex[i]['location_id'], 'modified_form = ',proteinComplex[i]['Modified form'], 'comment_fk_id =', proteinComplex[i]['comment_id'], 'reaction_direction =', direction
		###3p .sa ve)

	print len(proteinComplex)

def input_modified(filename):

	###modified_csv = read_csv(filename,10)	#proteincomplex_modified
	#modified_csv = read_csv(filename,11)	#rna_modified
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
		locations[i.location_id]= i.id
			
	#comment 
	comments = {}
	all_comments = Comment.objects.all()
	for i in all_comments:
		comments[i.comment_str]= i.id	

	#unmodified
	unmodifie_form = {}

	#all_unmodified_forms = ProteinComplex.objects.all() #proteincomplex_modified
	#all_unmodified_forms = Rna.objects.all() #rna_modified
	all_unmodified_forms = ProteinMonomers.objects.all() #monomer_modified

	for i in all_unmodified_forms:
		#unmodifie_form[i.protein_complex_id] = i.id #proteincomplex_modified
		#unmodifie_form[i.frame_id_id] = i.id #rna_modified
		unmodifie_form[i.frame_id_id] = i.id #monomer_modified

	#check
	for i in modified_csv:
		frame_id = modified_csv[i]['Frame ID']	
		if frame_id not in molecules:
			print frame_id, 'not in molecules'
			exit(0) 
		modified_csv[i]['molecule_id'] = molecules[frame_id]
		##
		unmodified_id = modified_csv[i]['Unmodified Form']	
		if unmodified_id not in molecules:
			print unmodified_id, 'absent'
			exit(0) 
		if molecules[unmodified_id] not in unmodifie_form:
			print unmodified_id, 'absent'
			exit(0) 
		modified_csv[i]['unmodified_id'] = unmodifie_form[molecules[unmodified_id]]		
		##
		location_id = modified_csv[i]['Location']
		location_id = location_id.replace('"','')
		location_id = location_id.replace(']','')
		location_id = location_id.replace('[','')

		if location_id not in locations:
			print location_id, 'not in locations', locations
			exit(0) 
		modified_csv[i]['location_id'] = locations[location_id]
		##
		comment_id = modified_csv[i]['Comments']	
		if comment_id not in comments:
			print comment_id, 'not in comments'
			exit(0) 
		modified_csv[i]['comment_id'] = comments[comment_id]
		##
		
	for i in modified_csv:
		## proteinComplex modified		
		##p = ProteinComplexModified(protein_complex_mod_id = modified_csv[i]['molecule_id'] ,name = modified_csv[i]['Name'], unmodified_protein_complex_fk_id = modified_csv[i]['unmodified_id'], location_fk_id = modified_csv[i]['location_id'],comment_fk_id = modified_csv[i]['comment_id'])
		##print 'protein_complex_mod_id =', modified_csv[i]['molecule_id'] ,'name = ',modified_csv[i]['Name'], 'unmodified_protein_complex_fk_id = ',modified_csv[i]['unmodified_id'], 'location_fk_id = ',modified_csv[i]['location_id'],'comment_fk_id =', modified_csv[i]['comment_id']

		## RNA modified		
		##p = RnaModified(rna_mod_id = modified_csv[i]['molecule_id'] ,name = modified_csv[i]['Name'], unmodified_rna_fk_id = modified_csv[i]['unmodified_id'], location_fk_id = modified_csv[i]['location_id'],comment_fk_id = modified_csv[i]['comment_id'])
		##print 'rna_mod_id =', modified_csv[i]['molecule_id'] ,'name = ',modified_csv[i]['Name'], 'unmodified_rna_fk_id = ',modified_csv[i]['unmodified_id'], 'location_fk_id = ',modified_csv[i]['location_id'],'comment_fk_id =', modified_csv[i]['comment_id']
		
		#monomer modified
		p = ProteinMonomerModified(protein_monomer_mod_id = modified_csv[i]['molecule_id'] ,name = modified_csv[i]['Name'], unmodified_protein_monomer_fk_id = modified_csv[i]['unmodified_id'], location_fk_id = modified_csv[i]['location_id'],comment_fk_id = modified_csv[i]['comment_id'])
		print 'protein_monomer_mod_id =', modified_csv[i]['molecule_id'] ,'name = ',modified_csv[i]['Name'], 'unmodified_protein_monomer_fk_id = ',modified_csv[i]['unmodified_id'], 'location_fk_id = ',modified_csv[i]['location_id'],'comment_fk_id =', modified_csv[i]['comment_id']
	
		##########p . s  a ve()

	print len(modified_csv)

###############################################################################3
#######input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers.csv',6,'Comments')
#######input_monomers(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers.csv')

#######input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna.csv',6,'Comments')
#######input_rna(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna.csv')

##########input_molecules(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes.csv', 8)
#####input_proteinComplexes(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes.csv')

###########input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes_modified.csv',10,'Comments')
###########input_molecules(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes_modified.csv', 10)
###########input_modified(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinComplexes_modified.csv')

##########input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna_modified.csv',11,'Comments')
##########input_molecules(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna_modified.csv', 11)
##########input_modified('/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/rna_modified.csv')

##########input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers_modified.csv',10,'Comments')
##########input_molecules(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers_modified.csv', 10)
##########input_modified(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/proteinMonomers_modified.csv')

####input_comments(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/metabolites.csv',12,'Comments')
###########input_molecules(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/metabolites.csv', 12)
################input_metabolites(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/metabolites.csv')


