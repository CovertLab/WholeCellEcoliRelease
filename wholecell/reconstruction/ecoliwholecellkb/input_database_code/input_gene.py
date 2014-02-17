import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, GeneAbsolutentPosition, GeneSplices

import re

GENOME_SIZE = 4639675

def read_csv(filen, columns):
	f = open(filen,'r')
	a = {}
	id_a = 1
	header = re.split('\t',f.readline())
	for line in f:
		t = re.split('\t',line)
		a[id_a] = {}
		for i in range(0,columns):			
			a[id_a][header[i]] = t[i]
		id_a = id_a + 1
	
	f.close()
	return a

def validate_genes(genes):
	#validate
	global GENOME_SIZE

	for i in genes:
		#direction
		if (genes[i]['Direction'] != '-' and genes[i]['Direction'] != '+'):
			print 'Invalid input for gene direction'
			exit(0)
		
		#length
		if (int(genes[i]['Length']) <= 0 ):
			print 'Invalid input for gene ',i,', length = ',genes[i]['Length']
			exit(0)

		#coordinate
		if (int(genes[i]['Coordinate']) < 0 or int(genes[i]['Coordinate'])> GENOME_SIZE ):
			print 'Invalid input for gene ',i,', Coordinate = ',genes[i]['Coordinate'], GENOME_SIZE
			exit(0)

		if (genes[i]['Direction'] == '-'):
			if (int(genes[i]['Coordinate'])- int(genes[i]['Length']) < 0 ):
				print 'Invalid input for gene ',i,', direction = -, Coordinate = ',genes[i]['Coordinate']
				exit(0)

		if (genes[i]['Direction'] == '+'):
			if (int(genes[i]['Coordinate']) + int(genes[i]['Length']) > GENOME_SIZE ):
				print 'Invalid input for gene ',i,', direction = +, Coordinate = ',genes[i]['Coordinate']
				exit(0)

def input_genes(): 
	genes = read_csv('/home/wholecell/WholeCellKB_ecoli/public/fixtures/cvs_ecoli/genes.csv',12)	
		
	validate_genes(genes)
		
	#input for other tables	
	gene_diff_position = {} # pr_key, gene_id, abs, old, new
	splices = {} # pr_key, gene_id, splices (boolean value in gene table)
	type_gene = {} #id_type, type
	id_type = 1
	product = {} #id_product, product
	id_product = 1	
	entry_pos_float_data = {}
	id_entry = 1

	for i in genes:
		#splices
		if genes[i]['Splices'] != '[]':
			splices[i] = genes[i]['Splices']
			splices[i] = splices[i].replace('[','')
			splices[i] = splices[i].replace(']','')
			splices[i] = splices[i].replace(' ','')
		#gene_diff_position
		if genes[i]['(absolute nt position, old, new)'] != '[]':
			t1 = genes[i]['(absolute nt position, old, new)']
			t1 = t1.replace('[','')
			t1 = t1.replace(']','')
			t1 = t1.replace('"','')
			t = re.split(',',t1.strip()) 
			gene_diff_position[i] = {}
			gene_diff_position[i]['abs'] = t[0].strip()
			gene_diff_position[i]['old'] = t[1].strip()
			gene_diff_position[i]['new'] = t[2].strip()
	
		#type
		temp = genes[i]['Type']
		flag = -1	
		for x in type_gene:
			if type_gene[x] == temp: #already stored
				flag = x
				genes[i]['Type_id'] = x
				break

		if flag == -1: #new type
			type_gene[id_type] = temp
			genes[i]['Type_id'] = id_type
			id_type = id_type + 1			
			
		#product
		temp = genes[i]['Product']
		flag = -1	
		for x in product:
			if product[x] == temp: #already stored
				flag = x
				genes[i]['Product_id'] = x
				break

		if flag == -1: #new product
			product[id_product] = temp
			genes[i]['Product_id'] = id_product
			id_product = id_product + 1			
		
		#entry_pos_float_data
		temp = genes[i]['Half life (s)']
		flag = -1	
		for x in entry_pos_float_data:
			if entry_pos_float_data[x]['data'] == temp: #already stored
				flag = x
				genes[i]['Entry_id_life'] = x
				break

		if flag == -1: #new entry
			entry_pos_float_data[id_entry] = {}
			entry_pos_float_data[id_entry]['data'] = temp
			entry_pos_float_data[id_entry]['unit'] = 's'
			genes[i]['Entry_id_life'] = id_entry
			id_entry = id_entry + 1			
	
		temp = genes[i]['Expression']
		flag = -1	
		for x in entry_pos_float_data:
			if entry_pos_float_data[x]['data'] == temp: #already stored
				flag = x
				genes[i]['Entry_id_exp'] = x
				break

		if flag == -1: #new entry
			entry_pos_float_data[id_entry] = {}
			entry_pos_float_data[id_entry]['data'] = temp
			entry_pos_float_data[id_entry]['unit'] = '[]'
			genes[i]['Entry_id_exp'] = id_entry
			id_entry = id_entry + 1			
	
	######## input for tables ################
	
	for i in genes:

		genes[i]['Name'] = genes[i]['Name'].replace('&','\&')
		genes[i]['Name'] = genes[i]['Name'].replace(';','\;')
		genes[i]['Name'] = genes[i]['Name'].replace(':',' ')
		genes[i]['Name'] = genes[i]['Name'].replace('\'','')
		genes[i]['Name'] = genes[i]['Name'].replace('\"','')
		genes[i]['Name'] = genes[i]['Name'].replace('[','')
		genes[i]['Name'] = genes[i]['Name'].replace(']','')
		namet = genes[i]['Name']		

		direct = ''
		if (genes[i]['Direction'] == '-'):
			direct = 'r'
		else:
			direct = 'f'

		spl = ''
		if (genes[i]['Splices']) == '[]':
			spl ='0'
		else:
			spl ='1'

		abpos = ''
		if (genes[i]['(absolute nt position, old, new)']) == '[]':
			abpos ='0'
		else:
			abpos ='1'
 
		p = Gene(frame_id=genes[i]['Frame ID'], symbol=genes[i]['Symbol'], name = namet, chromosome_id=1,coordinate=genes[i]['Coordinate'],length=genes[i]['Length'],direction=direct, expression_id=genes[i]['Entry_id_exp'], half_life_id=genes[i]['Entry_id_life'], productname_id=genes[i]['Product_id'], typegene_id=genes[i]['Type_id'], splices= spl,absolute_nt_position=abpos)
		#p.save() uncomment when need to insert table 
	

	for i in splices:
		temp = re.split(',', splices[i])
		p = GeneSplices(gene_id=i, start1=temp[0], stop1=temp[1], start2=temp[2], stop2=temp[3])
		#p.save()
	
	for i in gene_diff_position:
		p = GeneAbsolutentPosition(gene_id=i, abs_nt_pos=gene_diff_position[i]['abs'], old=gene_diff_position[i]['old'], new=gene_diff_position[i]['new'])
		#p.save()


# already take input for the gene table input_genes(); not run it again for same data as some table has some keys with no unique keys 

