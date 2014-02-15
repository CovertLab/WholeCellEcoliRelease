import re
import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Chromosome
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core import validators

GENOME_SIZE = 4639675

def validate_dna_sequence(seq):
	temp = validators.RegexValidator(regex=r'^[ACGT]+$', message='Enter a valid DNA sequence consisting of only the letters A, C, G, and T')(seq)
	
def readFasta(seqfile):
	global GENOME_SIZE

	f = open(seqfile, 'r')	
	data = f.readlines()
	f.close()
	
	wid = None
	sequences = {}
	for i in range(len(data)):
		if data[i][0] == '>':
			wid = data[i].strip()[1:]
			sequences[wid] = ''
		else:
			if wid is None:
				print 'File does not match FASTA format.'
				exit(0)
			sequences[wid] += data[i].strip()	
	
	#validate
	for wid, sequence in sequences.iteritems():
		temp = validate_dna_sequence(sequence)
		#if (temp != 'None'):
		#	print 'Error in DNA sequence\n'
		#	exit(0)
		if len(sequence) != GENOME_SIZE:
			print 'DNA sequence size doesn\'t match with predifined size = ', GENOME_SIZE
			exit(0)

	return sequences		


def input_chromosome(fw):
	global all_tables
	
	sequences = readFasta('/home/wholecell/WholeCellKB_ecoli/public/fixtures/cvs_ecoli/sequence.txt')

	chr_id = 1
	#fw.write('- model: public.molecule\n')
	#fw.write('  pk: '+str(chr_id)+'\n')
	#fw.write('  fields:\n')		
	
	for wid in sequences:
		fw.write('- model: public.chromosome\n')
		fw.write('  pk: '+str(chr_id)+'\n')
		chr_id = chr_id + 1
		fw.write('  fields:\n')
		fw.write('    sequence: '+sequences[wid]+'\n')
		fw.write('    length: '+str(GENOME_SIZE)+'\n')

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

def read_csv_bad_formating(filen, columns):
	f = open(filen,'r')
	a = {}
	id_a = 1
	header = re.split('\t',f.readline())
	for line in f:
		t = re.split('\t',line)
		a[id_a] = {}
		if(len(t) != len(header)):
			print id_a,t
			continue
		for i in range(0,columns):			
			a[id_a][header[i]] = t[i]
		id_a = id_a + 1
	f.close()
	return a

def input_genes(fw): 
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
			if gene_diff_position[i]['old'] == '*':
				gene_diff_position[i]['old'] = '\\*'

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
	'''
    for i in entry_pos_float_data:
		fw.write('- model: public.entrypositivefloatdata\n')
		fw.write('  pk: '+str(i)+'\n')
		fw.write('  fields:\n')
		fw.write('    value: '+entry_pos_float_data[i]['data']+'\n')
		fw.write('    units: '+entry_pos_float_data[i]['unit']+'\n')
	
	for i in type_gene:
		fw.write('- model: public.genetype\n')
		fw.write('  pk: '+str(i)+'\n')
		fw.write('  fields:\n')
		fw.write('    type_gene: '+type_gene[i]+'\n')
	
	for i in product:
		fw.write('- model: public.geneproduct\n')
		fw.write('  pk: '+str(i)+'\n')
		fw.write('  fields:\n')
		fw.write('    product: '+product[i]+'\n')
	
	
	for i in genes:
		fw.write('- model: public.gene\n')
		fw.write('  pk: '+str(i)+'\n') 
		fw.write('  fields:\n')
		fw.write('    frame_id: '+genes[i]['Frame ID']+'\n')
		fw.write('    symbol: '+genes[i]['Symbol']+'\n')

		genes[i]['Name'] = genes[i]['Name'].replace('&','\&')
		genes[i]['Name'] = genes[i]['Name'].replace(';','\;')
		genes[i]['Name'] = genes[i]['Name'].replace(':',' ')
		genes[i]['Name'] = genes[i]['Name'].replace('\'','')
		genes[i]['Name'] = genes[i]['Name'].replace('\"','')
		genes[i]['Name'] = genes[i]['Name'].replace('[','')
		genes[i]['Name'] = genes[i]['Name'].replace(']','')
		fw.write('    name: '+genes[i]['Name']+'\n')

		fw.write('    chromosome_id_refs_id_d8e6c528: 1\n')
		fw.write('    coordinate: '+genes[i]['Coordinate']+'\n')
		fw.write('    length: '+genes[i]['Length']+'\n')

		if (genes[i]['Direction'] == '-'):
			fw.write('    direction: r\n')
		else:
			fw.write('    direction: f\n')
		
		fw.write('    expression_id: '+str(genes[i]['Entry_id_exp'])+'\n')
		fw.write('    half_life_id_refs_id_946e5e57: '+str(genes[i]['Entry_id_life'])+'\n')
		fw.write('    productname_id: '+str(genes[i]['Product_id'])+'\n')
		fw.write('    typegene_id: '+str(genes[i]['Type_id'])+'\n')

		if (genes[i]['Splices']) == '[]':
			fw.write('    splices: 0\n')
		else:
			fw.write('    splices: 1\n')

		if (genes[i]['(absolute nt position, old, new)']) == '[]':
			fw.write('    absolute_nt_position: 0\n')
		else:
			fw.write('    absolute_nt_position: 1\n')
	'''
	x = 1
	for i in splices:
		temp = re.split(',', splices[i])
		fw.write('- model: public.genesplices\n')
		fw.write('  pk: '+str(x)+'\n')
		fw.write('  fields:\n')
		fw.write('    : '+str(i)+'\n')
		fw.write('    start1: '+temp[0]+'\n')
		fw.write('    stop1: '+temp[1]+'\n')
		fw.write('    start2: '+temp[2]+'\n')
		fw.write('    stop2: '+temp[3]+'\n')
		x = x + 1

	x = 1
	for i in gene_diff_position:
		fw.write('- model: public.geneabsolutentposition\n')
		fw.write('  pk: '+str(x)+'\n')
		fw.write('  fields:\n')
		fw.write('    gene_id: '+str(i)+'\n')
		fw.write('    abs_nt_pos: '+gene_diff_position[i]['abs']+'\n')
		fw.write('    old: '+gene_diff_position[i]['old']+'\n')
		fw.write('    new: '+gene_diff_position[i]['new']+'\n')
		x = x + 1


def create_genes(fw):
	# tables order: chromosome, entrypositiveflotedata, type, product, gene, splices,              geneabsolutentposition 	

	#input_chromosome(fw)
	input_genes(fw) 

def create_ymal_file(filename):
	
	fw = open(filename,'w')
	genes = create_genes(fw)
	
	fw.close()

#create_ymal_file('public/fixtures/gene_data.yaml')


