import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, Molecule, Location, Comment, ProteinMonomers, Rna, Metabolite, MetaboliteEquivalentEnzyme, ProteinComplex, ProteinComplexModified, ProteinMonomerModified, RnaModified, MetaboliteBiomass

import re

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


def input_metabolites(filename):
	
	#molecule
	molecules = {}
	all_molecules = Molecule.objects.all()
	for i in all_molecules:
		molecules[i.product] = i.id

	##
	metabolites = read_csv(filename,12)	
	x = 0
	y = 0
	biomass = {}
	for i in metabolites:

		if '""core"": [], ""wildtype"": []}' in metabolites[i]['Biomass information']:
			biomass[molecules[metabolites[i]['Frame ID']]] = False
		else:
			biomass[molecules[metabolites[i]['Frame ID']]] = True
	

	
	#### input metabolites id 
	all_metab = Metabolite.objects.all()

	for p in all_metab:
		print biomass[p.metabolite_id_id], p.metabolite_id_id	
		###p.has_biomass = biomass[p.metabolite_id_id]
		#p .   save)
		
###input_metabolites(/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/metabolites.csv')

