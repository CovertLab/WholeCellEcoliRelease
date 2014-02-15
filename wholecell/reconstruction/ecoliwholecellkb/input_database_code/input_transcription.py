import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, GeneAbsolutentPosition, GeneSplices, Promoter, Terminator, TranscriptionUnit, TranscriptionUnitGene, TranscriptionUnitTerminator, Location

import re

GENOME_SIZE = 4639675

def my_split(string_to_split,ignorelist,splitchar):
	
	for i in ignorelist:
		string_to_split = string_to_split.replace(i,'')
	
	temp = re.split(splitchar,string_to_split)
	i = 0
	while i < len(temp):
		temp[i] = temp[i].strip()
		i = i + 1			
	return temp

def validate_length(data_dict,header):

	global GENOME_SIZE

	for i in data_dict:
		if ((int(data_dict[i][header]) > GENOME_SIZE) and (int(data_dict[i][header]) < 0) ):
			print 'Invalid input for ',promoter,i,' = ',data_dict[i][header]
			exit(0)
	
def read_csv(filen, columns):
	f = open(filen,'r')
	a = {}
	id_a = 1
	header = re.split('\t',f.readline())
	for line in f:
		t = re.split('\t',line.strip())
		a[id_a] = {}
		for i in range(0,columns):			
			a[id_a][header[i].strip()] = t[i]
		id_a = id_a + 1
	
	f.close()
	return a

class InputTranscriptionUnitData():

	def input_promoter(self,filen):
		promoters = read_csv(filen,6)	
		validate_length(promoters,'Absolute +1 Position')
			
		for i in promoters:
			direct= promoters[i]['Direction']
			if direct == '-':
				direct = 'r'
			else:
				direct = 'f'

			p = Promoter(promoter_id=promoters[i]['Frame ID'], name=promoters[i]['Name'], position=promoters[i]['Absolute +1 Position'], direction=direct) 

			#print 'promoter_id=',promoters[i]['Frame ID'], 'name=',promoters[i]['Name'], 'position=',promoters[i]['Absolute +1 Position'], 'direction=',direct
			#### data already saves p.save()


	def input_terminator(self,filen):
		terminators = read_csv(filen,6)	
		validate_length(terminators,'Left')
		validate_length(terminators,'Right')
			
		for i in terminators:
			rho = '0'
			if terminators[i]['Rho Dependent']=='True':
				rho = '1'

			p = Terminator(terminator_id=terminators[i]['Frame ID'], name=terminators[i]['Name'], left=terminators[i]['Left'], right=terminators[i]['Right'], rho_dependent=rho) 

			#print 'terminator_id=',terminators[i]['Frame ID'], 'name=',terminators[i]['Name'], 'left=',terminators[i]['Left'], 'right=',terminators[i]['Right'], 'rho_dependent=',rho
			####### data already saved p.save()
	
	def input_transcriptionUnit(self,filen):
		trn_unit = read_csv(filen,10)	
		
		for i in trn_unit:
			if trn_unit[i]['Direction'] == '+':
				trn_unit[i]['Direction'] = 'f'
			else:
				trn_unit[i]['Direction'] = 'r'

		#validate
		validate_length(trn_unit,'Left')
		validate_length(trn_unit,'Right')

		#check promoter, terminators, genes foreign key

		#take input from database
		#promoter
		frame_id_promoter = {}
		all_promoter = Promoter.objects.all()
		for i in all_promoter:
			frame_id_promoter[i.promoter_id] = i
		#terminator
		frame_id_terminator = {}
		all_terminators = Terminator.objects.all()
		for i in all_terminators:
			frame_id_terminator[i.terminator_id]= i
		#genes		
		frame_id_gene = {}
		all_genes = Gene.objects.all()
		for i in all_genes:
			frame_id_gene[i.frame_id]= i

		#check for forein key
		for i in trn_unit:
			#promoter
			if trn_unit[i]['Promoter'] not in frame_id_promoter:
				print 'foreign key is missing for promoter', trn_unit[i]['Promoter']
				exit(0)

			#check treminator
			temp_terminator = my_split(trn_unit[i]['Terminators'],'\"[]',',')
			for t in temp_terminator:
				if t not in frame_id_terminator:
					print 'foreign key is missing for terminator', t
					exit(0)
				
			#check gene foreign key and direction
			temp_gene = my_split(trn_unit[i]['Genes'],'\"[]',',')
			for g in temp_gene:
				if g not in frame_id_gene:
					print 'foreign key is missing for gene', g
					exit(0)					
				for j in all_genes:					
					if j.frame_id == g:
						if trn_unit[i]['Direction'] != j.direction:
							print 'gene direction and trn_unit direction dont match ', g,i
							exit(0)	
						

		#input for transcriptional unit
		for i in trn_unit:
			p = TranscriptionUnit(transcription_unit_id=trn_unit[i]['Frame ID'], name=trn_unit[i]['Name'], left=trn_unit[i]['Left'], right=trn_unit[i]['Right'], direction=trn_unit[i]['Direction'], promoter_id_fk =frame_id_promoter[trn_unit[i]['Promoter']],degradation_rate=trn_unit[i]['Degradation rate (1/min)'],expression_rate=trn_unit[i]['Expression rate (a.u./min)'])
 
			##print 'transcription_unit_id=',trn_unit[i]['Frame ID'],' name=',trn_unit[i]['Name'],' left=',trn_unit[i]['Left'],' right=',trn_unit[i]['Right'],', direction=',trn_unit[i]['Direction'],', promoter_id_fk =',frame_id_promoter[trn_unit[i]['Promoter']],'degradation_rate=',trn_unit[i]['Degradation rate (1/min)'],'expression_rate=',trn_unit[i]['Expression rate (a.u./min)']
			## already saved data p.save()


		print 'Done transcriptional unit input\n*********************************************'
		
		#input for TranscriptionUnitTerminator
		all_trnunit = TranscriptionUnit.objects.all()
		frame_id_trnunit = {}
		
		for i in all_trnunit:
			frame_id_trnunit[i.transcription_unit_id]= i

		for i in trn_unit:
			temp_terminator = my_split(trn_unit[i]['Terminators'],'\"[]',',')
			for t in temp_terminator:
				p = TranscriptionUnitTerminator(transcription_unit_id_fk=frame_id_trnunit[trn_unit[i]['Frame ID']], terminator_id_fk=frame_id_terminator[t])
				#print 'transcription_unit_id_fk=',frame_id_trnunit[trn_unit[i]['Frame ID']], 'terminator_id_fk=',frame_id_terminator[t]
				## already saved data p.save()............		
			

		print 'Done transcriptional unit terminator input \n**********************************'	    

		#input for TranscriptionUnitGene
		for i in trn_unit:
			temp_gene = my_split(trn_unit[i]['Genes'],'\"[]',',')
			for t in temp_gene:					
				p = TranscriptionUnitGene(transcriptionunit_id_fk=frame_id_trnunit[trn_unit[i]['Frame ID']], gene_id_fk=frame_id_gene[t])
				# print 'transcriptionunit_id_fk=',frame_id_trnunit[trn_unit[i]['Frame ID']], 'gene_id_fk=',frame_id_gene[t]
				## already save data p.save()..........		

		print 'Done transcriptional unit GENE input \n**********************************'


####################################### done transcriptional Unit ###################################

class InputLocation():

	def input_location_data(self,filen):
	
		location = read_csv(filen,2)	
			
		for i in location:
			p = Location(location_id=location[i]['Frame ID'], abbreviation=location[i]['Abbreviation'])

			##print 'location_id=',location[i]['Frame ID'],' abbreviation=',location[i]['Abbreviation']
			#### data already saves p.save()............



	

#####################################################################################################
################################      Call functions to input data ##################################
#####################################################################################################

############################## Done this tables #####################################################

## x = InputTranscriptionUnitData()
## already saved data x.input_promoter('/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/promoters.csv')
## already saved data x.input_terminator('/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/terminators.csv')
## already save data x.input_transcriptionUnit('/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/transcriptionUnits_with_rates.csv')

##x = InputLocation()
## already save data x.input_location_data('/home/wholecell/ecoliwholecellkb/public/fixtures/cvs_ecoli/locations.csv')

######################################################################################################
######################################################################################################








	
