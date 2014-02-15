from django.db import models
from django.core import validators
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.db.models import Model, OneToOneField, CharField, IntegerField, URLField, PositiveIntegerField, FloatField, ForeignKey, BooleanField, SlugField, ManyToManyField, TextField, DateTimeField, options, permalink, SET_NULL, Min, Max

# Create your models here.
''' BEGIN: choices '''

CHOICES_DIRECTION = (
	('f', 'Forward'),
	('r', 'Reverse'),
)


''' BEGIN: validators '''
def validate_dna_sequence(seq):
	validators.RegexValidator(regex=r'^[ACGT]+$', message='Enter a valid DNA sequence consisting of only the letters A, C, G, and T')(seq)



class Chromosome(models.Model):
	#additional fields
	sequence = TextField(blank=True, default='', verbose_name='Sequence', validators=[validate_dna_sequence])
	length = PositiveIntegerField(verbose_name='Length (nt)')

	def get_by_natural_key(self):
		return (self.length)
	def __unicode__(self):
		return self.length
		

		
class EntryPositiveFloatData(models.Model):
	value = FloatField(verbose_name='Value', validators=[validators.MinValueValidator(0)])
	units = CharField(max_length = 255, blank=True, default='', verbose_name='Units')
	def get_by_natural_key(self,value):
		return (self.value)
	def __unicode__(self,value):
		return self.value

class GeneType(models.Model):
	type_gene = CharField(max_length = 40, blank=True, default='', verbose_name='Type of Gene',unique=True)
	def get_by_natural_key(self):
		return (self.type_gene)
	def __unicode__(self):
		return self.type_gene

'''
class MoleculeType(models.Model):
	molecule_type = CharField(max_length = 255, blank=True, default='', verbose_name='Type of Molecules',unique=True)
	def get_by_natural_key(self):
		return (self.molecule_type)
	def __unicode__(self):
		return self.molecule_type
'''

class Molecule(models.Model):
	product = CharField(max_length = 255, blank=True, default='', verbose_name='Molecules',unique=True)
	def get_by_natural_key(self):
		return (self.product)
	def __unicode__(self):
		return self.product

class Gene(models.Model):
	
	#additional fields
	frame_id = CharField(max_length=40, blank=True, default='', verbose_name='frame id',unique=True)
	symbol = CharField(max_length=255, blank=True, default='', verbose_name='Symbol')
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	chromosome = ForeignKey(Chromosome, related_name='genes', verbose_name='Chromosome')
	coordinate = PositiveIntegerField(verbose_name='Coordinate (nt)')
	length = PositiveIntegerField(verbose_name='Length (nt)')
	direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')
	expression = ForeignKey(EntryPositiveFloatData, verbose_name='Relative expression', related_name='+')
	half_life = ForeignKey(EntryPositiveFloatData, verbose_name='Half life', related_name='+')
	productname = ForeignKey(Molecule, related_name='productname', verbose_name='molecule')
	typegene = ForeignKey(GeneType, related_name='typegene', verbose_name='GeneType')
	splices = CharField(max_length=10, blank=True, default='', verbose_name='Splices')
	absolute_nt_position = CharField(max_length=10, blank=True, default='', verbose_name='absolute position, old, new')
	
class GeneSplices(models.Model):
	gene = ForeignKey(Gene, verbose_name='Gene Id')
	start1 = PositiveIntegerField(verbose_name='Start position 1')	
	stop1 = PositiveIntegerField(verbose_name='Stop position 1')	
	start2 = PositiveIntegerField(verbose_name='Start position 2')	
	stop2 = PositiveIntegerField(verbose_name='Stop position 2')	

class GeneAbsolutentPosition(models.Model):
	gene = ForeignKey(Gene, verbose_name='Gene Id')
	abs_nt_pos = PositiveIntegerField(verbose_name='Absolute nt position')
	old = CharField(max_length=10, blank=True, default='', verbose_name='Old nt')
	new = CharField(max_length=10, blank=True, default='', verbose_name='New nt')


###########################################################################################
########################## done all tables related with Gene.csv ##########################
###########################################################################################

class Promoter(models.Model):
	promoter_id = CharField(max_length=40, blank=True, default='', verbose_name='frame id',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	position = PositiveIntegerField(verbose_name='Absolute + 1 position')
	direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')
	#need to create another table for ''''signma''''

class Terminator(models.Model):
	terminator_id = CharField(max_length=40, blank=True, default='', verbose_name='frame id',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	left = PositiveIntegerField(verbose_name='left position')
	right = PositiveIntegerField(verbose_name='right position')
	rho_dependent = CharField(max_length=10, verbose_name='Rho Dependent')

class TranscriptionUnit(models.Model):
	transcription_unit_id = CharField(max_length=40, blank=True, default='', verbose_name='frame id',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	left = PositiveIntegerField(verbose_name='left position')
	right = PositiveIntegerField(verbose_name='right position')
	direction = CharField(max_length=10, choices=CHOICES_DIRECTION, verbose_name='Direction')
	degradation_rate = FloatField(verbose_name='degradation rate (1/min)')
	expression_rate = FloatField(verbose_name='expression rate (a.u./min)')
	promoter_id_fk = ForeignKey(Promoter, verbose_name='Promoter Id',related_name ='+')

class TranscriptionUnitGene(models.Model):
	transcriptionunit_id_fk = ForeignKey(TranscriptionUnit, verbose_name='Transcription unit Id',related_name ='+')
	gene_id_fk = ForeignKey(Gene, verbose_name='Gene Id',related_name ='+')
	
	class Meta:
		unique_together = ('transcriptionunit_id_fk','gene_id_fk')

class TranscriptionUnitTerminator(models.Model):
	transcription_unit_id_fk = ForeignKey(TranscriptionUnit, verbose_name='Transcription unit Id',related_name ='+')
	terminator_id_fk = ForeignKey(Terminator, verbose_name='Gene Id',related_name ='+')
	
	class Meta:
		unique_together = ('transcription_unit_id_fk','terminator_id_fk')

#########################################################################################################

class Location(models.Model):
	location_id = CharField(max_length=255, blank=True, default='', verbose_name='location',unique=True)
	abbreviation = CharField(max_length=20, blank=True, default='', verbose_name='abbreviation', unique=True)

############################## monomer, RNA, protein Complex #############################################

class Comment(models.Model):
	comment_str = CharField(max_length=255, blank=True, default='', verbose_name='Comment',unique=True)

class ProteinMonomers(models.Model):
	frame_id = ForeignKey(Molecule, related_name='monomerframeid', verbose_name='Gene_product',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	gene_fk = ForeignKey(Gene, related_name='monomergenefk', verbose_name='Gene')
	location_fk = ForeignKey(Location, related_name='monomerlocationfk', verbose_name='location_fk')
	is_modified = CharField(max_length=4, blank=True, default='0', verbose_name='Name')
	comment_fk = ForeignKey(Comment, related_name='monomercommentfk', verbose_name='comment_fk')
	
class Rna(models.Model):
	frame_id = ForeignKey(Molecule, related_name='rnaframeid', verbose_name='Gene_product',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	gene_fk = ForeignKey(Gene, related_name='rnagenefk', verbose_name='Gene')
	location_fk = ForeignKey(Location, related_name='rnalocation_fk', verbose_name='location_fk')
	is_modified = CharField(max_length=4, blank=True, default='0', verbose_name='Name')
	comment_fk = ForeignKey(Comment, related_name='rnacommentfk', verbose_name='comment_fk')
	
class Metabolite(models.Model):
	metabolite_id = ForeignKey(Molecule, related_name='metaboliteid', verbose_name='metabolites',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	feist_formula = CharField(max_length=255, blank=True, default='', verbose_name='Natural formula')		
	ph_formula = CharField(max_length=255, blank=True, default='', verbose_name='ph 7.2 formula')
	ph_charge =  IntegerField(verbose_name='ph 7.2 charge')
	ph_weight =  FloatField(verbose_name='ph 7.2 weight')
	media_concentration =  FloatField(verbose_name='media concentration (mM)',default = 0)
	has_biomass = BooleanField(verbose_name='Fake Metabolite',default = False)	
	maximum_exchange_rate = FloatField(verbose_name='maximum exchange rate (mmol/gDSW/hr)',default = 0)
	fake_metabolite = BooleanField(verbose_name='Fake Metabolite',default = False)	
	has_equivalent_enzyme = BooleanField(verbose_name='equivalent_enzyme',default = False)	
	comment_fk = ForeignKey(Comment, related_name='metabolitescommentfk', verbose_name='comment fk')
	
class MetaboliteBiomass(models.Model):
	metabolite_id_fk = ForeignKey(Metabolite, related_name='metaboliteidfk', verbose_name='metabolites')	
	biomass_concentration =  FloatField(verbose_name='biomass concentration (molecules/cell)',default = 0)
	biomass_location_fk = ForeignKey(Location, related_name='biomass_locationfk', verbose_name='location_fk')
	is_core = BooleanField(verbose_name='core',default = False)	
	is_wildtype = BooleanField(verbose_name='wildtype',default = False)	

class MetaboliteEquivalentEnzyme(models.Model):
	metabolite_id_fk = ForeignKey(Metabolite, related_name='metabolite_idfk', verbose_name='metabolites')	
	location_fk = ForeignKey(Location, related_name='equ_enzyme_locationfk', verbose_name='location_fk')
	equivalent_enzyme_id_fk = ForeignKey(Molecule, related_name='equivalent_enzymeid_fk', verbose_name='equivalent enzyme')

	
class ProteinComplex(models.Model):
	protein_complex = ForeignKey(Molecule, related_name='proteincomplex', verbose_name='protein complex',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	location_fk = ForeignKey(Location, related_name='location_protein_complex', verbose_name='location_fk')
	modified_form = BooleanField(verbose_name='Modified Form',default = False)	
	comment_fk = ForeignKey(Comment, related_name='comment_proteincomplex', verbose_name='comment_fk', default = '')
	reaction_direction = CharField(max_length=4, blank=True, default='0', verbose_name='reaction direction')

###############
class RelationStoichiometry(models.Model):
	reactant_fk = ForeignKey(Molecule, related_name='%(app_label)s_%(class)s_related_molecule_id', verbose_name='reactant')
	coefficient =  FloatField(verbose_name='coefficient of reactant')
	location_fk = ForeignKey(Location, related_name='%(app_label)s_%(class)s_related_location_fk_id', verbose_name='location_fk')

class ProteinComplexReactionRelation(models.Model):
	protein_complex_fk = ForeignKey(ProteinComplex, related_name='%(app_label)s_%(class)s_related', verbose_name='protein complex')
	reactant_relation = ForeignKey(RelationStoichiometry, related_name='reactant_relation_stoichiometry', verbose_name='relation stoichiometry')
###############

class ProteinComplexModified(models.Model):
	protein_complex_mod = ForeignKey(Molecule, related_name='protein_complex_mod_id', verbose_name='protein complex mod',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	unmodified_protein_complex_fk = ForeignKey(ProteinComplex, related_name='protein_complex_fk_id', verbose_name='unmodified protein complex')
	location_fk = ForeignKey(Location, related_name='location_proteincomplexfk', verbose_name='location_fk')
	comment_fk = ForeignKey(Comment, related_name='comment_proteincomplexmod', verbose_name='comment_fk', default = '')

class ProteinComplexModifiedReaction(models.Model):
	protein_complex_mod_fk = ForeignKey(ProteinComplexModified, related_name='%(app_label)s_%(class)s_related_protein_complex_mod', verbose_name='protein complex mod')
	reaction_id = CharField(max_length=255, blank=True, default='', verbose_name='Reaction ID' ,unique=True)
	ec = CharField(max_length=80, blank=True, default='', verbose_name='EC')
	reaction_direction = CharField(max_length=4, blank=True, default='0', verbose_name='reaction direction')

class ProteinComplexModReactionEnzyme(models.Model):
	complex_mod_reaction_fk = ForeignKey(ProteinComplexModifiedReaction, related_name='%(class)s_reaction', verbose_name='protein complex_mod')
	reaction_enzyme_fk = ForeignKey(Molecule, related_name='%(app_label)s_%(class)s_related_reaction_enz', verbose_name='reaction_enzyme')

class ProteinComplexModReactionRelation(models.Model):
	complex_mod_reaction_fk = ForeignKey(ProteinComplexModifiedReaction, related_name='%(class)s_reaction', verbose_name='protein complex_mod')
	reactant_relation = ForeignKey(RelationStoichiometry, related_name='%(class)s_sch', verbose_name='stoichiometry')

################

class ProteinMonomerModified(models.Model):
	protein_monomer_mod = ForeignKey(Molecule, related_name='protein_monomer_mod_fk', verbose_name='protein monomer mod',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	unmodified_protein_monomer_fk = ForeignKey(ProteinMonomers, related_name='monomer_id', verbose_name='unmodified protein monomer')
	location_fk = ForeignKey(Location, related_name='location_monomer_mod', verbose_name='location_monomer_mod')
	comment_fk = ForeignKey(Comment, related_name='comment_monomer_mod', verbose_name='comment_monomer_mod', default = '')

class ProteinMonomerModifiedReaction(models.Model):
	protein_monomer_mod_fk = ForeignKey(ProteinMonomerModified, related_name='%(class)s_monomer_mod', verbose_name='protein monomer mod')
	reaction_id = CharField(max_length=255, blank=True, default='', verbose_name='Reaction ID' ,unique=True)
	ec = CharField(max_length=80, blank=True, default='', verbose_name='EC')
	reaction_direction = CharField(max_length=4, blank=True, default='0', verbose_name='reaction direction')

class ProteinMonomerModReactionEnzyme(models.Model):
	reaction_fk = ForeignKey(ProteinMonomerModifiedReaction, related_name='%(class)s_reaction', verbose_name='protein monomer_mod')
	reaction_enzyme_fk = ForeignKey(Molecule, related_name='%(class)s_related_reaction_enz', verbose_name='reaction_enzyme')

class ProteinMonomerModReactionRelation(models.Model):
	reaction_fk = ForeignKey(ProteinMonomerModifiedReaction, related_name='%(class)s_reaction', verbose_name='protein monomer mod reaction')
	reactant_relation = ForeignKey(RelationStoichiometry, related_name='%(class)s_stoichiometry', verbose_name='relation stoichiometry')
	
#################

class RnaModified(models.Model):
	rna_mod = ForeignKey(Molecule, related_name='molecule_rna_mod', verbose_name='rna mod',unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	unmodified_rna_fk = ForeignKey(Rna, related_name='rna_fk', verbose_name='unmodified rna')
	location_fk = ForeignKey(Location, related_name='location_rna_mod', verbose_name='location_fk')
	comment_fk = ForeignKey(Comment, related_name='comment_rna_mod', verbose_name='comment_fk', default = '')

class RnaModifiedReaction(models.Model):
	rna_mod_fk = ForeignKey(RnaModified, related_name='%(class)s_rna_mod', verbose_name='rna mod fk')
	reaction_id = CharField(max_length=255, blank=True, default='', verbose_name='Reaction ID' ,unique=True)
	ec = CharField(max_length=80, blank=True, default='', verbose_name='EC')
	reaction_direction = CharField(max_length=4, blank=True, default='0', verbose_name='reaction direction')

class RnaModReactionEnzyme(models.Model):
	reaction_fk = ForeignKey(RnaModifiedReaction, related_name='%(class)s_reaction', verbose_name='protein complex_mod')
	reaction_enzyme_fk = ForeignKey(Molecule, related_name='%(class)s_related_reaction_enz', verbose_name='reaction_enzyme')

class RnaModifiedReactionRelation(models.Model):
	rna_mod_reaction_fk = ForeignKey(RnaModifiedReaction, related_name='%(class)s_reaction', verbose_name='rna mod reaction')
	reactant_relation = ForeignKey(RelationStoichiometry, related_name='%(class)s_stoichiometry', verbose_name='relation stoichiometry')

####################

class MetaboliteReaction(models.Model):
	frame_id = CharField(max_length=255, blank=True, default='', verbose_name='Frame ID' ,unique=True)
	name = CharField(max_length=255, blank=True, default='', verbose_name='Name')
	ec = CharField(max_length=80, blank=True, default='', verbose_name='EC')
	reaction_direction = CharField(max_length=4, blank=True, default='0', verbose_name='reaction direction')
	upper_bound = FloatField(verbose_name='upper bound')
	lower_bound = FloatField(verbose_name='lower bound')	
	comment_fk = ForeignKey(Comment, related_name='%(class)s_comment_id', verbose_name='comment_fk', default = '')

class MetaboliteReactionEnzyme(models.Model):
	metabolite_reaction_fk = ForeignKey(MetaboliteReaction, related_name='%(class)s_reaction_id', verbose_name='metabolite reaction')
	enzyme_fk = ForeignKey(Molecule, related_name='%(class)s_enzyme', verbose_name='enzyme')

class MetaboliteReactionRelation(models.Model):
	metabolite_reaction_fk = ForeignKey(MetaboliteReaction, related_name='%(class)s_metabolite_reaction', verbose_name='metabolite reaction')
	reactant_relation = ForeignKey(RelationStoichiometry, related_name='%(class)s_stoichiometry', verbose_name='relation stoichiometry')

###################################################

