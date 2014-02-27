# Create your views here.

from django.http import HttpResponse


import os
import re
os.environ['DJANGO_SETTINGS_MODULE'] = 'ecoliwholecellkb.settings'
import ecoliwholecellkb.settings

from public.models import Gene, Molecule, Location, Comment, ProteinMonomers, Rna, Metabolite, ProteinComplex, ProteinComplexModified, ProteinMonomerModified, RnaModified, RelationStoichiometry, ProteinComplexReactionRelation,ProteinComplexModifiedReaction,ProteinComplexModReactionRelation, ProteinComplexModReactionEnzyme, ProteinMonomerModifiedReaction, ProteinMonomerModReactionEnzyme, ProteinMonomerModReactionRelation, RnaModifiedReaction, RnaModReactionEnzyme, RnaModifiedReactionRelation, MetaboliteReaction, MetaboliteReactionEnzyme, MetaboliteReactionRelation, MetaboliteBiomass, MetaboliteEquivalentEnzyme, Chromosome, GeneSplices, GeneAbsolutentPosition, EntryPositiveFloatData, GeneType, TranscriptionUnit,Promoter, Terminator

from copy import deepcopy
from django.shortcuts import render, render_to_response, get_object_or_404
from django.db.models.query import EmptyQuerySet
from django.db import models
import datetime
from django.core.urlresolvers import reverse
from django.template import RequestContext, loader
from django.db.models.fields.related import OneToOneField, RelatedObject, ManyToManyField, ForeignKey
from django.db.models import Count, Sum, Avg
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField

##
from django.db.models import get_app, get_models, get_model

def index(request):
	species = {}
	species['name'] = 'Escherichia coli'
	#species['wid'] = '4'
	species['comments'] = 'E. coli is a gram-negative, facultative anaerobic, rod-shaped bacterium that is commonly found in the lower intestine of warm-blooded organisms. E. coli is frequently used in modern biological engineering because of its long history of laboratory culture and ease of manupulation. This PGDB provides a comprehensive description of E. coli molecular biology.'
	
	species['genetic_code'] = 11

	## Content
	content = []
	
	
	content.append([
			[0, 'Compartments', Location.objects.count(),None]
	])
	#
	all_seq = Chromosome.objects.all()
	genome = ''
	genome_length = 0	
	gc_content = 0.0

	for i in all_seq:
		genome = i.sequence
		genome_length = i.length
		gc_content = i.gc_content		
		break
	
	content.append([
			[0, 'Chromosomes', len(all_seq),None],
			[1, 'Length', genome_length, 'nt',None],
			[1, 'GC-content', ('%0.1f' % gc_content), '%'],
	])		
	
	content.append([
			[0, 'Transcription units', TranscriptionUnit.objects.count(),None, reverse('public.views.list', kwargs={'model_type': 'Transcriptionunit'})],
			[1, 'Promoters', Promoter.objects.count(),None, reverse('public.views.list', kwargs={'model_type': 'Promoter'})],
			[1, 'Terminators', Terminator.objects.count(),None, reverse('public.views.list', kwargs={'model_type': 'Terminator'})]
		])

	content.append([
			[0, 'Genes', Gene.objects.count(), None, reverse('public.views.list', kwargs={'model_type': 'Gene'}) ],
			[1, 'mRNA', Gene.objects.filter(typegene_id = 1).count(), None, reverse('public.views.list', kwargs={'model_type': 'Gene'}) + '?typegene=1'],
			[1, 'tRNA', Gene.objects.filter(typegene_id = 2).count(), None, reverse('public.views.list', kwargs={'model_type': 'Gene'}) + '?typegene=2'],
			[1, 'rRNA', Gene.objects.filter(typegene_id = 4).count(), None, reverse('public.views.list', kwargs={'model_type': 'Gene'}) + '?typegene=4'],
			[1, 'miscRNA', Gene.objects.filter(typegene_id = 3).count(), None, reverse('public.views.list', kwargs={'model_type': 'Gene'}) + '?typegene=3']
			])
	monomer_count = ProteinMonomers.objects.count()
	complex_count = ProteinComplex.objects.count()
	monomermod_count = ProteinMonomerModified.objects.count()
	complexmod_count = ProteinComplexModified.objects.count()
	rna_count = Rna.objects.count()
	rnamod_count = RnaModified.objects.count()
	

	content.append([
			[0, 'RNA', rna_count + rnamod_count, None],
			[1, 'RNA', rna_count, None, reverse('public.views.list', kwargs={'model_type': 'Rna'})],
			[1, 'RNA Modified', rnamod_count, None, reverse('public.views.list', kwargs={'model_type': 'RnaModified'})]
			])


	content.append([
			[0, 'Proteins', monomer_count + complex_count + monomermod_count + complexmod_count, None],
			[1, 'Protein Monomers', monomer_count, None, reverse('public.views.list', kwargs={'model_type': 'ProteinMonomers'})],
			[1, 'Protein Monomer Modified', monomermod_count, None, reverse('public.views.list', kwargs={'model_type': 'ProteinMonomerModified'})],
			[1, 'Protein Complexes', complex_count, None, reverse('public.views.list', kwargs={'model_type': 'ProteinComplex'})],
			[1, 'Protein Complex Modified', complexmod_count, None, reverse('public.views.list', kwargs={'model_type': 'ProteinComplexModified'})]
			])

	content.append([
			[0, 'Metabolites', Metabolite.objects.count(),None, reverse('public.views.list', kwargs={'model_type': 'Metabolite'})]
	])


	#
	nContent = [len(x) for x in content]
	totContent = sum(nContent)
	cum = 0
	idx = 0
	breakIdxs = [0, 0]
	for x in nContent:
		cum += x
		idx += 1
		if cum > totContent * 1/ 3 and breakIdxs[0] == 0:
			breakIdxs[0] = idx
		if cum > totContent * 2 / 3 and breakIdxs[1] == 0:
			breakIdxs[1] = idx		
			
	contentCol1 = []
	contentCol2 = []
	contentCol3 = []
	i = 0
	for x in content[:breakIdxs[0]]:
		i += 1
		for y in x:
			contentCol1.append([i] + y)	
	i = 0
	for x in content[breakIdxs[0]:breakIdxs[1]]:
		i += 1
		for y in x:
			contentCol2.append([i] + y)
	i = 0
	for x in content[breakIdxs[1]:]:
		i += 1
		for y in x:
			contentCol3.append([i] + y)

			
	##

	context =  {'species': species, 
				'content': [contentCol1, contentCol2, contentCol3],
				'contentRows': range(max(len(contentCol1), len(contentCol2), len(contentCol3))),				
				}

	return render(request, 'public/index.html',context)


def filter_fields_list(model_type):

	filters = {'Gene': {'typegene' :{'formal':'Type','name':'type_gene'}, 
						'direction':{'formal':'Direction','name':'name'}
					   },
				'Transcriptionunit': {'direction':{'formal':'Direction','name':'name'}},
				'Rna': {'is_modified': {'formal':'Has modified form?', 'name':'name'}},
				'ProteinComplex': {'modified_form': {'formal':'Has modified form?', 'name':'name'}},
				'ProteinMonomers': {'is_modified': {'formal':'Has modified form?', 'name':'name'}},					   
			  }
	
	if model_type in filters:
		return filters[model_type]
	return {}

def list(request, model_type):

	fields_to_show = filter_fields_list(model_type)

	model = get_model( 'public',model_type)
	objects = model.objects.all()
	
	facet_fields = []	
	for field in model._meta.fields:
		#facet
		field_name = field.name
		
		if field_name not in fields_to_show: 
			continue

		tmp_model = model

		if field.get_internal_type() == "ForeignKey":
			tmp_model = field.rel.to

		tmp = model.objects.order_by(field_name).values(field_name).annotate(count=Count(field_name))
		
		facets = []

		for facet in tmp:
			value = facet[field_name]	
			
			if value is None or unicode(value) == '':
				continue
 
			if field.get_internal_type() == "ForeignKey":
				
				tmp2 = tmp_model.objects.values('id', fields_to_show[field_name]['name']).get(id=value)
				id = tmp2['id']
				name = tmp2[fields_to_show[field_name]['name']]

			elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):	
				
				id = value
				choices = [x[0] for x in field.choices]
				if id in choices:
					name = field.choices[choices.index(id)][1]
				else:
					name = value
			else:
				id = value
				name = value
			
			if name == '0':
				name = 'No'
			elif name == '1':
				name = 'Yes'
		
			if value is not None and unicode(value) != '':
				facets.append({
					'id': unicode(id), 
					'name': unicode(name),
					'count': facet['count']
				})
		
		if len(facets) > 1:
			facet_fields.append({ 
				'name': field_name,
				'verbose_name': fields_to_show[field_name]['formal'],
				'facets': facets,
				})
	
		#filter
		val = request.GET.get(field_name)		
		if val:
			if isinstance(field, (ForeignKey, ManyToManyField)):
				kwargs = {field_name +'_id': val}
			elif isinstance(field, (BooleanField, NullBooleanField)):
				kwargs = {field_name: val == 'True'}
			elif isinstance(field, (AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField)):
				kwargs = {field_name: float(val)}
			else:
				kwargs = {field_name: val}
			objects = objects.filter(**kwargs)

	
	##########

	data = {}
	for key, val in request.GET.iterlists():
		data[key] = val
	

	context = {
			'model_type': model_type,
			'facet_fields': facet_fields,
			'queryset' : objects,
			'queryargs' : data
			}

	return render(request, 'public/list.html',context)

def fieldset_lists(model_type):

	fieldsets = {}
	fieldsets['Gene'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'frame_id'}, 
					{'verbose_name':'Name','name':'name'}, 
					{'verbose_name':'Symbol','name':'symbol'}, 
					{'verbose_name': 'Product Name','name':'productname','values':['product']}
					]),
			('Classification',  
					[{'verbose_name':'Type','name':'typegene', 'values' :['type_gene']}
					]), 
			('Structure', 
				[{'verbose_name': 'Structure', 'name': 'structure'}, 
				{'verbose_name': 'Sequence', 'name': 'sequence'}, 
				{'verbose_name': 'Transcription unit', 'name': 'transcription_units'},
				{'verbose_name': 'Empirical formula (pH 7.5)', 'name': 'empirical_formula'},
				{'verbose_name': 'Molecular weight (pH 7.5; Da)', 'name': 'molecular_weight'}
				]), 
			('Functional genomics', 
				[{'verbose_name': 'Expression', 'name': 'expression', 'values' :['value']},
				{'verbose_name': 'Half Life', 'name': 'half_life', 'values' :['value', 'units']},
				]), 
			('Metadata', [{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]
	
	fieldsets['Metabolite'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'metabolite_id', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}, 
					]),
			('Structure', 
				[
				{'verbose_name': 'Feist Formula', 'name': 'feist_formula'},
				{'verbose_name': 'pH Formula', 'name': 'ph_formula'},
				{'verbose_name': 'pH Charge', 'name': 'ph_charge'},
				{'verbose_name': 'pH Weight', 'name': 'ph_weight'},
				{'verbose_name': 'Media Concentration', 'name': 'media_concentration'},
				{'verbose_name': 'MAximum Exchange Rate', 'name': 'maximum_exchange_rate'}
				]), 
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['TranscriptionUnit'] = [
			('Name',  
					[{'verbose_name':'Frame_ID','name':'transcription_unit_id'}, 
					{'verbose_name':'Name','name':'name'}
					]), 
			('Structure', 
				[{'verbose_name': 'Structure', 'name': 'structure'}, 
				{'verbose_name': 'Sequence', 'name': 'sequence'}, 
				{'verbose_name': 'Direction', 'name': 'direction'},
				{'verbose_name': 'Degradation Rate', 'name': 'degradation_rate'},
				{'verbose_name': 'Expression Rate', 'name': 'expression_rate'},
				{'verbose_name': 'Promoter(s)', 'name': 'promoter_id_fk', 'values' :['name']}
				])			
	]

	fieldsets['Promoter'] = [
			('Name',  
					[{'verbose_name':'Frame_ID','name':'promoter_id'}, 
					{'verbose_name':'Name','name':'name'}
					]), 
			('Structure', 
				[{'verbose_name': 'Direction', 'name': 'direction'},
				{'verbose_name': 'Position', 'name': 'position'}
				])			
	]


	fieldsets['Terminator'] = [
			('Name',  
					[{'verbose_name':'Frame_ID','name':'terminator_id'}, 
					{'verbose_name':'Name','name':'name'}
					]), 
			('Structure', 
				[{'verbose_name': 'Left', 'name': 'left'},
				{'verbose_name': 'Right', 'name': 'right'},
				{'verbose_name': 'Rho dependent', 'name': 'rho_dependent'}
				])			
	]

	fieldsets['Rna'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'frame_id', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[{'verbose_name': 'Gene','name':'gene_fk','values':['frame_id']},
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Modified Form', 'name': 'is_modified','values':['rna_mod_id']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['RnaModified'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'rna_mod', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Unmodified Form', 'name': 'unmodified_rna_fk','values':['name']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['ProteinMonomers'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'frame_id', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[{'verbose_name': 'Gene','name':'gene_fk','values':['frame_id']},
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Modified Form', 'name': 'is_modified','values':['rna_mod_id']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['ProteinMonomerModified'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'protein_monomer_mod', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Unmodified Form', 'name': 'unmodified_protein_monomer_fk','values':['name']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['ProteinComplex'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'protein_complex', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Modified Form', 'name': 'modified_from','values':['rna_mod_id']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]

	fieldsets['ProteinComplexModified'] = [
			('Name',  
					[{'verbose_name':'Frame ID','name':'protein_complex_mod', 'values':['product']}, 
					{'verbose_name':'Name','name':'name'}
					]),
			('Structure', 
				[
				{'verbose_name': 'Location', 'name': 'location_fk','values':['location_id']},
				{'verbose_name': 'Unmodified Form', 'name': 'unmodified_protein_complex_fk','values':['name']},
				 ]),
			('Metadata', [{'verbose_name': 'Comment', 'name': 'comment_fk','values':['comment_str']},{'verbose_name': 'Created', 'name': 'created_user'}, {'verbose_name': 'Last updated', 'name': 'last_updated_user'}
				])
	]


	if model_type in fieldsets:
		return fieldsets[model_type]
	return {}

def detail(request, model_type, frame_id):
		
	model = get_model( 'public',model_type)
	
	if model_type == 'TranscriptionUnit':
		obj = model.objects.filter(transcription_unit_id=frame_id)
	elif model_type == 'Promoter':
		obj = model.objects.filter(promoter_id=frame_id)
	elif model_type == 'Terminator':
		obj = model.objects.filter(terminator_id=frame_id)	
	elif model_type == 'Rna':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(frame_id=x)	
	elif model_type == 'RnaModified':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(rna_mod=x)	
	elif model_type == 'ProteinMonomers':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(frame_id=x)	
	elif model_type == 'ProteinMonomerModified':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(protein_monomer_mod=x)	
	elif model_type == 'ProteinComplex':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(protein_complex=x)	
	elif model_type == 'ProteinComplexModified':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(protein_complex_mod_id=x)	
	elif model_type == 'Metabolite':
		x = Molecule.objects.filter(product=frame_id).values('id')
		x = x[0]['id']
		obj = model.objects.filter(metabolite_id=x)	

	else:
		obj = model.objects.filter(frame_id=frame_id)

	fields = model._meta.fields
	field_names = [x.name for x in fields]

	fieldsets = fieldset_lists(model_type)		

	#form query set
	qs = objectToQuerySet(obj, model = model)

	q_dict = obj.values()
	q_dict = q_dict[0]
 	f_key = {}
	
	for field in fields:
		if field.name == 'chromosome':
			continue
		if field.get_internal_type() == "ForeignKey":
			tmp_model = field.rel.to
			tmp2 = tmp_model.objects.values().get(id=q_dict[field.name+'_id'])
			f_key[field.name]= tmp2


	for i in fieldsets:
		for field in i[1]:		
			if field['name'] in q_dict:
				field['data'] = q_dict[field['name']] 
			elif field['name'] in f_key:
				for j in field['values']:
					if 'data' in field:
						field['data'] = field['data'] + ' ' + str(f_key[field['name']][j]) 
					else:
						field['data'] = str(f_key[field['name']][j] )
			
			
	one_object = ''
	for i in obj:
		one_object = i
		break

	data = {'queryset' : qs,
			'queryargs' : {},
			'model_type': model_type,
			'object': one_object,			
			'fieldsets': fieldsets,
			'Name': q_dict['name'],
			'Frame_ID': one_object.get_id(),
			'message': request.GET.get('message', ''),
			}

	for key, val in request.GET.iterlists():
		data['queryargs'][key] = val


	context = data

	return render(request, 'public/detail.html',context)
	

	##############
	'''
	obj = getEntry(species_wid = species_wid, wid = wid)
	if obj is None:
		raise Http404
	
	model = obj.__class__
	model_type = model.__name__
	fieldsets = deepcopy(model._meta.fieldsets)
	
	#filter out type, metadata
	fieldset_names = [x[0] for x in fieldsets]
	if 'Type' in fieldset_names:
		idx = fieldset_names.index('Type')
		del fieldsets[idx]
		
	#filter out empty fields
	rmfieldsets = []
	for idx in range(len(fieldsets)):
		rmfields = []
		for idx2 in range(len(fieldsets[idx][1]['fields'])):
			if isinstance(fieldsets[idx][1]['fields'][idx2], dict):
				field = fieldsets[idx][1]['fields'][idx2]
				field_name = field['name']
				verbose_name = field['verbose_name']
			else:
				field_name = fieldsets[idx][1]['fields'][idx2]
				field = model._meta.get_field_by_name(field_name)[0]
				if isinstance(field, RelatedObject):
					verbose_name = capfirst(field.get_accessor_name())
				else:
					verbose_name = field.verbose_name
				
			data = format_field_detail_view(obj, field_name, request.user.is_anonymous())
			if (data is None) or (data == ''):
				rmfields = [idx2] + rmfields
				
			if issubclass(model, Parameter) and field_name == 'value':
				is_modeled_url = reverse('public.views.viewParameterInSimulation', kwargs = {
					'species_wid':species_wid, 
					'wid': wid,
					})
			elif (isinstance(field, (ForeignKey, ManyToManyField)) and issubclass(field.rel.to, Parameter)) or 			   (isinstance(field, (RelatedObject)) and issubclass(field.model, Parameter)):
				is_modeled_url = reverse('public.views.viewParametersInSimulation', kwargs = {
					'species_wid':species_wid, 
					'wid': wid,
					})
			else:
				try:
					if ModelProperty.objects.get(species__wid = species_wid, class_name = model_type, property_name = field_name).simulation_properties.exists():
						is_modeled_url = reverse('public.views.viewPropertyInSimulation', kwargs = {
							'species_wid':species_wid, 
							'class_name': model_type, 
							'property_name': field_name,
							})
				except ObjectDoesNotExist:
					is_modeled_url = ''
			
			fieldsets[idx][1]['fields'][idx2] = {
				'name': field_name,
				'verbose_name': verbose_name.replace(" ", '&nbsp;').replace("-", "&#8209;"), 
				'data': data,
				'is_modeled_url': is_modeled_url,
				}
		for idx2 in rmfields:
			del fieldsets[idx][1]['fields'][idx2]
		if len(fieldsets[idx][1]['fields']) == 0:
			rmfieldsets = [idx] + rmfieldsets
	for idx in rmfieldsets:
		del fieldsets[idx]
	'''



######### healper
def objectToQuerySet(obj, model = None):
	if model is None:
		model = obj.__class__
	qs = model.objects.none()
	qs._result_cache.append(obj)
	return qs

