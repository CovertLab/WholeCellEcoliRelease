import importlib


variants = [
	'aa_uptake_sensitivity',
	'add_one_aa',
	'add_one_aa_shift',
	'all_shuffle_params',
	'condition',
	'gene_knockout',
	'mene_params',
	'metabolism_kinetic_objective_weight',
	'metabolism_secretion_penalty',
	'monomer_deg_rate_shuffle_params',
	'param_sensitivity',
	'remove_aa_inhibition',
	'remove_one_aa',
	'remove_one_aa_shift',
	'rna_deg_rate_shuffle_params',
	'rrna_orientation',
	'tf_activity',
	'time_step',
	'timelines',
	'transcription_initiation_shuffle_params',
	'translation_efficiencies_shuffle_params',
	'wildtype',
	]

def get_function(variant):
	module = importlib.import_module(f'models.ecoli.sim.variants.{variant}')
	return getattr(module, variant)

nameToFunctionMapping = {v: get_function(v) for v in variants}

# Support the old names for compatibility with existing shell scripts.
nameToFunctionMapping.update({
	'allShuffleParams': get_function('all_shuffle_params'),
	'geneKnockout': get_function('gene_knockout'),
	'meneParams': get_function('mene_params'),
	'monomerDegRateShuffleParams': get_function('monomer_deg_rate_shuffle_params'),
	'nutrientTimeSeries': get_function('timelines'),
	'tfActivity': get_function('tf_activity'),
	'transcriptionInitiationShuffleParams': get_function('transcription_initiation_shuffle_params'),
	'translationEfficienciesShuffleParams': get_function('translation_efficiencies_shuffle_params'),
	})
