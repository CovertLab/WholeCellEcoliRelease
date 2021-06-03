from .add_one_aa import add_one_aa
from .all_shuffle_params import all_shuffle_params
from .condition import condition
from .gene_knockout import gene_knockout
from .mene_params import mene_params
from .metabolism_kinetic_objective_weight import metabolism_kinetic_objective_weight
from .metabolism_secretion_penalty import metabolism_secretion_penalty
from .monomer_deg_rate_shuffle_params import monomer_deg_rate_shuffle_params
from .param_sensitivity import param_sensitivity
from .remove_aa_inhibition import remove_aa_inhibition
from .remove_one_aa import remove_one_aa
from .rna_deg_rate_shuffle_params import rna_deg_rate_shuffle_params
from .rrna_orientation import rrna_orientation
from .tf_activity import tf_activity
from .time_step import time_step
from .timelines import timelines
from .transcription_initiation_shuffle_params import transcription_initiation_shuffle_params
from .translation_efficiencies_shuffle_params import translation_efficiencies_shuffle_params
from .wildtype import wildtype


nameToFunctionMapping = {
	"add_one_aa": add_one_aa,
	"all_shuffle_params": all_shuffle_params,
	"condition": condition,
	"gene_knockout": gene_knockout,
	"mene_params": mene_params,
	"metabolism_kinetic_objective_weight": metabolism_kinetic_objective_weight,
	"metabolism_secretion_penalty": metabolism_secretion_penalty,
	"monomer_deg_rate_shuffle_params": monomer_deg_rate_shuffle_params,
	"param_sensitivity": param_sensitivity,
	"remove_aa_inhibition": remove_aa_inhibition,
	"remove_one_aa": remove_one_aa,
	"rna_deg_rate_shuffle_params": rna_deg_rate_shuffle_params,
	"rrna_orientation": rrna_orientation,
	"tf_activity": tf_activity,
	"time_step": time_step,
	"timelines": timelines,
	"transcription_initiation_shuffle_params": transcription_initiation_shuffle_params,
	"translation_efficiencies_shuffle_params": translation_efficiencies_shuffle_params,
	"wildtype": wildtype,

	# Support the old names for compatibility with existing shell scripts.
	"allShuffleParams": all_shuffle_params,
	"geneKnockout": gene_knockout,
	"meneParams": mene_params,
	"monomerDegRateShuffleParams": monomer_deg_rate_shuffle_params,
	"nutrientTimeSeries": timelines,
	"tfActivity": tf_activity,
	"transcriptionInitiationShuffleParams": transcription_initiation_shuffle_params,
	"translationEfficienciesShuffleParams": translation_efficiencies_shuffle_params,
	}
