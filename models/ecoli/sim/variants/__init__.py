from __future__ import absolute_import, division, print_function

from .gene_knockout import gene_knockout
from .wildtype import wildtype
from .timelines import timelines
from .tf_activity import tf_activity
from .condition import condition
from .transcription_initiation_shuffle_params import transcription_initiation_shuffle_params
from .kinetic_target_shuffle_params import kinetic_target_shuffle_params
from .catalyst_shuffle_params import catalyst_shuffle_params
from .translation_efficiencies_shuffle_params import translation_efficiencies_shuffle_params
from .monomer_deg_rate_shuffle_params import monomer_deg_rate_shuffle_params
from .kinetic_catalyst_shuffle_params import kinetic_catalyst_shuffle_params
from .all_shuffle_params import all_shuffle_params
from .mene_params import mene_params
from .metabolism_kinetic_objective_weight import metabolism_kinetic_objective_weight
from .rna_deg_rate_shuffle_params import rna_deg_rate_shuffle_params


nameToFunctionMapping = {
	"all_shuffle_params": all_shuffle_params,
	"catalyst_shuffle_params": catalyst_shuffle_params,
	"condition": condition,
	"gene_knockout": gene_knockout,
	"kinetic_catalyst_shuffle_params": kinetic_catalyst_shuffle_params,
	"kinetic_target_shuffle_params": kinetic_target_shuffle_params,
	"mene_params": mene_params,
	"metabolism_kinetic_objective_weight": metabolism_kinetic_objective_weight,
	"monomer_deg_rate_shuffle_params": monomer_deg_rate_shuffle_params,
	"timelines": timelines,
	"rna_deg_rate_shuffle_params": rna_deg_rate_shuffle_params,
	"tf_activity": tf_activity,
	"transcription_initiation_shuffle_params": transcription_initiation_shuffle_params,
	"translation_efficiencies_shuffle_params": translation_efficiencies_shuffle_params,
	"wildtype": wildtype,

	# Support the old names for compatibility with existing shell scripts.
	"allShuffleParams": all_shuffle_params,
	"catalystShuffleParams": catalyst_shuffle_params,
	"geneKnockout": gene_knockout,
	"kineticCatalystShuffleParams": kinetic_catalyst_shuffle_params,
	"kineticTargetShuffleParams": kinetic_target_shuffle_params,
	"meneParams": mene_params,
	"monomerDegRateShuffleParams": monomer_deg_rate_shuffle_params,
	"nutrientTimeSeries": timelines,
	"tfActivity": tf_activity,
	"transcriptionInitiationShuffleParams": transcription_initiation_shuffle_params,
	"translationEfficienciesShuffleParams": translation_efficiencies_shuffle_params,
	}
