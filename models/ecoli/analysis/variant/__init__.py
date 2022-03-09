from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"aa_synthesis_sensitivity.py",
	"aa_uptake_sensitivity.py",
	"adder_sizer.py",
	"cell_growth.py",
	"doubling_time_histogram.py",
	"growthConditionComparison.py",
	"growth_condition_comparison_validation.py",
	"growth_correlations.py",
	"growth_expression_comparison.py",
	"growth_rate_time_series.py",
	"growth_trajectory.py",
	"massFractionSummary.py",
	"meneSensitivity.py",
	"metabolism_kinetic_objective_weight.py",
	"metabolism_secretion_penalty.py",
	"param_sensitivity.py",
	"ppgpp_conc.py",
	"remove_aa_inhibition.py",
	"tfFit.py",
	"tfFitComparison.py",
	"time_step.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"growthConditionComparison.py",
		],
	'ENVIRONMENT_SHIFT': [
		'growth_trajectory',
		],
	'PAPER': [
		"adder_sizer.py",
		"doubling_time_histogram.py",
		"meneSensitivity.py",
		"param_sensitivity.py",
		"growth_condition_comparison_validation.py",
		],
	'VALIDATION': [
		'cell_growth.py',
		'growth_expression_comparison',
		'growth_rate_time_series',
		],
	# Variant analyses to run with a specific simulation variant (key)
	'AA_SYNTHESIS_SENSITIVITY': [
		'aa_synthesis_sensitivity',
		],
	'AA_UPTAKE_SENSITIVITY': [
		'aa_uptake_sensitivity',
		],
	'ADD_ONE_AA': [
		'cell_growth.py',
		'growth_correlations',
		'growth_rate_time_series',
		],
	'ADD_ONE_AA_SHIFT': [
		'cell_growth.py',
		'growth_correlations',
		'growth_rate_time_series',
		'growth_trajectory',
		],
	'METABOLISM_KINETIC_OBJECTIVE_WEIGHT': [
		'metabolism_kinetic_objective_weight.py',
		],
	'METABOLISM_SECRETION_PENALTY': [
		'metabolism_secretion_penalty.py',
		],
	'PARAM_SENSITIVITY': [
		'param_sensitivity.py',
		],
	'PPGPP_CONC': [
		'ppgpp_conc',
		],
	'PPGPP_LIMITATIONS': [
		'ppgpp_conc',
		],
	'PPGPP_LIMITATIONS_RIBOSOME': [
		'ppgpp_conc',
		],
	'REMOVE_AA_INHIBITION': [
		'remove_aa_inhibition.py',
		],
	'REMOVE_AAS_SHIFT': [
		'cell_growth.py',
		'growth_correlations',
		'growth_rate_time_series',
		'growth_trajectory',
		],
	'REMOVE_ONE_AA': [
		'cell_growth.py',
		'growth_correlations',
		'growth_rate_time_series',
		],
	'REMOVE_ONE_AA_SHIFT': [
		'cell_growth.py',
		'growth_correlations',
		'growth_rate_time_series',
		'growth_trajectory',
		],
	'TF_ACTIVITY': [
		'tfFit.py',
		'tfFitComparison.py',
		],
	'TIME_STEP': [
		'time_step.py',
		],
	'TIMELINES': [
		'growth_trajectory',
		],
	}
