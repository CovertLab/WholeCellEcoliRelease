from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"adder_sizer.py",
	"cell_growth.py",
	"doubling_time_histogram.py",
	"growthConditionComparison.py",
	"growth_condition_comparison_validation.py",
	"growth_rate_time_series.py",
	"massFractionSummary.py",
	"meneSensitivity.py",
	"metabolism_kinetic_objective_weight.py",
	"metabolism_secretion_penalty.py",
	"param_sensitivity.py",
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
	'PAPER': [
		"adder_sizer.py",
		"doubling_time_histogram.py",
		"meneSensitivity.py",
		"param_sensitivity.py",
		"growth_condition_comparison_validation.py",
		],
	'VALIDATION': [
		'cell_growth.py',
		'growth_rate_time_series',
		],
	# Variant analyses to run with a specific simulation variant (key)
	'ADD_ONE_AA': [
		'cell_growth.py',
		'growth_rate_time_series',
		],
	'ADD_ONE_AA_SHIFT': [
		'cell_growth.py',
		'growth_rate_time_series',
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
	'REMOVE_AA_INHIBITION': [
		'remove_aa_inhibition.py',
		],
	'REMOVE_ONE_AA': [
		'cell_growth.py',
		'growth_rate_time_series',
		],
	'REMOVE_ONE_AA_SHIFT': [
		'cell_growth.py',
		'growth_rate_time_series',
		],
	'TF_ACTIVITY': [
		'tfFit.py',
		'tfFitComparison.py',
		],
	'TIME_STEP': [
		'time_step.py',
		],
	}
