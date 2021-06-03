from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"adder_sizer.py",
	"cell_growth.py",
	"doubling_time_histogram.py",
	"growthConditionComparison.py",
	"growth_condition_comparison_validation.py",
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
		],
	# Variant analyses to run with a specific simulation variant (key)
	'add_one_aa': [
		'cell_growth.py',
		],
	'remove_aa_inhibition': [
		'remove_aa_inhibition.py',
		],
	'remove_one_aa': [
		'cell_growth.py',
		],
	}
