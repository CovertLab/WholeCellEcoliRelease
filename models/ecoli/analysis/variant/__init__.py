# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"adder_sizer.py",
	"doubling_time_histogram.py",
	"param_sensitivity.py",
	"meneSensitivity.py",
	"growth_condition_comparison_validation.py",
	"growthConditionComparison.py",
	"massFractionSummary.py",
	"metabolism_kinetic_objective_weight.py",
	"tfFit.py",
	"tfFitComparison.py",
	"tfTarget.py",
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
	}
