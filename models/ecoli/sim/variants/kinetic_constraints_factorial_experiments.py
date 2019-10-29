"""
Used for analysis kinetic constraint interaction effects on succinate and
isocitrate dehydrogenase fluxes as well as glucose uptake rates.

Useful with variant analysis scripts:
	kinetic_objective_interactions.py
	kinetic_objective_range_violin.py

2 reactions have high simulation flux compared to Toya et al. 2010. fluxes when
CONSTRAINTS_TO_DISABLE (including commented out constraints) are disabled:
	SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.
	ISOCITDEH-RXN

Modifies:
	sim_data.process.metabolism.constraintsToDisable

Expected variant indices (depends on length of FACTORIAL_DESIGN_CONSTRAINTS):
	0: control simulation with the old list of constraints.
	1-256: different combinations of constraints are disabled,
		2**len(FACTORIAL_DESIGN_CONSTRAINTS) combinations
"""


# Commented out constraints are in the factorial design so removed here
CONSTRAINTS_TO_DISABLE = [
	# 'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.',
	# 'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'ALANINE--TRNA-LIGASE-RXN-ALA-tRNAs/L-ALPHA-ALANINE/ATP/PROTON//Charged-ALA-tRNAs/AMP/PPI.64.',
	'ARGININE--TRNA-LIGASE-RXN-ARG-tRNAs/ARG/ATP/PROTON//Charged-ARG-tRNAs/AMP/PPI.52.',
	'ASPARAGINE--TRNA-LIGASE-RXN-ASN-tRNAs/ASN/ATP/PROTON//Charged-ASN-tRNAs/AMP/PPI.52.',
	'ASPARTATE--TRNA-LIGASE-RXN-ASP-tRNAs/L-ASPARTATE/ATP/PROTON//Charged-ASP-tRNAs/AMP/PPI.60.',
	'CYSTEINE--TRNA-LIGASE-RXN-CYS-tRNAs/CYS/ATP/PROTON//Charged-CYS-tRNAs/AMP/PPI.52.',
	'GLUTAMINE--TRNA-LIGASE-RXN-GLN-tRNAs/GLN/ATP/PROTON//Charged-GLN-tRNAs/AMP/PPI.52.',
	'HISTIDINE--TRNA-LIGASE-RXN-HIS-tRNAs/HIS/ATP/PROTON//Charged-HIS-tRNAs/AMP/PPI.52.',
	'ISOLEUCINE--TRNA-LIGASE-RXN-ILE-tRNAs/ILE/ATP/PROTON//Charged-ILE-tRNAs/AMP/PPI.52.',
	'LEUCINE--TRNA-LIGASE-RXN-LEU-tRNAs/LEU/ATP/PROTON//Charged-LEU-tRNAs/AMP/PPI.52.',
	'LYSINE--TRNA-LIGASE-RXN-LYS/LYS-tRNAs/ATP/PROTON//Charged-LYS-tRNAs/AMP/PPI.52.__LYSS-CPLX',
	'METHIONINE--TRNA-LIGASE-RXN-Elongation-tRNAMet/MET/ATP/PROTON//Charged-MET-tRNAs/AMP/PPI.61.',
	'PHENYLALANINE--TRNA-LIGASE-RXN-PHE-tRNAs/PHE/ATP/PROTON//Charged-PHE-tRNAs/AMP/PPI.52.',
	'PROLINE--TRNA-LIGASE-RXN-PRO-tRNAs/PRO/ATP/PROTON//Charged-PRO-tRNAs/AMP/PPI.52.',
	'SERINE--TRNA-LIGASE-RXN-SER-tRNAs/SER/ATP/PROTON//Charged-SER-tRNAs/AMP/PPI.52.',
	'THREONINE--TRNA-LIGASE-RXN-THR-tRNAs/THR/ATP/PROTON//Charged-THR-tRNAs/AMP/PPI.52.',
	'TRYPTOPHAN--TRNA-LIGASE-RXN-TRP/TRP-tRNAs/ATP/PROTON//Charged-TRP-tRNAs/AMP/PPI.52.',
	'TYROSINE--TRNA-LIGASE-RXN-TYR-tRNAs/TYR/ATP/PROTON//Charged-TYR-tRNAs/AMP/PPI.52.',
	'VALINE--TRNA-LIGASE-RXN-VAL-tRNAs/VAL/ATP/PROTON//Charged-VAL-tRNAs/AMP/PPI.52.',
	]

# Reactions that were previously disabled and are now enabled
NEWLY_ENABLED = [
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.',
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	]

# Reactions identified with flux_sensitivity but disabled due to complex regulation
NEWLY_DISABLED = [
	'GLYOXYLATE-REDUCTASE-NADP+-RXN__CPLX0-235',
	'ISOCITDEH-RXN',
	]

# Previously disabled constraints and constraints identified with flux_sensitivity variant
FACTORIAL_DESIGN_CONSTRAINTS = [
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.',
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'CYTDEAM-RXN',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
	'PSERTRANSAM-RXN',
	'CITSYN-RXN__CITRATE-SI-SYNTHASE',
	]


def get_disabled_constraints(index):
	"""
	Determines constraints in the factorial design list to disable for a given index.

	Args:
		index (int): variant index to use to generate the disabled reactions

	Returns:
		list[int]: mask corresponding to each constraint in the factorial design
			(1 if constraint is disabled, 0 if not)
		list[str]: reaction IDs for constraints in the factorial design that
			will be disabled
	"""
	# index == 0 is for running the old constraints.
	if index == 0:
		return None, None
	new_index = index - 1

	disable_constraints = [new_index // 2**i % 2 for i in range(len(FACTORIAL_DESIGN_CONSTRAINTS))]
	additional_disabled = [rxn for rxn, disable in zip(FACTORIAL_DESIGN_CONSTRAINTS, disable_constraints) if disable]

	return disable_constraints, additional_disabled

def kinetic_constraints_factorial_experiments_indices(sim_data):
	return 0

def kinetic_constraints_factorial_experiments(sim_data, index):
	if index == 0:
		sim_data.process.metabolism.constraintsToDisable = CONSTRAINTS_TO_DISABLE + NEWLY_ENABLED
		return dict(
			shortName="control",
			desc="Simulation with old constraints list."
			), sim_data
	else:
		disable_constraints, additional_disabled = get_disabled_constraints(index)
		sim_data.process.metabolism.constraintsToDisable = CONSTRAINTS_TO_DISABLE + NEWLY_DISABLED + additional_disabled

	return dict(
		shortName="reactions disabled: {}".format(disable_constraints),
		desc="Simulation with disabled kinetic targets for {}.".format(additional_disabled)
		), sim_data
