"""
Removes amino acid inhibition feedback for single amino acids to match data in
Sander et al. Allosteric Feedback Inhibition Enables Robust Amino Acid
Biosynthesis in E. coli by Enforcing Enzyme Overabundance. 2019.

Associated variant analysis:
	remove_aa_inhibition: to produce plot similar to Fig 1B from Sander et al.

Modifies:
	sim_data.process.metabolism.aa_kis

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: wildtype
	1-7: enzyme mutants that show no feedback control (determined by order of AA_TO_ENZYME
"""

import numpy as np


AA_TO_ENZYME = {
	'ARG[c]': 'argA',
	'TRP[c]': 'trpE',
	'HIS[c]': 'hisG',
	'LEU[c]': 'leuA',
	'THR[c]': 'thrA',
	'ILE[c]': 'ilvA',
	'PRO[c]': 'proB',
	}


def remove_aa_inhibition(sim_data, index):
	if index > 0:
		aa = list(AA_TO_ENZYME.keys())[index-1]
		aa_index = sim_data.molecule_groups.amino_acids.index(aa)
		sim_data.process.metabolism.aa_kis[aa_index] *= np.inf
		short = aa
		desc = f'remove {aa} inhibition on {AA_TO_ENZYME[aa]} activity'
	else:
		short = 'wt'
		desc = 'wildtype with no enzyme adjustment'

	return dict(shortName=short, desc=desc), sim_data
