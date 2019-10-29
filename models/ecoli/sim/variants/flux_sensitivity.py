"""
Used for sensitivity analysis of succinate and isocitrate dehydrogenase fluxes
for each metabolic kinetic constraint.

Useful with variant analysis script flux_sensitivity.py.

2 reactions have high simulation flux compared to Toya et al. 2010. fluxes when
CONSTRAINTS_TO_DISABLE are disabled:
	SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.
	ISOCITDEH-RXN

Modifies:
	sim_data.process.metabolism.constraintsToDisable
	sim_data.process.metabolism.run_flux_sensitivity

Expected variant indices:
	0: run flux sensitivity
"""


CONSTRAINTS_TO_DISABLE = [
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.',
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
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


def flux_sensitivity_indices(sim_data):
	return 0

def flux_sensitivity(sim_data, index):
	sim_data.process.metabolism.constraintsToDisable = CONSTRAINTS_TO_DISABLE
	sim_data.process.metabolism.run_flux_sensitivity = True

	return dict(
		shortName='flux_sensitivity',
		desc='Disable kinetics constraints to determine flux impact'
		), sim_data
