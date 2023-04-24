"""
tRNA Synthetase Kinetics variant for simulations using the curated k_cat
measurements.

Modifies:
	sim_data.relation.synthetase_to_k_cat
	sim_data.trna_synthetase_kinetics_variant

Expected variant indices:
	0: control (KineticTrnaChargingModel)
	1: previous (TranslationSupplyElongationModel)
	2: max measured (CoarseKineticTrnaChargingModel)
	3: HisRS kcat drop (KineticTrnaChargingModel)
	4: HisRS kcat, intermediate low (KineticTrnaChargingModel)
	5: HisRS kcat, intermediate high (KineticTrnaChargingModel)
	6: ArgRS kcat drop (KineticTrnaChargingModel)

All variants occur in the basal condition (M9 + glucose, aerobic)
"""

import numpy as np
from wholecell.utils import units

synthetases_to_abbreviation = {
	'ALAS-CPLX[c]': 'AlaRS',
	'ARGS-MONOMER[c]': 'ArgRS',
	'ASNS-CPLX[c]': 'AsnRS',
	'ASPS-CPLX[c]': 'AspRS',
	'CYSS-MONOMER[c]': 'CysRS',
	'GLURS-MONOMER[c]': 'GluRS',
	'GLNS-MONOMER[c]': 'GlnRS',
	'GLYS-CPLX[c]': 'GlyRS',
	'HISS-CPLX[c]': 'HisRS',
	'ILES-MONOMER[c]': 'IleRS',
	'LEUS-MONOMER[c]': 'LeuRS',
	'LYSS-CPLX[c]': 'LysRS',
	'METG-CPLX[c]': 'MetRS',
	'PHES-CPLX[c]': 'PheRS',
	'PROS-CPLX[c]': 'ProRS',
	'SERS-CPLX[c]': 'SerRS',
	'THRS-CPLX[c]': 'ThrRS',
	'TRPS-CPLX[c]': 'TrpRS',
	'TYRS-CPLX[c]': 'TyrRS',
	'VALS-MONOMER[c]': 'ValRS',
	}

def trna_synthetase_kinetics(sim_data, index):
	# Polypeptide Elongation Model is set in wholecell/sim/simulation.py
	# using the following attribute:
	sim_data.trna_synthetase_kinetics_variant = index

	# Set kcat drops and short names
	if index == 0:
		short_name = 'control (KineticTrnaChargingModel)'

	elif index == 1:
		short_name = 'previous (TranslationSupplyElongationModel)'

	elif index == 2:
		short_name = 'max measured (CoarseKineticTrnaChargingModel)'

	####################################################################
	# HisRS kcat drop
	elif index in [3, 4, 5]:
		synthetase = 'HISS-CPLX[c]'
		optimized_k_cat = sim_data.relation.\
			synthetase_to_k_cat[synthetase].asNumber(1/units.s)
		measured_k_cat = sim_data.relation.\
			synthetase_to_max_curated_k_cats[synthetase].asNumber(1/units.s)
		# k_cats = 1 / units.s * np.linspace(measured_k_cat, optimized_k_cat, 4)
		k_cats = 1 / units.s * np.linspace(measured_k_cat, optimized_k_cat, 7)

		# Measured (142)
		if index == 3:
			k_cat = k_cats[0]

		# Intermediate low 183
		if index == 4:
			k_cat = k_cats[1]

		# Intermediate high 223
		elif index == 5:
			k_cat = k_cats[2]

		sim_data.relation.synthetase_to_k_cat[synthetase] = k_cat
		fold_change = optimized_k_cat / k_cat.asNumber(1/units.s)
		short_name = f'{synthetases_to_abbreviation[synthetase]} kcat drop ({fold_change:.1f}x drop)'

	####################################################################
	# ArgRS kcat drop
	elif index == 6:
		synthetase = 'ARGS-MONOMER[c]'
		optimized_k_cat = sim_data.relation.synthetase_to_k_cat[synthetase]
		k_cat = sim_data.relation.synthetase_to_max_curated_k_cats[synthetase]
		sim_data.relation.synthetase_to_k_cat[synthetase] = k_cat
		fold_change = optimized_k_cat.asNumber(1/units.s) / k_cat.asNumber(1/units.s)
		short_name = f'{synthetases_to_abbreviation[synthetase]} kcat drop ({fold_change:.1f}x drop)'

	####################################################################
	else:
		# Unrecognized indexes run the control variant
		short_name = 'control (KineticTrnaChargingModel)'
		sim_data.trna_synthetase_kinetics_variant = 0

	return dict(
		shortName=short_name,
		desc=f'tRNA Synthetase Kinetics Variant is set at index {index}.'
		), sim_data
