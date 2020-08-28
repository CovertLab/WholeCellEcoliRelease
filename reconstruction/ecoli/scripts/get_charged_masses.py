'''
Print id, mw, location for charged tRNAs to add to modifiedForms.tsv
Output file is charged_data.tsv in the same directory as the script

Notes:
- Masses of small molecules are added to the tRNA mass so that the charging
reactions are mass balanced
'''

from __future__ import absolute_import, division, print_function

import os
from typing import Any

import numpy as np
from six.moves import zip

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from reconstruction.spreadsheets import tsv_writer

# file paths
file_loc = os.path.dirname(__file__)
output_filename = os.path.join(file_loc, 'charged_data.tsv')

# suppress scientific notation output
np.set_printoptions(suppress = True)

# get raw and sim data
raw_data = KnowledgeBaseEcoli()  # type: Any
sim_data = SimulationDataEcoli()
sim_data.initialize(raw_data)

# determine masses and write to output file
with tsv_writer(output_filename, ["id", "mw", "location"]) as writer:
	trnas = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']]
	charged = [x['modifiedForms'] for x in raw_data.rnas if x['id']+'[c]' in trnas]
	filtered_charged = []
	for c1 in charged:
		for c2 in c1:
			if 'FMET' in c2 or 'modified' in c2:
				continue
			filtered_charged += [c2 + '[c]']

	mol_names = sim_data.internal_state.bulkMolecules.bulkData['id']
	mws = sim_data.internal_state.bulkMolecules.bulkData['mass']
	for rxn in raw_data.modificationReactions:
		reactants = []
		products = []
		for mol in rxn['stoichiometry']:
			if mol['coeff'] == -1:
				reactants += ['%s[%s]' % (mol['molecule'], mol['location'])]
			else:
				products += ['%s[%s]' % (mol['molecule'], mol['location'])]

		for trna, ctrna in zip(trnas, filtered_charged):
			if trna in reactants and ctrna in products:
				mass = 0
				for reactant in reactants:
					if reactant in mol_names:
						mass += mws[np.where(mol_names == reactant)[0][0]].asNumber()
					else:
						print('could not get mass for %s' % (reactant,))
				for product in products:
					if product == ctrna:
						continue

					if product in mol_names:
						mass -= mws[np.where(mol_names == product)[0][0]].asNumber()
					else:
						print('could not get mass for %s' % (product,))

				writer.writerow({
					'id': ctrna[:-3],
					'mw': mass,
					'location': [ctrna[-2]],
					})
