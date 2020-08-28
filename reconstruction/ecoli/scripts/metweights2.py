from __future__ import absolute_import, division, print_function

import os
import re
from typing import Dict

from Bio.Data.IUPACData import atom_weights
import six

from reconstruction.spreadsheets import tsv_writer

# Constants

SOURCE = os.path.join("reconstruction", "ecoli", "flat", "metabolism", "ecocyc-compound-formulas.dat")
DELIMITER = "$"
REGEX = "([A-Z][a-z]?)([0-9]*)"
NDIGITS = 3 # number of decimals in atomic weights
COMPARTMENTS = ["n", "j", "w", "c", "e", "m", "o", "p", "l", "i"]
ADDED_SPECIES = {
	"glycogen-monomer":{"GLC":+1, "PROTON":+1, "WATER":-1},
	}
OUTPUT_METS = os.path.join("reconstruction", "ecoli", "flat", "metabolites.tsv")
OUTPUT_WATER = os.path.join("reconstruction", "ecoli", "flat", "water.tsv")

# Compute weights

no_formula = set()
atoms = set()
formulas = {}  # type: Dict[str, Dict[str, int]]

for line in open(SOURCE):
	molecule_name, formula_string = line.strip().split(DELIMITER)

	if formula_string:
		formula = {}  # type: Dict[str, int]
		for atom_name, count_string in re.findall(REGEX, formula_string):
			atoms.add(atom_name)

			assert atom_name not in formula
			formula[atom_name] = int(count_string) if count_string else 1

		if molecule_name in formulas:
			assert formulas[molecule_name] == formula

		else:
			formulas[molecule_name] = formula

	else:
		no_formula.add(molecule_name)

weights = {
	molecule_name:sum(count * round(atom_weights[atom_name], NDIGITS) for atom_name, count in six.viewitems(formula))
	for molecule_name, formula in six.viewitems(formulas)
	}

# Add misc. species
for molecule_name, stoich in six.viewitems(ADDED_SPECIES):
	assert molecule_name not in weights
	weights[molecule_name] = sum(
		coeff * weights[mol] for mol, coeff in six.viewitems(stoich)
		)

# Write out metabolites
with tsv_writer(OUTPUT_METS, ["id", "mw", "location"]) as writer:
	for molecule_name, weight in six.viewitems(weights):
		if molecule_name == "WATER":
			continue

		writer.writerow({
			"id":molecule_name,
			"mw":weight,
			"location":COMPARTMENTS
			})

# Write out water
with tsv_writer(OUTPUT_WATER, ["id", "mw", "location"]) as writer:
	for molecule_name, weight in six.viewitems(weights):
		if molecule_name != "WATER":
			continue

		writer.writerow({
			"id":molecule_name,
			"mw":weight,
			"location":COMPARTMENTS
			})
