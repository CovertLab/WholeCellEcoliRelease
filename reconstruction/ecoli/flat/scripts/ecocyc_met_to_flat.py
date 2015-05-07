
from __future__ import division

import os
import re

from collections import defaultdict

from reconstruction.spreadsheets import JsonWriter

_DIR = os.path.join("reconstruction", "ecoli", "flat", "metabolism")
_FBA_FILE = os.path.join(_DIR, "ecocyc-full-biomass.fba") # biomass stoich, nutrients, secretions
_LP_FILE = os.path.join(_DIR, "ecocyc-full-biomass.lp") # reaction stoich, constraints, objective
_LPLOG_FILE = os.path.join(_DIR, "ecocyc-full-biomass.lplog") # rxn/compound mapping

_COMPARTMENT_CHAR = {
	"CCO-CYTOSOL":"c",
	"CCO-PM-BAC-NEG":"i", # inner membrance
	"CCO-PERI-BAC":"p",
	"CCO-EXTRACELLULAR":"e",
	# TODO: resolve these w/ compartments.tsv
	"CCO-OUT":"0", # not in CCO, related to N/GAM?
	"CCO-IN":"1", # not in CCO, related to N/GAM?
	"CCO-UNKNOWN-SPACE":"2", # not in CCO, related to cell envelope
	"CCI-PERI-BAC-GN":"3", # not in CCO, related to cell envelope
	}

def reduce_compartment(molecule_id):
	begin_compartment = molecule_id.find("[")
	compartment = molecule_id[begin_compartment+1:-1]
	return "{}[{}]".format(
		molecule_id[:begin_compartment],
		_COMPARTMENT_CHAR[compartment]
		)

# Load reaction/compound ID mapping

reaction_id_to_name = {}
reversible_reactions = set()
exchange_reactions = set()

compound_id_to_name = {}

for line in open(_LPLOG_FILE):
	try:
		lp_id, name = line.strip("\n").split(",")

	except ValueError:
		continue

	if lp_id.startswith("r"):
		reaction_id_to_name[lp_id] = name

		if lp_id.endswith("a") or lp_id.endswith("b"):
			reversible_reactions.add(lp_id)

		if lp_id.endswith("_e"):
			exchange_reactions.add(lp_id)

	elif lp_id.startswith("c"):
		compound_id_to_name[lp_id] = reduce_compartment(name)

	else:
		raise Exception("Unparsable line:\n{}".format(line))

# Load reaction stoich

started_constraints = False
compound_id = None

reaction_stoich = defaultdict(dict)

for line in open(_LP_FILE):
	if not started_constraints:
		if line == "Subject To\n":
			started_constraints = True

		else:
			continue

	if compound_id is None:
		match = re.match("c[0-9]+", line)

		if match is not None:
			compound_id = match.group()

	else:
		stoich_raw = re.findall("([+-][0-9.]*) (r[0-9]+[ab\_e]*)", line)

		if not stoich_raw:
			continue

		for stoich_str, reaction_id in stoich_raw:
			stoich = 1.0 if (stoich_str == "+") else float(stoich_str)

			if int(stoich) != stoich:
				raise Exception("non-float reaction stoich: {} {} ({})".format(
					stoich,
					compound_id,
					reaction_id
					))

			stoich = int(stoich)

			reaction_dict = reaction_stoich[reaction_id]

			if compound_id in reaction_dict:
				assert reaction_dict[compound_id] == stoich, "{}[{}]: {} != {}".format(reaction_id, compound_id, reaction_dict[compound_id], stoich)

			else:
				reaction_dict[compound_id] = stoich

		if re.findall("= 0", line):
			compound_id = None

# Check that all singleton "reactions" are exchange reactions
assert all(
	(rxn_id in exchange_reactions)
	for rxn_id, stoich in reaction_stoich.viewitems()
	if len(stoich) == 1
	)

# Transform stoich from LP IDs to Ecocyc IDs

reaction_stoich = {
	reaction_id_to_name[reaction_id]: {
		compound_id_to_name[compound_id]:coeff
		for compound_id, coeff in stoich.viewitems()
		}
	for reaction_id, stoich in reaction_stoich.viewitems()
	if reaction_id not in exchange_reactions
	}

# Write a table of true reactions:
# reaction name, stoich, is_reversible

with open(os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv"), "w") as outfile:
	writer = JsonWriter(
		outfile,
		["reaction id", "stoichiometry", "is reversible"],
		dialect = "excel-tab"
		)

	writer.writeheader()

	for reaction_id, stoich in reaction_stoich.viewitems():
		writer.writerow({
			"reaction id":reaction_id,
			"stoichiometry":stoich,
			"is reversible":(reaction_id in reversible_reactions)
			})

# Parse the .fba file

current_group = None

data = {}

with open(_FBA_FILE) as fba_file:
	for line in fba_file:
		split = line.split()

		if len(split) == 0:
			continue

		if line.startswith("#"):
			pass

		elif line.startswith(":group"):
			assert current_group is None
			assert len(split) == 2

			current_group = split[-1]

		elif line.startswith(":end-group"):
			current_group = None

		elif line[0].isalpha():
			if split[0].endswith(":"): # section heading
				current_section = split[0][:-1]

				assert current_section not in data.viewkeys()

				data[current_section] = []

				if len(split) > 1:
					data[current_section].append(" ".join(split[1:]))

			elif current_section in {"reactions", "remove-reactions"}:
				data[current_section].append(split[0])

			elif current_section == "biomass": # biomass component
				data[current_section].append({
					"group id": current_group if current_group is not None else "",
					"molecule id": reduce_compartment(split[0]),
					"coefficient": float(split[1])
					})

			elif current_section in {"nutrients", "secretions"}:
				try:
					lb_at = split.index(":lower-bound")

				except ValueError:
					lower_bound = None

				else:
					lower_bound = float(split[lb_at+1])

				try:
					ub_at = split.index(":upper-bound")

				except ValueError:
					upper_bound = None

				else:
					upper_bound = float(split[ub_at+1])

				data[current_section].append({
					"molecule id":reduce_compartment(split[0]),
					"lower bound":lower_bound,
					"upper bound":upper_bound
					})

			else:
				print "could not parse line:", line.strip()

# Write a table of biomass components:
# group name, metabolite name, coeff

with open(os.path.join("reconstruction", "ecoli", "flat", "biomass.tsv"), "w") as outfile:
	writer = JsonWriter(
		outfile,
		["group id", "molecule id", "coefficient"],
		dialect = "excel-tab"
		)

	writer.writeheader()
	writer.writerows(data["biomass"])

# Write a table of exchanges and constraints:
# metabolite name, is_ngam, lb, ub
for section_name in ("nutrients", "secretions"):
	with open(os.path.join("reconstruction", "ecoli", "flat", section_name + ".tsv"), "w") as outfile:
		writer = JsonWriter(
			outfile,
			["molecule id", "lower bound", "upper bound"],
			dialect = "excel-tab"
			)

		writer.writeheader()
		writer.writerows(data[section_name])
