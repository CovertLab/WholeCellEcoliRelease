
from __future__ import division

import os
import re

from collections import defaultdict

from reconstruction.spreadsheets import JsonWriter

_DIR = os.path.join("reconstruction", "ecoli", "flat", "metabolism")
_FBA_FILE = os.path.join(_DIR, "ecocyc-full-biomass.fba") # biomass stoich, nutrients, secretions
_LP_FILE = os.path.join(_DIR, "ecocyc-full-biomass.lp") # reaction stoich, constraints, objective
_LPLOG_FILE = os.path.join(_DIR, "ecocyc-full-biomass.lplog") # rxn/compound mapping

# TODO: check this
_COMPARTMENT_CHAR = {
	"CCO-CYTOSOL":"c",
	"CCO-PM-BAC-NEG":"i", # inner membrance
	"CCO-PERI-BAC":"p",
	"CCO-EXTRACELLULAR":"e",
	"CCO-OUT":"0", # not in CCO, related to N/GAM?
	"CCO-IN":"1", # not in CCO, related to N/GAM?
	"CCO-UNKNOWN-SPACE":"2", # not in CCO, related to cell envelope
	"CCI-PERI-BAC-GN":"3", # not in CCO, related to cell envelope
	}

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
		begin_compartment = name.find("[")
		compartment = name[begin_compartment+1:-1]
		compound_id_to_name[lp_id] = "{}[{}]".format(
			name[:begin_compartment],
			_COMPARTMENT_CHAR[compartment]
			)

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
			stoich = 1.0 if stoich_str == "+" else float(stoich_str)

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

# Write a table of biomass components:
# group name, metabolite name, coeff

# Write a table of exchanges and constraints:
# metabolite name, is_ngam, lb, ub

import ipdb; ipdb.set_trace()
