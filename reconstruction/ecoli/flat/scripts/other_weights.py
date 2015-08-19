"""
called to update various molecular weights
- polymerized monomers
- proteins
- RNAs
- complexes

should be idempotent (i.e. no change if called more than once)
"""

import os
from functools import partial
from collections import Counter, OrderedDict

import numpy as np

from reconstruction import spreadsheets

def flip_dict(dct):
	return {value:key for key, value in dct.viewitems()}

DIALECT = "excel-tab"
KEY = "id"

FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
MET_FILE = os.path.join(FLAT_DIR, "metabolites.tsv")
WATER_FILE = os.path.join(FLAT_DIR, "water.tsv")
POLY_FILE = os.path.join(FLAT_DIR, "polymerized.tsv")
PROT_FILE = os.path.join(FLAT_DIR, "proteins.tsv")
RNA_FILE = os.path.join(FLAT_DIR, "rnas.tsv")
COMP_FILE = os.path.join(FLAT_DIR, "proteinComplexes_large.tsv")
COMP_RXN_FILE = os.path.join(FLAT_DIR, "complexationReactions_large.tsv")

ID_DIR = os.path.join(FLAT_DIR, "ids")
NTPS_FILE = os.path.join(ID_DIR, "ntps.txt")
DNTPS_FILE = os.path.join(ID_DIR, "dntps.txt")
AAS_FILE = os.path.join(ID_DIR, "amino_acids.txt")

SYMBOL_TO_ID = OrderedDict((
	("A", "L-ALPHA-ALANINE[c]"), ("R", "ARG[c]"), ("N", "ASN[c]"), ("D", "L-ASPARTATE[c]"),
	("C", "CYS[c]"), ("E", "GLT[c]"), ("Q", "GLN[c]"), ("G", "GLY[c]"),
	("H", "HIS[c]"), ("I", "ILE[c]"), ("L", "LEU[c]"), ("K", "LYS[c]"),
	("M", "MET[c]"), ("F", "PHE[c]"), ("P", "PRO[c]"), ("S", "SER[c]"),
	("T", "THR[c]"), ("W", "TRP[c]"), ("Y", "TYR[c]"), ("U", "L-SELENOCYSTEINE[c]"),
	("V", "VAL[c]")
	))

AA_SYM_ORDER = {
	key:ind
	for ind, key in enumerate(SYMBOL_TO_ID.keys())
	}

AA_ORDER = {
	value[:-3]:ind
	for ind, value in enumerate(SYMBOL_TO_ID.values())
	}

NTP_ORDER = {
	"ATP":0,
	"CTP":1,
	"GTP":2,
	"UTP":3
	}

MW_KEYS = [ # TODO: add flat file for this, and load here/in sim data
	'23srRNA',
	'16srRNA',
	'5srRNA',
	'tRNA',
	'mRNA',
	'miscRNA',
	'protein',
	'metabolite',
	'water',
	'DNA',
	'RNA'
	]

N_MW = len(MW_KEYS)

COMPARTMENTS = ["n", "j", "w", "c", "e", "m", "o", "p", "l", "i"] # TODO: also to flat file

JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

def load_ids(id_file):
	return [s.strip() for s in open(id_file) if s.strip()]


def lod_to_dod(list_of_dicts, dict_key):
	# function name of the day award
	entries = [
		(dct[dict_key], dct) for dct in list_of_dicts
		]

	out = dict(entries)

	if len(entries) != len(out):
		raise Exception("Repeated keys (non-unique)")

	return out

ntp_ids = load_ids(NTPS_FILE)
dntp_ids = load_ids(DNTPS_FILE)
aa_ids = load_ids(AAS_FILE)

met = lod_to_dod(JsonReader(open(MET_FILE)), KEY)
water = lod_to_dod(JsonReader(open(WATER_FILE)), KEY)

# RNA: NTP - diphosphate = polymerized nucleotide (consistent w/ process)
# DNA: same as RNA
# protein: unclear, check process (AA - H2O = poly'd)

poly = []
ntp_weights = np.zeros(len(NTP_ORDER))

ppi_weight = met["PPI"]["mw7.2"]
rna_index = MW_KEYS.index("RNA")
for ntp_id in ntp_ids:
	weight = met[ntp_id]["mw7.2"] - ppi_weight

	mw = np.zeros(N_MW)
	mw[rna_index] = weight

	poly_id = "polymerized {}".format(ntp_id)
	poly.append({
		KEY:poly_id,
		"location":COMPARTMENTS,
		"mw":mw.tolist(),
		"is_aa":False,
		"is_ntp":True,
		"is_dntp":False,
		"is_end":False
		})

	ntp_weights[NTP_ORDER[ntp_id]] = weight

dna_index = MW_KEYS.index("DNA")
for dntp_id in dntp_ids:
	mw = np.zeros(N_MW)
	mw[dna_index] = met[dntp_id]["mw7.2"] - ppi_weight

	poly_id = "polymerized {}".format(dntp_id)
	poly.append({
		KEY:poly_id,
		"location":COMPARTMENTS,
		"mw":mw.tolist(),
		"is_aa":False,
		"is_ntp":False,
		"is_dntp":True,
		"is_end":False
		})

ntp_terminal_weight = ppi_weight

mw = np.zeros(N_MW)
mw[rna_index] = ntp_terminal_weight
poly.append({
	KEY:"nucleic acid terminal pyrophosphate",
	"location":COMPARTMENTS,
	"mw":mw.tolist(),
	"is_aa":False,
	"is_ntp":True,
	"is_dntp":True,
	"is_end":True
	})

aa_weights = np.zeros(len(AA_ORDER))

water_weight = water["WATER"]["mw7.2"]
protein_index = MW_KEYS.index("protein")
for aa_id in aa_ids:
	weight = met[aa_id]["mw7.2"] - water_weight

	mw = np.zeros(N_MW)
	mw[protein_index] = weight

	poly_id = "polymerized {}".format(aa_id)
	poly.append({
		KEY:poly_id,
		"location":COMPARTMENTS,
		"mw":mw.tolist(),
		"is_aa":True,
		"is_ntp":False,
		"is_dntp":False,
		"is_end":False
		})

	aa_weights[AA_ORDER[aa_id]] = weight

aa_terminal_weight = water_weight

mw = np.zeros(N_MW)
mw[protein_index] = aa_terminal_weight
poly.append({
	KEY:"polypeptide terminal hydroxyl",
	"location":COMPARTMENTS,
	"mw":mw.tolist(),
	"is_aa":True,
	"is_ntp":False,
	"is_dntp":False,
	"is_end":True
	})

species_weights = {}

met_index = MW_KEYS.index("metabolite")
for key, value in met.viewitems():
	w = np.zeros(N_MW)
	w[met_index] = value["mw7.2"]
	species_weights[key] = np.array(w)

water_index = MW_KEYS.index("water")
for key, value in water.viewitems():
	w = np.zeros(N_MW)
	w[water_index] = value["mw7.2"]
	species_weights[key] = np.array(value["mw7.2"])

del met
del water

with open(POLY_FILE, "r") as f:
	reader = JsonReader(f)
	fieldnames = reader.fieldnames

with open(POLY_FILE, "w") as f:
	writer = JsonWriter(f, fieldnames)

	writer.writeheader()
	writer.writerows(poly)

del poly

# update RNAs

with open(RNA_FILE, "r") as f:
	reader = JsonReader(f)
	fieldnames = reader.fieldnames

	rna_data = lod_to_dod(reader, KEY)

for rna_id, rna_entry in rna_data.viewitems():
	new_weight = ntp_terminal_weight + np.dot(ntp_weights, rna_entry["ntCount"])

	mw = rna_entry["mw"]
	weight_index = np.where(mw)[0]

	assert weight_index.size == 1

	mw[weight_index] = new_weight
	species_weights[rna_id] = np.array(mw)

with open(RNA_FILE, "w") as f:
	writer = JsonWriter(f, fieldnames)
	writer.writeheader()

	for key in sorted(rna_data.keys()):
		writer.writerow(rna_data[key])

del rna_data

# update proteins

with open(PROT_FILE, "r") as f:
	reader = JsonReader(f)
	fieldnames = reader.fieldnames

	prot_data = lod_to_dod(reader, KEY)

for prot_id, prot_entry in prot_data.viewitems():
	aa_counts = np.zeros(aa_weights.size, np.int64)
	for c in prot_entry["seq"]:
		aa_counts[AA_SYM_ORDER[c]] += 1

	prot_entry["aaCount"] = aa_counts.tolist()
	new_weight = aa_terminal_weight + np.dot(aa_weights, aa_counts)

	mw = prot_entry["mw"]
	weight_index = np.where(mw)[0]

	assert weight_index.size == 1

	mw[weight_index] = new_weight
	species_weights[prot_id] = np.array(mw)

with open(PROT_FILE, "w") as f:
	writer = JsonWriter(f, fieldnames)
	writer.writeheader()

	for key in sorted(prot_data.keys()):
		writer.writerow(prot_data[key])

del prot_data

# update complexes

with open(COMP_FILE, "r") as f:
	reader = JsonReader(f)
	fieldnames = reader.fieldnames

	comp_data = lod_to_dod(reader, KEY)

with open(COMP_RXN_FILE, "r") as f:
	reader = JsonReader(f)
	comp_rxns = lod_to_dod(reader, KEY)

while comp_rxns:
	to_remove = set()

	for comp_rxn_id, comp_rxn in comp_rxns.viewitems():
		stoich = {s["molecule"]:s["coeff"] for s in comp_rxn["stoichiometry"]}

		subunits = set(mid for mid, c in stoich.viewitems() if c < 0)

		(comp_id,) = [mid for mid, c in stoich.viewitems() if c > 0]

		if subunits <= set(species_weights.viewkeys()):
			weight = np.zeros(N_MW)
			for subunit in subunits:
				weight += -stoich[subunit] * species_weights[subunit]

			species_weights[comp_id] = weight
			comp_data[comp_id]["mw"] = weight.tolist()
			to_remove.add(comp_rxn_id)

	for rxn_id in to_remove:
		del comp_rxns[rxn_id]

	if len(to_remove) == 0:
		unrecognized = {
			s["molecule"]
			for comp_rxn in comp_rxns.viewvalues()
			for s in comp_rxn["stoichiometry"]
			} - species_weights.viewkeys() - comp_data.viewkeys()

		raise Exception("{} unrecognized subunits: {}".format(len(unrecognized), ",".join(unrecognized)))

with open(COMP_FILE, "w") as f:
	writer = JsonWriter(f, fieldnames)
	writer.writeheader()

	for key in sorted(comp_data.keys()):
		writer.writerow(comp_data[key])

del comp_data
