"""
This file compares the Moffitt 2016 mRNA half lives dataset with the Bernstein
2002 dataset which is used in this study (correlation plot reported in Figure S2).

Optional: This file also prepares the Moffitt 2016 mRNA half lives dataset for
use in the whole-cell framework. The output of this file
("rnas_alternate_half_lives_without_kas.tsv") can be relocated to
reconstruction/ecoli/flat/rnas_alternate_half_lives_without_kas.tsv, and the
alternate-rna-half-life option can be set when running the fitter by either:

python runscripts/manual/runFitter.py --alternate-rna-half-life

or

DESC="Simulation description" \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_RNA_HALF_LIFE=1 python runscripts/fw_queue.py

Sources:
Moffitt et al. "Spatial organization shapes the turnover of a bacterial
transcriptome". 2016. Figure 4 - source data 1. Column "Decay rate (1/min)".

Bernstein et al. "Global analysis of mRNA decay and abundance in Escherichia
coli at single-gene resolution using two-color fluorescent DNA microarrays".
2002. Table S5. Column "Hlaf-lives in M9 medium, min".
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

COLOR_HIGHLIGHT = "tab:orange"
PLOT_ALL_GENES = False

GENES_RNAP = ["rpoA", "rpoB", "rpoC"]
GENES_RIBO = ["rpsM", "rpsQ", "rpsU", "rpsO", "rpsJ", "rpsD", "rpsL", "rpsE",
              "rpsA", "rpsB", "rpsF", "rpsP", "rpsS", "rpsT", "rpsH", "sra",
              "rpsI", "rpsG", "rpsK", "rpsC", "rpsR", "rpsN", "rplK", "rplR",
              "rpmJ", "rplP", "rplO", "rpmH", "rplM", "rplU", "rplN", "rplX",
              "rpmI", "rpmC", "rplJ", "rplQ", "rpmB", "rplI", "rpmE", "rpmG",
              "rpmD", "rpmA", "rplF", "rplV", "rplW", "rplY", "rpmF", "rplA",
              "rplT", "rplB", "rplE", "rplS", "rplD", "rplL", "rplC"]

SYNONYMS = {
	'rapA': 'hepA',
	'thiB': 'tbpA',
	'yacF': 'zapD',
	'skp': 'hlpA',
	'yaeJ': 'arfB',
	'insX': 'ykgN',
	'ykgT': 'yagP',
	'yagV': 'ecpE',
	'yagW': 'ecpD',
	'yagX': 'ecpC',
	'matC': 'ecpB',
	'matB': 'ecpA',
	'matA': 'ecpR',
	'ybaZ': 'atl',
	'kefA': 'mscK',
	'chiX': 'micM',
	'selU': 'ybbB',
	'ylbA': 'allE',
	'exoD': 'ybcC',
	'peaD': 'ybcD',
	'quuD': 'ybcQ',
	'nohD': 'nohB',
	'ybeB': 'rsfS',
	'chiQ': 'ybfN',
	'ybhT': 'acrZ',
	'ybhO': 'clsB',
	'ybjK': 'rcdA',
	'ymdC': 'clsC',
	'opgC': 'mdoC',
	'opgG': 'mdoG',
	'opgH': 'mdoH',
	'ycfQ': 'comR',
	'cohE': 'ymfK',
	'croE': 'ymfT',
	'beeE': 'ymfO',
	'jayE': 'ymfP',
	'stfP': 'ycfK',
	'tfaP': 'ymfS',
	'ycgE': 'bluR',
	'ycgF': 'bluF',
	'ychM': 'dauA',
	# 'insZ': ['ychG_1', 'ychG_2'],
	'ycjZ': 'pgrR',
	'isrA': 'mcaS',
	'ydaC': 'rcbA',
	'ralR': 'lar',
	'opgD': 'mdoD',
	# 'insP': ['yncK_1', 'yncK_2'],
	'insQ': 'ydcM',
	'prr': 'patD',
	'yncA': 'mnaT',
	'yneO': 'ydeK',
	'stfQ': 'ydfN',
	'nohQ': 'nohA',
	'quuQ': 'ydfT',
	'intK': 'ydfW',
	'ydfC': 'rzpQ',
	'dgsA': 'mlc',
	'dtpA': 'tppB',
	'ydiI': 'menI',
	'cdgR': 'ydiV',
	'yeaZ': 'tsaB',
	'adrB': 'yoaD',
	'yebN': 'mntP',
	'ryeB': 'sdsR',
	'hiuH': 'yedX',
	'yeeU': 'cbeA',
	'cld': 'wzzB',
	'rfc': 'wbbH',
	# 'yehH': ['molR_1', 'molR_2', 'molR_3'],
	'preB': 'preT',
	'rhmD': 'yfaW',
	'yfcM': 'epmC',
	'gtrA': 'hemA',
	'gtrB': 'yfdH',
	'gtrS': 'yfdI',
	'oweS': 'yfdO',
	'yfiD': 'grcA',
	'yfiQ': 'pka',
	# 'yfjV': ['ypjM_3', 'ypjM_2', 'ypjM_1'],
	# 'ygaD': ['pncC', 'nrdF'],
	'ygbF': 'cas2',
	'ygbT': 'cas1',
	'ygcB': 'cas3',
	'csdL': 'tcdA',
	'ygfU': 'uacT',
	'ygfX': 'cptA',
	'ygfY': 'cptB',
	'yggG': 'loiP',
	'gpr': 'yghZ',
	'plsY': 'ygiH',
	'ygjD': 'tsaD',
	'yrbA': 'ibaG',
	'yhbJ': 'rapZ',
	'rimN': 'tsaC',
	'yhhK': 'panZ',
	'yhiQ': 'rsmJ',
	'yhiR': 'rlmJ',
	'yiaQ': 'sgbH',
	'yiaR': 'sgbU',
	'yiaS': 'sgbE',
	'yibD': 'waaH',
	'yibB': 'htrL',
	'rfaF': 'waaF',
	'rfaC': 'waaC',
	'rfaL': 'waaL',
	'rfaZ': 'waaZ',
	'rfaY': 'waaY',
	'rfaJ': 'waaJ',
	'rfaI': 'waaI',
	'rfaB': 'waaB',
	'rfaS': 'waaS',
	'rfaP': 'waaP',
	'rfaG': 'waaG',
	'rfaQ': 'waaQ',
	'yijP': 'eptC',
	'yjbD': 'pagB',
	'qorA': 'qor',
	'dsbD': 'dipZ',
	'yjeK': 'epmB',
	'yjeA': 'epmA',
	'yjeP': 'mscM',
	'yjeE': 'tsaE',
	'qorB': 'ytfG',
	'ytfM': 'tamA',
	'ytfN': 'tamB',
	'yjgF': 'ridA',
	'yjgI': 'bdcA',
	'yjgJ': 'bdcR',
	'yjgB': 'ahr',
	'qseD': 'hypT',
	'opgB': 'mdoB',
}

def openfile(filename):
	with open(filename, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		data = np.array([line for line in reader])
	return data

# Describe paths
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
this_dir = os.path.dirname(os.path.abspath(__file__))

moffitt_2016_file = os.path.join(this_dir, "Moffitt_2016.tsv")
bernstein_2002_file = os.path.join(this_dir, "Bernstein_2002.tsv")
rnas_file = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "rnas.tsv")
genes_file = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "genes.tsv")
synonyms_file = os.path.join(root_dir, "validation", "ecoli", "flat", "geneIDs.tsv")
output_plot_file = os.path.join(this_dir, "compare_mrna_half_lives{}.{}")
output_data_file = os.path.join(this_dir, "rnas_alternate_half_lives_without_kas.tsv")

# Get data
moffitt_2016 = openfile(moffitt_2016_file)
bernstein_2002 = openfile(bernstein_2002_file)

# Moffitt 2016: Average replicates
R1_LABEL = 'WT -kas replicate 1'
R2_LABEL = 'WT -kas replicate 2'
r1_indexes = np.where(moffitt_2016 == R1_LABEL)[0]
r2_indexes = np.where(moffitt_2016 == R2_LABEL)[0]

gene_ids_moffitt_2016 = moffitt_2016[r1_indexes, 1]
decay_rate_moffitt_2016_raw = np.zeros([len(gene_ids_moffitt_2016), 2])

for i, (i_r1, i_r2) in enumerate(zip(r1_indexes, r2_indexes)):
	assert moffitt_2016[i_r1, 1] == moffitt_2016[i_r2, 1]
	decay_rate_moffitt_2016_raw[i, 0] = moffitt_2016[i_r1, -2] if moffitt_2016[i_r1, -2] != '' else np.nan
	decay_rate_moffitt_2016_raw[i, 1] = moffitt_2016[i_r2, -2] if moffitt_2016[i_r2, -2] != '' else np.nan

decay_rate_moffitt_2016 = np.nanmean(decay_rate_moffitt_2016_raw, axis=1)
half_lives_moffitt_2016 = np.log(2) / decay_rate_moffitt_2016

# Bernstein
gene_ids_bernstein_2002 = bernstein_2002[9:, 1]
half_lives_bernstein_2002_raw = bernstein_2002[9:, np.where(bernstein_2002 == 'Half-lives in M9 ')[1][0]]
half_lives_bernstein_2002 = np.array([float(x) if x != '' else np.nan for x in half_lives_bernstein_2002_raw])

# Unify genes
synonyms_raw = openfile(synonyms_file)
synonyms = np.full([synonyms_raw.shape[0] -1, 23], None)
for i, row in enumerate(synonyms_raw[1:]):
	names = [x[1: -1] for x in row[1][1: -1].split(" ")]
	for j, name in enumerate(names):
		synonyms[i, j] = name

replace = {
	"insL": "insL-1",
	"insB": "insB-1",
	"insA": "insA-1",
	"ykfN": "b4626",
	"insI": "insI-1",
	"insH": "insH-1",
	"ykgS": "b4688",
	# "ptwF": "b4629", # discontinued
	# "ykgT": "b4695", # discontinued
	"ykgP": "b4630",
	"insE": "insE-1",
	"insF": "insF-1",
	"insC": "insC-1",
	"ylbI": "b4632",
	"selU": "b0503",
	"xisD": "b4633",
	"exoD": "b0539",
	"peaD": "b4508",
	"nohD": "b0560",
	"aaaD": "b4634",
	# "pauD": "b4635", # discontinued
	"sokE": "b4700",
	"ybfI": "b4636",
	"efeU": "efeU_1",
	"oweE": "b4692",
	"aaaE": "b4693",
	"ymiB": "b4672",
	"ttcC": "b4638",
	"insP": "yncK_1",
	"rzoQ": "b4689",
	"yeeH": "b4639",
	"yoeG": "b4640",
	"yoeH": "b4641",
	"yoeD": "b4642",
	"pawZ": "b4643",
	"ypjI": "b4644",
	"psaA": "b4645",
	"yrdF": "b4697",
	"yrdE": "b4646",
	"mokA": "b4647",
	"ysaC": "b4648",
	"ysaD": "b4649",
	"yibS": "b4650",
	"yibW": "b4651",
	"yibU": "b4652",
	"yicT": "b4653",
	"istR": "istR-1",
	"yjdQ": "b4654",
	"yjhY": "b4656",
	"yjhZ": "b4657",
}

for i, gene_id in enumerate(gene_ids_moffitt_2016):
	if gene_id in gene_ids_bernstein_2002:
		continue

	if gene_id in replace:
		gene_id = replace.get(gene_id)

	if gene_id in synonyms:
		for synonym in synonyms[np.where(synonyms == gene_id)[0][0]]:
			if synonym in gene_ids_bernstein_2002:
				gene_ids_moffitt_2016[i] = synonym
				break

gene_ids = list(set(gene_ids_moffitt_2016).intersection(set(gene_ids_bernstein_2002)))

# Compare datasets
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
for gene_id in gene_ids:
	moffitt_2016_val = half_lives_moffitt_2016[np.where(gene_ids_moffitt_2016 == gene_id)[0][0]]
	bernstein_2002_val = half_lives_bernstein_2002[np.where(gene_ids_bernstein_2002 == gene_id)[0][0]]

	if np.isnan(moffitt_2016_val) or np.isnan(bernstein_2002_val):
		continue

	if gene_id in GENES_RNAP:
		ax.plot(np.log10(moffitt_2016_val), np.log10(bernstein_2002_val),
		        color=COLOR_HIGHLIGHT, marker="o", markersize=10)
		print(gene_id)
	elif gene_id in GENES_RIBO:
		ax.plot(np.log10(moffitt_2016_val), np.log10(bernstein_2002_val),
		        color="tab:blue", marker="o", markersize=10)
	else:
		if PLOT_ALL_GENES:
			ax.plot(np.log10(moffitt_2016_val), np.log10(bernstein_2002_val),
			        color="tab:green", marker="o", alpha=0.2, markersize=10)
		else:
			pass

range_min = np.floor(min(ax.get_xlim()[0], ax.get_ylim()[0]) * 10)/10.
range_max = np.ceil(max(ax.get_xlim()[1], ax.get_ylim()[1]) * 10)/10.
ax.plot(
		[range_min, range_max],
		[range_min, range_max], "k")
ax.set_xlim(range_min, range_max)
ax.set_ylim(range_min, range_max)
ax.set_xticks([range_min, range_max])
ax.set_yticks([range_min, range_max])
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
plt.savefig(output_plot_file.format("__clean", "pdf"))

ax.set_xticklabels(ax.get_xticks())
ax.set_yticklabels(ax.get_yticks())
from matplotlib.ticker import FormatStrFormatter
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.set_xlabel("log10 Moffitt et al. 2016 (min)")
ax.set_ylabel("log10 Bernstein et al. 2002 (min)")
plt.legend(
	handles=[
		mlines.Line2D([], [], color="tab:red", linewidth=0., marker=".", label="RNA polymerase"),
		mlines.Line2D([], [], color="tab:blue", linewidth=0., marker=".", label="Ribosome")],
	loc="best")
plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75)
plt.savefig(output_plot_file.format("", "pdf"))
plt.close("all")

# Prepare Moffitt 2016 mRNA half lives dataset for wcEcoli
# Note: stable RNA half lives are unchanged (ie. from Bernstein 2002)

genes = openfile(genes_file)
rnas = openfile(rnas_file)

eg_to_name = {x[-2]: x[-4] for x in genes[1:]}
eg_to_type = {x[-2]: x[3] for x in rnas[1:]}
name_to_type = {x[-4]: eg_to_type.get(x[-2]) for x in genes[1:]}

half_lives_moffitt_2016_sec = half_lives_moffitt_2016 * 60.
gene_types_moffitt_2016 = np.array([name_to_type.get(x) for x in gene_ids_moffitt_2016])
half_lives_moffitt_2016_avg = np.nanmean(half_lives_moffitt_2016_sec[gene_types_moffitt_2016 == "mRNA"])

out = ['"id"\t"halfLife"']

for row in rnas[1:, :]:
	gene_id = row[-2]
	rna_id = row[-3]
	rna_type = row[3]

	if rna_type == "mRNA":
		gene_name = eg_to_name.get(gene_id)

		if gene_name not in gene_ids or np.isnan(half_lives_moffitt_2016_sec[np.where(gene_ids_moffitt_2016 == gene_name)[0][0]]):
			# Assign average half life for missing half lives
			half_life = half_lives_moffitt_2016_avg
		else:
			half_life = half_lives_moffitt_2016_sec[np.where(gene_ids_moffitt_2016 == gene_name)[0][0]]
	else:
		half_life = row[0]

	out.append('"{}"\t{}'.format(rna_id, half_life))

# Save output
with open(output_data_file, "w") as f:
	f.write("\n".join(out))

