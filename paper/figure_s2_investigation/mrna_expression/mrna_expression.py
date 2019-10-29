"""
This file compares the Covert 2004 RNA-seq dataset with the one
generated in this study (correlation plot reported in Figure S2).

Optional: This file also prepares the Covert 2004 RNA-seq dataset for use in
the whole-cell framework. The output of this file ("alternate_rna_seq.tsv") can
be relocated to reconstruction/ecoli/flat/rna_seq_data/alternate_rna_seq.tsv,
and the alternate-rna-seq option can be set when running the fitter by either:

python runscripts/manual/runFitter.py --alternate-rna-seq

or

DESC="Simulation description" \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_RNA_SEQ=1 python runscripts/fw_queue.py

Source:
Covert et al. "Integrating high-throughput and computational data elucidates
bacterial networks". 2004. Supplementary Data 7.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import pearsonr

from reconstruction.ecoli.dataclasses.process.transcription import RNA_SEQ_ANALYSIS

PLOT_ALL_GENES = False

BASAL_CONDITION = "M9 Glucose minus AAs"
GENES_RNAP = ["EG10893", "EG10894", "EG10895"]
GENES_RIBO = ["EG10912", "EG10916", "EG10920", "EG10914", "EG10909", "EG10903",
			"EG10911", "EG10904", "EG10900", "EG10901", "EG10905", "EG10915",
			"EG10918", "EG10919", "EG10907", "EG11508", "EG10908", "EG10906",
			"EG10910", "EG10902", "EG10917", "EG10913", "EG10872", "EG10879",
			"EG11232", "EG10877", "EG10876", "EG10892", "EG10874", "EG50001",
			"EG10875", "EG10884", "EG11231", "EG10887", "EG10871", "EG10878",
			"EG10886", "EG10870", "EG10889", "EG10891", "EG10888", "EG50002",
			"EG10869", "EG10882", "EG10883", "EG10885", "EG10890", "EG10864",
			"EG10881", "EG10865", "EG10868", "EG10880", "EG10867", "EG10873",
			"EG10866"]

COLOR_RIBOSOME = "tab:blue"
COLOR_RNAP = "tab:orange"

def openfile(filename):
	with open(filename, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		data = np.array([line for line in reader])
	return data

# Describe paths
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
this_dir = os.path.dirname(os.path.abspath(__file__))

rna_seq_this_study_file = os.path.join(root_dir, "reconstruction", "ecoli", "flat",
	"rna_seq_data", "rnaseq_{}_mean.tsv".format(RNA_SEQ_ANALYSIS))
rna_seq_covert_2004_file = os.path.join(this_dir, "Covert_2004.tsv")
gene_ids_file = os.path.join(this_dir, "gene_ids.tsv") #todo: use validation synonyms
rnas_file = os.path.join(root_dir, "reconstruction", "ecoli", "flat", "rnas.tsv")
output_plot_file = os.path.join(this_dir, "compare_mrna_expression{}.{}")
output_data_file = os.path.join(this_dir, "alternate_rna_seq.tsv")

# Get data
rna_seq_this_study = openfile(rna_seq_this_study_file)
rna_seq_covert_2004 = openfile(rna_seq_covert_2004_file)

# Unify gene IDs between data sets
# (this study uses EG numbers, Covert 2004 uses b numbers)
gene_ids = openfile(gene_ids_file)
b_to_eg = {gene_ids[1:, 1][i]: gene_ids[1:, 0][i] for i in range(len(gene_ids) -1) if gene_ids[1:, 1][i] != ""}
eg_to_b = {gene_ids[1:, 0][i]: gene_ids[1:, 1][i] for i in range(len(gene_ids) -1)if gene_ids[1:, 1][i] != ""}


# Prepare Covert 2004 RNA-seq dataset for wcEcoli
rnas = openfile(rnas_file)
eg_to_type = {rnas[1:, -2][i]: rnas[1:, 3][i] for i in range(len(rnas) -1)}
b_to_type = {b: eg_to_type.get(b_to_eg.get(b)) for b in b_to_eg.keys()}

covert_2004_rna_type = np.array([b_to_type.get(b) for b in rna_seq_covert_2004[:, 0]])

covert_2004_avgs = {
	"mRNA": rna_seq_covert_2004[covert_2004_rna_type == "mRNA", 1:].astype(float).mean(axis=1).mean(axis=0),
	"tRNA": rna_seq_covert_2004[covert_2004_rna_type == "tRNA", 1:].astype(float).mean(axis=1).mean(axis=0),
	"rRNA": rna_seq_covert_2004[covert_2004_rna_type == "rRNA", 1:].astype(float).mean(axis=1).mean(axis=0),
	"miscRNA": rna_seq_covert_2004[covert_2004_rna_type == "miscRNA", 1:].astype(float).mean(axis=1).mean(axis=0)}

out = ['"Gene"\t"{}"'.format(BASAL_CONDITION)]
x_values = []
y_values = []
fig, ax = plt.subplots(1, 1, figsize=(5, 5))

for i, gene_eg in enumerate(rna_seq_this_study[1:, 0]):
	gene_b = eg_to_b.get(gene_eg)

	if gene_b is None or gene_b not in rna_seq_covert_2004:
		# Assign (type-based) average expression for missing genes
		covert_2004_avg_exp = covert_2004_avgs.get(eg_to_type.get(gene_eg))

	else:
		covert_2004_index = np.where(rna_seq_covert_2004 == gene_b)[0][0]
		covert_2004_avg_exp = rna_seq_covert_2004[covert_2004_index, 1:].astype(float).mean()

		if PLOT_ALL_GENES:
			if eg_to_type.get(gene_eg) == "mRNA":
				if gene_eg in GENES_RNAP:
					ax.plot(
					np.log10(rna_seq_this_study[i +1, 1].astype(float)),
					np.log10(covert_2004_avg_exp),
					color=COLOR_RNAP, marker="o")

				elif gene_eg in GENES_RIBO:
					ax.plot(
					np.log10(rna_seq_this_study[i +1, 1].astype(float)),
					np.log10(covert_2004_avg_exp),
					color=COLOR_RIBOSOME, marker="o")

				else:
					ax.plot(
					np.log10(rna_seq_this_study[i +1, 1].astype(float)),
					np.log10(covert_2004_avg_exp),
					color="k", marker="o", alpha=0.2)


		if gene_eg in GENES_RNAP or gene_eg in GENES_RIBO:
			ax.plot(
				np.log10(rna_seq_this_study[i +1, 1].astype(float)),
				np.log10(covert_2004_avg_exp),
				color=COLOR_RNAP if gene_eg in GENES_RNAP else COLOR_RIBOSOME,
				marker="o", markersize=10)
			x_values.append(np.log10(rna_seq_this_study[i +1, 1].astype(float)))
			y_values.append(np.log10(covert_2004_avg_exp))

	out.append('"{}"\t{}'.format(gene_eg, covert_2004_avg_exp))

# Compute pearsonr
r_value, p_value = pearsonr(x_values, y_values)

# Save output
with open(output_data_file, "w") as f:
	f.write("\n".join(out))

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
ax.set_xlabel("log10 This Study (Macklin et al. 2019)")
ax.set_ylabel("log10 Covert et al. 2004")
plt.legend(	handles=[
	mlines.Line2D([], [], color=COLOR_RIBOSOME, linewidth=0., marker="o", label="Ribosome"),
	mlines.Line2D([], [], color=COLOR_RNAP, linewidth=0., marker="o", label="RNA Polymerase")],
	loc="best")
plt.title("r = {}\nR^2 = {}\np-value = {}".format(r_value, r_value**2, p_value))
plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75)
plt.savefig(output_plot_file.format("", "pdf"))
plt.close("all")
