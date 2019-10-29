"""
This file compares the Mohammad 2019 ribosome density dataset with the Li 2014
translation efficiency dataset which is used in this study (correlation plot
reported in Figure S2).

Optional: This file also prepares the Mohammad 2019 ribosome density dataset
for use in the whole-cell framework. The output of this file ("translationEfficiency_alternate.tsv") can be relocated to
reconstruction/ecoli/flat/translationEfficiency_alternate.tsv, and the
alternate-translation-efficiency option can be set when running the fitter by
either:

python runscripts/manual/runFitter.py --alternate-translation-efficiency

or

DESC="Simulation description" \
DISABLE_RIBOSOME_CAPACITY_FITTING=1 DISABLE_RNAPOLY_CAPACITY_FITTING=1 \
ALTERNATE_TRANSLATION_EFFICIENCY=1 python runscripts/fw_queue.py

Sources:
Mohammad et al. "A systematically-revised ribosome profiling method for bacteria
reveals pauses at single-codon resolution". 2019.
https://github.com/greenlabjhmi/2018_Bacterial_Pipeline_riboseq

Li et al. "Quantifying Absolute Protein Synthesis Rates Reveals Principles
Underlying Allocation of Cellular Resources". 2014. Table S4.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

COLOR_HIGHLIGHT = "tab:orange"

def open_file(filename, delimiter):
	with open(filename, "r") as f:
		reader = csv.reader(f, delimiter=delimiter)
		data = np.array([x for x in reader])
	return data

# Describe paths
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
this_dir = os.path.dirname(os.path.abspath(__file__))
mohammad_data_file = os.path.join(this_dir, "SRR7759814_pass.fastq.gz_genelist.csv")
li_data_file = os.path.join(this_dir, "li_table_s4.tsv")
wcm_data_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "rna_seq_data", "rnaseq_seal_rpkm_mean.tsv")
wcm_genes_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "genes.tsv")
ribosome_genes_file = os.path.join(this_dir, "ribosome_genes.tsv")
output_plot = os.path.join(this_dir, "compare_translation_efficiency{}.pdf")
output_file = os.path.join(
	root_dir, "reconstruction", "ecoli", "flat", "translationEfficiency_alternate.tsv")

# Identify genes of interest to highlight
rnap_genes = ["rpoA", "rpoB", "rpoC"]
ribosome_genes = [x[0] for x in open_file(ribosome_genes_file, "\n")]

# Mohammad et al. 2019 Ribosome Density (RPKM)
mohammad_data = np.array([x[0].split(",") for x in open_file(mohammad_data_file, ",")])
mohammad_data_header = mohammad_data[0].tolist()
mohammad_genes = mohammad_data[1:, mohammad_data_header.index("Alias")]
mohammad_ribosome_density = mohammad_data[1:, mohammad_data_header.index("RPKM")].astype(float)

# Li et al. 2014 Translation Efficiency (AU), mRNA (RPKM)
li_data = np.array(open_file(li_data_file, delimiter="\t"))
li_data_header = li_data[0].tolist()
li_genes = li_data[1:, li_data_header.index("Gene")]
li_mrna = np.array([x[1:-1] if "[" in x else x for x in li_data[1:, li_data_header.index("mRNA level (RPKM)")]], dtype=float)
li_translation_efficiency_raw = li_data[1:, li_data_header.index("Translation efficiency (AU)")]

# Exclude "NA" genes
mask = li_translation_efficiency_raw != "NA"
li_genes = li_genes[mask]
li_mrna = li_mrna[mask]
li_translation_efficiency = li_translation_efficiency_raw[mask].astype(float)

# Compute ribosome density
# Note: Translation Efficiency = Ribosome Density / mRNA
li_ribosome_density = li_translation_efficiency * li_mrna


# This study RNA (RPKM)
wcm_data = open_file(wcm_data_file, "\t")
wcm_data_header = wcm_data[0].tolist()
wcm_genes_egnumber = wcm_data[1:, wcm_data_header.index("Gene")]
wcm_rna = wcm_data[1:, wcm_data_header.index("M9 Glucose minus AAs")].astype(float)
wcm_genes_data = open_file(wcm_genes_file, "\t")
wcm_genes_data_header = wcm_genes_data[0].tolist()
egnumber_to_name = {x[wcm_genes_data_header.index("id")]: x[wcm_genes_data_header.index("symbol")] for x in wcm_genes_data[1:]}
name_to_egnumber = {x[wcm_genes_data_header.index("symbol")]: x[wcm_genes_data_header.index("id")] for x in wcm_genes_data[1:]}
wcm_genes = [egnumber_to_name.get(x) for x in wcm_genes_egnumber]

# Compute translation efficiency by scaling Mohammad et al. 2019 by this study or Li et al. 2014
def derive_translation_efficiency(ribosome_density_genes, ribosome_density, rna_genes, rna):
	derived_genes = []
	derived_translation_efficiency = []
	for gene, ribosome_density_value in zip(ribosome_density_genes, ribosome_density):
		if gene not in rna_genes:
			continue

		rna_value = rna[rna_genes.index(gene)]
		if rna_value == 0:
			continue
		derived_genes.append(gene)
		derived_translation_efficiency.append(ribosome_density_value / rna_value)
	derived_translation_efficiency = np.array(derived_translation_efficiency)
	return derived_genes, derived_translation_efficiency

derived_with_this_study_genes, derived_with_this_study_translation_efficiency = derive_translation_efficiency(
	mohammad_genes, mohammad_ribosome_density, wcm_genes, wcm_rna)
derived_with_li_genes, derived_with_li_translation_efficiency = derive_translation_efficiency(
	mohammad_genes, mohammad_ribosome_density, li_genes.tolist(), li_mrna)


# Analyze
def plot(ax, x_genes, x_data, y_genes, y_data):
	for x_gene, x_value in zip(x_genes, x_data):
		if x_gene not in y_genes:
			continue

		y_value = y_data[y_genes.index(x_gene)]

		color = COLOR_HIGHLIGHT if x_gene in rnap_genes else "tab:blue" if x_gene in ribosome_genes else "none"
		if color == "none":
			continue
		ax.scatter(
			np.log10(x_value),
			np.log10(y_value),
			s=90, c=color)
			# , edgecolors="tab:blue" if color is "none" else color)
	return


fig, ax = plt.subplots(1, 1, figsize=(5, 5))

plot(ax, li_genes, li_ribosome_density, mohammad_genes.tolist(), mohammad_ribosome_density)
range_min = np.floor(min(ax.get_xlim()[0], ax.get_ylim()[0]) * 10)/10.
range_max = np.ceil(max(ax.get_xlim()[1], ax.get_ylim()[1]) * 10)/10.
ax.set_xlim([range_min, range_max])
ax.set_ylim([range_min, range_max])
ax.plot(
		[range_min, range_max],
		[range_min, range_max], "k")
ax.set_xticks([range_min, range_max])
ax.set_yticks([range_min, range_max])
ax.set_xticklabels([])
ax.set_yticklabels([])
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
plt.savefig(output_plot.format("__clean"))

ax.set_xticklabels(ax.get_xticks())
ax.set_yticklabels(ax.get_yticks())
from matplotlib.ticker import FormatStrFormatter
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.set_xlabel("log10 Li ribosome density")
ax.set_ylabel("log10 Mohammad ribosome density")
plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75)
plt.savefig(output_plot.format(""))
plt.close("all")


# Write alternate translation efficiency
out = ['"geneId"\t"name"\t"translationEfficiency"']
for gene, value in zip(derived_with_li_genes, derived_with_li_translation_efficiency):
	if value > 10:
		continue
	egnumber = name_to_egnumber.get(gene)
	if egnumber is None:
		continue

	out.append('"{}"\t"{}"\t{}'.format(egnumber, gene, value))

with open(output_file, "w") as f:
	f.write("\n".join(out))
