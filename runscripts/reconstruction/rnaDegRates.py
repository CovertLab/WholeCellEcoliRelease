#!/usr/bin/env python

# Compares model RNA degradation rates to rate data from Moffitt et al. 2016
# Requires 2 files in validation/ecoli/flat: geneIDs.tsv and moffitt2016_mrna_deg_rates.tsv
# Outputs rnaDegRates.pdf plot to directory that script is run from

from __future__ import absolute_import, division, print_function

import io
import os
from typing import Any, Dict, Union

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.io import tsv

GENE_IDS = os.path.join("validation","ecoli","flat","geneIDs.tsv")
DEG_RATES = os.path.join("validation","ecoli","flat","moffitt2016_mrna_deg_rates.tsv")

geneDict = {}
with io.open(GENE_IDS, "rb") as csvfile:
	reader = tsv.dict_reader(csvfile)
	for row in reader:
		genes = row["Names"].replace("\"","").replace("(","").replace(")","").split(" ")
		geneDict[row["FrameID"]] = genes

rateDict = {}
with io.open(DEG_RATES, "rb") as csvfile:
	reader = tsv.dict_reader(csvfile)
	for row in reader:
		if row["Sample"] == "WT -kas replicate 1" and row["Rate"] != "":
			rateDict[row["Name"]] = float(row["Rate"])
		elif row["Sample"] == "WT -kas replicate 2" and row["Rate"] != "":
			if row["Name"] in rateDict:
				rateDict[row["Name"]] = (rateDict[row["Name"]] + float(row["Rate"])) / 2
			else:
				rateDict[row["Name"]] = float(row["Rate"])
		elif row["Sample"] == "WT +kas replicate 1" and row["Name"] not in rateDict:
			rateDict[row["Name"]] = -1

raw_data = KnowledgeBaseEcoli()  # type: Any

modelRates = {}
paperRates = {}  # type: Dict[str, Union[int, float]]

for rna in raw_data.rnas:
	geneID = rna["geneId"]
	modelRates[geneID] = 60. / rna["halfLife"]
	paperRates[geneID] = 0
	if geneID in geneDict:
		for gene in rateDict:
			if gene in geneDict[geneID]:
				paperRates[geneID] = rateDict[gene]
				break

model = []
paper = []
for key in modelRates.keys():
	if paperRates[key] > 0:
		model.append(modelRates[key])
		paper.append(paperRates[key])

plt.figure(figsize = (8, 8))
maxLine = 1.1 * max(max(paper), max(model))

plt.plot([0, maxLine], [0, maxLine], '--r')	
plt.plot(model, paper, 'o', markeredgecolor = 'k', markerfacecolor = 'none')
plt.axis([0, 1, 0, maxLine])
Correlation_ExpPred = np.corrcoef(model, paper)[0][1]

plt.xlabel("RNA decay rate expected from model [1/min]")
plt.ylabel("RNA decay rate from paper (Moffitt et al. 2016) [1/min]")

plt.savefig("rnaDegRates.pdf")

# print("no match in data:")
# count = 0
# for gene, rate in paperRates.items():
# 	if rate == 0:
# 		print(gene)
# 		count += 1
