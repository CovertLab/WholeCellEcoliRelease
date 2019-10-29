"""
@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Updated 7/14/16, Created 1/1/2015
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Path to files
dirname = os.path.dirname(__file__)
PATH_TO_RNASEQ_DATA = os.path.join(dirname, "data/RNAseq/seal_rpkm_without_withAAs-bis.tsv")
PATH_TO_RNASEQ_MAP = os.path.join(dirname, "data/RNAseq/mapRNAseq2EcoMAC.tsv")
PATH_TO_PROTEOME = os.path.join(dirname, "data/Proteome/Proteome_Schmidt.tsv")
PATH_TO_TRANSLATION_EFFICIENCY = os.path.join(dirname, "data/TranslationEfficiency.tsv")
PATH_TO_GENE = os.path.join(dirname, "data/EcoMAC/geneName.tsv")
PATH_TO_MRNA = os.path.join(dirname, "data/EcoMAC/length_mRNA.tsv")
PATH_TO_PROTEIN_DECAY = os.path.join(dirname, "data/proteinFastDecay.tsv")
PATH_TO_OUTPUT_TSV = os.path.join(dirname, "steadyStateOutput.tsv")
PATH_TO_OUTPUT_PDF = os.path.join(dirname, "steadyStateOutput.pdf")

# Constants
doublingTime = 45. # minutes
growthRate = np.log(2.) / (doublingTime / 60.) # doublings / hour
growthRatePerSec = growthRate / 3600. # doublings / sec
phi = 16. # amino acids / sec, Bremer & Dennis 1996
fRibActive = 0.80 # Bremer & Dennis 1996
nAvogadro = 6.e23 # molecules/mole
mrnaMassFraction = 0.041 # Neidhardt, Physiology of the Bacterial Cell, 1931, p4
nucleotideMass = 324. # representative nucleotide, g/mole
aminoAcidMass = 118.9 # representative amino acid, g/mole
shortProteinHalfLife = 2. / 60. # hours
proteinHalfLife = 10. # hours


def main():
	"""
	Performs steady state analysis, writes output to PATH_TO_OUTPUT_TSV, and
	produces plot shown in Fig. 5A.

	The rate of change of protein concentration was described as:
	dp/dt = phi / protein length * mRNA * translation efficiency * fRibActive
			* ribosomes - ((ln(2) / protein decay) - (ln(2) / doubling time))
			* protein

	where:
		phi is the ribosomal elongation rate (aa/s)
		fRibActive is the fraction of active ribosomes
	"""

	# Estimate macromolecular composition
	# (functions were fit to values reported by Bremer & Dennis, 1996, Table 2)
	cellMass = 6.e8 * growthRate**(-1.243) # cells/OD460
	rnaMass = 4.e16 * np.exp(growthRate * 0.3101) # nt/OD460
	proteinMass = 6.e17 * growthRate**(-0.19) # aa/OD460
	ribosome = 14576. * (growthRate**1.6473) # ribosome/cell

	# Get lengths of mRNA transcripts (nucleotides) and proteins (amino acids)
	mrnaLength = np.array(openfile(PATH_TO_MRNA)[:, 1], dtype = float)
	proteinLength = mrnaLength / 3

	# Compute mRNA
	mrna = computeMrna(cellMass, rnaMass, mrnaLength)

	# Compute translational efficiency
	translationEff = computeTranslationEff()

	# Assign protein decay according to N-end rule
	proteinDecay = computeProteinDecay()

	# Compute protein
	protein = computeProtein(cellMass, proteinMass, proteinLength)

	# Identify genes for which all data are available
	candidates = np.logical_and(np.logical_and(mrna > 0, protein > 0), translationEff > 0)
	proteinLength = proteinLength[candidates]
	mrna = mrna[candidates]
	translationEff = translationEff[candidates]
	protein = protein[candidates]
	proteinDecay = proteinDecay[candidates]

	# Compute production term
	x_vals = np.log10(phi / proteinLength * mrna * translationEff * fRibActive * ribosome)

	# Compute loss term
	y_vals = np.log10(protein * (proteinDecay + growthRatePerSec))

	# Save output
	with open(PATH_TO_OUTPUT_TSV, "w") as f:
		for i, x, y in zip(1 + np.where(candidates)[0], x_vals, y_vals):
			f.write("{0}\t{1}\t{2}\n".format(i, x, y))

	# Generate plot
	fig, ax = plt.subplots(1, 1, figsize=(4, 4))
	ax.scatter(x_vals, y_vals, c = "k", alpha = 0.1, s=8)
	ax_min = np.floor(min(ax.get_xlim()[0], ax.get_ylim()[0]))
	ax_max = np.ceil(max(ax.get_xlim()[1], ax.get_ylim()[1]))
	ax.plot([ax_min, ax_max], [ax_min, ax_max], "k", linewidth=1)
	ax.plot([ax_min, ax_max - 1], [ax_min + 1, ax_max], "k")
	ax.plot([ax_min + 1, ax_max], [ax_min, ax_max - 1], "k")

	# Highlight genes that were investigated further experimentally
	genesTested = ["gshA", "pnp", "carA", "cdsA", "dcuR", "bioD", "rph"]
	colors = ["b", "b", "b", "b", "r", "r", "r"]
	geneName2Id = makeGeneName2EcoMacId()
	for geneId, color in zip(genesTested, colors):
		ecomacId = geneName2Id[geneId]
		candidateIndex = np.where(1 + np.where(candidates)[0] == ecomacId)[0][0]
		ax.scatter(x_vals[candidateIndex], y_vals[candidateIndex], c=color, s=16)
		ax.text(x_vals[candidateIndex] + 0.05, y_vals[candidateIndex] - 0.05, geneId)

	# Format and save
	ax.set_xlim([ax_min, ax_max])
	ax.set_ylim([ax_min, ax_max])
	ax.set_xlabel("Log10 production rate (protein/s)")
	ax.set_ylabel("Log10 loss rate (protein/s)")
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.spines["left"].set_position(("outward", 10))
	ax.spines["bottom"].set_position(("outward", 10))
	plt.tight_layout()
	plt.savefig(PATH_TO_OUTPUT_PDF)
	plt.close("all")
	return

def computeMrna(cellMass, rnaMass, mrnaLength):
	# Load data
	rnaSeq2EcoMAC = np.array(openfile(PATH_TO_RNASEQ_MAP)[:, 0], dtype = int)
	rnaSeqData = np.array(openfile(PATH_TO_RNASEQ_DATA)[:, 1:], dtype = float)

	# Take log2 average of RNA seq data
	mask = rnaSeqData >= 1
	rnaSeqDataLog2 = np.zeros(rnaSeqData.shape)
	rnaSeqDataLog2[mask] = np.log2(rnaSeqData[mask])
	rnaSeqEcoMACFormat = rnaSeqDataLog2[rnaSeq2EcoMAC - 1, :]
	rnaLog2Averaged = 2. ** np.mean(rnaSeqEcoMACFormat, axis = 1)

	# Compute relative expression of each mRNA
	relativeExp = rnaLog2Averaged / np.sum(rnaLog2Averaged)

	# Compute counts of mRNA
	totalRnaMass = nucleotideMass / nAvogadro * rnaMass / cellMass # g/cell
	mrnaMass = mrnaLength / nAvogadro * nucleotideMass # g/mRNA
	mrnaTotalCounts = totalRnaMass * mrnaMassFraction / np.dot(relativeExp, mrnaMass) # counts/cell
	mrna = relativeExp * mrnaTotalCounts
	return mrna

def computeTranslationEff():
	data = openfile(PATH_TO_TRANSLATION_EFFICIENCY)
	geneNames = data[:, 0]
	translationEff = np.array(data[:, 1], dtype = float)
	translationEff = mapData2EcoMAC(geneNames, translationEff)

	# Compute normalized translational efficiency
	mask = translationEff > 0
	translationEffAvg = 10.**np.mean(np.log(translationEff[mask]))
	translationEff[np.logical_not(mask)] = -translationEffAvg
	translationEffNormalized = translationEff / sum(abs(translationEff))
	return translationEffNormalized

def computeProteinDecay():
	fastDecay = openfile(PATH_TO_PROTEIN_DECAY)[:, 1]
	geneName2Id = makeGeneName2EcoMacId()
	fastDecayIndices = np.array([geneName2Id[x] for x in fastDecay if x in geneName2Id])

	# Compute protein decay rates (1/s)
	fastProteinDecay = np.log(2.) / (shortProteinHalfLife * 3600.)
	proteinDecay = np.log(2.) / (proteinHalfLife * 3600.) * np.ones(len(geneName2Id))
	proteinDecay[fastDecayIndices - 1] = fastProteinDecay
	return proteinDecay

def computeProtein(cellMass, proteinMass, proteinLength):
	# Load data
	rawdata = openfile(PATH_TO_PROTEOME)
	data = np.array(rawdata[:, 1], dtype = float)
	proteome = mapData2EcoMAC(rawdata[:, 0], data)

	# Compute total counts of protein
	mask = proteome <= 0
	proteomeMinCount1 = np.copy(proteome)
	proteomeMinCount1[mask] = 1.
	relativeExp = proteomeMinCount1 / np.sum(proteomeMinCount1)
	totalProteinMass = aminoAcidMass / nAvogadro * proteinMass / cellMass # g/cell
	proteinMass = proteinLength / nAvogadro * aminoAcidMass # g/protein
	proteinTotalCounts = totalProteinMass / np.dot(relativeExp, proteinMass) # counts/cell

	# Compute counts of each protein type
	proteomeMinCount0 = np.copy(proteome)
	proteomeMinCount0[mask] = 0
	protein = proteomeMinCount0 / np.sum(proteomeMinCount0) * proteinTotalCounts
	return protein

def mapData2EcoMAC(geneNames, data):
	geneName2Id = makeGeneName2EcoMacId()
	dataMappedToEcoMAC = -1 * np.ones(len(geneName2Id))
	for geneName, value in zip(geneNames, data):
		if geneName in geneName2Id:
			dataMappedToEcoMAC[geneName2Id[geneName] - 1] = value
	return dataMappedToEcoMAC

def makeGeneName2EcoMacId():
	data = openfile(PATH_TO_GENE)
	geneName2Id = dict(zip(data[:, 2], np.array(data[:, 0], dtype = int)))
	return geneName2Id

def openfile(filename):
	with open(filename, "r") as f:
		rawdata = f.readlines()
	data = []
	for line in rawdata:
		data.append(line.split())
	return np.array(data)


if __name__ == "__main__":
	main()
