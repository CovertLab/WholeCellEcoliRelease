#!/usr/bin/env python

"""
Hard coded data file
"""
import collections

AMINO_ACID_1_TO_3_ORDERED = collections.OrderedDict(( # TOKB
	("A", "ALA-L[c]"), ("R", "ARG-L[c]"), ("N", "ASN-L[c]"), ("D", "ASP-L[c]"),
	("C", "CYS-L[c]"), ("E", "GLU-L[c]"), ("Q", "GLN-L[c]"), ("G", "GLY[c]"),
	("H", "HIS-L[c]"), ("I", "ILE-L[c]"), ("L", "LEU-L[c]"), ("K", "LYS-L[c]"),
	("M", "MET-L[c]"), ("F", "PHE-L[c]"), ("P", "PRO-L[c]"), ("S", "SER-L[c]"),
	("T", "THR-L[c]"), ("W", "TRP-L[c]"), ("Y", "TYR-L[c]"), ("U", "SEC-L[c]"),
	("V", "VAL-L[c]")
	))

MOLECULAR_WEIGHT_KEYS = [
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
	'RNA' # nonspecific RNA
	]

MOLECULAR_WEIGHT_ORDER = {
	key:index for index, key in enumerate(MOLECULAR_WEIGHT_KEYS)
	}

COMPLEXES_REQUIRE_MODIFIED = ['ACETYL-COA-CARBOXYLMULTI-CPLX', 'BCCP-CPLX',
	'CPLX0-263', 'CPLX0-2901','CPLX0-7721', 'CPLX0-7748', 'CPLX0-7754',
	'CPLX0-7795', 'CPLX0-7884','CPLX0-7885', 'ENTMULTI-CPLX', 'GCVMULTI-CPLX', 
	'PHOSPHASERDECARB-CPLX', 'PHOSPHASERDECARB-DIMER', 'PHOSPHO-OMPR', 
	'PROTEIN-NRIP', 'SAMDECARB-CPLX']

COMPLEXES_NOT_FORMED = [
	# RNA poly + sigma factor
	"RNAPE-CPLX", "CPLX0-221", "CPLX0-222", "RNAPS-CPLX", "RNAP32-CPLX",
	"RNAP54-CPLX", "RNAP70-CPLX",
	]

REACTION_ENZYME_ASSOCIATIONS = {
	# problem: multiple associated enzymes
	## PTS system
	"FEIST_ACMANAptspp":None,
	"FEIST_ACMUMptspp":None,
	"FEIST_ASCBptspp":None,
	"FEIST_DHAPT":None,
	"FEIST_FRUpts2pp":None,
	"FEIST_FRUptspp":None,
	"FEIST_GALTptspp":None,
	"FEIST_GAMptspp":None,
	"FEIST_GTHRDHpp":None,
	"FEIST_MALTptspp":None,
	"FEIST_MANGLYCptspp":None,
	"FEIST_MANptspp":None,
	"FEIST_MNLptspp":None,
	"FEIST_SBTptspp":None,
	"FEIST_TREptspp":None,

	## proenzymes
	"FEIST_ADMDC":["SPED-MONOMER"],
	"FEIST_PSD120":["PSD-MONOMER"],
	"FEIST_PSD140":["PSD-MONOMER"],
	"FEIST_PSD141":["PSD-MONOMER"],
	"FEIST_PSD160":["PSD-MONOMER"],
	"FEIST_PSD161":["PSD-MONOMER"],
	"FEIST_PSD180":["PSD-MONOMER"],
	"FEIST_PSD181":["PSD-MONOMER"],
	"FEIST_ASP1DC":["ASPDECARBOX-MONOMER"],

	## acyl carrier protein
	"FEIST_CITL":None, # ['ACECITLY-CPLX', 'G6340-MONOMER']

	## multiple enzyme associations, not complexed
	"FEIST_NO3R2pp":None, # ['NARW-MONOMER', 'NITRATREDUCTZ-CPLX']
	"FEIST_O16AP1pp":None, # ['EG11982-MONOMER', 'G7090-MONOMER']
	"FEIST_O16AP2pp":None, # ['EG11982-MONOMER', 'G7090-MONOMER']
	"FEIST_O16AP3pp":None, # ['EG11982-MONOMER', 'G7090-MONOMER']
	"FEIST_PDX5PS":None, # ['CPLX0-7847', 'CPLX0-321', 'CPLX0-743']
	"FEIST_RIBabcpp":None, # ['ABC-28-CPLX', 'CPLX0-7646']
	"FEIST_THZPSN":None, # ['CPLX0-248', 'THIF-MONOMER', 'CPLX-8029', 'THII-MONOMER']

	}

METABOLITE_CONCENTRATIONS = { # mol / L # TODO: move to SQL
	"glu-L": 9.60e-2,
	"gthrd": 1.70e-2,
	"fdp": 1.50e-2,
	"atp": 9.60e-3,
	"u3aga": 9.20e-3,
	"utp": 8.30e-3,
	"gtp": 4.90e-3,
	"dttp": 4.60e-3,
	"asp-L": 4.20e-3,
	"val-L": 4.00e-3,
	"6pgc": 3.80e-3,
	"gln-L": 3.80e-3,
	"ctp": 2.70e-3,
	"ala-L": 2.60e-3,
	"nad": 2.60e-3,
	"udpg": 2.50e-3,
	"uri": 2.10e-3,
	"cit": 2.00e-3,
	"udp": 1.80e-3,
	"mal-L": 1.70e-3,
	"3pg": 1.50e-3,
	"citr-L": 1.40e-3,
	"coa": 1.40e-3,
	"glyc-R": 1.40e-3,
	"gam6p": 1.20e-3,
	"actp": 1.10e-3,
	"6pgl": 1.00e-3,
	"gdp": 6.80e-4,
	"accoa": 6.10e-4,
	"cbasp": 5.90e-4,
	"arg-L": 5.70e-4,
	"succ": 5.70e-4,
	"udpglcur": 5.70e-4,
	"adp": 5.60e-4,
	"asn-L": 5.10e-4,
	"akg": 4.40e-4,
	"lys-L": 4.10e-4,
	"pro-L": 3.90e-4,
	"dtdp": 3.80e-4,
	"dhap": 3.70e-4,
	"hcys-L": 3.70e-4,
	"cmp": 3.60e-4,
	"amp": 2.80e-4,
	"succoa": 2.30e-4,
	"gua": 1.90e-4,
	"pep": 1.80e-4,
	"amet": 1.80e-4,
	"thr-L": 1.80e-4,
	"fad": 1.70e-4,
	"met-L": 1.50e-4,
	"23dhb": 1.40e-4,
	"fum": 1.20e-4,
	"nadph": 1.20e-4,
	"phpyr": 9.00e-5,
	"nadh": 8.30e-5,
	"acgam1p": 8.20e-5,
	"his-L": 6.80e-5,
	"ser-L": 6.80e-5,
	"4hbz": 5.20e-5,
	"dgmp": 5.10e-5,
	"glyc3p": 4.90e-5,
	"acorn": 4.30e-5,
	"glcn": 4.20e-5,
	# "23camp": 3.50e-5, # can't be formed by the reaction network
	"dctp": 3.50e-5,
	"malcoa": 3.50e-5,
	"tyr-L": 2.90e-5,
	"gmp": 2.40e-5,
	"aacoa": 2.20e-5,
	"ribflv": 1.90e-5,
	"phe-L": 1.80e-5,
	"acon-C": 1.60e-5,
	"datp": 1.60e-5,
	"csn": 1.40e-5,
	"skm": 1.40e-5,
	"histd": 1.30e-5,
	"dhor-S": 1.20e-5,
	"quln": 1.20e-5,
	"trp-L": 1.20e-5,
	"orn": 1.00e-5,
	"damp": 8.80e-6,
	"aps": 6.60e-6,
	# "inost": 5.70e-6, # can't be formed by the reaction network
	"ppcoa": 5.30e-6,
	"adpglc": 4.30e-6,
	"anth": 3.50e-6,
	"dad-2": 2.80e-6,
	"cytd": 2.60e-6,
	"nadp": 2.10e-6,
	"gsn": 1.60e-6,
	"ade": 1.50e-6,
	"dgsn": 5.20e-7,
	"adn": 1.30e-7,
	}

POLYMERIZED_AMINO_ACID_WEIGHTS = [
	{"base molecule":"ALA-L",	"frame id":"Polymerized ALA-L",	"mw":71.0785},
	{"base molecule":"ARG-L",	"frame id":"Polymerized ARG-L",	"mw":157.1957},
	{"base molecule":"ASN-L",	"frame id":"Polymerized ASN-L",	"mw":114.1034},
	{"base molecule":"ASP-L",	"frame id":"Polymerized ASP-L",	"mw":114.0796},
	{"base molecule":"CYS-L",	"frame id":"Polymerized CYS-L",	"mw":103.1385},
	{"base molecule":"GLU-L",	"frame id":"Polymerized GLU-L",	"mw":128.1064},
	{"base molecule":"GLN-L",	"frame id":"Polymerized GLN-L",	"mw":128.1302},
	{"base molecule":"GLY",	"frame id":"Polymerized GLY",	"mw":57.0517},
	{"base molecule":"HIS-L",	"frame id":"Polymerized HIS-L",	"mw":137.1413},
	{"base molecule":"ILE-L",	"frame id":"Polymerized ILE-L",	"mw":113.1589},
	{"base molecule":"LEU-L",	"frame id":"Polymerized LEU-L",	"mw":113.1589},
	{"base molecule":"LYS-L",	"frame id":"Polymerized LYS-L",	"mw":129.1817},
	{"base molecule":"MET-L",	"frame id":"Polymerized MET-L",	"mw":131.1921},
	{"base molecule":"PHE-L",	"frame id":"Polymerized PHE-L",	"mw":147.1761},
	{"base molecule":"PRO-L",	"frame id":"Polymerized PRO-L",	"mw":97.1163},
	{"base molecule":"SER-L",	"frame id":"Polymerized SER-L",	"mw":87.0775},
	{"base molecule":"THR-L",	"frame id":"Polymerized THR-L",	"mw":101.1043},
	{"base molecule":"TRP-L",	"frame id":"Polymerized TRP-L",	"mw":186.213},
	{"base molecule":"TYR-L",	"frame id":"Polymerized TYR-L",	"mw":163.1751},
	{"base molecule":"SEC-L",	"frame id":"Polymerized SEC-L",	"mw":149.0252},
	{"base molecule":"VAL-L",	"frame id":"Polymerized VAL-L",	"mw":99.1321},
	]

POLYPEPTIDE_END_WEIGHT = {"base molecule":"H2O",	"frame id":"Polypeptide terminal hydroxyl",	"mw":18.0148}

POLYMERIZED_NUCLEOTIDE_WEIGHTS = [
	{"base molecule":"ATP",	"frame id":"Polymerized ADN",	"mw":328.1999},
	{"base molecule":"CTP",	"frame id":"Polymerized CYTD",	"mw":304.1739},
	{"base molecule":"GTP",	"frame id":"Polymerized GSN",	"mw":344.1989},
	{"base molecule":"UTP",	"frame id":"Polymerized URI",	"mw":305.158},
	]

POLYMERIZED_DEOXY_NUCLEOTIDE_WEIGHTS = [
	{"base molecule":"DATP",	"frame id":"Polymerized DAD-2",	"mw":312.2009},
	{"base molecule":"DCTP",	"frame id":"Polymerized DCYT",	"mw":288.1749},
	{"base molecule":"DGTP",	"frame id":"Polymerized DGSN",	"mw":328.1999},
	{"base molecule":"DTTP",	"frame id":"Polymerized THYMD",	"mw":303.1858},
	]

RNA_END_WEIGHT = {"base molecule":"PPI",	"frame id":"Nucleic acid terminal pyrophosphate",	"mw":174.9489}