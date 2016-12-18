"""
SimulationData for rna decay process

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 08/18/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np
import theano.tensor as T
import theano

class RnaDecay(object):
	""" RnaDecay """

	def __init__(self, raw_data, sim_data):
		self._buildRnaDecayData(raw_data, sim_data)

	def _buildRnaDecayData(self, raw_data, sim_data):
		self.mrna_index = 2
		self.rrna_index = 3
		self.trna_index = 4
		self.rtrna_index = 5

		self.endoRnaseIds = [x["endoRnase"].encode("utf-8") for x in raw_data.endoRnases]
		self.kcats = (1 / units.s) * np.array([x["kcat"].asNumber(1 / units.s) for x in raw_data.endoRnases])
		self.StatsFit = {
				'LossKm': 0.,
				'LossKmOpt': 0.,
				'RnegKmOpt': 0.,
				'ResKm': 0.,
				'ResKmOpt': 0.,
				'ResEndoRNKm': 0.,
				'ResEndoRNKmOpt': 0.,
				'ResScaledKm': 0.,
				'ResScaledKmOpt': 0.,
			}

		# store Residuals re-scaled (sensitivity analysis "alpha")
		self.SensitivityAnalysisAlphaResidual = {}
		self.SensitivityAnalysisAlphaRegulariNeg = {}

		# store Km's and Residuals re-scaled (sensitivity analysis "kcat EndoRNases")
		self.SensitivityAnalysisKcat = {}
		self.SensitivityAnalysisKcat_ResIni = {}
		self.SensitivityAnalysisKcat_ResOpt = {}

		# store Km's from first-order RNA decay
		self.KmFirstOrderDecay = []

		# store convergence of non-linear Km's (G(km))
		self.KmConvergence = []


		self.TargetEndoRNasesFullMRNA = np.zeros(len(self.endoRnaseIds))
		self.TargetEndoRNasesFullTRNA = np.zeros(len(self.endoRnaseIds))
		self.TargetEndoRNasesFullRRNA = np.zeros(len(self.endoRnaseIds))

		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10856-MONOMER[p]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10857-MONOMER[c]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("G7175-MONOMER[c]")] = 0
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10859-MONOMER[i]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG11299-MONOMER[c]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10860-MONOMER[c]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10861-MONOMER[c]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("G7365-MONOMER[c]")] = self.mrna_index
		self.TargetEndoRNasesFullMRNA[self.endoRnaseIds.index("EG10862-MONOMER[c]")] = self.mrna_index

		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10856-MONOMER[p]")] = self.trna_index
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10857-MONOMER[c]")] = 0
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("G7175-MONOMER[c]")] = 1
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10859-MONOMER[i]")] = self.trna_index
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG11299-MONOMER[c]")] = 0
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10860-MONOMER[c]")] = self.trna_index
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10861-MONOMER[c]")] = self.trna_index
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("G7365-MONOMER[c]")] = self.trna_index
		self.TargetEndoRNasesFullTRNA[self.endoRnaseIds.index("EG10862-MONOMER[c]")] = self.trna_index

		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10856-MONOMER[p]")] = self.rrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10857-MONOMER[c]")] = self.rtrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("G7175-MONOMER[c]")] = 0
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10859-MONOMER[i]")] = self.rrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG11299-MONOMER[c]")] = self.rtrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10860-MONOMER[c]")] = self.rrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10861-MONOMER[c]")] = self.rrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("G7365-MONOMER[c]")] = self.rrna_index
		self.TargetEndoRNasesFullRRNA[self.endoRnaseIds.index("EG10862-MONOMER[c]")] = self.rrna_index

	def kmLossFunction(self, vMax, rnaConc, kDeg, isEndoRnase, alpha):
		N = rnaConc.size
		km = T.dvector()

		# Residuals of non-linear optimization
		residual = (vMax / km / kDeg) / (1 + (rnaConc / km).sum()) - np.ones(N)
		residual_aux = (vMax * rnaConc / km) / (1 + (rnaConc / km).sum()) - (kDeg * rnaConc)

		# Counting negative Km's (first regularization term) 
		regularizationNegativeNumbers = (np.ones(N) - km / np.abs(km)).sum() / N

		# Penalties for EndoR Km's, which might be potentially nonf-fitted
		regularizationEndoR = (isEndoRnase * np.abs(residual)).sum()
		
		# Multi objective-based regularization
		WFendoR = 0.1 # weighting factor to protect Km optimized of EndoRNases
		regularization = regularizationNegativeNumbers + (WFendoR * regularizationEndoR)

		# Loss function
		LossFunction = T.log(T.exp(residual) + T.exp(alpha * regularization)) - T.log(2)
		LossFunction_aux = T.log(T.exp(residual_aux) + T.exp(alpha * regularization)) - T.log(2)

		J = theano.gradient.jacobian(LossFunction, km)
		J_aux = theano.gradient.jacobian(LossFunction_aux, km)
		Jacob = theano.function([km], J)
		Jacob_aux = theano.function([km], J_aux)
		L = theano.function([km], LossFunction)
		L_aux = theano.function([km], LossFunction_aux)
		Rneg = theano.function([km], regularizationNegativeNumbers)
		R = theano.function([km], residual)
		Lp = theano.function([km], J)
		Lp_aux = theano.function([km], J_aux)
		R_aux = theano.function([km], residual_aux)

		return L, Rneg, R, Lp, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux
