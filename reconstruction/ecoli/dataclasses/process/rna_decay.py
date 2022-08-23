"""
SimulationData for rna decay process
"""

from wholecell.utils import units

import aesara.tensor as T
from aesara import function, gradient
import numpy as np


class RnaDecay(object):
	""" RnaDecay """

	def __init__(self, raw_data, sim_data):
		self._buildRnaDecayData(raw_data, sim_data)

	def _buildRnaDecayData(self, raw_data, sim_data):
		_ = sim_data  # unused
		self.endoRNase_ids = [
			x['endoRnase'] + sim_data.getter.get_compartment_tag(x['endoRnase'])
			for x in raw_data.endoRNases]
		self.kcats = (1 / units.s) * np.array([x["kcat"].asNumber(1 / units.s) for x in raw_data.endoRNases])
		self.stats_fit = {
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
		self.sensitivity_analysis_alpha_residual = {}
		self.sensitivity_analysis_alpha_regular_i_neg = {}

		# store Km's and Residuals re-scaled (sensitivity analysis "kcat EndoRNases")
		self.sensitivity_analysis_kcat = {}
		self.sensitivity_analysis_kcat_res_ini = {}
		self.sensitivity_analysis_kcat_res_opt = {}

		# store Km's from first-order RNA decay
		self.Km_first_order_decay = []

		# store convergence of non-linear Km's (G(km))
		self.Km_convergence = []


	def km_loss_function(self, vMax, rnaConc, kDeg, isEndoRnase, alpha):
		'''
		Generates the functions used for estimating the per-RNA affinities (Michaelis-Menten
		constants) to the endoRNAses.

		The optimization problem is formulated as a multidimensional root-finding problem; the goal
		is to find a set of Michaelis-Menten constants such that the endoRNAse-mediated degradation
		under basal concentrations is consistent with the experimentally observed half-lives, thus

		(nonlinear rate) = (linear rate)

		where the nonlinear rate is the rate as predicted from some kinetic rate law, and the
		linear rate is proportional to the inverse of the observed half-life.  Then, reordering,

		0 = (nonlinear rate) - (linear rate)

		is (for the moment) the root we wish to find, for each RNA species, giving us the
		multidimensional function

		R_aux = (nonlinear rate) - (linear rate)

		This is the unnormalized residual function; the normalized residuals are

		R = (nonlinear rate)/(linear rate) - 1

		In addition to matching our half-lives we also desire the Michaelis-Menten constants to be
		non-negative (negative values have no physical meaning).  Thus we introduce a penalty term
		for negative values.  TODO (John): explain penalty term

		The two terms (the residuals R and the negative value penalty Rneg) are combined into one
		'loss' function L (alpha is the weighting on the negative value penalty):

		L = ln((exp(R) + exp(alpha*Rneg))/2)
		  = ln(exp(R) + exp(alpha*Rneg)) - ln(2)

		The loss function has one element for each RNA.  This functional form is a soft
		(continuous and differentiable) approximation to

		L = max(R, alpha*Rneg)

		The root finder, provided with L, will attempt to make each element of L as close to zero
		as possible, and therefore minimize both R and Rneg.

		The third-party package Aesara (formerly Theano) creates the functions and finds an analytic
		expression for the Jacobian.

		Parameters
		----------
		vMax: scalar
			The total endoRNAse capacity, in dimensions of amount per volume per time.
		rnaConc: 1-D array, float
			Concentrations of RNAs (that will be degraded), in dimensions of amount per volume.
		kDeg: 1-D array, float
			Experimentally observed degradation rates (computed from half-lives), in dimensions of
			per unit time.
		isEndoRnase: 1-D array, bool
			A vector that is True everywhere that an RNA corresponds to an endoRNAse; that is, an
			endoRNAse (or endoRNAse subunit) mRNA.
		alpha: scalar, >0
			Regularization weight, used to penalize for negative Michaelis-Menten value predictions
			during the course of the optimization.  Typical value is 0.5.

		Returns
		-------
		L: function
			The 'loss' function.
		Rneg: function
			The negative Michaelis-Menten constant penalty terms.
		R: function
			The residual error (deviation from steady-state).
		Lp: function
			The Jacobian of the loss function L with respect to the Michaelis-Menten constants.
		R_aux: function
			Unnormalized 'residual' function.
		L_aux: function
			Unnormalized 'loss' function.
		Lp_aux: function
			Jacobian of the unnormalized 'loss' function.
		Jacob: function
			Duplicate with Lp.
		Jacob_aux: function
			Duplicate with Lp_aux.

		Notes
		-----
		The regularization term also includes a penalty for the endoRNAse residuals, as well as a
		fixed weighting (WFendoR = 0.1).
		TODO (John): Why is this needed?  It seems redundant.
		TODO (John): How do we know this weight is sufficient?

		All of the outputs are Aesara functions, and take a 1-D array of Michaelis-Menten constants
		as their sole inputs.  All of the functions return a 1-D array, with the exception of the
		Jacobians, which return matrices.

		TODO (John): Remove the redundant outputs.

		TODO (John): Consider redesigning this as an objective minimization problem rather than a
		root finding problem.

		TODO (John): Consider replacing the Michaelis-Menten constants with logarithmic
		equivalents, thereby eliminating the requirement for the negative value penalty.

		TODO (John): Consider moving this method out of this class, as it is, in fact, a static
		method, and isn't utilized anywhere within this class.
		'''

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

		J = gradient.jacobian(LossFunction, km)
		J_aux = gradient.jacobian(LossFunction_aux, km)
		Jacob = function([km], J)
		Jacob_aux = function([km], J_aux)
		L = function([km], LossFunction)
		L_aux = function([km], LossFunction_aux)
		Rneg = function([km], regularizationNegativeNumbers)
		R = function([km], residual)
		Lp = function([km], J)
		Lp_aux = function([km], J_aux)
		R_aux = function([km], residual_aux)

		return L, Rneg, R, Lp, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux
