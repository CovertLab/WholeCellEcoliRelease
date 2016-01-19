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

	def km1(self, km, vMax, rnaConc, kDeg):
		print "hello"
		return  vMax * rnaConc / km / (1 + (rnaConc / km)) - kDeg * rnaConc

	def km2(self, km, vMax, rnaConc, kDeg):
		return  vMax * rnaConc / km / (1 + (rnaConc / km).sum()) - kDeg * rnaConc

	def kmFunctions(self, vMax, rnaConc, kDeg):
		assert rnaConc.size == kDeg.size
		N = rnaConc.size
		km = T.dvector()
		denominator = 1 + (rnaConc / km).sum()
		residual = []
		print "start"
		residual = (vMax * rnaConc / km) / denominator - kDeg * rnaConc
		residual2 = (vMax / km) / denominator - kDeg
		# for i in xrange(N):
		# 	numerator = vMax * rnaConc[i] / km[i]
		# 	denominator = 1 + rnaConc[i] / km[i]
		# 	residual.append(
		# 		numerator / denominator - (kDeg[i] * rnaConc[i])
		# 		)
		J = theano.gradient.jacobian(residual, km)
		J2 = theano.gradient.jacobian(residual2, km)
		print "end"
		# f = theano.function([km], T.stack(*residual))
		f = theano.function([km], residual)
		f2 = theano.function([km], residual2)

		print "done"
		fp = theano.function([km], J)
		fp2 = theano.function([km], J2)

		return f, fp, f2, fp2
