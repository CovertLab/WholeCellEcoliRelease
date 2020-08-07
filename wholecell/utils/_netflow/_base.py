'''
Base class for netflow interfaces with different solver backends.
All functions required for implementation with modular_fba.py are listed.
'''

from __future__ import absolute_import
from __future__ import division


class NetworkFlowProblemBase(object):
    _maximize = True
    quadratic_objective = False

    def setFlowMaterialCoeff(self, flow, material, coefficient):
        raise NotImplementedError()

    def setFlowBounds(selfs, flow, ub=None, lb=None):
        raise NotImplementedError()

    def setFlowObjectiveCoeff(self, flow, coefficient):
        raise NotImplementedError()

    def getFlowObjectiveCoeff(self, flow):
        raise NotImplementedError()

    def getFlowRates(self, flows):
        raise NotImplementedError()

    def getShadowPrices(self, materials):
        raise NotImplementedError()

    def getReducedCosts(self, materials):
        raise NotImplementedError()

    def getObjectiveValue(self):
        raise NotImplementedError()

    def getSMatrix(self):
        raise NotImplementedError()

    def getFlowNames(self):
        raise NotImplementedError()

    def getMaterialNames(self):
        raise NotImplementedError()

    def getUpperBounds(self):
        raise NotImplementedError()

    def getLowerBounds(self):
        raise NotImplementedError()

    def getObjective(self):
        raise NotImplementedError()

    def buildEqConst(self):
        raise NotImplementedError()

    def maximizeObjective(self, doMax):
        self._maximize = doMax
