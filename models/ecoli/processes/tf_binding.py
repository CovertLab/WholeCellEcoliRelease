#!/usr/bin/env python

"""
TfBinding

Bind transcription factors to DNA

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/14/16
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound


class TfBinding(wholecell.processes.process.Process):
    """ TfBinding """

    _name = "TfBinding"

    # Constructor
    def __init__(self):
        super(TfBinding, self).__init__()

    # Construct object graph
    def initialize(self, sim, sim_data):
        super(TfBinding, self).initialize(sim, sim_data)

        # Get IDs of free DNA targets (alpha), and transcription factor bound targets (tfNames)
        recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
        self.tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames if x.split("__")[-1] != "alpha"]))
        alphaNames = [x for x in recruitmentColNames if x.endswith("__alpha")]
        tfNames = {}
        for tf in self.tfs:
            tfNames[tf] = [x for x in recruitmentColNames if x.endswith("__" + tf)]

        # Get constants
        self.nAvogadro = sim_data.constants.nAvogadro
        self.cellDensity = sim_data.constants.cellDensity

        # Create dictionaries and method
        self.tfNTargets = sim_data.process.transcription_regulation.tfNTargets
        self.pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF
        self.tfToTfType = sim_data.process.transcription_regulation.tfToTfType

        # Build views
        self.alphaView = self.bulkMoleculesView(alphaNames)
        self.tfDnaBoundViews = {}
        self.tfMoleculeActiveView = {}
        self.tfMoleculeInactiveView = {}
        for tf in self.tfs:
            self.tfDnaBoundViews[tf] = self.bulkMoleculesView(tfNames[tf])
            self.tfMoleculeActiveView[tf] = self.bulkMoleculeView(tf + "[c]")
            if self.tfToTfType[tf] == "1CS":
                if tf == sim_data.process.transcription_regulation.activeToBound[tf]:
                    self.tfMoleculeInactiveView[tf] = self.bulkMoleculeView(sim_data.process.equilibrium.getUnbound(tf + "[c]"))
                else:
                    self.tfMoleculeInactiveView[tf] = self.bulkMoleculeView(sim_data.process.transcription_regulation.activeToBound[tf] + "[c]")
            elif self.tfToTfType[tf] == "2CS":
                self.tfMoleculeInactiveView[tf] = self.bulkMoleculeView(sim_data.process.two_component_system.activeToInactiveTF[tf + "[c]"])


    def calculateRequest(self):
        # Request all counts of free DNA targets
        self.alphaView.requestAll()

        # Request all counts of active transcription factors
        for view in self.tfMoleculeActiveView.itervalues():
            view.requestAll()

        # Request all counts of DNA bound transcription factors
        for view in self.tfDnaBoundViews.itervalues():
            view.requestAll()


    def evolveState(self):
        # Set counts of all free DNA targets to 1
        self.alphaView.countsIs(1)

        # Create vectors for storing values
        nTfs = len(self.tfs)
        pPromotersBound = np.zeros(nTfs, np.float64)
        nPromotersBound = np.zeros(nTfs, np.float64)
        nActualBound = np.zeros(nTfs, np.float64)

        for i, tf in enumerate(self.tfs):
            # Get counts of transcription factors
            tfActiveCounts = self.tfMoleculeActiveView[tf].count()
            tfInactiveCounts = None
            if self.tfToTfType[tf] != "0CS":
                tfInactiveCounts = self.tfMoleculeInactiveView[tf].total()
            tfBoundCounts = self.tfDnaBoundViews[tf].counts()
            promoterCounts = self.tfNTargets[tf]

            # Free all DNA-bound transcription factors into free active transcription factors
            self.tfDnaBoundViews[tf].countsIs(0)
            self.tfMoleculeActiveView[tf].countInc(tfBoundCounts.sum())

            # If there are no active transcription factors, continue to the next transcription factor
            if tfActiveCounts == 0:
                continue

            # Compute probability of binding the promoter
            pPromoterBound = None
            if self.tfToTfType[tf] == "0CS":
                if tfActiveCounts + tfBoundCounts.sum() > 0:
                    pPromoterBound = 1.
                else:
                    pPromoterBound = 0.
            else:
                pPromoterBound = self.pPromoterBoundTF(tfActiveCounts, tfInactiveCounts)

            # Determine the number of available promoter sites to bind
            nToBind = int(stochasticRound(self.randomState, promoterCounts * pPromoterBound))

            # If there are no promoter sites to bind, continue to the next transcription factor
            if nToBind == 0:
                continue

            # Determine randomly which DNA targets to bind based on which of the following is most limiting:
            # number of promoter sites to bind, number of DNA targets, or number of active transcription factors
            boundLocs = np.zeros_like(tfBoundCounts)
            boundLocs[
                self.randomState.choice(tfBoundCounts.size, size = np.min((nToBind, tfBoundCounts.size, self.tfMoleculeActiveView[tf].count())), replace = False)
                ] = 1

            # Update counts of free and DNA-bound transcription factors
            self.tfMoleculeActiveView[tf].countDec(boundLocs.sum())
            self.tfDnaBoundViews[tf].countsIs(boundLocs)

            # Record values
            pPromotersBound[i] = pPromoterBound
            nPromotersBound[i] = pPromoterBound * self.tfNTargets[tf]
            nActualBound[i] = boundLocs.sum()

        # Write values to listeners
        self.writeToListener("RnaSynthProb", "pPromoterBound", pPromotersBound)
        self.writeToListener("RnaSynthProb", "nPromoterBound", nPromotersBound)
        self.writeToListener("RnaSynthProb", "nActualBound", nActualBound)
