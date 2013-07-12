#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

import numpy

class Fitter(object):
	""" Fitter """

	@staticmethod
	def FitSimulation(sim, kb):
		mc = sim.getState("MoleculeCounts")

		idx = {}

		rnap_ids = ["EG10893_RNA", "EG10894_RNA", "EG10895_RNA", "EG10896_RNA"]
		idx["rnaExp"] = {"rnap_70": numpy.array([i for i, x in enumerate(kb.rnas) if x["id"] in rnap_ids])}

		mc.rnaExp[idx["rnaExp"]["rnap_70"]] = numpy.mean(mc.rnaExp[idx["rnaExp"]["rnap_70"]])
		mc.rnaExp /= numpy.sum(mc.rnaExp)

		idx["monExp"] = {"rnap_70": numpy.array([i for i, x in enumerate(kb.rnas) if x["monomerId"] != None and x["id"] in rnap_ids])}
		mc.monExp[idx["monExp"]["rnap_70"]] = numpy.mean(mc.monExp[idx["monExp"]["rnap_70"]])
		mc.monExp /= numpy.sum(mc.monExp)