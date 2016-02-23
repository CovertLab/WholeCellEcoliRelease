#!/usr/bin/env python

from __future__ import division

from models.ecoli.sim.variants.gene_knockout import geneKnockout
from models.ecoli.sim.variants.gene_knockout import geneKnockoutTotalIndices

from models.ecoli.sim.variants.wildtype import wildtype
from models.ecoli.sim.variants.wildtype import wildtypeTotalIndices

from models.ecoli.sim.variants.time_step import timeStep
from models.ecoli.sim.variants.time_step import timeStepTotalIndices

from models.ecoli.sim.variants.kinetics_flux_coeff import kineticsFluxCoeff
from models.ecoli.sim.variants.kinetics_flux_coeff import kineticsFluxCoeffTotalIndices

from models.ecoli.sim.variants.growth_rate import growthRate
from models.ecoli.sim.variants.growth_rate import growthRateTotalIndices

from models.ecoli.sim.variants.environment import environment
from models.ecoli.sim.variants.environment import environmentTotalIndices


nameToFunctionMapping = {
	"geneKnockout": geneKnockout,
	"wildtype": wildtype,
	"timeStep": timeStep,
	"kineticsFluxCoeff": kineticsFluxCoeff,
	"growthRate": growthRate,
	"environment": environment,
}

nameToNumIndicesMapping = {
	"geneKnockout": geneKnockoutTotalIndices,
	"wildtype": wildtypeTotalIndices,
	"timeStep": timeStepTotalIndices,
	"kineticsFluxCoeff": kineticsFluxCoeffTotalIndices,
	"growthRate": growthRateTotalIndices,
	"environment": environmentTotalIndices,
}