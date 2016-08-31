#!/usr/bin/env python

from __future__ import division

from models.ecoli.sim.variants.gene_knockout import geneKnockout
from models.ecoli.sim.variants.gene_knockout import geneKnockoutTotalIndices

from models.ecoli.sim.variants.wildtype import wildtype
from models.ecoli.sim.variants.wildtype import wildtypeTotalIndices

from models.ecoli.sim.variants.time_step import timeStep
from models.ecoli.sim.variants.time_step import timeStepTotalIndices

from models.ecoli.sim.variants.starvation_variant import starvationVariant
from models.ecoli.sim.variants.starvation_variant import starvationVariantTotalIndices
from models.ecoli.sim.variants.kinetics_flux_coeff import kineticsFluxCoeff
from models.ecoli.sim.variants.kinetics_flux_coeff import kineticsFluxCoeffTotalIndices

from models.ecoli.sim.variants.growth_rate import growthRate
from models.ecoli.sim.variants.growth_rate import growthRateTotalIndices

from models.ecoli.sim.variants.nutrientTimeSeries import nutrientTimeSeries
from models.ecoli.sim.variants.nutrientTimeSeries import nutrientTimeSeriesTotalIndices

from models.ecoli.sim.variants.scaling_factor import scalingFactor
from models.ecoli.sim.variants.scaling_factor import scalingFactorTotalIndices

from models.ecoli.sim.variants.tf_activity import tfActivity
from models.ecoli.sim.variants.tf_activity import tfActivityTotalIndices

from models.ecoli.sim.variants.condition import condition
from models.ecoli.sim.variants.condition import conditionIndices

nameToFunctionMapping = {
	"geneKnockout": geneKnockout,
	"wildtype": wildtype,
	"timeStep": timeStep,
	"starvationVariant": starvationVariant,
	"kineticsFluxCoeff": kineticsFluxCoeff,
	"growthRate": growthRate,
	"scalingFactor": scalingFactor,
	"nutrientTimeSeries": nutrientTimeSeries,
	"tfActivity": tfActivity,
	"condition": condition,
}

nameToNumIndicesMapping = {
	"geneKnockout": geneKnockoutTotalIndices,
	"wildtype": wildtypeTotalIndices,
	"timeStep": timeStepTotalIndices,
	"starvationVariant": starvationVariantTotalIndices,
	"kineticsFluxCoeff": kineticsFluxCoeffTotalIndices,
	"growthRate": growthRateTotalIndices,
	"scalingFactor": scalingFactorTotalIndices,
	"nutrientTimeSeries": nutrientTimeSeriesTotalIndices,
	"tfActivity": tfActivityTotalIndices,
	"condition": conditionIndices,
}