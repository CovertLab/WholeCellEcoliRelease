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

from models.ecoli.sim.variants.metabolism_target_range import metabolismTargetRange
from models.ecoli.sim.variants.metabolism_target_range import metabolismTargetRangeTotalIndices

from models.ecoli.sim.variants.metabolism_objective_kinetic_homeostatic_ratio import metabolismKineticHomeostaticRatio
from models.ecoli.sim.variants.metabolism_objective_kinetic_homeostatic_ratio import metabolismKineticHomeostaticRatioTotalIndices

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

from models.ecoli.sim.variants.anaerobic import anaerobic
from models.ecoli.sim.variants.anaerobic import anaerobicIndices

nameToFunctionMapping = {
	"geneKnockout": geneKnockout,
	"wildtype": wildtype,
	"timeStep": timeStep,
	"starvationVariant": starvationVariant,
	"metabolismTargetRange":metabolismTargetRange,
	"growthRate": growthRate,
	"scalingFactor": scalingFactor,
	"nutrientTimeSeries": nutrientTimeSeries,
	"tfActivity": tfActivity,
	"condition": condition,
	"anaerobic": anaerobic,
	"metabolismKineticHomeostaticRatio": metabolismKineticHomeostaticRatio,
}

nameToNumIndicesMapping = {
	"geneKnockout": geneKnockoutTotalIndices,
	"wildtype": wildtypeTotalIndices,
	"timeStep": timeStepTotalIndices,
	"starvationVariant": starvationVariantTotalIndices,
	"metabolismTargetRangeTotalIndices":metabolismTargetRangeTotalIndices,
	"growthRate": growthRateTotalIndices,
	"scalingFactor": scalingFactorTotalIndices,
	"nutrientTimeSeries": nutrientTimeSeriesTotalIndices,
	"tfActivity": tfActivityTotalIndices,
	"condition": conditionIndices,
	"anaerobic": anaerobicIndices,
	"metabolismKineticHomeostaticRatio": metabolismKineticHomeostaticRatioTotalIndices,
}