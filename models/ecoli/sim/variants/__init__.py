#!/usr/bin/env python

from __future__ import division

from models.ecoli.sim.variants.gene_knockout import geneKnockout
from models.ecoli.sim.variants.gene_knockout import geneKnockoutTotalIndices

from models.ecoli.sim.variants.wildtype import wildtype
from models.ecoli.sim.variants.wildtype import wildtypeTotalIndices


nameToFunctionMapping = {
	"geneKnockout": geneKnockout,
	"wildtype": wildtype,
}

nameToNumIndicesMapping = {
	"geneKnockout": geneKnockoutTotalIndices,
	"wildtype": wildtypeTotalIndices,
}