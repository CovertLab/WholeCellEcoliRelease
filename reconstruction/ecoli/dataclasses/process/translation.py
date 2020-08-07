"""
SimulationData for translation process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/09/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates, make_elongation_rates_flat

class Translation(object):
    """ Translation """

    def __init__(self, raw_data, sim_data, options):
        self._buildMonomerData(raw_data, sim_data, options['disable_measured_protein_deg'])
        self._buildTranslation(raw_data, sim_data)
        self._buildTranslationEfficiency(raw_data, sim_data, options['alternate_translation_efficiency'])

        self.ribosomal_proteins = self._build_ribosomal_proteins(raw_data, sim_data)
        self.elongation_rates = self._build_elongation_rates(raw_data, sim_data)

    def _buildMonomerData(self, raw_data, sim_data, disable_measured_protein_deg):
        assert all([len(protein['location']) == 1 for protein in raw_data.proteins])
        ids = ['{}[{}]'.format(protein['id'], protein['location'][0]) for protein in raw_data.proteins]

        rnaIds = []

        for protein in raw_data.proteins:
            rnaId = protein['rnaId']

            rnaLocation = None
            for rna in raw_data.rnas:
                if rna['id'] == rnaId:
                    assert len(rna['location']) == 1
                    rnaLocation = rna['location'][0]
                    break

            rnaIds.append('{}[{}]'.format(
                rnaId,
                rnaLocation
                ))

        lengths = []
        aaCounts = []
        sequences = []

        for protein in raw_data.proteins:
            sequence = protein['seq']

            counts = []

            for aa in sim_data.amino_acid_1_to_3_ordered.viewkeys():
                counts.append(
                    sequence.count(aa)
                    )

            lengths.append(len(sequence))
            aaCounts.append(counts)
            sequences.append(sequence)

        maxSequenceLength = max(len(seq) for seq in sequences)

        mws = np.array([protein['mw'] for protein in raw_data.proteins]).sum(axis = 1)

        size = len(rnaIds)

        nAAs = len(aaCounts[0])

        # Calculate degradation rates based on N-rule
        # TODO: citation
        deg_rate_units = 1 / units.s
        fastRate = (np.log(2) / (2*units.min)).asNumber(deg_rate_units)
        slowRate = (np.log(2) / (10*60*units.min)).asNumber(deg_rate_units)
        self.fastRate = fastRate

        fastAAs = ["R", "K", "F", "L", "W", "Y"]
        slowAAs = ["H", "I", "D", "E", "N", "Q", "C", "A", "S", "T", "G", "V", "M"]
        noDataAAs = ["P", "U"]

        NruleDegRate = {}
        NruleDegRate.update(
            (fastAA, fastRate) for fastAA in fastAAs
            )
        NruleDegRate.update(
            (slowAA, slowRate) for slowAA in slowAAs
            )
        NruleDegRate.update(
            (noDataAA, slowRate) for noDataAA in noDataAAs
            ) # Assumed slow rate because of no data

        # Build list of ribosomal proteins
        # Give all ribosomal proteins the slowAA rule
        ribosomalProteins = []
        ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s30_proteins])
        ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s50_proteins])

        # Get degradation rates from measured protein half lives
        measured_deg_rates = {}
        if not disable_measured_protein_deg:
            measured_deg_rates.update({
                p['id']: (np.log(2) / p['half life']).asNumber(deg_rate_units)
                for p in raw_data.protein_half_lives
                })

        degRate = np.zeros(len(raw_data.proteins))
        isRProtein = []
        for i,m in enumerate(raw_data.proteins):
            if m['id'] in measured_deg_rates:
                degRate[i] = measured_deg_rates[m['id']]
                isRProtein.append(m['id'] in ribosomalProteins)
            elif m['id'] in ribosomalProteins:
                degRate[i] = slowRate
                isRProtein.append(True)
            else:
                degRate[i] = NruleDegRate[m['seq'][0]]
                isRProtein.append(False)
        monomerData = np.zeros(
            size,
            dtype = [
                ('id', 'a50'),
                ('rnaId', 'a50'),
                ('degRate', 'f8'),
                ('length', 'i8'),
                ('aaCounts', '{}i8'.format(nAAs)),
                ('mw', 'f8'),
                ('sequence', 'a{}'.format(maxSequenceLength)),
                ('isRProtein', 'bool'),
                ]
            )

        monomerData['id'] = ids
        monomerData['rnaId'] = rnaIds
        monomerData['degRate'] = degRate
        monomerData['length'] = lengths
        monomerData['aaCounts'] = aaCounts
        monomerData['mw'] = mws
        monomerData['sequence'] = sequences
        monomerData['isRProtein'] = isRProtein

        field_units = {
            'id'        :    None,
            'rnaId'        :    None,
            'degRate'    :    deg_rate_units,
            'length'    :    units.aa,
            'aaCounts'    :    units.aa,
            'mw'        :    units.g / units.mol,
            'sequence'  :   None,
            'isRProtein':   None,
            }

        self.monomerData = UnitStructArray(monomerData, field_units)

    def _buildTranslation(self, raw_data, sim_data):
        sequences = self.monomerData["sequence"] # TODO: consider removing sequences

        maxLen = np.int64(
            self.monomerData["length"].asNumber().max()
            + sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
            )

        self.translationSequences = np.empty((sequences.shape[0], maxLen), np.int8)
        self.translationSequences.fill(polymerize.PAD_VALUE)

        aaIDs_singleLetter = sim_data.amino_acid_1_to_3_ordered.keys()

        aaMapping = {aa:i for i, aa in enumerate(aaIDs_singleLetter)}

        for i, sequence in enumerate(sequences):
            for j, letter in enumerate(sequence):
                self.translationSequences[i, j] = aaMapping[letter]

        aaIDs = sim_data.amino_acid_1_to_3_ordered.values()

        self.translationMonomerWeights = (
            (
                sim_data.getter.getMass(aaIDs)
                - sim_data.getter.getMass(["WATER[c]"])
                )
            / raw_data.constants['nAvogadro']
            ).asNumber(units.fg)

        self.translationEndWeight = (sim_data.getter.getMass(["WATER[c]"]) / raw_data.constants['nAvogadro']).asNumber(units.fg)

    def _buildTranslationEfficiency(self, raw_data, sim_data, alternate_translation_efficiency):
        """
        Since the translation efficiency data set from Li et al. 2014 does not
        report a measurement for all genes, genes that do not have a measurement
        are assigned the average translation efficiency.

        If alternate_translation_efficiency is set to True, the translation
        efficiency described by Mohammad et al. 2019 is used instead.
        """
        monomerIds = [x["id"].encode("utf-8") + "[" + sim_data.getter.getLocation([x["id"]])[0][0] + "]" for x in raw_data.proteins]
        monomerIdToGeneId = dict([(x["id"].encode("utf-8") + "[" + sim_data.getter.getLocation([x["id"]])[0][0] + "]", x["geneId"].encode("utf-8")) for x in raw_data.proteins])

        if alternate_translation_efficiency:
            geneIdToTrEff = dict(
                [(x["geneId"].encode("utf-8"), x["translationEfficiency"]) for x in raw_data.translationEfficiency_alternate])
        else:
            geneIdToTrEff = dict([(x["geneId"].encode("utf-8"), x["translationEfficiency"]) for x in raw_data.translationEfficiency if type(x["translationEfficiency"]) == float])

        trEffs = []
        for monomerId in monomerIds:
            geneId = monomerIdToGeneId[monomerId]
            if geneId in geneIdToTrEff:
                trEffs.append(geneIdToTrEff[geneId])
            else:
                trEffs.append(np.nan)

        self.translationEfficienciesByMonomer = np.array(trEffs)
        self.translationEfficienciesByMonomer[np.isnan(self.translationEfficienciesByMonomer)] = np.nanmean(self.translationEfficienciesByMonomer)

    def _build_ribosomal_proteins(self, raw_data, sim_data):
        self.ribosomal_proteins = [
            rprotein['id'] + rprotein['location']
            for rprotein in raw_data.ribosomal_protein_transcripts]

    def _build_elongation_rates(self, raw_data, sim_data):
        self.protein_ids = self.monomerData['id']
        self.ribosomal_protein_ids = sim_data.moleculeGroups.rProteins

        self.protein_indexes = {
            protein: index
            for index, protein in enumerate(self.protein_ids)}

        self.ribosomal_proteins = {
            rprotein: self.protein_indexes.get(rprotein, -1)
            for rprotein in self.ribosomal_protein_ids}

        self.rprotein_indexes = np.array([
            index
            for index in self.ribosomal_proteins.values()
            if index >= 0], dtype=np.int64)

        self.base_elongation_rate = sim_data.constants.ribosomeElongationRateBase.asNumber(units.aa / units.s)
        self.max_elongation_rate = sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
        self.elongation_rates = np.full(
            self.protein_ids.shape,
            self.base_elongation_rate,
            dtype=np.int64)

        self.elongation_rates[self.rprotein_indexes] = self.max_elongation_rate

    def make_elongation_rates_flat(self, base, flat_elongation=False):
        return make_elongation_rates_flat(
            self.protein_ids.shape,
            base,
            self.rprotein_indexes,
            self.max_elongation_rate,
            flat_elongation)

    def make_elongation_rates(
            self,
            random,
            base,
            time_step,
            flat_elongation=False):

        return make_elongation_rates(
            random,
            self.protein_ids.shape,
            base,
            self.rprotein_indexes,
            self.max_elongation_rate,
            time_step,
            flat_elongation)
