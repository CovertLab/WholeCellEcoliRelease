"""
Extracts protein half-lives from Supplementary Table S2 of Nagar et. al.,
"Harnessing Machine Learning to Unravel Protein Degradation in Escherichia coli"
(2021).
"""

import io
import os
import time
from typing import Dict

from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH

# Arbitrarily-determined higher half-life time cutoff in minutes. Proteins with
# half-lives above MAX_HALF_LIFE are all given a half-life of MAX_HALF_LIFE.
MAX_HALF_LIFE = 5000.

# Number of decimals when rounding half-life values (in units of min)
ROUND_N_DECIMALS = 1

# Directories
FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))

# Files
INPUT = os.path.join(FILE_LOCATION, 'msystems.csv')
RNAS_FILE = os.path.join(
    ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'rnas.tsv')
OUTPUT_FLAT_FILE = os.path.join(
    ROOT_PATH, 'reconstruction', 'ecoli', 'flat',
    'protein_half_lives_pulsed_silac.tsv')


def get_symbols_to_monomer_ids():
    # type: () -> Dict[str, str]
    """
    Builds a mapping from gene symbols to the protein monomer IDs, used in later
    data import functions.

    Returns: symbols_to_monomer_ids: Dictionary from gene symbol to first
    associated protein monomer ID.
    """
    symbols_to_monomer_ids = {}

    # Map gene/protein symbols to protein monomer IDs
    with io.open(RNAS_FILE, 'rb') as f:
        reader = tsv.reader(f, delimiter='\t')

        headers = next(reader)
        while headers[0].startswith('#'):
            headers = next(reader)

        gene_symbol_index = headers.index('common_name')
        protein_id_index = headers.index('monomer_ids')

        # TODO(ahzhang): Handle cistrons with two or more associated proteins
        # Currently, each gene symbol in data is associated with the first
        # associated protein id, which is given a half-life in translation.py.
        # Translation.py also only accounts for first associated protein id when
        # assigning half-lives.
        for line in reader:
            gene_symbol = line[gene_symbol_index]
            protein_id = list(line[protein_id_index][2:-2].split('", "'))[0]

            symbols_to_monomer_ids[gene_symbol] = protein_id

        return symbols_to_monomer_ids


def get_half_lives(symbols_to_monomer_ids):
    # type: (Dict[str, str]) -> Dict[str, float]
    """
    Reads protein half-lives from Nagar et. al. (2021) and stores them in a
    dictionary. Any protein listed that does not exist in the model is skipped.

    Args:
        symbols_to_monomer_ids: Dictionary that maps gene symbols to protein
        monomer ID's used by the model.

    Returns:
        id_to_half_life: Dictionary that maps protein monomer ids to their
        half-lives as reported in Nagar et. al. (2021).
    """

    id_to_half_life = {}

    with io.open(INPUT, 'rb') as f:
        reader = tsv.reader(f, delimiter=',')

        headers = next(reader)
        # TODO(ahzhang): change to something of the form
        #  headers.index('Gene names'). Right now this doesn't work because
        #  the raw csv file is encoded in utf-8-sig, which can't be decoded
        #  correctly by binary reading, and gives a \ufeff character at the
        #  beginning of the file.
        gene_symbol_index = 0
        half_life_index = headers.index('half_life')

        for line in reader:
            gene_symbol = line[gene_symbol_index]
            half_life = float(line[half_life_index])

            if gene_symbol in symbols_to_monomer_ids:
                id_to_half_life[symbols_to_monomer_ids[gene_symbol]] = half_life

    return id_to_half_life

def build_half_life_table(raw_half_lives):
    # type: (Dict[str, float]) -> None
    """
    Builds the protein half-life flat file that is used as raw data for
    simulations using them. Assigns all proteins with half-lives greater than
    the arbitrarily chosen MAX_HALF_LIFE as having a half-life of MAX_HALF-LIFE.

    Args:
        raw_half_lives: Dictionary that maps monomer ID's to their half-lives
        as reported in Nagar et. al. (2021).
    """
    # Builds list of tuples that stores half-life data
    half_life_data = [
        (monomer_id, round(min(half_life, MAX_HALF_LIFE), ROUND_N_DECIMALS))
        for (monomer_id, half_life) in raw_half_lives.items()
        ]

    # Sort by monomer ID
    half_life_data.sort(key=lambda v: v[0])

    # Write to flat file
    with io.open(OUTPUT_FLAT_FILE, 'wb') as f:
        print('Writing to {}'.format(f.name))
        writer = tsv.writer(f, quotechar="'", lineterminator='\n')
        writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
        writer.writerow(['"id"', '"half_life (units.min)"'])

        for row in half_life_data:
            writer.writerow([f'"{row[0]}"', f'{row[1]}'])


if __name__ == '__main__':
    symbols_to_monomer_ids = get_symbols_to_monomer_ids()
    raw_half_lives = get_half_lives(symbols_to_monomer_ids)

    build_half_life_table(raw_half_lives)
