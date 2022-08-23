#! /usr/bin/env python

"""
Use Network Component Analysis (NCA) to determine TF regulatory impacts on
each gene based on expression in all conditions.

Note: many variables are named assuming a transcrition factor (TF) is the
regulator but with new features and additional regulation added, this is now
not always the case and could be refactored for clarity.
"""

import argparse
import csv
import os
import re
import sys
import time
from typing import Dict, List, Optional, Set, Tuple, cast

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

import nca


BASE_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')

# Output related
OUTPUT_DIR = os.path.join(BASE_DIR, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)
REGULATION_FILE = 'regulation.tsv'
ACTIVITY_FILE = 'activity.tsv'
HISTOGRAM_FILE = 'histogram.png'
FOLD_CHANGE_FILE = 'fold_changes.tsv'

# Cached results
NETWORK_CACHE_FILE = 'network.npy'
TF_CACHE_FILE = 'tfs.npy'

# RegulonDB related
REGULON_DB_DIR = os.path.join(DATA_DIR, 'regulon-db')
REGULON_DB_SCRIPT = os.path.join(REGULON_DB_DIR, 'download_regulondb.sh')
TF_GENE_FILE = 'tf_genes.tsv'

# EcoCyc regulation
ECOCYC_DIR = os.path.join(DATA_DIR, 'ecocyc')
ECOCYC_REGULATION_FILE = os.path.join(ECOCYC_DIR, 'regulation.tsv')
LEU_LRP_REGULATION_FILE = os.path.join(ECOCYC_DIR, 'leu-lrp-regulation.tsv')

# Sequencing related
WCM_FILE = os.path.join(DATA_DIR, 'wcm_fold_changes.tsv')
COMPENDIUM_DIR = os.path.join(DATA_DIR, 'compendium')
SAMPLES_FILE = os.path.join(COMPENDIUM_DIR, 'samples.tsv')
GENE_NAMES_FILE = os.path.join(COMPENDIUM_DIR, 'gene_names.tsv')
GENE_SYNONYMS_FILE = os.path.join(COMPENDIUM_DIR, 'gene_synonyms.tsv')
SEQ_DIR = os.path.join(COMPENDIUM_DIR, 'seq')
SEQ_FILES = [
    'EcoMAC.tsv',
    # 'RNAseqEcoMACFormat.tsv',
    # 'GSE29076.tsv',
    # 'GSE72525.tsv',
    # 'GSE55662.tsv',
    # 'GSE50529.tsv',
    # 'GSE55365aerobic.tsv',
    # 'GSE55365anaerobic.tsv',
    ]

# Convert some whole-cell model IDs to RegulonDB IDs
TF_MAPPING = {
    'glnG': 'NtrC',
    'ihfA': 'IHF',
    'hns': 'H-NS',
    'yjiE': 'HypT',
    'bglJ': 'RcsB-BglJ',
    }

# Convert some lowercase gene IDs to whole-cell model IDs
WCM_MAPPING = {
    'btsr': 'yehT',
    'cecr': 'ybiH',
    'csqr': 'yihW',
    'glar': 'csiR',
    'pyrr': 'ypdB',
    'sutr': 'ydcN',
    'xynr': 'yagI',
    'alsr': 'rpiR',
    'pgrr': 'ycjZ',
    'rcdA': 'ybjK',
    }

TRNA_MAPPING = {
    'an l-histidyl-[trnahis]': 'tRNA-His',
    'an l-isoleucyl-[trnaile]': 'tRNA-Ile',
    'an l-leucyl-[trnaleu]': 'tRNA-Leu',
    'an l-phenylalanyl-[trnaphe]': 'tRNA-Phe',
    'an l-threonyl-[trnathr]': 'tRNA-Thr',
    'an l-tryptophanyl-[trnatrp]': 'tRNA-Trp',
    'an l-valyl-[trnaval]': 'tRNA-Val',
    }


def load_regulon_db_file(filename: str) -> List[List[str]]:
    """
    Perform check to see if regulonDB data needs to be downloaded before loading
    a tsv file from the regulonDB data directory.

    Args:
        filename: name of file in regulonDB data directory to load

    Returns:
        data: data from a loaded tsv file
    """

    # Check if data exists or if the script needs to be run
    path = os.path.join(REGULON_DB_DIR, filename)
    if not os.path.exists(path):
        raise IOError(f'{path} does not exist. Run {REGULON_DB_SCRIPT} to download regulonDB data.')

    # Load data
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        data = [line for line in reader if not line[0].startswith('#')]

    return data

def load_gene_names() -> Tuple[np.ndarray, np.ndarray, Dict[str, str]]:
    """
    Loads genes names associated with sequencing data rows.

    Returns:
        b_numbers: b number of each gene, ordered to match sequencing data rows
        symbols: gene symbol of each gene, ordered to match sequencing data rows
        synonyms: mapping of gene symbols to b number for gene synonyms
    """

    # Load sequencing data related genes
    b_numbers = []
    symbols = []
    with open(GENE_NAMES_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            b_numbers.append(line[1])
            symbols.append(line[2])

    # Load gene synonyms
    synonyms = {}
    with open(GENE_SYNONYMS_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            b = line[0]
            for synonym in line[1:]:
                synonyms[synonym] = b

    return np.array(b_numbers), np.array(symbols), synonyms

def load_seq_data(linearize: bool, average: bool) -> Tuple[np.ndarray, Dict[int, int]]:
    """
    Load sequencing data from the compendium.

    Args:
        linearize: if set, counts will be transformed from log2 space to linear space
        average: if True, average sequencing in replicate samples

    Returns:
        data: matrix of normalized sequencing data representing counts (n genes, m samples)
            linear if linearize, log2 otherwise
        idx_mapping: mapping of loaded sample column index to new sample column index
            which will be a 1:1 mapping if average is False

    TODO:
    - normalize data to include other files instead of just EcoMAC in SEQ_FILES
    """

    idx_mapping = {}
    seq_data = []
    for filename in SEQ_FILES:
        path = os.path.join(SEQ_DIR, filename)

        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            seq_data.append(list(reader))

    data = np.hstack(seq_data).astype(np.float64)

    if average:
        with open(SAMPLES_FILE) as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)

            replicate_col = header.index('Replicate')
            mask = np.ones(len(header), bool)
            mask[:2] = False
            mask[replicate_col] = False
            old_line = np.array(header)[mask]
            averaged_data = []
            start_idx = None
            old_replicate = 1
            new_idx = -1
            for i, line in enumerate(reader):
                replicate = int(line[replicate_col][-1])
                new_line = np.array(line)[mask]
                if replicate != old_replicate + 1 or np.any(new_line != old_line):
                    if start_idx is not None:
                        averaged_data.append(data[:, start_idx:i].mean(1))

                    start_idx = i
                    new_idx += 1

                old_replicate = replicate
                old_line = new_line
                idx_mapping[i] = new_idx
            averaged_data.append(data[:, start_idx:].mean(1))

            data = np.hstack((np.vstack(averaged_data).T, data[:, i+1:]))

    if linearize:
        data = 2**data

    return data, idx_mapping

def load_sigma_gene_interactions() -> Dict[str, Dict[str, int]]:
    """
    Load regulonDB Sigma factor-gene interactions.

    Returns:
        sigma_genes: relationship between sigma factors and genes
            {sigma factor: {gene: regulatory direction}}

    TODO:
        - add option for EcoCyc sigma factors
    """

    data = load_regulon_db_file('sigma_genes.tsv')

    sigma_genes = {}  # type: Dict[str, Dict[str, int]]
    sigma_idx = 0
    gene_idx = 1
    dir_idx = 2
    # evidence_idx = 4  # TODO: filter on evidence?
    for line in data:
        sigma_factors = line[sigma_idx].split(', ')
        gene = line[gene_idx]
        effect = line[dir_idx]

        if effect != 'activator':
            raise ValueError(f'Unrecognized sigma factor effect: {effect}')

        for sigma_factor in sigma_factors:
            if sigma_factor not in sigma_genes:
                sigma_genes[sigma_factor] = {}
            sigma_genes[sigma_factor][gene] = 1

    return sigma_genes

def add_tf_gene_data(
        tf_genes: Dict[str, Dict[str, int]],
        tf: str,
        gene: str,
        effect: str,
        activator_key: str,
        repressor_key: str,
        ambiguous_key: str,
        split: bool,
        verbose: bool,
        ) -> None:
    """
    Add TF-gene interaction data to a data structure from multiple sources.

    Args:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
            regulatory direction is 1 for positive regulation
            regulatory direction is -1 for negative regulation
            regulatory direction is 0 for unknown or conflicting regulation
        tf: TODO
        gene: TODO
        effect: TODO
        activator_key: TODO
        repressor_key: TODO
        ambiguous_key: TODO
        split: if True, splits TFs into activator and repressor forms
        verbose: if True, prints warnings about loaded data
    """

    # Check type of regulation
    if effect == activator_key:
        direction = 1
    elif effect == repressor_key:
        direction = -1
    elif effect == ambiguous_key:
        direction = 0
    else:
        raise ValueError(f'Unrecognized TF effect: {effect}')

    # Store new data
    if split:
        if effect == activator_key or effect == ambiguous_key:
            split_tf = f'{tf}:activator'
            genes = tf_genes.get(split_tf, {})
            genes[gene] = 1
            tf_genes[split_tf] = genes

        if effect == repressor_key or effect == ambiguous_key:
            split_tf = f'{tf}:repressor'
            genes = tf_genes.get(split_tf, {})
            genes[gene] = -1
            tf_genes[split_tf] = genes
    else:
        genes = tf_genes.get(tf, {})
        if genes.get(gene, direction) != direction:
            if verbose:
                print(f'Inconsistent regulatory effects for {tf} on {gene}')
            genes[gene] = 0
        else:
            genes[gene] = direction
        tf_genes[tf] = genes

def load_ecocyc_tf_gene_interactions(
        valid_types: Set[str],
        split: bool = False,
        verbose: bool = True,
        ) -> Dict[str, Dict[str, int]]:
    """
    Load EcoCyc TF-gene interactions.

    Args:
        valid_types: types of regulation that should be included in the NCA problem
        split: if True, splits TFs into activator and repressor forms
        verbose: if True, prints warnings about loaded data

    Returns:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
            regulatory direction is 1 for positive regulation
            regulatory direction is -1 for negative regulation
            regulatory direction is 0 for unknown or conflicting regulation
    """

    tf_genes = {}  # type: Dict[str, Dict[str, int]]
    positive = '+'
    negative = '-'
    unknown = 'NIL'

    # Load Leu-Lrp regulation to differentiate from just Lrp regulation.
    # Leu-Lrp will have roughly the opposite effect of Lrp but they are not
    # distinct transcription factors in ECOCYC_REGULATION_FILE.
    leu_lrp_regulation = {}
    with open(LEU_LRP_REGULATION_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            if line[0].startswith('#'):
                continue

            if line[1] == 'x' and line[2] == 'x':
                reg = unknown
            elif line[1] == 'x':
                reg = positive
            else:
                reg = negative

            leu_lrp_regulation[line[0]] = reg

    with open(ECOCYC_REGULATION_FILE) as f:
        reader = csv.reader(f, delimiter='\t')

        type_idx = 2
        dir_idx = 3
        common_name_idx = 4
        for line in reader:
            # Skip lines that are not needed
            if line[0].startswith('#'):
                continue
            if line[type_idx] not in valid_types:
                continue

            # Extract data from the line
            tf, gene = re.findall('(.*) [-+]*> (.*)', line[common_name_idx])[0]
            effect = line[dir_idx]

            # Handle special case for Leu-Lrp vs just Lrp regulation
            if tf == 'lrp':
                if gene in leu_lrp_regulation:
                    effect = leu_lrp_regulation[gene]
                elif effect == positive:
                    effect = negative
                elif effect == negative:
                    effect = positive

            add_tf_gene_data(tf_genes, tf, gene, effect, positive, negative, unknown, split, verbose)

    return tf_genes

def load_regulondb_tf_gene_interactions(
        split: bool = False,
        verbose: bool = True,
        ) -> Dict[str, Dict[str, int]]:
    """
    Load regulonDB TF-gene interactions.

    Args:
        split: if True, splits TFs into activator and repressor forms
        verbose: if True, prints warnings about loaded data

    Returns:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
            regulatory direction is 1 for positive regulation
            regulatory direction is -1 for negative regulation
            regulatory direction is 0 for unknown or conflicting regulation
    """

    data = load_regulon_db_file(TF_GENE_FILE)

    tf_genes = {}  # type: Dict[str, Dict[str, int]]
    tf_idx = 0
    gene_idx = 1
    dir_idx = 2
    # evidence_idx = 4  # TODO: filter on evidence?
    for line in data:
        # Extract columns of interest
        tf = line[tf_idx]
        gene = line[gene_idx]
        effect = line[dir_idx]
        add_tf_gene_data(tf_genes, tf, gene, effect, 'activator', 'repressor', 'unknown', split, verbose)

    return tf_genes

def load_tf_gene_interactions(
        valid_types: Set[str],
        ecocyc: bool = True,
        regulondb: bool = False,
        split: bool = False,
        verbose: bool = True,
        ) -> Dict[str, Dict[str, int]]:
    """
    Load TF-gene interactions.

    Args:
        valid_types: types of regulation that should be included in the NCA problem
            (only used for EcoCyc)
        ecocyc: if True, uses regulation interactions from EcoCyc
        regulondb: if True, uses regulation interactions from RegulonDB
        split: if True, splits TFs into activator and repressor forms
        verbose: if True, prints warnings about loaded data

    Returns:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
            regulatory direction is 1 for positive regulation
            regulatory direction is -1 for negative regulation
            regulatory direction is 0 for unknown or conflicting regulation

    TODO:
        - create consistent naming convention from EcoCyc and RegulonDB (eg EcoCyc TF arcA vs RegulonDB ArcA)
        - check for clashes in regulation direction if using both ecocyc and regulondb
    """

    if ecocyc and regulondb:
        raise NotImplementedError('Using both EcoCyc and RegulonDB requires reconciling differences in TF/gene IDs.')

    tf_genes = {}
    if ecocyc:
        tf_genes.update(load_ecocyc_tf_gene_interactions(valid_types, split=split, verbose=verbose))
        if verbose:
            print('Loaded EcoCyc gene regulation')

    if regulondb:
        tf_genes.update(load_regulondb_tf_gene_interactions(split=split, verbose=verbose))
        if verbose:
            print('Loaded RegulonDB gene regulation')

    return tf_genes

def load_wcm_fold_changes() -> Dict[str, Dict[str, Tuple[float, int]]]:
    """
    Load fold changes currently implemented in the whole-cell model.

    Returns:
        fold_changes: regulatory pairs for TF-gene interactions
            {tf: {gene: (log2 FC, direction)}}
    """

    fold_changes = {}  # type: Dict[str, Dict[str, Tuple[float, int]]]

    with open(WCM_FILE) as f:
        reader = csv.reader(f, delimiter='\t')

        next(reader)  # remove comment line
        headers = next(reader)
        tf_idx = headers.index('TF')
        gene_idx = headers.index('Target')
        fc_idx = headers.index('F_avg')
        direction_idx = headers.index('Regulation_direct')

        for line in reader:
            if line[0].startswith('#'):
                continue

            tf = line[tf_idx]
            gene = line[gene_idx]
            fc = float(line[fc_idx])
            direction = int(line[direction_idx])

            data = fold_changes.get(tf, {})
            data.update({gene: (fc, direction)})
            fold_changes[tf] = data

    return fold_changes

def create_tf_map(
        gene_names: np.ndarray,
        synonyms: Dict[str, str],
        tf_genes: Dict[str, Dict[str, int]],
        verbose: bool = True,
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create an initial map reflecting the known regulatory network topology between
    transcription factors and genes. Also includes a column of ones for constitutive
    expression of each gene.

    Args:
        gene_names: gene IDs corresponding to each row of expression matrix
        synonyms: mapping of gene symbols to b number for gene synonyms
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        verbose: if True, prints warnings about loaded data

    Returns:
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            positive regulation: positive number
            negative regulation: negative number
            unknown/ambiguous regulation: NaN
            no regulation: 0
        tfs: IDs of TFs associated with each column of mapping
    """

    np.random.seed(0)

    gene_idx = {name: idx for idx, name in enumerate(gene_names)}
    n_genes = len(gene_names)
    n_tfs = len(tf_genes)
    mapping = np.zeros((n_genes, n_tfs))

    # Populate matrix with links between TFs and genes
    tfs = []
    for j, (tf, genes) in enumerate(sorted(tf_genes.items())):
        for gene, direction in genes.items():
            b = synonyms.get(gene)
            if b is None or b not in gene_idx:
                if verbose:
                    print(f'Unknown transcription factor gene: {gene}')
                continue

            if direction == 0:
                mapping[gene_idx[b], j] = np.nan
            else:
                mapping[gene_idx[b], j] = direction * np.random.rand()
        tfs.append(tf)

    return mapping, np.array(tfs)

def save_regulation(
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        output_dir: str,
        ) -> None:
    """
    Save the results of NCA to a file showing each TF/gene pair.

    Args:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: TF activity for each condition
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        output_dir: path to directory to save the output
    """

    # Regulation data (TF-gene pairs)
    relation_idx = np.where(A.T)
    with open(os.path.join(output_dir, REGULATION_FILE), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['TF', 'Gene', 'FC'])
        for tf_idx, gene_idx in zip(*relation_idx):
            writer.writerow([tfs[tf_idx], genes[gene_idx], A[gene_idx, tf_idx]])

    # Activity data (TF in each condition)
    with open(os.path.join(output_dir, ACTIVITY_FILE), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Condition:'] + list(range(P.shape[1])))  # type: ignore
        for tf, activity in zip(tfs, P):
            writer.writerow([tf] + list(activity))

def load_regulation(directory: str, genes: np.ndarray, tfs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load the results of a previous NCA run saved to file with save_regulation.
    Return values should match the arguments passed to save_regulation.

    Args:
        directory: path to folder containing saved NCA results
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)

    Returns:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: TF activity for each condition
    """

    with open(os.path.join(directory, REGULATION_FILE)) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = np.array(list(reader))

    fc_genes = data[:, headers.index('Gene')]
    fc_tfs = data[:, headers.index('TF')]
    fcs = data[:, headers.index('FC')].astype(float)

    # Recreate A matrix
    A = np.zeros((len(genes), len(tfs)))
    gene_idx = {gene: idx for idx, gene in enumerate(genes)}
    tf_idx = {tf: idx for idx, tf in enumerate(tfs)}
    for gene, tf, fc in zip(fc_genes, fc_tfs, fcs):
        A[gene_idx[gene], tf_idx[tf]] = fc

    # Recreate P matrix
    with open(os.path.join(directory, ACTIVITY_FILE)) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # strip headers
        P = np.array(list(reader))[:, 1:].astype(float)

    return A, P

def add_global_expression(
        tfs: np.ndarray,
        mapping: Optional[np.ndarray] = None,
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Expand out TFs to include global expression (consituitive and regulated).
    This will capture genes not captured by TFs.

    Args:
        tfs: IDs of TFs associated with each column of mapping
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            if None, no update is performed

    Returns:
        tfs: updated IDs with global expression IDs
        mapping: updated matrix with global expression columns

    TODO:
        - add ppGpp
    """

    tfs = np.hstack((np.array(['constituitive', 'regulated']), tfs))

    if mapping is None:
        mapping = np.array([])
    else:
        n_genes = mapping.shape[0]
        no_regulation = np.sum(mapping, axis=1) == 0

        # Add constituitive expression for any gene not regulated
        constituitive = np.zeros((n_genes, 1))
        constituitive[no_regulation] = 1

        # Add regulated expression for general expression level of regulated genes
        regulated = np.zeros((n_genes, 1))
        regulated[~no_regulation] = 1

        mapping = np.hstack((constituitive, regulated, mapping))

    return tfs, cast(np.ndarray, mapping)

def add_sigma_factors(
        tfs: np.ndarray,
        genes: np.ndarray,
        synonyms: Dict[str, str],
        mapping: Optional[np.ndarray] = None,
        verbose: bool = False,
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Expand out TFs to include sigma factors. This will capture some genes not
    regulated by TFs.

    Args:
        tfs: IDs of TFs associated with each column of mapping
        genes: IDs of genes associated with each row of mapping
        synonyms: mapping of gene symbols to b number for gene synonyms
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            if None, no update is performed
        verbose: if True, prints warnings about unmatched genes

    Returns:
        tfs: updated IDs with sigma factor IDs
        mapping: updated matrix with sigma factor columns
    """

    sigma_genes = load_sigma_gene_interactions()
    sigma_factors = list(sigma_genes.keys())
    tfs = np.hstack((np.array(sigma_factors), tfs))

    if mapping is None:
        mapping = np.array([])
    else:
        n_genes = mapping.shape[0]
        sigma_regulation = np.zeros((n_genes, len(sigma_factors)))
        gene_idx = {gene: idx for idx, gene in enumerate(genes)}
        sigma_idx = {sigma: idx for idx, sigma in enumerate(sigma_factors)}
        for sigma_factor, regulation in sigma_genes.items():
            for gene, direction in regulation.items():
                b = synonyms.get(gene)
                if b in gene_idx:
                    sigma_regulation[gene_idx[b], sigma_idx[sigma_factor]] = direction
                elif verbose:
                    print(f'Unknown sigma factor gene: {gene}')
        mapping = np.hstack((sigma_regulation, mapping))

    return tfs, cast(np.ndarray, mapping)

def add_noisy_expression(
        tfs: np.ndarray,
        noise: int,
        repeats: int,
        seed: int,
        mapping: Optional[np.ndarray] = None,
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Expand out TFs to include random expression but assigning genes to multiple
    "transcription factors".

    Args:
        tfs: IDs of TFs associated with each column of mapping
        noise: number of TFs to include for noise
        repeats: number of times to repeat a gene in noise TFs
        seed: random seed for assigning genes
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            if None, no update is performed

    Returns:
        tfs: updated IDs with global expression IDs
        mapping: updated matrix with global expression columns

    TODO:
        - verify this leads to reasonable results and doesn't cause issues
    """

    if mapping is None:
        mapping = np.array([])
    else:
        np.random.seed(seed)
        n_genes = mapping.shape[0]
        n_noisy = noise * repeats
        tfs = np.hstack(([f'noise-{i}' for i in range(n_noisy)], tfs))
        for repeat in range(repeats):
            idx = np.arange(n_genes)
            np.random.shuffle(idx)
            for n in range(noise):
                new = np.zeros((n_genes, 1))
                new[idx[n_genes * n // noise:n_genes * (n + 1) // noise]] = 1
                mapping = np.hstack((new, mapping))

    return tfs, cast(np.ndarray, mapping)

def calculate_fold_changes(
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        n_std: float = 1.,
        min_conditions: int = 10,
        ) -> Tuple[np.ndarray, Dict[str, Dict[str, float]]]:
    """
    Calculate expected expression fold changes from NCA solutions.

    Args:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF/condition relationship (m TFs, o conditions)
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        n_std: TODO
        min_conditions: TODO

    Returns:
        fcs: fold changes corresponding to entries in A (n genes, m TFs)
        regulatory_pairs: fold changes for TFs and genes from NCA solution
            {TF: {gene: fold change}}
    """

    # Find the upper and lower bounds based on number of standard deviations
    # while maintaining a minimum number of conditions to keep data for
    sort = np.sort(P, 1)
    lower = np.fmax(P.mean(1) - n_std*P.std(1), sort[:, min_conditions])
    upper = np.fmin(P.mean(1) + n_std*P.std(1), sort[:, -min_conditions])

    # Exclude data that does not fall outside the limits
    low_P = P.copy().T
    low_P[low_P > lower] = 0
    high_P = P.copy().T
    high_P[high_P < upper] = 0

    # Determine fold change as TF effect on gene multiplied by the difference
    # in activity between average high and low conditions
    fcs = A * (high_P.sum(0) / np.sum(high_P != 0, 0) - low_P.sum(0) / np.sum(low_P != 0, 0))
    regulatory_pairs = {}  # type: Dict[str, Dict[str, float]]
    for i, gene in enumerate(genes):
        for j, tf in enumerate(tfs):
            processed_tf = tf.split(':')[0].lower()
            data = regulatory_pairs.get(processed_tf, {})
            data[gene] = data.get(gene, 0) + fcs[i, j]
            if data[gene] == 0:
                del data[gene]
            regulatory_pairs[processed_tf] = data

    return fcs, regulatory_pairs

def compare_wcm_nca(nca_pairs: Dict[str, Dict[str, float]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get matching arrays of fold changes from the whole-cell model and NCA results
    for comparison and plotting.

    Args:
        nca_pairs: fold changes for TFs and genes from NCA solution
            {TF: {gene: fold change}}

    Returns:
        wcm_fcs: fold change for each TF/gene pair in the whole-cell model
        nca_fcs: fold change for each TF/gene pair from the NCA solution
        wcm_consistent: True if whole-cell model fold change direction is
            consistent with annotations
    """
    wcm_fcs = []
    nca_fcs = []
    wcm_consistent = []
    for tf, regulation in load_wcm_fold_changes().items():
        for gene, (fc, direction) in regulation.items():
            wcm_fcs.append(fc * np.sign(direction))
            nca_fc = nca_pairs.get(TF_MAPPING.get(tf, tf).lower(), {}).get(gene, 0)
            nca_fcs.append(nca_fc)
            wcm_consistent.append(np.abs(direction) == 1)

    return np.array(wcm_fcs), np.array(nca_fcs), np.array(wcm_consistent)

def match_statistics(
        E: np.ndarray,
        A: np.ndarray,
        P: np.ndarray,
        tfs: np.ndarray,
        tf_genes: Dict[str, Dict[str, int]],
        genes: np.ndarray,
        ) -> None:
    """
    Assess the accuracy of NCA results by printing statistics about how well
    the regulation direction and overall expression is captured.

    Args:
        E: original expression data (n genes, m conditions)
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF activity in each condition (o TFs, m conditions)
        tfs: names of each TF corresponding to columns in A (m TFs)
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        genes: IDs for each gene corresponding to rows in A (n genes)

    TODO:
        - stats for each TF
        - track skipped TFs and genes needed to satisfy NCA constraints
    """

    # Variable to track stats
    total_neg = 0
    correct_neg = 0
    dropped_neg = 0
    total_pos = 0
    correct_pos = 0
    dropped_pos = 0
    ambiguous = 0

    gene_idx = {gene: i for i, gene in enumerate(genes)}
    tf_idx = {tf: i for i, tf in enumerate(tfs)}

    # Check each entry in the mapping matrix against the annotated data
    for tf, regulation in tf_genes.items():
        if tf not in tf_idx:
            continue

        for gene, annotated in regulation.items():
            if gene not in gene_idx:
                continue

            predicted = A[gene_idx[gene], tf_idx[tf]]
            if annotated == 0:
                ambiguous += 1
            elif annotated > 0:
                total_pos += 1
                if predicted > 0:
                    correct_pos += 1
                elif predicted == 0:
                    dropped_pos += 1
            else:
                total_neg += 1
                if predicted < 0:
                    correct_neg += 1
                elif predicted == 0:
                    dropped_neg += 1

    # Print statistics
    annotated_neg = total_neg - dropped_neg
    annotated_pos = total_pos - dropped_pos
    correct = correct_neg + correct_pos
    dropped = dropped_neg + dropped_pos
    matched = annotated_neg + annotated_pos
    total = matched + dropped + ambiguous
    percent_dropped = 0 if total == 0 else 100*dropped/total
    percent_match = 0 if matched == 0 else 100*correct/matched
    percent_match_neg = 0 if annotated_neg == 0 else 100*correct_neg/annotated_neg
    percent_match_pos = 0 if annotated_pos == 0 else 100*correct_pos/annotated_pos
    print(f'Dropped regulation: {dropped}/{total} ({percent_dropped:.1f}%)')
    print(f'Overall matches: {correct}/{matched} ({percent_match:.1f}%)')
    print(f'Negative regulation matches: {correct_neg}/{annotated_neg} ({percent_match_neg:.1f}%)')
    print(f'Positive regulation matches: {correct_pos}/{annotated_pos} ({percent_match_pos:.1f}%)')
    print(f'Number of ambiguous regulatory interactions: {ambiguous}')

    # E prediction from results
    E_est = A.dot(P)
    no_prediction = E_est == 0
    n_zeros = np.sum(no_prediction)
    n_negatives = np.sum(E_est < 0)
    n_samples = E.shape[0] * E.shape[1]
    error = np.sqrt(np.mean(((E - E_est) / E)**2))
    error_predictions = np.sqrt(np.mean(((E[~no_prediction] - E_est[~no_prediction]) / E[~no_prediction])**2))
    print(f'\nError fitting original data: {error:.3f}')
    print(f'Error fitting original data (excluding unpredicted samples): {error_predictions:.3f}')
    print(f'No prediction made: {n_zeros}/{n_samples} ({100 * n_zeros / n_samples:.1f}%)')
    print(f'Negative prediction made: {n_negatives}/{n_samples} ({100 * n_negatives / n_samples:.1f}%)')

    # WCM fold change correlations
    fcs, nca_pairs = calculate_fold_changes(A, P, genes, tfs)
    wcm_fcs, nca_fcs, wcm_consistent = compare_wcm_nca(nca_pairs)
    nonzero = nca_fcs != 0
    pearson_all = stats.pearsonr(wcm_fcs, nca_fcs)
    pearson_consistent = stats.pearsonr(wcm_fcs[wcm_consistent], nca_fcs[wcm_consistent])
    if np.any(nonzero):
        pearson_no_zeros = stats.pearsonr(wcm_fcs[nonzero], nca_fcs[nonzero])
    else:
        pearson_no_zeros = [np.nan, np.nan]
    print('\nFold change correlation with WCM:')
    print(f'All: r={pearson_all[0]:.3f} (p={pearson_all[1]:.0e}, n={len(wcm_fcs)})')
    print(f'Only if WCM consistent: r={pearson_consistent[0]:.3f} (p={pearson_consistent[1]:.0e}, n={np.sum(wcm_consistent)})')
    print(f'Excluding unmatched: r={pearson_no_zeros[0]:.3f} (p={pearson_no_zeros[1]:.0e}, n={np.sum(nonzero)})')

def plot_results(
        tf_genes: Dict[str, Dict[str, int]],
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        output_dir: str,
        max_bins: int = 200,
        ) -> None:
    """
    Plot NCA results for easier inspection.

    Args:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF/condition relationship (m TFs, o conditions)
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        output_dir: path to directory to save the plot
        max_bins: maximum number of bins for histogram plots
    """

    def plot(ax, series, label, hist_range, n_bins, color):
        mean = sum(series) / len(series) if len(series) else 0
        ax.hist(series, color=color, bins=n_bins, range=hist_range, alpha=0.5, label=label)
        ax.axvline(mean, color=color, linestyle='--', label=f'{label} mean: {mean:.2f}')

    def no_spines(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    negative = np.zeros(A.shape, bool)
    positive = np.zeros(A.shape, bool)
    ambiguous = np.zeros(A.shape, bool)

    # Check each entry in the mapping matrix against the annotated data
    for i, j in zip(*np.where(A)):
        gene = genes[i]
        tf = tfs[j]
        annotated = tf_genes.get(tf, {}).get(gene)

        if annotated is not None:
            if annotated < 0:
                negative[i, j] = True
            elif annotated > 0:
                positive[i, j] = True
            else:
                ambiguous[i, j] = True

    cmap = plt.get_cmap('tab10')

    # TF regulation types
    positive_tfs = np.array([np.all([v > 0 for v in tf_genes.get(tf, {}).values()]) for tf in tfs])
    negative_tfs = np.array([np.all([v < 0 for v in tf_genes.get(tf, {}).values()]) for tf in tfs])
    combined_tfs = (~positive_tfs) & (~negative_tfs)

    # Get fold changes from NCA and WCM
    fcs, nca_pairs = calculate_fold_changes(A, P, genes, tfs)
    wcm_fcs, nca_fcs, wcm_consistent = compare_wcm_nca(nca_pairs)

    # Plot figure
    plt.figure(figsize=(12, 18))
    gs = gridspec.GridSpec(3, 2)
    ax_A = plt.subplot(gs[0, 0])
    ax_fc = plt.subplot(gs[1, 0])
    ax_wcm = plt.subplot(gs[2, 0])
    ax_P1 = plt.subplot(gs[0, 1])
    ax_P2 = plt.subplot(gs[1, 1])
    ax_P3 = plt.subplot(gs[2, 1])

    ## Plot A results
    hist_range = (np.floor(A.min()), np.ceil(A.max()))
    n_bins = min(5 * int(hist_range[1] - hist_range[0]), max_bins)
    plot(ax_A, A[negative], 'Negative', hist_range, n_bins, cmap(0))
    plot(ax_A, A[positive], 'Positive', hist_range, n_bins, cmap(1))
    plot(ax_A, A[ambiguous], 'Ambiguous', hist_range, n_bins, cmap(2))
    ax_A.legend(fontsize=8, frameon=False)
    ax_A.set_xlabel('TF-Gene Interation\n(A entries)', fontsize=10)
    ax_A.set_ylabel('Count', fontsize=10)
    no_spines(ax_A)

    ## Plot fold changes
    hist_range = (np.floor(fcs.min()), np.ceil(fcs.max()))
    n_bins = min(5 * int(hist_range[1] - hist_range[0]), max_bins)
    plot(ax_fc, fcs[negative], 'Negative', hist_range, n_bins, cmap(0))
    plot(ax_fc, fcs[positive], 'Positive', hist_range, n_bins, cmap(1))
    plot(ax_fc, fcs[ambiguous], 'Ambiguous', hist_range, n_bins, cmap(2))
    ax_fc.legend(fontsize=8, frameon=False)
    ax_fc.set_xlabel('TF-Gene Fold Change', fontsize=10)
    ax_fc.set_ylabel('Count', fontsize=10)
    no_spines(ax_fc)

    ## Compare to whole-cell model fold changes
    nonzero = nca_fcs != 0
    pearson_all = stats.pearsonr(wcm_fcs, nca_fcs)
    pearson_consistent = stats.pearsonr(wcm_fcs[wcm_consistent], nca_fcs[wcm_consistent])
    if np.any(nonzero):
        pearson_no_zeros = stats.pearsonr(wcm_fcs[nonzero], nca_fcs[nonzero])
    else:
        pearson_no_zeros = [np.nan, np.nan]
    ax_wcm.plot(wcm_fcs[wcm_consistent], nca_fcs[wcm_consistent], 'x', label='WCM consistent')
    ax_wcm.plot(wcm_fcs[~wcm_consistent], nca_fcs[~wcm_consistent], 'x', label='WCM inconsistent')
    xlim = ax_wcm.get_xlim()
    ylim = ax_wcm.get_ylim()
    min_val = min(wcm_fcs.min(), nca_fcs.min())
    max_val = max(wcm_fcs.max(), nca_fcs.max())
    ax_wcm.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')
    ax_wcm.set_xlim(xlim)
    ax_wcm.set_ylim(ylim)
    ax_wcm.set_xlabel('Whole-cell model fold change')
    ax_wcm.set_ylabel('NCA predicted fold change')
    ax_wcm.legend(fontsize=8, frameon=False)
    ax_wcm.set_title(f'All: r={pearson_all[0]:.3f} (p={pearson_all[1]:.0e}, n={len(wcm_fcs)})\n'
        f'Consistent: r={pearson_consistent[0]:.3f} (p={pearson_consistent[1]:.0e}, n={np.sum(wcm_consistent)})\n'
        f'Exclude zeros: r={pearson_no_zeros[0]:.3f} (p={pearson_no_zeros[1]:.0e}, n={np.sum(nonzero)})',
        fontsize=10)
    no_spines(ax_wcm)

    ## Plot P results
    P_ave = P.mean(1)
    hist_range = (np.floor(P_ave.min()), np.ceil(P_ave.max()))
    n_bins = min(5 * int(hist_range[1] - hist_range[0]), max_bins)
    plot(ax_P1, P_ave[negative_tfs], 'All negative', hist_range, n_bins, cmap(0))
    plot(ax_P1, P_ave[positive_tfs], 'All positive', hist_range, n_bins, cmap(1))
    plot(ax_P1, P_ave[combined_tfs], 'Multiple', hist_range, n_bins, cmap(2))
    ax_P1.legend(fontsize=8, frameon=False)
    ax_P1.set_xlabel('TF Average Activity\n(mean of P rows)', fontsize=10)
    ax_P1.set_ylabel('Count', fontsize=10)
    no_spines(ax_P1)

    ## Plot P range results
    P_range = P.max(1) - P.min(1)
    hist_range = (np.floor(P_range.min()), np.ceil(P_range.max()))
    n_bins = min(2 * int(hist_range[1] - hist_range[0]), max_bins)
    plot(ax_P2, P_range[negative_tfs], 'All negative', hist_range, n_bins, cmap(0))
    plot(ax_P2, P_range[positive_tfs], 'All positive', hist_range, n_bins, cmap(1))
    plot(ax_P2, P_range[combined_tfs], 'Multiple', hist_range, n_bins, cmap(2))
    ax_P2.legend(fontsize=8, frameon=False)
    ax_P2.set_xlabel('TF Activity Range\n(range of P rows)', fontsize=10)
    ax_P2.set_ylabel('Count', fontsize=10)
    no_spines(ax_P2)

    ## Plot P distributions
    for tf_activity in P:
        ax_P3.hist(tf_activity, linewidth=1, alpha=0.2, histtype='step')
    ax_P3.set_xlabel('TF Activity\n(P rows)', fontsize=10)
    ax_P3.set_ylabel('Count', fontsize=10)
    no_spines(ax_P3)

    ## Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, HISTOGRAM_FILE))
    plt.close('all')

def save_fold_changes(
        valid_types: Set[str],
        ecocyc: bool,
        regulondb: bool,
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        output_dir: str,
        attenuation: bool = False,
        ) -> None:
    """
    Save fold changes to a file to use in the whole-cell model.

    Args:
        valid_types: types of regulation that should be included in the NCA problem
        ecocyc: if True, uses regulation interactions from EcoCyc
        regulondb: if True, uses regulation interactions from RegulonDB
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF/condition relationship (m TFs, o conditions)
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        output_dir: path to directory to save the data
        attenuation: if True, sets different headers for output file

    TODO:
        - strip columns that are not needed after updating WCM load
    """

    curated_data = load_tf_gene_interactions(
        valid_types, ecocyc=ecocyc, regulondb=regulondb, split=False, verbose=False)
    output_file = os.path.join(output_dir, FOLD_CHANGE_FILE)
    _, nca_pairs = calculate_fold_changes(A, P, genes, tfs)
    lowercase_mapping = {k.lower(): k for k in load_gene_names()[2]}
    lowercase_mapping.update(WCM_MAPPING)

    with open(output_file, 'w') as f:
        # Write file generation comments
        metadata = f'# Generated by {__file__} on {time.ctime()}\n'
        f.write(metadata)
        command = f'# Command: {" ".join(sys.argv)}\n'
        f.write(command)

        writer = csv.writer(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)

        # Write header
        if attenuation:
            header = ['tRNA', 'Target', 'log2 FC']
        else:
            header = ['TF', 'Target', 'log2 FC mean', 'log2 FC std', 'Regulation_direct']
        writer.writerow(header)

        # Write row for each regulatory interaction
        for tf, nca_genes in sorted(nca_pairs.items()):
            if tf not in lowercase_mapping:
                print(f'Warning: did not find a match for {tf}')
            tf = lowercase_mapping.get(tf, tf)

            for gene, fc in sorted(nca_genes.items()):
                curated_dir = curated_data.get(tf, {}).get(gene, 0)

                # Skip fold changes that do not agree with curated directions
                if curated_dir * fc < 0:
                    continue

                tf = TRNA_MAPPING.get(tf, tf)
                writer.writerow([tf, gene, round(fc, 2)])

def parse_args() -> argparse.Namespace:
    """Parse command line args for options to run."""

    parser = argparse.ArgumentParser()

    default_nca = nca.METHODS[0]
    default_label = 'nca-results'
    default_iterative_iterations = 100
    default_iterative_split = 5
    default_robust_iterations = 5
    default_status_update = 0.1

    # General options
    parser.add_argument('-l', '--label',
        default=default_label,
        help=f'Label for output directory to save results to (default: {default_label}).')
    parser.add_argument('-v', '--verbose',
        action='store_true',
        help='If set, prints status updates for creating initial network map.')

    # Options for efficient analysis by reusing saved results
    parser.add_argument('-a', '--analysis',
        help='Path to directory containing saved regulation data to just run analysis on. Will skip NCA.')
    parser.add_argument('-c', '--cache',
        help='Path to cache directory to load network files. Defaults to output directory if not specified.')
    parser.add_argument('-f', '--force',
        action='store_true',
        help='If set, force a rerun of identifying the initial TF map, otherwise use cached values.')

    # Network options
    parser.add_argument('--no-ecocyc', dest='ecocyc',
        action='store_false',
        help='If set, does not use EcoCyc regulation data.')
    parser.add_argument('--regulondb',
        action='store_true',
        help='If set, uses RegulonDB regulation data.')

    # Data options
    parser.add_argument('--attenuation-only',
        action='store_true',
        help='If set, only looks at tRNA attenuation fold changes.')
    parser.add_argument('--average-seq',
        action='store_true',
        help='If set, averages sequencing data for replicate samples and reduces the number of samples.')
    parser.add_argument('--linear',
        action='store_true',
        help='If set, use linear counts from sequencing data otherwise keep log2 counts.')
    parser.add_argument('--global-expression',
        action='store_true',
        help='If set, creates pseudo transcription factors to capture global expression capacity (basal rates).')
    parser.add_argument('--split',
        action='store_true',
        help='If set, split transcription factors into positive and negative regulation.')
    parser.add_argument('--sigma-factors',
        action='store_true',
        help='If set, loads sigma factor-gene interactions into network.')
    parser.add_argument('--noise',
        nargs=3,
        type=int,
        metavar=('NOISE_ELEMENTS', 'REPEATS', 'SEED'),
        help='Number of noisy transcription factors to add, number of times to repeat with different genes and random seed for assigning genes.')

    # NCA options
    parser.add_argument('-m', '--method',
        choices=nca.METHODS,
        default=default_nca,
        help=f'NCA method to use, defined in nca.py (default: {default_nca}).')
    parser.add_argument('--status',
        type=float,
        default=default_status_update,
        help=f'Print status update every time this fraction of steps has been completed (default: {default_status_update}).')

    ## ISNCA specific options
    parser.add_argument('-i', '--iterative',
        action='store_true',
        help='If set, performs iterative sub-network component analysis.')
    parser.add_argument('--parallel',
        action='store_true',
        help='If set, runs sub-network NCA solvers in parallel.')
    parser.add_argument('--iterative-iterations',
        type=int,
        default=default_iterative_iterations,
        help=f'Maximum iterations for iterative sub-network component analysis (default: {default_iterative_iterations}).')
    parser.add_argument('--iterative-splits',
        type=int,
        default=default_iterative_split,
        help=f'Maximum splits for iterative sub-network component analysis (default: {default_iterative_split}).')

    ## ROBNCA specific options
    parser.add_argument('--robust-iterations',
        type=int,
        default=default_robust_iterations,
        help=f'Maximum iterations for robust_nca (default: {default_robust_iterations}).')

    return parser.parse_args()


if __name__ == '__main__':
    start = time.time()
    args = parse_args()

    # TODO: allow combinations of valid regulation types
    # TODO: add other types of interest ('Transcription-Attenuation', 'Small-Molecule-Mediated-Attenuation')
    if args.attenuation_only:
        valid_types = {'Ribosome-Mediated-Attenuation'}
    else:
        valid_types = {'Transcription-Factor-Binding'}

    # TODO: normalize seq data or only use EcoMAC data for now
    # TODO: handle options better when loading analysis or cache
    #   - seq_data can be different with args.linear
    #   - tf_genes can be different with args.split
    #   - tfs might not match saved regulation with args.global_expression
    print('Loading data from files...')
    seq_data, idx_mapping = load_seq_data(args.linear, args.average_seq)
    b_numbers, gene_symbols, synonyms = load_gene_names()
    tf_genes = load_tf_gene_interactions(
        valid_types, ecocyc=args.ecocyc, regulondb=args.regulondb,
        split=args.split, verbose=args.verbose,
        )

    if args.analysis is None:
        # Check input cached files
        output_dir = os.path.join(OUTPUT_DIR, args.label)
        os.makedirs(output_dir, exist_ok=True)
        cache_dir = args.cache if args.cache else os.path.join(output_dir, 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        network_cache_file = os.path.join(cache_dir, NETWORK_CACHE_FILE)
        tf_cache_file = os.path.join(cache_dir, TF_CACHE_FILE)
        no_cache = not (os.path.exists(network_cache_file) and os.path.exists(tf_cache_file))

        # Create or load network mapping and TF IDs
        if args.force or no_cache or args.iterative:
            print('Creating initial network mapping...')
            initial_tf_map, tfs = create_tf_map(b_numbers, synonyms, tf_genes, verbose=args.verbose)

            if not args.iterative:
                initial_tf_map, tfs = nca.nca_criteria_check(initial_tf_map, tfs, verbose=args.verbose)
        else:
            print('Loading cached initial network mapping...')
            initial_tf_map = np.load(network_cache_file)
            tfs = np.load(tf_cache_file)

        # Output cached files
        cache_dir = os.path.join(output_dir, 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        network_cache_file = os.path.join(cache_dir, NETWORK_CACHE_FILE)
        tf_cache_file = os.path.join(cache_dir, TF_CACHE_FILE)
        np.save(network_cache_file, initial_tf_map)
        np.save(tf_cache_file, tfs)

        if args.global_expression:
            tfs, initial_tf_map = add_global_expression(tfs, mapping=initial_tf_map)
        if args.sigma_factors:
            tfs, initial_tf_map = add_sigma_factors(tfs, b_numbers, synonyms, mapping=initial_tf_map, verbose=args.verbose)
        if args.noise:
            tfs, initial_tf_map = add_noisy_expression(tfs, *args.noise, mapping=initial_tf_map)  # type: ignore

        # Solve NCA problem
        nca_method = getattr(nca, args.method)
        if args.iterative:
            A, P, tfs = nca.iterative_sub_nca(nca_method, seq_data, initial_tf_map, tfs,
                statistics=match_statistics, statistics_args=(tf_genes, gene_symbols),
                n_iters=args.iterative_iterations, splits=args.iterative_splits, parallel=args.parallel,
                robust_iters=args.robust_iterations, status_step=args.status, verbose=args.verbose)
            np.save(tf_cache_file, tfs)
        else:
            A, P = nca_method(seq_data, initial_tf_map, n_iters=args.robust_iterations, status_step=args.status)

        # Save results
        save_regulation(A, P, gene_symbols, tfs, output_dir)
    else:
        print('Loading regulation results without running NCA...')
        output_dir = args.analysis
        if not os.path.exists(output_dir):
            raise IOError(f'Directory does not exist: {output_dir}')
        cache_dir = args.cache if args.cache else os.path.join(args.analysis, 'cache')

        tfs = np.load(os.path.join(cache_dir, TF_CACHE_FILE))
        if not args.iterative:
            if args.global_expression:
                tfs, _ = add_global_expression(tfs)
            if args.sigma_factors:
                tfs, _ = add_sigma_factors(tfs, gene_symbols, synonyms)
            if args.noise:
                tfs, _ = add_noisy_expression(tfs, *args.noise)

        A, P = load_regulation(args.analysis, gene_symbols, tfs)

    match_statistics(seq_data, A, P, tfs, tf_genes, gene_symbols)
    plot_results(tf_genes, A, P, gene_symbols, tfs, output_dir)
    save_fold_changes(valid_types, args.ecocyc, args.regulondb, A, P,
        gene_symbols, tfs, output_dir, args.attenuation_only)

    print(f'Completed in {(time.time() - start) / 60:.1f} min')
