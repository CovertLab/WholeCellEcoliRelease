#!/usr/bin/env python
"""
Plot the growth rate vs R/P ratio in multiple conditions to show the new
environmental capabilities of the model.

Requires the directory of all the sims run for the paper as SIM_DIR with
certain descriptions assumed to select sets of runs.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# Input paths
SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
BASE_SIM_DIR = 'Conditions_without_regulation_or_charging'
NEW_C_SOURCE_DIR = 'Conditions_with_regulation'
ADD_ONE_DIR = 'Add_one_amino_acid_shift'
REMOVE_ONE_DIR = 'Remove_one_amino_acid_shift'
PPGPP_DIR = 'ppGpp_sensitivity'
PPGPP_LIMITATION_LOW_DIR = 'ppGpp_limitations_-_low_ppGpp'
PPGPP_LIMITATION_HIGH_DIR = 'ppGpp_limitations_-_high_ppGpp'
NEW_AA_SOURCE_DIR = 'Amino_acid_combinations_in_media'
INHIBITION_NO_PPGPP_DIR = 'Remove_amino_acid_inhibition_-_constant_ppgpp'
INHIBITION_DIR = 'Remove_amino_acid_inhibition'
FILE_PATH = 'plotOut/{}.tsv'

# Output paths
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_FILE = 'combined-growth-rp.pdf'

# Cache paths
CACHE_DIR = os.path.join(FILE_LOCATION, 'cache')
os.makedirs(CACHE_DIR, exist_ok=True)
USE_CACHE = True
SAVE_CACHE = True

CONTROL_IDX = 19
GROWTH_HEADER = 'Growth'
RP_HEADER = 'R/P ratio'
RPA_HEADER = 'R/(P+A) ratio'
MEAN_HEADER = ' mean'
STD_HEADER = ' std'

DENNIS_BREMER_2021 = np.array([
    [0.1691, 0.415888308335967],
    [0.2056, 0.693147180559945],
    [0.2576, 1.03972077083992],
    [0.3307, 1.38629436111989],
    [0.4176, 1.73286795139986],
    [0.5023, 2.07944154167984],
])

# Growth rate (1/hr) and RNA/protein from Zhu and Dai. Growth suppression by altered ppGpp. 2019. NAR.
ZHU_GROWTH = 0
ZHU_RP = 1
## Fig 3H SpoT OE
ZHU_LOW_PPGPP = np.array([
    [0.841111111111112, 0.241406177962194],
    [0.557777777777778, 0.330790686952513],
    [0.518888888888889, 0.34746887966805],
    [0.44, 0.378337943752882],
    [0.407777777777778, 0.407450437989857],
    [0.235555555555556, 0.484156293222684],
])
## Fig 2D RelA OE
ZHU_HIGH_PPGPP = np.array([
    [0.950099800399202, 0.275764143355069],
    [0.766467065868263, 0.228071469002285],
    [0.497005988023952, 0.17593320821044],
    [0.345309381237525, 0.149870408436854],
])

ONE_AA_OPTIONS = dict(alpha=0.5, markersize=4, markeredgewidth=0)
PPGPP_OPTIONS = dict(alpha=0.5, markersize=6, markeredgewidth=0)
FADE_OPTIONS = dict(alpha=0.2, markersize=4, color='black', markeredgewidth=0)

FIG_SIZE = (4, 4)

ORANGE = '#D55E00'
BLUE = '#0072B2'
GREEN = '#009E73'


def load_data(desc, filename='growth_trajectory', save_cache=SAVE_CACHE):
    filepath = FILE_PATH.format(filename)
    cached_filepath = os.path.join(CACHE_DIR, f'{desc}-{os.path.basename(filepath)}')

    if USE_CACHE and os.path.exists(cached_filepath):
        path = cached_filepath
        save_cache = False  # don't need to save cache if loaded from cache
    else:
        dirs = os.listdir(SIM_DIR)
        for d in dirs:
            if d.endswith(desc):
                path = os.path.join(SIM_DIR, d, filepath)
                break
        else:
            raise ValueError(f'{desc} not found in sim directory')

    data = {}
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')

        headers = next(reader)
        for row in reader:
            data[float(row[0])] = dict(zip(headers[1:], np.array(row[1:], float)))

    if save_cache:
        with open(path) as f:
            with open(cached_filepath, 'w') as fc:
                fc.write(f.read())

    return data

def plot(data, variants=None, exclude=None, std=True, label=None, options=None, x_header=RP_HEADER, y_header=GROWTH_HEADER):
    def extract_data(header, data_type):
        return np.array([data[v][header + data_type] for v in variants if v not in exclude])

    if variants is None:
        variants = list(data)
    if exclude is None:
        exclude = set()
    if options is None:
        options = {}

    rp_ratio = extract_data(x_header, MEAN_HEADER)
    growth = extract_data(y_header, MEAN_HEADER)

    if std:
        rp_ratio_std = extract_data(x_header, STD_HEADER)
        growth_std = extract_data(y_header, STD_HEADER)
    else:
        rp_ratio_std = None
        growth_std = None

    plt.errorbar(rp_ratio, growth, xerr=rp_ratio_std, yerr=growth_std, fmt='o', label=label, linewidth=1, **options)

    return rp_ratio, growth

def plot_conditions(fade=False, grouping=False, options=None, std=True, label=True):
    if fade:
        options = FADE_OPTIONS
        std = False
        label = False

    new_std = std

    original_options = dict(options) if options else {}
    parameterized_options = dict(options) if options else {}
    unparameterized_options = dict(options) if options else {}
    if grouping:
        # Colors from seaborn-colorblind prop cycle
        original_options.update(dict(color=ORANGE, markersize=10, alpha=0.5, markeredgewidth=0))
        parameterized_options.update(dict(color=BLUE, markersize=10, alpha=0.5, markeredgewidth=0))
        unparameterized_options.update(dict(color=GREEN, markersize=7, alpha=0.3, markeredgewidth=0))
        new_std = False
        label = False

    original_xy = plot(no_regulation, std=std, label='Original conditions' if label else '', variants=np.arange(3),
        options=original_options)
    new_xy = []
    new_xy.append(plot(regulation, std=std, label='New carbon sources with growth regulation' if label else '', variants=np.arange(3, 5),
        options=parameterized_options))
    new_xy.append(plot(add_one, std=std, variants=[CONTROL_IDX], label='Minimal + glc with growth regulation' if label else '',
        options=parameterized_options))
    new_xy.append(plot(remove_one, std=std, variants=[CONTROL_IDX], label='Rich + glc with growth regulation' if label else '',
        options=parameterized_options))
    new_xy.append(plot(new_aa, std=new_std, variants=[1, 2], label='New amino acid media with growth regulation' if label else '',
        options=unparameterized_options))
    new_xy.append(plot(add_one, std=False, exclude=[CONTROL_IDX], label='Add one AA to minimal with growth regulation' if label else '',
        options=unparameterized_options if unparameterized_options else ONE_AA_OPTIONS))
    new_xy.append(plot(remove_one, std=False, exclude=[CONTROL_IDX], label='Remove one AA from rich with growth regulation' if label else '',
        options=unparameterized_options if unparameterized_options else ONE_AA_OPTIONS))

    def print_r2(label, data):
        r, p = stats.pearsonr(*data)
        print(f'{label}: r2 = {r**2:.4f} ({p=:.2g})')

    print_r2('original', original_xy)
    print_r2('new', np.hstack(new_xy))

def plot_ppgpp():
    plot(add_one, variants=[CONTROL_IDX], label='Minimal + glc', options=dict(color='k'))
    plot(remove_one, variants=[CONTROL_IDX], label='Rich + glc', options=FADE_OPTIONS)
    plot(ppgpp, std=False, variants=range(10, 11), options=FADE_OPTIONS)
    plot(ppgpp, std=False, variants=range(12, 15), options=FADE_OPTIONS)  # 12, 20 for all
    plot(ppgpp, std=False, label='Minimal + glc,\nppGpp perturbed', variants=range(0, 4),
        options=dict(color=BLUE, marker='s', **PPGPP_OPTIONS))  # 0, 4 for all
    plot(ppgpp, std=False, variants=range(5, 10),
        options=dict(color=ORANGE, marker='s', **PPGPP_OPTIONS))  # 5, 10 for all

def plot_zhu():
    plt.plot(ZHU_LOW_PPGPP[:, ZHU_RP], ZHU_LOW_PPGPP[:, ZHU_GROWTH], 'X', label='Literature (Zhu et al.)', color=BLUE, **PPGPP_OPTIONS)
    plt.plot(ZHU_HIGH_PPGPP[:, ZHU_RP], ZHU_HIGH_PPGPP[:, ZHU_GROWTH], 'X', color=ORANGE, **PPGPP_OPTIONS)

def plot_inhibition():
    plot(inhib_no_ppgpp, variants=range(1, 8), std=False, label='Removed allosteric inhibition without ppGpp', options=PPGPP_OPTIONS)
    plot(inhib, variants=range(1, 8), std=False, label='Removed allosteric inhibition with ppGpp', options=PPGPP_OPTIONS)
    plot(add_one, variants=[CONTROL_IDX], label='Minimal + glc')
    plot(remove_one, variants=[CONTROL_IDX], label='Rich + glc')

def plot_trends(all=False):
    options = dict(linewidth=1, alpha=0.5)
    plt.plot([-0.3, 0.93], [-2, 4], '--k', **options)  # Dennis and Bremer. 1996. (dry_mass_composition.tsv)

    if all:
        plt.plot([-0.35, 0.91], [-2, 4], '--k', **options)  # Zhu et al. Growth suppression by altered (p)ppGpp levels... 2019.
        plt.plot(DENNIS_BREMER_2021[:, 0], DENNIS_BREMER_2021[:, 1], 'X', **options)

def format_plot(legend=True, ylim=(0, 2)):
    # Show legend
    if legend:
        plt.legend(fontsize=6, frameon=False)

    # Set axes
    plt.xlabel('RNA/protein mass ratio', fontsize=8)
    plt.ylabel('Growth rate (1/hr)', fontsize=8)
    plt.xlim([0, 0.6])
    plt.ylim(ylim)

    # Remove axes borders
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ax formatting
    ax.tick_params(labelsize=8)

    plt.tight_layout()

def save_fig(output_file):
    filename = os.path.join(OUTPUT_DIR, output_file)
    plt.savefig(filename)
    print(f'Saved to {filename}')
    plt.close('all')


if __name__ == '__main__':
    no_regulation = load_data(BASE_SIM_DIR)
    regulation = load_data(NEW_C_SOURCE_DIR)
    add_one = load_data(ADD_ONE_DIR)
    remove_one = load_data(REMOVE_ONE_DIR)
    ppgpp = load_data(PPGPP_DIR)
    ppgpp_low = load_data(PPGPP_LIMITATION_LOW_DIR)
    ppgpp_high = load_data(PPGPP_LIMITATION_HIGH_DIR)
    new_aa = load_data(NEW_AA_SOURCE_DIR)
    inhib_no_ppgpp = load_data(INHIBITION_NO_PPGPP_DIR)
    inhib = load_data(INHIBITION_DIR)

    plt.figure(figsize=FIG_SIZE)
    plot_conditions()
    plot_trends()
    format_plot()
    save_fig(OUTPUT_FILE)

    plt.figure(figsize=FIG_SIZE)
    plot_conditions(grouping=True)
    plot_trends()
    format_plot(legend=False)
    save_fig('groups-' + OUTPUT_FILE)

    plt.figure(figsize=(3.2, 3.2))
    plot_ppgpp()
    plot_zhu()
    plot_trends()
    format_plot()
    save_fig('ppgpp-' + OUTPUT_FILE)

    plt.figure(figsize=FIG_SIZE)
    plot_inhibition()
    plot_trends()
    format_plot()
    save_fig('inhibition-' + OUTPUT_FILE)

    # Could also plot regulation with variants 2 and 0
    plt.figure(figsize=(2.7, 2.7))
    plot(new_aa, variants=[0], options=dict(color=GREEN, markersize=4))
    plot(new_aa, variants=[3], options=dict(color=ORANGE, markersize=4))
    plot_trends()
    format_plot(legend=False, ylim=[-0.5, 3.5])
    save_fig('shifts-' + OUTPUT_FILE)
