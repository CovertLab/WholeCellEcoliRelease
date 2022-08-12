#!/usr/bin/env python
"""
Plot the growth rate from fixed ppGpp runs with expression adjustments to
amino acid synthesis enzymes or ribosomes to show limitations.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


# Input paths
SIM_DIR = '/home/travis/scratch/wcEcoli_out/'
SENSITIVITY_DIR = 'ppGpp_sensitivity'
LIMITATIONS_LOW_DIR = 'ppGpp_limitations_-_low_ppGpp'
LIMITATIONS_HIGH_DIR = 'ppGpp_limitations_-_high_ppGpp'
RIBOSOME_LIMITATIONS_INHIBITION_DIR = 'ppGpp_limitations_with_ribosomes_at_high_ppGpp'
RIBOSOME_LIMITATIONS_NO_INHIBITION_DIR = 'ppGpp_limitations_with_ribosomes_at_high_ppGpp,_no_ppGpp_translation_inhibition'
FILE_PATH = 'plotOut/{}.tsv'

# Output paths
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cache paths
CACHE_DIR = os.path.join(FILE_LOCATION, 'cache')
os.makedirs(CACHE_DIR, exist_ok=True)
USE_CACHE = True
SAVE_CACHE = True

GROWTH_HEADER = 'Growth'
MEAN_HEADER = ' mean'


def load_data(desc, filename='2+gen-growth_trajectory', save_cache=SAVE_CACHE):
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

def plot_data(data, label=True):
    y = list(data.values())
    x = np.arange(len(y))
    if label:
        labels = list(data.keys())
        height = 2
    else:
        labels = [''] * len(x)
        height = 1

    plt.figure(figsize=(1 + 0.25 * len(x), height))
    plt.bar(x, y, color='k', alpha=0.8)

    plt.xticks(x, labels, rotation=30, ha='right', fontsize=8)
    plt.ylabel('Growth rate (1/hr)', fontsize=8)
    plt.ylim([0.4, 0.8])

    # Remove axes borders
    ax = plt.gca()
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def save_fig(output_file):
    plt.tight_layout()
    filename = os.path.join(OUTPUT_DIR, output_file + '.pdf')
    plt.savefig(filename)
    print(f'Saved to {filename}')
    plt.close('all')


if __name__ == '__main__':
    # Load data from disk
    sensitivity = load_data(SENSITIVITY_DIR)
    limitations_low = load_data(LIMITATIONS_LOW_DIR)
    limitations_high = load_data(LIMITATIONS_HIGH_DIR)
    ribosome_limit_inhibition = load_data(RIBOSOME_LIMITATIONS_INHIBITION_DIR)
    ribosome_limit_no_inhibition = load_data(RIBOSOME_LIMITATIONS_NO_INHIBITION_DIR)

    # Compile data to plot
    growth_key = GROWTH_HEADER + MEAN_HEADER
    low_data = {
        'Control': sensitivity[4][growth_key],  # Optimal conc (50 uM)
        'Low ppGpp': limitations_low[0][growth_key],  # No modifications (20 uM)
        'Increase enzymes': limitations_low[25][growth_key],  # 25% increase
        'Increase ribosomes': limitations_low[34][growth_key],  # 25% increase
        }
    high_data = {
        'Control': sensitivity[4][growth_key],  # Optimal conc (50 uM)
        'High ppGpp': limitations_high[74][growth_key],  # No modifications (90 uM)
        'Increase enzymes': limitations_high[99][growth_key],  # 25% increase
        'Increase ribosomes': ribosome_limit_inhibition[44][growth_key],  # 50% increase rProtein, 100% increase rRNA, 45 also works ok (10% rProtein, 50% rRNA)
        'No GTPase inhibition': ribosome_limit_no_inhibition[32][growth_key],  # No modifications, with no inhibition
        'Increase ribosomes,\nno GTPase inhibition': ribosome_limit_no_inhibition[44][growth_key],  # 50% increase rProtein, 100% increase rRNA, 45 also works ok (10% rProtein, 50% rRNA)
        }

    # Plot low ppGpp sim results
    plot_data(low_data)
    save_fig('low-ppgpp')
    plot_data(low_data, label=False)
    save_fig('low-ppgpp-clean')

    # Plot high ppGpp sim results
    plot_data(high_data)
    save_fig('high-ppgpp')
    plot_data(high_data, label=False)
    save_fig('high-ppgpp-clean')
