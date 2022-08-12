#!/usr/bin/env python
"""
Compare Sander et al data to output from the model at different effective KI
values.

Requires data from remove_aa_inhibition variant and saving numpy arrays from the
variant analysis plot (original files will have underscores and start with
remove_aa_inhibition):
    control_conc: remove-inhib-conc-control.npy
    conc: remove-inhib-conc.npy
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
VALIDATION_DATA = os.path.join(FILE_LOCATION, 'sander-allosteric-aa-conc.tsv')
MODEL_DATA = os.path.join(FILE_LOCATION, 'remove-inhib-conc')
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# From variant analysis plot
AMINO_ACIDS = [
	'GLN',
	'GLT',
	'ARG',
	'PRO',
	'L-ASPARTATE',
	'ASN',
	'LYS',
	'MET',
	'THR',
	'VAL',
	'L-ALPHA-ALANINE',
	'SER',
	'GLY',
	'HIS',
	'PHE',
	'TRP',
	'TYR',
	'ILE',
	'LEU',
    ]
ENZYMES = [
	'argA',
	'trpE',
	'hisG',
	'leuA',
	'thrA',
	'ilvA',
	'proB',
    ]
AA_IDX = {aa: i for i, aa in enumerate(sorted(AMINO_ACIDS))}
ENZ_IDX = {enz: i for i, enz in enumerate(ENZYMES)}
KI_FACTORS = sorted([np.inf, 2, 5, 10, 100])[::-1]  # from sim variant file
KIS = np.array([0.15, 0.17, 0.07327, 0.28, 0.167, 0.15, 0.06])  # KIs (mM) from sim_data.process.metabolism.aa_kis corresponding to AA in COMPARISONS

# Data to plot
COMPARISONS = [
    ('ARG', 'argA'),
    ('TRP', 'trpE'),
    ('HIS', 'hisG'),
    ('ILE', 'ilvA'),
    ('LEU', 'leuA'),
    ('THR', 'thrA'),
    ('PRO', 'proB'),
    ]
BAR_AAS = [
    ['ARG'],
    ['TRP'],
    ['HIS'],
    ['ILE', 'LEU'],
    ['THR'],
    ['PRO'],
    ]

FIG_WIDTH = 7
FIG_HEIGHT = 2.5


def load_validation():
    data = {}
    with open(VALIDATION_DATA) as f:
        reader = csv.reader(f, delimiter='\t')

        for _ in range(3):
            next(reader)
        conditions = next(reader)[1:]
        for row in reader:
            aa = row[0]
            conc = np.array(row[1:], float)
            data[aa] = dict(zip(conditions, conc))

    return data

def load_model(appended=''):
    control = np.load(f'{MODEL_DATA}-control{appended}.npy')
    conc = np.load(f'{MODEL_DATA}{appended}.npy')
    return control, conc

def get_validation(validation, aa, enz):
    return validation[aa]['WT'], np.array([validation[aa][enz]])

def get_model(model, aa, enz):
    aa_idx = AA_IDX[aa]
    control = model[0][aa_idx]
    conc = model[1][AA_IDX[aa], ENZ_IDX[enz], :]
    return control, conc

def plot_ki_prediction(validation, model, show_stats=True):
    # TODO: decide on metric (need to do 1 - metric?)
    # TODO: fit curve to points and drop prediction vline
    cols = 4
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols, figsize=(FIG_WIDTH, FIG_HEIGHT), constrained_layout=True)
    hide_axes = np.ones_like(axes, dtype=bool)

    ki_factors = np.array(KI_FACTORS)
    inverse_factors = 1 / ki_factors
    predicted_kis = []
    for i, ((aa, enz), ki) in enumerate(zip(COMPARISONS, KIS)):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        val = get_validation(validation, aa, enz)
        wcm = get_model(model, aa, enz)

        val_increase = val[1][0] / val[0]
        wcm_increase = wcm[1] / wcm[0]

        # Original data
        ax.plot(inverse_factors, wcm_increase, 'o', alpha=0.5)

        # Interp function to calculate expected value based on validation increase
        x_interp = np.log(wcm_increase[1:])  # Fit in log space for smoother fit
        y_interp = np.log(inverse_factors[1:])
        interp = interp1d(x_interp, y_interp)
        val_match = np.exp(interp(np.log(val_increase)))
        predicted_ki = 1 / val_match * ki
        predicted_kis.append(predicted_ki)
        y_fit = np.linspace(x_interp.min(), x_interp.max(), 1000)
        ax.plot(np.exp(interp(y_fit)), np.exp(y_fit))

        ax.axhline(val_increase, linestyle='--', color='k', linewidth=0.5, alpha=0.5)
        ax.axvline(val_match, linestyle='--', color='k', linewidth=0.5, alpha=0.5)

        ax.set_xscale('symlog', linthreshx=0.01)
        ax.set_yscale('log')

        # Consistent x range across amino acids
        ax.set_xlim([-0.001, 1])

        # Get proper tick labels (more than one)
        if ax.get_ylim()[0] > 1:
            ax.set_ylim([1, ax.get_ylim()[1]])
        if ax.get_ylim()[1] < 10:
            ax.set_ylim([ax.get_ylim()[0], 10])

        if show_stats:
            xlabel = f'Original KI / adjusted KI\n(ki={predicted_ki:.2f}, val={val_match:.2f})'
            ax.set_xlabel(xlabel, fontsize=8, labelpad=2)
            if col == 0:
                ax.set_ylabel('Conc fold change', fontsize=8, labelpad=2)
        ax.set_title(aa, fontsize=8, pad=0)
        ax.tick_params(labelsize=8, pad=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

    return predicted_kis

def plot_kis(predicted_kis):
    labels = [c[1] for c in COMPARISONS]

    # Plot bar of predicted KI vs baseline KI
    plt.figure(figsize=(1.6, 1.4), constrained_layout=True)
    width = 0.4
    x = np.arange(len(KIS))
    plt.bar(x, KIS, width=-width, align='edge')
    plt.bar(x, predicted_kis, width=width, align='edge')
    plt.xticks(x, labels, rotation=45, fontsize=8)
    plt.ylabel('KI (mM)', fontsize=8, labelpad=0)

    ax = plt.gca()
    ax.tick_params(labelsize=8, pad=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_fig('kis.pdf')

def plot_ki_range(validation, model):
    _, axes = plt.subplots(len(COMPARISONS), 2, figsize=(FIG_WIDTH, 20), constrained_layout=True)
    for i, (aa, enz) in enumerate(COMPARISONS):
        val = get_validation(validation, aa, enz)
        wcm = get_model(model, aa, enz)

        ax = axes[i, 0]
        data = [val[1][0] / val[0]] + list(wcm[1] / wcm[0])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'Normalized {aa} conc to WT')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val'] + KI_FACTORS, rotation=45, fontsize=8)

        ax = axes[i, 1]
        data = [val[0]] + list(val[1]) + [wcm[0]] + list(wcm[1])
        x = np.arange(len(data))
        ax.bar(x, data)
        ax.set_yscale('log')
        ax.set_ylabel(f'{aa} conc')
        ax.set_xticks(x)
        ax.set_xticklabels(['Val WT', 'Val mutant', 'WCM WT'] + KI_FACTORS, rotation=45, fontsize=8)

def plot_bars(datasets, functions, log=False, normalize=False, bottom=None, single_first=True):
    # TODO: label x and y axes
    # TODO: highlight enzyme
    # TODO: option to normalize based on control
    cols = 3
    rows = int(np.ceil(len(BAR_AAS) / cols))
    _, axes = plt.subplots(rows, cols, figsize=(FIG_WIDTH / 1.4, FIG_HEIGHT), constrained_layout=True)
    hide_axes = np.ones_like(axes, dtype=bool)

    for i, aas in enumerate(BAR_AAS):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        all_conc = []
        for j, (data, fun) in enumerate(zip(datasets, functions)):
            for aa in aas:
                conc = np.array([fun(data, aa, ENZYMES[0])[0]] + [fun(data, aa, enz)[1][0] for enz in ENZYMES])
                if normalize:
                    conc = np.log10(conc / conc[0])
                all_conc.append(conc)

                if j == 0 and single_first:
                    break

        if len(all_conc) == 0:
            continue

        n_mutants = len(all_conc[0])
        n_datasets = len(all_conc)
        x = np.arange(n_mutants)
        width = 0.8 / n_datasets
        offsets = np.arange(n_datasets) * width - 0.4 + width/2

        for j, offset in enumerate(offsets):
            ax.bar(x + offset, all_conc[j], width, bottom=bottom)

        if log:
            ax.set_yscale('log')

        if normalize:
            ax.axhline(0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
        if i == 0:
            ax.set_ylabel('log10 increase\nover WT' if normalize else 'Concentration (mM)', fontsize=8, labelpad=2)

        ax.set_xticks(x)
        ax.set_xticklabels(['WT'] + ENZYMES, rotation=45, fontsize=8)
        ax.set_title(', '.join(aas), fontsize=8, pad=0)
        ax.tick_params(labelsize=8, pad=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

def plot_sub_scatter(validation, model):
    cols = 4
    rows = int(np.ceil(len(COMPARISONS) / cols))
    _, axes = plt.subplots(rows, cols, figsize=(2*FIG_WIDTH, 2*FIG_HEIGHT), constrained_layout=True)
    hide_axes = np.ones_like(axes, dtype=bool)

    for i, (aa, _) in enumerate(COMPARISONS):
        i += 1  # skip first ax
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        hide_axes[row, col] = False

        plot_scatter(ax, validation, model, label=aa, amino_acids=[aa], plot_all=True)

    # Hide unused axes
    for ax in axes[hide_axes]:
        ax.set_visible(False)

def plot_scatter(ax, validation, model, label='Amino acid', amino_acids=None,
        enzymes=None, legend=True, plot_all=False):
    if amino_acids is None:
        amino_acids = AMINO_ACIDS
    if enzymes is None:
        enzymes = ENZYMES

    val_control = []
    model_control = []
    val_mutants = []
    model_mutants = []
    val_other = []
    model_other = []

    for aa, enz in COMPARISONS:
        if aa not in amino_acids:
            continue

        val = get_validation(validation, aa, enz)
        val_control.append(val[0])
        val_mutants.append(val[1][0])

        mod = get_model(model, aa, enz)
        model_control.append(mod[0])
        model_mutants.append(mod[1][0])

        if plot_all:
            for m in mod[1][1:]:
                val_mutants.append(val[1][0])
                model_mutants.append(m)

    for aa in amino_acids:
        for enz in enzymes:
            if (aa, enz) in COMPARISONS:
                continue
            val_other.append(get_validation(validation, aa, enz)[1][0])
            model_other.append(get_model(model, aa, enz)[1][0])

    ax.loglog(val_control, model_control, 'or', alpha=0.5, label='Allosteric AA in WT')
    ax.loglog(val_mutants, model_mutants, 'ob', alpha=0.5, label='Allosteric AA in mutant')
    ax.loglog(val_other, model_other, 'ok', alpha=0.2, markersize=2, label='Other AA')

    if legend:
        ax.legend(fontsize=8, frameon=False)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_ax = min(xlim[0], ylim[0])
    max_ax = max(xlim[1], ylim[1])
    xy_line = [min_ax, max_ax]
    ax.loglog(xy_line, xy_line, '--k', alpha=0.2, linewidth=1)

    ax.set_xlabel(f'{label} conc\nin validation (mM)')
    ax.set_ylabel(f'{label} conc\nin model (mM)')

def save_fig(filename):
    file = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(file)
    plt.close('all')
    print(f'Saved to {file}')

def save_model_table(mean, var, filename='model.tsv'):
    file = os.path.join(OUTPUT_DIR, filename)
    with open(file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Amino acid', 'Wildtype'] + ENZYMES)

        for aa, enz in COMPARISONS:
            wt_mean = get_model(mean, aa, enz)[0]
            wt_std = np.sqrt(get_model(var, aa, enz)[0])
            mutant_data = [f'{get_model(mean, aa, mutant)[1][0]:.2g} +/- {np.sqrt(get_model(var, aa, mutant)[1][0]):.2g}' for mutant in ENZYMES]
            writer.writerow([aa, f'{wt_mean:.2g} +/- {wt_std:.2g}'] + mutant_data)

    print(f'Saved to {file}')

if __name__ == '__main__':
    validation = load_validation()
    model = load_model()
    model_var = load_model(appended='-var')

    # Compare validation data AA increase to level of inhibition reduction in the model
    predicted_kis = plot_ki_prediction(validation, model)
    save_fig('aa-ki-prediction.pdf')
    plot_ki_prediction(validation, model, show_stats=False)
    save_fig('aa-ki-prediction-clean.pdf')
    plot_kis(predicted_kis)

    # Compare validation data to range of KI values produced in the model
    plot_ki_range(validation, model)
    save_fig('aa-ki-range.pdf')

    # Validation bar plot to reproduce Fig 1B from paper
    plot_bars([validation], [get_validation])
    save_fig('validation-bar.pdf')

    # Comparable bar plot from model output
    plot_bars([model], [get_model], log=True, bottom=1e-2, single_first=False)
    save_fig('model-bar.pdf')

    # Side by side bar plot with validation and model data
    plot_bars([validation, model], [get_validation, get_model], log=True, bottom=1e-2)
    save_fig('side-by-side-bar.pdf')

    # Side by side bar plot with validation and model data
    plot_bars([validation, model], [get_validation, get_model], normalize=True)
    save_fig('normalized-side-by-side-bar.pdf')

    # Scatter plot between validation and model
    plt.figure()
    plot_scatter(plt.gca(), validation, model)
    plt.tight_layout()
    save_fig('scatter.pdf')

    # Scatter for each amino acid
    plot_sub_scatter(validation, model)
    save_fig('sub-scatter.pdf')

    # Table of model concentrations and variance
    save_model_table(model, model_var)
