#! /usr/bin/env python

"""
Shows output from a simple toy data example to compare different data and
modeling approaches. Plots outputs from NCA on a test data set with figure
comparison from https://www.eee.hku.hk/~cqchang/gNCA-fig.pdf for each method.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import nca


BASE_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
NCA_DATA_DIR = os.path.join(DATA_DIR, 'nca')
NCA_TEST_TOPOLOGY_FILE = os.path.join(NCA_DATA_DIR, 'subnet1_top.tsv')
NCA_TEST_EXPRESSION_FILE = os.path.join(NCA_DATA_DIR, 'subnet1_data.tsv')

OUTPUT_DIR = os.path.join(BASE_DIR, 'out', 'test-plots')
os.makedirs(OUTPUT_DIR, exist_ok=True)


# Randomly generated data to show the importance of a general expression column
# 4 conditions: wt, TF1 high, TF2 high, TF1 and TF2 low
# Total expression is roughly normalized in each condition
TEST_EXPRESSION = np.array([
    [8, 8.5, 7.7, 7.1],
    [8, 8.8, 7.5, 6.9],
    [8, 6.5, 9, 8],
    [8, 7.5, 7, 9],
    ])
TEST_TOPOLOGY = np.array([
    [1, 0, 1],
    [1, 0, 1],
    [-1, 1, 1],
    [0, -1, 1],
    ])


def chang_nca(method: str) -> None:
    """
    Plots TFAs to match results from Chang et al. Bioinformatics. 2008. Fig 5.

    Args:
        method: NCA method to use (function must be defined in nca.py)
    """

    def plot_output(P, start, samples, label, method):
        tfs = np.array([1, 8, 9, 14, 15, 21, 30, 32, 33, 34, 35])
        data = P[tfs, start:start+samples]
        n_tfs = data.shape[0]

        plt.figure(figsize=(5, 15))
        for i, tf in enumerate(data):
            ax = plt.subplot(n_tfs, 1, i+1)
            ax.plot(range(len(tf)), tf)

        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f'{method}_{label}.png'))
        plt.close('all')

    with open(NCA_TEST_EXPRESSION_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        expression = np.array(list(reader), float)

    with open(NCA_TEST_TOPOLOGY_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        topology = np.array(list(reader), float)

    A, P = getattr(nca, method)(expression, topology)

    plot_output(P, 0, 14, 'elutriation', method)
    plot_output(P, 14, 18, 'alpha', method)
    plot_output(P, 32, 24, 'cdc', method)
    plot_output(P, 56, 13, 'cycle', method)


if __name__ == '__main__':
    # TODO: create toy model to generate data to have exact fold change targets
    # Compare choice between linear, log2 and relative log2 input expression data
    # and adding a general expression "TF" vs not - log2 and general are best
    ## No general expression column with log2 expression input
    E = TEST_EXPRESSION
    A = TEST_TOPOLOGY[:, :-1]
    Ae, Pe = nca.robust_nca(E, A, verbose=False)
    with np.printoptions(precision=2, suppress=True):
        print('Results with log2 expression input without general expression')
        print('- high fold changes (in A)')
        print('- does not capture correct sign for TF2 (second column in A)')
        print('- TF activity does not match expectations (mostly second row in P)')
        print(f'A =\n{Ae}')
        print(f'P =\n{Pe}')
        print(f'error: {np.linalg.norm(2**TEST_EXPRESSION - 2**(Ae.dot(Pe))):.3f}')

    ## No general expression column with linear expression input
    E = 2**TEST_EXPRESSION
    A = TEST_TOPOLOGY[:, :-1]
    Ae, Pe = nca.robust_nca(E, A, verbose=False)
    with np.printoptions(precision=2, suppress=True):
        print('\nResults with linear expression input without general expression')
        print('- reasonable fold changes')
        print('- TF activity is consistent')
        print(f'A =\n{Ae}')
        print(f'P =\n{Pe}')
        print(f'error: {np.linalg.norm(2**TEST_EXPRESSION - Ae.dot(Pe)):.3f}')

    ## Added general expression column (all 1s) with linear expression input
    E = 2**TEST_EXPRESSION
    A = TEST_TOPOLOGY
    Ae, Pe = nca.robust_nca(E, A, verbose=False)
    with np.printoptions(precision=2, suppress=True):
        print('\nResults with linear expression input and general expression')
        print('- better fold changes')
        print('- TF activity is consistent')
        print('- better error')
        print(f'A =\n{Ae}')
        print(f'P =\n{Pe}')
        print(f'error: {np.linalg.norm(2**TEST_EXPRESSION - Ae.dot(Pe)):.3f}')

    ## Added general expression column (all 1s) with relative expression input
    ## Similar to Liao PNAS 2003 eq 7 implementation
    E = TEST_EXPRESSION - TEST_EXPRESSION.mean(1)
    A = TEST_TOPOLOGY
    Ae, Pe = nca.robust_nca(E, A, verbose=False)
    with np.printoptions(precision=2, suppress=True):
        print('\nResults with relative expression input and general expression:')
        print('- similar results as absolute log2 input')
        print(f'A =\n{Ae}')
        print(f'P =\n{Pe}')
        print(f'error: {np.linalg.norm(2**TEST_EXPRESSION - 2**(Ae.dot(Pe) + TEST_EXPRESSION.mean(1))):.3f}')

    ## Added general expression column (all 1s) with log2 expression input
    E = TEST_EXPRESSION
    A = TEST_TOPOLOGY
    Ae, Pe = nca.robust_nca(E, A, verbose=False)
    with np.printoptions(precision=2, suppress=True):
        print('\nResults with log2 expression input and general expression:')
        print('- best error and consistency')
        print(f'A =\n{Ae}')
        print(f'P =\n{Pe}')
        print(f'error: {np.linalg.norm(2**TEST_EXPRESSION - 2**(Ae.dot(Pe))):.3f}')

    # Test all methods on example data from Chang to compare output plots to paper
    print('\nTesting all methods on Chang data')
    for method in nca.METHODS:
        chang_nca(method)
