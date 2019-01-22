#! /usr/bin/python

'''
Direct port of Bosdriesz from mathematica file
Uses ppGpp kinetics to model up/down shifts in nutrients

Working version of the mathematica script supplied in the supplement.
Naming convention mostly mirrors the mathematica script for easy comparison.
Shifts in conditions are simulated as a change in AA production rate and produce
similar results to their paper figure although not explicitly shown in their
mathematica file.
'''

from __future__ import division

import argparse
import csv
import multiprocessing as mp
import os
import time

from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import ode


file_location = os.path.dirname(os.path.realpath(__file__))

nAA = 20
cellVolume = 2.5e-15
nAvogadro = 6e23
proteinContent = 15e8

# indices for concentrations
aa_indices = range(nAA)
ta_indices = range(nAA, 2 * nAA)
ppgpp_index = ta_indices[-1] + 1
r_index = ppgpp_index + 1

np.random.seed(10)

# pars
params = {
	'e': 0.05,
	'kn': 0.2 * np.random.lognormal(np.log(1), 0.2, 20),
	'kIa': 100,
	'sTot': 1,
	'ks': 100,
	'kMaa': 100,
	'kMtf': 1,
	'krib': 20,
	'krta': 1,
	'krt': 500,
	'kRelA': 75,
	'RelAtot': 100 / (nAvogadro * cellVolume / 1e6),
	'kDRelA': 0.26,
	'vSpoTSynt': 0.001,
	'kSpoTdeg': np.log(2) / 30,
	'vInitMax': 2000,
	'rnapF': 1,
	'kMrrn': 20,
	'kIppGpp': 1,
	'nppGpp': 1,
	'gammamax': 1,
	'tau': 0.5,
	'f': 0.05,
	'bm': proteinContent / (cellVolume*nAvogadro/1e6),
	'nARib': 7459 * 1.65,
	'nAmet': 300,
	'rmax': proteinContent / (cellVolume*nAvogadro/1e6) / (7459*1.65),
	}

# Parameters for AA noise
aa_mean = 1
aa_std = 0.05


def dcdt_ode(t, c, args):
	'''
	Derivatives function that is called by ode from scipy.integrate

	Switches order of arguments to fit the expectation of ode but sharing
	dcdt code.

	Args:
		t (float): time of integration step
		c (array[float]): concentrations at integration step
		args (tuple): tuple to match the additional arguments in dcdt
			(params (dict), shift (float), single_shift (bool),
			shift_time (float), f_aa (float or array[float])

	Returns:
		array[float]: rates of change of each concentration
	'''

	return dcdt(c, t, *args)

def dcdt(c, t, params, shift=0, single_shift=False, shift_time=2000, f_aa=None):
	'''
	Derivatives function that is called by odeint from scipy.integrate

	Args:
		t (float): time of integration step
		c (array[float]): concentrations at integration step
		params (dict): dictionary of model parameters with keys as defined at the top of the file
		shift (float): indicator of nutrient shift direction
			0 (default): no shift
			positive: upshift
			negative: downshit
		single_shift (bool): whether to shift all amino acids (False, default)
			or a single AA (True)
		shift_time (float): time when nutrient shift (if any) occurs
		f_aa (float or array[float]): fraction of each amino acid present,
			if float it is constant for each, otherwise an array should contain
			a value for each amino acid

	Returns:
		array[float]: rates of change of each concentration
	'''

	dc = np.zeros_like(c)

	if f_aa is None:
		f_aa = params['f']

	# shift - not in mathematica file
	shift_magnitude = np.ones(nAA)
	if t > shift_time:
		# downshifts
		if shift < 0:
			if single_shift:
				shift_magnitude[-1] = 0.5  # single nutrient downshift
			else:
				shift_magnitude = 0.5*np.ones(nAA)  # all downshift - increase ppGpp, decrease ribosomes
		# upshifts
		elif shift > 0:
			if single_shift:
				shift_magnitude [-1] = 2  # single nutrient upshift
			else:
				shift_magnitude = 2*np.ones(nAA)  # all upshift - decrease ppGpp, increase ribosomes

	aa = c[aa_indices]
	taa = c[ta_indices]
	ppGpp = c[ppgpp_index]
	r = c[r_index]

	tf = params['tau'] * r - taa

	vAAsynt = shift_magnitude * params['bm'] * params['e'] * params['kn'] * (1 - r/params['rmax']) / (params['nAmet'] * (1 + aa / params['kIa']))
	vtRNAchar = params['ks'] * params['sTot'] * tf * aa / (params['kMaa'] * params['kMtf'] * (1 + tf / params['kMtf'] + aa / params['kMaa'] + tf * aa / params['kMaa'] / params['kMtf']))  # modified with `1 +`
	numeratorRibosome = 1 + np.sum(f_aa * (params['krta']/taa + tf/taa*params['krta']/params['krt']))
	vR = params['krib']*r / numeratorRibosome
	mu = vR / params['bm']
	vrrnInit = params['vInitMax'] * params['rnapF'] / (params['kMrrn'] + params['rnapF']) / (1 + (ppGpp / params['kIppGpp'])**params['nppGpp']) / (cellVolume * nAvogadro / 1e6)
	vribosome = vR
	vRsynt = min(vrrnInit, params['gammamax']*vR/params['nARib'])
	vRdilution = r * mu
	rtfSolutions = r * (f_aa * tf / taa * params['krta'] / params['krt']) / numeratorRibosome
	rtfTot = np.sum(rtfSolutions)
	vRelA = params['kRelA'] * params['RelAtot'] / (1 + params['kDRelA'] / rtfTot)
	vSpoTdeg = params['kSpoTdeg'] * ppGpp

	odesAA = vAAsynt - vtRNAchar
	odesTAA = vtRNAchar - f_aa*vribosome
	odesppGpp = vRelA + params['vSpoTSynt'] - vSpoTdeg
	odesMacromolComp = vRsynt - vRdilution

	# derivatives
	dc[aa_indices] = odesAA
	dc[ta_indices] = odesTAA
	dc[ppgpp_index] = odesppGpp
	dc[r_index] = odesMacromolComp

	return dc

def simulate(args, params, output_file):
	'''
	Simulate the ODE system and plot the results

	Args:
		args (argparse.Namespace): arguments parsed from the command line
		params (dict): dictionary of model parameters with keys as defined at the top of the file
		output_file (str): path to plot output

	Returns:
		array[float]: final values of ppGpp, ribosomes, elongation rate, AA, tRNA
	'''

	f_aa = params['f'] * np.ones(nAA)
	f_params = (params, args.shift, args.single_shift, args.shift_time, f_aa)
	f_all = [f_aa]

	# initial conditions
	co = np.zeros(2*nAA + 2)
	co[aa_indices] = params['kIa']  # aa (100)
	co[ta_indices] = 0.1 * params['tau'] * params['rmax']  # charged tRNA (4.063)
	co[ppgpp_index] = params['kIppGpp']  # ppGpp (1)
	co[r_index] = 0.2 * params['rmax']  # ribosome (16.25)
	tmax = args.tmax
	to = 0
	t = np.linspace(to,tmax,tmax)

	# solve ode with odeint (lsoda)
	if args.method == 'lsoda':
		sol = odeint(dcdt, co, t, args=f_params)
	# solve ode with ode (vode, forward or backward)
	else:
		solver = ode(dcdt_ode)
		solver.set_f_params(f_params)
		solver.set_integrator('vode', method=args.method, order=args.order)
		solver.set_initial_value(co, to)

		sol = [co]
		t = [to]
		while solver.successful() and solver.t < tmax:
			if args.noise:
				f_aa = np.random.normal(aa_mean, aa_std, nAA)
				f_aa /= f_aa.sum()
				f_params = (params, args.shift, args.single_shift, f_aa)
			solver.set_f_params(f_params)
			solver.integrate(solver.t + args.timestep)
			sol.append(solver.y)
			t.append(solver.t)
			f_all.append(f_aa)

	t = np.array(t)
	sol = np.array(sol)

	# solution timeseries
	aa = sol[:, aa_indices]
	taa = sol[:, ta_indices]
	ppgpp = sol[:, ppgpp_index]
	r = sol[:, r_index]

	# derived timeseries
	tf = params['tau'] * r.reshape(-1,1) - taa
	numeratorRibosome = 1 + np.sum(f_all * (params['krta'] / taa + tf / taa * params['krta'] / params['krt']), axis=1)
	vElongation = params['krib'] / numeratorRibosome

	# plot results
	n_subplots = 5
	plt.figure(figsize=(6, 9))
	plt.subplot(n_subplots, 1, 1)
	plt.plot(t, ppgpp)
	plt.ylabel('[ppGpp]')

	plt.subplot(n_subplots, 1, 2)
	plt.plot(t, r)
	plt.ylabel('[ribosomes]')

	plt.subplot(n_subplots, 1, 3)
	plt.plot(t, aa)
	plt.ylabel('[AA]')

	plt.subplot(n_subplots, 1, 4)
	plt.plot(t, taa)
	plt.ylabel('[charged tRNA]')

	plt.subplot(n_subplots, 1, 5)
	plt.plot(t, vElongation)
	plt.ylabel('Elongation Rate (AA/s)')

	plt.savefig(output_file)
	plt.close('all')

	return np.hstack((ppgpp[-1], r[-1], vElongation[-1], aa[-1, :].mean(), taa[-1, :].mean(), aa[-1, :], taa[-1, :]))

def update_params(params, key, value):
	'''
	Replaces a value in the params dictionary and returns the new dictionary.

	Args:
		params (dict): dictionary of model parameters with keys as defined at the top of the file
		key (str): key in params to update the value of
		value (float): value to update params to

	Returns:
		dict: copy of params input with the updated value for the specified key
	'''

	d = params.copy()
	d[key] = value
	return d

def parse():
	'''
	Parses arguments from the command line.

	Returns:
		argparse.Namespace: parsed arguments and values
	'''

	parser = argparse.ArgumentParser(description='Simulate growth control with ppGpp dynamics')

	parser.add_argument('-o', '--output', default='ppgpp',
		help='Output filename for plot and sensitivity analysis (default: ppgpp)')
	parser.add_argument('-s', '--shift', type=int, default=0,
		help='Shift direction (0 (default): no shift, positive: upshift, negative: downshift)')
	parser.add_argument('-m', '--method', default='adams',
		help='ODE solver method (adams (default, forward), bdf (backward), lsoda (adaptive)')
	parser.add_argument('--order', type=int, default=2,
		help='Order of solver method (default: 2), <= 12 for adams, <= 5 for bdf, not implemented for lsoda')
	parser.add_argument('-t', '--timestep', type=float, default=1,
		help='Timestep to advance for each integration (default: 1), not implemented for lsoda')
	parser.add_argument('--tmax', type=float, default=5000,
		help='Max time to integrate to (default: 5000)')
	parser.add_argument('--shift-time', type=float, default=2000,
		help='Time when shift occurs, only works with -s set (default: 2000)')
	parser.add_argument('--noise', action='store_true',
		help='Add noise to AA usage if set, not implemented for lsoda')
	parser.add_argument('--single-shift', action='store_true',
		help='Shift only one amino acid if set')
	parser.add_argument('--sensitivity', action='store_true',
		help='Perform sensitivity analysis if set, otherwise run one simulation')
	parser.add_argument('--parallel', action='store_true',
		help='Perform sensitivity analysis in parallel if set, only works with --sensitivity')

	return parser.parse_args()


if __name__ == '__main__':
	start = time.time()

	args = parse()

	# Analyze sensitivity to parameters
	if args.sensitivity:
		# Store output in sensitivity directory
		output_dir = os.path.join(file_location, 'sensitivity')
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)

		# Factors to vary parameters by
		variations = [0.1, 0.2, 0.5, 0.9, 1.1, 2, 5, 10]

		print('Running baseline')
		baseline = simulate(args, params, os.path.join(output_dir, args.output))

		# Perform sensitivity for each parameter in params
		with open(os.path.join(output_dir, '{}.tsv'.format(args.output)), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['Parameter', 'Factor', 'ppGpp', 'Ribosomes', 'Elongation Rate',
				'Average AA', 'Average tRNA'] + ['AA_{}'.format(i) for i in range(nAA)]
				+ ['tRNA_{}'.format(i) for i in range(nAA)])
			writer.writerow(['Baseline', ''] + list(baseline))

			for key, value in params.items():
				print('Running sensitivity for {}'.format(key))

				sim_args = [(args, update_params(params, key, value*factor),
					os.path.join(output_dir,'{}_{}_{}.png'.format(args.output, key, factor)))
					for factor in variations]

				if args.parallel:
					pool = mp.Pool(processes=mp.cpu_count())
					results = [pool.apply_async(simulate, sa) for sa in sim_args]
					pool.close()
					pool.join()

					for result, factor in zip(results, variations):
						if result.successful():
							sol = result.get()
							writer.writerow([key, factor] + list(sol / baseline))
						else:
							print('*** Error in multiprocessing for {} x{} ***'.format(key, factor))
							writer.writerow([key, factor])
				else:
					for sa, factor in zip(sim_args, variations):
						sol = simulate(*sa)
						writer.writerow([key, factor] + list(sol / baseline))
	# Run one simulation
	else:
		output_file = os.path.join(file_location, 'output', args.output)
		simulate(args, params, output_file)

	print('Completed in {:.1f} min'.format((time.time() - start) / 60))
