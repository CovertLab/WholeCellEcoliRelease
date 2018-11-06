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
import os

from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.integrate import ode

file_location = os.path.dirname(os.path.realpath(__file__))

nAA = 20
tmax = 10000
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
e = 0.05
kn = 0.2 * np.random.lognormal(np.log(1), 0.2, 20)
kIa = 100
sTot = 1
ks = 100
kMaa = 100
kMtf = 1
krib = 20
krta = 1
krt = 500
kRelA = 75
RelAtot = 100 / (nAvogadro * cellVolume / 1e6)
kDRelA = 0.26
vSpoTSynt = 0.001
kSpoTdeg = np.log(2) / 30
vInitMax = 2000
rnapF = 1
kMrrn = 20
kIppGpp = 1
nppGpp = 1
gammamax = 1
tau = 0.5
f = 0.05
bm = proteinContent / (cellVolume*nAvogadro/1e6)
nARib = 7459 * 1.65
nAmet = 300
rmax = proteinContent / (cellVolume*nAvogadro/1e6) / (7459*1.65)


def dcdt_ode(t, c, args):
	'''
	Derivatives function that is called by ode from scipy.integrate

	Switches order of arguments to fit the expectation of ode but sharing
	dcdt code.

	Args:
		t (float): time of integration step
		c (array[floats]): concentrations at integration step
		args (tuple): tuple to match the additional arguments in dcdt
			(shift (float), single_shift (bool))

	Returns:
		array[floats]: rates of change of each concentration
	'''

	return dcdt(c, t, *args)

def dcdt(c, t, shift=0, single_shift=False):
	'''
	Derivatives function that is called by odeint from scipy.integrate

	Args:
		t (float): time of integration step
		c (array[floats]): concentrations at integration step
		shift (float): indicator of nutrient shift direction
			0 (default): no shift
			positive: upshift
			negative: downshit
		single_shift (bool): whether to shift all amino acids (False, default)
			or a single AA (True)

	Returns:
		array[floats]: rates of change of each concentration
	'''

	dc = np.zeros_like(c)

	# shift - not in mathematica file
	shift_magnitude = np.ones(nAA)
	if t > 2000:
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

	tf = tau * r - taa

	vAAsynt = shift_magnitude * bm * e * kn * (1 - r/rmax) / (nAmet * (1 + aa / kIa))
	vtRNAchar = ks * sTot * tf * aa / (kMaa * kMtf * (1 + tf / kMtf + aa / kMaa + tf * aa / kMaa / kMtf))  # modified with `1 +`
	numeratorRibosome = 1 + np.sum(f * (krta/taa + tf/taa*krta/krt))
	vR = krib*r / numeratorRibosome
	mu = vR / bm
	vrrnInit = vInitMax*rnapF/(kMrrn + rnapF) / (1 + ppGpp / kIppGpp * nppGpp) / (cellVolume*nAvogadro/1e6)
	vribosome = vR
	vRsynt = min(vrrnInit, gammamax*vR/nARib)
	vRdilution = r * mu
	rtfSolutions = r*(f*tf/taa*krta/krt)/numeratorRibosome
	rtfTot = np.sum(rtfSolutions)
	vRelA = kRelA * RelAtot / (1 + kDRelA/rtfTot)
	vSpoTdeg = kSpoTdeg * ppGpp

	odesAA = vAAsynt - vtRNAchar
	odesTAA = vtRNAchar - f*vribosome
	odesppGpp = vRelA + vSpoTSynt - vSpoTdeg
	odesMacromolComp = vRsynt - vRdilution

	# derivatives
	dc[aa_indices] = odesAA
	dc[ta_indices] = odesTAA
	dc[ppgpp_index] = odesppGpp
	dc[r_index] = odesMacromolComp

	return dc

def simulate(args):
	'''
	Simulate the ODE system and plot the results

	Args:
		args: arguments parsed from the command line
	'''

	# Handle arguments from argparse
	shift = args.shift
	single_shift = args.single_shift
	method = args.method  # adams (forward), bdf (backward), lsoda
	order = args.order
	dt = args.timestep
	output_file = os.path.join(file_location, args.output)
	f_params = (shift, single_shift)

	# initial conditions
	co = np.zeros(2*nAA + 2)
	co[aa_indices] = kIa  # aa (100)
	co[ta_indices] = 0.1*tau*rmax  # charged tRNA (4.063)
	co[ppgpp_index] = kIppGpp  # ppGpp (1)
	co[r_index] = 0.2*rmax  # ribosome (16.25)
	tmax = 5000
	to = 0
	t = np.linspace(to,tmax,tmax)

	# solve ode with odeint (lsoda)
	if method == 'lsoda':
		sol = odeint(dcdt, co, t, args=f_params)
	# solve ode with ode (vode, forward or backward)
	else:
		sol = np.zeros((tmax, len(co)))
		solver = ode(dcdt_ode)
		solver.set_f_params(f_params)
		solver.set_integrator('vode', method=method, order=order)
		solver.set_initial_value(co, to)

		sol = [co]
		t = [to]
		while solver.successful() and solver.t < tmax:
			solver.integrate(solver.t + dt)
			sol.append(solver.y)
			t.append(solver.t)

	t = np.array(t)
	sol = np.array(sol)

	# solution timeseries
	aa = sol[:,aa_indices]
	taa = sol[:,ta_indices]
	ppgpp = sol[:,ppgpp_index]
	r = sol[:,r_index]

	# derived timeseries
	tf = tau * r.reshape(-1,1) - taa
	numeratorRibosome = 1 + np.sum(f * (krta / taa + tf / taa * krta / krt), axis=1)
	vElongation = krib / numeratorRibosome

	# plot results
	n_subplots = 5
	plt.figure(figsize=(6,9))
	plt.subplot(n_subplots,1,1)
	plt.plot(t, ppgpp)
	plt.ylabel('[ppGpp]')

	plt.subplot(n_subplots,1,2)
	plt.plot(t, r)
	plt.ylabel('[ribosomes]')

	plt.subplot(n_subplots,1,3)
	plt.plot(t, aa)
	plt.ylabel('[AA]')

	plt.subplot(n_subplots,1,4)
	plt.plot(t, taa)
	plt.ylabel('[charged tRNA]')

	plt.subplot(n_subplots,1,5)
	plt.plot(t, vElongation)
	plt.ylabel('Elongation Rate (AA/s)')

	plt.savefig(output_file)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Simulate growth control with ppGpp dynamics')

	parser.add_argument('-o', '--output', default='ppgpp',
		help='Output filename for plot (default: ppgpp)')
	parser.add_argument('-s', '--shift', type=int, default=0,
		help='Shift direction (0 (default): no shift, positive: upshift, negative: downshift)')
	parser.add_argument('-m', '--method', default='adams',
		help='ODE solver method (adams (default, forward), bdf (backward), lsoda (adaptive)')
	parser.add_argument('--order', type=int, default=2,
		help='Order of solver method (default: 2), <= 12 for adams, <= 5 for bdf, not implemented for lsoda')
	parser.add_argument('-t', '--timestep', type=float, default=1,
		help='Timestep to advance for each integration (default: 1), not implemented for lsoda')
	parser.add_argument('--single-shift', action='store_true',
		help='Shift only one amino acid if set')

	args = parser.parse_args()

	simulate(args)