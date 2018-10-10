'''
Direct port of Bosdriesz from mathematica file
Uses ppGpp kinetics to model up/down shifts

Working version of the mathematica script supplied in the supplement.
Shifts in conditions are simulated as a change in AA production rate and produce
similar results to their paper figure although not explicitly shown in their
mathematica file.
'''

from __future__ import division

from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import numpy as np

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


def dcdt_ode(t, c):
	return dcdt(c, t)

def dcdt(c, t):
	dc = np.zeros_like(c)

	# shift - not in mathematica file
	shift = np.ones(nAA)
	if t > 2000:
		# downshifts
		shift = 0.5*np.ones(nAA)  # all downshift - increase ppGpp, decrease ribosomes
		# shift[-1] = 0.5  # single nutrient downshift

		# upshifts
		# shift = 2*np.ones(nAA)  # all upshift - decrease ppGpp, increase ribosomes
		# shift [-1] = 2  # single nutrient upshift

	aa = c[aa_indices]
	taa = c[ta_indices]
	ppGpp = c[ppgpp_index]
	r = c[r_index]

	tf = tau * r - taa

	vAAsynt = shift * bm * e * kn * (1 - r/rmax) / (nAmet * (1 + aa / kIa))
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


# initial conditions
co = np.zeros(2*nAA + 2)
co[aa_indices] = kIa  # aa (100)
co[ta_indices] = 0.1*tau*rmax  # charged tRNA (4.063)
co[ppgpp_index] = kIppGpp  # ppGpp (1)
co[r_index] = 0.2*rmax  # ribosome (16.25)
tmax = 5000
to = 0
t = np.linspace(to,tmax,tmax)

## solve ode with odeint (lsoda)
# sol = odeint(dcdt, co, t)

## solve ode with ode (vode, forward or backward)
sol = np.zeros((tmax, len(co)))
solver = ode(dcdt_ode)
solver.set_integrator('vode', method='adams', order=2)  # forward
# solver.set_integrator('vode', method='bdf', order=2)  # backward
solver.set_initial_value(co, to)

sol = [co]
t = [to]
dt = 1
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

plt.show()
