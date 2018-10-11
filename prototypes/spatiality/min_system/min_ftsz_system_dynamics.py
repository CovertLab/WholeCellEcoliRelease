"""
PDE simulation of MinE/MinD spatial dynamics.

The cytoplasm is modeled as a 2-dimensional field that contains three types of molecules: MinD-ADP, MinD-ATP, MinE.
The membrane is a 1-dimensional field that wraps from the midpoint of one cap to the midpoint of the other cap. It
contains MinD-ATP and MinE-MinD-ATP.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, print_function
# from __future__ import division
import numpy as np
from scipy.ndimage import convolve

import matplotlib

ANIMATE = True
SAVE_PLOT = True
INIT_IN_HALF = False

if ANIMATE:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

if ANIMATE:
	plt.ion()
	fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1)

## Simulation parameters
# Time
TIME_TOTAL = 50.0  # total time
DT = .0001  # time step
N_ITERATIONS = int(TIME_TOTAL / DT)  # number of iterations

# discretization of lattice
BIN_SIZE = 0.1  # micrometers

# laplacian kernel for diffusion
LAPLACIAN_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])

## Animation parameters
N_ANIMATE = 50
N_PLOT = 20
ANIMATION_STEP_SIZE = N_ITERATIONS // N_ANIMATE
PLOT_STEP_SIZE = N_ITERATIONS // N_PLOT

## Model parameters
# cell size
length = 10.0  # micrometers
radius = 0.5  # micrometers
diameter = 2 * radius

# chemical parameters
# capital letters for cytoplasmic molecules, lower case for membrane molecules
params = {
	'diffusion' : 2.5,  # micrometer**2/sec
	'k_ADP_ATP' : 1,  # sec**-1. Conversion rate of MinD:ADP to MinD:ATP
	'k_D' : 0.025,  # micrometer/sec. Spontaneous attachment rate of MinD:ATP to membrane
	'k_dD' : 0.0015,  # micrometer**3/sec. Recruitment of MinD:ATP to membrane by attached MinD
	'k_de' : 0.7*3,  # sec**-1. Rate of ATP hydrolysis in MinE:MinD:ATP complex, breaks apart complex and releases phosphate
	'k_E' : 0.093*3,  # micrometer**3/sec. Rate of MinE attachment to membrane-associated MinD:ATP complex
	'k_fE' : 0.01*3,  # micrometer**3/sec. Rate of MinE inhibition of FtsZ membrane-binding
	'k_fF' : 0.05,  # micrometer**3/sec. Rate of membrane-bound FtsZ recruitment of cytoplasmic FtsZ
	'k_F' : 0.0005,  # micrometer**3/sec. Spontaneous attachment rate of FtsZ to membrane
	}

# molecule indices
index = {
	'MinD-ADP[c]' : 0,
	'MinD-ATP[c]' : 1,
	'MinE[c]' : 2,
	'FtsZ[c]' : 3,
	'MinD-ATP[m]' : 0,
	'MinE-MinD-ATP[m]' : 1,
	'FtsZ[m]' : 2,
	}

# initial concentrations
MinD_initial_conc = 1000 / (np.pi * radius ** 2)  # reported: 1000/micrometer
MinE_initial_conc = 350 / (np.pi * radius ** 2)  # reported: 350/micrometer
FtsZ_initial_conc = 150 / (np.pi * radius ** 2)

# get number of bins and bin size
bins_x = int(length / BIN_SIZE)  # number of bins along x
bins_y = int(diameter / BIN_SIZE)  # number of bins along y
bins_mem = bins_x + bins_y - 2  # membrane goes along x and wraps around to each cap's midpoint. subtract 2 for the corners
dx = length / bins_x

# make templates
edges_template = np.pad(np.zeros((bins_y-2, bins_x-2)), (1, 1), 'constant', constant_values=(1, 1))

# indices for the membrane, specifying contact sites along the body of the cylinder, and the caps.
cap_pos1 = bins_y/2 - 1
cap_pos2 = bins_y/2 + bins_x - 2

top_membrane = [(i, 0) for i in range(1, bins_y/2)]
top_membrane.extend([(0, i) for i in range(bins_x)])
top_membrane.extend([(i, -1) for i in range(1, bins_y/2)])

bottom_membrane = [(i, 0) for i in range(bins_y/2, bins_y-1)]
bottom_membrane.extend([(-1, i) for i in range(bins_x)])
bottom_membrane.extend([(i, -1) for i in range(bins_y/2, bins_y-1)])


def run_step(cyto, mem):
	## Initialize deltas and diffusion arrays
	d_cyto = np.zeros_like(cyto)
	d_mem = np.zeros_like(mem)
	diffusion_c = np.zeros_like(cyto)

	## Save current state
	cyto_0 = np.copy(cyto)
	mem_0 = np.copy(mem)

	## Get diffusion
	diffusion_c[index['MinD-ADP[c]'], :, :] = convolve(
		cyto_0[index['MinD-ADP[c]'], :, :], LAPLACIAN_2D, mode='reflect'
		) / dx**2 * params['diffusion']
	diffusion_c[index['MinD-ATP[c]'], :, :] = convolve(
		cyto_0[index['MinD-ATP[c]'], :, :], LAPLACIAN_2D, mode='reflect'
		) / dx**2 * params['diffusion']
	diffusion_c[index['MinE[c]'], :, :] = convolve(
		cyto_0[index['MinE[c]'], :, :], LAPLACIAN_2D, mode='reflect'
		) / dx**2 * params['diffusion']
	diffusion_c[index['FtsZ[c]'], :, :] = convolve(
		cyto_0[index['FtsZ[c]'], :, :], LAPLACIAN_2D, mode='reflect'
		) / dx**2 * params['diffusion']

	## Get values at cytoplasm-membrane contact sites
	# cytoplasm to membrane contacts
	DT_c_top = np.array([cyto_0[index['MinD-ATP[c]'], :, :][idx] for idx in top_membrane])
	DT_c_bottom = np.array([cyto_0[index['MinD-ATP[c]'], :, :][idx] for idx in bottom_membrane])
	exchange_DT_c = (DT_c_top + DT_c_bottom) / 2

	E_c_top = np.array([cyto_0[index['MinE[c]'], :, :][idx] for idx in top_membrane])
	E_c_bottom = np.array([cyto_0[index['MinE[c]'], :, :][idx] for idx in bottom_membrane])
	exchange_E_c = (E_c_top + E_c_bottom) / 2

	FtsZ_c_top = np.array([cyto_0[index['FtsZ[c]'], :, :][idx] for idx in top_membrane])
	FtsZ_c_bottom = np.array([cyto_0[index['FtsZ[c]'], :, :][idx] for idx in bottom_membrane])
	exchange_FtsZ_c = (FtsZ_c_top + FtsZ_c_bottom) / 2

	# membrane to cytoplasm contacts
	exchange_DT_m = np.zeros((bins_y, bins_x))
	for mem_idx, cyto_idx in enumerate(top_membrane):
		exchange_DT_m[cyto_idx] = mem_0[index['MinD-ATP[m]'], :][mem_idx]
	for mem_idx, cyto_idx in enumerate(bottom_membrane):
		exchange_DT_m[cyto_idx] = mem_0[index['MinD-ATP[m]'], :][mem_idx]

	exchange_EDT_m = np.zeros((bins_y, bins_x))
	for mem_idx, cyto_idx in enumerate(top_membrane):
		exchange_EDT_m[cyto_idx] = mem_0[index['MinE-MinD-ATP[m]'], :][mem_idx]
	for mem_idx, cyto_idx in enumerate(bottom_membrane):
		exchange_EDT_m[cyto_idx] = mem_0[index['MinE-MinD-ATP[m]'], :][mem_idx]

	exchange_FtsZ_m = np.zeros((bins_y, bins_x))
	for mem_idx, cyto_idx in enumerate(top_membrane):
		exchange_FtsZ_m[cyto_idx] = mem_0[index['FtsZ[m]'], :][mem_idx]
	for mem_idx, cyto_idx in enumerate(bottom_membrane):
		exchange_FtsZ_m[cyto_idx] = mem_0[index['FtsZ[m]'], :][mem_idx]

	## Calculate reaction rates
	# rates for cytoplasm
	rxn_1_c = (params['k_D'] + params['k_dD'] * (exchange_DT_m + exchange_EDT_m)
		) * edges_template * cyto_0[index['MinD-ATP[c]'], :, :]
	rxn_2_c = params['k_E'] * exchange_DT_m * cyto_0[index['MinE[c]'], :, :]
	rxn_3_c = params['k_de'] * exchange_EDT_m
	rxn_4_c = params['k_ADP_ATP'] * cyto_0[index['MinD-ADP[c]'], :, :]
	rxn_5_c = params['k_fE'] * exchange_FtsZ_m * cyto_0[index['MinE[c]'], :, :]
	rxn_6_c = (params['k_F'] + params['k_fF'] * (exchange_FtsZ_m)
		) * edges_template * cyto_0[index['FtsZ[c]'], :, :]

	# rates for membrane
	rxn_1_m = (params['k_D'] + params['k_dD'] * (mem_0[index['MinD-ATP[m]'], :] + mem_0[index['MinE-MinD-ATP[m]'], :])
		) * exchange_DT_c
	rxn_2_m = params['k_E'] * mem_0[index['MinD-ATP[m]'], :] * exchange_E_c
	rxn_3_m = params['k_de'] * mem_0[index['MinE-MinD-ATP[m]'], :]
	rxn_5_m = params['k_fE'] * mem_0[index['FtsZ[m]'], :] * exchange_E_c
	rxn_6_m = (params['k_F'] + params['k_fF'] * (mem_0[index['FtsZ[m]'], :])
			   ) * exchange_FtsZ_c

	## Get derivatives
	d_cyto[index['MinD-ADP[c]'], :, :] = (diffusion_c[index['MinD-ADP[c]'], :, :] - rxn_4_c + rxn_3_c)
	d_cyto[index['MinD-ATP[c]'], :, :] = (diffusion_c[index['MinD-ATP[c]'], :, :] + rxn_4_c - rxn_1_c)
	d_cyto[index['MinE[c]'], :, :] = (diffusion_c[index['MinE[c]'], :, :] + rxn_3_c - rxn_2_c)
	d_cyto[index['FtsZ[c]'], :, :] = (diffusion_c[index['FtsZ[c]'], :, :] - rxn_6_c + rxn_5_c)

	d_mem[index['MinD-ATP[m]'], :] = (-rxn_2_m + rxn_1_m)
	d_mem[index['MinE-MinD-ATP[m]'], :] = (-rxn_3_m + rxn_2_m)
	d_mem[index['FtsZ[m]'], :] = (-rxn_5_m + rxn_6_m)

	return d_cyto, d_mem

def animate_step(DT_m, EDT_m, FtsZ_m, MinDD_c, MinDT_c, MinE_c, FtsZ_c, time):
	# clear non-imshow plots
	ax1.clear()
	ax2.clear()
	ax3.clear()

	## Plot
	# plot cytoplasm
	ax4.imshow(MinDD_c, cmap='YlGnBu', aspect="auto")
	ax5.imshow(MinDT_c, cmap='YlGnBu', aspect="auto")
	ax6.imshow(MinE_c, cmap='YlOrRd', aspect="auto")
	ax7.imshow(FtsZ_c, cmap='BuGn', aspect="auto")

	# plot membrane
	ax1.plot(DT_m.T, 'b')
	ax2.plot(EDT_m.T, 'r')
	ax3.plot(FtsZ_m.T, 'g')

	# add vertical lines for cap location
	ax1.axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
	ax1.axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')
	ax2.axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
	ax2.axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')
	ax3.axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
	ax3.axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')

	# turn off axes
	ax1.xaxis.set_visible(False)
	ax2.xaxis.set_visible(False)
	ax3.xaxis.set_visible(False)
	ax4.xaxis.set_visible(False)
	ax5.xaxis.set_visible(False)
	ax6.xaxis.set_visible(False)
	ax7.xaxis.set_visible(False)
	ax4.yaxis.set_visible(False)
	ax5.yaxis.set_visible(False)
	ax6.yaxis.set_visible(False)
	ax7.yaxis.set_visible(False)

	ax1.set_title('[MinD:ATP] membrane')
	ax2.set_title('[MinE:MinD:ATP] membrane')
	ax3.set_title('[FtsZ] membrane')
	ax4.set_title('[MinD:ADP] cytoplasm')
	ax5.set_title('[MinD:ATP] cytoplasm')
	ax6.set_title('[MinE] cytoplasm')
	ax7.set_title('[FtsZ] cytoplasm')

	plt.subplots_adjust(hspace=0.9)
	fig.suptitle('t = ' + str(time) + ' (s)', fontsize=16)
	plt.pause(0.00000001)

def save_plot(cyto_out, mem_out):
	fig2, axes = plt.subplots(N_PLOT, 8, figsize=(8,0.2*N_PLOT))

	for slice in xrange(N_PLOT):
		axes[slice, 0].text(0.5, 0.5, str(slice*PLOT_STEP_SIZE*DT)+'s')

		# plot cytoplasm
		axes[slice, 1].imshow(cyto_out[index['MinD-ATP[c]'],:, :, slice], cmap='YlGnBu', aspect="auto")
		axes[slice, 2].imshow(cyto_out[index['MinD-ADP[c]'],:, :, slice], cmap='YlGnBu', aspect="auto")
		axes[slice, 3].imshow(cyto_out[index['MinE[c]'],:, :, slice], cmap='YlOrRd', aspect="auto")
		axes[slice, 4].imshow(cyto_out[index['FtsZ[c]'], :, :, slice], cmap='BuGn', aspect="auto")

		# plot membrane
		axes[slice, 5].plot(mem_out[index['MinD-ATP[m]'],:, slice].T)
		axes[slice, 6].plot(mem_out[index['MinE-MinD-ATP[m]'],:, slice].T)
		axes[slice, 7].plot(mem_out[index['FtsZ[m]'], :, slice].T)

		# add vertical lines for cap location
		axes[slice, 5].axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
		axes[slice, 5].axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')
		axes[slice, 6].axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
		axes[slice, 6].axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')
		axes[slice, 7].axvline(x=cap_pos1, linestyle='--', linewidth=1, color='k')
		axes[slice, 7].axvline(x=cap_pos2, linestyle='--', linewidth=1, color='k')

		axes[slice, 0].axis('off')

		for x in range(1,8):
			axes[slice, x].set_xticks([])
			axes[slice, x].set_yticks([])

	axes[0, 1].set_title('[MinD:ATP] cytoplasm', fontsize=6)
	axes[0, 2].set_title('[MinD:ADP] cytoplasm', fontsize=6)
	axes[0, 3].set_title('[MinE] cytoplasm', fontsize=6)
	axes[0, 4].set_title('[FtsZ] cytoplasm', fontsize=6)
	axes[0, 5].set_title('[MinD:ATP] membrane', fontsize=6)
	axes[0, 6].set_title('[MinE:MinD:ATP] membrane', fontsize=6)
	axes[0, 7].set_title('[FtsZ] membrane', fontsize=6)
	plt.subplots_adjust(wspace=0.7, hspace=0.3)
	plt.savefig('prototypes/spatiality/min_system/output/min_ftsz_dynamics.pdf', bbox_inches='tight')


## Initialize molecular fields
# TODO (eran) make cytoplasm and membrane size automaticall from number of molecules in index.
# cytoplasm
cytoplasm = np.zeros((4, bins_y, bins_x))
cytoplasm[index['MinD-ADP[c]'], :, :] = np.random.normal(MinD_initial_conc, 1, (bins_y, bins_x))
cytoplasm[index['MinD-ATP[c]'], :, :] = np.random.normal(MinD_initial_conc, 1, (bins_y, bins_x))
cytoplasm[index['MinE[c]'], :, :] = np.random.normal(MinE_initial_conc, 1, (bins_y, bins_x))
cytoplasm[index['FtsZ[c]'], :, :] = np.random.normal(FtsZ_initial_conc, 1, (bins_y, bins_x))

if INIT_IN_HALF:
	# put all MinD in one half of the cytoplasm to speed up time to oscillations
	cytoplasm[index['MinD-ADP[c]'], :, 0:bins_x/2] *= 2.0
	cytoplasm[index['MinD-ADP[c]'], :, bins_x/2:bins_x] *= 0.0
	cytoplasm[index['MinD-ATP[c]'], :, 0:bins_x/2] *= 2.0
	cytoplasm[index['MinD-ATP[c]'], :, bins_x/2:bins_x] *= 0.0

# membrane
membrane = np.zeros((3, bins_mem))
membrane[index['MinD-ATP[m]'], :] = np.random.normal(MinD_initial_conc, 1, (bins_mem,))
membrane[index['MinE-MinD-ATP[m]'], :] = np.random.normal(MinE_initial_conc, 1, (bins_mem,))
membrane[index['FtsZ[m]'], :] = np.random.normal(FtsZ_initial_conc, 1, (bins_mem,))

## Initialize arrays for output, for use in save_plot
cytoplasm_out = np.empty((4, bins_y, bins_x, N_PLOT))
membrane_out = np.empty((3, bins_mem, N_PLOT))

# Simulate the PDEs with the finite difference method.
for i in range(N_ITERATIONS):

	d_cytoplasm, d_membrane = run_step(cytoplasm, membrane)
	cytoplasm += d_cytoplasm * DT
	membrane += d_membrane * DT

	# ## Disallow negatives
	# TODO (Eran) throw an exception if negative
	# cytoplasm[index['MinD-ADP[c]'], :, :][cytoplasm[index['MinD-ADP[c]'], :, :] < 0] = 0.0
	# cytoplasm[index['MinD-ATP[c]'], :, :][cytoplasm[index['MinD-ATP[c]'], :, :] < 0] = 0.0
	# cytoplasm[index['MinE[c]'], :, :][cytoplasm[index['MinE[c]'], :, :] < 0] = 0.0
	# membrane[index['MinD-ATP[m]'], :][membrane[index['MinD-ATP[m]'], :] < 0] = 0.0
	# membrane[index['MinE-MinD-ATP[m]'], :][membrane[index['MinE-MinD-ATP[m]'], :] < 0] = 0.0

	# Plot the state of the system.
	if ANIMATE and i % ANIMATION_STEP_SIZE == 0 and i < N_ANIMATE * ANIMATION_STEP_SIZE:
		animate_step(
			membrane[index['MinD-ATP[m]'], :],
			membrane[index['MinE-MinD-ATP[m]'], :],
			membrane[index['FtsZ[m]'], :],
			cytoplasm[index['MinD-ADP[c]'], :, :],
			cytoplasm[index['MinD-ATP[c]'], :, :],
			cytoplasm[index['MinE[c]'], :, :],
			cytoplasm[index['FtsZ[c]'], :, :],
			i * DT
			)

	if SAVE_PLOT and i % PLOT_STEP_SIZE == 0 and i < N_PLOT * PLOT_STEP_SIZE:
		cytoplasm_out[:, :, :, i / PLOT_STEP_SIZE] = cytoplasm[:, :, :]
		membrane_out[:, :, i / PLOT_STEP_SIZE] = membrane[:, :]

if SAVE_PLOT:
	save_plot(cytoplasm_out, membrane_out)
