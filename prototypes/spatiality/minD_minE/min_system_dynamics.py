"""
PDE simulation of MinE/MinD spatial dynamics.

The cytoplasm is modeled as a 2-dimensional field that contains three types of molecules: MinD-ADP, MinD-ATP, MinE.
The membrane is a 1-dimensional field that wraps from the midpoint of one cap to the midpoint of the other cap. It
contains MinD-ATP and MinE-MinD-ATP.

based on:
	Huang, K.C., Meir, Y., & Wingreen, N.S. (2003).
	Dynamic structures in Escherichia coli: spontaneous formation of MinE rings and MinD polar zones.
	Proceedings of the National Academy of Sciences, 100(22), 12724-12728.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
# from __future__ import division
# from __future__ import print_function
import numpy as np
from scipy.ndimage import convolve

import matplotlib

ANIMATE = True
SAVE_PLOT = True

if ANIMATE:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

if ANIMATE:
	plt.ion()
	fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1)

## Simulation parameters
# cell size
length = 4.0  # micrometers
radius = 0.5  # micrometers
diameter = 2 * radius

# discretization of lattice
bin_size = 0.05  # micrometers

# time
T = 50.0  # total time
dt = .0001  # time step
n = int(T / dt)  # number of iterations

## Animation parameters
n_animate = 50
n_plot = 20
animate_step = n // n_animate
plot_step = n // n_plot

# get number of bins and bin size
bins_x = int(length / bin_size)  # number of bins along x
bins_y = int(diameter / bin_size)  # number of bins along y
bins_mem = bins_x + bins_y - 2  # membrane goes along x and wraps around to each cap's midpoint. subtract 2 for the corners
dx = length / bins_x

# make templates
edges_template = np.pad(np.zeros((bins_y-2, bins_x-2)), (1,1), 'constant', constant_values=(1,1))

# indices for the membrane, specifying contact sites along the body of the cylinder, and the caps.
edge_length = bins_y/2 - 1 # TODO (Eran) -- get rid of this

top_membrane = [(i,0) for i in range(1,bins_y/2)]
top_membrane.extend([(0,i) for i in range(bins_x)])
top_membrane.extend([(i,-1) for i in range(1,bins_y/2)])

bottom_membrane = [(i,0) for i in range(bins_y/2,bins_y-1)]
bottom_membrane.extend([(-1,i) for i in range(bins_x)])
bottom_membrane.extend([(i,-1) for i in range(bins_y/2,bins_y-1)])

# laplacian kernels for diffusion
laplace_kernel_2D = np.array([[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]])
# laplace_kernel_2D = np.array([[0.5, 1.0, 0.5], [1.0, -6., 1.0], [0.5, 1.0, 0.5]])
laplace_kernel_1D = np.array([1.0, -2.0, 1.0])

# chemical parameters
diffusion = 2.5  # micrometer**2/sec

k_ADP_ATP = 1  # sec**-1
k_D = 0.025  # micrometer/sec
k_dD = 0.0015  # micrometer**3/sec
k_de = 0.7  # sec**-1
k_E = 0.093  # micrometer**3/sec

## Initialize molecular fields
MinD_conc = 110 / (np.pi * radius ** 2)  # reported: 1000/micrometer
MinE_conc = 100 / (np.pi * radius ** 2)  # reported: 350/micrometer
# cytoplasm
DD_c = abs(np.random.normal(MinD_conc, 50, (bins_y, bins_x)))
DT_c = abs(np.random.normal(MinD_conc, 50, (bins_y, bins_x)))
E_c = abs(np.random.normal(MinE_conc, 50, (bins_y, bins_x)))
# membrane
DT_m = abs(np.random.normal(0, 50, (bins_mem,)))
EDT_m = abs(np.random.normal(0, 50, (bins_mem,)))

## Initialize arrays for plotting
# cytoplasm
DD_c_out = np.empty((bins_y, bins_x, n_plot))
DT_c_out = np.empty((bins_y, bins_x, n_plot))
E_c_out = np.empty((bins_y, bins_x, n_plot))
# membrane
DT_m_out = np.empty((bins_mem, n_plot))
EDT_m_out = np.empty((bins_mem, n_plot))

save_concentrations = np.zeros((5,n_animate))

def show_patterns(DT_m, EDT_m, MinDD_c, MinDT_c, MinE_c, concentrations, time):
	# clear non-imshow plots
	ax1.clear()
	ax2.clear()
	ax6.clear()

	# plot
	ax1.plot(DT_m.T, 'b')
	# ax1.plot(EDT_m.T, 'r')
	ax2.plot(EDT_m.T, 'r')
	ax1.axvline(x=edge_length, linestyle='--', linewidth=1, color='k')
	ax1.axvline(x=edge_length + bins_x - 2, linestyle='--', linewidth=1, color='k')
	ax2.axvline(x=edge_length, linestyle='--', linewidth=1, color='k')
	ax2.axvline(x=edge_length + bins_x - 2, linestyle='--', linewidth=1, color='k')
	ax3.imshow(MinDD_c, cmap='YlGnBu', aspect="auto")
	ax4.imshow(MinDT_c, cmap='YlGnBu', aspect="auto")
	ax5.imshow(MinE_c, cmap='YlOrRd', aspect="auto")
	ax6.plot(concentrations.T)

	# turn off axes
	ax1.xaxis.set_visible(False)
	ax2.xaxis.set_visible(False)
	ax3.xaxis.set_visible(False)
	ax4.xaxis.set_visible(False)
	ax5.xaxis.set_visible(False)
	ax3.yaxis.set_visible(False)
	ax4.yaxis.set_visible(False)
	ax5.yaxis.set_visible(False)

	ax1.set_title('[MinD:ATP] membrane')
	ax2.set_title('[MinE:MinD:ATP] membrane')
	ax3.set_title('[MinD:ADP] cytoplasm')
	ax4.set_title('[MinD:ATP] cytoplasm')
	ax5.set_title('[MinE] cytoplasm')
	ax6.set_title('avg concentrations')

	plt.subplots_adjust(hspace=0.9)
	fig.suptitle('t = ' + str(time) + ' (s)', fontsize=16)
	plt.pause(0.00000001)


# Simulate the PDEs with the finite difference method.
for i in range(n):

	## Save current values
	DD_c_0 = np.copy(DD_c)
	DT_c_0 = np.copy(DT_c)
	E_c_0 = np.copy(E_c)
	DT_m_0 = np.copy(DT_m)
	EDT_m_0 = np.copy(EDT_m)

	## Diffusion
	d_DD_c = convolve(DD_c, laplace_kernel_2D, mode='reflect') / dx**2 * diffusion
	d_DT_c = convolve(DT_c, laplace_kernel_2D, mode='reflect') / dx**2 * diffusion
	d_E_c = convolve(E_c, laplace_kernel_2D, mode='reflect') / dx**2 * diffusion
	d_ET_m = convolve(DT_m, laplace_kernel_1D, mode='reflect') / dx**2 * diffusion #**0.5
	d_EDT_m = convolve(EDT_m, laplace_kernel_1D, mode='reflect') / dx**2 * diffusion #**0.5

	## Get values at cytoplasm-membrane contact sites
	# cytoplasm to membrane contacts
	E_c_top = np.array([E_c_0[idx] for idx in top_membrane])
	E_c_bottom = np.array([E_c_0[idx] for idx in bottom_membrane])
	E_c_exchange = (E_c_top + E_c_bottom)/2

	DT_c_top = np.array([DT_c_0[idx] for idx in top_membrane])
	DT_c_bottom = np.array([DT_c_0[idx] for idx in bottom_membrane])
	DT_c_exchange = (DT_c_top + DT_c_bottom)/2

	# membrane to cytoplasm contacts
	DT_m_exchange = np.zeros((bins_y, bins_x))
	for mem_idx, cyto_idx in enumerate(top_membrane):
		DT_m_exchange[cyto_idx] = DT_m_0[mem_idx]
	for mem_idx, cyto_idx in enumerate(bottom_membrane):
		DT_m_exchange[cyto_idx] = DT_m_0[mem_idx]

	EDT_m_exchange = np.zeros((bins_y, bins_x))
	for mem_idx, cyto_idx in enumerate(top_membrane):
		EDT_m_exchange[cyto_idx] = EDT_m_0[mem_idx]
	for mem_idx, cyto_idx in enumerate(bottom_membrane):
		EDT_m_exchange[cyto_idx] = EDT_m_0[mem_idx]

	## Calculate reaction rates
	# get rates for cytoplasm
	rxn_1_c = (k_D + k_dD * (DT_m_exchange + EDT_m_exchange)) * edges_template * DT_c_0
	rxn_2_c = k_E * DT_m_exchange * E_c_0
	rxn_3_c = k_de * EDT_m_exchange
	rxn_4_c = k_ADP_ATP * DD_c_0

	# get rates for membrane
	rxn_1_m = (k_D + k_dD * (DT_m_0 + EDT_m_0)) * DT_c_exchange
	rxn_2_m = k_E * DT_m_0 * E_c_exchange
	rxn_3_m = k_de * EDT_m_0

	## Update the variables
	DD_c = DD_c_0 + (d_DD_c - rxn_4_c + rxn_3_c) * dt
	DT_c = DT_c_0 + (d_DT_c + rxn_4_c - rxn_1_c) * dt
	E_c = E_c_0 + (d_E_c + rxn_3_c - rxn_2_c) * dt
	DT_m = DT_m_0 + (d_ET_m - rxn_2_m + rxn_1_m) * dt
	EDT_m = EDT_m_0 + (d_EDT_m - rxn_3_m + rxn_2_m) * dt

	## Disallow negatives
	DD_c[DD_c < 0] = 0.0
	DT_c[DT_c < 0] = 0.0
	E_c[E_c < 0] = 0.0
	DT_m[DT_m < 0] = 0.0
	EDT_m[EDT_m < 0] = 0.0

	# Plot the state of the system.
	if ANIMATE and i % animate_step == 0 and i < n_animate * animate_step:
		# check reaction balance
		diff_rxn_1 = np.sum(rxn_1_c) - 2 * np.sum(rxn_1_m)
		diff_rxn_2 = np.sum(rxn_2_c) - 2 * np.sum(rxn_2_m)
		diff_rxn_3 = np.sum(rxn_3_c) - 2 * np.sum(rxn_3_m)

		# mean concentrations
		mean_DD_c = np.mean(DD_c)
		mean_DT_c = np.mean(DT_c)
		mean_E_c = np.mean(E_c)
		mean_DT_m = np.mean(DT_m)
		mean_EDT_m = np.mean(EDT_m)

		print "diff_rxn_1 = {:.2f} | diff_rxn_2 = {:.2f} | diff_rxn_3 = {:.2f}".format(diff_rxn_1,diff_rxn_2,diff_rxn_3)

		save_concentrations[:,i/animate_step] = [mean_DD_c,mean_DT_c,mean_E_c,mean_DT_m,mean_EDT_m]
		show_patterns(DT_m, EDT_m, DD_c, DT_c, E_c, save_concentrations, i * dt)


	if SAVE_PLOT and i % plot_step == 0 and i < n_plot * plot_step:
		## Save state for plot
		# cytoplasm
		DD_c_out[:,:,i/plot_step] = DD_c
		DT_c_out[:,:,i/plot_step] = DT_c
		E_c_out[:,:,i/plot_step] = E_c
		# membrane
		DT_m_out[:,i/plot_step] = DT_m
		EDT_m_out[:,i/plot_step] = EDT_m


if SAVE_PLOT:
	fig, axes = plt.subplots(n_plot, 6, figsize=(8,0.5*n_plot))

	for slice in xrange(n_plot):
		axes[slice, 0].text(0.5, 0.5, str(slice*plot_step*dt)+'s')

		axes[slice, 1].imshow(DT_c_out[:,:,slice], cmap='YlGnBu', aspect="auto")
		axes[slice, 2].imshow(DD_c_out[:, :, slice], cmap='YlGnBu', aspect="auto")
		axes[slice, 3].imshow(E_c_out[:,:,slice], cmap='YlOrRd', aspect="auto")
		axes[slice, 4].plot(DT_m_out[:,slice].T)
		axes[slice, 5].plot(EDT_m_out[:,slice].T)

		axes[slice, 4].axvline(x=edge_length, linestyle='--', linewidth=1, color='k')
		axes[slice, 4].axvline(x=edge_length+bins_x-2, linestyle='--', linewidth=1, color='k')
		axes[slice, 5].axvline(x=edge_length, linestyle='--', linewidth=1, color='k')
		axes[slice, 5].axvline(x=edge_length+bins_x-2, linestyle='--', linewidth=1, color='k')

		axes[slice, 0].axis('off')
		axes[slice, 1].set_xticks([])
		axes[slice, 2].set_xticks([])
		axes[slice, 3].set_xticks([])
		axes[slice, 4].set_xticks([])
		axes[slice, 5].set_xticks([])
		axes[slice, 1].set_yticks([])
		axes[slice, 2].set_yticks([])
		axes[slice, 3].set_yticks([])
		axes[slice, 4].set_yticks([])
		axes[slice, 5].set_yticks([])

	axes[0, 1].set_title('[MinD:ATP] cytoplasm', fontsize=6)
	axes[0, 2].set_title('[MinD:ADP] cytoplasm', fontsize=6)
	axes[0, 3].set_title('[MinE] cytoplasm', fontsize=6)
	axes[0, 4].set_title('[MinD:ATP] membrane', fontsize=6)
	axes[0, 5].set_title('[MinE:MinD:ATP] membrane', fontsize=6)
	plt.subplots_adjust(hspace=0.5)
	plt.savefig('prototypes/spatiality/minD_minE/out/min_dynamics.pdf', bbox_inches='tight')
