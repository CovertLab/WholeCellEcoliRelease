"""
PDE simulation of MinE/MinD spatial dynamics.

based on:
	Huang, K.C., Meir, Y., & Wingreen, N.S. (2003).
	Dynamic structures in Escherichia coli: spontaneous formation of MinE rings and MinD polar zones.
	Proceedings of the National Academy of Sciences, 100(22), 12724-12728.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np
from scipy.ndimage import convolve

import matplotlib

ANIMATE = False
SAVE_PLOT = True

if ANIMATE:
	matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

if ANIMATE:
	plt.ion()
	fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1)

# simulation parameters
length = 4.0  # micrometers
radius = 0.5  # micrometers
diameter = 2 * radius

bin_size = 0.05  # micrometers
bins_x = int(length / bin_size)  # number of bins along x
bins_y = int(diameter / bin_size)  # number of bins along y
bins_mem = bins_x + bins_y - 2  # membrane goes along x and wraps around to each cap's midpoint. subtract 2 for the corners
dx = length / bins_x

zero_cyto = np.zeros((bins_y, bins_x))
zero_mem = np.zeros((bins_mem,))
edges_template = np.pad(np.zeros((bins_y-2, bins_x-2)), (1,1), 'constant', constant_values=(1,1))

# indices for the membrane specifying sites along the cylinder, and the caps.
edge_length = bins_y/2 - 1
cyl_indices = range(edge_length, edge_length + bins_x)
cap1_indices = range(0, edge_length)
cap2_indices = range(edge_length + bins_x, bins_mem)

T = 40.0  # total time
dt = .0001  # time step
n = int(T / dt)  # number of iterations

# animation parameters
n_animate = 50
n_plot = 20
animate_step = n // n_animate
plot_step = n // n_plot

# laplacian kernels for diffusion
laplace_kernel_2D = np.array([[0.5, 1.0, 0.5], [1.0, -6., 1.0], [0.5, 1.0, 0.5]])
laplace_kernel_1D = np.array([1.0, -2.0, 1.0])

# constants
PI = np.pi

# chemical parameters
diffusion = 2.5  # micrometer**2/sec

k_ADP_ATP = 1  # sec**-1
k_D = 0.025  # micrometer/sec
k_dD = 0.0015  # micrometer**3/sec
k_de = 0.7  # sec**-1
k_E = 0.093  # micrometer**3/sec

## Initialize molecular fields
MinD_conc = 110 / (PI * radius ** 2)  # 1000/micrometer
MinE_conc = 100 / (PI * radius ** 2)  # 350/micrometer

# cytoplasm
DD_c = abs(np.random.normal(MinD_conc, 50, (bins_y, bins_x)))
DT_c = abs(np.random.normal(MinD_conc, 50, (bins_y, bins_x)))
E_c = abs(np.random.normal(MinE_conc, 50, (bins_y, bins_x)))

# membrane
DT_m = abs(np.random.normal(0, 50, (bins_mem,)))
EDT_m = abs(np.random.normal(0, 50, (bins_mem,)))

# Initialize arrays for plot
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
	ax1.plot(DT_m.T)
	ax2.plot(EDT_m.T)
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

def c_to_m(Z):
	# membrane contact
	contact = np.copy(zero_mem)
	cap1 = Z[1:-1,0]
	cap2 = Z[1:-1,-1]

	# membrane contact along cylinder body and both caps
	contact[cyl_indices] += (Z[0,:] + Z[-1,:])
	contact[cap1_indices] += np.flipud(cap1[:edge_length]) + cap1[edge_length:]
	contact[cap2_indices] += np.flipud(cap2[:edge_length]) + cap2[edge_length:]

	return contact / 2

def m_to_c(Z):
	# cytoplasm contact along cylinder body and both caps
	contact = np.copy(zero_cyto)
	body = Z[cyl_indices]
	cap1 = np.concatenate((np.flipud(Z[cap1_indices]), Z[cap1_indices]))
	cap2 = np.concatenate((np.flipud(Z[cap2_indices]), Z[cap2_indices]))

	contact[0,:] += body
	contact[-1,:] += body
	contact[1:-1,0] += cap1
	contact[1:-1,-1] += cap2

	return contact


# Simulate the PDE with the finite difference method.
for i in range(n):

	## Diffusion
	# cytoplasm
	d_DD_c = convolve(DD_c, laplace_kernel_2D, mode='reflect') // dx**2 * diffusion
	d_DT_c = convolve(DT_c, laplace_kernel_2D, mode='reflect') // dx**2 * diffusion
	d_E_c = convolve(E_c, laplace_kernel_2D, mode='reflect') // dx**2 * diffusion

	# membrane
	d_ET_m = convolve(DT_m, laplace_kernel_1D, mode='reflect') // dx * diffusion #**0.5
	d_EDT_m = convolve(EDT_m, laplace_kernel_1D, mode='reflect') // dx * diffusion #**0.5

	## Save current values
	# cytoplasm
	DD_c_0 = np.copy(DD_c)
	DT_c_0 = np.copy(DT_c)
	E_c_0 = np.copy(E_c)

	# membrane
	DT_m_0 = np.copy(DT_m)
	EDT_m_0 = np.copy(EDT_m)

	## Calculate reaction rates
	# get rates for cytoplasm
	m_to_c_DT_m_0 = m_to_c(DT_m_0)
	m_to_c_EDT_m_0 = m_to_c(EDT_m_0)

	rxn_1_c = (k_D + k_dD * (m_to_c_DT_m_0 + m_to_c_EDT_m_0)) * edges_template * DT_c_0
	rxn_2_c = k_E * m_to_c_DT_m_0 * E_c_0
	rxn_3_c = k_de * m_to_c_EDT_m_0
	rxn_4_c = k_ADP_ATP * DD_c_0

	# get rates for membrane
	rxn_1_m = (k_D + k_dD * (DT_m_0 + EDT_m_0)) * c_to_m(DT_c_0)
	rxn_2_m = k_E * DT_m_0 * c_to_m(E_c_0)
	rxn_3_m = k_de * EDT_m_0

	## Update the variables
	# update cytoplasm
	DD_c = DD_c_0 + (d_DD_c - rxn_4_c + rxn_3_c) * dt
	DT_c = DT_c_0 + (d_DT_c + rxn_4_c - rxn_1_c) * dt
	E_c = E_c_0 + (d_E_c + rxn_3_c - rxn_2_c) * dt

	# update membrane
	DT_m = DT_m_0 + (d_ET_m - rxn_2_m + rxn_1_m) * dt
	EDT_m = EDT_m_0 + (d_EDT_m - rxn_3_m + rxn_2_m) * dt

	# disallow negatives
	DD_c[DD_c < 0] = 0.0
	DT_c[DT_c < 0] = 0.0
	E_c[E_c < 0] = 0.0
	DT_m[DT_m < 0] = 0.0
	EDT_m[EDT_m < 0] = 0.0

	# plot the state of the system.
	if ANIMATE and i % animate_step == 0 and i < n_animate * animate_step:
		# reaction balance
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
	plt.savefig('user/min_system.pdf', bbox_inches='tight')
