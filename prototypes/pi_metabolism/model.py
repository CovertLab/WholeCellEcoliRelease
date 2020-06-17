
from __future__ import absolute_import, division, print_function

import numpy as np

import matplotlib.pyplot as plt
from six.moves import range

# Biological constants (basic)

# These are used to establish a model with numbers reasonably close to what is
# actually modeled.

TRANSLATION_RATE = 16 # AA/ribosome-sec
AA_RELATIVE_ABUNDANCE = 0.10 # fraction alanine per protein
ATP_PER_TRANSLATION = 4 # ATP per amino acid translated

CONC_ATP_BASAL = 9.6e-3 # M ATP
CONC_AA_BASAL = 2.6e-3 # M alanine

SAT_ATP = 1e-5 # M
SAT_AA = 1e-5 # M

VOLUME = 1e-15 # L - note that I don't simulate growth/dilution

N_RIBOSOME = 20000 # molecules

N_AVOGADRO = 6.022e23 # molecules per mol

# Biological constants (derived)

N_ATP_BASAL = N_AVOGADRO * CONC_ATP_BASAL * VOLUME
N_AA_BASAL = N_AVOGADRO * CONC_AA_BASAL * VOLUME

MAX_POLYMERIZATION_RATE = TRANSLATION_RATE * N_RIBOSOME

# Simulation constants

DURATION = 3600.0 # seconds
SHORT_STEPSIZE = 1.0 # seconds
LONG_STEPSIZE = 10.0 # seconds

# The 'cut' is a discrete shift partway through the simulation that constrains
# the maximum rate of alanine production to some fraction of its anticipated
# consumption rate, given the number of ribosomes and their translation rate.

CUT_TO = 0.9
CUT_AT = 1800.0 # seconds

# Controller constants

K_PROPORTIONAL = 1.0 # 1.0 is most equivalent to current model
K_INTEGRAL = 0.1

# The bootstrap value is used to initialize the accumulated error to something
# other than zero

# This was found *empirically* and will change with other constants/parameters,
# so be careful!
BOOTSTRAP_VALUE = -0.0053

def simulate(
		save_as = None,
		do_longstep = False,
		do_integral = False,
		do_bootstrap = False,
		do_cut = False,
		do_saturate = False
		):

	# Parse arguments
	stepsize = LONG_STEPSIZE if do_longstep else SHORT_STEPSIZE
	k_integral = K_INTEGRAL if do_integral else 0.0
	bootstrap_value = BOOTSTRAP_VALUE if do_bootstrap else 0.0
	sat_atp = SAT_ATP if do_saturate else 0.0
	sat_aa = SAT_AA if do_saturate else 0.0

	# Set up simulation

	# The effective translation rate adjusts for the fact that we need to boost
	# the max translation rate for sufficient production if we include the
	# effects of saturation (which necessarily lower the rate of reaction).
	effective_translation_rate = (
		TRANSLATION_RATE
		* (1 + sat_atp/CONC_ATP_BASAL)
		* (1 + sat_aa/CONC_AA_BASAL)
		)

	t = np.arange(0, DURATION, stepsize)

	timepoints = t.size

	n_atp = np.empty(timepoints)
	n_aa = np.empty(timepoints)
	n_polymerized = np.empty(timepoints)

	conc_atp_accumulated_error = bootstrap_value*stepsize * ATP_PER_TRANSLATION
	conc_aa_accumulated_error = bootstrap_value*stepsize * AA_RELATIVE_ABUNDANCE

	for i in range(timepoints):
		# Compute concentrations

		if i == 0:
			n_atp_init = N_ATP_BASAL
			n_aa_init = N_AA_BASAL

		else:
			n_atp_init = n_atp[i-1]
			n_aa_init = n_aa[i-1]

		conc_atp = n_atp_init / N_AVOGADRO / VOLUME
		conc_aa = n_aa_init / N_AVOGADRO / VOLUME

		# Simulate translation

		# This doesn't totally emulate partitioning, so an improper model can
		# allow for abundances to go negative.  This isn't too much of an issue
		# since I don't let negative numbers propagate in dangerous ways.
		if conc_atp > 0 and conc_aa > 0:
			poly_rate = (
				effective_translation_rate
				* conc_atp/(conc_atp+sat_atp)
				* conc_aa/(conc_aa+sat_aa)
				)

		else:
			poly_rate = 0

		n_poly_events = stepsize * N_RIBOSOME * poly_rate

		change_atp_translation = -ATP_PER_TRANSLATION * n_poly_events
		change_aa_translation = -AA_RELATIVE_ABUNDANCE * n_poly_events

		n_polymerized[i] = n_poly_events

		# Simulate metabolism

		conc_atp_error = conc_atp - CONC_ATP_BASAL
		conc_aa_error = conc_aa - CONC_AA_BASAL

		flux_atp = - (
			K_PROPORTIONAL * conc_atp_error
			+ k_integral * conc_atp_accumulated_error
			) / stepsize
		flux_aa = - (
			K_PROPORTIONAL * conc_aa_error
			+ k_integral * conc_aa_accumulated_error
			) / stepsize

		conc_atp_accumulated_error += conc_atp_error
		conc_aa_accumulated_error += conc_aa_error

		if do_cut and t[i] > CUT_AT:
			flux_aa = min(
				CUT_TO * MAX_POLYMERIZATION_RATE * AA_RELATIVE_ABUNDANCE / VOLUME / N_AVOGADRO,
				flux_aa
				)

		change_atp_metabolism = flux_atp * stepsize * VOLUME * N_AVOGADRO
		change_aa_metabolism = flux_aa * stepsize * VOLUME * N_AVOGADRO

		# Finalize updates

		n_atp_final = n_atp_init + change_atp_translation + change_atp_metabolism
		n_aa_final = n_aa_init + change_aa_translation + change_aa_metabolism

		n_atp[i] = n_atp_final
		n_aa[i] = n_aa_final


	# Plot

	# noinspection PyTypeChecker
	(fig, axes) = plt.subplots(
		figsize = (10, 5),
		ncols = 2,
		sharex = True,
		constrained_layout = True
		)

	# Left panel - metabolite abundances

	axes[0].step(t/60, n_atp, color = 'r', label = 'ATP')
	axes[0].axhline(N_ATP_BASAL, color = 'r', lw = 0.5, ls = '--', label = 'ATP target')

	axes[0].step(t/60, n_aa, color = 'b', label = 'AA')
	axes[0].axhline(N_AA_BASAL, color = 'b', lw = 0.5, ls = '--', label = 'AA target')

	if do_cut:
		axes[0].axvline(CUT_AT/60, color = 'k', lw = 1, alpha = 0.5)

	axes[0].set_ylim(0, max(N_ATP_BASAL, N_AA_BASAL)*1.2)
	axes[0].legend(loc = 'best')

	axes[0].set_ylabel('Molecule abundances')
	axes[0].set_xlabel('Time (min)')

	# Right panel - protein translation rate

	axes[1].step(t/60, n_polymerized / stepsize, color = 'g', label = 'AA/s')
	axes[1].axhline(MAX_POLYMERIZATION_RATE, color = 'g', ls = '--', lw = 0.5)
	if do_cut:
		axes[1].axvline(CUT_AT/60, color = 'k', lw = 1, alpha = 0.5)

	axes[1].legend(loc = 'best')

	axes[1].set_ylabel('Incorporation rate')

	axes[1].set_ylim(0, MAX_POLYMERIZATION_RATE * 1.2)

	# Finalize plots

	if save_as is None:
		plt.show()

	else:
		fig.savefig(save_as)
		plt.close(fig)

if __name__ == '__main__':
	simulate(
		do_longstep = True,
		do_integral = True,
		do_bootstrap = True,
		do_cut = True,
		do_saturate = True
		)
