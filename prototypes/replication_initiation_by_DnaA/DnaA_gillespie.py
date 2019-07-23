#! /usr/bin/python
"""
Runs a stochastic Gillespie simulation of DnaA-mediated replication initiation
in E. coli.
"""

# Imports
from __future__ import division, print_function, absolute_import
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Seed for numpy
SEED = 21

# Fixed parameters
"""
N_CONDITIONS (int): Number of conditions to simulate (3 by default)
CONDITION_NAMES (list of strs): Names of the conditions to simulate
DOUBLING_TIMES (list of floats, min): Expected doubling times under each
	condition being simulated
CRITICAL_MASSES (list of floats, fg): Critical mass per origin at which 
	replication is initiated, empirically measured (source?)
INITIAL_BOX_COUNTS (list of ints): Initial counts of nonfunctional DNA boxes on
	the chromosome (from model outputs)
INITIAL_VOLUMES (list of floats, um^3): Initial volumes of simulated cells 
	under each condition being simulated
N_AVOGADRO (float): Avogadro's number
"""
N_CONDITIONS = 3
CONDITION_NAMES = ["basal", "anaerobic", "AA"]
DOUBLING_TIMES = [44., 100., 25.]  # minutes
CRITICAL_MASSES = [975., 600., 975.]  # fg
INITIAL_BOX_COUNTS = [450, 300, 750]
INITIAL_VOLUMES = [1.1, 0.4, 2.0]  # um^3
N_AVOGADRO = 6.02e23  # 1/mol


# Adjustable parameters
"""
INITIAL_DNAA_COUNTS (list of ints): Initial counts of DnaA-ATP complexes
BINDING_RATE_CONSTANT (float): Value for the rate constant of the binding 
	reaction between DnaA boxes and the DnaA-ATP complex
DISSOCIATION_RATE_CONSTANT (float): Value for the rate constant of the
	unbinding reaction of the bound DnaA box
"""
INITIAL_DNAA_COUNTS = [100, 60, 325]
BINDING_RATE_CONSTANT = 0.00022  # 1/(s*nM)
DISSOCIATION_RATE_CONSTANT = 0.0067  # 1/s


class DnaA_gillespie(object):
	"""
	Stochastic simulator of DnaA-mediated replication initiation in E. coli.
	"""
	def __init__(self):
		"""
		Loads all global variables as instance variables.
		"""
		# Load all fixed parameters
		self.n_conditions = N_CONDITIONS
		self.condition_names = CONDITION_NAMES
		self.doubling_times = np.array(DOUBLING_TIMES)*60  # minutes -> seconds
		self.critical_masses = np.array(CRITICAL_MASSES)*1e-15  # fg -> g
		self.initial_box_counts = np.array(INITIAL_BOX_COUNTS)
		self.initial_volumes = np.array(INITIAL_VOLUMES)*1e-15  # um^3 -> L
		self.n_avogadro = N_AVOGADRO

		# Load all adjustable parameters
		self.initial_DnaA_counts = np.array(INITIAL_DNAA_COUNTS)
		self.k1 = BINDING_RATE_CONSTANT/(self.n_avogadro*1e-9)  # 1/(s*nM) -> L/s
		self.k2 = DISSOCIATION_RATE_CONSTANT

		# Compute all derived parameters
		self.Kd = self.k2/self.k1

		# Initialize output values
		self.time = {}
		self.free_box_counts = {}
		self.bound_box_counts = {}
		self.free_DnaA_counts = {}

		# Seed random number generator
		np.random.seed(SEED)


	def run(self):
		"""
		Run Gillespie simulation, and store simulation results.
		"""
		# Loop through each condition
		for cond_idx in range(self.n_conditions):
			condition_name = self.condition_names[cond_idx]
			initial_volume = self.initial_volumes[cond_idx]
			doubling_time = self.doubling_times[cond_idx]
			initial_box_count = self.initial_box_counts[cond_idx]
			initial_DnaA_count = self.initial_DnaA_counts[cond_idx]

			# Get initial counts of each molecule
			init_free_box, init_bound_box, init_free_DnaA = self._get_initial_counts(
				cond_idx)

			# Initialize output arrays
			time = [0.0]
			free_box_counts = [init_free_box]
			bound_box_counts = [init_bound_box]
			free_DnaA_counts = [init_free_DnaA]

			while time[-1] < self.doubling_times[cond_idx]:
				time_old = time[-1]
				free_box_old = free_box_counts[-1]
				bound_box_old =	bound_box_counts[-1]
				free_DnaA_old = free_DnaA_counts[-1]

				# Calculate values of c's for both reactions
				c1 = self.k1/self._exponential_curve(
					initial_volume, time_old, doubling_time)
				c2 = self.k2

				# Calculate values of a's for both reactions
				a1 = c1*free_DnaA_old*free_box_old
				a2 = c2*bound_box_old
				a_total = a1 + a2
				p1 = a1/a_total

				# Calculate time until next reaction
				r1 = np.random.rand()
				tau = (1./a_total) * np.log(1./r1)
				time_new = time_old + tau

				# Compute which reaction will happen
				free_box_new = free_box_old
				bound_box_new = bound_box_old
				free_DnaA_new = free_DnaA_old

				r2 = np.random.rand()
				if r2 < p1:
					free_box_new -= 1
					bound_box_new += 1
					free_DnaA_new -= 1
				else:
					free_box_new += 1
					bound_box_new -= 1
					free_DnaA_new += 1

				# Add counts of newly synthesized boxes and DnaA complexes
				free_box_added = int(self._exponential_curve(
					initial_box_count, time_new, doubling_time)) - int(
					self._exponential_curve(
						initial_box_count, time_old, doubling_time))
				free_DnaA_added = int(self._exponential_curve(
					initial_DnaA_count, time_new, doubling_time)) - int(
					self._exponential_curve(
						initial_DnaA_count, time_old, doubling_time))

				free_box_new += free_box_added
				free_DnaA_new += free_DnaA_added

				# Append new values to list
				time.append(time_new)
				free_box_counts.append(free_box_new)
				bound_box_counts.append(bound_box_new)
				free_DnaA_counts.append(free_DnaA_new)

			# Finalize output values
			self.time[condition_name] = np.array(time)
			self.free_box_counts[condition_name] = np.array(free_box_counts)
			self.bound_box_counts[condition_name] = np.array(bound_box_counts)
			self.free_DnaA_counts[condition_name] = np.array(free_DnaA_counts)


	def plot_results(self):
		plotOutDir = os.path.join(
			os.path.dirname(os.path.abspath(__file__)), "output")

		# Loop through each condition
		for cond_idx in range(self.n_conditions):
			condition_name = self.condition_names[cond_idx]
			time = self.time[condition_name]

			fig = plt.figure()
			fig.set_size_inches(8, 11.5)

			gs = gridspec.GridSpec(4, 1)

			# Plot number of free boxes
			ax = plt.subplot(gs[0, 0])
			ax.plot(time, self.free_box_counts[condition_name])
			ax.set_ylabel("Number of free DnaA boxes")

			# Plot number of bound boxes
			ax = plt.subplot(gs[1, 0])
			ax.plot(time, self.bound_box_counts[condition_name])
			ax.set_ylabel("Number of bound DnaA boxes")

			# Plot probability that a box will be bound
			total_box_counts = self.bound_box_counts[condition_name] + self.free_box_counts[condition_name]
			ax = plt.subplot(gs[2, 0])
			ax.plot(time,
				self.bound_box_counts[condition_name].astype(np.float64)/total_box_counts)
			ax.set_ylabel("Probability of DnaA box being bound")

			# Plot number of free DnaA
			ax = plt.subplot(gs[3, 0])
			ax.plot(time, self.free_DnaA_counts[condition_name])
			ax.set_ylabel("Number of free DnaA-ATP complexes")
			ax.set_xlabel("Time [s]")

			fig.tight_layout()
			plt.savefig(os.path.join(plotOutDir, condition_name + '.png'))
			plt.close('all')


	def _get_initial_counts(self, cond_idx):
		"""
		Computes the initial counts of all molecules involved in the
		simulation, assuming that the reaction is at equilibrium.
		"""
		total_DnaA = self.initial_DnaA_counts[cond_idx]
		total_boxes = self.initial_box_counts[cond_idx]
		initial_volume = self.initial_volumes[cond_idx]

		# Computes coefficients of the quadratic equation resulting from
		# equilibrium
		coeffs = [
			1.0,
			-(total_DnaA + total_boxes + self.Kd*initial_volume),
			total_DnaA*total_boxes
			]

		# Solve quadratic equation
		roots = np.roots(coeffs)

		# Get mask for physically sensible roots, assert there's only one
		root_mask = np.logical_and(roots < min(total_DnaA, total_boxes), roots > 0)
		assert root_mask.sum() == 1

		# Get physically sensible root as counts of bound boxes
		bound_box = int(np.round(roots[np.where(root_mask)[0][0]]))
		free_box = total_boxes - bound_box
		free_DnaA = total_DnaA - bound_box

		return free_box, bound_box, free_DnaA


	def _exponential_curve(self, initial_value, t, doubling_time):
		"""
		Returns the value of a certain parameter at time t, assuming
		exponential growth from a given initial value and doubling time.
		"""
		return initial_value*np.exp(np.log(2)/doubling_time * t)


if __name__ == '__main__':
	sim = DnaA_gillespie()
	sim.run()
	sim.plot_results()
