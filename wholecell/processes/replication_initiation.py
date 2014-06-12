#!/usr/bin/env python

"""
ReplicationInitiation

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process


class ReplicationInitiation(wholecell.processes.process.Process):
	""" ReplicationInitiation """

	_name = "ReplicationInitiation"

	# Constructor
	def __init__(self):


		super(ReplicationInitiation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ReplicationInitiation, self).initialize(sim, kb)

		# Load parameters
		


		# Create views
		self.pi = self.bulkMoleculeView(['PI[c]'])
		self.atp = self.bulkMoleculeView(['ATP[c]'])
		self.adp = self.bulkMoleculeView(['ADP[c]'])

		self.dnaA_box = self.bulkChromosomesView(['R1_dnaA',
												'R2_dnaA',
												'R3_dnaA',
												'R4_dnaA',
												'R5_dnaA'])

		self.dnaA_box_atp_polymer = self.bulkChromosomesView(['R1_dnaA_atp_polymer',
																'R2_dnaA_atp_polymer',
																'R3_dnaA_atp_polymer',
																'R4_dnaA_atp_polymer',
																'R5_dnaA_atp_polymer'])

		self.dnaA_box_adp_polymer = self.bulkChromosomesView(['R1_dnaA_adp_polymer',
																'R2_dnaA_adp_polymer',
																'R3_dnaA_adp_polymer',
																'R4_dnaA_adp_polymer',
																'R5_dnaA_adp_polymer'])

	def calculateRequest(self):
		pass

	# Calculate temporal evolution
	def evolveState(self):
		pass