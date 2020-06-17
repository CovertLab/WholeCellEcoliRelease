'''
2019/12/16
Min system modeled by compartmentalized ODEs/toast-slice model

@organization: Covert Lab, Department of Bioengineering, Stanford University
'''

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from wholecell.utils import units
from wholecell.utils import spatial_tool
from six.moves import range


# From colorbrewer2.org, qualitative 8-class set 1
COLORS_256 = [
    [55, 126, 184],
    [77, 175, 74],
    [255, 127, 0],
    [152, 78, 163],
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
    [228, 26, 28],
    ]
COLORS = np.asarray(
    [[colorValue/255. for colorValue in color] for color in COLORS_256])
R = 0.5 # radius of E coli, unit: um
A = np.pi*R**2 # cross-sectional area of E coli, unit: um^2
D = 2.5 # diffusion constant of MinD and MinE, unit: um^2/sec
CHANGE_SIZE = False
NUCLEOID = False
SAVE_PLOT = True

class MinODEMaster(object):
    def __init__(self, n, l_0):
        '''
        Args:
            n: the number of compartment in the system.
            l_0: the initial length of E coli.
        '''
        self.n = n
        self.l_0 = l_0
        self.min_d_conc = 1000 # unit: 1000/um
        self.min_e_conc = 350 # unit: 350/um

    def min_master(self, dl_dt, t_max):
        '''
        Initiate simulation and generate barcode graph for visualizing results.
        The time step is default to be 1 sec.

        Args:
            dl_dt: cell length growth rate, unit: um/sec
            t_max: simulation time.
        '''
        t_test = 100  # unit: sec
        t = np.linspace(0, t_test, 2*t_test)
        parms = [0]
        y0 = self._initial_count()
        sol_test = odeint(self._ode_main, y0.flatten(), t, args=(parms,))

        t = np.linspace(0, t_max, t_max)
        parms = [dl_dt]
        sol = odeint(self._ode_main, sol_test[-1, :], t, args=(parms,))
        l = self.l_0 + dl_dt * t
        self._generate_barcode_graph(t, sol, l, sampling = 5)

        if SAVE_PLOT:
            plt.savefig('prototypes/spatiality/min_system/output/min_ode.pdf')

        return sol

    def _compute_cytoplasm_volume(self, l):
        '''
        This function computes how each toast slice will be affected by the
        presence of nucleoid. the presence of nucleoid will decrease the volume
        of cytoplasm in each slice and decrease the area that is available for
        diffusion.
        Args:
            l: cell length

        Returns:
            dv_vec: the cytoplasm volume of each slice, in vector form
            a_diffusion_vec: area available for diffusion of each slice, in
                vector form.
        '''
        l_nuc, d_nuc = spatial_tool.compute_nucleoid_size(self.l_0, 2 * R)
        l_nuc = l_nuc.asNumber(units.um)
        d_nuc = d_nuc.asNumber(units.um)
        l_cap = (l - l_nuc) / 2
        l_slice = l / self.n
        a_nuc = A*(d_nuc/(2*R))**2

        n_affected_each = int(np.ceil(((l - l_nuc) / l_slice) / 2))
        dv_cap = A * l_slice
        dv_nuc = (A - a_nuc) * l_slice
        dv_transition = (n_affected_each * l_slice - l_cap) * (A - a_nuc) + (
                    l_cap - (n_affected_each - 1) * l_slice) * A
        dv_vec = np.ones(self.n) * dv_nuc
        dv_vec[(n_affected_each - 1, -n_affected_each),] = dv_transition
        dv_vec[0:(n_affected_each - 1)] = dv_cap

        a_diffusion_vec = np.ones((n, 2)) * (A - a_nuc)
        a_diffusion_vec[(n_affected_each - 1), :] = [A, (A - a_nuc)]
        a_diffusion_vec[-n_affected_each, :] = [(A - a_nuc), A]
        a_diffusion_vec[0:(n_affected_each - 1), :] = [A, A]

        if n_affected_each >= 2:
            dv_vec[(-n_affected_each + 1):] = dv_cap
            a_diffusion_vec[(-n_affected_each + 1):, :] = [A, A]

        return dv_vec, a_diffusion_vec

    def _initial_count(self):
        '''
        Find the initial condition of the toast-slice model for Min system.
        The number of compartments = n, and the initial length of the cell = l_0.

        Returns:
            y0: the initial condition of ODE, in molecular counts.
        '''
        # total number of MinD & MinE
        n_min_d_total = self.min_d_conc * self.l_0
        n_min_e_total = self.min_e_conc * self.l_0

        # divide into each compartment. no membrane components in the beginning.
        y0 = np.random.rand(self.n + 2, 5)
        if NUCLEOID:
            dv_vec, _ = self._compute_cytoplasm_volume(l_0)
            y0[1:-1, :] *= dv_vec.reshape((-1, 1))

        y0[(0, -1), :] = 0
        y0[:, (2, 4)] = 0
        y0[:, (0, 1)] *= n_min_d_total / sum(sum(y0[:, (0, 1)]))
        y0[:, 3] *= n_min_e_total / sum(y0[:, 3])

        return y0.flatten()

    def _extra_count(self, extra_min_d, extra_min_e, l = None):
        '''
        Stochastically add extra MinD and MinE into the system.
        The newly added MinD & MinE always in cytosol.

        Args:
            extra_min_d: extra number of MinD that should be distributed
            extra_min_e: extra number of MinE that should be distributed
            l: cell length

        Returns:
            y_extra: Stochastically distributed new MinD and MinE
        '''
        if l is not None and NUCLEOID:
            dv_vec, _ = self._compute_cytoplasm_volume(l)
            y_extra = np.ones((self.n + 2, 5))
            y_extra[1:-1, :] *= dv_vec.reshape((-1, 1))
            y_extra[(0, -1), :] = 0
            y_extra[:, (2, 4)] = 0
            y_extra[:, (0, 1)] *= extra_min_d / sum(sum(y_extra[:, (0, 1)]))
            y_extra[:, 3] *= extra_min_e / sum(y_extra[:, 3])
        else:
            y_extra = np.zeros((self.n + 2, 5))
            y_extra[:, (0, 1)] = extra_min_d / (2 * self.n)
            y_extra[:, 3] = extra_min_e / self.n
            y_extra[(0, -1), :] = 0

        return y_extra.flatten()

    def _ode_main(self, y, t, parms):
        '''
        The main ODE of the Min system.

        units: [D_ADP], [D_ATP], [E_c] = 1/um^3; [D_m], [DE_m] = 1/um^2;
        units: k1, k5: 1/sec; k3, k4: um^3/sec; k2: um/sec

        d[D_ADP]/dt = - k1[D_ADP] + "k5[DE_m]"
        d[D_ATP]/dt = k1[D_ADP] - "k2[D_ATP]" - "k3[D_ATP]([D_m] + [DE_m])"
        d[D_m]/dt = k2[D_ATP] + k3[D_ATP]([D_m] + [DE_m]) - k4[E_c][D_m]
        d[E_c]/dt = - "k4[E_c][D_m]" + "k5[DE_m]"
        d[DE_m]/dt = k4[E_c][D_m] - k5[DE_m]

        Ref:
        Proc. Natl. Acad. Sci. U. S. A. (2003). doi:10.1073/pnas.2135445100
        '''
        [dl_dt] = parms
        l = self.l_0 + dl_dt * t

        # 0-1. kinetic parameters:
        k1 = 1.0  # unit: 1/sec
        k2 = 0.025  # unit: um/sec
        k3 = 0.0015  # unit: um^3/sec
        k4 = 0.093  # unit: um^3/sec
        k5 = 0.7  # unit: 1/sec

        # 0-2. size & bin related parameters
        dx = l / self.n  # unit: um, the width of the toast slice
        da = 2 * np.pi * R * dx  # unit: um^2, the area of toast crust

        if NUCLEOID:
            dv_vec, a_diffusion_vec = self._compute_cytoplasm_volume(l)
        else:
            dv = A * dx  # unit: um^3, the volume of each toast slice

        # 0-3. reshape and divide the initial conditions
        y = np.reshape(y, (self.n + 2, 5))
        y_left_cap = y[0, :]
        y_middle = y[1:(n + 1), :]
        y_right_cap = y[-1, :]

        # 1-1. kinetic reactions of middle part
        if NUCLEOID:
            d_d_adp_c_dt_middle = (- k1 * y_middle[:, 0]
                                   + k5 * y_middle[:, 4]).reshape((-1, 1))
            d_d_atp_c_dt_middle = (+ k1 * y_middle[:, 0]
                                   - k2 * y_middle[:, 1] * da/dv_vec
                                   - k3 * y_middle[:, 1]
                                   * (y_middle[:, 2] + y_middle[:, 4])/dv_vec).reshape((-1, 1))
            d_d_m_dt_middle = (+ k2 * y_middle[:, 1] * da/dv_vec
                               + k3 * y_middle[:, 1] * (y_middle[:, 2] + y_middle[:, 4]) / dv_vec
                               - k4 * y_middle[:, 3] * y_middle[:, 2]/dv_vec).reshape((-1, 1))
            d_e_c_dt_middle = (- k4 * y_middle[:, 3] * y_middle[:, 2]/dv_vec
                               + k5 * y_middle[:, 4]).reshape((-1, 1))
            d_de_m_dt_middle = (+ k4 * y_middle[:, 3] * y_middle[:, 2]/dv_vec
                                - k5 * y_middle[:, 4]).reshape((-1, 1))
        else:
            d_d_adp_c_dt_middle = (- k1 * y_middle[:, 0]
                                   + k5 * y_middle[:, 4]).reshape((-1, 1))
            d_d_atp_c_dt_middle = (+ k1 * y_middle[:, 0]
                                   - k2 * y_middle[:, 1] * da/dv
                                   - k3 * y_middle[:, 1] * (
                                           y_middle[:, 2]
                                           + y_middle[:, 4])/dv).reshape((-1, 1))
            d_d_m_dt_middle = (+ k2 * y_middle[:, 1] * da/dv
                               + k3 * y_middle[:, 1] * (y_middle[:, 2] + y_middle[:, 4])/dv
                               - k4 * y_middle[:, 3] * y_middle[:, 2]/dv).reshape((-1, 1))
            d_e_c_dt_middle = (- k4 * y_middle[:, 3] * y_middle[:, 2]/dv
                               + k5 * y_middle[:, 4]).reshape((-1, 1))
            d_de_m_dt_middle = (+ k4 * y_middle[:, 3] * y_middle[:, 2]/dv
                                - k5 * y_middle[:, 4]).reshape((-1, 1))

        # 1-2. combine the kinetic reactions of middle part together
        dydt_middle_kinetic = np.hstack((
            d_d_adp_c_dt_middle,
            d_d_atp_c_dt_middle,
            d_d_m_dt_middle,
            d_e_c_dt_middle,
            d_de_m_dt_middle))

        # 2-1. consider the kinetic reactions with membrane
        def kinetic_rxn_with_membrane(dydt_middle_kinetic,
                                      left_right, y_cap, dv):
            # left or right membrane
            if left_right == 'left':
                index = 0
            else:
                index = -1

            # kinetic reactions:
            if NUCLEOID:
                d_d_m_dt_cap = (+ k2 * y_middle[index, 1]*A/dv_vec[index]
                                + k3 * y_middle[index, 1]*(y_cap[2] + y_cap[4])/dv_vec[index]
                                - k4 * y_middle[index, 3]*y_cap[2]/dv_vec[index])
                d_de_m_dt_cap = (+ k4 * y_middle[index, 3]*y_cap[2]/dv_vec[index]
                                 - k5 * y_cap[4])
                dydt_cap_kinetic = np.array([0, 0, d_d_m_dt_cap, 0, d_de_m_dt_cap])

                dydt_middle_kinetic[index, 0] += k5 * y_cap[4]
                dydt_middle_kinetic[index, 1] -= (k2*y_middle[index, 1]*A/dv_vec[index]
                                                  + k3*y_middle[index, 1]*(y_cap[2] + y_cap[4])/dv_vec[index])
                dydt_middle_kinetic[index, 3] += (k5*y_cap[4]
                                                  - k4*y_middle[index, 3]* y_cap[2]/dv_vec[index])
            else:
                d_d_m_dt_cap = (+ k2 * y_middle[index, 1] * A/dv
                                + k3 * y_middle[index, 1] * (y_cap[2] + y_cap[4])/dv
                                - k4 * y_middle[index, 3] * y_cap[2]/dv)
                d_de_m_dt_cap = (+ k4 * y_middle[index, 3] * y_cap[2]/dv
                                 - k5 * y_cap[4])
                dydt_cap_kinetic = np.array([0, 0, d_d_m_dt_cap, 0, d_de_m_dt_cap])

                dydt_middle_kinetic[index, 0] += k5 * y_cap[4]
                dydt_middle_kinetic[index, 1] -= (k2 * y_middle[index, 1] * A/dv
                                                  + k3 * y_middle[index, 1] * (
                                                          y_cap[2] + y_cap[4])/dv)
                dydt_middle_kinetic[index, 3] += (k5 * y_cap[4]
                                                  - k4 * y_middle[index, 3]*y_cap[2]/dv)

            return dydt_middle_kinetic, dydt_cap_kinetic

        if NUCLEOID:
            dydt_middle_kinetic, dydt_left_cap_kinetic = kinetic_rxn_with_membrane(
                dydt_middle_kinetic, 'left', y_left_cap, dv_vec)
            dydt_middle_kinetic, dydt_right_cap_kinetic = kinetic_rxn_with_membrane(
                dydt_middle_kinetic, 'right', y_right_cap, dv_vec)
        else:
            dydt_middle_kinetic, dydt_left_cap_kinetic = kinetic_rxn_with_membrane(
                dydt_middle_kinetic, 'left', y_left_cap, dv)
            dydt_middle_kinetic, dydt_right_cap_kinetic = kinetic_rxn_with_membrane(
                dydt_middle_kinetic, 'right', y_right_cap, dv)

        # diffusion. Proteins on membrane do not diffuse.
        dydt_middle_diffusion = np.zeros((self.n, 5))
        if NUCLEOID:
            for i in range(self.n):
                if i == 0:
                    diffusion = a_diffusion_vec[i, 1] * D * (
                            y_middle[i, :] - y_middle[i + 1, :]) / dx
                    dydt_middle_diffusion[i, :] = - diffusion
                elif i == (n - 1):
                    diffusion = a_diffusion_vec[i, 0] * D * (
                            y_middle[i - 1, :] - y_middle[i, :]) / dx
                    dydt_middle_diffusion[i, :] = + diffusion
                else:
                    diffusion_1 = a_diffusion_vec[i, 1] * D * (
                            y_middle[i, :] - y_middle[i + 1, :]) / dx
                    diffusion_2 = a_diffusion_vec[i, 0] * D * (
                            y_middle[i - 1, :] - y_middle[i, :]) / dx
                    dydt_middle_diffusion[i, :] = - diffusion_1 + diffusion_2
        else:
            for i in range(self.n):
                if i == 0:
                    diffusion = A * D * (y_middle[i, :] - y_middle[i + 1, :])/dx
                    dydt_middle_diffusion[i, :] = - diffusion
                elif i == (self.n - 1):
                    diffusion = A * D * (y_middle[i - 1, :] - y_middle[i, :])/dx
                    dydt_middle_diffusion[i, :] = + diffusion
                else:
                    diffusion_1 = A * D * (y_middle[i, :] - y_middle[i + 1, :])/dx
                    diffusion_2 = A * D * (y_middle[i - 1, :] - y_middle[i, :])/dx
                    dydt_middle_diffusion[i, :] = - diffusion_1 + diffusion_2
        dydt_middle_diffusion[:, (2, 4)] = 0

        # add extra component
        extra_min_d_dt = dl_dt * self.min_d_conc  # unit: molecules
        extra_min_e_dt = dl_dt * self.min_e_conc  # unit: molecules
        if NUCLEOID:
            dy_extra_dt = self._extra_count(extra_min_d_dt, extra_min_e_dt, l)
        else:
            dy_extra_dt = self._extra_count(extra_min_d_dt, extra_min_e_dt)

        # combine together
        dydt = np.vstack((dydt_left_cap_kinetic,
                          dydt_middle_kinetic + dydt_middle_diffusion,
                          dydt_right_cap_kinetic))
        dydt = dydt.flatten()
        dydt += dy_extra_dt
        return dydt

    def _generate_barcode_graph(self, t, sol, l, sampling = None):
        '''
        Generate the barcode graph of the solution.

        Args:
            t: the time points
            sol: the solution of _ode_main
            l: the time series of cell length.
            sampling: sampling interval of generating graph.
        '''
        if sampling is None:
            sampling = 200 # sampling rate for generating plot.

        # separate solutions into 5 species
        min_d_adp = sol[:, 0::5]
        min_d_atp = sol[:, 1::5]
        min_d_m = sol[:, 2::5]
        min_e_c = sol[:, 3::5]
        min_de_m = sol[:, 4::5]

        # plot
        def plot_component_barcode(
                min_element, title_label, n, l, t, c_code, ax):
            min_element_normalized = np.round(min_element/np.max(min_element), 7)
            m_thickness = 0.3  # membrane thickness on plot
            for i in range(0, len(t), sampling):
                tp = t[i]
                l_current = l[i]
                for j in range(n + 2):
                    if j == 0:
                        ax.plot([0, m_thickness],
                                [tp, tp],
                                lw=2,
                                Color=COLORS[c_code]*min_element_normalized[i, j])
                    elif j == n + 1:
                        ax.plot([l_current + m_thickness, l_current + 2*m_thickness],
                                [tp, tp],
                                lw=2,
                                Color=COLORS[c_code]*min_element_normalized[i, j])
                    else:
                        ax.plot([m_thickness + l_current/n*(j - 1), m_thickness + l_current/n*j],
                                [tp, tp],
                                lw=2,
                                Color=COLORS[c_code]*min_element_normalized[i, j])
            ax.set_xlabel('length($\mu$m)')
            ax.set_ylabel('time(s)')
            ax.invert_yaxis()
            ax.set_title(title_label)
            return

        fig, axes = plt.subplots(figsize = (8, 20), nrows = 5, ncols = 1)
        plot_component_barcode(
            min_d_adp, r'$MinD_{ADP}$', self.n, l, t, 0, axes[0])
        plot_component_barcode(
            min_d_atp, r'$MinD_{ATP}$', self.n, l, t, 1, axes[1])
        plot_component_barcode(
            min_d_m, r'$MinD_{m}$', self.n, l, t, 3, axes[2])
        plot_component_barcode(
            min_e_c, r'$MinE_{c}$', self.n, l, t, 5, axes[3])
        plot_component_barcode(
            min_de_m, r'$MinDE_{m}$', self.n, l, t, 6, axes[4])
        plt.tight_layout()

        return

if __name__ == "__main__":
    l_0 = 4 # unit: um
    n = 6 # number of compartments in the beginning
    min_ode_master = MinODEMaster(n, l_0)
    if CHANGE_SIZE:
        dl_dt = 6 / 2400
    else:
        dl_dt = 0
    t_max = 500
    sol = min_ode_master.min_master(dl_dt, t_max)
