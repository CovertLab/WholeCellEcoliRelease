"""
Functions that may be useful for future investigation in spatial model.
Some of the functions contain default values specific to E coli. Please read
before you use. References for the numbers are included in each function.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/09/2019
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import scipy.constants
from unum import Unum

from wholecell.utils import units


K_B = scipy.constants.Boltzmann*units.J/units.K  # Boltzmann constant, unit: J/K

def compute_hydrodynamic_radius(mw, mtype = None):
    '''
    This function compute the hydrodynamic diameter of a macromolecules from
    its molecular weight. It is important to note that the hydrodynamic
    diameter is mainly used for computation of diffusion constant, and can
    be different from the observed diameter under microscopes or the radius
    of gyration, especially for loose polymers such as RNAs. This function
    is not E coli specific.

    References: Bioinformatics (2012). doi:10.1093/bioinformatics/bts537

    Args:
        mw: molecular weight of the macromolecules, units: Daltons.
        mtype: There are 5 possible mtype options: protein, RNA, linear_DNA,
            circular_DNA, and supercoiled_DNA.

    Returns: the hydrodynamic radius (in unit of nm) of the macromolecules
        using the following formula
        - rp = 0.0515*MW^(0.392) nm (Hong & Lei 2008) (protein)
        - rp = 0.0566*MW^(0.38) nm (Werner 2011) (RNA)
        - rp = 0.024*MW^(0.57) nm (Robertson et al 2006) (linear DNA)
        - rp = 0.0125*MW^(0.59) nm (Robertson et al 2006) (circular DNA)
        - rp = 0.0145*MW^(0.57) nm (Robertson et al 2006) (supercoiled DNA)
    '''
    if mtype is None:
        mtype = 'protein'

    dic_rp = {'protein': (0.0515, 0.392),
                'RNA': (0.0566, 0.38),
                'linear_DNA': (0.024, 0.57),
                'circular_DNA': (0.0125, 0.59),
                'supercoiled_DNA': (0.0145, 0.57),
                }

    if isinstance(mw, Unum):
        mw_unitless = mw.asNumber(units.g / units.mol)
    else:
        mw_unitless = mw

    if mtype in dic_rp:
        rp_0, rp_power = dic_rp[mtype]
    else:
        raise KeyError("The input 'mtype' should be one of the following 5 "
                        "options: 'protein', 'RNA', 'linear_DNA', "
                        "'circular_DNA', 'supercoiled_DNA'.")

    rp = rp_0*mw_unitless**rp_power
    rp = units.nm*rp
    return rp

def compute_rg_rna(n_nt):
    '''
    This function computes the radius of gyration (in units of nm) of an RNA
    with length n_nt(in units of nt). This function is not E coli specific.
    It is important to note that radius of gyration can be very different
    from the hydrodynamic radius or the radius you can observe under
    microscopes, especially for RNAs.

    References: Nucleic Acids Research, 2017, vol 45, 5: 2919-2934.
        doi: 10.1093/nar/gkx023

    Formula: Rg = a*N**nu. a = 0.366 +/- 0.146 nm, nu = 0.50 +/- 0.05

    Args:
        n_nt: the number of nucleotide in an RNA.

    Returns:
        The radius of gyration of an RNA in the unit of nm.
    '''
    if isinstance(n_nt, Unum):
        n_nt = n_nt.asNumber(units.nt)

    a = units.nm*0.366  # unit: nm
    nu = 0.50
    return a*n_nt**nu

def compute_diffusion_constant_from_mw(mw, mtype = None,
                                       loc = None,
                                       temp = None,
                                       parameters = None):
    '''
    Warning: The default values of the 'parameters' are E coli specific.

    This function computes the hypothesized diffusion constant of
    macromolecules within the nucleoid and the cytoplasm region.
    In literature, there is no known differentiation between the diffusion
    constant of a molecule in the nucleoid and in the cytoplasm up to the
    best of our knowledge in 2019. However, there is a good reason why we
    can assume that previously reported diffusion constant are in fact the
    diffusion constant of a protein in the nucleoid region:
    (1) The image traces of a protein within a bacteria usually cross the
        nucleoid regions.
    (2) The nucleoid region, compared to the cytoplasm, should be the main
    limiting factor restricting the magnitude of diffusion constant.
    (3) The same theory of diffusion constant has been implemented to
    mammalian cells, and the term 'rh', the average hydrodynamic radius of
    the biggest crowders, are different in mammalian cytoplasm, and it seems
    to reflect the hydrodynamic radius of the actin filament (note: the
    hydrodynamic radius of actin filament should be computed based on the
    average length of actin fiber, and is not equal to the radius of the
    actin filament itself.) (ref: Nano Lett. 2011, 11, 2157-2163).
    As for E coli, the 'rh' term = 40nm, which may correspond to the 80nm
    DNA fiber. On the other hand, for the diffusion constant of E coli in
    the true cytoplasm, we will expect the value of 'rh' term to be
    approximately 10 nm, which correspond to the radius of active ribosomes.

    However, all the above statements are just hypothesis. If you want to
    compute the diffusion constant of a macromolecule in the whole E coli
    cell, you should set loc = 'nucleoid': the formula for this is obtained
    from actual experimental data set. When you set loc = 'cytoplasm', the
    entire results are merely hypothesis.

    Ref: Kalwarczyk, T., Tabaka, M. & Holyst, R.
    Bioinformatics (2012). doi:10.1093/bioinformatics/bts537

    D_0 = K_B*T/(6*pi*eta_0*rp)
    ln(D_0/D_cyto) = ln(eta/eta_0) = (xi^2/Rh^2 + xi^2/rp^2)^(-a/2)
    D_0 = the diffusion constant of a macromolecule in pure solvent
    eta_0 = the viscosity of pure solvent, in this case, water
    eta = the size-dependent viscosity experienced by the macromolecule.
    xi = average distance between the surface of proteins
    rh = average hydrodynamic radius of the biggest crowders
    a = some constant of the order of 1
    rp = hydrodynamic radius of probed molecule

    In this formula, since we allow the changes in temperature, we also
    consider the viscosity changes of water under different temperature:
    Ref: Dortmund Data Bank
    eta_0 = A*10^(B/(T-C))
    A = 2.414*10^(-5) Pa*sec
    B = 247.8 K
    C = 140 K

    Args:
        mw: molecular weight(unit: Da) of the macromolecule
        mtype: There are 5 possible mtype options: protein, RNA, linear_DNA,
        circular_DNA, and supercoiled_DNA.
        loc: The location of the molecule: 'nucleoid' or 'cytoplasm'.
        temp: The temperature of interest. unit: K.
        parameters: The 4 parameters required to compute the diffusion
            constant: xi, a, rh_nuc, rh_cyto.
            The default values are E coli specific.

    Returns:
        dc: the diffusion constant of the macromolecule, units: um**2/sec
    '''
    if mtype is None:
        mtype = 'protein'
    if loc is None:
        loc = 'nucleoid'
    if temp is None:
        temp = 300 * units.K
    if parameters is None:
        parameters = (0.51, 0.53, 42, 10)

    rp = compute_hydrodynamic_radius(mw, mtype = mtype)
    dc = compute_diffusion_constant_from_rp(rp,
                                       loc = loc,
                                       temp = temp,
                                       parameters = parameters)
    return dc

def compute_diffusion_constant_from_rp(rp,
                                       loc = None,
                                       temp = None,
                                       parameters = None):
    '''
    Warning: The default values of the 'parameters' are E coli specific.

    This is the same function as 'compute_diffusion_constant_from_mw'
    except that it accepts the hydrodynamic radius of the macromolecules
    as the input. All hypothesis and statements in
    'compute_diffusion_constant_from_mw' apply here as well.

    Args:
        rp: the hydrodynamic radius of the macromolecule, units: nm.
        loc: The location of the molecule: 'nucleoid' or 'cytoplasm'.
        temp: The temperature of interest. unit: K.
        parameters: The 4 parameters required to compute the diffusion
            constant: xi, a, rh_nuc, rh_cyto.
            The default values are E coli specific.

    Returns:
        dc: the diffusion constant of the macromolecule, units: um**2/sec
    '''
    if loc is None:
        loc = 'nucleoid'
    if temp is None:
        temp = 300 * units.K
    if parameters is None:
        parameters = (0.51, 0.53, 42, 10)

    # unpack constants required for the calculation
    xi, a, rh_nuc, rh_cyto = parameters  # unit: nm, 1, nm, nm

    # viscosity of water
    a_visc = 2.414*10**(-5)*units.Pa*units.s  # unit: Pa*sec
    b_visc = 247.8*units.K  # unit: K
    c_visc = 140*units.K  # unit: K
    eta_0 = a_visc*10**(b_visc/(temp - c_visc))  # unit: Pa*sec

    # determine Rh
    if loc == 'nucleoid':
        rh = rh_nuc  # unit: nm
    elif loc == 'cytoplasm':
        rh = rh_cyto  # unit: nm
    else:
        raise NameError(
            "The location can only be 'nucleoid' or 'cytoplasm'.")

    if not isinstance(rp, Unum):
        rp = units.nm * rp

    # compute DC(diffusion constant)
    dc_0 = K_B*temp/(6*np.pi*eta_0*rp)
    dc = dc_0*np.exp(
        -(xi**2/rh**2 + xi**2/rp.asNumber(units.nm)**2)**(-a/2))
    return dc

def compute_nucleoid_size(l_cell, d_cell,
                          length_scaling_parameters = None,
                          nucleoid_area_ratio = None):
    '''
    Warning 1: This function contains default values that are specific to
    E coli grew on M9 media under aerobic condition only. The shape and size
    of nucleoid of a bacteria can be very different across species.
    Warning 2: This function is not suitable to compute the nucleoid size of
    E coli when its shape turn filamentaous. This is because the scaling
    formula of the length of the nucleoid is obtained from a dnaC mutant
    E coli strain. According to our reference, when the cells are allowed
    to replicate their DNA normally, a constant N/C area ratio is maintained
    even for filamentous variants (treated with cephalexin). However, if the
    DNA is not allowed to replicate, the N/C area ratio will decrease as the
    cell elongate. It is therefore recommended to carefully examine the
    condition when the length of the cells grow beyond 3 um.
    Warning 3: the default values of length_scaling_parameters and
    nucleoid_area_ratio are set to be E coli specific, grew on M9 media
    under aerobic condition. The N/C area ratio for E coli grew on LB under
    aerobic condition is close to 0.4. For E coli grew under anaerobic
    condition it is close to 0.5.
    Warning 4: It is also important to note that a single cell can contain
    more than 1 nucleoid, and this formula may not be suitable for these
    cases.

    Reference on nucleoid length:
    Wu, F. et al. Curr. Biol. (2019). doi:10.1016/j.cub.2019.05.015
    Reference on nucleoid/cytoplasm area ratio:
    Gray, W. T. et al. Cell (2019). doi:10.1016/j.cell.2019.05.017

    Args:
        l_cell: the length of the cell, in units of um
        d_cell: the width of the cell, in units of um
        length_scaling_parameters: the parameters used for the scaling
            formula of the length of the nucleoid with respect to the length
            of the whole cell.
        nucleoid_area_ratio: the nucleoid/cytoplasm area ratio measured
            under microscope.
    Returns:
        l_nuc: the length of the nucleoid
        d_nuc: the diameter of the nucleoid
    '''
    if length_scaling_parameters is None:
        length_scaling_parameters = (6.6, 8.3)
    if nucleoid_area_ratio is None:
        nucleoid_area_ratio = 0.60

    l_sat, l_c = length_scaling_parameters # unit: um
    if isinstance(l_cell, Unum):
        l_cell = l_cell.asNumber(units.um)
    if isinstance(d_cell, Unum):
        d_cell = d_cell.asNumber(units.um)
    if isinstance(l_sat, Unum):
        l_sat = l_sat.asNumber(units.um)
    if isinstance(l_c, Unum):
        l_c = l_c.asNumber(units.um)

    l_nuc = l_sat * (1 - np.exp(-l_cell / l_c))
    d_nuc = nucleoid_area_ratio * l_cell * d_cell / l_nuc
    return units.um*l_nuc, units.um*d_nuc

def compute_free_volume_ratio(q, eta):
    '''
    This is a function that computes the free volume ratio of existing
    macromolecules with respect to a new molecule based on the scaled
    particle theory. This is not E coli specific. It is important to note
    that the results may not be accurate for large q or eta > 0.50.

    References:
    Malloggi, F. Soft Matter at Aqueous Interfaces.
    Lecture Notes in Physics (2016). doi:10.1007/978-3-319-24502-7.
    Chapter 3-4-2.

    Args:
        q: the size ratio between the new particle and existing particles.
            For example, if the radius of the new particle = delta, and the
            radius of the existing particle = R, then q = delta/R.
        eta: the compaction ratio, or the volume occupancy of the existing
            particle. eta = n*v_particle/v_total

    Returns:
        alpha: the free volume ratio of the space = v_free/v_total.
    '''
    a = 3*q + 3*q**2 + q**3
    b = 9/2*q**2 + 3*q**3
    c = 3*q**3
    y = eta/(1 - eta)
    q_capital = a*y + b*y**2 + c*y**3
    alpha = (1 - eta) * np.exp(-q_capital)
    return alpha

def compute_n_blob(choice_model, choice_unit, bp_dna = None):
    '''
    This function computes the number of DNA blobs within a bacteria. This
    functions contain multiple assumptions that demand careful examination.

    Previous research (Ref: Skoko D, Wong B, Johnson R, Marko J (2004)
    Biochemistry 43: 13867-13874.) directly observed the structure of
    chromosomal DNA of E coli under AFM. They observed that the DNA are
    composed of 40nm & 80nm fibers. The 80nm fibers are the main type within
    the cells, but it is hypothesized to be folded from 2 threads of 40nm
    fibers. Therefore, we hypothesized that chromosomal DNA are composed of
    DNA blobs with diameter of 40nm.

    For the model choice, there are 2 possible options. In the Flory chain
    option, the DNA within a blob is assumed to be self-avoidant. In the
    ideal chain option, the DNA within a blob is assumed to be in a melted
    state.

    For the unit choice, there are 3 possible options. You can choose single
    base pair of DNA as the smallest unit in the formation of DNA blob.
    In conventional physics, 2*persistence length is regarded as the
    standard choice of the smallest unit in a polymer. It is important to
    note that, the persistence length of DNA can be drastically different
    under in vitro & in vivo conditions. The persistence length of DNA in
    vitro is ~ 50nm, while the persistence length of DNA in vivo can be
    down to 20nm. This is because the binding of DNA binding proteins &
    DNA supercoiling can decrease the energy of DNA and makes it softer.
    Another possible choice of the smallest unit of DNA is 10bp, which
    corresponds to the average binding distance between HU protein, a
    structural DNA binding protein.

    Choosing different model setting and different smallest unit of DNA can
    result in drastically different results.

    The default value of genomic size is set to be E coli specific.

    Args:
        choice_model: the physical model of DNA within a blob. The 2
            possible options are 'ideal' & 'flory'.
        choice_unit: the choice of the smallest unit within a blob. The 3
            possible options are 'single', 'lp', or an integer number
            indicating the number of base pairs within a single unit.
        bp_dna: the number of base pairs in the chromosome. The default
            value is set to be the chromosome size of E coli.

    Returns:
        n_blob: the number of DNA blobs in a bacteria.
    '''
    if bp_dna is None:
        bp_dna = 4639221

    lp = 20  # persistence length of the DNA in vivo, unit: nm
    r_blob = 20  # radius of the DNA blob deduced from AFM, unit: nm
    length_bp_conv = 0.34  # unit: nm/bp
    bp_per_lp = lp / length_bp_conv  # unit: bp
    d_dna = 2  # unit: nm

    if isinstance(choice_unit, Unum):
        choice_unit = choice_unit.asNumber(units.nt)

    if isinstance(bp_dna, Unum):
        bp_dna = bp_dna.asNumber(units.nt)

    if choice_model == 'ideal':
        if choice_unit == 'single':
            n_bp_per_blob = 6 * (r_blob/length_bp_conv)**2
            n_blob = bp_dna / n_bp_per_blob
        elif choice_unit == 'lp':
            n_lp_per_blob = 6 * (r_blob/(2*lp))**2
            n_bp_per_blob = n_lp_per_blob * bp_per_lp * 2
            n_blob = bp_dna / n_bp_per_blob
        elif isinstance(choice_unit, int):
            length_per_unit = choice_unit * length_bp_conv
            n_unit_per_blob = 6 * (r_blob/length_per_unit)**2
            n_bp_per_blob = choice_unit * n_unit_per_blob
            n_blob = bp_dna / n_bp_per_blob
        else:
            raise NameError("Your choice_unit should either be 'single', "
                            "'lp', or an integer indicating"
                            "the number of base pair within a unit.")
    elif choice_model == 'flory':
        if choice_unit == 'single':
            n_bp_per_blob = (r_blob/
                             (3/8*length_bp_conv**4 * d_dna)**(1/5))**(5/3)
            n_blob = bp_dna / n_bp_per_blob
        elif choice_unit == 'lp':
            n_lp_per_blob = (r_blob/(3/8*lp**4*d_dna)**(1/5))**(5/3)
            n_bp_per_blob = n_lp_per_blob * bp_per_lp
            n_blob = bp_dna / n_bp_per_blob
        elif isinstance(choice_unit, int):
            length_per_unit = choice_unit * length_bp_conv
            n_unit_per_blob = (r_blob/
                               (3/8*length_per_unit**4*d_dna)**(1/5))**(5/3)
            n_bp_per_blob = choice_unit * n_unit_per_blob
            n_blob = bp_dna / n_bp_per_blob
        else:
            raise NameError("Your choice_unit should either be 'single', "
                            "'lp', or an integer indicating"
                            "the number of base pair within a unit.")
    else:
        raise NameError(
            "Your choice_model should either be 'ideal' or 'flory'.")
    return n_blob
