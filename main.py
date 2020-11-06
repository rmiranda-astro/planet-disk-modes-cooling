import numpy as np
from math import pi
import scipy.integrate as integ
import matplotlib.pyplot as plt
import cmath
import time
import multiprocessing
from master_equation_cool import *
from parameters import *

# Define the logarithmic radial grid
r = np.logspace(np.log10(r_in), np.log10(r_out), int(r_steps))

# Get the grid index of the planet's orbital radius r = 1
i1 = np.argmin(np.abs(r - 1.0))

# Surface density \Sigma at radius s, assuming a s^-p power law
# (p is specified in the parameter file)
def Sigma_func(s):
    Sigma = 1.0 / s ** p
    return Sigma


# Derivative of the surface density
def dSigmadr_func(s):
    dSigmadr = -p / s ** (p + 1.0)
    return dSigmadr


# Keplerian frequency \Omega_K at radius s
def OmegaK_func(s):
    OmegaK = 1.0 / s ** (3.0 / 2.0)
    return OmegaK


# Orbital frequency \Omega (Keplerian modified by pressure support) at radius s
def Omega_func(s):
    Omega = OmegaK_func(s)
    Omega *= np.sqrt(1.0 - h0 ** 2.0 * (p + q) * s ** (1.0 - q))
    return Omega


# Derivative of the orbital frequency
def dOmegadr_func(s):
    dOmegadr = -3.0 / 2.0 * OmegaK_func(s) / s
    dOmegadr *= (
        -(s ** (-q))
        * (h0 ** 2.0 * (2.0 + q) * (p + q) * s - 3.0 * s ** q)
        / 3.0
        / np.sqrt(1.0 - h0 ** 2.0 * (p + q) * s ** (1.0 - q))
    )
    return dOmegadr


# Radial epicyclic frequency at radius s
def kappa_func(s):
    kappa = OmegaK_func(s)
    kappa *= np.sqrt(1.0 + h0 ** 2.0 * (q - 2.0) * (p + q) * s ** (1.0 - q))
    return kappa


# Derivative of the radial epicyclic frequency
def dkappadr_func(s):
    dkappadr = -3.0 / 2.0 * OmegaK_func(s) / s
    dkappadr *= (
        s ** (-q)
        * (h0 ** 2.0 * (p + q) * (q ** 2.0 - 4.0) * s + 3.0 * s ** q)
        / 3.0
        / np.sqrt(1.0 + h0 ** 2.0 * (q - 2.0) * (p + q) * s ** (1.0 - q))
    )
    return dkappadr


# Perturbation frequency for azimuthal number m
def omega_func(m):
    omega = m
    return omega


# The squared frequency D = \kappa^2 - \tilde\omega^2 at radius s for azimuthal number m
def D_func(s, m):
    omega = omega_func(m)
    D = kappa_func(s) ** 2.0 - (m * Omega_func(s) - omega) ** 2.0
    return D


# The derivative of D
def dDdr_func(s, m):
    omega = omega_func(m)
    dDdr = 2.0 * kappa_func(s) * dkappadr_func(s) - 2.0 * m * (
        m * Omega_func(s) - omega
    ) * dOmegadr_func(s)
    return dDdr


# Isothermal sound speed at radius s, assuming s^-q temperature profile
# (q is specified in the parameter file)
def csiso_func(s):
    csiso = h0 / s ** (q / 2.0)
    return csiso


# Derivative of isothermal sound speed at radius s
def dcsisodr_func(s):
    dcsisodr = -(q / 2.0) * h0 / s ** (q / 2.0 + 1.0)
    return dcsisodr


# Adiabatic sound speed at radius s
def csadi_func(s):
    csadi = np.sqrt(gamma) * csiso_func(s)
    return csadi


# Derivative of the adiabatic sound speed
def dcsadidr_func(s):
    dcsadidr = np.sqrt(gamma) * dcsisodr_func(s)
    return dcsadidr


# Squared Brunt-Vaisala N_r^2 at radius s
def Nr2_func(s):
    csadi = csadi_func(s)
    zeta = (q + (1.0 - gamma) * p) * (q + p) / gamma ** 2.0
    Nr2 = -zeta * csadi ** 2.0 / s ** 2.0
    return Nr2


# Derivative of Brunt-Vaisala frequency
def dNr2dr_func(s):
    csadi = csadi_func(s)
    zeta = (q + (1.0 - gamma) * p) * (q + p) / gamma ** 2.0
    dNr2dr = zeta * (q + 2.0) * csadi ** 2.0 / s ** 3.0
    return dNr2dr


# Inverse entropy lengthscale 1/L_S at radius s
def LSinv_func(s):
    eta = (q + (1.0 - gamma) * p) / gamma
    LSinv = -eta / s
    return LSinv


# Derivative of the inverse entropy lengthscale
def dLSinvdr_func(s):
    eta = (q + (1.0 - gamma) * p) / gamma
    dLSinvdr = eta / s ** 2.0
    return dLSinvdr


# Inverse temperature lengthscale 1/L_T at radius s
def LTinv_func(s):
    LTinv = -q / s
    return LTinv


# Derivative of the inverse temperature lengthscale
def dLTinvdr_func(s):
    dLTinvdr = q / s ** 2.0
    return dLTinvdr


# Cooling timescale t_c, assuming constant dimensionless cooling timescale \beta
# (\beta = \Omega t_c)
def tc_func(s):
    tc = beta / Omega_func(s)
    return tc


# Derivative of the cooling timescale
def dtcdr_func(s):
    dtcdr = -beta / Omega_func(s) ** 2.0 * dOmegadr_func(s)
    return dtcdr


# Real part of the WKB wavenumber k at radius s for mode with azimuthal number m
def k_func(s, m):
    omega_tilde = omega_func(m) - m * Omega_func(s)
    tc = tc_func(s)
    cseff = (
        (1.0 + (gamma * omega_tilde * tc) ** 2.0) / (1.0 + (omega_tilde * tc) ** 2.0)
    ) ** (1.0 / 4.0) * csiso_func(s)
    Deff = D_func(s, m) + (gamma * omega_tilde * tc) ** 2.0 / (
        1.0 + (gamma * omega_tilde * tc) ** 2.0
    ) * Nr2_func(s)
    k = cmath.sqrt(-Deff) / cseff
    return k


# The squared frequency D_S = D + N_r^2 (buoyancy-modified analogue of D)
def DS_func(s, m):
    D = D_func(s, m)
    Nr2 = Nr2_func(s)
    DS = D + Nr2
    return DS


# Derivative of D_S
def dDSdr_func(s, m):
    dDdr = dDdr_func(s, m)
    dNr2dr = dNr2dr_func(s)
    dDSdr = dDdr + dNr2dr
    return dDSdr


# The coefficient C in the relation dh' = C*dh that is satisfied by the enthalpy
# perturbation dh for a WKB outgoing wave (taking into account the variation of the
# background disk properties)
def wkb_func(s, m):
    omega_tilde = omega_func(m) - m * Omega_func(s)
    tc = tc_func(s)
    rek = k_func(s, m)
    imk = (
        1.0
        / 2.0
        * (gamma - 1.0)
        * omega_tilde
        * tc
        / (1.0 + gamma * (omega_tilde * tc) ** 2.0)
        * rek
    )
    k = rek + 1.0j * imk
    wkb = 1.0j * k
    wkb += (
        dDdr_func(s, m) / (4.0 * D_func(s, m))
        - dSigmadr_func(s) / (2.0 * Sigma_func(s))
        + dcsisodr_func(s) / (2.0 * csiso_func(s))
        - 1.0 / (2.0 * s)
    )
    wkb += (
        1.0
        / 2.0
        * (1.0 + 1.0j * gamma * omega_tilde * tc)
        / (1.0 + (gamma * omega_tilde * tc) ** 2.0)
        * LTinv_func(s)
    )
    return wkb


# Derivative of the phase for mode f ([dh,dh'])
def dphidr_func(f):
    dphidr = (
        np.real(f[:, 0]) * np.imag(f[:, 1]) - np.real(f[:, 1]) * np.imag(f[:, 0])
    ) / (np.real(f[:, 0]) ** 2.0 + np.imag(f[:, 0]) ** 2.0)
    return dphidr


# Smoothed Laplace coefficient b_s^m with smoothing parameter epsilon at radius \alpha
def laplace_func(s, m, alpha, epsilon):
    phi = np.linspace(0.0, 2.0 * pi, 1001, endpoint=True)
    integrand = (
        np.cos(m * phi)
        * (1.0 + alpha ** 2.0 - 2.0 * alpha * np.cos(phi) + epsilon ** 2.0) ** (-s)
        / pi
    )
    integral = integ.simps(integrand, phi)
    return integral


# Derivative of smoothed Laplace coefficient b_s^m
def laplace_deriv_func(s, m, alpha, epsilon):
    phi = np.linspace(0.0, 2.0 * pi, 1001, endpoint=True)
    integrand = (
        2.0
        * s
        * np.cos(m * phi)
        * (np.cos(phi) - alpha)
        * (1.0 + alpha ** 2.0 - 2.0 * alpha * np.cos(phi) + epsilon ** 2.0)
        ** (-(s + 1.0))
        / pi
    )
    integral = integ.simps(integrand, phi)
    return integral


# Second derivative of smoothed Laplace coefficient b_s^m
def laplace_2deriv_func(s, m, alpha, epsilon):
    phi = np.linspace(0.0, 2.0 * pi, 1001, endpoint=True)
    integrand = (
        2.0
        * s
        * np.cos(m * phi)
        * (
            s
            + (1.0 + 2.0 * s) * (alpha ** 2.0 - 2.0 * alpha * np.cos(phi))
            + (1.0 + s) * np.cos(2.0 * phi)
            - epsilon ** 2.0
        )
        * (1.0 + alpha ** 2.0 - 2.0 * alpha * np.cos(phi) + epsilon ** 2.0)
        ** (-(s + 2.0))
        / pi
    )
    integral = integ.simps(integrand, phi)
    return integral


# Planetary potential at radius s for harmonic m
def Phi_m_func(s, m):
    Phi = -laplace_func(1.0 / 2.0, m, s, soft)
    if m == 1 and indirect == 1:
        Phi += s
    return Phi


# Derivative of the planetary potential harmonic
def dPhidr_m_func(s, m):
    dPhidr = -laplace_deriv_func(1.0 / 2.0, m, s, soft)
    if m == 1 and indirect == 1:
        dPhidr += 1.0
    return dPhidr


# Second derivative of the planetary potential harmonic
def d2Phidr2_m_func(s, m):
    d2Phidr2 = -laplace_2deriv_func(1.0 / 2.0, m, s, soft)
    return d2Phidr2


def f_deriv_func(r, f, m, ih):

    """
    Get the first and second derivatives of the enthalpy perturbation dh using
    "f" = [Re(dh),Im(dh),Re(dh'),Im(dh')], at radius "r" for azimuthal number "m".
    Returns "dfdr" = [Re(dh'),Im(dh'),Re(dh''),Im(dh'')], as needed for the integrators
    in the function "integrate_ode". The first derivative dh' is simply copied from "f",
    and the second derivative dh'' is computed according the master equation using dh
    and dh' (with coefficients calculated by the function "master_equation_cool"). The
    flag "ih" controls whether or not to include the planetary potential terms.
    """

    Omega = Omega_func(r)
    dOmegadr = dOmegadr_func(r)

    Sigma = Sigma_func(r)
    dSigmadr = dSigmadr_func(r)

    cs = csadi_func(r)

    omega_tilde = omega_func(m) - m * Omega

    D = D_func(r, m)
    dDdr = dDdr_func(r, m)

    LSinv = LSinv_func(r)
    dLSinvdr = dLSinvdr_func(r)

    LTinv = LTinv_func(r)
    dLTinvdr = dLTinvdr_func(r)

    Nr2 = Nr2_func(r)
    dNr2dr = dNr2dr_func(r)

    tc = tc_func(r)
    dtcdr = dtcdr_func(r)

    if ih == 1:

        Phi = Phi_m_func(r, m)
        dPhidr = dPhidr_m_func(r, m)
        d2Phidr2 = d2Phidr2_m_func(r, m)

    else:

        Phi, dPhidr, d2Phidr2 = (0.0 for i in range(3))

    dfdr = np.empty(4)

    C1R, C1I, C0R, C0I, PsiR, PsiI = master_equation_cool(
        r,
        gamma,
        cs,
        Omega,
        dOmegadr,
        Sigma,
        dSigmadr,
        LTinv,
        dLTinvdr,
        LSinv,
        dLSinvdr,
        omega_tilde,
        tc,
        dtcdr,
        D,
        dDdr,
        Nr2,
        dNr2dr,
        pole_disp,
        Phi,
        dPhidr,
        d2Phidr2,
        m,
    )

    dfdr[0] = f[2]
    dfdr[1] = f[3]
    dfdr[2] = PsiR - C1R * f[2] + C1I * f[3] - C0R * f[0] + C0I * f[1]
    dfdr[3] = PsiI - C1R * f[3] - C1I * f[2] - C0R * f[1] - C0I * f[0]

    return dfdr


def integrate_ode(x, bc, m, ih):

    """
    Integrate the master equation on the radial grid "x", starting from the corotation
    radius and integrating both inwards and outwards. The array "bc" gives the value of
    the (complex) enthalpy perturbation dh and its derivative dh' at the corotation
    radius. The flag "ih" turns the planetary potential terms on or off. Returns
    "f" = [dh,dh'].
    """

    i_co = np.argmin(np.abs(m * Omega_func(x) - omega_func(m)))

    f = np.zeros((len(x), 2), complex)
    f[i_co] = bc

    # Initialize the integrator. Use the flexible lsoda integrator unless dealing with
    # the m = 1 mode with "m1_switch" set to 1 in the parameter file, in which case
    # switch to an explicit Runge-Kutta integrator that is much slower but more robust
    # against some spurious behavior that can arise in the inner disk for m = 1 when
    # the lsoda integrator is used.

    if m1_switch == 1 and m == 1:

        solver = integ.ode(f_deriv_func).set_integrator("dopri5")

    else:

        solver = integ.ode(f_deriv_func).set_integrator(
            "lsoda", max_order_ns=4, max_order_s=4
        )

    # Integrate outwards from corotation.

    solver.set_initial_value(
        [np.real(bc[0]), np.imag(bc[0]), np.real(bc[1]), np.imag(bc[1])], x[i_co]
    ).set_f_params(m, ih)

    for i in range(i_co + 1, len(x)):

        solver.integrate(x[i])
        f[i, 0] = solver.y[0] + 1.0j * solver.y[1]
        f[i, 1] = solver.y[2] + 1.0j * solver.y[3]

    solver.set_initial_value(
        [np.real(bc[0]), np.imag(bc[0]), np.real(bc[1]), np.imag(bc[1])], x[i_co]
    ).set_f_params(m, ih)

    # Integrate inwards from corotation.

    for i in range(i_co - 1, -1, -1):

        solver.integrate(x[i])
        f[i, 0] = solver.y[0] + 1.0j * solver.y[1]
        f[i, 1] = solver.y[2] + 1.0j * solver.y[3]

    return f


def homogeneous_coefficients(yh1, yh2, yih, Ci, Co, i_in, i_out):

    """
    Given two linearly independent homogeneous mode solutions "yh1" and "yh2", an
    inhomogeneous solution "yih" (all three are of the form [dh,dh']), the WKB
    coefficients at the inner/outer boundaries "Ci" and "Co" (see the function
    "wkb_func"), and the grid indices of the inner and outer boundaries "i_in" and
    "i_out", compute the coefficients a1 and a2 such that a1*yh1 + a2*yh2 + yih
    satisfies the WKB outgoing boundary condition at the inner and outer boundaries.
    """

    mat = np.empty((2, 2), complex)
    mat[0, 0] = yh1[i_in, 1] - Ci * yh1[i_in, 0]
    mat[0, 1] = yh2[i_in, 1] - Ci * yh2[i_in, 0]
    mat[1, 0] = yh1[i_out, 1] - Co * yh1[i_out, 0]
    mat[1, 1] = yh2[i_out, 1] - Co * yh2[i_out, 0]

    det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]

    vec = np.empty(2, complex)
    vec[0] = Ci * yih[i_in, 0] - yih[i_in, 1]
    vec[1] = Co * yih[i_out, 0] - yih[i_out, 1]

    a1 = (mat[1, 1] * vec[0] - mat[0, 1] * vec[1]) / det
    a2 = (-mat[1, 0] * vec[0] + mat[0, 0] * vec[1]) / det

    return a1, a2


def phase_error_diagnostic(f, m, i_in, i_out):

    """
    Given the mode "f", azimuthal number "m", and inner/outer boundary grid indices
    "i_in" and "i_out", return the peak-to-peak amplitude of the phase gradient error
    (the difference between the derivative of the mode phase and the WKB wavenumber k)
    within one scale height of each boundary. The mode solving algorithm attempts to
    minimize this quantity as part of one of its refinement steps.
    """

    H_in = csiso_func(r[i_in]) / Omega_func(r[i_in])
    H_out = csiso_func(r[i_out]) / Omega_func(r[i_out])

    i_cut_in = np.argmin(np.abs(r - (r[i_in] + H_in)))
    i_cut_out = np.argmin(np.abs(r - (r[i_out] - H_out)))

    i_bounds = np.concatenate((np.arange(i_in, i_cut_in), np.arange(i_cut_out, i_out)))

    dphidr = dphidr_func(f)

    phase_grad_error = np.empty(len(r))

    for i in i_bounds:

        phase_grad_error[i] = np.real(dphidr[i] - k_func(r[i], m))

    error_in = np.max(phase_grad_error[i_in:i_cut_in]) - np.min(
        phase_grad_error[i_in:i_cut_in]
    )
    error_out = np.max(phase_grad_error[i_cut_out:i_out]) - np.min(
        phase_grad_error[i_cut_out:i_out]
    )

    return error_in, error_out


def get_boundaries(f, thresh=1.0e-4):

    """
    Find new boundary locations for mode solution "f", by finding where its magnitude
    has dropped below "thresh" (default value 1e-4) of its value at r = 1. If the
    magnitude never drops below this threshold, the boundaries are simply the first and
    last points of the grid.
    """

    max = np.abs(f[i1, 0])

    i_in = 0

    for i in range(i1, -1, -1):

        if np.abs(f[i, 0]) < thresh * max:

            i_in = i
            break

    i_out = len(f) - 1

    for i in range(i1, len(f)):

        if np.abs(f[i, 0]) < thresh * max:

            i_out = i
            break

    return i_in, i_out


def get_wkb_bc(i_in, i_out, m):

    """
    Get the WKB outgoing wave coefficients (see the function "wkb_func") at the inner
    and outer boundaries, specified by the radial grid indices "i_in" and "i_out", for
    azimuthal number "m".
    """

    C_i = wkb_func(r[i_in], m)

    if m == 1:

        C_i -= 2.0j * k_func(r[i_in], m)

    C_o = wkb_func(r[i_out], m)

    return C_i, C_o


def solve_m_mode(m):

    """
    The main function for solving for planet-disk modes. Takes the azimuthal number "m"
    of the mode to be solved as its single argument. Returns the mode solution f,
    an array of complex numbers with shape ("r_steps",2) ("r_steps" is the grid size,
    defined in the parameter file), where f[i,0] and f[i,1] correspond to dh and its
    radial derivative dh' (both are needed to calculate the other fluid perturbation
    variables), respectively, at the grid point r[i].
    """

    start = time.time()

    i_in = 0
    i_out = len(r) - 1

    # Find two linearly independent homogeneous solutions and the inhomogeneous solution

    f_h1 = integrate_ode(r, [1.0, 0.0], m, 0)
    f_h2 = integrate_ode(r, [0.0, 1.0], m, 0)
    f_ih = integrate_ode(r, [0.0, 0.0], m, 1)

    # Find the coefficients of the homoegenous solutions needed to satisfy the outgoing
    # wave boundary conditions at r_in and r_out.

    C_i, C_o = get_wkb_bc(i_in, i_out, m)
    a1, a2 = homogeneous_coefficients(f_h1, f_h2, f_ih, C_i, C_o, i_in, i_out)

    f_m = a1 * f_h1 + a2 * f_h2 + f_ih

    # Get new boundary locations (see the function "get_boundaries") and find a new
    # solution satisfying outgoing wave conditions at the new boundaries. This prevents
    # rapidly decaying solutions from blowing up (due to contamination by an
    # exponentially growing solution).

    i_in, i_out = get_boundaries(f_m)
    C_i, C_o = get_wkb_bc(i_in, i_out, m)
    a1, a2 = homogeneous_coefficients(f_h1, f_h2, f_ih, C_i, C_o, i_in, i_out)

    # Reintegrate the the inhomogeneous solution starting with the (more or less)
    # correct values of dh and dh' at the corotation radius (rather than random
    # guesses). Find a new combination of the two homogeneous solutions and this new
    # inhomogeneous solution that satisfies the boundary conditions. This leads to a
    # more accurate solution. Perform this refinement step "reint_steps" times.

    for i in range(reint_steps):

        i_co = np.argmin(np.abs(m * Omega_func(r) - omega_func(m)))
        f_ih = integrate_ode(r, f_m[i_co], m, 1)
        a1, a2 = homogeneous_coefficients(f_h1, f_h2, f_ih, C_i, C_o, i_in, i_out)

        f_m = a1 * f_h1 + a2 * f_h2 + f_ih

        i_in, i_out = get_boundaries(f_m)
        C_i, C_o = get_wkb_bc(i_in, i_out, m)
        a1, a2 = homogeneous_coefficients(f_h1, f_h2, f_ih, C_i, C_o, i_in, i_out)

    # Compute the phase error (difference between the mode phase derivative and the WKB
    # wavenumber) near the boundaries, and attempt to minimize it by finding new,
    # slightly modified outgoing wave boundary conditions (i.e., coefficients C in
    # dh' = C*dh) and finding a new (full homogeneous + inhomogeneous) solution
    # satisfying these conditions. Repeat this process "phgrad_steps" times.

    error = phase_error_diagnostic(f_m, m, i_in, i_out)

    for i in range(phgrad_steps):

        delta = 1.0e-6

        C_ip = C_i * (1.0 + delta)
        C_im = C_i * (1.0 - delta)
        C_op = C_o * (1.0 + delta)
        C_om = C_o * (1.0 - delta)

        a1_ip, a2_ip = homogeneous_coefficients(
            f_h1, f_h2, f_ih, C_ip, C_o, i_in, i_out
        )
        f_m_ip = a1_ip * f_h1 + a2_ip * f_h2 + f_ih

        a1_im, a2_im = homogeneous_coefficients(
            f_h1, f_h2, f_ih, C_im, C_o, i_in, i_out
        )
        f_m_im = a1_im * f_h1 + a2_im * f_h2 + f_ih

        a1_op, a2_op = homogeneous_coefficients(
            f_h1, f_h2, f_ih, C_i, C_op, i_in, i_out
        )
        f_m_op = a1_op * f_h1 + a2_op * f_h2 + f_ih

        a1_om, a2_om = homogeneous_coefficients(
            f_h1, f_h2, f_ih, C_i, C_om, i_in, i_out
        )
        f_m_om = a1_om * f_h1 + a2_om * f_h2 + f_ih

        error_ip = phase_error_diagnostic(f_m_ip, m, i_in, i_out)
        error_im = phase_error_diagnostic(f_m_im, m, i_in, i_out)
        error_op = phase_error_diagnostic(f_m_op, m, i_in, i_out)
        error_om = phase_error_diagnostic(f_m_om, m, i_in, i_out)

        jac = np.empty((2, 2), complex)
        jac[0, 0] = (error_ip[0] - error_im[0]) / (C_ip - C_im)
        jac[0, 1] = (error_op[0] - error_om[0]) / (C_op - C_om)
        jac[1, 0] = (error_ip[1] - error_im[1]) / (C_ip - C_im)
        jac[1, 1] = (error_op[1] - error_om[1]) / (C_op - C_om)
        jac_det = jac[0, 0] * jac[1, 1] - jac[0, 1] * jac[1, 0]

        if m == 1:

            C_o -= error[1] / jac[1, 1]

        else:

            C_i -= (jac[1, 1] * error[0] - jac[0, 1] * error[1]) / jac_det
            C_o -= (-jac[1, 0] * error[0] + jac[0, 0] * error[1]) / jac_det

        a1, a2 = homogeneous_coefficients(f_h1, f_h2, f_ih, C_i, C_o, i_in, i_out)
        f_m = a1 * f_h1 + a2 * f_h2 + f_ih

    # Set the solution outside of the new boundaries to zero.

    f_m[:i_in, :] = 0.0
    f_m[i_out + 1 :, :] = 0.0

    end = time.time()

    print("mode", m, "solution found, time taken = ", end - start)

    return f_m


def vars_from_enthalpy(f, r, m):

    """
    Compute the perturbations \delta u_r, \delta u_\phi and \delta\Sigma (the radial
    velocity, azimuthal velocity, and surface density perturbation), given "f", the mode
    solution for the enthalpy described produced by "solve_m_mode", the radial grid "r",
    and the azimuthal mode number "m". Returns the three complex perturbation variables,
    each is a complex array of length "r_steps".
    """

    omega = omega_func(m)

    u_r, u_phi, sigma = (np.empty(len(r), complex) for i in range(3))

    for i in range(len(r)):

        Omega = Omega_func(r[i])
        kappa = kappa_func(r[i])

        omega_tilde = omega - m * Omega
        omega_tilde += 1.0j * pole_disp
        tc = tc_func(r[i])
        beta_tilde = omega_tilde * tc

        Nr2 = Nr2_func(r[i])

        LSinv = LSinv_func(r[i])
        LTinv = LTinv_func(r[i])

        Leffinv = (LTinv - 1.0j * gamma * beta_tilde * LSinv) / (
            1.0 - 1.0j * gamma * beta_tilde
        )

        D = D_func(r[i], m)
        Dc = D - 1.0j * gamma * beta_tilde / (1.0 - 1.0j * gamma * beta_tilde) * Nr2

        Phi = Phi_m_func(r[i], m)
        dPhidr = dPhidr_m_func(r[i], m)

        csadi = csadi_func(r[i])

        u_r[i] = (
            1.0j
            / Dc
            * (
                omega_tilde * (f[i, 1] + dPhidr)
                - 2.0 * m * Omega / r[i] * (f[i, 0] + Phi)
                - Leffinv * omega_tilde * f[i, 0]
            )
        )
        u_phi[i] = (
            1.0
            / Dc
            * (
                kappa ** 2.0 / (2.0 * Omega) * (f[i, 1] + dPhidr)
                - m * (kappa ** 2.0 - Dc) / (r[i] * omega_tilde) * (f[i, 0] + Phi)
                - Leffinv * kappa ** 2.0 / (2.0 * Omega) * f[i, 0]
            )
        )

        sigma[i] = (
            Sigma_func(r[i])
            / (1.0 / gamma - 1.0j * beta_tilde)
            * ((1.0 - 1.0j * beta_tilde) * f[i, 0] / csadi ** 2.0 + LSinv * tc * u_r[i])
        )

    return u_r, u_phi, sigma
