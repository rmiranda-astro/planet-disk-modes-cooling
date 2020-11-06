import numpy as np
from math import pi
from main import *
import sys

# Read the number of radial grid points "nr_new" for the torque density/angular momentum
# flux profiles, the maximum azimuthal number "m_max" to be used to construct them, and
# a file extension "file_ext" to append to the output file from the command line
# arguments.

nr_new = int(sys.argv[1])
m_max = int(sys.argv[2])
file_ext = sys.argv[3]

scale = 1.0 / h0 ** 3.0

r = np.loadtxt("rmodes.out")
nr = len(r)

# Specify the new radial grid.

r_2 = np.logspace(np.log10(r[0]), np.log10(r[-1]), nr_new)

hr, hi, hpr, hpi = (np.empty((m_max, nr)) for i in range(4))

f_m = np.empty((m_max, nr, 2), complex)
f_m_2 = np.empty((m_max, len(r_2), 2), complex)
ur_m_2, uphi_m_2, sigma_m_2 = (np.empty((m_max, len(r_2)), complex) for i in range(3))

# Read the enthalpy perturbation from the mode files for each azimuthal number,
# interpolate them to the new radial grid, and use them to compute the surface density
# and velocity perturbations that will be used to compute the torque density and angular
# momentum flux.

for i in range(m_max):

    m = i + 1

    print("reading/synthesizing m =", m, "...")

    hr[i], hi[i], hpr[i], hpi[i] = np.loadtxt("m_%d.out" % m, unpack=True)

    f_m[i, :, 0] = hr[i] + 1.0j * hi[i]
    f_m[i, :, 1] = hpr[i] + 1.0j * hpi[i]

    f_m_2[i, :, 0] = np.interp(r_2, r, np.real(f_m[i, :, 0])) + 1.0j * np.interp(
        r_2, r, np.imag(f_m[i, :, 0])
    )
    f_m_2[i, :, 1] = np.interp(r_2, r, np.real(f_m[i, :, 1])) + 1.0j * np.interp(
        r_2, r, np.imag(f_m[i, :, 1])
    )

    ur_m_2[i], uphi_m_2[i], sigma_m_2[i] = vars_from_enthalpy(f_m_2[i], r_2, m)

dTdr = np.zeros(len(r_2))
amf = np.zeros(len(r_2))

T_tot = 0.0

# Compute the torque density "dTdr" and angular momentum flux "amf" using the surface
# density and velocity perturbations.

for i in range(len(r_2)):

    Sigma = Sigma_func(r_2[i])

    Omega = Omega_func(r_2[i])
    Omega_p = 1.0

    for j in range(m_max):

        m = j + 1

        omega_tilde = m * (Omega_p - Omega)

        Phi = Phi_m_func(r_2[i], m)

        dTdr[i] += -pi * r_2[i] * m * Phi * np.imag(sigma_m_2[j, i])

        amf[i] += (
            pi * Sigma * r_2[i] ** 2 * np.real(ur_m_2[j, i] * np.conj(uphi_m_2[j, i]))
        )

    if i > 0:
        T_tot += dTdr[i] * (r_2[i] - r_2[i - 1])

# Compute the integrated torque "T_int" (the torque density integrated from the
# corotation radius).

T_int = np.zeros(len(r_2))
i1 = np.argmin(np.abs(Omega_func(r_2) - 1.0))

for i in range(i1 + 1, len(r_2)):

    T_int[i] = T_int[i - 1] + dTdr[i] * (r_2[i] - r_2[i - 1])

for i in range(i1 - 1, -1, -1):

    T_int[i] = T_int[i + 1] - dTdr[i] * (r_2[i + 1] - r_2[i])

dTdr_amf_file = open("dTdr-amf_%s.out" % file_ext, "w")

# Write the output file, which has four columns: the disk radius, torque density,
# integrated torque, and angular momentum flux (scaled by their natural units).

for i in range(len(r_2)):

    dTdr_amf_file.write(
        "%e  %e  %e  %e\n" % (r_2[i], dTdr[i] / scale, T_int[i] / scale, amf[i] / scale)
    )

dTdr_amf_file.close()
