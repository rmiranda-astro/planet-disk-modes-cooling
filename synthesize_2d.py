import numpy as np
from math import pi
from main import *
import sys

# Read the number of radial and azimuthal grid points for the 2D maps "nr_new" and
# "nphi", the maximum azimuthal number "m_max" to use in the 2D synthesis, and a
# file extension "file_ext" to append to the output files from the command line
# arguments.

nr_new = int(sys.argv[1])
nphi = int(sys.argv[2])
m_max = int(sys.argv[3])
file_ext = sys.argv[4]

r = np.loadtxt("rmodes.out")
nr = len(r)

r_in_new = r[0]
r_out_new = r[-1]

# Specify the new radial grid and the azimuthal grid.

r_new = np.logspace(np.log10(r_in_new), np.log10(r_out_new), nr_new)

phi = np.linspace(-pi, pi, nphi)

hr, hi, hpr, hpi = (np.empty((m_max, nr)) for i in range(4))

f_m = np.empty((m_max, nr, 2), complex)

# Read the enthalpy perturbations from the mode files for each azimuthal number.

for i in range(m_max):

    m = i + 1

    print("reading m =", m, "...")

    hr[i], hi[i], hpr[i], hpi[i] = np.loadtxt("m_%d.out" % m, unpack=True)

    f_m[i, :, 0] = hr[i] + 1.0j * hi[i]
    f_m[i, :, 1] = hpr[i] + 1.0j * hpi[i]

sigma2d, ur2d, uphi2d = (np.zeros((nr_new, nphi), complex) for i in range(3))

# Interpolate the enthalpy perturbations from the original grid that the modes were
# solved on to the new grid to be used for the 2D maps, compute the velocity and surface
# density perturbations from the enthalpy perturbation, and construct their 2D maps.

for i in range(m_max):

    m = i + 1

    print("synthesizing m = ", m, "...")

    htemp = np.empty((nr_new, 2), complex)

    htemp[:, 0] = np.interp(r_new, r, np.real(f_m[i, :, 0])) + 1.0j * np.interp(
        r_new, r, np.imag(f_m[i, :, 0])
    )
    htemp[:, 1] = np.interp(r_new, r, np.real(f_m[i, :, 1])) + 1.0j * np.interp(
        r_new, r, np.imag(f_m[i, :, 1])
    )

    ur_m, uphi_m, sigma_m = vars_from_enthalpy(htemp, r_new, m)

    sigma_m /= Sigma_func(r_new)

    for j in range(nr_new):

        sigma2d[j] += sigma_m[j] * np.exp(1.0j * m * phi)
        ur2d[j] += ur_m[j] * np.exp(1.0j * m * phi)
        uphi2d[j] += uphi_m[j] * np.exp(1.0j * m * phi)

# Scale the velocity and surface density maps by their natural units.

sigma2d *= h0 ** 3.0
ur2d *= h0 ** 2.0
uphi2d *= h0 ** 2.0

# Write the output files, the radial and azimiuthal grids for the 2D maps, and the maps
# themselves for the velocity and surface density perturbations.

print("writing output...")

r_out_file = open("r2d_%s.out" % file_ext, "w")
phi_out_file = open("phi2d_%s.out" % file_ext, "w")
sigma2d_out_file = open("sigma2d_%s.out" % file_ext, "w")
ur2d_out_file = open("ur2d_%s.out" % file_ext, "w")
uphi2d_out_file = open("uphi2d_%s.out" % file_ext, "w")

for i in range(nr_new):

    r_out_file.write("%e\n" % r_new[i])

    for j in range(len(phi)):

        if i == 0:

            phi_out_file.write("%e\n" % phi[j])

        sigma2d_out_file.write("%e\n" % np.real(sigma2d[i, j]))
        ur2d_out_file.write("%e\n" % np.real(ur2d[i, j]))
        uphi2d_out_file.write("%e\n" % np.real(uphi2d[i, j]))

r_out_file.close()
phi_out_file.close()
sigma2d_out_file.close()
ur2d_out_file.close()
uphi2d_out_file.close()
