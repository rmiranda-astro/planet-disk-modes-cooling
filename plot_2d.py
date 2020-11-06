import numpy as np
from math import pi
import matplotlib.pyplot as plt
import sys

# Read the file extension (as specified in "synthesize_2d_map_cool") of the 2D maps/grid
# files to be read.

file_ext = sys.argv[1]

# Read the grid files and a 2D surface density map given a file extension. Return the 2D
# radial and azimuthal grids "rr" and "phiphi", and the surface density map "f2d" as
# 2D arrays.

def read_2d(ext):

    r = np.loadtxt("r2d_%s.out" % ext)
    phi = np.loadtxt("phi2d_%s.out" % ext)
    f2d = np.loadtxt("sigma2d_%s.out" % ext).reshape(len(r), len(phi))

    rr, phiphi = np.meshgrid(r, phi)

    return rr, phiphi, f2d

# Read and plot a 2D surface density map.

rr_test, phiphi_test, f2d_test = read_2d(file_ext)

cblim = 2.0

plt.axes().set_aspect("equal")
plt.xlim(-1, 1)
plt.ylim(np.log10(rr_test[0, 0]), np.log10(rr_test[-1, -1]))
map = plt.pcolormesh(
    phiphi_test / pi,
    np.log10(rr_test),
    f2d_test.T,
    vmin=-cblim,
    vmax=cblim,
    rasterized=True,
)
plt.colorbar(map)
plt.savefig("out.pdf", dpi=300)
