from main import *
import numpy as np
from math import pi
import multiprocessing
import sys

# Read the maximum azimuthal number "m_max" and the number of processors "nproc" from
# the command line arguments.

m_max = int(sys.argv[1])
nproc = int(sys.argv[2])

f_m = np.empty((m_max, len(r), 2), complex)

# If only using one processor, solve for each mode from m = 1 to m = m_max sequentially.

if nproc == 1:

    for i in range(m_max):
    
        m = i + 1
        f_m[i] = solve_m_mode(m)
        
# Use multiprocessing to solve for the modes in parallel if the specified number of
# processors is greater than one.
        
else:

    pool = multiprocessing.Pool(processes=nproc)

    m_array = np.arange(1, m_max + 1)
    modes = pool.map(solve_m_mode, m_array)

    for i in range(m_max):
    
        f_m[i] = modes[i]
        
# Write the mode output files, each of which has four columns, giving Re(dh), Im(dh),
# Re(dh') and Im(dh') at each grid point.

for i in range(m_max):

    m = i + 1

    print("writing m =", m, "file...")

    out_file = open("m_%d.out" % m, "w")
    
    for j in range(len(r)):
    
        out_file.write("%e  " % np.real(f_m[i, j, 0]))
        out_file.write("%e  " % np.imag(f_m[i, j, 0]))
        out_file.write("%e  " % np.real(f_m[i, j, 1]))
        out_file.write("%e  " % np.imag(f_m[i, j, 1]))
        out_file.write("\n")
        
    out_file.close()

out_file = open("rmodes.out", "w")

# Write the grid file that gives the radial coordinates of each grid point.

for i in range(len(r)):

    out_file.write("%e\n" % r[i])
    
out_file.close()
