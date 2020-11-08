This code is used to solve the linear perturbation equations for planet-disk interaction with constant-beta cooling, as presented in [Miranda & Rafikov (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...892...65M/abstract), hereafter MR20, using the numerical solution method described in the Appendix of [Miranda & Rafikov (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...875...37M/abstract), hereafter MR19.

# Mode Solutions

[main.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/main.py) contains the main machinery for finding mode solutions of the master equation (equations (9) & (24)-(26) of MR20) for the enthalpy perturbation \delta h. It contains two key functions:

* `solve_m_mode(m)`, which takes the azimuthal number m of the mode to be solved as its single argument, computes and returns the mode solution f, an array of complex numbers with shape (r_steps, 2) (where r_steps is the grid size defined in the parameter file), where f[i,0] and f[i,1] correspond to \delta h and its radial derivative (which is needed to calculate the other fluid perturbation variables), respectively, at the grid point r[i].

* `vars_from_enthalpy(f,r,m)` computes the surface density perturbation \delta\Sigma and the velocity perturbations \delta u_r and \delta u_\phi, given a mode solution f for the enthalpy produced by `solve_m_mode`, its corresponding grid array r, and the azimuthal mode number m. The three fluid variables are returned as complex arrays of length r_steps.

The parameters for the disk model, numerical grid, and some options for the mode solver are specified in [parameters.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/parameters.py):

* Disk Model Parameters

	* h0: Disk aspect ratio H/r at the planetary orbital radius r = r_p (or r = 1 in code units)

	* q: Temperature power law index (T follows an r^-q power law)

	* p: Surface density power law index (\Sigma follows an r^-p power law)

	* gamma: The adiabatic index

	* beta: Dimensionless cooling timescale \Omega t_c. Note that beta can safely be set to exactly zero, which will result in the correct behavior for a locally isothermal disk.

	* soft: Softening length for the planetary potential (as a fraction of r_p)

* Grid Parameters

	* r_in: Inner grid radius (in terms of r_p)

	* r_out: Outer grid radiuss

	* r_steps: Size of the logarithmic grid (number of grid points)

* Solver Options

	* indirect: Turns on (1) or off (0) the indirect potential associated with the motion of the star around the center of mass of the star + planet system. Including the indirect potential only slightly modifies the m = 1 mode, which typically has a small amplitude compared to modes with larger azimuthal numbers. 

	* pole_disp: Amount by which the corotation singularity is displaced from the real line. Must be positive and should be as small as possible. A good fiducial value is 10^-6. If the integrator used in the mode solver returns an error (step size too small), it may be because pole_disp is too small, in which case a somewhat larger value should be used.

	* reint_steps: Number of times to re-integrate the inhomogeneous solution to refine the the coefficients of the homogeneous solutions in the full mode solution (see MR19). This should always be at least 1 to ensure accurate solutions. Increasing it to 2 may marginally improve the solution accuracy, and using a larger value has a negligible effect.

	* phgrad_steps: Number of iterations for the phase gradient refinement step (see Fig. 15 of MR19). This should usually be set to 1.

	* m1_switch: Switches on (1)/off (0) the use of a more accurate, but much slower integrator when solving for the m = 1 mode only. This solves an issue in which the m = 1 mode blows up in the inner disk, although this problem is rare in the current version of the code, due to the use of adaptive boundary locations. This can usually be set to 0, but if the resulting m = 1 mode doesn't look right (i.e., doesn't properly decay towards the inner disk), then setting it to 1 should solve the problem.

Running [find_m_modes.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/find_m_modes.py) calls `solve_m_mode` to produce mode solutions for values of the azimuthal number m up to a specified maximum value m_max. The maximum azimuthal number and number of parallel processors to use for the calculation are specified as command line arguments. For example, to solve for m_max = 80 modes using 28 processors, use

`python find_m_modes.py 80 28`

When a solution for each mode is found, message is printed to the terminal, indicating the time taken to compute it. When the script is finished, the following files are written:

* rmodes.out: Contains the radial coordinates of the grid points used for the mode solutions

* m_1.out, m_2.out, ... m_{m_max}.out, containing the enthalpy perturbation for each mode number, written as four columns

  * Real part of \delta h
  * Imaginary part of \delta h
  * Real part of the derivative of \delta h
  * Imaginary part of the derivative of \delta h

[master_equation_cool.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/master_equation_cool.py), which is called during the mode solution routine, computes the real and imaginary parts of the ODE coefficients C_1, C_0, and \Psi, defined in equation (9) of MR20 (taking into account the displacement of the corotation singularity from the real line). The master equation can be modified by modifying this file.

# 2D Maps of the Fluid Variables

Running [synthesize_2d.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/synthesize_2d.py) produces 2D maps (in polar coordinates) of the fractional surface density perturbation \delta\Sigma/\Sigma and the velocity perturbations \delta u_r and \delta u_\phi, using the enthalpy mode solutions. The mode files (m_1.out, m_2.out, etc), found using find_m_modes.py, are required to produce the maps.

Several parameters are specified as command line arguments: the size of the grid used for the maps, the number of modes to be used in their construction, and a file extension that is added to the end of the output filenames (containing, e.g., the disk parameters). For example, to produce maps with N_r x N_phi = 1024 x 2048 grid points, using m_max = 80 modes with the extension "file_ext", use

`python synthesize_2d.py 1024 2048 80 file_ext`

The following files are produced:

* r2d_(file_ext).out: The radial grid for the 2D maps

* phi2d_(file_ext).out: The azimuthal grid for the maps

* sigma2d_(file_ext).out: 2D map of \delta\Sigma/\Sigma, scaled by M_p/M_th (where M_th is the thermal mass, see equation (35) of MR20)

* ur2d_(file_ext).out: 2D map of \delta u_r, scaled by c_s,p M_p/M_th (where c_s,p is the isothermal sound speed at r_p)

* uphi2d_(file_ext).out: 2D map of \delta u_\phi (scaled the same as the above)

[plot_2d.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/plot_2d.py) gives an example of a simple routine to read and plot these files.

# Torque Density & Angular Momentum Flux

[dTdr_amf.py](https://github.com/rmiranda-astro/planet-disk-modes-cooling/blob/main/dTdr_amf.py) computes the torque density (see equation (36) of MR20), integrated torque, and angular momentum flux (see equation (39) of MR20), using the enthalpy mode solutions. Here again the mode files (m_1.out, m_2.out, etc) are needed. The size of the new radial grid for these quantities, the maximum mode number, and a file extension (as in the 2D maps) are specified as command line arguments, e.g., 

`python dTdr_amf.py 10000 80 file_ext`

This produces the file dTdr-amf_(file_ext).out, which has four columns:

* Radial position
* Torque density dT/dr (scaled by F_J,0/r_p, see equation (64) of MR20)
* Integrated torque (dT/dr integrated from the corotation radius; scaled by F_J,0)
* Angular momentum flux F_J (scaled by F_J,0)
