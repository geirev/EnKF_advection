#EnKF_advection
EnKF implementation with advection equation as used in Evensen (2009).

To compile you must also clone the EnKF_Analysis and the EnKF_sampling repositories.

In EnKF_analysis/lib/makefile define the build catalog to be the same as for EnKF_advection and make.

In EnKF_sampling/lib/makefile define the build catalog to be the same as for EnKF_advection and make.

In EnKF_advection/source run make

Copy the infile.in file to a catalog where you wish to run advect.lin.

