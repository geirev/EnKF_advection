EnKF implementation with advection equation as used in Evensen (2009).

To compile you must also clone the EnKF_Analysis repository.
Then in lib/makefile define the build catalog to be the same as for EnKF_advection.
make to install analysis library and .mod files.

Now you should be able to compile the EnKF_advection distribution in source.
Copy the infile.in file to a catalog where you wish to run advect.lin.

