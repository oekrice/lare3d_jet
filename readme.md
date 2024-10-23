<<<<<<< HEAD
# lare3d_jet
lare3d for the Pariat jet model

For the full LARE readme look for readme_lare.md or equivalent, which explains more about how the code actually works. This readme should explain how to get the jet model to run.

Steps to install and run:

Navigate to a suitable place and clone this repository:
```
  git clone https://github.com/oekrice/lare3d_jet/
```
which should download the files and things. The base folder is mainly full of auxiliary python scripts to move things around and plot rather than the base code itself, which is in the folder 'src'. 

# To compile the code 

This requires the use of one of the Makefiles. Two are provided: one which works on my local linux machine using gfortran, and the other which works on an . They are really quite different and it is likely that they would require some tweaking to be used elsewhere. Pick whichever one is more appropriate and rename:

```
  scp Makefile_gfortran Makefile
```
or 

```
  scp Makefile_intel Makefile
```
Then update the module locations for MPIF90 and the NETCDF libraries (near the top of the file). On the HPCs I've used the locations of such things can be obtained using 

```
  module show netcdf-fortran/intel-2020.4/4.5.4
```
or equivalent, but finding them can be tricky.

If this has been done correctly, the code can then be compiled merely with 

```
make
```

# To run the code 

Before actually running the Fortran code, there is some setup to be done. The python wrapper 'run.py' is used to define various variables (resolutions, boundaries, number of plots etc.), which are then saved to a 'parameters' file and read in to Fortran in real time. The initial condition (PFSS) is also called in this file. I recommend running this on a login node before running the bash script to get Fortran actually going. As is, the code is set up on a 64^3 grid which should produce an eruption and not take too long. 

```
python run.py 0
```

where the 0 is the `run number' -  a functionality that doesn't yet work but would be used for multiple runs.

LARE is then run by calling mpiexec or equivalent. Some examples:

```
/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/lare3d 0
```
```
mpirun ./bin/lare3d 0
```
depending on where your MPI is stored. `-np` is the number of processes, which is automatically determined by the machine on an HPC. The `0` is again the run number.

This should print out various information, and save files to `./Data/', or wherever else has been specified.

A common error is ` No such file or directory NetCDF: Not a valid ID', which is usually because the output filenames have gone wrong.

# To plot things

The plotting wrapper is called `plotslice.py', which can be run just by using 
```
python plotslice.py
```
This plots some colourmaps of cuts of things like the density and magnetic field, and also uses the field line tracer from fltrace.py for some visualisation. These are saved in `./plots'. There is an optional pyvista equivalent, which can be activated by setting `plot_vista = True' on line 163, and uncommenting the relevant imports at the top of fltrace.py. But pyvista doesn't work on every machine so it's left off by default.




