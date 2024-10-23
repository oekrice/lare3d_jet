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

Then update the module locations for MPIF90 and the NETCDF libraries (near the top of the file). On the HPCs I've used the locations of such things can be obtained using 

```
  module show netcdf-fortran/intel-2020.4/4.5.4
```
or equivalent, but finding them can be tricky.
