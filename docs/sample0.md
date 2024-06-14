# Quick Start

This document describes a working prototype of the Surface Flux
Transport model (SFT). The SFT was developed to provide a computationaly
inexpensive solver to explore the parameter space to better understand
the evolution of solar surface magnetic field evolution in response to
the changes in source functions and transport parameters. In its current
form the SFT is numerically stable and provides a user-friendly test-bed
to understand how the properties of sunspots, advection profiles and
diffusivity impact the global solar surface magnetic field evolution.

The core of the SFT are written in Fortran 90. As the equations are
simplified to evolve in one dimensions, presently the code does not
implement parallel communications libraries. The SFT creates a single
executable which can be compiled with netCDF libraries and run on any
computer as long as the GNU compiler and netCDF libraries are properly
installed.

## Acknowledgments

This numerical implementation of the model follows [Yeates (2020)](https://doi.org/10.1007/s11207-020-01688-y). The
first version of the SFT1D was developed at the Institute for Astronomy
(IFA) of the University of Hawaii.

## System Requirements

In order to install and run the SFT 1D the following minimum system
requirements apply.

-   The SFT runs only under the UNIX/Linux operating systems. This now
    includes Macintosh system 10.x because it is based on BSD UNIX. It
    does not run under any Microsoft Windows operating system.

-   A GNU FORTRAN compiler must be installed.

-   The file writing subroutine use netCDF output. For these the
    serial/parallel version of the netCDF library has to be installed.

-   In order to generate the documentation, LaTex has to be installed.

As the code solves the magnetic field evolution in only one dimension,
it does not need huge computing power or higher memory load. It can run
on a single processor with a nominal RAM memory.

In addition to the above requirements, the SFT output is typically
visualized using Python. Other visualization packages may also be used,
but the output file formats should be suported by those visualization
softwares.

## A Brief Description of the SFT 1D code

The distribution in the form of the compressed tar image includes the
SFT source code. The top level directory contains the following
subdirectories:

-   `doc` - the documentation directory

-   `bipole_file` - the bipole properties to generate the source terms
    for the model

-   `input_files` - initial magnetic field configuration for the SFT
    code, typically a synoptic magnetogram interpolated into the
    targeted sine-latitude -- longitude grid

-   `src` - all the fortran routines for building the executable

-   `plots` - some example plots of a standard run

and the following files

-   `initial_Parameters.nml` - user input of the initial parameter
    values for the variables

-   `Makefile` - the main makefile

-   `makefile.git` - makefile to build and test the code on github with
    continous integration (CI)

-   `hmi_polar_field.p` - pickle file contaiing HMI observed polar field
    obtained from JSOC webpage

-   `install_netcdf.sh` - shell script to download and install necessary
    netCDF libraries on a linux system

There are eight fortran files in the `src` directory which cotains the
source code.

-   `evolSFT.f90` - a short instruction on installation and usage

-   `flows.f90` - the analytical form of the flow profile

-   `gridSFT.f90` - to construct the grid of the solver

-   `init_condition.f90` - code to set-up the initial conditions

-   `output.f90` - subroutine to write the files in ASCII format

-   `main.f90` - the main code

-   `variables.f90` - definition of the variables to be used in the code

-   `write_data.f90` - code to write the output data in netCDF format.

## General Hints

#### Getting help with the Makefile

You can find all the possible targets that can be built by typing

    make help

#### Compiling the code

Before going to compile the code it is essential to modify a few things
in the `Makefile`. Currently, the fortran compiler and the paths of the
netCDF libraries are set to this path.

    # Set FORTRAN90 compiler
    FC = gfortran

    # Location of files for netcdf library
    NETCDF = -I/usr/local/netcdf/include/
    NETCDFLIB = -L/usr/local/netcdf/lib/ -lnetcdff

You need to change the paths and the compiler options according to your
system. Once this is correctly specified, compile the code by typing
`make`. This will create two additional sub-directories; `bin, obj`. The
executable `SFT_1D` will be created in the bin directory by compiling
codes from `src` directory, if there are no errors in the compilation
process.

#### Running the code

Surface flux transport model requires two major input from the user.
First being the meridional flow profile and magnetic diffusivity
($\eta$). And the second being the source functions i.e., bipole
properties which are added to the evolution at the time when they appear
on the solar photosphere.

-   Bipolar Magnetic Region (BMRs) properties are pre-written to an
    ASCII file and saved in `bipole_file/all_bmrs.txt`.

-   Meridional flow profile is defined in the subroutine `MC_flow`
    contained in `flows.f90` file.

-   Magnetic diffusivity ($\eta$) is defined in units of km$^2$/s in the
    `initial_Parameters.nml` file.

These parameters are passed to the code using a `namelist` file
`initial_Parameters.nml`. This will create the necessary subdirectories
for saving the output and execute the code based on the selected
parameters. After compiling the code, type

    >./bin/SFT\_1D initial\_Parameters.nml

to run the simulation.

There is a sample plotting file (`plot_results.py`) is also distributed
with this version of the code along with some supporting files. If all
the necessary python packages are available, this script will generate
the butterfly diagram, comparative plot between the HMI polar field and
the SFT polar field and the dipole moment. With the current set of
initial conditions, the example of these plots are provided with this
version of the code.
