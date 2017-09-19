## General purpose 3D lattice Boltzmann code

SunlightLB is an open-source 3D lattice Boltzmann code which can be
used to solve a variety of hydrodynamics problems, including some
passive scalar transport problems. Features include: 

* C-based code ;
* scriptable ;
* D3Q19 lattice ;
* twin-relaxation-time (stress relaxation) scheme ;
* moving objects ;
* link-bounce-back boundary conditions ;
* zero and non-zero Reynolds number flows.

SunlightLB is implemented as a library of C functions.  Scripting
language support is enabled by a [SWIG](http://www.swig.org/)
interface file, which allows, for example, SunlightLB to be used as a
perl extension module.  Examples of the use of SunlightLB are provided
in both C, perl, and python.  SunlightLB implements a standard lattice
Boltzmann algorithm for three dimensional simulations. It uses a D3Q19
lattice with a twin relaxation time scheme [1-3].  Objects, possibly
moving, are included by a link bounce-back method [2].  This enables
SunlightLB to solve a variety of hydrodynamics problems such as the
computation of flows through pore spaces, the computation of
resistance matrices for colloidal hydrodynamics problems, and so
on. Both zero Reynolds number flows, and non-zero Reynolds number
flows, can be solved.  In addition, passive scalar transport is
implemented on top of the lattice Boltzmann scheme via a
tagged-particle propagation algorithm, with a variety of boundary
conditions. This allows simulation of a variety of
reaction-advection-diffusion problems, such as a passive scalar
adsorbing in a porous material in the presence of a flow (deep-bed
filtration).  The algorithms underlying SunlightLB are published in
the following literature :

[1] Qian, d'Humieres and Lallemand, Europhys. Lett. **17**, 479
(1992).

[2] Ladd, [J. Fluid Mech. **271**, 285 (1994)](https://doi.org/10.1017/S0022112094001771);
[J. Fluid Mech. **271**, 311 (1994)](https://doi.org/10.1017/S0022112094001783). 

[3] Behrend, Harris and Warren, [Phys. Rev. E **50**, 4586 (1994)](https://doi.org/10.1103/PhysRevE.50.4586).

### Installation Notes

These notes are in the process of being updated!

The basic lattice Boltzmann code is implemented as a set of library
routines, coded in vanilla C.  The Makefile is configured to build
these into a shared object (so) library.  A
[SWIG](http://www.swig.org/) (Simplified Wrapper and Interface
Generator) interface file allows access for scripting languages like
perl or python, by building extension modules.  Note that the
installation steps below may require root priviledges.

To build and install SunlightLB as a shared object library.

* Go into the `src` directory.
* Edit the top part of the `Makefile` (don't worry about the perl- or python-related variables at this point).
* Type `make` to build the `libsunlightlb.so*` library files.
* Type `make install` to install the header and library files.
* Go into the `examples` directory.
* Edit the top part of the `Makefile` as above.
* Type `make sphere` to compile the simplest example.</li>
* Typing `./sphere` runs the code.

Note that the Makefile for the examples assumes that the header file
and library will be found by the compiler (eg they are in
`/usr/local/include` and `/usr/local/lib` respectively). If this is
not the case, you may have to specify the locations using the `-I` and
`-L` compiler options.  In addition, you may have to set
`LD_LIBRARY_PATH` so that the linker can find the directory containing
the library.

For example, it is possible to compile the examples with purely local
header and library files.  The library files `libsunlightlb.so*`
should first be built in the `src` directory, as in the first three
steps above, then :

* Go into the `examples` directory.
* Edit the top part of the `Makefile`, in particular add `-I../src -L../src` to the `CC_FLAGS` variable.
* Type `make sphere` (it should work now).
* Typing `LD_LIBRARY_PATH=../src/ ./sphere` should run the code. 

If you don't want to work with a shared object library, the source
code can be copied into the same directory as the driver code and all
compiled together.

#### To build and install SunlightLB as a perl extension module.

Firstly, build and install the shared object library as before.

* Go into the `src` directory.
* Edit the perl-related parts of the `Makefile.`
* Type `make 4perl` to use SWIG to make the perl module files.
* Type `make install4perl` to install the perl module files.
* Typing `perl sphere.pl` in the `examples` directory runs the code. 

#### To build and install SunlightLB as a python module.

This is very similar to building a perl module.

* Go into the `src` directory.
* Edit the python-related parts of the `Makefile`.
* Type `make 4python` to use SWIG to make the python module files.
* Type `make install4python` to install the python module files.
* Typing `python2.2 sphere.py` in the `examples` directory runs the code. 

#### Monitor utility

The Tcl/Tk script `lbmonitor` has been included in the `src`
directory.  This provides a GUI to monitor and check-point lattice
Boltzmann simulations that generate `*.mon` files.  To run this, you
can copy it to the directory in which the LB executables / scripts
will be run, and run it with `wish lbmonitor`. Alternatively, find out
where the Tcl/Tk `wish` command is and edit the first line of
`lbmonitor`, then put it somewhere that will be found on the path.

### Copying

SunlightLB is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SunlightLB is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<http://www.gnu.org/licenses/>.

### Copyright

SunlightLB is based on an original code copyright &copy; 2005-2017
Unilever UK Central Resources Ltd (Registered in London number 29140;
Registered Office: Unilever House, 100 Victoria Embankment, London
EC4Y 0DY, UK).

### Contact

Send email to patrick{dot}warren{at}unilever{dot}com

Send paper mail to Dr Patrick B Warren, Unilever R&D Port Sunlight,
Quarry Road East, Bebington, Wirral, CH63 3JW, UK.
