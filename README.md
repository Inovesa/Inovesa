Inovesa
=======

Inovesa (Inovesa Numerical Optimized Vlasov Equation Solver Application) is
a tool developed to simulate the dynamics of an electron bunch in a storage
ring, including the self-interaction with its own wake field.
To do so, it  uses the well established method to numerically solve the
Vlasov-Fokker-Planck equation.

Inovesa is modularly extensible and uses OpenCL to massively parallelize the
computation. It was designed with standard desktop PCs and usability in mind.
The working principle and example numerical studies can be found in the
[paper describing Inovesa (PRAB)][1].

Inovesa and the experience of any user greatly profits from
contributions of any form -- starting from questions or bug reports
going all the way to implementation of new features.
If you consider contributing to Inovesa,
please read the [contibuting file][2].


[1]: https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.20.030704 "Parallelized Vlasov-Fokker-Planck solver for desktop personal computers"

[2]: (CONTRIBUTING.md)
"Information on contributing to Inovesa"


Build Requirenments
-------------------

CMake will check whether build requirenmens are met.
On Debian 8 the following packages have to be installed.

Mandatory:
* cmake
* g++
* libboost-dev
* libboost-system-dev
* libboost-filesystem-dev
* libboost-program-options-dev
* libfftw3-dev

Optional:
* libglew-dev (for GUI support)
* libglfw3-dev (for GUI support)
* libhdf5-dev (for HDF5 support)
* libpng++-dev (for PNG support)
* libxrandr-dev (for GUI support)
* opencl-dev (optional for parallelization)
* libclfft-dev (for faster FFT using parallelization)

