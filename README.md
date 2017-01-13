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
[paper describing Inovesa][1].


[1]: https://arxiv.org/abs/1611.05293 "A Parallelized Vlasov-Fokker-Planck-Solver for Desktop PCs"


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

