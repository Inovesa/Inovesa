INOVESA
=======

INOVESA Numerical OpenCL Vlasov Equation Solver Application


Build Requirenments
-------------------

CMake will check whether build requirenmens are met.
On Debian 8 the following packages have to be installed:

* libglew-dev
* libglfw3-dev (own cmake script)
* opencl-dev (if in doubt: amd-opencl-dev runs on all x86 CPU, own cmake script)
* libxrandr-dev
* libglm-dev (own cmake script)
* libgsl0-dev (for comparission with legacy code, own cmake script)
* libfftw3-dev (currently needed additionally to GSL, own cmake script)
* libhdf5-dev
* libpng++-dev
* libboost-dev
* libboost-system-dev
* libboost-filesystem-dev
* libboost-program-options-dev

