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

[2]: CONTRIBUTING.md
"Information on contributing to Inovesa"

Installation
------------

Compiling inovesa follows the standard build procedure of CMake.
Details can be found at the corresponding [Inovesa Wiki page](https://github.com/Inovesa/Inovesa/wiki/Installation).

Build Requirenments
-------------------

CMake will check whether build requirenmens are met.
The following is needed:

### Mandatory Build Requirenments
* [cmake](https://cmake.org/) (Version 3.1 or later)
* C++ 14 compatible compiler (e.g. [g++](https://gcc.gnu.org/) 5.0 or later, [clang](http://clang.llvm.org/) 3.5 or later)
* [boost](http://www.boost.org/) (boost-system, boost-filesystem, boost-program-options)
* [FFTW](http://fftw.org/) (Version 3)

### Optional Build Requirenments
* [GLEW](https://www.opengl.org/sdk/libs/GLEW/) (for GUI support)
* [GLFW](http://www.glfw.org/) (for GUI support)
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) (to write out data)
* [PNG++](http://www.nongnu.org/pngpp/) (for PNG support)
* [OpenCL](https://www.khronos.org/opencl/) (>=1.1, optional for parallelization)
* [clFFT](https://github.com/clMathLibraries/clFFT) (for faster FFT using parallelization)

