Inovesa
=======

Inovesa (Inovesa Numerical Optimized Vlasov Equation Solver Application) is
a tool developed to simulate the dynamics of an electron bunch in a storage
ring, including the self-interaction with its own wake field.
To do so, it  uses the well established method to numerically solve the
Vlasov-Fokker-Planck equation.

Inovesa is modularly extensible and uses OpenCL to massively parallelize the
computation. It was designed with standard desktop PCs and usability in mind.
Inovesa and the experience of any user greatly profits from contributions of
any form -- starting from questions or bug reports going all the way to
implementation of new features. We warmly welcome input to Inovesa.
If you consider contributing, the [contibuting file](CONTRIBUTING.md)
is a good starting point.

[![DOI](https://zenodo.org/badge/73905339.svg)](https://zenodo.org/badge/latestdoi/73905339)
[![Inovesa build status (branch: Master)](https://travis-ci.org/Inovesa/Inovesa.svg)](https://travis-ci.org/Inovesa/Inovesa/branches)
[![Coverage Status](https://coveralls.io/repos/github/Inovesa/Inovesa/badge.svg)](https://coveralls.io/github/Inovesa/Inovesa)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/Inovesa/Inovesa.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Inovesa/Inovesa/context:cpp)


Installation
------------

Compiling inovesa follows the standard build procedure of CMake.
In general, CMake will check whether build requirenmens are met.
Details can be found at the corresponding [Inovesa Wiki page](https://github.com/Inovesa/Inovesa/wiki/Installation).

### Mandatory Dependencies
* [cmake](https://cmake.org/) (Version 3.1 or later)
* C++ 14 compatible compiler (e.g. [g++](https://gcc.gnu.org/) 5.0 or later, [clang](http://clang.llvm.org/) 3.5 or later)
* [boost](http://www.boost.org/) (boost-system, boost-filesystem, boost-program-options)
* [FFTW](http://fftw.org/) (Version 3)

### Optional Dependencies
* [GLEW](https://www.opengl.org/sdk/libs/GLEW/) (for GUI support)
* [GLFW](http://www.glfw.org/) (for GUI support)
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/) (to write out data)
* [PNG++](http://www.nongnu.org/pngpp/) (for PNG support)
* [OpenCL](https://www.khronos.org/opencl/) (>=1.1, for parallelization)
* [clFFT](https://github.com/clMathLibraries/clFFT) (for faster FFT using parallelization)



Publications
------------

### Background on the implementation

* The fundamental working principle (v1.0) and example numerical studies can be found in the [first paper on Inovesa][2].
* The features added in v1.1 are described in a [paper on modeling synchrotron motion][3]


### Example results

* Analysis of [bunch profiles using machine learning][4] (Inovesa v1.0)
* Studies of [phase shift][5] (Inovesa v1.0)
* A wide [range of applications][6] (Inovesa v1.0)
* [Emitted synchrotron radiation][7] of different wavelengths (Inovesa v1.0)


### Citing Inovesa

We use the zenodo project to get a DOI for each version. So, when you use
Inovesa to obtain your simulation results, you can
[search zenodo for the right citation of your Inovesa version](https://zenodo.org/search?page=1&size=20&q=conceptrecid:597356&all_versions&sort=-version)
to directly point to that specific version.



[1]: CONTRIBUTING.md
"Information on contributing to Inovesa"

[2]: https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.20.030704 "Parallelized Vlasov-Fokker-Planck solver for desktop personal computers"

[3]: http://iopscience.iop.org/article/10.1088/1742-6596/1067/6/062025/meta "Elaborated Modeling of Synchrotron Motion in Vlasov-Fokker-Planck Solvers"

[4]: https://doi.org/10.18429/JACoW-IPAC2018-THPAK030 "Studies of Longitudinal Dynamics in the Micro-Bunching Instability Using Machine Learning"

[5]: https://doi.org/10.18429/JACoW-IPAC2018-WEPAL028 "Study of the Influence of the CSR Impedance on the Synchronous Phase Shift at KARA"

[6]: https://doi.org/10.5445/ir/1000084466 "Simulation and measurement of the dynamics of ultra-short electron bunch profiles for the generation of coherent THz radiation"

[7]: https://arxiv.org/abs/1710.09568 "Continuous bunch-by-bunch spectroscopic investigation of the micro-bunching instability"
