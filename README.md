# RcppPlanc

## About

RcppPlanc is, as the name implies, an Rcpp wrapper for [planc](https://github.com/ramkikannan/planc). Currently, it
implements wrappers for the vanilla NMF and NNLS algorithms. It also implements integrative NMF as described in
[Welch, J. D. et al, 2019](https://doi.org/10.1016/j.cell.2019.05.006), online integrative NMF as described in
[Gao, C. et al, 2021](https://doi.org/10.1038/s41587-021-00867-x) and unshared integrative NMF as described in
[Kriebel and Welch, 2022](https://doi.org/10.1038/s41467-022-28431-4). This is, at present, the fastest NMF implementation on CRAN, assuming you set nCores on call.

## Requirements

- A fast BLAS library that provides a cblas header ([OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) works well. Mac
OS on ARM64 will always use the built in [vecLib](https://developer.apple.com/documentation/accelerate/veclib?language=objc) framework.)
- A modern R toolchain (Rtools on Windows, the system toolchain on UNIX-alikes)
- A modern R (tested on 4.2+)
- CMake >=3.20
- HDF5 (tested against 1.12.0)
- OpenMP implementation that uses v4.0 semantics (if using OpenMP). For the GNU stack, this means gcc 9.

## Installation

### From Source:
1. Install the requirements above.
2. Ensure your libraries can be found by CMake.
3. `devtools::install_github(welch-lab/RcppPlanc)`
### Binary (R-universe, CRAN pending)
#### MacOS and Windows
`install.packages('RcppPlanc', repos = c('https://welch-lab.r-universe.dev', 'https://cloud.r-project.org'))`
#### Ubuntu 22.04 (YMMV on other linuxes)
1. Ensure libhdf5, libhwloc, and libopenblas are installed.

```
install.packages("RcppPlanc", repos = c(
linux = 'https://welch-lab.r-universe.dev/bin/linux/jammy/4.3',
sources = 'https://welch-lab.r-universe.dev',
cran = 'https://cloud.r-project.org'))
```

## Caveats

Non-fatal errors during the CMake run are to be expected as the compiler flags are tailored to each system.
Currently, the cblas locator logic only finds headers at the root of the given include directories or in the openblas
and flexiblas subdirectories. If CMake cannot find a cblas header, its directory must be specified manually in CMakeCache.txt.
This software will not use OpenMP on MacOS and performance on Intel Macs will suffer accordingly. Blame Apple. ARM64 Macs
more than make up for the threaded performance loss with their inbuilt matrix coprocessor. If you're on an Intel mac and feeling
adventurous, you can always try following the instructions at <https://mac.r-project.org/openmp/> but your mileage may vary.

## Citations and Prior Works

A huge shoutout goes to [Ramakrishnan Kannan](https://github.com/ramkikannan) for the original work on which this program is based.
Relevant citations can be found in his repository at <https://github.com/ramkikannan/planc/blob/master/papers.md>.

When and if a citation for this package becomes available, it will be posted here.

Parts of this code source are licensed under the 3-clause BSD (c) 2016 , UT-Battelle, LLC and (c) 2017 Ramakrishnan Kannan.
The compiled binary and sources not marked otherwise are licensed under the GPLv2 or later.

Licenses for the various CMake scripts (FindHWLOC, FindR, and FindRModule, and PatchFile) are included at the top
of each file.

## Contributing and Support

Please reach out to us at the Welch Lab using the issue tracker in this repository *before* contacting upstream.
Much of the inherited code is heavily modified and half of the algorithms have not yet been upstreamed.

Pull requests are welcome.

## Documentation

Vignettes and R documentation files are available in the standard folders.
