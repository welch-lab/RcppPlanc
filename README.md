# RcppPlanc

## About

RcppPlanc is, as the name implies, an Rcpp wrapper for [planc](https://github.com/ramkikannan/planc). Currently, it
implements wrappers for the vanilla NMF and NNLS algorithms. It also implements integrative NMF as described in
[Yang and Michailidis, 2016](https://doi.org/10.1093/bioinformatics/btv544), online integrative NMF as described in
[Welch, J. D. et al, 2019](https://doi.org/10.1038/s41587-021-00867-x) and unshared integrative NMF as described in
[Kriebel and Welch, 2022](https://doi.org/10.1038/s41467-022-28431-4).

## Requirements

- A fast BLAS library that provides a cblas header ([OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) works well. Mac
OS will always use the built in veclib framework.)
- A modern R toolchain (Rtools on Windows, the system toolchain on UNIX-alikes)
- A modern R (tested on 4.2+)
- CMake >=3.20
- HDF5 (tested against 1.12.0)

## Caveats

Non-fatal errors during the CMake run are to be expected as the compiler flags are tailored to each system.
Currently, the cblas locator logic only finds headers at the root of the given include directories or in the openblas
and flexiblas subdirectories. If CMake cannot find a cblas header, its directory must be specified manually in CMakeCache.txt.
This software will not use OpenMP on MacOS and performance on Intel Macs will suffer accordingly. Blame Apple. ARM64 Macs
more than make up for the threaded performance loss with their inbuilt matrix coprocessor.

## Citations and Prior Works

A huge shoutout goes to [Ramakrishnan Kannan](https://github.com/ramkikannan) for the original work on which this program is based.
Relevant citations can be found in his repository at <https://github.com/ramkikannan/planc/blob/master/papers.md>.

When and if a citation for this package becomes available, it will be posted here.

## Contributing and Support

Please reach out to us at the Welch Lab using the issue tracker in this repository *before* contacting upstream.
Much of the inherited code is heavily modified and half of the algorithms have not yet been upstreamed.

Pull requests are welcome.

## Documentation

Vignettes and R documentation files are available in the standard folders.
