Addressed all NOTES that I believe fixable.
_______________________________________________________________________________________________________________
# Version 1.0.0

## Test Environments

* Local:
  * Windows 11, R 4.3.2 (x86_64-w64-mingw32)
  * r-hub/containers/ubuntu-release (R 4.3.2)
  * r-hub/containers/ubuntu-next (R 4.3.1 Patched r85316)
  * r-hub/containers/ubuntu-gcc12 (R-Devel r85316) -- additional packages: libdeflate-dev
  * r-hub/containers/ubuntu-clang (R-Devel r85316) -- additional packages: libdeflate-dev, libstdc++-12-dev
  * r-hub/containers/atlas (Fedora 38, R-Devel r85316)
  * r-hub/containers/gcc13 (Fedora 38, R-Devel r85316)
  * r-hub/containers/intel (Fedora 38, R-Devel r85316)
  * r-hub/containers/mkl (Fedora 38, R-Devel r85316)
  * r-hub/containers/nold  (Ubuntu 22.04, R-Devel r85316) -- additional packages: libdeflate-dev
  * r-hub/containers/clang-asan (Ubuntu 22.04, R-Devel r85316) -- additional packages: libdeflate-dev
* win-builder:
  * Windows Server 2022, Devel
* R-hub Builder (<https://builder.r-hub.io>)
  * Windows Server 2022 R-Devel
  * Ubuntu Linux 20.04.1 LTS, R-devel with rchk


## R CMD check results

``` <!-- language: lang-none -->

/root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:88:10: warning: explicit comparison with infinity in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:98:10: warning: explicit comparison with infinity in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:134:10: warning: explicit comparison with NaN in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:144:10: warning: explicit comparison with NaN in fast floating point mode [-Wtautological-constant-compare]
See ‘/check/RcppPlanc.Rcheck/00install.out’ for details.
```

This warning is exclusive to the intel compiler and is an issue with RcppArmadillo.

``` <!-- language: lang-none -->
* checking R files for syntax errors ... WARNING
Warning in Sys.setlocale("LC_CTYPE", "en_US.UTF-8") :
  OS reports request to set locale to "en_US.UTF-8" cannot be honored
```

The R-Hub ubuntu containers appear to not have locales installed. At least, some of them.


  ``` <!-- language: lang-none -->

❯ checking installed package size ... NOTE
    installed size is 33.0Mb
    sub-directories of 1Mb or more:
      libs  30.9Mb
```

Static linkage of quite a few libraries will do this. Reported size is from my Windows machine. Unix-alike installs
use dynamic linkage and embedded hwloc and will be smaller. Linux builds are typically ~5MB even with statically linked
hwloc.

``` <!-- language: lang-none -->
❯ checking line endings in C/C++/Fortran sources/headers ... NOTE
 Found the following sources/headers with CR or CRLF line endings:
  src/CMakeFiles/3.24.3/CompilerIdC/CMakeCCompilerId.c
  src/CMakeFiles/3.24.3/CompilerIdCXX/CMakeCXXCompilerId.cpp
  src/CMakeFiles/CheckTypeSize/SIZE_OF_VOIDP.c
  src/CMakeFiles/FindOpenMP/OpenMPCheckVersion.c
  src/CMakeFiles/FindOpenMP/OpenMPCheckVersion.cpp
  src/CMakeFiles/FindOpenMP/OpenMPTryFlag.c
  src/CMakeFiles/FindOpenMP/OpenMPTryFlag.cpp
  src/_deps/highfive-src/include/highfive/H5Version.hpp
Some Unix compilers require LF line endings.
```

  These are not included in the tarball, they are downloaded at configure time by git. On UNIX-alikes this note is absent.

``` <!-- language: lang-none -->
Found the following Makefile(s) with CR or CRLF line endings:
  src/Makefile
  src/_deps/highfive-build/Makefile
  src/_deps/highfive-subbuild/Makefile
  src/_deps/hwloc-subbuild/Makefile
Some Unix 'make' programs require LF line endings.

```
These are not included in the tarball, they are generated at configure time. On UNIX-alikes this note is absent.


``` <!-- language: lang-none -->
❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-freciprocal-math' '-funsafe-math-optimizations'
```

These compilation flags are checked at configure time and will not be used if they are not available.


  ``` <!-- language: lang-none -->
❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

```

Per <https://github.com/r-hub/rhub/issues/560> this appears to be a bug.

  ``` <!-- language: lang-none -->
❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException''

```

Per <https://github.com/r-hub/rhub/issues/503> this appears to be a bug.

  ``` <!-- language: lang-none -->
❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable

```

This is an Rhub issue. It works locally.

``` <!-- language: lang-none -->
❯ checking for future file timestamps ... NOTE
  unable to verify current time

```

Per <https://github.com/r-hub/rhub/issues/440> this appears to be a bug.


```c++
  Direct leak of 16 byte(s) in 1 object(s) allocated from:
#0 0x555ffa2c6318 in __interceptor_calloc (/opt/R/devel-asan/lib/R/bin/exec/R+0xb8318) (BuildId: f9ace1cd1e557acb45d82ef6d66e646229ab6542)
#1 0x7f81f807ea0f  (/usr/lib/x86_64-linux-gnu/libc++abi.so.1+0x25a0f) (BuildId: 47782717d7fe8d90800ba4d11684a17cfa249698)
#2 0x7f815f03f2a4 in planc::INMF<arma::Mat<double>>::INMF(std::__1::vector<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>, std::__1::allocator<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>>>&, unsigned int, double, bool) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/common/inmf.hpp:125:19
#3 0x7f815f035237 in planc::BPPINMF<arma::Mat<double>>::BPPINMF(std::__1::vector<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>, std::__1::allocator<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>>>&, unsigned int, double) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/nmf/bppinmf.hpp:137:82
#4 0x7f815f035237 in Rcpp::Vector<19, Rcpp::PreserveStorage> runINMF<arma::Mat<double>>(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>, unsigned int, double, unsigned int, bool) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:572:23
#5 0x7f815ed901f5 in bppinmf_dense(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>> const&, unsigned int, double, unsigned int, bool, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<arma::Mat<double>>) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:632:16
#6 0x7f815eda1d36 in bppinmf(Rcpp::Vector<19, Rcpp::PreserveStorage>, unsigned int, double, unsigned int, bool, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<arma::Mat<double>>) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:666:16
#7 0x7f815f1cc172 in _RcppPlanc_bppinmf /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/RcppExports.cpp:134:34
#8 0x7f81f8283262 in R_doDotCall /tmp/R-devel/src/main/dotcode.c
```
This appears to be a libc++ issue with exception throwing per https://stackoverflow.com/questions/61015622/throw-with-clang-shows-possible-memory-leak.

```<!-- language: lang-none -->
Function Rcpp::Vector<19, Rcpp::PreserveStorage> onlineINMF_S1_mem<arma::Mat<double> >(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double> > >, unsigned int, double, unsigned int, unsigned int, unsigned int, bool)
  [UP] ignoring variable <unnamed var:   %36 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %40 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %44 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %49 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %53 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %57 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  and similar rchk errors
```
These appear to be issues with rchk's handling of Armadillo's template metaprogramming.


0 errors ✔ | 0 warnings ✖ | 4 notes ✖

## revdepcheck results

None as this is a new package.
