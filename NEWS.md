# RcppPlanc 2.0.12
* Added a few symbol removals to hwloc patches
* Added an error message to address the real issue in [17](https://github.com/welch-lab/RcppPlanc/issues/17/)
* Silenced confusing CMake find_package calls
* Removed unnecessary CMakePresets from source.

# RcppPlanc 2.0.11
* Re-enabled static library builds for CRAN on *nix to silence (questionable) ODR violation from libASAN.
* Replaced Rcpp::Nullable\<Rcpp::NumericMatrix\> with Rcpp::Nullable\<arma::mat\>, addressing  issue [17](https://github.com/welch-lab/RcppPlanc/issues/17/).
* Re-enabled roxygen in-tree build.
* Added mechanism for finding OPENMP_SHLIB_CXXFLAGS in static RCPP case

# RcppPlanc 2.0.10
* Removed masked cblas_sgemm as it was never being called and causing all sorts of linkage issues
* Remove cblas library linkage and some includes.

# RcppPlanc 2.0.9
* Add check for path specifically used by Prof. Ripley's CRAN builders to disable MKL
* Fix generic linkage for BLAS/LAPACK in BUILD_RCPP case

# RcppPlanc 2.0.8
* Add runtime check to make sure MKL actually works if detected

# RcppPlanc 2.0.7
* Allow for empty LIBR_STRING on certain Linux machines.

# RcppPlanc 2.0.6
* Sanitized regex from Rscript call in src/planc/CMakeLists/FindR.cmake

# RcppPlanc 2.0.5
* Fixed typo in FindR.cmake (oops!)
* Use HWLOC_LIBRARIES instead of HWLOC_LDFLAGS
* Fixed mistake in H5SpMat constructor example
* Added copyright indicators for HighFive and Armadillo
* Fixed MKL detection
* Fix error in UINMF logic

# RcppPlanc 2.0.4
* Fixed data.cpp opening H5SpMat as ReadWrite instead of ReadOnly
* CMake now checks R_RHOME before path.
* Added fallback to guess R_LDFLAGS when Makeconf parsing fails.
* Removed bad check_symbol_exists test in FindHWLOC.
* Use HWLOC_LDFLAGS instead of pkgcfg_lib_HWLOC_hwloc to capture both pkgconf and find_library situations.

# RcppPlanc 2.0.3
* 3rd attempt to refactor the BLAS detection of the CMake system. (0.1 and 0.2 were failed attempts at same)

# RcppPlanc 2.0.0

* Initial CRAN release
* Fixed critical algorithmic bugs in INMF/UINMF with intial seeding (API Breaking)
* Used libplanc/nmflib 1.0.0 for backend

# RcppPlanc 1.0.0-rc3

* Non-CRAN release
* Created wrapper functions for NMF with ANLS-BPP, ADMM, HALS and MU algorithms supported
* Created wrapper functions for symNMF with ANLS-BPP and GNSYM algortihms supported
* Implemented iNMF, onlineINMF and UINMF with C++
