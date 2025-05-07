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
