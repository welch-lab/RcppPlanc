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
