//
// Created by andrew on 11/21/2023.
//
#include "detect_blas.h"
#include <dlfcn.h>
#ifdef _OPENMP
#include <omp.h>
#endif

openblas_init_t get_openblas_parallel(void* libloc) {
    const openblas_init_t parallel_address = (openblas_init_t)dlsym(libloc, "openblas_get_parallel");
    return parallel_address;
}

openblas_set_t get_openblas_set(void* libloc) {
    const openblas_set_t set_address = (openblas_set_t)dlsym(libloc, "openblas_set_num_threads");
    return set_address;
}

bool is_openmp(void) {
#ifdef _OPENMP
    const int threads = omp_get_num_threads();
    if (threads == 1) return false;
    else return true;
#endif
    return false;
}