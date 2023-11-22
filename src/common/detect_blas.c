//
// Created by andrew on 11/21/2023.
//
#include "detect_blas.h"
#include <dlfcn.h>

openblas_init_t get_openblas_parallel(void* libloc) {
    const openblas_init_t parallel_address = (openblas_init_t)dlsym(libloc, "openblas_get_parallel");
    return parallel_address;
}

openblas_set_t get_openblas_set(void* libloc) {
    const openblas_set_t set_address = (openblas_set_t)dlsym(libloc, "openblas_set_num_threads");
    return set_address;
}