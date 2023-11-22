//
// Created by andrew on 11/21/2023.
//

#ifndef DETECT_BLAS_H
#define DETECT_BLAS_H
#include <stdbool.h>

typedef int(*openblas_init_t)(void);
openblas_init_t get_openblas_parallel(void* libloc);

typedef void(*openblas_set_t)(int);
openblas_set_t get_openblas_set(void* libloc);

bool is_openmp(void);

#endif //DETECT_BLAS_H
