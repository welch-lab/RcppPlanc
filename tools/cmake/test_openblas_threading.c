#include "../../src/common/config.h"
#include <stdio.h>
int main(int argc, char *argv[]){
    int outnum = openblas_get_parallel();
    printf("%d\n", outnum);
    return 0;
}