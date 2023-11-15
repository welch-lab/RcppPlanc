function(OPENBLAS_USES_OPENMP OUT_ARG)
    try_run(runSuccess compileSuccess ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_openblas_threading.c
            COMPILE_DEFINITIONS "-DHAVE_OPENBLAS_CBLAS_H"
            LINK_LIBRARIES openblas RUN_OUTPUT_VARIABLE runOut COMPILE_OUTPUT_VARIABLE compileOut)
    set(${OUT_ARG} ${runOut} PARENT_SCOPE)

endfunction()