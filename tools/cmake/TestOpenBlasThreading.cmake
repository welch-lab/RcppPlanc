function(OPENBLAS_USES_OPENMP OUT_ARG HEADER_DEFINED)
    if(NOT HEADER_DEFINED MATCHES "^%")
        string(PREPEND HEADER_DEFINED "-D")
    endif()
    try_run(runSuccess compileSuccess ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_openblas_threading.c
            COMPILE_DEFINITIONS ${HEADER_DEFINED}
            LINK_LIBRARIES openblas RUN_OUTPUT_VARIABLE runOut COMPILE_OUTPUT_VARIABLE compileOut)
    message("OpenBLAS parallelism of type ${runOut}")
    set(${OUT_ARG} ${runOut} PARENT_SCOPE)
endfunction()