cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
set(CMAKE_OSX_DEPLOYMENT_TARGET "11" CACHE STRING "Target MacOS Version")
mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
if(WIN32)
SET(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
endif()
project(RcppPlanc VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CheckIncludeFile)
include(CheckCXXCompilerFlag)
include(FetchContent)
include(ExternalProject)
include(CheckSymbolExists)
include(CheckCXXSymbolExists)
include(cmake/FindRModule.cmake)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

if(APPLE)
    set(CMAKE_PROGRAM_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/bin)
    set(CMAKE_LIBRARY_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/lib)
    set(CMAKE_INCLUDE_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/include)
endif()

separate_arguments(PASSTHROUGH_CONFIGURE_ARGS)

find_package(R REQUIRED)

add_library(RcppPlanc SHARED "${PROJECT_SOURCE_DIR}/../src/rcppplanc_nmf.cpp" "${PROJECT_SOURCE_DIR}/../src/RcppExports.cpp" "${PROJECT_SOURCE_DIR}/../src/conversion.cpp" "${PROJECT_SOURCE_DIR}/../src/common/detect_blas.c" $<$<BOOL:${WIN32}>:${PROJECT_SOURCE_DIR}/../src/RcppPlanc-win.def>)
if(WIN32)
target_link_options(RcppPlanc PRIVATE -static-libgcc $<$<CONFIG:RELEASE>:-s> $<$<CONFIG:RELWITHDEBINFO>:-s>)
else()
    target_link_options(RcppPlanc PRIVATE "-shared")
endif()
set_target_properties(RcppPlanc PROPERTIES PREFIX "")
if(APPLE AND R_SHLIB_EXT MATCHES ".so")
set_target_properties(RcppPlanc PROPERTIES SUFFIX ${R_SHLIB_EXT})
endif()
target_include_directories(RcppPlanc
PRIVATE "${PROJECT_SOURCE_DIR}/../src/common;${PROJECT_SOURCE_DIR}/../src/nnls;${PROJECT_SOURCE_DIR}/../src/nmf")
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
execute_process(COMMAND clang -v ERROR_VARIABLE CLANGID)
string(REGEX MATCH "\\(.*\\)" CLANGPARSE ${CLANGID})
if(CLANGPARSE MATCHES "\\(CLANG\\: AOCC.*")
set(BLA_VENDOR AOCL_mt)
target_compile_options(RcppPlanc PRIVATE -zopt -O3 -flto -fstrip-mining -freciprocal-math -finline-recursion=3
                                     -fremap-arrays -fnt-store=aggressive -mrecip=all -optimize-strided-mem-cost
                                     -fenable-aggressive-gather -favoid-fpe-causing-opt)
endif()
endif()

set(MKL_DIR "${PROJECT_SOURCE_DIR}/cmake")
find_package(MKL QUIET)
if(MKL_FOUND AND NOT DEFINED BLA_VENDOR)
    message(STATUS "found MKL")
    set(USE_MKL)
    set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
    target_compile_definitions(RcppPlanc PRIVATE "MKL_FOUND")
    set(CBLAS_H_DIR "${MKL_INCLUDE_DIR}")
elseif(WIN32 AND NOT DEFINED ENV{CONDA_BUILD})
    find_package(OpenBLAS CONFIG REQUIRED)
    set(BLAS_LIBRARIES "${OpenBLAS_LIBRARIES}")
    set(LAPACK_LIBRARIES "")
    set(CBLAS_H_DIR "${OpenBLAS_INCLUDE_DIR}")
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

if(NOT(APPLE AND(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES arm64)))
    find_package(OpenMP)
endif()

find_r_module(Rcpp)
find_r_module(RcppArmadillo)
find_r_module(RcppProgress)

if(R_INCLUDE_DIR)
target_include_directories(RcppPlanc PRIVATE "${R_INCLUDE_DIR}")
else()
message(FATAL_ERROR "missing R headers")
endif()
if(R_RCPP)
target_include_directories(RcppPlanc PRIVATE "${R_RCPP}/include")
else()
message(FATAL_ERROR "missing RcppArmadillo")
endif()
if(R_RCPPARMADILLO)
target_include_directories(RcppPlanc PRIVATE "${R_RCPPARMADILLO}/include")
else()
message(FATAL_ERROR "missing RcppArmadillo")
endif()
if(R_RCPPPROGRESS)
target_include_directories(RcppPlanc PRIVATE "${R_RCPPPROGRESS}/include")
else()
message(FATAL_ERROR "missing RcppProgress")
endif()

# from scivision/patchfile, MIT license c. Michael

if(WIN32)
    # prioritize Git Patch on Windows as other Patches may be very old and incompatible.
    find_package(Git)

    if(Git_FOUND)
        get_filename_component(GIT_DIR ${GIT_EXECUTABLE} DIRECTORY)
        get_filename_component(GIT_DIR ${GIT_DIR} DIRECTORY)
    endif()
endif()

find_program(PATCH
    NAMES patch
    HINTS ${GIT_DIR}
    PATH_SUFFIXES usr/bin
)

add_library(hw_detect OBJECT "${PROJECT_SOURCE_DIR}/../src/common/hw_detect.c")
    find_package(HWLOC)
    if(HWLOC_FOUND)
    target_link_libraries(hw_detect PRIVATE "${HWLOC_LDFLAGS}")
    target_include_directories(hw_detect PRIVATE "${HWLOC_INCLUDE_DIRS}")
    else()
    file(GLOB hwloc_patch_list LIST_DIRECTORIES FALSE "${CMAKE_CURRENT_SOURCE_DIR}/patches/*.patch")
    set(RTOOLS_LD "ld")
    set(cmd "sh")
    list(APPEND cmd ../hwloc/configure --with-pic --prefix=${hwloc_SOURCE_DIR}/out/ --disable-picky --enable-shared=no --enable-static=yes)
    list(APPEND cmd --disable-cairo --disable-libxml2 --disable-io --disable-plugin-dlopen --disable-plugin-ltdl --enable-embedded-mode)
    if(WIN32)
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "AMD64")
        list(APPEND cmd --build=x86_64-w64-mingw32)
    elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "ARM64")
        list(APPEND cmd --build=aarch64-w64-mingw32)
        set(RTOOLS_LD "ld.lld")
    endif()
    elseif(PASSTHROUGH_CONFIGURE_ARGS)
        list(APPEND cmd ${PASSTHROUGH_CONFIGURE_ARGS})
    endif()
    ExternalProject_Add(hwloc
        URL https://download.open-mpi.org/release/hwloc/v2.9/hwloc-2.9.3.tar.gz
        URL_HASH SHA256=5985db3a30bbe51234c2cd26ebe4ae9b4c3352ab788b1a464c40c0483bf4de59
        PATCH_COMMAND ${CMAKE_COMMAND} -E cat ${hwloc_patch_list} | ${PATCH} -p1 -l
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env LD=${RTOOLS_LD} "CFLAGS=-DNDEBUG $ENV{CFLAGS} $<$<BOOL:${APPLE}>:-mmacosx-version-min=11>" ${cmd}
    )
        ExternalProject_Get_Property(hwloc SOURCE_DIR BINARY_DIR)
    target_include_directories(hw_detect PRIVATE "${BINARY_DIR}/include" "${SOURCE_DIR}/include")
        add_custom_command(TARGET hwloc POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E remove_directory ${SOURCE_DIR}/utils/lstopo
                VERBATIM)
    add_dependencies(hw_detect hwloc)
    target_link_libraries(hw_detect PRIVATE "${BINARY_DIR}/hwloc/.libs/libhwloc_embedded.a")
endif()

if(DEFINED ENV{CONDA_BUILD})
    find_package(HDF5 REQUIRED)
    set(CMAKE_REQUIRED_INCLUDES highfive)
    check_include_file_cxx("highfive/highfive.hpp" HighFive_Header)
    if(HighFive_Header)
        target_include_directories(RcppPlanc PRIVATE highfive)
    endif()
else()
    if(APPLE)
        target_link_libraries(hw_detect INTERFACE "-framework Foundation" "-framework IOKit" "-framework OpenCL")
    endif()
set(HIGHFIVE_EXAMPLES OFF CACHE INTERNAL BOOL)
set(HIGHFIVE_BUILD_DOCS OFF CACHE INTERNAL BOOL)
set(HIGHFIVE_USE_BOOST OFF CACHE INTERNAL BOOL)
set(HIGHFIVE_UNIT_TESTS OFF CACHE INTERNAL BOOL)
set(HIGHFIVE_HAS_CONCEPTS OFF CACHE INTERNAL BOOL)
if((WIN32 OR APPLE))
    set(HDF5_USE_STATIC_LIBRARIES ON CACHE INTERNAL BOOL)
    set(HIGHFIVE_STATIC_HDF5 ON CACHE INTERNAL BOOL)
endif()
FetchContent_Declare(
    HighFive
    URL https://github.com/BlueBrain/HighFive/archive/refs/tags/v2.9.0.tar.gz
    URL_HASH SHA256=6301def8ceb9f4d7a595988612db288b448a3c0546f6c83417dab38c64994d7e
)
FetchContent_MakeAvailable(HighFive)
# from https://gitlab.kitware.com/cmake/cmake/-/issues/18872#note_1090297
# Since there is no compiler line to help us, add the additional required
# libraries manually.
set(_additional_libs sz z dl m)
foreach(_additional_lib IN LISTS _additional_libs)
    # If both static and shared are available, prefer static to avoid libdl
    # annoyances ("Using 'dlopen' in statically linked applications requires
    # at runtime the shared libraries from the glibc version used for
    # linking")
    if(HDF5_USE_STATIC_LIBRARIES)
        set(_libnames ${_additional_lib} lib${_additional_lib}.a)
    else()
        set(_libnames ${_additional_lib})
    endif(HDF5_USE_STATIC_LIBRARIES)
    set(_libvar "LIB_${_additional_lib}")
    find_library(${_libvar}
            NAMES ${_libnames}
            HINTS ${HDF5_ROOT}
            PATH_SUFFIXES lib Lib lib64 Lib64
            REQUIRED)
    if(HDF5_FIND_DEBUG)
        message(STATUS "${_additional_lib} (${_libvar}): ${${_libvar}}")
    endif(HDF5_FIND_DEBUG)
    if(${${_libvar}} STREQUAL "${_libvar}-NOTFOUND")
        message(FATAL_ERROR "${_additional_lib} not found, required for HDF5 install")
    endif()
    list(APPEND HDF5_LIBRARIES ${${_libvar}})
endforeach()
list(REMOVE_DUPLICATES HDF5_LIBRARIES)
endif()
string(REGEX MATCH "-lR" R_LINK_LINE ${R_LDFLAGS})
string(REGEX MATCH "-framework R" R_FRAMEWORK_LINE ${R_LDFLAGS})
string(REPLACE " -lR" "" R_LDFLAGS ${R_LDFLAGS})
string(REPLACE " -framework R" "" R_LDFLAGS ${R_LDFLAGS})
target_link_options(RcppPlanc PRIVATE ${R_LDFLAGS})
target_link_libraries(RcppPlanc PRIVATE ${R_LINK_LINE} ${R_FRAMEWORK_LINE})
if(LINUX)
    #rhel 8 libstdc++ unexported symbol workaround
    check_cxx_symbol_exists(std::filesystem::exists "filesystem" HAVE_FS)
    if(NOT HAVE_FS)
    target_link_libraries(RcppPlanc PRIVATE stdc++fs)
    endif()
endif()
target_link_libraries(RcppPlanc PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} HighFive hw_detect)
#target_include_directories(RcppPlanc BEFORE PRIVATE ${HighFive_BINARY_DIR}/include)
if(OpenMP_CXX_FOUND)
target_compile_options(RcppPlanc PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(RcppPlanc PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
target_link_options(RcppPlanc PRIVATE ${OpenMP_CXX_FLAGS})
endif()
endif()

#check_cxx_compiler_flag("-funsafe-math-optimizations" UNSAFE_MATH)
#if(UNSAFE_MATH)
#    target_compile_options(RcppPlanc PRIVATE "-funsafe-math-optimizations")
#endif()
check_cxx_compiler_flag("-ffinite-math-only" FINITE_MATH)
if(FINITE_MATH)
    target_compile_options(RcppPlanc PRIVATE "-ffinite-math-only")
endif()
check_cxx_compiler_flag("-fno-signaling-nans" SIGNALING_NANS)
if(SIGNALING_NANS)
    target_compile_options(RcppPlanc PRIVATE "-fno-signaling-nans")
endif()
check_cxx_compiler_flag("-fno-rounding-math" ROUNDING_MATH)
if(SIGNALING_NANS)
    target_compile_options(RcppPlanc PRIVATE "-fno-rounding-math")
endif()
check_cxx_compiler_flag("-fno-trapping-math" TRAPPING_MATH)
if(TRAPPING_MATH)
    target_compile_options(RcppPlanc PRIVATE "-fno-trapping-math")
endif()
check_cxx_compiler_flag("-fno-math-errno" ERRNO_MATH)
if(ERRNO_MATH)
    target_compile_options(RcppPlanc PRIVATE "-fno-math-errno")
endif()
check_cxx_compiler_flag("-fno-signed-zeros" UNSIGNED_ZERO)
if(UNSIGNED_ZERO)
    target_compile_options(RcppPlanc PRIVATE "-fno-signed-zeros")
endif()
#check_cxx_compiler_flag("-freciprocal-math" RECIPROCAL_MATH)
#if(RECIPROCAL_MATH)
#    target_compile_options(RcppPlanc PRIVATE "-freciprocal-math")
#endif()
list(JOIN BLAS_LIBRARIES " " BLAS_LIBRARY_STRING)
string(TOUPPER ${BLAS_LIBRARY_STRING} BLAS_LIBRARY_STRING)
#if(BLAS_LIBRARIES MATCHES "")
#target_compile_definitions(RcppPlanc PRIVATE BLAS_IMPLICIT)
#endif()
if(BLAS_LIBRARY_STRING MATCHES ^.*OPENBLAS.*)
    check_include_file("openblas/cblas.h" OPENBLAS_HEADER)
    if(OPENBLAS_HEADER)
        target_compile_definitions(RcppPlanc PRIVATE "HAVE_OPENBLAS_CBLAS_H")
        set(USING_OPENBLAS TRUE)
        set(OPENBLAS_DEFINED "HAVE_OPENBLAS_CBLAS_H")
    endif()
endif()
if(APPLE)
target_include_directories(RcppPlanc PRIVATE "${BLAS_LIBRARIES}/Frameworks/vecLib.framework/Headers")
endif()
if(NOT APPLE AND NOT USE_MKL AND NOT USING_OPENBLAS)
        set(CBLAS_H_DIR CACHE PATH "Directory containing cblas.h")
if(CBLAS_H_DIR)
target_include_directories(RcppPlanc PRIVATE ${CBLAS_H_DIR})
else()
    check_include_file("flexiblas/cblas.h" FLEXIBLAS_HEADER)
    if(FLEXIBLAS_HEADER)
        target_compile_definitions(RcppPlanc PRIVATE "HAVE_FLEXIBLAS_CBLAS_H")
    else()
        check_include_file("cblas.h" CBLAS_HEADER)
        if(NOT CBLAS_HEADER AND NOT DEFINED BLAS_IMPLICIT)
            message(FATAL_ERROR "no cblas header, please set include directory with CBLAS_H_DIR variable")
        elseif(BLAS_LIBRARY_STRING MATCHES ^.*LIBBLAS.*)
            set(CMAKE_REQUIRED_LIBRARIES blas)
            check_symbol_exists(cblas_sgemm "cblas.h" HAS_CBLAS_SYMBOLS)
            if(NOT HAS_CBLAS_SYMBOLS)
                find_library(CBLAS_REFERENCE_LIB cblas libcblas REQUIRED)
                target_link_libraries(RcppPlanc PRIVATE ${CBLAS_REFERENCE_LIB})
            endif()
        endif()
    endif()
endif()
endif()

#target_compile_definitions(RcppPlanc PRIVATE ARMA_64BIT_WORD)

add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${HighFive_BINARY_DIR}/CMakeFiles/hdf5/cmake_hdf5_test.c
        VERBATIM)
if(WIN32)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${HighFive_SOURCE_DIR}/Makefile
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${HighFive_BINARY_DIR}/CMakeFiles/hdf5/cmake_hdf5_test.c
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${HighFive_SOURCE_DIR}/include/highfive/H5Version.hpp
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${HighFive_BINARY_DIR}/include/highfive/H5Version.hpp
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/CMakeFiles/CheckTypeSize/SIZE_OF_VOIDP.c
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/CMakeFiles/FindOpenMP/OpenMPCheckVersion.c
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/FindOpenMP/OpenMPCheckVersion.cpp
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/FindOpenMP/OpenMPTryFlag.c
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/FindOpenMP/OpenMPTryFlag.cpp
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/${CMAKE_VERSION}/CompilerIdC/CMakeCCompilerId.c
        VERBATIM)
add_custom_command(TARGET RcppPlanc POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_SOURCE_DIR}/../src/${CMAKE_VERSION}/CompilerIdCXX/CMakeCXXCompilerId.cpp
        VERBATIM)
endif()
